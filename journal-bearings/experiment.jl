# Shared, versioned experiment geometry manifest. The manifest is the source
# of truth for the mesh and bearing geometry; output VTK is not parsed for
# simulation inputs. Caller loads bearing.jl before these manifest helpers.
using TOML, SHA

const DEFAULT_CASES_ROOT = joinpath(@__DIR__, "out", "cases")

function safe_case_name(name::AbstractString)
    occursin(r"^[a-z0-9][a-z0-9_-]*$", name) ||
        throw(ArgumentError("name must match [a-z0-9][a-z0-9_-]* (no paths or dots): $name"))
    return String(name)
end

function case_directory(name::AbstractString; cases_root=DEFAULT_CASES_ROOT)
    dir = joinpath(cases_root, safe_case_name(name))
    islink(dir) && throw(ArgumentError("case directory must not be a symlink: $dir"))
    return dir
end

mesh_hash(path) = bytes2hex(open(sha256, path))

function experiment_manifest(name, variant, p, ntheta, nz, mesh_path;
                             feed_diameter=nothing, feed_z=nothing, mesh_grid=nothing)
    data = Dict{String,Any}(
        "schema_version" => 1, "name" => safe_case_name(name), "variant" => variant,
        "bearing" => Dict("journal_radius_mm" => p.journal_midpoint_radius,
                          "length_mm" => p.bearing_length,
                          "clearance_mm" => p.radial_clearance,
                          "cone_angle_deg" => p.semicone_angle_deg,
                          "eccentricity_ratio" => p.eccentricity_ratio),
        "grid" => Dict("n_theta" => ntheta, "n_axial" => nz,
                       "mesh_sha256" => mesh_hash(mesh_path)))
    if variant == "surface-feed"
        data["feed"] = Dict{String,Any}("diameter_mm" => feed_diameter, "z_mm" => feed_z)
    end
    if !isnothing(mesh_grid)
        fitted = hasproperty(mesh_grid, :feed_edges)
        data["grid"]["layout"] = fitted ? mesh_grid.layout :
            variant == "surface-feed" ? "cell-centre-mask" : "periodic-rectangular"
        data["grid"]["surface_cells"] = length(mesh_grid.quads)
        data["grid"]["nodes"] = length(mesh_grid.nodes)
        if fitted
            variant == "surface-feed" || throw(ArgumentError("fitted inlet requires surface-feed variant"))
            data["grid"]["feed_rim_edges"] = length(mesh_grid.feed_edges)
            data["grid"]["quality"] = Dict{String,Any}(string(k)=>v for (k,v) in pairs(mesh_grid.quality))
            data["feed"]["representation"] = "rim-boundary; reservoir interior omitted"
        end
    end
    return data
end

function save_experiment_manifest(dir, name, variant, p, ntheta, nz;
                                  feed_diameter=nothing, feed_z=nothing, mesh_grid=nothing)
    path = joinpath(dir, "experiment.toml")
    ispath(path) && throw(ArgumentError("experiment manifest already exists: $path"))
    data = experiment_manifest(name, variant, p, ntheta, nz, joinpath(dir, "mesh.vtk");
                               feed_diameter, feed_z, mesh_grid)
    tmp, io = mktemp(dir)
    try
        TOML.print(io, data; sorted=true)
        close(io)
        mv(tmp, path)  # no replacement if another manifest appeared
    finally
        isopen(io) && close(io)
        ispath(tmp) && rm(tmp)
    end
    return path
end

function load_experiment(name; cases_root=DEFAULT_CASES_ROOT)
    dir = case_directory(name; cases_root)
    path = joinpath(dir, "experiment.toml")
    isfile(path) || throw(ArgumentError("no experiment manifest at $path"))
    data = TOML.parsefile(path)
    version = get(data, "schema_version", nothing)
    version isa Int && !(version isa Bool) && version == 1 ||
        throw(ArgumentError("unsupported experiment schema version in $path"))
    get(data, "name", nothing) == name || throw(ArgumentError("experiment name mismatch: $path"))
    variant = get(data, "variant", nothing)
    variant in ("unported", "surface-feed") || throw(ArgumentError("unknown geometry variant"))
    bearing = get(data, "bearing", nothing)
    grid = get(data, "grid", nothing)
    bearing isa AbstractDict && grid isa AbstractDict ||
        throw(ArgumentError("missing bearing/grid in $path"))
    number(d, k) = begin
        v = get(d, k, nothing)
        v isa Real && !(v isa Bool) && isfinite(v) ||
            throw(ArgumentError("missing or invalid $k in $path"))
        Float64(v)
    end
    p = BearingParameters(journal_midpoint_radius=number(bearing, "journal_radius_mm"),
                          bearing_length=number(bearing, "length_mm"),
                          radial_clearance=number(bearing, "clearance_mm"),
                          semicone_angle_deg=number(bearing, "cone_angle_deg"),
                          eccentricity_ratio=number(bearing, "eccentricity_ratio"))
    ntheta, nz = get(grid, "n_theta", nothing), get(grid, "n_axial", nothing)
    ntheta isa Int && !(ntheta isa Bool) && ntheta >= 3 &&
        nz isa Int && !(nz isa Bool) && nz >= 1 ||
        throw(ArgumentError("invalid grid resolution in $path"))
    if variant == "surface-feed"
        feed = get(data, "feed", nothing)
        feed isa AbstractDict || throw(ArgumentError("missing feed parameters in $path"))
        diameter, z = number(feed, "diameter_mm"), number(feed, "z_mm")
        diameter > 0 && 0 < z < p.bearing_length ||
            throw(ArgumentError("invalid feed geometry in $path"))
    elseif haskey(data, "feed")
        throw(ArgumentError("unported experiment has feed parameters in $path"))
    end
    # Schema 1 defaults to a rectangular grid or cell-centre feed mask.
    layout = get(grid, "layout", variant == "surface-feed" ? "cell-centre-mask" : "periodic-rectangular")
    layout in ("periodic-rectangular", "cell-centre-mask", "controlled-fitted-land-v1") ||
        throw(ArgumentError("unsupported mesh layout in $path"))
    (variant == "unported") == (layout == "periodic-rectangular") ||
        throw(ArgumentError("mesh layout and geometry variant disagree in $path"))
    if variant == "surface-feed" && layout == "cell-centre-mask"
        !haskey(data["feed"], "representation") ||
            throw(ArgumentError("cell-centre masks must not declare a fitted rim representation"))
    end
    if layout == "controlled-fitted-land-v1"
        for key in ("surface_cells", "nodes", "feed_rim_edges")
            v = get(grid, key, nothing)
            v isa Int && !(v isa Bool) && v > 0 || throw(ArgumentError("invalid fitted $key in $path"))
        end
        get(data["feed"], "representation", nothing) == "rim-boundary; reservoir interior omitted" ||
            throw(ArgumentError("fitted inlet must declare its rim-boundary representation"))
        (ntheta,nz) in ((128,64),(256,128),(512,256)) ||
            throw(ArgumentError("unsupported controlled inlet background resolution"))
        level=ntheta÷128
        (grid["surface_cells"],grid["nodes"],grid["feed_rim_edges"]) ==
            (8416level^2,8416level^2+144level-1,32level) ||
            throw(ArgumentError("counts do not match controlled fitted layout"))
        (p.journal_midpoint_radius,p.bearing_length,p.radial_clearance,p.semicone_angle_deg,
         data["feed"]["diameter_mm"],data["feed"]["z_mm"]) == (50.,100.,.05,10.,5.,50.) ||
            throw(ArgumentError("parameters do not match controlled fitted reference geometry"))
    end
    mesh = joinpath(dir, "mesh.vtk")
    isfile(mesh) || throw(ArgumentError("missing case mesh: $mesh"))
    expected = get(grid, "mesh_sha256", nothing)
    expected isa String && mesh_hash(mesh) == expected ||
        throw(ArgumentError("mesh does not match experiment manifest: $mesh"))
    return (; dir, data, variant, layout, p, ntheta, nz, mesh)
end
