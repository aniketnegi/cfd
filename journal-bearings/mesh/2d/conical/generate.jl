# Command-line entry point for the conical housing surface grids (lengths in mm).
include(joinpath(@__DIR__, "../../../bearing.jl"))
include(joinpath(@__DIR__, "grid.jl"))
include(joinpath(@__DIR__, "surface_feed.jl"))
include(joinpath(@__DIR__, "../../vtk.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "experiment.jl"))

const USAGE = """NAME
    generate.jl - generate a conical housing surface grid

SYNOPSIS
    julia --project=. mesh/2d/conical/generate.jl [OPTION]...

DESCRIPTION
    Write a periodic Quad4 housing surface as VTK. All lengths are in mm.
    The feed variant fits the circular-drill rim with controlled quad blocks;
    it exports land cells plus boundary lines, not a stepped cell mask or 3D
    passage. It prescribes no supply pressure and solves no fluid equation.
    Fitted inlets use the supplied bearing dimensions, diameter and axial
    position, with drilling along +Y. The opening and its surrounding transition
    must fit inside the surface. Counts directly size the background;
    optimization runs at that resolution. Final counts differ after replacing
    the inlet patch. Requests failing the quality checks are rejected.
    Use --output FILE for fitted inlets. Named fitted-case manifests are not
    supported by the current external case loader.
    For unported meshes, --experiment NAME creates
    out/cases/NAME/{mesh.vtk,experiment.toml}; the name must be new.
    Options accept --key VALUE or --key=VALUE.

REQUIRED OPTIONS
    --geometry=GEOMETRY        unported or surface-feed
    --journal-radius-mm=MM     journal radius at the axial midpoint
    --length-mm=MM             axial bearing length
    --clearance-mm=MM          radial clearance
    --cone-angle-deg=DEG       cone semi-angle; 0 gives a cylinder
    --eccentricity-ratio=RATIO journal offset divided by radial clearance
    --experiment=NAME         new unported case name; exclusive with --output
    --output=FILE             standalone .vtk path; exclusive with --experiment

OPTIONS
    --n-theta=N               background circumferential count (default: 512)
    --n-axial=N               background axial count (default: 256)
                               fitted cells replace a patch in this background
    --feed-diameter-mm=MM      required for surface-feed only
    --feed-z-mm=MM             feed centre from the large axial end;
                               required for surface-feed only
    -h, --help                 display this help and exit

EXIT STATUS
    0    Success.
    2    Invalid arguments.
"""

const OPTIONS = Set(("--geometry", "--journal-radius-mm", "--length-mm",
                     "--clearance-mm", "--cone-angle-deg", "--eccentricity-ratio",
                     "--experiment", "--output", "--n-theta", "--n-axial",
                     "--feed-diameter-mm", "--feed-z-mm"))
const REQUIRED = ("--geometry", "--journal-radius-mm", "--length-mm",
                  "--clearance-mm", "--cone-angle-deg", "--eccentricity-ratio")

function parse_grid_args(args; cases_root=DEFAULT_CASES_ROOT)
    opts = Dict{String,String}()
    i = 1
    while i <= length(args)
        argument = args[i]
        key, value = if startswith(argument, "--") && occursin('=', argument)
            i += 1
            split(argument, '='; limit=2)
        else
            i < length(args) && !startswith(args[i + 1], "--") ||
                throw(ArgumentError("missing value for $argument"))
            i += 2
            (argument, args[i - 1])
        end
        key in OPTIONS || throw(ArgumentError("unknown option: $key"))
        haskey(opts, key) && throw(ArgumentError("repeated option: $key"))
        isempty(value) && throw(ArgumentError("missing value for $key"))
        opts[key] = value
    end
    for key in REQUIRED
        haskey(opts, key) || throw(ArgumentError("missing required option: $key"))
    end

    number(key) = begin
        value = tryparse(Float64, opts[key])
        isnothing(value) && throw(ArgumentError("invalid number for $key: $(opts[key])"))
        value
    end
    count(key, fallback) = begin
        value = tryparse(Int, get(opts, key, string(fallback)))
        isnothing(value) && throw(ArgumentError("invalid integer for $key"))
        value
    end

    geometry = opts["--geometry"]
    geometry in ("unported", "surface-feed") ||
        throw(ArgumentError("--geometry must be unported or surface-feed"))
    has_feed = haskey(opts, "--feed-diameter-mm") || haskey(opts, "--feed-z-mm")
    if geometry == "surface-feed"
        all(haskey(opts, key) for key in ("--feed-diameter-mm", "--feed-z-mm")) ||
            throw(ArgumentError("surface-feed requires --feed-diameter-mm and --feed-z-mm"))
    elseif has_feed
        throw(ArgumentError("feed options are only valid with surface-feed"))
    end

    haskey(opts, "--output") != haskey(opts, "--experiment") ||
        throw(ArgumentError("specify exactly one of --experiment or --output"))
    experiment = get(opts, "--experiment", nothing)
    geometry == "surface-feed" && !isnothing(experiment) && throw(ArgumentError(
        "named fitted cases require an updated case loader; use --output FILE",
    ))
    dir = isnothing(experiment) ? nothing : case_directory(experiment; cases_root)
    path = isnothing(experiment) ? opts["--output"] : joinpath(dir, "mesh.vtk")
    isnothing(experiment) && !endswith(lowercase(path), ".vtk") &&
        throw(ArgumentError("--output must end in .vtk"))
    existing = isnothing(dir) ? path : dir
    ispath(existing) && throw(ArgumentError("output already exists: $existing"))
    p = BearingParameters(
        journal_midpoint_radius=number("--journal-radius-mm"),
        bearing_length=number("--length-mm"),
        radial_clearance=number("--clearance-mm"),
        semicone_angle_deg=number("--cone-angle-deg"),
        eccentricity_ratio=number("--eccentricity-ratio"),
    )
    n_theta = count("--n-theta", 512)
    n_axial = count("--n-axial", 256)
    n_theta >= 3 || throw(ArgumentError("--n-theta must be at least 3"))
    n_axial >= 1 || throw(ArgumentError("--n-axial must be positive"))
    feed_diameter = geometry == "surface-feed" ? number("--feed-diameter-mm") : nothing
    feed_z = geometry == "surface-feed" ? number("--feed-z-mm") : nothing
    return (; geometry, p, n_theta, n_axial, feed_diameter, feed_z, path, experiment, dir)
end

function main(args=ARGS; cases_root=DEFAULT_CASES_ROOT)
    if "--help" in args || "-h" in args
        print(USAGE)
        return
    end
    cfg = parse_grid_args(args; cases_root)
    grid = if cfg.geometry == "surface-feed"
        surface_feed_grid(cfg.p, cfg.n_theta, cfg.n_axial;
                          feed_diameter=cfg.feed_diameter, feed_z=cfg.feed_z)
    else
        housing_grid(cfg.p, cfg.n_theta, cfg.n_axial)
    end
    if !isnothing(cfg.dir)
        mkpath(dirname(cfg.dir))
        mkdir(cfg.dir)  # exclusive; never overwrite a case
    end
    write_vtk(cfg.path, grid)
    if !isnothing(cfg.experiment)
        save_experiment_manifest(cfg.dir, cfg.experiment, cfg.geometry, cfg.p,
                                 cfg.n_theta, cfg.n_axial;
                                 feed_diameter=cfg.feed_diameter, feed_z=cfg.feed_z, mesh_grid=grid)
        println("Saved experiment: $(abspath(joinpath(cfg.dir, "experiment.toml")))")
    end
    if cfg.geometry == "surface-feed"
        println("Fitted inlet rim: $(length(grid.feed_edges)) edges; reservoir interior omitted")
    end
    println("Surface quality (NOT Hex8 OQ): $(grid.quality)")
    println("Mesh acceptance: $(grid.acceptance)")
    println("Wrote $(length(grid.quads)) quads to $(abspath(cfg.path))")
    println("Cell-centre radial gap: $(minimum(grid.radial_gaps)) to $(maximum(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch err
        if err isa ArgumentError
            println(stderr, "generate.jl: $(sprint(showerror, err))")
            println(stderr, "Try 'julia --project=. mesh/2d/conical/generate.jl --help' for more information.")
            exit(2)
        end
        rethrow()
    end
end
