# Case-specific conical full-film refinement run.
include(joinpath(@__DIR__,"../../reynolds/src/conical.jl"))
include(joinpath(@__DIR__,"../../../experiment.jl"))

function run_conical_verification()
    case = load_experiment("conical-unported-e0p6")
    p = case.p
    case.variant == "unported" || error("verification requires the unported case")
    # This script documents one fixed geometry and must not silently reuse another.
    (p.journal_midpoint_radius,p.bearing_length,p.radial_clearance,
     p.semicone_angle_deg,p.eccentricity_ratio,case.ntheta,case.nz) ==
        (50.0,100.0,0.05,10.0,0.6,512,256) || error("unexpected verification geometry/grid")
    dir = joinpath(case.dir,"runs","full-film-verification-vg32-40c")
    ispath(dir) && error("refusing to overwrite existing verification run: $dir")
    mkpath(dir)
    metadata = Dict{String,Any}(
        "completed" => false, "julia_version" => string(VERSION),
        "model" => "rigid-conical-signed-full-film-verification",
        "case" => "conical-unported-e0p6", "mesh_sha256" => mesh_hash(case.mesh),
        "rpm" => 2000.0, "viscosity_pa_s" => 32e-6*860.0,
        "nominal_temperature_c" => 40.0, "assumed_density_kg_m3" => 860.0,
        "nominal_kinematic_viscosity_m2_s" => 32e-6,
        "end_pressure_gauge_pa" => [0.0,0.0], "cavitation" => false,
        "load_control" => false, "elasticity" => false,
        "gap" => "radial-ray-gap times cos(cone angle), thin-film projection",
        "force" => "pressure force on housing conical surface; excludes ends and shear",
        "warning" => "Not a physical bearing prediction or Gangrade reproduction",
        "source_sha256" => Dict(relpath(path,joinpath(@__DIR__,"../../..")) => mesh_hash(path)
            for path in (joinpath(@__DIR__,"../../../cli.jl"),
                         joinpath(@__DIR__,"../../reynolds/src/conical.jl"), @__FILE__,
                         joinpath(@__DIR__,"../../reynolds/src/cylindrical.jl"),
                         joinpath(@__DIR__,"../../reynolds/src/full_film.jl"),
                         joinpath(@__DIR__,"../../reynolds/src/face_fluxes.jl"),
                         joinpath(@__DIR__,"../../reynolds/src/fitted_full_film.jl"),
                         joinpath(@__DIR__,"../../reynolds/src/pressure_vtk.jl"),
                         joinpath(@__DIR__,"../../../bearing.jl"),
                         joinpath(@__DIR__,"../../../mesh/2d/grid.jl"),
                         joinpath(@__DIR__,"../../../mesh/2d/quality.jl"),
                         joinpath(@__DIR__,"../../../mesh/2d/metrics.jl"),
                         joinpath(@__DIR__,"../../../mesh/vtk.jl"))))
    save_metadata() = open(joinpath(dir,"verification.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    save_metadata()
    previous = nothing
    open(joinpath(dir,"convergence.tsv"),"w") do io
        println(io,join(("ntheta","nz","min_pressure_pa","max_pressure_pa",
            "housing_fx_n","housing_fy_n","housing_fz_n","absolute_end_flow_m3_s",
            "relative_l1_cell_imbalance","relative_global_imbalance",
            "relative_pressure_mismatch"),'\t'))
        for (nt,nz) in ((32,16),(64,32),(128,64),(256,128),(512,256))
            r = solve_conical_full_film(p,nt,nz;rpm=2000.0,viscosity_pa_s=32e-6*860.0)
            mismatch = NaN
            if !isnothing(previous)
                pc = reshape(previous.pressure_pa,nt÷2,nz÷2)
                pf,af = reshape(r.pressure_pa,nt,nz),reshape(r.grid.areas,nt,nz)
                average = [sum(pf[2i-1:2i,2j-1:2j].*af[2i-1:2i,2j-1:2j]) /
                           sum(af[2i-1:2i,2j-1:2j]) for i in 1:nt÷2,j in 1:nz÷2]
                mismatch = norm(pc-average)/norm(average)
            end
            d = r.diagnostics
            println(io,join((nt,nz,minimum(r.pressure_pa),maximum(r.pressure_pa),
                r.housing_pressure_force_n...,d.absolute_end_flow_m3_s,
                d.relative_l1_cell_imbalance,d.relative_global_imbalance,mismatch),'\t'))
            println("$(nt)×$(nz): pressure mismatch=$mismatch; continuity=$(d.relative_l1_cell_imbalance)")
            previous = r
        end
    end
    write_conical_pressure_vtk(joinpath(dir,"full-film.vtk"),previous)
    metadata["diagnostics"] = Dict(string(k) => (v isa Tuple ? collect(v) : v)
                                    for (k,v) in pairs(previous.diagnostics))
    metadata["housing_pressure_force_n"] = collect(previous.housing_pressure_force_n)
    metadata["completed"] = true
    save_metadata()
    println("Saved signed-full-film numerical verification to $dir")
    return dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "../../../cli.jl"))
    exit(BearingCLI.study_main(; output=false,
        description="Refine the saved conical-unported-e0p6 case; output stays under its runs directory.") do cfg
        run_conical_verification()
    end)
end
