# Numerical verification with a broad artificial reservoir, not a rig port.
using TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"cases.jl"))

function verify_conical_jfo(out)
    ispath(out) && error("refusing to overwrite existing evidence: $out")
    grids = [32,64,128,256]
    runs = [conical_jfo_numerical_case(n,n÷2) for n in grids]
    mms = [conical_jfo_mms(n,n÷2;rpm) for rpm in (-2000.0,2000.0) for n in grids]
    mkpath(out)
    open(joinpath(out,"refinement.tsv"),"w") do io
        println(io,"nt\tnz\tpeak_pa\tminimum_fill\tfeed_m3_s\tpressure_mismatch_previous\tfill_mismatch_previous\tmax_cell_relative_imbalance\tglobal_relative_imbalance\titerations")
        for (k,n) in enumerate(grids)
            r = runs[k]
            e = k==1 ? (;pressure_relative_l2=NaN,fill_area_mean_absolute=NaN) :
                conical_jfo_nested_error(runs[k-1],r,n÷2,n÷4)
            println(io,join((n,n÷2,maximum(r.pressure_pa),minimum(r.fill_fraction),
                r.diagnostics.net_feed_flow_m3_s,e.pressure_relative_l2,e.fill_area_mean_absolute,
                r.diagnostics.relative_max_cell_imbalance,r.diagnostics.relative_global_imbalance,
                r.diagnostics.iterations),'\t'))
        end
    end
    open(joinpath(out,"manufactured.tsv"),"w") do io
        println(io,join(string.(keys(mms[1])),'\t'))
        for r in mms
            println(io,join(values(r),'\t'))
        end
    end
    n = grids[end]; nz = n÷2; r = runs[end]
    write_conical_jfo_vtk(joinpath(out,"cavitation.vtk"),r)
    open(joinpath(out,"fields.tsv"),"w") do io
        println(io,"theta_rad\tz_m\tpressure_pa\tfill_fraction\tfeed_cell\tfeed_source_m3_s")
        for j in 1:nz,i in 1:n
            k = i+n*(j-1)
            println(io,join(((i-0.5)*2pi/n,(j-0.5)*0.1/nz,r.pressure_pa[k],
                r.fill_fraction[k],Int(r.feed_cells[k]),r.feed_source_m3_s[k]),'\t'))
        end
    end
    root = normpath(joinpath(@__DIR__,"../../../.."))
    paths = ["simulation/reynolds/src/"*f for f in
        ("conical_jfo.jl","jfo_active_set.jl","conical.jl","cylindrical.jl",
         "full_film.jl","pressure_vtk.jl")]
    append!(paths,["simulation/verification/jfo/conical/"*f for f in
        ("cases.jl","conical_jfo_mms.jl","conical_jfo_verification.jl")])
    append!(paths,["simulation/tests/conical_jfo_checks.jl",
        "bearing.jl","mesh/2d/grid.jl","mesh/2d/quality.jl",
        "mesh/2d/metrics.jl","mesh/vtk.jl","Project.toml","Manifest.toml"])
    hashes = Dict(p=>bytes2hex(sha256(read(joinpath(root,p)))) for p in paths)
    metadata = Dict("scope"=>"2D periodic cone JFO numerical verification; artificial reservoir, NOT rig prediction or EHD",
        "julia_version"=>string(VERSION),"source_sha256"=>hashes,
        "inputs"=>Dict("journal_midpoint_radius_m"=>0.05,"axial_length_m"=>0.1,
            "radial_clearance_m"=>50e-6,"semi_cone_angle_deg"=>10.0,"eccentricity_ratio"=>0.6,
            "rpm"=>2000.0,"viscosity_pa_s"=>0.02752,"end_pressure_pa"=>[0.0,0.0],
            "cavitation_pressure_pa"=>0.0,"feed_pressure_pa"=>5e5,
            "feed_theta_range_rad"=>[pi/4,3pi/4],"feed_z_range_m"=>[0.0375,0.0625],
            "feed_model"=>"fully liquid bidirectional prescribed-pressure cells; no passage or resistance",
            "pressure_reference"=>"arbitrary numerical zero, not a physical VG32 cavitation choice",
            "n_theta"=>grids,"n_axial"=>grids.÷2),
        "finest_diagnostics"=>Dict(string(k)=>(v isa Tuple ? collect(v) : v) for (k,v) in pairs(r.diagnostics)))
    open(joinpath(out,"verification.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    println("Saved 2D numerical cone evidence to $out (not a physical bearing prediction)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("expected at most one fresh output directory")
    out = isempty(ARGS) ? joinpath(@__DIR__,"../../../../out/cases/conical-jfo-numerical/runs/transport-verification") : ARGS[1]
    verify_conical_jfo(normpath(out))
end
