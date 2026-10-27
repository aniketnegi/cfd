# Numerical verification with a broad artificial reservoir, not a rig port.
using TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__,"conical_jfo.jl"))
include(joinpath(@__DIR__,"conical_jfo_mms.jl"))

function conical_jfo_numerical_case(nt,nz)
    p = BearingParameters(journal_midpoint_radius=50.0,bearing_length=100.0,
        radial_clearance=0.05,semicone_angle_deg=10.0,eccentricity_ratio=0.6)
    mask = vec([pi/4<(i-0.5)*2pi/nt<3pi/4 && 0.375<(j-0.5)/nz<0.625
                for i in 1:nt,j in 1:nz])
    return solve_conical_jfo(p,nt,nz;rpm=2000.0,viscosity_pa_s=0.02752,
        end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0,
        feed_cells=mask,feed_pressure_pa=5e5)
end

function conical_jfo_nested_error(coarse,fine,nt,nz)
    weights = reshape(fine.grid.areas,2nt,2nz)
    function average(field)
        f = reshape(field,2nt,2nz)
        [sum(f[2i-1:2i,2j-1:2j].*weights[2i-1:2i,2j-1:2j]) /
         sum(weights[2i-1:2i,2j-1:2j]) for i in 1:nt,j in 1:nz]
    end
    p,f = vec(average(fine.pressure_pa)),vec(average(fine.fill_fraction))
    return (; pressure_relative_l2=norm(coarse.pressure_pa-p)/norm(p),
        fill_area_mean_absolute=sum(abs.(coarse.fill_fraction-f).*coarse.grid.areas)/sum(coarse.grid.areas))
end

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
    root = normpath(joinpath(@__DIR__,"../.."))
    paths = ["simulation/reynolds/"*f for f in
        ("conical_jfo.jl","jfo_active_set.jl","conical_jfo_mms.jl","conical_jfo_checks.jl",
         "conical_jfo_verification.jl","conical.jl","cylindrical.jl","jfo_slider.jl",
         "jfo_slider_reference.jl","jfo_slider_checks.jl")]
    append!(paths,["bearing.jl","mesh/2d/grid.jl","expts/conical-jfo-verification/plot.py"])
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
    out = isempty(ARGS) ? joinpath(@__DIR__,"../../out/cases/conical-jfo-numerical/runs/transport-verification") : ARGS[1]
    verify_conical_jfo(normpath(out))
end
