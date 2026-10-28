# Case-specific evidence runner, deliberately separate from the bearing run CLI.
# Optional sole argument: a fresh output directory. Existing runs are refused.
using TOML, SHA, Printf, LinearAlgebra
include(joinpath(@__DIR__,"../../../reynolds/src/jfo_slider.jl"))
include(joinpath(@__DIR__,"jfo_slider_reference.jl"))

function verify_jfo_slider(out)
    ispath(out) && error("refusing to overwrite existing evidence: $out")
    ref = slider_reference()
    resolutions = [40,80,160,320,640,1280]
    results = [solve_jfo_slider(slider_gap,n;length_m=0.125,viscosity_pa_s=0.015,
        entrainment_m_s=2.0,end_pressure_pa=(1e6,1e6),cavitation_pressure_pa=0.0,
        inflow_fill_fraction=1.0) for n in resolutions]
    mkpath(out)
    open(joinpath(out,"convergence.tsv"),"w") do io
        println(io,"cells\tdx_m\tpressure_relative_l2\tfill_mean_absolute_error\tpeak_pa\tflux_relative_error\trupture_m\treformation_m\tmax_cell_relative_imbalance\tboundary_relative_imbalance\tflux_relative_spread\titerations")
        for (n,r) in zip(resolutions,results)
            p,f = ref.pressure.(r.x_m),ref.fill.(r.x_m)
            cavity = findall(r.fill_fraction .< 1-1e-8)
            row = [n,0.125/n,norm(r.pressure_pa-p)/norm(p),sum(abs,r.fill_fraction-f)/n,
                maximum(r.pressure_pa),abs(r.face_flux_m2_s[1]/ref.constant_flux_m2_s-1),
                r.face_x_m[first(cavity)],r.face_x_m[last(cavity)+1],
                r.diagnostics.relative_max_cell_imbalance,r.diagnostics.relative_boundary_imbalance,
                r.diagnostics.relative_flux_spread,r.diagnostics.iterations]
            println(io,join(row,'\t'))
        end
    end
    r = results[end]
    open(joinpath(out,"profiles.tsv"),"w") do io
        println(io,"x_m\tpressure_pa\treference_pressure_pa\tfill_fraction\treference_fill_fraction")
        for i in eachindex(r.x_m)
            x = r.x_m[i]
            println(io,join((x,r.pressure_pa[i],ref.pressure(x),r.fill_fraction[i],ref.fill(x)),'\t'))
        end
    end
    open(joinpath(out,"face-flux.tsv"),"w") do io
        println(io,"x_m\tflux_m2_s\treference_flux_m2_s")
        for (x,q) in zip(r.face_x_m,r.face_flux_m2_s)
            println(io,join((x,q,ref.constant_flux_m2_s),'\t'))
        end
    end
    root = normpath(joinpath(@__DIR__,"../../../.."))
    source_paths = ["simulation/reynolds/src/jfo_slider.jl",
        "simulation/reynolds/src/jfo_active_set.jl",
        "simulation/verification/jfo/slider/jfo_slider_reference.jl",
        "simulation/verification/jfo/slider/jfo_slider_verification.jl",
        "simulation/tests/jfo_slider_checks.jl","Project.toml","Manifest.toml"]
    hashes = Dict(p=>bytes2hex(sha256(read(joinpath(root,p)))) for p in source_paths)
    metadata = Dict(
        "scope"=>"Steady rigid 1D mass-conserving cavitation verification; NOT conical/EHD/experimental validation",
        "julia_version"=>string(VERSION),"source_sha256"=>hashes,
        "source"=>Dict("doi"=>"10.1016/j.triboint.2013.12.012",
            "section"=>"4.1, with author supplement version 1.3 (2014-02-01)",
            "paper_url"=>"https://bluebox.ippt.pan.pl/~sstupkie/papers/Lengiewicz_etal_TribInt_2014_accepted.pdf",
            "supplement_url"=>"https://bluebox.ippt.pan.pl/~sstupkie/files/bearing.pdf"),
        "inputs"=>Dict("length_m"=>0.125,"mean_gap_m"=>20e-6,"gap_amplitude_m"=>5e-6,
            "viscosity_pa_s"=>0.015,"entrainment_m_s"=>2.0,
            "end_pressure_above_cavitation_pa"=>[1e6,1e6],"cavitation_pressure_pa"=>0.0,
            "inflow_fill_fraction"=>1.0,"pressure_reference"=>"relative to benchmark cavitation pressure; not an oil vapour-pressure claim",
            "cells"=>resolutions),
        "continuum_reference"=>Dict("rupture_m"=>ref.rupture_m,"reformation_m"=>ref.reformation_m,
            "peak_pressure_pa"=>ref.peak_pressure_pa,"pressure_at_zero_pa"=>ref.pressure_at_zero_pa,
            "constant_flux_m2_s"=>ref.constant_flux_m2_s,"hstar_m"=>ref.hstar_m,
            "method"=>"Constant-flux integration, 64-point Gauss-Legendre quadrature and bracketed bisection; checked against printed supplement outputs"),
        "finest_diagnostics"=>Dict(string(k)=>v for (k,v) in pairs(r.diagnostics)))
    open(joinpath(out,"verification.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    @printf("Saved slider evidence to %s\nReference peak %.9f MPa, rupture %.9f mm, reformation %.9f mm\n",
        out,ref.peak_pressure_pa/1e6,ref.rupture_m*1e3,ref.reformation_m*1e3)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("expected at most one fresh output directory")
    out = isempty(ARGS) ? joinpath(@__DIR__,"../../../../out/cases/lengiewicz-sinusoidal-slider/runs/jfo-verification") : ARGS[1]
    verify_jfo_slider(normpath(out))
end
