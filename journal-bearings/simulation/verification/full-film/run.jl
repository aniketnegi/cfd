include(joinpath(@__DIR__,"benchmark.jl"))
include(joinpath(@__DIR__,"donaldson_series.jl"))

function run_benchmark(output;full=false)
    table,header=readdlm(joinpath(@__DIR__,"table-6-5.csv"),',';header=true)
    eccentricities=parse.(Float64,String.(vec(header)[2:end]))
    points=[(row[1],ecc,row[column+1]) for row in eachrow(table) for (column,ecc) in enumerate(eccentricities)]
    full || filter!(point->(point[1],point[2]) in ((.7,.2),(1.,.8),(.1,.99)),points)
    fresh_output(output;suite=full ? "full-film-all320" : "full-film-representative")
    open(joinpath(output,"results.tsv"),"w") do io
        first=true
        for (aspect,ecc,target) in points
            basis=spectral_basis(ecc,320)
            reference=spectral_load(basis,aspect)
            @assert abs(reference-spectral_load(spectral_basis(ecc,160),aspect))<1e-9
            historical=series_observables(ecc,donaldson_basis(ecc,15),aspect).ratio
            for model in ("classical","production"),nt in (128,256,512)
                nz=refinement_axial_count(aspect,nt)
                p=model=="classical" ? DonaldsonCylinder(50.,100aspect,.05,ecc,0.,0.) :
                    BearingParameters(journal_midpoint_radius=50.,bearing_length=100aspect,
                        radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=ecc)
                grid=model=="classical" ? benchmark_grid(p,nt,nz) : housing_grid(p,nt,nz)
                elapsed=@elapsed result=solve_full_film(grid;rpm=RPM,viscosity_pa_s=MU)
                infinity=model=="classical" ? infinite_load(p) : production_infinite_load(p)
                force=result.housing_pressure_force_n
                ratio=hypot(force[1],force[2])/infinity
                d=result.diagnostics
                @assert d.relative_l1_cell_imbalance<1e-10 && d.relative_global_imbalance<1e-10
                row=(;model,L_D=aspect,epsilon=ecc,nt,nz,published=target,load_ratio=ratio,
                    classical_spectral_load=reference,donaldson_N15_load=historical,
                    relative_cell_imbalance=d.relative_l1_cell_imbalance,
                    relative_global_imbalance=d.relative_global_imbalance,seconds=elapsed)
                first && println(io,join(keys(row),'\t'));first=false
                println(io,join(values(row),'\t'));flush(io)
                model=="classical" && nt==512 && save_profiles(p,result,"LD$(aspect)-eps$(ecc)";directory=joinpath(output,"profiles"))
                @printf("%s LD=%g eps=%g %dx%d ratio=%.10f\n",model,aspect,ecc,nt,nz,ratio);flush(stdout)
                GC.gc()
            end
        end
    end
end
include(joinpath(@__DIR__, "../../../cli.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    exit(BearingCLI.study_main(; full=true,
        description="Compare signed full-film loads with Donaldson and independent spectral references.") do cfg
        run_benchmark(cfg.output; full=cfg.full)
    end)
end
