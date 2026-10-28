include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"independent_fd.jl"))
function compare_fd(output;full=false)
    fresh_output(output;suite=full ? "independent-FD-FV-1024" : "independent-FD-FV-representative")
    grids=full ? (64,128,256,512,1024) : (64,128,256)
    warm=nothing;last_fd=nothing;errors=Float64[]
    @testset "Independent finite-width pressure, leakage and rupture comparison" begin
        open(joinpath(output,"fd-comparison.tsv"),"w") do io
            first=true
            for nt in grids
                nz=nt÷2
                fd=solve_fd(nt,nz;initial=warm,tolerance=1e-9)
                @test fd.residual.kkt<=1e-9
                obs=fd_observables(fd)
                native,field=cylindrical_case(1.,.8,"classical",nt)
                scale=.02752*(2000*2pi/60)*.05^2/(50e-6)^2
                fv=reshape(field.pressure_pa,nt,nz)/scale
                sampled=(fd.pressure[2:end,1:end-1]+fd.pressure[2:end,2:end])/2
                difference=fv-sampled
                error=norm(difference)/norm(sampled);push!(errors,error)
                theta=field.theta_centres;h=1 .+ .8 .* cos.(theta)
                @test maximum(abs,theta-fd.coefficients.theta[2:end])<1e-12
                gradient_fd,_=end_gradients(fd.pressure,2/nz,2)
                gradient_fv=fv[:,1]*nz
                density_fd=h.^3 .* gradient_fd[2:end]/12
                density_fv=h.^3 .* gradient_fv/12
                row=(;nt,nz,FD_Q0=obs.Q0_second,FD_Q0_first=obs.Q0_first,FD_Q0_third=obs.Q0_third,
                    FV_Q0=native.flow_Q0,FD_lower_Q0=obs.Q0_lower,FD_upper_Q0=obs.Q0_upper,
                    pressure_relative_L2=error,pressure_relative_Linf=maximum(abs,difference)/maximum(abs,sampled),
                    end_gradient_relative_L2=norm(gradient_fv-gradient_fd[2:end])/norm(gradient_fd[2:end]),
                    end_flux_density_relative_L2=norm(density_fv-density_fd)/norm(density_fd),
                    FD_kkt=fd.residual.kkt,FD_sweeps=fd.sweeps)
                row!(io,row;header=first);first=false
                threshold=1e-8*max(maximum(fv),maximum(sampled))
                open(joinpath(output,"rupture-n$(nt).tsv"),"w") do profile
                    println(profile,"z_m\tFV_last_positive_rad\tFD_interpolated_last_positive_rad")
                    for j in 1:nz
                        println(profile,join((field.z_centres[j],theta[findlast(>(threshold),fv[:,j])],
                            theta[findlast(>(threshold),sampled[:,j])]),'\t'))
                    end
                end
                open(joinpath(output,"profiles-n$(nt).tsv"),"w") do profile
                    println(profile,"theta_rad\tFV_Pi_midplane\tFD_Pi_midplane\tFV_lower_flux_density\tFD_lower_flux_density")
                    for i in 1:nt
                        println(profile,join((theta[i],(fv[i,nz÷2]+fv[i,nz÷2+1])/2,
                            fd.pressure[i+1,nz÷2+1],density_fv[i],density_fd[i]),'\t'))
                    end
                end
                serialize(joinpath(output,"fd-n$(nt).jls"),fd.pressure)
                last_fd=fd
                warm=nt==last(grids) ? nothing : doubled_initial(fd.pressure)
                @printf("nt=%d FDQ0=%.9f FVQ0=%.9f pressureL2=%.6g\n",nt,obs.Q0_second,native.flow_Q0,error);flush(stdout)
                GC.gc()
            end
            @test all(errors[2:end] .< .3 .* errors[1:end-1])
            @test errors[end]<.001
        end
        if full
            original=fd_observables(last_fd).Q0_second
            open(joinpath(output,"tolerance.tsv"),"w") do io
                println(io,"relaxation\tkkt\tflow_relative_change")
                for relaxation in (1.8,1.95)
                    check=solve_fd(1024,512;initial=last_fd.pressure,relaxation,tolerance=1e-10)
                    change=abs(fd_observables(check).Q0_second/original-1)
                    @test check.residual.kkt<1e-10
                    @test change<1e-7
                    println(io,join((relaxation,check.residual.kkt,change),'\t'))
                end
            end
        end
    end
end
include(joinpath(@__DIR__, "../../../cli.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    exit(BearingCLI.study_main(; full=true,
        description="Compare independent finite differences with finite-volume pressure and end flows.") do cfg
        compare_fd(cfg.output; full=cfg.full)
    end)
end
