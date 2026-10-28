include(joinpath(@__DIR__,"independent_fd.jl"))
using Test, TOML, SHA
function run_checks(output)
    OUTPUT=output
    ispath(OUTPUT) && error("refusing overwrite")
    mkpath(OUTPUT)
    open(joinpath(OUTPUT,"inputs.toml"),"w") do io
        root=normpath(joinpath(@__DIR__,"../../.."))
        files=[joinpath(@__DIR__,name) for name in ("fd_checks.jl","independent_fd.jl","../../../cli.jl")]
        TOML.print(io,Dict("suite"=>"independent-FD-MMS","julia_version"=>string(VERSION),
            "epsilon"=>.8,"beta"=>1.,"physical_validation"=>false,
            "source_sha256"=>Dict(relpath(path,root)=>bytes2hex(sha256(read(path))) for path in files));sorted=true)
    end
    @test !isdefined(@__MODULE__,:assemble_full_film)
    @test !isdefined(@__MODULE__,:full_film_fluxes)

    @testset "Independent analytical end-flux reconstruction" begin
        ecc=.8;beta=1.
        exact=pi/(3beta)*(1+1.5ecc^2-1.5ecc-.375ecc^3)
        for nt in (32,64,128)
            nz=nt÷2;c=fd_coefficients(nt,nz;ecc,beta)
            p=[sin(t/2)^2*(1-(w/beta)^2) for t in c.theta,w in c.axial]
            result=(;pressure=p,coefficients=c)
            obs=fd_observables(result;beta)
            @test obs.Q0_second ≈ exact rtol=1e-12
            @test obs.Q0_third ≈ exact rtol=1e-12
            @test obs.Q0_lower ≈ exact/2 rtol=1e-12
            @test obs.Q0_upper ≈ exact/2 rtol=1e-12
            @test obs.Q0_first ≈ exact*(1-1/nz) rtol=1e-12
        end
    end

    @testset "Independent finite-width manufactured PDE and side flow" begin
        ecc=.8
        forcing(t,w)=-(cos(t)/2-3ecc*sin(t)^2/(2*(1+ecc*cos(t))))*(1-w^2)+2sin(t/2)^2
        errors=Float64[];flow_errors=Float64[]
        exact_flow=pi/3*(1+1.5ecc^2-1.5ecc-.375ecc^3)
        open(joinpath(OUTPUT,"fd-mms.tsv"),"w") do io
            println(io,"nt\tnz\tpressure_relative_L2\tQ0\texact_Q0\tflow_relative_error")
            for nt in (32,64,128,256)
                result=solve_fd(nt,nt÷2;ecc,forcing,obstacle=false,tolerance=1e-11)
                c=result.coefficients
                exact=[sin(t/2)^2*(1-w^2) for t in c.theta,w in c.axial]
                error=norm(result.pressure-exact)/norm(exact)
                obs=fd_observables(result)
                push!(errors,error);push!(flow_errors,abs(obs.Q0_second/exact_flow-1))
                println(io,join((nt,nt÷2,error,obs.Q0_second,exact_flow,last(flow_errors)),'\t'))
                @test result.residual.kkt<=1e-11
            end
        end
        @test all(errors[2:end] .< .3 .* errors[1:end-1])
        @test all(flow_errors[2:end] .< .3 .* flow_errors[1:end-1])
        @test errors[end]<1e-4
        @test flow_errors[end]<1e-4
    end

    @testset "Projected SOR tolerance and relaxation control" begin
        first=solve_fd(64,32;tolerance=1e-11,relaxation=1.)
        second=solve_fd(64,32;tolerance=1e-11,relaxation=1.7)
        @test minimum(first.pressure)>=0
        @test norm(first.pressure-second.pressure)/norm(first.pressure)<1e-8
        @test abs(fd_observables(first).Q0_second/fd_observables(second).Q0_second-1)<1e-8
        finer=solve_fd(128,64;tolerance=1e-11)
        warm=solve_fd(128,64;initial=doubled_initial(second.pressure),tolerance=1e-11)
        @test norm(finer.pressure-warm.pressure)/norm(finer.pressure)<1e-8
        @test all(warm.pressure[[1,end],:].==0)
        @test all(warm.pressure[:,[1,end]].==0)
    end
end

include(joinpath(@__DIR__, "../../../cli.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    exit(BearingCLI.study_main(; description="Run analytical end-flow checks.") do cfg
        run_checks(cfg.output)
    end)
end
