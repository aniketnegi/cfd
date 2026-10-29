using Test
include(joinpath(@__DIR__,"run.jl"))

function check_references()
    @testset "Cylindrical continuous forcing and independent FD" begin
        for direction in (-1,1),t in (.2,1.4,2.8,4.2),z in (-.7,.3)
            delta=1e-4
            P(t,z)=cylinder_exact(t,z,1.,.6;direction).P
            F(t,z)=cylinder_exact(t,z,1.,.6;direction).F
            H(t)=1-.6cos(t)
            qt(t,z)=direction*H(t)*F(t,z)/2-H(t)^3*(P(t+delta,z)-P(t-delta,z))/(24delta)
            qz(t,z)=-H(t)^3*(P(t,z+delta)-P(t,z-delta))/(24delta)
            numerical=(qt(t+delta,z)-qt(t-delta,z)+qz(t,z+delta)-qz(t,z-delta))/(2delta)
            exact=cylinder_exact(t,z,1.,.6;direction).source
            @test numerical≈exact atol=2e-6 rtol=2e-6
        end
        errors=Float64[]
        for n in (32,64,128)
            a=CylindricalJFOReference.assemble(1.,.6,-1.,n,n÷2)
            # Signed continuous manufactured field; only the FD pressure operator.
            rhs=[cos(t)*(1-z^2+2)+1.8sin(t)^2/(1-.6cos(t))*(1-z^2)
                for t in a.theta,z in a.zeta[2:end-1]]
            exact=[cos(t)*(1-z^2) for t in a.theta,z in a.zeta[2:end-1]]
            computed=a.A\vec(rhs)
            push!(errors,norm(computed-vec(exact))/norm(exact))
        end
        @test 1.95<log2(errors[1]/errors[2])<2.05
        @test 1.95<log2(errors[2]/errors[3])<2.05
        for pv in (-1.,-.1),direction in (-1,1)
            a=CylindricalJFOReference.assemble(1.,.6,pv,32,16;direction)
            sor=CylindricalJFOReference.solve(a;tolerance=1e-11)
            state=solve_jfo_active_set(a.A,a.B,a.b,1.)
            @test norm(state.pressure-vec(sor.P[:,2:end-1]))/norm(state.pressure)<1e-7
            @test maximum(abs,state.void-vec(1 .-sor.F[:,2:end-1]))<1e-7
            @test 0<=minimum(sor.F)<=maximum(sor.F)<=1
            @test sor.residual<1e-11
        end
        for graded in (false,true),rpm in (-RPM,RPM)
            coarse=manufactured_case(32,16;graded,rpm)
            fine=manufactured_case(64,32;graded,rpm)
            @test fine.pressure_relative_L2<coarse.pressure_relative_L2
            @test fine.fill_mean_absolute<coarse.fill_mean_absolute
            @test fine.local_balance<1e-10
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;output=false,
        description="Check cylindrical JFO forcing, FD consistency and independent algebra.") do cfg
        check_references()
    end)
end
