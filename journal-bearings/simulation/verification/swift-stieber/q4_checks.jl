include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"q4.jl"))

function check_q4()
    @testset "Public mesh, nodal topology and Q4 energy" begin
        grid=housing_grid(cone(),12,4)
        a=assemble_q4(grid)
        @test a.grid===grid
        @test length(a.grid.nodes)==60
        @test length(a.grid.quads)==48
        @test a.grid.quads[1]==(1,2,14,13)
        @test count(a.fixed)==24
        @test size(a.A)==(60,60)
        @test norm(a.A-a.A')<1e-14*norm(a.A)
        @test norm(a.A*ones(60))<1e-14*norm(a.A)
        @test isposdef(Symmetric(Matrix(a.A[.!a.fixed,.!a.fixed])))
        r=solve_q4(grid)
        doubled=solve_q4(grid;viscosity_pa_s=2MU)
        @test doubled.pressure_pa≈2r.pressure_pa rtol=1e-11
        @test all(>=(0),r.pressure_pa)
        @test r.complementarity<1e-10
        @test q4_periodic_sample([0.,pi],[2.,4.],[0.,pi,2pi,pi/2])≈[2.,4.,2.,3.]
        for z in (.02,.05,.08)
            @test q4_plane(a,3 .+2a.theta.+5a.z;z_m=z)≈3 .+2grid.theta_nodes.+5z
        end
        reverse=solve_q4(grid;rpm=-2000.)
        for j in 0:4,i in 0:11
            @test r.pressure_pa[12j+i+1]≈reverse.pressure_pa[12j+mod(18-i,12)+1] rtol=1e-10 atol=1e-6
        end
        flooded=solve_q4(grid;end_pressure_pa=(1e8,1e8))
        signed=solve_q4(grid;closure=:signed,end_pressure_pa=(1e8,1e8))
        @test flooded.pressure_pa≈signed.pressure_pa rtol=1e-11
        # Differentiate the interpolated Couette flux rather than integrating N_theta*C.
        strong_b=zeros(60);ca=cosd(10.)
        for (cell,ids) in enumerate(grid.quads)
            i,j=mod1(cell,12),(cell-1)÷12+1
            dt,dz=a.axes.theta_widths[i],a.axes.z_widths[j]
            h=a.hnode[collect(ids)]
            for xi in (-1/sqrt(3),1/sqrt(3)),eta in (-1/sqrt(3),1/sqrt(3))
                N=[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)]/4
                htheta=((1-eta)*(h[2]-h[1])+(1+eta)*(h[3]-h[4]))/(2dt)
                z=a.axes.z_edges[j]+(1+eta)*dz/2
                rj=.05+(.05-z)*tand(10.)
                for q in 1:4
                    strong_b[ids[q]]-=dt*dz/4*N[q]*(2000*2pi/60*rj/(2ca))*htheta
                end
            end
        end
        @test norm(a.b-strong_b)/norm(strong_b)<1e-12
    end
    @testset "Q4 stationary cylinder and cone references" begin
        for angle in (0.,10.)
            p=cone(;angle,ecc=0.)
            errors=Float64[]
            for nz in (4,8,16)
                r=solve_q4(housing_grid(p,12,nz);rpm=0.,end_pressure_pa=(5e5,1e5))
                z=r.a.z
                r0=.05005+.05tand(angle);rL=.05005-.05tand(angle)
                exact=iszero(angle) ? 5e5 .-4e5.*z/.1 :
                    [5e5-4e5*log((.05005+(.05-x)*tand(angle))/r0)/log(rL/r0) for x in z]
                push!(errors,norm(r.pressure_pa-exact)/norm(exact))
            end
            @test iszero(angle) ? maximum(errors)<1e-11 : all(errors[2:end].<.3 .*errors[1:end-1])
            println("Q4 stationary angle=$angle errors=$errors");yield()
        end
    end
    @testset "Independent continuous cone forcing on the public mesh" begin
        p=cone(;ecc=0.);errors=Float64[]
        K=(50e-6*cosd(10.))^3/(12MU);ca=cosd(10.);L=.1
        for nt in (12,24,48,96)
            grid=housing_grid(p,nt,nt÷3);a=assemble_q4(grid;rpm=0.)
            source=zeros(length(a.b));points=(-1/sqrt(3),1/sqrt(3))
            for (cell,ids) in enumerate(grid.quads)
                i,j=mod1(cell,nt),(cell-1)÷nt+1
                dt,dz=a.axes.theta_widths[i],a.axes.z_widths[j]
                for xi in points,eta in points
                    theta=a.axes.theta_edges[i]+(1+xi)*dt/2
                    z=a.axes.z_edges[j]+(1+eta)*dz/2
                    r=.05005+(.05-z)*tand(10.)
                    exact=sin(theta)*sin(pi*z/L)
                    force=(K/(r*ca)+K*r*ca*(pi/L)^2)*exact+
                        K*tand(10.)*ca*pi/L*sin(theta)*cos(pi*z/L)
                    N=[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)]/4
                    for q in 1:4
                        source[ids[q]]+=dt*dz/4*N[q]*force
                    end
                end
            end
            free=findall(.!a.fixed);P=zeros(length(a.b));P[free]=a.A[free,free]\source[free]
            exact=sin.(a.theta).*sin.(pi.*a.z/L)
            push!(errors,norm(P-exact)/norm(exact))
            yield()
        end
        @test all(errors[2:end].<.3 .*errors[1:end-1])
        @test errors[end]<.001
        println("Q4 signed manufactured field errors=$errors")
    end
end
include(joinpath(@__DIR__,"../../../cli.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;output=false,description="Check the verification-only Q4 assembly on public unported meshes.") do cfg
        check_q4()
    end)
end
