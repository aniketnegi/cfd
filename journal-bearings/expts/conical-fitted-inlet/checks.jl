using Test
include(joinpath(@__DIR__,"operator.jl"))
reference=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
    bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
pose=floating_bush_pose(reference,(-9.152501266813789e-7,1.7844427760798323e-7,0.))

# Continuum derivatives, evaluated independently of the discrete matrices.
function fitted_mms_point(x,y,m;rpm=2000.)
    theta=x/m.R;z=y*m.ca;L=.1;c=50e-6;e=30e-6;mu=.0277;ca=m.ca;k=tand(10.)
    st,ct=sincos(theta);rj=.05+(.05-z)*k;r=rj+c;D=sqrt(rj^2-e^2*st^2)
    h=ca*(c-e*ct+e^2*st^2/(rj+D));ht=ca*(e*st+e^2*st*ct/D);hz=ca*k*(rj/D-1)
    K=h^3/(12mu);Kt=h^2*ht/(4mu);Kz=h^2*hz/(4mu)
    t=ct+.25*cos(2pi*z/L);tt=-st;ttt=-ct
    tz=-.25*(2pi/L)*sin(2pi*z/L);tzz=-.25*(2pi/L)^2*cos(2pi*z/L)
    f=z*(L-z)/L^2;df=(L-2z)/L^2;ddf=-2/L^2;P0=1e7
    P=0.;F=1.;pt=0.;ptt=0.;pz=0.;pzz=0.;Ft=0.
    if t>0
        P=P0*f*t^4;pt=P0*f*4t^3*tt;ptt=P0*f*(12t^2*tt^2+4t^3*ttt)
        pz=P0*(df*t^4+f*4t^3*tz)
        pzz=P0*(ddf*t^4+8df*t^3*tz+f*(12t^2*tz^2+4t^3*tzz))
    else
        F=1-.2*(-t)^4;Ft=.8*(-t)^3*tt
    end
    diffusion=(Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
    transport=rpm*2pi/60*rj*(ht*F+h*Ft)/(2ca)
    return (;P,F,source=(transport-diffusion)*ca/m.R)
end

function run_fitted_checks()
    @testset "Fitted surface topology, geometry and forces" begin
        for level in (1,2,4)
            m=fitted_mesh(level)
            @test m.quality.rim_edges==32level
            @test count(==(:large),m.boundary)==128level
            @test count(==(:small),m.boundary)==128level
            @test m.quality.minimum_scaled_corner>0
            @test m.quality.minimum_face_alignment>0
            rimids=unique(vcat([collect(m.quads[m.owner[f]]) for f in eachindex(m.owner) if m.boundary[f]==:feed]...))
            # Test actual boundary edge endpoints, not the other cell corners.
            edges=Dict{Tuple{Int,Int},Int}()
            for q in m.quads,i in 1:4
                edge=Tuple(sort([q[i],q[mod1(i+1,4)]]));edges[edge]=get(edges,edge,0)+1
            end
            boundarynodes=unique(vcat([collect(e) for (e,count) in edges if count==1 &&
                all(id->1e-10<m.nodes[id][2]<m.L/m.ca-1e-10,e)]...))
            @test length(boundarynodes)==32level
            @test all(boundarynodes) do id
                x,y=m.nodes[id];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                abs((rh*cos(x/m.R))^2+(z-.05)^2-.0025^2)<1e-17
            end
            @test all(v->v in (1,2),values(edges))
            uniform=fitted_force(m,fill(5e5,length(m.quads))).force
            @test abs(uniform[1])<.002
            @test abs(uniform[2]/(-5e5*pi*.0025^2)-1)<.008/level^2
            expected_z=5e5*2pi*m.R*.1*tand(10.)-5e5*pi*.0025^2*tand(10.)
            @test abs(uniform[3]-expected_z)<.02/level^2+.003
            println("mesh level=$level $(m.quality)")
        end
    end
    @testset "Regular-grid reduction and exact linear/constant fields" begin
        m=fitted_mesh(1;regular=true)
        a=fitted_operator(m,pose)
        old=assemble_conical_full_film(pose.parameters,128,64;rpm=2000.,viscosity_pa_s=.0277,
            eccentricity_angle_rad=pose.eccentricity_angle_rad)
        @test norm(a.A-old.A)/norm(old.A)<1e-9
        @test norm(a.b-old.b)/norm(old.b)<1e-7
        m=fitted_mesh()
        a=fitted_operator(m,pose;rpm=0.,boundary_pressure=(x,y,t)->3e5)
        @test maximum(abs,a.A\a.b.-3e5)/3e5<1e-10
        a=fitted_operator(m,pose;unit_mobility=true,boundary_pressure=(x,y,t)->y)
        exact=m.centre[:,2]
        residual=a.A*exact-a.b-sind(10.)/m.R*m.areas
        @test maximum(abs,residual)<1e-11
        @test maximum(abs,a.Gx*exact+a.g0[:,1])<1e-10
        @test maximum(abs,a.Gy*exact+a.g0[:,2].-1)<1e-10
    end
    @testset "Independent smooth pressure manufactured solution" begin
        errors=Float64[]
        for level in (1,2,4)
            m=fitted_mesh(level);kx=2pi/m.period;ky=pi/(m.L/m.ca)
            exact(x,y)=1+.2*cos(kx*x)*sin(ky*y)
            a=fitted_operator(m,pose;unit_mobility=true,boundary_pressure=(x,y,t)->exact(x,y))
            P=[exact(m.centre[k,:]...) for k in eachindex(m.quads)]
            source=[begin
                x,y=m.centre[k,:];rh=m.R+(.05-y*m.ca)*tand(10.)
                .2*cos(kx*x)*((m.R/rh*kx^2+rh/m.R*ky^2)*sin(ky*y)+sind(10.)/m.R*ky*cos(ky*y))
            end for k in eachindex(m.quads)]
            result=a.A\(a.b+source.*m.areas)
            err=sqrt(sum((result-P).^2 .*m.areas)/sum((P.-1).^2 .*m.areas));push!(errors,err)
            @test minimum(result)>0
            println("pressure MMS level=$level error=$err")
        end
        @test errors[1]/errors[2]>2.5
        @test errors[2]/errors[3]>2.5
        @test errors[3]<.002
    end
    @testset "Independent curved pressure/content manufactured fields" begin
        manufactured_pose=floating_bush_pose(reference,(-30e-6,0.,0.))
        for rpm in (2000.,-2000.)
            pe=Float64[];fe=Float64[]
            for level in (1,2)
                m=fitted_mesh(level)
                a=fitted_operator(m,manufactured_pose;rpm,
                    boundary_pressure=(x,y,t)->fitted_mms_point(x,y,m;rpm).P,
                    incoming_fill=(x,y,t)->fitted_mms_point(x,y,m;rpm).F)
                exact=[fitted_mms_point(m.centre[k,:]...,m;rpm) for k in eachindex(m.quads)]
                source=[f.source for f in exact].*m.areas
                s=solve_jfo_active_set(a.A,a.B,a.b+source,1e8;max_iterations=250)
                P=[f.P for f in exact];F=[f.F for f in exact]
                ep=sqrt(sum((s.pressure-P).^2 .*m.areas)/sum(P.^2 .*m.areas))
                ef=sum(abs.(1 .-s.void-F).*m.areas)/sum(m.areas)
                push!(pe,ep);push!(fe,ef)
                @test minimum(s.pressure)>=-1e-5
                @test minimum(1 .-s.void)>=0 && maximum(1 .-s.void)<=1+1e-10
                flux=a.pressure_flux*s.pressure+a.constant+a.adv-a.transport*s.void
                @test maximum(abs,a.divergence*flux-source)/maximum(abs,a.adv)<1e-8
                println("JFO MMS rpm=$rpm level=$level pressure=$ep fill=$ef iterations=$(s.iterations)")
            end
            @test pe[2]<pe[1]/1.5
            @test fe[2]<fe[1]/1.5
            @test pe[2]<.01 && fe[2]<.005
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;run_fitted_checks();end
