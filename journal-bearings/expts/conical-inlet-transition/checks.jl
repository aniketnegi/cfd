using Test
include(joinpath(@__DIR__,"../conical-fitted-inlet/checks.jl"))
include(joinpath(@__DIR__,"layout.jl"))

function run_transition_checks()
    base=fitted_mesh();data=macro_data(base)
    x=vec(copy(data.nodes[data.free,:]));f,g=layout_energy(base,data,x)
    result=controlled_layout(base);controlled=result.mesh
    @testset "Coupled placement and unchanged physical geometry" begin
        indices=unique(round.(Int,range(1,length(x),length=12)))
        for i in indices
            dx=zeros(length(x));dx[i]=1e-8
            numerical=(layout_energy(base,data,x+dx)[1]-layout_energy(base,data,x-dx)[1])/2e-8
            @test abs(numerical-g[i])/max(1.,abs(g[i]))<1e-5
        end
        @test result.final_energy<result.initial_energy
        @test all(diff([h.energy for h in result.history]).<=0)
        @test controlled.quads==base.quads
        fixed=setdiff(1:length(base.nodes),data.free)
        @test controlled.nodes[fixed]==base.nodes[fixed]
        @test maximum(abs,controlled.areas-base.areas)>0
        @test abs(sum(controlled.areas)-sum(base.areas))<1e-12
        @test controlled.quality.maximum_neighbour_area_ratio<2
        @test controlled.quality.minimum_scaled_corner>.9
    end
    @testset "Shared-edge refinement, exact rim and worst-cell quality" begin
        for level in (1,2,4)
            m=subdivide_layout(controlled,level);old=subdivide_layout(base,level)
            d=macro_data(m)
            @test length(m.quads)==length(base.quads)*level^2
            @test m.quality.rim_edges==32level
            @test count(==(5),d.valence)==4
            @test all(v->v in (2,4,5),d.valence)
            @test all(v->d.valence[v]==2,d.boundary)
            @test m.quality.minimum_scaled_corner>.9
            @test m.quality.minimum_face_alignment>.9
            @test m.quality.maximum_neighbour_area_ratio<2
            @test abs(sum(m.areas)-sum(old.areas))<1e-12
            @test all(vcat(collect.(collect(d.rim))...)) do v
                x,y=m.nodes[v];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                abs(hypot(rh*cos(x/m.R),z-.05)-.0025)<1e-12
            end
            @test m.quality.minimum_scaled_corner>old.quality.minimum_scaled_corner
            @test m.quality.minimum_face_alignment>old.quality.minimum_face_alignment
            println("controlled level=$level $(m.quality)")
        end
    end
    @testset "Independent pressure on controlled macro maps" begin
        errors=Float64[]
        for level in (1,2,4)
            m=subdivide_layout(controlled,level);kx=2pi/m.period;ky=pi/(m.L/m.ca)
            exact(x,y)=1+.2*cos(kx*x)*sin(ky*y)
            a=fitted_operator(m,pose;unit_mobility=true,boundary_pressure=(x,y,t)->exact(x,y))
            P=[exact(m.centre[k,:]...) for k in eachindex(m.quads)]
            forcing=[begin
                x,y=m.centre[k,:];rh=m.R+(.05-y*m.ca)*tand(10.)
                .2*cos(kx*x)*((m.R/rh*kx^2+rh/m.R*ky^2)*sin(ky*y)+sind(10.)/m.R*ky*cos(ky*y))
            end for k in eachindex(m.quads)]
            computed=a.A\(a.b+forcing.*m.areas)
            error=sqrt(sum((computed-P).^2 .*m.areas)/sum((P.-1).^2 .*m.areas));push!(errors,error)
            @test minimum(computed)>0
            println("controlled pressure MMS level=$level error=$error")
        end
        @test errors[1]/errors[2]>2.5
        @test errors[2]/errors[3]>2.5
        @test errors[3]<.002
    end
    @testset "Curved cavity on controlled maps, both rotation signs" begin
        manufactured_pose=floating_bush_pose(reference,(-30e-6,0.,0.))
        for rpm in (2000.,-2000.)
            pe=Float64[];fe=Float64[]
            for level in (1,2)
                m=subdivide_layout(controlled,level)
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
                println("controlled JFO MMS rpm=$rpm level=$level pressure=$ep fill=$ef")
            end
            @test pe[2]<pe[1]/1.5
            @test fe[2]<fe[1]/1.5
            @test pe[2]<.01 && fe[2]<.005
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;run_transition_checks();end
