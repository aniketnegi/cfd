using Test
include(joinpath(@__DIR__,"solve.jl"))
BLAS.set_num_threads(1)
problem=fitted_load_problem(1);m=problem.mesh
nodes=copy(m.nodes);quads=copy(m.quads);centres=copy(m.centre)
@testset "Fitted load adapter: independent translated gaps and fixed housing" begin
    for d in ((0.,0.,0.),(12e-6,-8e-6,1e-4),(-5e-6,10e-6,-1e-4))
        pose=floating_bush_pose(problem.reference,d)
        a=fitted_operator(m,pose;rpm=0.,boundary_pressure=(x,y,t)->3e5)
        @test maximum(abs,a.A\a.b.-3e5)<1e-3
        for k in round.(Int,range(1,length(m.quads);length=12))
            theta=m.centre[k,1]/m.R;z=m.centre[k,2]*m.ca;bx,by,bz=d
            rj=.05+(.05-z-bz)*tand(10.);rh=.05005+(.05-z)*tand(10.)
            ray=-bx*cos(theta)-by*sin(theta)+sqrt(rj^2-(bx*sin(theta)-by*cos(theta))^2)
            @test a.h[k]≈(rh-ray)*cosd(10.) atol=1e-16
        end
        a=fitted_operator(m,pose)
        for f in round.(Int,range(1,length(m.owner);length=12))
            theta=m.mid[f][1]/m.R;z=m.mid[f][2]*m.ca;bx,by,bz=d
            rj=.05+(.05-z-bz)*tand(10.);rh=.05005+(.05-z)*tand(10.)
            ray=-bx*cos(theta)-by*sin(theta)+sqrt(rj^2-(bx*sin(theta)-by*cos(theta))^2)
            expected=2000*2pi/60*rj*(rh-ray)*cosd(10.)/2*m.normals[f][1]*m.lengths[f]
            @test a.adv[f]≈expected atol=1e-17 rtol=1e-10
        end
    end
    @test m.nodes==nodes && m.quads==quads && m.centre==centres
    @test_throws ArgumentError floating_bush_state(problem,(0.,0.,0.),128,64;feed_cells=falses(8192))
    @test_throws ArgumentError floating_bush_state(problem,(0.,0.,0.),256,128;feed_cells=nothing)
end
@testset "Known fitted state and transverse inverse recovery at shifted seating" begin
    d=(-9.152501266813789e-7,1.7844427760798323e-7,0.)
    state=floating_bush_state(problem,d,128,64;feed_cells=nothing)
    @test collect(state.bush_force_n)≈[-18.24857749,-1002.14633878,321.13450985] atol=1e-6
    @test state.film.net_feed*60000≈.0308873032652 rtol=1e-9
    @test state.film.imbalance<1e-9 && state.film.local_error<1e-9
    known=(-2e-6,.4e-6,1e-4)
    forward=floating_bush_state(problem,known,128,64;feed_cells=nothing)
    inverse=equilibrate_floating_bush(problem,128,64;feed_cells=nothing,fixed_axial=true,
        external_force_n=Tuple(-collect(forward.bush_force_n)),
        initial_displacement_m=(known[1]*1.03,known[2]*.98,known[3]),fd_step=1e-6)
    @test inverse.converged
    @test maximum(abs,collect(inverse.residual_n)[1:2])<.05
    @test norm(collect(inverse.state.pose.displacement_m).-collect(known))<1e-9
    @test inverse.state.pose.displacement_m[3]==known[3]
    @test m.nodes==nodes && m.quads==quads && m.centre==centres
end
@testset "Safeguarded scalar root and failure reporting" begin
    root=bisect_axial(x->10*(x-31.25),12.5,100.)
    @test root.converged && abs(root.c-31.25)<.002
    @test all(h->h.a<=h.c<=h.b,root.history)
    @test bisect_axial(x->x-12.5,12.5,100.).c==12.5
    @test_throws ArgumentError bisect_axial(x->x^2+1,12.5,100.)
    @test_throws ErrorException bisect_axial(x->NaN,12.5,100.)
    @test_throws ErrorException bisect_axial(x->error("radial solve failed"),12.5,100.)
    @test !bisect_axial(x->x-31.25,12.5,100.;max_iterations=1).converged
    @test_throws ArgumentError bisect_axial(identity,100.,12.5)
    @test_throws ArgumentError radial_equilibrium(problem,-1.,(0.,0.))
end
