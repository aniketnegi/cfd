using Test
include(joinpath(@__DIR__,"../verification/jfo/cylindrical/cases.jl"))

@testset "Supplied tensor JFO: compatibility and independent face balances" begin
    p=cylinder_parameters(1.,.6;model="production")
    grid=housing_grid(p,32,16)
    kwargs=(;rpm=RPM,viscosity_pa_s=MU,end_pressure_pa=(0.,0.),cavitation_pressure_pa=-PRESSURE_SCALE)
    counted=solve_conical_jfo(p,32,16;kwargs...)
    supplied=solve_conical_jfo(grid;kwargs...)
    @test counted.pressure_pa==supplied.pressure_pa
    @test counted.fill_fraction==supplied.fill_fraction
    @test counted.diagnostics==supplied.diagnostics
    for graded in (false,true),rpm in (-RPM,RPM)
        g=cylinder_grid(cylinder_parameters(1.,.6),16,8;graded)
        a=assemble_conical_jfo(g;kwargs...,rpm)
        n=length(a.b);P=[2e5*(1+sin(k)) for k in 1:n];F=[.7+.2cos(k) for k in 1:n]
        flux=conical_jfo_fluxes(a,P,F)
        @test flux.balance≈a.A*P+a.B*(1 .-F)-a.b rtol=1e-12
        @test maximum(abs,vec(sum(a.B;dims=1)))==0
        @test sum(flux.balance)≈sum(flux.qz[:,end])-sum(flux.qz[:,1]) atol=1e-18
        axes=a.axes
        for i in (1,5,16),j in (1,8)
            next=mod1(i+1,16);up=rpm>0 ? i : next
            t=axes.theta_edges[i+1];h=CLEARANCE*(1-.6cos(t));width=axes.z_widths[j]
            expected=rpm*2pi/60*R*h*width/2*F[up+16*(j-1)]+
                h^3*width/(12MU*R*((axes.theta_widths[i]+axes.theta_widths[next])/2))*(P[i+16*(j-1)]-P[next+16*(j-1)])
            @test flux.qt[i,j]≈expected rtol=1e-13
        end
    end
end

@testset "Cylindrical JFO analytical limits and submerged exchange" begin
    for graded in (false,true)
        g=cylinder_grid(cylinder_parameters(1.,0.),16,8;graded)
        r=solve_conical_jfo(g;rpm=0.,viscosity_pa_s=MU,end_pressure_pa=(2e5,1e5),cavitation_pressure_pa=0.)
        a=full_film_axes(g)
        exact=[2e5-1e5*z/.1 for z in a.z_centres for t in a.theta_centres]
        @test r.pressure_pa≈exact rtol=1e-12
        @test r.fill_fraction==ones(length(exact))
        Q=2pi*R*CLEARANCE^3/(12MU)*1e5/.1
        @test -sum(r.axial_face_flow_m3_s[:,1])≈-Q rtol=1e-12
        @test sum(r.axial_face_flow_m3_s[:,end])≈Q rtol=1e-12
    end
    for pv in (-1.,-.1)
        r=cylinder_solution(1.,.6,pv,64)
        s=observables(r)
        @test minimum(r.pressure_pa)>=pv*PRESSURE_SCALE
        @test 0<=minimum(r.fill_fraction)<1
        @test s.cavitated_cells>0
        @test s.gross_in_m3_s>0 && s.gross_out_m3_s>0
        @test s.gross_in_m3_s≈s.gross_out_m3_s rtol=1e-10
        @test s.local_balance<1e-10 && s.global_balance<1e-10 && s.complementarity<1e-10
    end
    g=cylinder_grid(cylinder_parameters(.25,.2),64,8)
    r=solve_conical_jfo(g;rpm=RPM,viscosity_pa_s=MU,end_pressure_pa=(0.,0.),cavitation_pressure_pa=-PRESSURE_SCALE)
    signed=solve_full_film(g;rpm=RPM,viscosity_pa_s=MU,end_pressure_pa=(0.,0.))
    @test r.fill_fraction==ones(64*8)
    @test r.pressure_pa≈signed.pressure_pa rtol=1e-11 atol=1e-8
end

@testset "Reference shifts, scaling, reversal and sampling" begin
    g=cylinder_grid(cylinder_parameters(1.,.6),32,16)
    kw=(;rpm=RPM,viscosity_pa_s=MU,end_pressure_pa=(0.,0.),cavitation_pressure_pa=-PRESSURE_SCALE)
    r=solve_conical_jfo(g;kw...)
    shifted=solve_conical_jfo(g;kw...,end_pressure_pa=(101325.,101325.),cavitation_pressure_pa=101325-PRESSURE_SCALE)
    @test shifted.pressure_pa≈r.pressure_pa.+101325 rtol=1e-11
    @test shifted.fill_fraction≈r.fill_fraction rtol=1e-11
    reverse=solve_conical_jfo(g;kw...,rpm=-RPM)
    @test reshape(reverse.pressure_pa,32,16)≈Base.reverse(reshape(r.pressure_pa,32,16);dims=1) rtol=1e-10 atol=1e-7
    @test reshape(reverse.fill_fraction,32,16)≈Base.reverse(reshape(r.fill_fraction,32,16);dims=1) rtol=1e-10
    viscous=solve_conical_jfo(g;kw...,viscosity_pa_s=2MU,cavitation_pressure_pa=-2PRESSURE_SCALE)
    @test viscous.pressure_pa≈2r.pressure_pa rtol=1e-11
    @test viscous.fill_fraction≈r.fill_fraction rtol=1e-11
    @test viscous.axial_face_flow_m3_s≈r.axial_face_flow_m3_s rtol=1e-10
    @test periodic_profile([0.,pi],[2.,4.],[0.,pi,2pi,pi/2])≈[2.,4.,2.,3.]
    axes=full_film_axes(g);field=[3+5z for z in axes.z_centres for t in axes.theta_centres]
    @test plane_profile(g,field)≈fill(3.25,32)
    matrix=reshape(field,32,16)
    for t in (-4.,0.,8.),z in (.02,.05,.08)
        @test field_at(axes,matrix,t,z)≈3+5z
    end
    @test_throws ArgumentError solve_conical_jfo(g;kw...,cavitation_pressure_pa=0.)
    @test_throws ArgumentError assemble_conical_jfo((;g...,feed_edges=[]);kw...)
    @test_throws ArgumentError assemble_conical_jfo(g;kw...,rpm=NaN)
    @test_throws ArgumentError assemble_conical_jfo(g;kw...,viscosity_pa_s=0.)
    @test_throws ArgumentError axial_count(1.,65)
end
