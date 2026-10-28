using Test
include(joinpath(@__DIR__,"../../simulation/reynolds/jfo_slider.jl"))
include(joinpath(@__DIR__,"reynolds_cavitation.jl"))

@testset "Independent Reynolds rupture reference (zero-pressure inlet/outlet slider)" begin
    L,mu,u = 0.125,0.015,2.0
    gap(x) = 20e-6-5e-6*cos(2pi*x/L)
    # Integrate the continuum ODE directly, independent of finite-volume A.
    order = 96
    eigenpair = eigen(SymTridiagonal(zeros(order),[k/sqrt(4k^2-1) for k in 1:order-1]))
    nodes,weights = eigenpair.values,2 .* eigenpair.vectors[1,:].^2
    integrate(f,a,b) = (b-a)/2*sum(weights .* f.((a+b)/2 .+ (b-a)/2 .* nodes))
    residual(xr) = integrate(x->(gap(x)-gap(xr))/gap(x)^3,-L/2,xr)
    lo,hi = 0.0,L/2
    @test residual(lo)>0>residual(hi)
    for _ in 1:60
        mid = (lo+hi)/2
        if residual(mid)>0
            lo=mid
        else
            hi=mid
        end
    end
    rupture = (lo+hi)/2
    gradient(x) = 12mu*u*(gap(x)-gap(rupture))/gap(x)^3
    pressure(x) = x<rupture ? integrate(gradient,-L/2,x) : 0.0
    @test abs(gradient(rupture)) < 1e-8
    @test abs(integrate(gradient,-L/2,rupture)) < 1e-5
    errors,front_errors,slopes = Float64[],Float64[],Float64[]
    for n in (80,160,320,640)
        a = assemble_jfo_slider(gap,n;length_m=L,viscosity_pa_s=mu,entrainment_m_s=u,
            end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0)
        s = pressure_obstacle(a.A,a.b,a.pressure_scale_pa,maximum(abs,a.adv))
        exact = pressure.(a.x_m)
        push!(errors,norm(s.pressure-exact)/norm(exact))
        last_full = findlast(s.pressure .> 1e-8)
        push!(front_errors,abs(a.face_x_m[last_full+1]-rupture))
        push!(slopes,abs(s.pressure[last_full]/(L/n))/(maximum(exact)/L))
        @test minimum(s.pressure) >= 0
        @test minimum(s.reaction) >= -1e-14
        @test front_errors[end] < 2L/n
        @test maximum(abs,(s.pressure/a.pressure_scale_pa).*(s.reaction/maximum(abs,a.adv))) < 1e-10
    end
    @test all(errors[2:end] .< 0.4 .* errors[1:end-1])
    @test errors[end] < 1e-4
    @test slopes[end] < 0.05 # one-sided rupture gradient tends to zero with refinement
    println("Reynolds slider rupture=$(rupture*1000) mm; pressure errors=$errors; front errors=$front_errors; normalized rupture slopes=$slopes")
end

@testset "Pressure-only cone: full film, KKT, and not a one-pass clamp" begin
    p = BearingParameters(journal_midpoint_radius=50.0,bearing_length=100.0,
        radial_clearance=0.05,semicone_angle_deg=10.0,eccentricity_ratio=0.6)
    flooded = fig6_reynolds(p,32,16;end_pressure_pa=(3e7,3e7))
    old = solve_full_film(housing_grid(p,32,16);rpm=2000.0,viscosity_pa_s=0.0277,end_pressure_pa=(3e7,3e7))
    @test flooded.pressure_pa ≈ old.pressure_pa rtol=1e-11
    for ends in ((0.0,0.0),(5e5,0.0),(0.0,5e5),(5e5,5e5))
        r = fig6_reynolds(p,64,32;end_pressure_pa=ends)
        @test minimum(r.pressure_pa) >= 0
        @test r.positive_pressure_relative_residual < 1e-10
        @test r.complementarity_relative < 1e-10
        @test r.cavity_full_film_flux_defect_m3_s > 0
        @test r.extrapolated_end_outflow_m3_s ≈ r.cavity_full_film_flux_defect_m3_s rtol=1e-10
        s = solve_full_film(housing_grid(p,64,32);rpm=2000.0,viscosity_pa_s=0.0277,end_pressure_pa=ends)
        @test norm(r.pressure_pa-max.(s.pressure_pa,0.0))/norm(r.pressure_pa) > 0.01
    end
    r = fig6_reynolds(p,32,16;end_pressure_pa=(5e5,5e5))
    shifted = fig6_reynolds(p,32,16;end_pressure_pa=(601325.0,601325.0),cavitation_pressure_pa=101325.0)
    @test shifted.pressure_above_cavitation_pa ≈ r.pressure_pa rtol=1e-12
    reversed = fig6_reynolds(p,32,16;rpm=-2000.0,end_pressure_pa=(5e5,5e5))
    @test reshape(reversed.pressure_pa,32,16) ≈ reverse(reshape(r.pressure_pa,32,16);dims=1) rtol=1e-11
    for column in (1,17)
        fixed = fig6_reynolds(p,32,16;pinned_column=column,end_pressure_pa=(5e5,5e5))
        @test all(fixed.pressure_pa[column:32:end] .== 0)
        @test fixed.positive_pressure_relative_residual < 1e-10
        @test fixed.extrapolated_end_outflow_m3_s ≈
              fixed.cavity_full_film_flux_defect_m3_s+fixed.pinned_line_exchange_m3_s rtol=1e-10
    end
    @test_throws ArgumentError fig6_reynolds(p,32,16;pinned_column=0)
    @test_throws ArgumentError fig6_reynolds(p,32,16;end_pressure_pa=(-1.0,0.0))
    @test_throws ErrorException fig6_reynolds(p,32,16;max_iterations=1)
end
