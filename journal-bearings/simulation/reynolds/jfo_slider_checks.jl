using Test, LinearAlgebra
include(joinpath(@__DIR__,"jfo_slider.jl"))
include(joinpath(@__DIR__,"jfo_slider_reference.jl"))

# All benchmark inputs are from Lengiewicz et al.; these are NOT VG32 bearing inputs.
benchmark(n; u=2.0,mu=0.015,ends=(1e6,1e6),pc=0.0,kwargs...) =
    solve_jfo_slider(slider_gap,n;length_m=0.125,viscosity_pa_s=mu,
        entrainment_m_s=u,end_pressure_pa=ends,cavitation_pressure_pa=pc,kwargs...)

@testset "Continuum reference and printed author-supplement values" begin
    ref = slider_reference()
    # Published rounded outputs in mm, MPa; tolerate only their printed precision.
    @test ref.rupture_m ≈ 20.8491e-3 atol=5.1e-8 rtol=0
    @test ref.reformation_m ≈ 56.7045e-3 atol=5.1e-8 rtol=0
    @test ref.hstar_m ≈ 0.0175034e-3 atol=5.1e-11 rtol=0
    @test ref.pressure_at_zero_pa ≈ 3.31208e6 atol=5.1 rtol=0
    @test ref.pressure(-0.0625) ≈ 1e6 atol=1e-7
    @test ref.pressure(0.0625) ≈ 1e6 atol=1e-7
    @test abs(ref.gradient(ref.rupture_m)) < 1e-5
    @test ref.gradient(ref.reformation_m) > 0
    @test ref.fill(ref.reformation_m-1e-9) < 0.8
    @test ref.fill(ref.reformation_m+1e-9) == 1
    # Integration accuracy is controlled independently of the FVM grid.
    finer = slider_reference(quadrature_order=96)
    @test ref.rupture_m ≈ finer.rupture_m atol=1e-13
    @test ref.reformation_m ≈ finer.reformation_m atol=1e-13
    xs = range(-0.0625,0.0625;length=201)
    @test maximum(abs,ref.pressure.(xs)-finer.pressure.(xs)) < 1e-5
    @test_throws ArgumentError slider_reference(quadrature_order=2)
end

@testset "Published slider: whole fields, fronts and conserved flux" begin
    ref = slider_reference()
    errors,fill_l1 = Float64[],Float64[]
    for n in (40,80,160,320,640,1280)
        r = benchmark(n)
        exact = ref.pressure.(r.x_m)
        push!(errors,norm(r.pressure_pa-exact)/norm(exact))
        push!(fill_l1,sum(abs,r.fill_fraction-ref.fill.(r.x_m))/n)
        d = r.diagnostics
        @test minimum(r.pressure_pa) >= 0
        @test minimum(r.fill_fraction) >= 0
        @test maximum(r.fill_fraction) <= 1
        @test d.complementarity_residual < 1e-12
        @test d.relative_max_cell_imbalance < 1e-10
        @test d.relative_flux_spread < 1e-10
        @test d.relative_boundary_imbalance < 1e-10
        cavity = findall(r.fill_fraction .< 1-1e-8)
        @test !isempty(cavity)
        @test cavity == collect(first(cavity):last(cavity))
        @test 1 < first(cavity) < last(cavity) < n # both rupture and reformation
        # Cell classification locates a front only to its grid neighbourhood.
        rupture = r.face_x_m[first(cavity)]
        reformation = r.face_x_m[last(cavity)+1]
        @test abs(rupture-ref.rupture_m) <= 2*0.125/n
        @test abs(reformation-ref.reformation_m) <= 2*0.125/n
        @test abs(r.face_flux_m2_s[1]/ref.constant_flux_m2_s-1) < 1e-3
        @test abs(maximum(r.pressure_pa)/ref.peak_pressure_pa-1) < 1e-3
    end
    @test all(errors[2:end] .< 0.4 .* errors[1:end-1])
    @test errors[end] < 1e-4
    @test fill_l1[end] < 1e-3 # jump alignment makes this non-monotonic
    @test fill_l1[end] < 0.1fill_l1[1]
    println("Pressure relative L2 errors: $errors")
    println("Fill mean absolute errors (jump alignment matters): $fill_l1")
end

@testset "Pressure reference, reverse motion and scaling" begin
    r = benchmark(160)
    shifted = benchmark(160;ends=(1e6+101325,1e6+101325),pc=101325)
    @test shifted.pressure_above_cavitation_pa ≈ r.pressure_pa rtol=1e-12
    @test shifted.fill_fraction ≈ r.fill_fraction rtol=1e-12
    @test shifted.face_flux_m2_s ≈ r.face_flux_m2_s rtol=1e-12
    reverse_run = benchmark(160;u=-2.0)
    @test reverse_run.pressure_pa ≈ reverse(r.pressure_pa) rtol=1e-11
    @test reverse_run.fill_fraction ≈ reverse(r.fill_fraction) rtol=1e-11
    @test reverse_run.face_flux_m2_s ≈ -reverse(r.face_flux_m2_s) rtol=1e-11
    faster = benchmark(160;u=4.0,ends=(2e6,2e6))
    @test faster.pressure_pa ≈ 2r.pressure_pa rtol=1e-11
    @test faster.fill_fraction ≈ r.fill_fraction rtol=1e-11
    @test faster.face_flux_m2_s ≈ 2r.face_flux_m2_s rtol=1e-11
    thicker = benchmark(160;mu=0.030,ends=(2e6,2e6))
    @test thicker.pressure_pa ≈ 2r.pressure_pa rtol=1e-11
    @test thicker.fill_fraction ≈ r.fill_fraction rtol=1e-11
    @test thicker.face_flux_m2_s ≈ r.face_flux_m2_s rtol=1e-11
end

@testset "Full-film and explicit starved/dry inflow limits" begin
    flat = solve_jfo_slider(x->20e-6,40;length_m=0.125,viscosity_pa_s=0.015,
        entrainment_m_s=2.0,end_pressure_pa=(2e6,1e6),cavitation_pressure_pa=0.0)
    exact = 2e6 .- (flat.x_m .+ 0.0625)/0.125 .* 1e6
    @test flat.pressure_pa ≈ exact rtol=1e-11
    @test flat.fill_fraction == ones(40)
    expected_flux = 2*20e-6 + (20e-6)^3/(12*0.015)*1e6/0.125
    @test flat.face_flux_m2_s ≈ fill(expected_flux,41) rtol=1e-11
    diverging(x) = 20e-6+5e-6*(x+0.0625)/0.125
    for fraction in (0.0,0.4,1.0)
        r = solve_jfo_slider(diverging,80;length_m=0.125,viscosity_pa_s=0.015,
            entrainment_m_s=2.0,end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0,
            inflow_fill_fraction=fraction)
        @test maximum(abs,r.pressure_pa) < 1e-8
        @test r.face_flux_m2_s ≈ fill(2*20e-6*fraction,81) atol=1e-16
        # Upwind fraction at a cell is transported through its east face.
        expected_fill = 20e-6*fraction ./ diverging.(r.face_x_m[2:end])
        @test r.fill_fraction ≈ expected_fill atol=1e-12
        @test r.fill_fraction[end] < 1 # do not impose flooded outflow
    end
end

@testset "Reformation stress case and rejection of pressure clipping" begin
    # Section 4.2 geometry; conservation/two-cavity check only, not a full-field reference.
    wavy(x) = slider_gap(x)-2.5e-6*cos(20pi*x/0.125)
    r = solve_jfo_slider(wavy,320;length_m=0.125,viscosity_pa_s=0.015,
        entrainment_m_s=2.0,end_pressure_pa=(1e6,1e6),cavitation_pressure_pa=0.0)
    c = r.fill_fraction .< 1-1e-8
    @test count(i -> c[i] && (i==1 || !c[i-1]),eachindex(c)) == 2
    @test r.diagnostics.relative_flux_spread < 1e-10
    a = assemble_jfo_slider(slider_gap,160;length_m=0.125,viscosity_pa_s=0.015,
        entrainment_m_s=2.0,end_pressure_pa=(1e6,1e6),cavitation_pressure_pa=0.0)
    clipped = max.(a.A \ a.b,0.0)
    # The tempting pressure-only clamp does not satisfy the same face mass balance.
    q = a.g .* diff(-vcat(1e6,clipped,1e6)) .+ a.adv
    @test (maximum(q)-minimum(q))/maximum(abs,a.adv) > 0.05
end

@testset "Unsupported or inconsistent input fails explicitly" begin
    @test_throws ArgumentError benchmark(1)
    @test_throws ArgumentError benchmark(40;u=0.0)
    @test_throws ArgumentError benchmark(40;u=Inf)
    @test_throws ArgumentError benchmark(40;mu=0.0)
    @test_throws ArgumentError benchmark(40;mu=NaN)
    @test_throws ArgumentError benchmark(40;ends=(-1.0,1e6))
    @test_throws ArgumentError benchmark(40;pc=2e6)
    @test_throws ArgumentError benchmark(40;inflow_fill_fraction=0.5) # p_in>pcav
    @test_throws ArgumentError benchmark(40;inflow_fill_fraction=-0.1)
    @test_throws ArgumentError benchmark(40;inflow_fill_fraction=1.1)
    @test_throws ArgumentError benchmark(40;max_iterations=0)
    @test_throws ArgumentError benchmark(40;tolerance=NaN)
    @test_throws ErrorException benchmark(40;max_iterations=1)
    for gap in (x->0.0,x->-1e-5,x->NaN,x->Inf)
        @test_throws ArgumentError solve_jfo_slider(gap,40;length_m=0.125,
            viscosity_pa_s=0.015,entrainment_m_s=2.0,
            end_pressure_pa=(1e6,1e6),cavitation_pressure_pa=0.0)
    end
end
