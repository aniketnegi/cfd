using Test, LinearAlgebra
include(joinpath(@__DIR__,"../reynolds/src/cylindrical.jl"))

function cylinder(; length_mm=100.0, ecc=0.6, angle=0.0)
    BearingParameters(journal_midpoint_radius=50.0, bearing_length=length_mm,
                      radial_clearance=0.050, semicone_angle_deg=angle,
                      eccentricity_ratio=ecc)
end
const MU = 32e-6 * 860.0  # nominal ISO VG 32 at 40 °C, illustrative density
const RPM = 2000.0
solve(p, nt=32, nz=16; rpm=RPM, mu=MU) =
    solve_cylindrical_full_film(p, nt, nz; rpm, viscosity_pa_s=mu)

@testset "Cylindrical signed full-film conservation and symmetries" begin
    p = cylinder()
    r = solve(p)
    @test minimum(r.pressure_pa) < 0 < maximum(r.pressure_pa)
    @test r.diagnostics.relative_l1_residual < 1e-11
    @test r.diagnostics.relative_max_cell_residual < 1e-10
    @test r.diagnostics.relative_global_imbalance < 1e-11
    @test r.diagnostics.absolute_side_flow_m3_s > 0
    # Equal ambient axial ends and no z-dependence of gap or speed.
    pressure = reshape(r.pressure_pa, 32, 16)
    @test pressure ≈ reverse(pressure; dims=2) rtol=1e-10
    @test solve(p; rpm=-RPM).pressure_pa ≈ -r.pressure_pa rtol=1e-10
    @test solve(p; mu=2MU).pressure_pa ≈ 2r.pressure_pa rtol=1e-10
    @test maximum(abs, solve(cylinder(ecc=0.0)).pressure_pa) < 1e-10
    @test_throws ArgumentError solve(cylinder(angle=10.0))
    @test_throws ArgumentError solve(p; mu=0)
    @test_throws ArgumentError solve(p; rpm=Inf)
    @test_throws ArgumentError solve_cylindrical_full_film(p, 3, 2; rpm=RPM, viscosity_pa_s=MU)
    mktempdir() do dir
        path = joinpath(dir, "pressure.vtk")
        write_pressure_vtk(path, r)
        output = read(path, String)
        @test occursin("SCALARS pressure_gauge_pa double 1\nLOOKUP_TABLE default", output)
        @test length(split(output, "SCALARS pressure_gauge_pa double 1\nLOOKUP_TABLE default\n")) == 2
    end
end

@testset "Finite L/D=1 grid convergence" begin
    p = cylinder()
    coarse, middle, fine = [solve(p, nt, nz) for (nt, nz) in
                            ((32, 16), (64, 32), (128, 64))]
    # Compare coarse cell pressures with area averages of four fine cells.
    function mismatch(a, b, nt, nz)
        ac = reshape(a.pressure_pa, nt, nz)
        bf = reshape(b.pressure_pa, 2nt, 2nz)
        avg = [(bf[2i-1, 2j-1] + bf[2i, 2j-1] +
                bf[2i-1, 2j] + bf[2i, 2j]) / 4 for i in 1:nt, j in 1:nz]
        norm(ac - avg) / norm(avg)
    end
    first = mismatch(coarse, middle, 32, 16)
    second = mismatch(middle, fine, 64, 32)
    @test 0 < second < 0.4first
    @test second < 0.005
end

@testset "Short-bearing analytical limit (San Andrés Eq. 4.10)" begin
    # L/D = 0.01: axial Poiseuille flow dominates, signed full-film
    # p ≈ 3μ(U/R)(dh/dθ)(z² - L²/4)/h³, with h=c-e*cosθ.
    p = cylinder(length_mm=1.0)
    errors = Float64[]
    for (nt, nz) in ((24, 12), (48, 24), (96, 48))
        r = solve(p, nt, nz)
        R, L = p.housing_midpoint_radius * 1e-3, p.bearing_length * 1e-3
        c, e = p.radial_clearance * 1e-3, p.eccentricity * 1e-3
        U = (RPM * 2pi / 60) * p.journal_midpoint_radius * 1e-3
        expected = [begin
            theta = (i - 0.5) * 2pi / nt
            z = (j - 0.5) * L / nz - L / 2
            h = c - e * cos(theta)
            3MU * U / R * e * sin(theta) * (z^2 - L^2/4) / h^3
        end for i in 1:nt, j in 1:nz]
        push!(errors, norm(r.pressure_pa - vec(expected)) / norm(expected))
        @test r.diagnostics.relative_global_imbalance < 1e-10
    end
    @test errors[1] > errors[2] > errors[3]
    @test errors[3] < 0.001  # limit also has O((L/R)^2) and radial-gap errors
end
