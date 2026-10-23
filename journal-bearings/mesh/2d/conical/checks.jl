using Test
include(joinpath(@__DIR__, "surface_feed.jl"))

p = BearingParameters(journal_midpoint_radius=50.0, radial_clearance=0.05,
                      bearing_length=100.0, semicone_angle_deg=10.0,
                      eccentricity_ratio=0.6)

@testset "Shared conical housing grid" begin
    plain = housing_grid(p, 512, 256)
    fed = surface_feed_grid(p, 512, 256; feed_diameter=5.0, feed_z=50.0)
    @test fed.nodes == plain.nodes
    @test fed.quads == plain.quads
    @test fed.radial_gaps == plain.radial_gaps
    @test fed.areas == plain.areas
    @test length(fed.nodes) == 512 * 257
    @test length(fed.quads) == length(fed.feed_cells) == 512 * 256
    @test fed.quads[512] == (512, 1, 513, 1024) # periodic seam, no duplicate nodes
    @test minimum(fed.radial_gaps) > 0
    @test sum(fed.areas) ≈ 2pi * p.housing_midpoint_radius * p.bearing_length /
                            cosd(p.semicone_angle_deg) rtol=1e-12

    # Independent projected mouth area: at each z the admissible angular band
    # on +Y has width 2asin(sqrt(r²-(z-z0)²)/Rh(z)). Integrate the cone metric.
    r, z0, n = 2.5, 50.0, 10_000
    dz = 2r / n
    exact_area = sum(1:n) do j
        z = z0 - r + (j - 0.5) * dz
        rh = p.housing_midpoint_radius + (p.half_bearing_length - z) * p.slope
        rh * 2asin(sqrt(r^2 - (z - z0)^2) / rh) * dz /
            cosd(p.semicone_angle_deg)
    end
    # Points on the drilled rim lie on the conical housing and have the
    # prescribed XZ-plane radius even though the grid does not follow the rim.
    for angle in range(0, 2pi; length=17)
        x = r * cos(angle)
        z = z0 + r * sin(angle)
        rh = p.housing_midpoint_radius + (p.half_bearing_length - z) * p.slope
        theta = atan(sqrt(rh^2 - x^2), x)
        xh, yh, _ = bearing_point(p, theta, z, 1.0)
        @test yh > 0
        @test xh^2 + (z - z0)^2 ≈ r^2 atol=1e-12
    end
    masked_area = sum(fed.areas[fed.feed_cells])
    @test 0 < count(fed.feed_cells) < length(fed.quads)
    # The 512×256 grid still stair-steps the exact circle; this is an area
    # diagnostic, not a claim that the mask has a conforming 5 mm boundary.
    @test masked_area ≈ exact_area rtol=0.03

    dtheta = 2pi / fed.n_theta
    axial_step = p.bearing_length / fed.n_axial
    opposite_intersection = false
    for j in 1:fed.n_axial, i in 1:fed.n_theta
        theta, z = (i - 0.5) * dtheta, (j - 0.5) * axial_step
        x, y, _ = bearing_point(p, theta, z, 1.0)
        inside_cylinder = x^2 + (z - z0)^2 <= r^2
        k = i + fed.n_theta * (j - 1)
        @test fed.feed_cells[k] == (y > 0 && inside_cylinder)
        opposite_intersection |= y < 0 && inside_cylinder && !fed.feed_cells[k]
    end
    @test opposite_intersection

    mktempdir() do dir
        path = joinpath(dir, "feed.vtk")
        write_vtk(path, fed; feed_cells=fed.feed_cells)
        output = read(path, String)
        @test occursin("POINTS $(length(fed.nodes)) double", output)
        @test occursin("POLYGONS $(length(fed.quads)) $(5 * length(fed.quads))", output)
        @test occursin("CELL_DATA $(length(fed.quads))", output)
        mask = split(output, "SCALARS feed_cell int 1\nLOOKUP_TABLE default\n")
        @test length(mask) == 2
        @test count(==("1"), split(chomp(mask[2]), '\n')) == count(fed.feed_cells)
        @test_throws ArgumentError write_vtk(joinpath(dir, "bad.vtk"), fed;
                                              feed_cells=[true])
        write_vtk(joinpath(dir, "unported.vtk"), plain)
        @test !occursin("SCALARS feed_cell", read(joinpath(dir, "unported.vtk"), String))
    end
end

@testset "Feed geometry limits" begin
    make_grid(; diameter=5.0, z=50.0, ntheta=256, nz=64) =
        surface_feed_grid(p, ntheta, nz; feed_diameter=diameter, feed_z=z)
    for diameter in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError make_grid(; diameter)
    end
    for z in (NaN, -1.0, 2.5, 97.5, 102.5)
        @test_throws ArgumentError make_grid(; z)
    end
    @test_throws ArgumentError make_grid(; diameter=101.0)
    @test_throws ArgumentError make_grid(; ntheta=8, nz=4) # unresolved at cell centres
    @test_throws ArgumentError make_grid(; ntheta=2)
    @test_throws ArgumentError make_grid(; nz=0)

    cylinder = BearingParameters(journal_midpoint_radius=50.0, radial_clearance=0.05,
                                 bearing_length=100.0, semicone_angle_deg=0.0,
                                 eccentricity_ratio=0.6)
    sample = surface_feed_grid(cylinder, 512, 256; feed_diameter=5.0, feed_z=50.0)
    @test sum(sample.areas[sample.feed_cells]) ≈ pi * 2.5^2 rtol=0.03
end
