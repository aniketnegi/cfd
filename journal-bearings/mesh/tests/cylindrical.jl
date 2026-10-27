# Explicit cylindrical coverage for the shared surface mesher and CLI.
using Test
include(joinpath(@__DIR__, "../2d/generate.jl"))

function check_cylinder(grid, p, theta_centres)
    radius = p.journal_midpoint_radius + p.radial_clearance
    @test grid.acceptance.accepted
    @test grid.acceptance.valid_geometry
    @test grid.acceptance.quality_pass
    @test grid.acceptance.placement_stable
    @test all(node -> isapprox(hypot(node[1], node[2]), radius; atol=1e-11, rtol=0), grid.nodes)
    @test extrema(node[3] for node in grid.nodes) == (0., p.bearing_length)
    @test all(area -> isfinite(area) && area > 0, grid.areas)
    # Independent ray/circle intersection, not the shared radius/gap helpers.
    @test all(eachindex(grid.quads)) do cell
        theta = theta_centres[cell]
        ray = p.eccentricity*cos(theta) +
            sqrt(p.journal_midpoint_radius^2 - (p.eccentricity*sin(theta))^2)
        isapprox(grid.radial_gaps[cell], radius-ray; atol=1e-11, rtol=0)
    end
end

@testset "Cylindrical unported and fitted surfaces" begin
    for (name, dimensions, counts, diameter, feed_z) in (
        ("reference", (journal_midpoint_radius=50., bearing_length=100., radial_clearance=.05), (128,64), 5., 50.),
        ("varied", (journal_midpoint_radius=40., bearing_length=80., radial_clearance=.08), (160,52), 4., 35.),
    )
        @testset "$name" begin
            p = BearingParameters(; dimensions..., semicone_angle_deg=0., eccentricity_ratio=.6)
            nt, nz = counts
            plain = housing_grid(p, nt, nz)
            theta = [(i-.5)*2pi/nt for j in 1:nz for i in 1:nt]
            check_cylinder(plain, p, theta)
            @test length(plain.quads) == nt*nz
            @test length(plain.nodes) == nt*(nz+1)
            # Accumulate in high precision so this tests cell areas, not the
            # rounding error of a long Float64 reduction.
            @test sum(BigFloat, plain.areas) ≈ 2pi*(p.journal_midpoint_radius+p.radial_clearance)*p.bearing_length rtol=1e-14
            @test plain.acceptance.termination_reason == :direct_construction
            @test plain.quads[nt] == (nt, 1, nt+1, 2nt)

            fitted = surface_feed_grid(p, nt, nz; feed_diameter=diameter, feed_z)
            surface = fitted.surface
            check_cylinder(fitted, p, surface.centre[:,1]./surface.R)
            @test surface.ca == 1.
            @test surface.slope == 0.
            @test fitted.areas ≈ surface.areas.*1e6 rtol=1e-14
            @test length(fitted.feed_edges) == count(==(:feed), surface.boundary) > 0
            @test all(fitted.feed_edges) do edge
                all(edge) do node
                    x, y, z = fitted.nodes[node]
                    y > 0 && isapprox(x^2+(z-feed_z)^2, (diameter/2)^2; atol=1e-10, rtol=0)
                end
            end
            # All exterior faces are either axial exits or the fitted rim;
            # the angular seam must remain an internal periodic connection.
            @test all(eachindex(surface.owner)) do face
                (surface.neighbour[face] == 0) == (surface.boundary[face] in (:large, :small, :feed))
            end
            # Independent circular-drill footprint integral on a cylinder.
            # z=feed_z+a*sin(t) removes the square-root endpoint singularity.
            radius = p.journal_midpoint_radius+p.radial_clearance
            a = diameter/2
            steps = 4096
            dt = pi/steps
            integrand(t) = 2radius*asin(a*cos(t)/radius)*a*cos(t)
            opening = dt/3 * sum(i -> (i in (0,steps) ? 1 : isodd(i) ? 4 : 2)*
                integrand(-pi/2+i*dt), 0:steps)
            missing = 2pi*radius*p.bearing_length-sum(fitted.areas)
            @test 0 < missing < 1.001opening
            @test abs(missing/opening-1) < .007*(8/surface.patch_cells)^2
            println("Cylindrical $name: unported=$(length(plain.quads)), fitted=$(length(fitted.quads)); acceptance=$(fitted.acceptance); quality=$(fitted.quality)")
        end
    end
end

@testset "Cylindrical CLI from an unrelated working directory" begin
    project = normpath(joinpath(@__DIR__, "../.."))
    generator = joinpath(project, "mesh/2d/generate.jl")
    p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
        radial_clearance=.05, semicone_angle_deg=0., eccentricity_ratio=.6)
    mktempdir() do directory
        for variant in ("unported", "surface-feed")
            path = joinpath(directory, "$variant.vtk")
            args = ["--geometry", variant, "--journal-radius-mm", "50", "--length-mm", "100",
                "--clearance-mm", ".05", "--cone-angle-deg", "0", "--eccentricity-ratio", ".6",
                "--n-theta", "128", "--n-axial", "64", "--output", path]
            if variant == "surface-feed"
                append!(args, ["--feed-diameter-mm", "5", "--feed-z-mm", "50"])
            end
            command = `$(Base.julia_cmd()) --project=$project --startup-file=no --threads=1 --gcthreads=1 --check-bounds=yes $generator $args`
            command = addenv(Cmd(command; dir=directory), "OPENBLAS_NUM_THREADS"=>"1", "OMP_NUM_THREADS"=>"1")
            @test success(run(ignorestatus(command)))
            @test isfile(path)
            direct = variant == "unported" ? housing_grid(p,128,64) :
                surface_feed_grid(p,128,64; feed_diameter=5., feed_z=50.)
            expected = joinpath(directory, "$variant-direct.vtk")
            write_vtk(expected, direct)
            @test read(path) == read(expected)
            @test occursin("LINES ", read(path,String)) == (variant == "surface-feed")
        end
    end
end
