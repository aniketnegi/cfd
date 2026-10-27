using Test
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))

@testset "Requested counts construct the actual inlet background" begin
    p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
        radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
    for (n_theta, n_axial) in ((127, 65), (160, 52), (192, 64))
        @testset "$n_theta×$n_axial" begin
            grid = surface_feed_grid(p, n_theta, n_axial; feed_diameter=5., feed_z=50.)
            surface = grid.surface
            patch_cells = surface.patch_cells
            layers = surface.transition_layers
            @test (grid.n_theta, grid.n_axial) == (n_theta, n_axial)
            @test (surface.n_theta, surface.n_axial) == (n_theta, n_axial)
            @test iseven(patch_cells) && patch_cells >= 2
            @test layers >= 1
            @test length(grid.quads) == n_theta*n_axial - patch_cells^2 + 4patch_cells*layers
            @test length(grid.feed_edges) == 4patch_cells
            @test length(grid.nodes) == length(grid.quads) + n_theta + 2patch_cells - 1
            @test count(==(:large), surface.boundary) == n_theta
            @test count(==(:small), surface.boundary) == n_theta
            # The fixed periodic cut has one vertex per axial position.
            @test count(node -> node[1] > 0 && abs(node[2]) < 1e-12, grid.nodes) == n_axial + 1
            @test length(grid.nodes) == length(unique(grid.nodes))
            @test minimum(grid.areas) > 0
            @test minimum(grid.radial_gaps) > 0
            @test grid.quality.minimum_scaled_corner > .9
            @test grid.quality.minimum_face_alignment > .9
            @test grid.quality.maximum_neighbour_area_ratio < 2

            edge_counts = Dict{Tuple{Int,Int},Int}()
            for quad in grid.quads, corner in 1:4
                edge = minmax(quad[corner], quad[mod1(corner+1,4)])
                edge_counts[edge] = get(edge_counts,edge,0) + 1
            end
            @test all(count -> count in (1,2), values(edge_counts))
            @test count(==(1), values(edge_counts)) == 2n_theta + length(grid.feed_edges)
            @test all(edge -> edge_counts[edge] == 1, grid.feed_edges)
            @test length(grid.nodes) - length(edge_counts) + length(grid.quads) == -1
            @test all(grid.feed_edges) do edge
                all(edge) do node
                    x,y,z = grid.nodes[node]
                    y > 0 && isapprox(x^2+(z-50.)^2, 2.5^2; atol=1e-10)
                end
            end
        end
    end
    @test_throws ArgumentError surface_feed_grid(p, 8, 4; feed_diameter=5., feed_z=50.)
    @test_throws r"does not meet the surface-quality limits" surface_feed_grid(p, 12, 6;
        feed_diameter=5., feed_z=50.)
end

@testset "Cached objective and analytic gradient" begin
    p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
        radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
    mesh = ControlledSurfaceInlet.reference_seed(p, 64, 32; feed_diameter=5., feed_z=50.)
    selection = ControlledSurfaceInlet.macro_data(mesh)
    original_nodes = copy(selection.nodes)
    variables = vec(copy(selection.nodes[selection.free,:]))
    energy, gradient = ControlledSurfaceInlet.layout_energy(mesh, selection, variables)
    @test isfinite(energy)
    @test all(isfinite, gradient)
    @test length(gradient) == length(variables)
    @test length(selection.initial_areas) == length(mesh.quads)
    @test all(selection.active_faces) do face
        mesh.neighbour[face] > 0
    end

    # Check the derivative without reusing its algebraic formula.
    delta = 1e-8
    for index in round.(Int, range(1,length(variables); length=8))
        plus = copy(variables)
        minus = copy(variables)
        plus[index] += delta
        minus[index] -= delta
        plus_energy = first(ControlledSurfaceInlet.layout_energy(mesh, selection, plus))
        minus_energy = first(ControlledSurfaceInlet.layout_energy(mesh, selection, minus))
        finite_difference = (plus_energy-minus_energy)/(2delta)
        @test isapprox(gradient[index], finite_difference; rtol=1e-5, atol=5e-5)
    end
    @test selection.nodes == original_nodes

    invalid = copy(variables)
    invalid[1] += 1.
    @test !isfinite(first(ControlledSurfaceInlet.layout_energy(mesh, selection, invalid)))
end
