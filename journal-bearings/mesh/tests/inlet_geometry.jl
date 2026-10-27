using Test
using LinearAlgebra: cross, dot, norm
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))

# Independently lift a chart displacement into a Cartesian tangent (mm).
# This checks the metric-based measurements without repeating their 2D formula.
function cone_tangent(p, surface, point, displacement)
    theta = point[1] / surface.R
    z = point[2] * surface.ca * 1e3
    radius = p.housing_midpoint_radius + (p.half_bearing_length - z) * p.slope
    dtheta = displacement[1] / surface.R
    dz = displacement[2] * surface.ca * 1e3
    return [-radius*sin(theta)*dtheta - p.slope*cos(theta)*dz,
             radius*cos(theta)*dtheta - p.slope*sin(theta)*dz, dz]
end

# Integrate the circular drilling footprint directly on the cone, independently
# of chart polygons and their centroid-based physical-area conversion.
function drilling_area(p, diameter, feed_z)
    radius = diameter / 2
    intervals = 10_000
    dz = diameter / intervals
    return sum(1:intervals) do interval
        z = feed_z - radius + (interval - .5) * dz
        housing_radius = p.housing_midpoint_radius + (p.half_bearing_length - z) * p.slope
        half_width = sqrt(radius^2 - (z-feed_z)^2)
        housing_radius * 2asin(half_width/housing_radius) * dz / cosd(p.semicone_angle_deg)
    end
end

@testset "Supplied geometry reaches every mesh stage" begin
    reference = (journal_midpoint_radius=50., bearing_length=100., radial_clearance=.05,
        semicone_angle_deg=10., eccentricity_ratio=.6)
    cases = (
        ("reference", reference, 5., 50., (128,64)),
        ("radius-40", (; reference..., journal_midpoint_radius=40.), 5., 50., (128,64)),
        ("radius-60", (; reference..., journal_midpoint_radius=60.), 5., 50., (128,64)),
        ("length-80", (; reference..., bearing_length=80.), 5., 40., (128,64)),
        ("length-120", (; reference..., bearing_length=120.), 5., 60., (128,64)),
        ("decimal-length", (; reference..., bearing_length=123.4), 5., 61.7, (128,64)),
        ("decimal-length-angle-15", (; reference..., bearing_length=123.4,
            semicone_angle_deg=15.), 5., 61.7, (128,64)),
        ("clearance-0p1", (; reference..., radial_clearance=.1), 5., 50., (128,64)),
        ("zero-angle limit", (; reference..., semicone_angle_deg=0.), 5., 50., (128,64)),
        ("angle-5", (; reference..., semicone_angle_deg=5.), 5., 50., (128,64)),
        ("angle-15", (; reference..., semicone_angle_deg=15.), 5., 50., (128,64)),
        ("diameter-3", reference, 3., 50., (128,64)),
        ("diameter-6", reference, 6., 50., (128,64)),
        ("diameter-8", reference, 8., 50., (128,64)),
        ("inlet-40", reference, 5., 40., (128,64)),
        ("inlet-60", reference, 5., 60., (128,64)),
        ("inlet-43p7", reference, 5., 43.7, (128,64)),
        ("combined-odd", (journal_midpoint_radius=45., bearing_length=90.,
            radial_clearance=.06, semicone_angle_deg=12., eccentricity_ratio=.6),
            4., 35.3, (127,65)),
        ("combined", (journal_midpoint_radius=60., bearing_length=120.,
            radial_clearance=.08, semicone_angle_deg=15., eccentricity_ratio=.6),
            6., 70.3, (160,64)),
    )
    for (name, dimensions, diameter, feed_z, counts) in cases
        @testset "$name" begin
            p = BearingParameters(; dimensions...)
            grid = surface_feed_grid(p, counts...; feed_diameter=diameter, feed_z)
            surface = grid.surface
            # Rebuilding geometry after optimization must not restore defaults.
            @test surface.R == p.housing_midpoint_radius * 1e-3
            @test surface.L == p.bearing_length * 1e-3
            @test surface.ca == cosd(p.semicone_angle_deg)
            @test surface.slope === p.slope
            @test surface.feed_z_m == feed_z * 1e-3
            @test surface.inlet_radius_m == diameter/2 * 1e-3
            @test surface.period == 2pi * surface.R
            @test (surface.n_theta, surface.n_axial) == counts
            @test grid.quality.minimum_scaled_corner > .9
            @test grid.quality.minimum_face_alignment > .9
            @test grid.quality.maximum_neighbour_area_ratio < 2
            @test grid.quality.minimum_scaled_corner == minimum(surface.corner_scores)
            @test grid.quality.minimum_face_alignment == minimum(surface.face_alignment)

            @test all(grid.nodes) do (x,y,z)
                radius = p.housing_midpoint_radius + (p.half_bearing_length-z)*p.slope
                isapprox(hypot(x,y), radius; atol=1e-10, rtol=0) &&
                    0 <= z <= p.bearing_length
            end
            @test minimum(node[3] for node in grid.nodes) == 0
            @test maximum(node[3] for node in grid.nodes) == p.bearing_length
            @test all(grid.feed_edges) do edge
                all(edge) do node
                    x,y,z = grid.nodes[node]
                    y > 0 && isapprox(x^2+(z-feed_z)^2, (diameter/2)^2; atol=1e-10, rtol=0)
                end
            end
            rim_z = [grid.nodes[node][3] for node in unique(collect(Iterators.flatten(grid.feed_edges)))]
            @test minimum(rim_z) ≈ feed_z-diameter/2 atol=1e-10
            @test maximum(rim_z) ≈ feed_z+diameter/2 atol=1e-10
            @test all(eachindex(grid.quads)) do cell
                theta = surface.centre[cell,1] / surface.R
                z = surface.centre[cell,2] * surface.ca * 1e3
                journal_radius = p.journal_midpoint_radius + (p.half_bearing_length-z)*p.slope
                journal_ray = p.eccentricity*cos(theta) +
                    sqrt(journal_radius^2-(p.eccentricity*sin(theta))^2)
                expected = journal_radius + p.radial_clearance - journal_ray
                isapprox(grid.radial_gaps[cell], expected; atol=1e-10, rtol=0)
            end
            whole_area = 2pi*p.housing_midpoint_radius*p.bearing_length/cosd(p.semicone_angle_deg)
            opening_area = whole_area - sum(grid.areas)
            true_area = drilling_area(p, diameter, feed_z)
            @test 0 < opening_area < 1.001true_area
            @test abs(opening_area/true_area-1) < .007*(8/surface.patch_cells)^2

            # Probe both the distant background and the inlet cells using the
            # actual 3D cone tangents, rather than the implementation's metric.
            feed_faces = findall(==(:feed), surface.boundary)
            cells = unique(vcat(round.(Int, range(1,length(grid.quads); length=12)),
                surface.owner[feed_faces]))
            for cell in cells
                coordinates = surface.localnodes[cell]
                scores = Float64[]
                for corner in 1:4
                    point = coordinates[corner,:]
                    theta = point[1] / surface.R
                    normal = [cos(theta), sin(theta), p.slope]
                    next = cone_tangent(p, surface, point,
                        coordinates[mod1(corner+1,4),:]-point)
                    previous = cone_tangent(p, surface, point,
                        coordinates[mod1(corner-1,4),:]-point)
                    push!(scores, dot(cross(next,previous),normal)/(norm(next)*norm(previous)*norm(normal)))
                end
                @test surface.corner_scores[cell] ≈ minimum(scores) atol=1e-12 rtol=0
            end
            faces = unique(vcat(round.(Int, range(1,length(surface.owner); length=20)), feed_faces))
            for face in faces
                point = surface.mid[face]
                theta = point[1] / surface.R
                normal = [cos(theta), sin(theta), p.slope]
                chart_normal = surface.normals[face]
                tangent = cone_tangent(p, surface, point, [-chart_normal[2], chart_normal[1]])
                physical_normal = cross(tangent,normal)
                neighbour = surface.neighbour[face]
                connection = (neighbour == 0 ? point : surface.centre[neighbour,:]) -
                    surface.centre[surface.owner[face],:]
                connection[1] = mod(connection[1]+surface.period/2,surface.period)-surface.period/2
                direction = cone_tangent(p, surface, point, connection)
                alignment = dot(direction,physical_normal)/(norm(direction)*norm(physical_normal))
                @test surface.face_alignment[face] ≈ alignment atol=1e-12 rtol=0
            end

            initial = ControlledSurfaceInlet.reference_seed(p, counts...;
                feed_diameter=diameter, feed_z)
            selection = ControlledSurfaceInlet.macro_data(initial)
            half_width_m = 4diameter * 1e-3
            expected_free = findall(eachindex(initial.nodes)) do node
                x,y = initial.nodes[node]
                selection.valence[node] > 0 && !(node in selection.boundary) &&
                    abs(x-p.housing_midpoint_radius*1e-3*pi/2) < half_width_m &&
                    abs(y-feed_z*1e-3/cosd(p.semicone_angle_deg)) < half_width_m/cosd(p.semicone_angle_deg)
            end
            @test selection.free == expected_free
            free_set = Set(expected_free)
            @test all(node -> node in free_set || surface.nodes[node] == initial.nodes[node],
                eachindex(initial.nodes))
        end
    end
end

@testset "Eccentricity changes only gaps for nonreference geometry" begin
    dimensions = (journal_midpoint_radius=45., bearing_length=90., radial_clearance=.06,
        semicone_angle_deg=12.)
    centred = BearingParameters(; dimensions..., eccentricity_ratio=0.)
    offset = BearingParameters(; dimensions..., eccentricity_ratio=.7)
    first = surface_feed_grid(centred,127,65; feed_diameter=4.,feed_z=35.3)
    second = surface_feed_grid(offset,127,65; feed_diameter=4.,feed_z=35.3)
    @test first.nodes == second.nodes
    @test first.quads == second.quads
    @test first.feed_edges == second.feed_edges
    @test first.surface == second.surface
    @test first.areas == second.areas
    @test first.quality == second.quality
    @test first.radial_gaps != second.radial_gaps
    @test all(isapprox.(first.radial_gaps,.06; atol=1e-11, rtol=0))
    @test all(offset.minimum_gap-1e-11 .<= second.radial_gaps .<= offset.maximum_gap+1e-11)
end
