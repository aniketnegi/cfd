# Shared surface metrics and acceptance for fitted and unported cone meshes.
module SurfaceMeshQuality
using LinearAlgebra
import ..housing_radius
include("metrics.jl")

function passes_quality_limits(quality)
    return quality.minimum_scaled_corner > 0.9 &&
        quality.minimum_face_alignment > 0.9 &&
        quality.maximum_neighbour_area_ratio < 2
end

function mesh_acceptance(quality; valid_geometry, placement_stable, termination_reason)
    quality_pass = passes_quality_limits(quality)
    return (; valid_geometry, quality_pass, placement_stable, termination_reason,
        accepted=valid_geometry && quality_pass && placement_stable)
end

# Exact cone chord, evaluated stably even for tiny movements and seam edges.
# Chart lengths/metres; the result is physical metres, not chart distance.
function physical_distance(mesh, first, second)
    dx, dy = second[1] - first[1], second[2] - first[2]
    abs(dx) > mesh.period / 2 && (dx -= copysign(mesh.period, dx))
    z0 = chart_axial_position_mm(mesh.geometry_mm, mesh.ca, first[2])
    z1 = chart_axial_position_mm(mesh.geometry_mm, mesh.ca, second[2])
    r0 = housing_radius(mesh.geometry_mm, z0) * 1e-3
    r1 = housing_radius(mesh.geometry_mm, z1) * 1e-3
    return hypot(dy, 2 * sqrt(r0 * r1) * sin(dx / (2 * mesh.R)))
end

function local_edge_lengths(mesh)
    lengths = fill(Inf, length(mesh.nodes))
    for (a, b) in mesh.edge_nodes
        distance = physical_distance(mesh, mesh.nodes[a], mesh.nodes[b])
        lengths[a] = min(lengths[a], distance)
        lengths[b] = min(lengths[b], distance)
    end
    return lengths
end

# Every vertex's path length matters: back-and-forth movement must not cancel.
function accumulate_motion!(motion, mesh, previous, current, local_lengths)
    count = length(motion)
    for node in 1:count
        first = (previous[node], previous[node + count])
        second = (current[node], current[node + count])
        motion[node] += physical_distance(mesh, first, second) / local_lengths[node]
    end
    return maximum(motion; init=0.0)
end

"""
    tensor_grid_quality(theta, z, geometry, ca)

Analytic form of the SAME surface metrics for an orthogonal tensor chart.
Corners and face connections are orthogonal. Chart areas are proportional to
angular width times axial width; the cone metric scales circumferential edges
by the housing radius at each endpoint. No optimizer or full face assembly.
"""
function tensor_grid_quality(theta, z, geometry, ca)
    widths = diff([theta; first(theta) + 2pi])
    all(w -> 0 < w < pi, widths) ||
        throw(ArgumentError("each angular cell must span less than pi radians"))
    heights = diff(z)
    ratios = [max(widths[i], widths[mod1(i+1,length(widths))]) /
        min(widths[i], widths[mod1(i+1,length(widths))]) for i in eachindex(widths)]
    for j in 1:length(heights)-1
        push!(ratios, max(heights[j], heights[j+1]) / min(heights[j], heights[j+1]))
    end
    edge_ratio = 1.0
    for j in eachindex(heights), width in widths
        edges = (housing_radius(geometry,z[j])*width,
            housing_radius(geometry,z[j+1])*width, heights[j]/ca)
        all(x -> isfinite(x) && x > 0, edges) || error("invalid physical cell edges")
        edge_ratio = max(edge_ratio, maximum(edges)/minimum(edges))
    end
    return (; minimum_scaled_corner=1.0, maximum_edge_ratio=edge_ratio,
        minimum_face_alignment=1.0, maximum_nonorthogonality_deg=0.0,
        maximum_neighbour_area_ratio=maximum(ratios),
        cells=length(widths)*length(heights), rim_edges=0)
end
end
