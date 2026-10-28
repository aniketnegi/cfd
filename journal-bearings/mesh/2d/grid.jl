# Periodic Quad4 grid on the housing surface.
# Lengths are mm, theta is in radians.
include("quality.jl")

"""
    background_grid(theta_nodes, z_nodes)

Build a periodic quadrilateral background from angular positions in radians and
axial positions in mm. Angles must increase strictly within [0, 2pi), without a
duplicate 2pi vertex; axial positions must be finite and increase strictly.

Returned nodes are (theta, z), not Cartesian points. This constructs only the
background vertices and cell connections, with angular index varying fastest.
"""
function background_grid(theta_nodes::AbstractVector{<:Real}, z_nodes::AbstractVector{<:Real})
    angles = collect(Float64, theta_nodes)
    axial_positions = collect(Float64, z_nodes)
    n_theta = length(angles)
    n_axial = length(axial_positions) - 1

    n_theta >= 3 || throw(ArgumentError("at least three angular positions are required"))
    n_axial >= 1 || throw(ArgumentError("at least two axial positions are required"))
    all(isfinite, angles) && all(theta -> 0 <= theta < 2pi, angles) &&
        all(>(0), diff(angles)) ||
        throw(ArgumentError("angular positions must be finite and strictly increasing within [0, 2pi)"))
    all(isfinite, axial_positions) && all(>(0), diff(axial_positions)) ||
        throw(ArgumentError("axial positions must be finite and strictly increasing"))

    node_id(i, j) = i + n_theta * (j - 1)
    nodes = Vector{NTuple{2,Float64}}(undef, n_theta * (n_axial + 1))
    quads = Vector{NTuple{4,Int}}(undef, n_theta * n_axial)

    for j in 1:(n_axial + 1), i in 1:n_theta
        nodes[node_id(i, j)] = (angles[i], axial_positions[j])
    end
    for j in 1:n_axial, i in 1:n_theta
        # Close circumference by connecting the last vertex to the first.
        next_i = i == n_theta ? 1 : i + 1
        cell = i + n_theta * (j - 1)
        quads[cell] = (
            node_id(i, j),
            node_id(next_i, j),
            node_id(next_i, j + 1),
            node_id(i, j + 1),
        )
    end

    return (; nodes, quads, n_theta, n_axial)
end

"""
    patch_axis(count, centre, extent, position, half_width, patch_cells)

Place a `patch_cells`-wide central patch inside position ± half_width, with a
surrounding support extending to position ± 2half_width. The remaining intervals
are uniform on either side. `centre` is a zero-based grid-node index.
"""
function patch_axis(count, centre, extent, position, half_width, patch_cells)
    support_start = centre - patch_cells
    support_end = centre + patch_cells
    0 < support_start < support_end < count ||
        throw(ArgumentError("support must lie inside grid"))
    0 < position - 2 * half_width < position + 2 * half_width < extent ||
        throw(ArgumentError("physical support crosses domain boundary"))

    axis = Vector{Float64}(undef, count + 1)
    for index in 0:count
        if index <= support_start
            coordinate = (position - 2 * half_width) * index / support_start
        elseif index <= support_end
            coordinate = position + 2 * half_width * (index - centre) / patch_cells
        else
            coordinate = position + 2 * half_width +
                (extent - position - 2 * half_width) * (index - support_end) / (count - support_end)
        end
        axis[index + 1] = coordinate
    end
    axis[1], axis[end] = 0.0, extent
    return axis
end

"""
    housing_grid(p, n_theta, n_axial; theta_nodes=nothing, z_nodes=nothing)

Map a periodic background onto the cylindrical/conical housing and calculate
cell radial gaps and physical areas. Spacing is uniform by default. Supply angular/axial position
vectors to adjust either direction: lengths must be n_theta and n_axial+1, with
z_nodes spanning exactly 0 to p.bearing_length.

The supplied positions are copied into the background, not modified. With
nonuniform spacing, each cell uses its own midpoint and angular/axial widths.
"""
function housing_grid(
    p::BearingParameters,
    n_theta::Int,
    n_axial::Int;
    theta_nodes=nothing,
    z_nodes=nothing,
)
    n_theta >= 3 || throw(ArgumentError("n_theta must be at least 3"))
    n_axial >= 1 || throw(ArgumentError("n_axial must be positive"))

    dtheta = 2pi / n_theta
    dz = p.bearing_length / n_axial
    ds = dz / cosd(p.semicone_angle_deg)

    # fn(i,j) => index for array
    node_id(i, j) = i + n_theta * (j - 1)

    uniform_theta = isnothing(theta_nodes)
    uniform_z = isnothing(z_nodes)
    if uniform_theta
        theta_nodes = [(i - 1) * dtheta for i in 1:n_theta]
    end
    if uniform_z
        z_nodes = [(j - 1) * dz for j in 1:(n_axial + 1)]
        # Anchor the ends even when multiplication rounds past bearing_length.
        z_nodes[end] = p.bearing_length
    end
    length(theta_nodes) == n_theta ||
        throw(DimensionMismatch("theta_nodes must contain n_theta positions"))
    length(z_nodes) == n_axial + 1 ||
        throw(DimensionMismatch("z_nodes must contain n_axial+1 positions"))

    background = background_grid(theta_nodes, z_nodes)
    first(background.nodes)[2] == 0 && last(background.nodes)[2] == p.bearing_length ||
        throw(ArgumentError("axial positions must span 0 to bearing_length"))

    # Map shared (theta, z) vertices onto the physical housing surface.
    nodes = [bearing_point(p, theta, z, 1.0) for (theta, z) in background.nodes]
    quads = background.quads
    radial_gaps = Vector{Float64}(undef, length(quads))
    areas = Vector{Float64}(undef, length(quads))

    for j in 1:n_axial
        z_low = background.nodes[node_id(1, j)][2]
        z_high = background.nodes[node_id(1, j + 1)][2]
        z = uniform_z ? (j - 0.5) * dz : (z_low + z_high) / 2
        cell_ds = uniform_z ? ds : (z_high - z_low) / cosd(p.semicone_angle_deg)
        for i in 1:n_theta
            k = i + n_theta * (j - 1)
            theta_low = background.nodes[node_id(i, j)][1]
            theta_high = i == n_theta ?
                background.nodes[node_id(1, j)][1] + 2pi :
                background.nodes[node_id(i + 1, j)][1]
            theta = uniform_theta ? (i - 0.5) * dtheta : (theta_low + theta_high) / 2
            cell_dtheta = uniform_theta ? dtheta : theta_high - theta_low
            geometry = bearing_radial_geometry(p, theta, z)
            radial_gaps[k] = geometry.radial_gap

            # Analytic cone area, not the area of the faceted output quad.
            areas[k] = geometry.housing_radius * cell_dtheta * cell_ds
        end
    end

    all(point -> all(isfinite, point), nodes) &&
        all(x -> isfinite(x) && x > 0, radial_gaps) &&
        all(x -> isfinite(x) && x > 0, areas) || error("invalid housing geometry, gaps or areas")
    quality = SurfaceMeshQuality.tensor_grid_quality(collect(Float64,theta_nodes),
        collect(Float64,z_nodes), bearing_geometry(p), cosd(p.semicone_angle_deg))
    acceptance = SurfaceMeshQuality.mesh_acceptance(quality; valid_geometry=true,
        placement_stable=true, termination_reason=:direct_construction)
    acceptance.accepted || error("unported mesh does not meet the surface-quality limits: $quality")
    # Retain the supplied tensor coordinates; solvers must not infer uniform spacing.
    theta_nodes = [background.nodes[i][1] for i in 1:n_theta]
    z_nodes = [background.nodes[node_id(1, j)][2] for j in 1:n_axial+1]
    return (; nodes, quads, radial_gaps, areas, n_theta, n_axial, quality, acceptance,
        parameters=p, theta_nodes, z_nodes)
end
