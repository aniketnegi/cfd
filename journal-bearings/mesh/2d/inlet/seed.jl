# Initial inlet layout: replace a background patch with rings of quadrilaterals.
# Construction uses (theta, z) with z in mm, then converts to metre-based chart
# coordinates for node placement.

# Recover the +Y housing point from its XZ projection; return (theta, z).
function housing_coordinates(parameters, x::Float64, z::Float64)
    radius = housing_radius(parameters, z)
    abs(x) < radius || throw(ArgumentError("inlet point lies outside housing"))
    y = sqrt(radius^2 - x^2)
    return (atan(y, x), z)
end

# Visit the rectangle's boundary once, without repeating the corner vertices.
function rectangle_perimeter(ids::Matrix{Int}, start_i::Int, start_j::Int, cells_per_side::Int)
    perimeter = Int[]
    append!(perimeter, (ids[start_i + step, start_j] for step in 0:(cells_per_side - 1)))
    append!(perimeter, (ids[start_i + cells_per_side, start_j + step] for step in 0:(cells_per_side - 1)))
    append!(perimeter, (ids[start_i + cells_per_side - step, start_j + cells_per_side] for step in 0:(cells_per_side - 1)))
    append!(perimeter, (ids[start_i, start_j + cells_per_side - step] for step in 0:(cells_per_side - 1)))
    return perimeter
end

"""
    harmonic_support!(nodes, grid_ids, parameters, centre_i, centre_j, patch_cells)

Smooth background nodes between the fixed patch and support perimeters. Each
unknown XZ position is the average of its four neighbours: 4x - sum(x_neighbour)=0.
Solve for projected x and z together, then map the points back onto the cone.
"""
function harmonic_support!(nodes, grid_ids, parameters, centre_i, centre_j, patch_cells)
    unknown_indices = [
        (i, j)
        for j in (centre_j - patch_cells + 1):(centre_j + patch_cells - 1),
            i in (centre_i - patch_cells + 1):(centre_i + patch_cells - 1)
        if !(abs(i - centre_i) <= patch_cells ÷ 2 && abs(j - centre_j) <= patch_cells ÷ 2)
    ]
    row_for_index = Dict(index => row for (row, index) in enumerate(unknown_indices))
    rows = Int[]
    columns = Int[]
    values = Float64[]
    boundary_positions = zeros(length(unknown_indices), 2)

    for (row, (i, j)) in enumerate(unknown_indices)
        push!(rows, row)
        push!(columns, row)
        push!(values, 4.0)
        for index in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
            neighbour_row = get(row_for_index, index, 0)
            if neighbour_row == 0
                theta, z = nodes[grid_ids[index...]]
                boundary_positions[row, 1] += housing_radius(parameters, z) * cos(theta)
                boundary_positions[row, 2] += z
            else
                push!(rows, row)
                push!(columns, neighbour_row)
                push!(values, -1.0)
            end
        end
    end

    matrix = sparse(rows, columns, values, length(unknown_indices), length(unknown_indices))
    positions = matrix \ boundary_positions
    for (row, index) in enumerate(unknown_indices)
        nodes[grid_ids[index...]] = housing_coordinates(parameters, positions[row, 1], positions[row, 2])
    end
    return nothing
end

"""
    reference_seed(p, n_theta, n_axial; feed_diameter, feed_z)

Build the initial cone mesh around a +Y circular drilling. Bearing dimensions,
diameter and axial position use mm. The patch count follows physical background
spacing; layers retain the chosen nine-to-eight layer/patch ratio.

The opening must fit between the ends and leave room for the unchanged ring
and smoothing construction. Invalid requests are rejected, not moved or resized.
"""
function reference_seed(parameters::BearingParameters, n_theta::Int, n_axial::Int;
    feed_diameter::Real, feed_z::Real)
    n_theta >= 3 || throw(ArgumentError("n_theta must be at least 3"))
    n_axial >= 1 || throw(ArgumentError("n_axial must be positive"))
    diameter = Float64(feed_diameter)
    inlet_z = Float64(feed_z)
    isfinite(diameter) && diameter > 0 ||
        throw(ArgumentError("feed_diameter must be positive and finite"))
    isfinite(inlet_z) || throw(ArgumentError("feed_z must be finite"))
    inlet_radius = diameter / 2
    inlet_radius < inlet_z < parameters.bearing_length - inlet_radius ||
        throw(ArgumentError("the complete inlet opening must lie strictly between the bearing ends"))
    patch_half_width = 1.5 * diameter
    control_radius_factor = 2.87
    control_square_blend = 0.36

    # The adjusted axes place these nodes exactly at the inlet centre, even
    # when the requested counts would not put a uniform-grid node there.
    centre_i = n_theta ÷ 4 + 1
    centre_j = round(Int, inlet_z * n_axial / parameters.bearing_length) + 1

    cos_angle = cosd(parameters.semicone_angle_deg)
    circumferential_step = 2pi * housing_radius(parameters, inlet_z) / n_theta
    meridional_step = parameters.bearing_length / (n_axial * cos_angle)
    background_cell_area = circumferential_step * meridional_step
    patch_area = (2 * patch_half_width)^2 / cos_angle
    desired_patch_cells = sqrt(patch_area / background_cell_area)

    # An even patch count keeps the inlet-centre node between opposite sides.
    # Limit it to leave background cells beyond the smoothing support.
    available_cells = minimum((centre_i - 2, n_theta - centre_i - 1,
        centre_j - 2, n_axial - centre_j))
    maximum_patch_cells = 2 * (available_cells ÷ 2)
    maximum_patch_cells >= 2 || throw(ArgumentError(
        "background $n_theta×$n_axial leaves insufficient room for the inlet transition",
    ))
    patch_cells = clamp(2 * round(Int, desired_patch_cells / 2), 2, maximum_patch_cells)
    transition_layers = max(1, round(Int, 9 * patch_cells / 8))

    start_i = centre_i - patch_cells ÷ 2
    start_j = centre_j - patch_cells ÷ 2
    1 <= start_i && start_i + patch_cells <= n_theta ||
        throw(ArgumentError("inlet patch crosses seam"))
    1 <= start_j && start_j + patch_cells <= n_axial + 1 ||
        throw(ArgumentError("inlet patch crosses end"))
    1 < centre_i - patch_cells && centre_i + patch_cells < n_theta ||
        throw(ArgumentError("support crosses seam"))
    1 < centre_j - patch_cells && centre_j + patch_cells < n_axial + 1 ||
        throw(ArgumentError("support crosses end"))
    isfinite(patch_half_width) && patch_half_width > inlet_radius ||
        throw(ArgumentError("patch_half_width must be finite and exceed inlet radius"))

    theta_axis = patch_axis(n_theta, centre_i - 1, 2pi, pi/2,
        patch_half_width / housing_radius(parameters, inlet_z), patch_cells)[1:end-1]
    z_axis = patch_axis(n_axial, centre_j - 1, parameters.bearing_length,
        inlet_z, patch_half_width, patch_cells)

    # Check the smaller-radius side too: a cone's mouth clearance varies with z.
    minimum(housing_radius(parameters, z) * sin(pi/2 - theta_axis[start_i])
        for z in (z_axis[start_j], z_axis[start_j + patch_cells])) > inlet_radius &&
        min(inlet_z - z_axis[start_j], z_axis[start_j + patch_cells] - inlet_z) > inlet_radius ||
        throw(ArgumentError("background patch does not enclose inlet"))

    background = background_grid(theta_axis, z_axis)
    surface_nodes = NTuple{2,Float64}[]
    retained_ids = zeros(Int, length(background.nodes))
    grid_ids = reshape(retained_ids, n_theta, n_axial + 1)
    for (node, position) in enumerate(background.nodes)
        i = mod1(node, n_theta)
        j = (node - 1) ÷ n_theta + 1
        # Keep patch perimeter vertices, but remove its old interior vertices.
        if !(start_i < i < start_i + patch_cells && start_j < j < start_j + patch_cells)
            push!(surface_nodes, position)
            retained_ids[node] = length(surface_nodes)
        end
    end

    surface_quads = NTuple{4,Int}[]
    for (cell, quad) in enumerate(background.quads)
        i = mod1(cell, n_theta)
        j = (cell - 1) ÷ n_theta + 1
        in_patch = start_i <= i < start_i + patch_cells && start_j <= j < start_j + patch_cells
        in_patch && continue
        push!(surface_quads, ntuple(corner -> retained_ids[quad[corner]], 4))
    end

    control_ring = rectangle_perimeter(grid_ids, start_i, start_j, patch_cells)
    rim_segments = 4 * patch_cells

    # Round the square patch boundary without changing its vertex connections.
    # Blend 0 would be circular in XZ; blend 1 would retain the original square.
    for (segment, node) in enumerate(control_ring)
        angle = 5pi/4 + 2pi * (segment - 1) / rim_segments
        theta, z = surface_nodes[node]
        x = housing_radius(parameters, z) * cos(theta)
        control_radius = control_radius_factor * inlet_radius
        x = (1 - control_square_blend) * (-control_radius * cos(angle)) + control_square_blend * x
        z = (1 - control_square_blend) * (inlet_z + control_radius * sin(angle)) + control_square_blend * z
        hypot(x, z - inlet_z) > inlet_radius ||
            throw(ArgumentError("control ring must lie outside inlet"))
        surface_nodes[node] = housing_coordinates(parameters, x, z)
    end
    harmonic_support!(surface_nodes, grid_ids, parameters, centre_i, centre_j, patch_cells)

    # Rim vertices satisfy x² + (z-inlet_z)² = inlet_radius² on the +Y cone.
    # The starting angle and traversal match the control ring's vertex order.
    rim_ring = Int[]
    for segment in 0:(rim_segments - 1)
        angle = 5pi/4 + 2pi * segment / rim_segments
        x = inlet_radius * cos(angle)
        z = inlet_z + inlet_radius * sin(angle)
        push!(surface_nodes, housing_coordinates(parameters, -x, z))
        push!(rim_ring, length(surface_nodes))
    end

    function interpolate_ring(inner::Vector{Int}, outer::Vector{Int}, fraction::Float64)
        ring = Int[]
        for (inner_node, outer_node) in zip(inner, outer)
            inner_theta, inner_z = surface_nodes[inner_node]
            outer_theta, outer_z = surface_nodes[outer_node]
            inner_x = housing_radius(parameters, inner_z) * cos(inner_theta)
            outer_x = housing_radius(parameters, outer_z) * cos(outer_theta)
            x = (1 - fraction) * inner_x + fraction * outer_x
            z = (1 - fraction) * inner_z + fraction * outer_z
            push!(surface_nodes, housing_coordinates(parameters, x, z))
            push!(ring, length(surface_nodes))
        end
        return ring
    end

    # Each pair of successive rings bounds one cell layer.
    rings = Vector{Int}[rim_ring]
    for layer in 1:(transition_layers - 1)
        push!(rings, interpolate_ring(rim_ring, control_ring, layer / transition_layers))
    end
    push!(rings, control_ring)
    for layer in 1:transition_layers, segment in 1:rim_segments
        inner = rings[layer]
        outer = rings[layer + 1]
        next_segment = mod1(segment + 1, rim_segments)
        push!(surface_quads, (inner[segment], outer[segment], outer[next_segment], inner[next_segment]))
    end

    rim = Set(minmax(rim_ring[segment], rim_ring[mod1(segment + 1, rim_segments)])
        for segment in 1:rim_segments)
    reference_radius = parameters.housing_midpoint_radius * 1e-3
    length_m = parameters.bearing_length * 1e-3
    nodes = [(reference_radius * theta, z * 1e-3 / cos_angle) for (theta, z) in surface_nodes]
    mesh = assemble_surface_mesh(nodes, surface_quads, rim;
        R=reference_radius, L=length_m, ca=cos_angle, slope=parameters.slope,
        geometry_mm=bearing_geometry(parameters),
        feed_z_m=inlet_z * 1e-3, inlet_radius_m=inlet_radius * 1e-3)
    return (; mesh..., n_theta, n_axial, patch_cells, transition_layers)
end
