# Improve node positions on the requested mesh.
# Coordinates are (R_ref * theta, z / cos(angle)), in metres.

# Translate each cell before measuring it, and unwrap cells crossing the seam.
function cell_coordinates(mesh, positions, quad)
    coordinates = Matrix(positions[collect(quad), :])
    if maximum(coordinates[:, 1]) - minimum(coordinates[:, 1]) > mesh.period / 2
        coordinates[coordinates[:, 1] .< mesh.period / 2, 1] .+= mesh.period
    end
    coordinates .-= coordinates[1:1, :]
    return coordinates
end

function cell_chart_area(coordinates)
    return sum(
        coordinates[corner, 1] * coordinates[mod1(corner + 1, 4), 2] -
        coordinates[corner, 2] * coordinates[mod1(corner + 1, 4), 1]
        for corner in 1:4
    ) / 2
end

"""
    macro_data(mesh)

Select movable nodes and the cells affected by them. Keep the inlet rim, axial
ends, periodic seam and distant background fixed. For the 5 mm inlet, interior
nodes may move within an approximately ±20 mm neighbourhood of its centre.
"""
function macro_data(mesh)
    nodes = hcat([collect(point) for point in mesh.nodes]...)'
    nodes = Matrix(nodes)
    edge_counts = Dict{Tuple{Int,Int},Int}()
    node_valence = zeros(Int, length(mesh.nodes))
    for quad in mesh.quads
        for corner in 1:4
            node_valence[quad[corner]] += 1
            edge = minmax(quad[corner], quad[mod1(corner + 1, 4)])
            edge_counts[edge] = get(edge_counts, edge, 0) + 1
        end
    end

    boundary_nodes = Set(node for (edge, count) in edge_counts if count == 1 for node in edge)
    rim_edges = Set(edge for (edge, count) in edge_counts if count == 1 &&
        all(node -> 1e-10 < nodes[node, 2] < mesh.L / mesh.ca - 1e-10, edge))
    movable_half_width = 0.020 * (mesh.inlet_radius_m / 0.0025)
    movable_nodes = [
        node for node in eachindex(node_valence)
        if node_valence[node] > 0 && !(node in boundary_nodes) &&
            abs(nodes[node, 1] - mesh.R * pi/2) < movable_half_width &&
            abs(nodes[node, 2] - mesh.feed_z_m / mesh.ca) < movable_half_width / mesh.ca
    ]
    movable_node_set = Set(movable_nodes)
    active_cells = findall(quad -> any(node -> node in movable_node_set, quad), mesh.quads)
    active_mask = falses(length(mesh.quads))
    active_mask[active_cells] .= true
    active_faces = [
        face for face in eachindex(mesh.owner)
        if mesh.neighbour[face] > 0 &&
            (active_mask[mesh.owner[face]] || active_mask[mesh.neighbour[face]])
    ]

    # Cells without movable vertices never change during optimization. Cache
    # their areas once instead of measuring the entire fine mesh at every trial.
    initial_areas = [cell_chart_area(cell_coordinates(mesh, nodes, quad)) for quad in mesh.quads]
    all(>(0), initial_areas) || error("nonpositive initial cell area")

    return (;
        nodes,
        free=movable_nodes,
        active=active_cells,
        active_faces,
        initial_areas,
        rim=rim_edges,
        valence=node_valence,
        boundary=boundary_nodes,
    )
end

"""
    layout_energy(mesh, selection, variables; area_weight=0.25, anchor_weight=0.01)

Return the mesh objective and its gradient with respect to movable coordinates.
Penalize corner cotangents squared, neighbouring log-area differences, and smaller
log-area changes relative to the initial mesh. Weights are construction choices.
An inverted trial cell returns infinite energy so the line search rejects it.
"""
function layout_energy(mesh, selection, variables; area_weight=0.25, anchor_weight=0.01)
    positions = copy(selection.nodes)
    positions[selection.free, :] = reshape(variables, :, 2)
    gradient = zeros(size(positions))
    areas = copy(selection.initial_areas)
    area_gradient = zeros(length(mesh.quads))
    energy = 0.0

    for cell in selection.active
        quad = mesh.quads[cell]
        coordinates = cell_coordinates(mesh, positions, quad)
        area = cell_chart_area(coordinates)
        area > 0 || return (Inf, zeros(length(variables)))
        areas[cell] = area

        for corner in 1:4
            next_corner = mod1(corner + 1, 4)
            previous_corner = mod1(corner - 1, 4)
            next_edge = coordinates[next_corner, :] - coordinates[corner, :]
            previous_edge = coordinates[previous_corner, :] - coordinates[corner, :]
            corner_cross = next_edge[1] * previous_edge[2] - next_edge[2] * previous_edge[1]
            edge_dot = dot(next_edge, previous_edge)
            corner_cross > 0 || return (Inf, zeros(length(variables)))
            energy += (edge_dot / corner_cross)^2

            # Differentiate cot(angle)^2 with respect to both incident edges;
            # the corner vertex receives the opposite of their combined change.
            next_edge_gradient = 2 * edge_dot / corner_cross^2 * previous_edge -
                2 * edge_dot^2 / corner_cross^3 * [previous_edge[2], -previous_edge[1]]
            previous_edge_gradient = 2 * edge_dot / corner_cross^2 * next_edge -
                2 * edge_dot^2 / corner_cross^3 * [-next_edge[2], next_edge[1]]
            gradient[quad[next_corner], :] .+= next_edge_gradient
            gradient[quad[previous_corner], :] .+= previous_edge_gradient
            gradient[quad[corner], :] .-= next_edge_gradient + previous_edge_gradient
        end
    end

    # Area ratios discourage an abrupt jump in cell size across a shared face.
    for face in selection.active_faces
        owner_cell = mesh.owner[face]
        neighbour_cell = mesh.neighbour[face]
        log_ratio = log(areas[owner_cell] / areas[neighbour_cell])
        energy += area_weight * log_ratio^2
        area_gradient[owner_cell] += 2 * area_weight * log_ratio / areas[owner_cell]
        area_gradient[neighbour_cell] -= 2 * area_weight * log_ratio / areas[neighbour_cell]
    end
    for cell in selection.active
        # The smaller anchor discourages collapse/drift without fixing cell sizes.
        log_ratio = log(areas[cell] / mesh.areas[cell])
        energy += anchor_weight * log_ratio^2
        area_gradient[cell] += 2 * anchor_weight * log_ratio / areas[cell]
    end

    # An inactive cell has no movable vertices, so its coordinate derivatives
    # cannot contribute to the returned gradient.
    for cell in selection.active
        quad = mesh.quads[cell]
        area_gradient[cell] == 0 && continue
        for corner in 1:4
            next_corner = mod1(corner + 1, 4)
            previous_corner = mod1(corner - 1, 4)
            dx = positions[quad[previous_corner], 1] - positions[quad[next_corner], 1]
            dx = mod(dx + mesh.period / 2, mesh.period) - mesh.period / 2
            gradient[quad[corner], 1] += area_gradient[cell] *
                (positions[quad[next_corner], 2] - positions[quad[previous_corner], 2]) / 2
            gradient[quad[corner], 2] += area_gradient[cell] * dx / 2
        end
    end
    return energy, vec(gradient[selection.free, :])
end

"""
    controlled_layout(mesh; max_iterations=180)

Improve movable node positions using limited-memory BFGS and backtracking.
Keep up to eight position/gradient-change pairs, and limit each node's trial
motion to 0.2 mm for the 5 mm inlet. Return the last accepted layout and history;
reaching the iteration cap is not a claim that the optimum was found.
"""
function controlled_layout(mesh; max_iterations=180)
    # The 5 mm reference normalizes the existing numerical step/tolerance
    # settings; it does not fix the actual inlet diameter.
    diameter_scale = mesh.inlet_radius_m / 0.0025
    selection = macro_data(mesh)
    variables = vec(copy(selection.nodes[selection.free, :]))
    energy, gradient = layout_energy(mesh, selection, variables)
    position_changes = Vector{Float64}[]
    gradient_changes = Vector{Float64}[]
    history = NamedTuple[]

    for iteration in 1:max_iterations
        push!(history, (; iteration, energy, gradient_norm=norm(gradient)))
        norm(gradient) * diameter_scale < 1e-5 && break

        # L-BFGS two-loop recursion: estimate inverse-curvature times gradient
        # from recent accepted steps, without forming a large Hessian matrix.
        working_gradient = copy(gradient)
        alpha = zeros(length(position_changes))
        for index in reverse(eachindex(position_changes))
            alpha[index] = dot(position_changes[index], working_gradient) /
                dot(gradient_changes[index], position_changes[index])
            working_gradient .-= alpha[index] .* gradient_changes[index]
        end
        inverse_curvature_scale = isempty(position_changes) ? 1e-7 * diameter_scale^2 :
            dot(position_changes[end], gradient_changes[end]) / dot(gradient_changes[end], gradient_changes[end])
        direction = inverse_curvature_scale .* working_gradient
        for index in eachindex(position_changes)
            beta = dot(gradient_changes[index], direction) /
                dot(gradient_changes[index], position_changes[index])
            direction .+= position_changes[index] .* (alpha[index] - beta)
        end
        direction .*= -1
        if !(dot(direction, gradient) < 0)
            direction = -gradient * (1e-7 * diameter_scale^2)
        end

        maximum_motion = maximum(sqrt.(sum(reshape(direction, :, 2).^2, dims=2)))
        step_limit = 0.0002 * diameter_scale
        if maximum_motion > step_limit
            direction .*= step_limit / maximum_motion
        end

        accepted = false
        for backtrack in 0:22
            trial = variables + 2.0^(-backtrack) * direction
            trial_energy, trial_gradient = layout_energy(mesh, selection, trial)
            if isfinite(trial_energy) && trial_energy <= energy + 1e-4 * dot(gradient, trial - variables)
                position_change = trial - variables
                gradient_change = trial_gradient - gradient
                # Retain only pairs with positive, sufficiently resolved curvature.
                if dot(position_change, gradient_change) > 1e-12 * norm(position_change) * norm(gradient_change)
                    push!(position_changes, position_change)
                    push!(gradient_changes, gradient_change)
                    if length(position_changes) > 8
                        popfirst!(position_changes)
                        popfirst!(gradient_changes)
                    end
                end
                variables, energy, gradient = trial, trial_energy, trial_gradient
                accepted = true
                break
            end
        end
        accepted || break
    end

    nodes = copy(selection.nodes)
    nodes[selection.free, :] = reshape(variables, :, 2)
    optimized_mesh = assemble_surface_mesh(
        [Tuple(nodes[node, :]) for node in axes(nodes, 1)], mesh.quads, selection.rim;
        R=mesh.R, L=mesh.L, ca=mesh.ca, slope=mesh.slope,
        feed_z_m=mesh.feed_z_m, inlet_radius_m=mesh.inlet_radius_m,
    )
    return (;
        mesh=optimized_mesh,
        history,
        initial_energy=first(history).energy,
        final_energy=energy,
        free_nodes=length(selection.free),
        maximum_move_m=maximum(sqrt.(sum((nodes - selection.nodes).^2, dims=2))),
    )
end
