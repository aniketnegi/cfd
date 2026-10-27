# Build cell/face geometry in coordinates (R_ref * theta, z / cos(angle)).
# Lengths here are metres. Physical surface measurements require the cone metric.

cross2(a, b) = a[1] * b[2] - a[2] * b[1]

# Keep the supplied end plane exact when converting the metre chart back to mm.
function chart_axial_position_mm(geometry_mm, ca, chart_y)
    end_y = geometry_mm.bearing_length * 1e-3 / ca
    return chart_y == end_y ? geometry_mm.bearing_length : chart_y * ca * 1e3
end

"""
    cell_chart_area(coordinates)

Signed polygon area in the coordinate chart. Subtract the first vertex before
forming cross products so a small cell does not lose accuracy far from zero.
The input must already be unwrapped across the periodic seam.
"""
function cell_chart_area(coordinates)
    origin_x, origin_y = coordinates[1, 1], coordinates[1, 2]
    corners = size(coordinates, 1)
    twice_area = sum(1:corners) do corner
        next_corner = mod1(corner + 1, corners)
        cross2(
            (coordinates[corner, 1] - origin_x, coordinates[corner, 2] - origin_y),
            (coordinates[next_corner, 1] - origin_x, coordinates[next_corner, 2] - origin_y),
        )
    end
    return twice_area / 2
end

"""
    cell_chart_centroid(coordinates, area=cell_chart_area(coordinates))

Polygon centroid in the same unwrapped chart. Form its first moment near the
first vertex, then restore the origin; area is signed and must be nonzero.
"""
function cell_chart_centroid(coordinates, area=cell_chart_area(coordinates))
    iszero(area) && throw(ArgumentError("a polygon centroid requires nonzero chart area"))
    origin = coordinates[1, :]
    relative = coordinates .- transpose(origin)
    corners = size(coordinates, 1)
    moment = sum(1:corners) do corner
        next_corner = mod1(corner + 1, corners)
        (relative[corner, :] + relative[next_corner, :]) *
            cross2(relative[corner, :], relative[next_corner, :])
    end
    return origin + moment / (6 * area)
end

"""
    assemble_surface_mesh(nodes, quads, rim; ...)

Build cell centres, coordinate areas, shared faces and surface-quality measures.
`rim` contains unordered vertex-ID pairs identifying the inlet boundary.

Faces have an owner cell and either a neighbour cell or neighbour=0 at a boundary.
Face normals point out of the owner in the coordinate chart. R/L/feed_z_m and
inlet_radius_m are metres, ca is cos(cone angle), and slope is tan(cone angle).
geometry_mm retains the original bearing dimensions for shared radius queries.
Returned field names are also consumed by the mesh exporter and fitted-grid
calculations.
"""
function assemble_surface_mesh(
    nodes,
    quads,
    rim;
    R,
    L,
    ca,
    slope,
    geometry_mm,
    feed_z_m,
    inlet_radius_m,
    period=2pi * R,
    regular=false,
)
    all(point -> all(isfinite, point) && 0 <= point[1] < period &&
        0 <= point[2] <= L / ca, nodes) || error("invalid surface node coordinates")
    centres = zeros(length(quads), 2)
    chart_areas = zeros(length(quads))
    corner_values = Float64[]
    edge_ratios = Float64[]
    local_nodes = Matrix{Float64}[]

    for (cell, quad) in enumerate(quads)
        coordinates = hcat([collect(nodes[node]) for node in quad]...)'
        coordinates = Matrix(coordinates)

        # A seam-crossing cell must be measured across its short angular span,
        # not across nearly the full circumference of the coordinate chart.
        if maximum(coordinates[:, 1]) - minimum(coordinates[:, 1]) > period / 2
            coordinates[coordinates[:, 1] .< period / 2, 1] .+= period
        end
        push!(local_nodes, coordinates)

        area = cell_chart_area(coordinates)
        area > 0 || error("nonpositive cell $cell")
        centres[cell, :] = cell_chart_centroid(coordinates, area)
        centres[cell, 1] = mod(centres[cell, 1], period)
        chart_areas[cell] = area

        edge_lengths = Float64[]
        for corner in 1:4
            z_mm = chart_axial_position_mm(geometry_mm, ca, coordinates[corner, 2])
            radius = housing_radius(geometry_mm, z_mm) * 1e-3
            surface_scale = [radius / R, 1.0]
            next_edge = (coordinates[mod1(corner + 1, 4), :] - coordinates[corner, :]) .* surface_scale
            previous_edge = (coordinates[mod1(corner - 1, 4), :] - coordinates[corner, :]) .* surface_scale

            # Positive sine of the corner angle: 1 at a right angle, 0 at collapse.
            push!(corner_values, cross2(next_edge, previous_edge) / (norm(next_edge) * norm(previous_edge)))
            push!(edge_lengths, norm(next_edge))
        end
        push!(edge_ratios, maximum(edge_lengths) / minimum(edge_lengths))
    end

    owner = Int[]
    neighbour = Int[]
    face_midpoints = Vector{Float64}[]
    face_normals = Vector{Float64}[]
    face_lengths = Float64[]
    boundary = Symbol[]
    edge_nodes = Tuple{Int,Int}[]
    face_for_edge = Dict{Tuple{Int,Int},Int}()

    for (cell, quad) in enumerate(quads), corner in 1:4
        edge = minmax(quad[corner], quad[mod1(corner + 1, 4)])
        if haskey(face_for_edge, edge)
            face = face_for_edge[edge]
            neighbour[face] == 0 || error("nonmanifold edge")
            neighbour[face] = cell
        else
            coordinates = local_nodes[cell]
            start = coordinates[corner, :]
            finish = coordinates[mod1(corner + 1, 4), :]
            tangent = finish - start
            edge_length = norm(tangent)
            midpoint = (start + finish) / 2
            midpoint[1] = mod(midpoint[1], period)

            push!(owner, cell)
            push!(neighbour, 0)
            push!(face_midpoints, midpoint)
            push!(face_normals, [tangent[2], -tangent[1]] / edge_length)
            push!(face_lengths, edge_length)
            push!(edge_nodes, edge)
            face_for_edge[edge] = length(owner)
            push!(boundary, :internal)
        end
    end

    wrap_distance(distance) = mod(distance + period / 2, period) - period / 2
    face_alignment = Float64[]
    neighbour_area_ratios = Float64[]
    for face in eachindex(owner)
        owner_cell = owner[face]
        neighbour_cell = neighbour[face]
        midpoint = face_midpoints[face]

        if neighbour_cell == 0
            if edge_nodes[face] in rim
                boundary[face] = :feed
            elseif abs(midpoint[2]) < 1e-12
                boundary[face] = :large
            elseif abs(midpoint[2] - L / ca) < 1e-12
                boundary[face] = :small
            else
                error("unclassified boundary")
            end
        end

        connection = neighbour_cell == 0 ?
            midpoint - centres[owner_cell, :] :
            centres[neighbour_cell, :] - centres[owner_cell, :]
        connection[1] = wrap_distance(connection[1])

        # Directions scale with the metric; normals use its inverse so their
        # dot product still describes the physical surface-face alignment.
        z_mm = chart_axial_position_mm(geometry_mm, ca, midpoint[2])
        radius = housing_radius(geometry_mm, z_mm) * 1e-3
        surface_scale = [radius / R, 1.0]
        normal = face_normals[face] ./ surface_scale
        push!(face_alignment,
            dot(connection .* surface_scale, normal) / (norm(connection .* surface_scale) * norm(normal)))
        if neighbour_cell > 0
            push!(neighbour_area_ratios,
                max(chart_areas[owner_cell], chart_areas[neighbour_cell]) /
                min(chart_areas[owner_cell], chart_areas[neighbour_cell]))
        end
    end

    minimum(corner_values) > 0 && minimum(face_alignment) > 0 ||
        error("invalid corner or face orientation")
    quality = (;
        minimum_scaled_corner=minimum(corner_values),
        maximum_edge_ratio=maximum(edge_ratios),
        minimum_face_alignment=minimum(face_alignment),
        maximum_nonorthogonality_deg=acosd(minimum(face_alignment)),
        maximum_neighbour_area_ratio=maximum(neighbour_area_ratios),
        cells=length(quads),
        rim_edges=count(==(:feed), boundary),
    )
    return (;
        nodes,
        quads,
        edge_nodes,
        centre=centres,
        areas=chart_areas,
        owner,
        neighbour,
        mid=face_midpoints,
        normals=face_normals,
        lengths=face_lengths,
        boundary,
        localnodes=local_nodes,
        R,
        L,
        ca,
        period,
        quality,
        regular,
        inlet_radius_m,
        corner_scores=vec(minimum(reshape(corner_values, 4, :), dims=1)),
        face_alignment,
        neighbour_area_ratios,
        slope,
        feed_z_m,
        geometry_mm,
    )
end
