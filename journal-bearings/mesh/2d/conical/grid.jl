# Shared periodic circumferential–axial Quad4 grid on the housing cone.
# Lengths are mm, theta is in radians. This grid defines no feed or solver conditions.
include(joinpath(@__DIR__, "..", "..", "..", "bearing.jl"))

function housing_grid(p::BearingParameters, n_theta::Int, n_axial::Int)
    n_theta >= 3 || throw(ArgumentError("n_theta must be at least 3"))
    n_axial >= 1 || throw(ArgumentError("n_axial must be positive"))
    p.eccentricity < p.journal_midpoint_radius - p.half_bearing_length * p.slope ||
        throw(ArgumentError("journal eccentricity exceeds the small-end radius"))

    dtheta = 2pi / n_theta
    dz = p.bearing_length / n_axial
    ds = dz / cosd(p.semicone_angle_deg)
    node_id(i, j) = i + n_theta * (j - 1)
    nodes = Vector{NTuple{3,Float64}}(undef, n_theta * (n_axial + 1))
    quads = Vector{NTuple{4,Int}}(undef, n_theta * n_axial)
    radial_gaps = Vector{Float64}(undef, length(quads))
    areas = Vector{Float64}(undef, length(quads))

    # z=0 is the large end. Unique circumferential nodes wrap at the seam.
    for j in 1:(n_axial + 1)
        z = (j - 1) * dz
        for i in 1:n_theta
            theta = (i - 1) * dtheta
            nodes[node_id(i, j)] = bearing_point(p, theta, z, 1.0)
        end
    end

    for j in 1:n_axial
        z = (j - 0.5) * dz
        for i in 1:n_theta
            next_i = i == n_theta ? 1 : i + 1
            k = i + n_theta * (j - 1)
            quads[k] = (node_id(i, j), node_id(next_i, j),
                        node_id(next_i, j + 1), node_id(i, j + 1))

            theta = (i - 0.5) * dtheta
            xh, yh, _ = bearing_point(p, theta, z, 1.0)
            xj, yj, _ = bearing_point(p, theta, z, 0.0)
            radial_gaps[k] = hypot(xh, yh) - hypot(xj, yj)
            # Analytic cone area, not the area of the faceted output quad.
            areas[k] = hypot(xh, yh) * dtheta * ds
        end
    end

    minimum(radial_gaps) > 0 || error("journal intersects housing")
    return (; nodes, quads, radial_gaps, areas, n_theta, n_axial)
end

function write_vtk(path::AbstractString, grid; feed_cells::Union{Nothing,AbstractVector{Bool}}=nothing)
    isnothing(feed_cells) || length(feed_cells) == length(grid.quads) ||
        throw(ArgumentError("feed_cells must have one value per quad"))
    feed_edges = hasproperty(grid, :feed_edges) ? grid.feed_edges : NTuple{2,Int}[]
    !isempty(feed_edges) && !isnothing(feed_cells) &&
        throw(ArgumentError("a fitted rim is a boundary, not a feed-cell mask"))
    # Keep the unported writer's existing temporary-file contract used by
    # atomic solver exports. New fitted exports never overwrite an existing file.
    !isempty(feed_edges) && ispath(path) && throw(ArgumentError("output already exists: $path"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Conical housing surface (mm)")
        println(io, "ASCII\nDATASET POLYDATA")
        println(io, "POINTS $(length(grid.nodes)) double")
        for (x, y, z) in grid.nodes
            println(io, "$x $y $z")
        end
        if !isempty(feed_edges)
            println(io, "LINES $(length(feed_edges)) $(3 * length(feed_edges))")
            for (a, b) in feed_edges
                println(io, "2 $(a - 1) $(b - 1)")
            end
        end
        println(io, "POLYGONS $(length(grid.quads)) $(5 * length(grid.quads))")
        for (a, b, c, d) in grid.quads
            println(io, "4 $(a - 1) $(b - 1) $(c - 1) $(d - 1)")
        end
        # POLYDATA cell ordering: lines first, then polygons. Rim lines are
        # boundary annotations, not additional surface/pressure cells.
        println(io, "CELL_DATA $(length(feed_edges) + length(grid.quads))")
        for (name, values) in (("radial_gap_mm", grid.radial_gaps),
                               ("analytic_area_mm2", grid.areas))
            println(io, "SCALARS $name double 1\nLOOKUP_TABLE default")
            for _ in feed_edges
                println(io, name == "analytic_area_mm2" ? 0.0 : NaN)
            end
            for value in values
                println(io, value)
            end
        end
        if !isempty(feed_edges)
            println(io, "SCALARS feed_rim int 1\nLOOKUP_TABLE default")
            for _ in feed_edges; println(io, 1); end
            for _ in grid.quads; println(io, 0); end
            println(io, "SCALARS surface_corner_score double 1\nLOOKUP_TABLE default")
            for _ in feed_edges; println(io, NaN); end
            for value in grid.surface.corner_scores; println(io, value); end
        end
        if !isnothing(feed_cells)
            println(io, "SCALARS feed_cell int 1\nLOOKUP_TABLE default")
            for selected in feed_cells
                println(io, Int(selected))
            end
        end
    end
end
