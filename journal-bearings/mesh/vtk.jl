function write_vtk(
    path::AbstractString,
    grid;
    feed_cells::Union{Nothing,AbstractVector{Bool}} = nothing
)
    (
        isnothing(feed_cells)
        || length(feed_cells) == length(grid.quads)
        || throw(ArgumentError("feed_cells must have one value per quad"))
    )

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
