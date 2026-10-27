# Compare requested inlet backgrounds. Outputs are generated evidence, not source.
using LinearAlgebra: norm
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))
include(joinpath(@__DIR__, "../vtk.jl"))

function area_weighted_edge_ratio(grid)
    surface = grid.surface
    ratios = Float64[]
    for coordinates in surface.localnodes
        edge_lengths = Float64[]
        for corner in 1:4
            z = coordinates[corner, 2] * surface.ca
            radius = surface.R + (surface.L / 2 - z) * surface.slope
            edge = coordinates[mod1(corner + 1, 4), :] - coordinates[corner, :]
            push!(edge_lengths, norm(edge .* [radius / surface.R, 1.0]))
        end
        push!(ratios, maximum(edge_lengths) / minimum(edge_lengths))
    end
    return sum(ratios .* grid.areas) / sum(grid.areas)
end

function compare_counts(output_directory)
    ispath(output_directory) && error("output directory already exists: $output_directory")
    mkpath(output_directory)
    p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
        radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
    counts = ((512, 256), (640, 208), (512, 160), (768, 256), (800, 256), (513, 257))
    columns = ("n_theta", "n_axial", "cells", "vertices", "rim_edges", "patch_cells",
        "transition_layers", "minimum_corner", "minimum_face_alignment", "maximum_area_ratio",
        "maximum_edge_ratio", "area_weighted_edge_ratio", "generation_seconds", "status", "error")

    open(joinpath(output_directory, "comparison.tsv"), "w") do table
        println(table, join(columns, '\t'))
        for (n_theta, n_axial) in counts
            grid = nothing
            seconds = 0.0
            try
                seconds = @elapsed grid = surface_feed_grid(p, n_theta, n_axial;
                    feed_diameter=5., feed_z=50.)
            catch err
                message = sprint(showerror, err)
                # Preserve rejected cases, but do not hide programming errors.
                err isa ErrorException && occursin("does not meet the surface-quality limits", message) || rethrow()
                println(table, join((n_theta, n_axial, fill(NaN, 11)...,
                    "quality rejected", replace(message, '\n' => ' ', '\t' => ' ')), '\t'))
                flush(table)
                println("$n_theta×$n_axial rejected: $message")
                continue
            end

            quality = grid.quality
            mean_ratio = area_weighted_edge_ratio(grid)
            values = (n_theta, n_axial, length(grid.quads), length(grid.nodes),
                length(grid.feed_edges), grid.surface.patch_cells, grid.surface.transition_layers,
                quality.minimum_scaled_corner, quality.minimum_face_alignment,
                quality.maximum_neighbour_area_ratio, quality.maximum_edge_ratio,
                mean_ratio, seconds, "accepted", "")
            println(table, join(values, '\t'))
            flush(table)
            write_vtk(joinpath(output_directory, "inlet-$n_theta-$n_axial.vtk"), grid)
            println("$n_theta×$n_axial: $(length(grid.quads)) cells; mean edge ratio $mean_ratio; $seconds s")
            flush(stdout)
        end
    end
    println("Wrote comparison.tsv and accepted meshes to $(abspath(output_directory))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("usage: compare_inlet_counts.jl [FRESH_OUTPUT_DIRECTORY]")
    directory = isempty(ARGS) ? joinpath(@__DIR__, "../../out/cases/conical-inlet-count-comparison") : only(ARGS)
    compare_counts(directory)
end
