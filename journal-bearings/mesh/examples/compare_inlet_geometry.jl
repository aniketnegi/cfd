# Vary geometry without changing the construction settings or quality limits.
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))
include(joinpath(@__DIR__, "../vtk.jl"))

function geometry_cases()
    reference = (radius=50., length=100., clearance=.05, angle=10.,
        diameter=5., feed_z=50., n_theta=128, n_axial=64)
    return (
        "reference" => reference,
        "radius-40" => (; reference..., radius=40.),
        "radius-60" => (; reference..., radius=60.),
        "length-80" => (; reference..., length=80., feed_z=40.),
        "length-120" => (; reference..., length=120., feed_z=60.),
        "length-123p4" => (; reference..., length=123.4, feed_z=61.7),
        "length-123p4-angle-15" => (; reference..., length=123.4, angle=15., feed_z=61.7),
        "clearance-0p1" => (; reference..., clearance=.1),
        "angle-0" => (; reference..., angle=0.),
        "angle-5" => (; reference..., angle=5.),
        "angle-15" => (; reference..., angle=15.),
        "diameter-3" => (; reference..., diameter=3.),
        "diameter-6" => (; reference..., diameter=6.),
        "diameter-8" => (; reference..., diameter=8.),
        "inlet-40" => (; reference..., feed_z=40.),
        "inlet-60" => (; reference..., feed_z=60.),
        "inlet-43p7" => (; reference..., feed_z=43.7),
        "combined-odd" => (radius=45., length=90., clearance=.06, angle=12.,
            diameter=4., feed_z=35.3, n_theta=127, n_axial=65),
        "combined-fine" => (radius=60., length=120., clearance=.08, angle=15.,
            diameter=6., feed_z=70.3, n_theta=640, n_axial=208),
        "opening-touches-end" => (; reference..., feed_z=2.5),
        "transition-crosses-end" => (; reference..., feed_z=5.),
        "oversized-opening" => (; reference..., diameter=100.),
    )
end

function compare_geometry(output_directory)
    ispath(output_directory) && error("output directory already exists: $output_directory")
    mkpath(output_directory)
    columns = ("name", "radius_mm", "length_mm", "clearance_mm", "angle_deg",
        "diameter_mm", "feed_z_mm", "n_theta", "n_axial", "status", "cells",
        "vertices", "rim_edges", "minimum_corner", "minimum_alignment",
        "maximum_chart_area_ratio", "maximum_edge_ratio", "seconds", "reason")

    open(joinpath(output_directory, "comparison.tsv"), "w") do table
        println(table, join(columns, '\t'))
        for (name, geometry) in geometry_cases()
            inputs = (name, geometry.radius, geometry.length, geometry.clearance,
                geometry.angle, geometry.diameter, geometry.feed_z, geometry.n_theta, geometry.n_axial)
            grid = nothing
            seconds = 0.0
            try
                p = BearingParameters(journal_midpoint_radius=geometry.radius,
                    bearing_length=geometry.length, radial_clearance=geometry.clearance,
                    semicone_angle_deg=geometry.angle, eccentricity_ratio=.6)
                seconds = @elapsed grid = surface_feed_grid(p, geometry.n_theta, geometry.n_axial;
                    feed_diameter=geometry.diameter, feed_z=geometry.feed_z)
            catch err
                message = sprint(showerror, err)
                expected_rejection = err isa ArgumentError ||
                    (err isa ErrorException && (occursin("does not meet the surface-quality limits", message) ||
                        occursin("not accepted:", message)))
                expected_rejection || rethrow()
                println(table, join((inputs..., "rejected", fill(NaN, 8)...,
                    replace(message, '\n' => ' ', '\t' => ' ')), '\t'))
                flush(table)
                println("$name rejected: $message")
                continue
            end

            quality = grid.quality
            values = (inputs..., "accepted", length(grid.quads), length(grid.nodes),
                length(grid.feed_edges), quality.minimum_scaled_corner, quality.minimum_face_alignment,
                quality.maximum_neighbour_area_ratio, quality.maximum_edge_ratio, seconds, "")
            println(table, join(values, '\t'))
            flush(table)
            write_vtk(joinpath(output_directory, "$name.vtk"), grid)
            println("$name accepted: $(length(grid.quads)) cells; $quality")
            flush(stdout)
        end
    end
    println("Wrote geometry comparison to $(abspath(output_directory))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("usage: compare_inlet_geometry.jl [FRESH_OUTPUT_DIRECTORY]")
    directory = isempty(ARGS) ? joinpath(@__DIR__, "../../out/cases/conical-inlet-geometry-comparison") : only(ARGS)
    compare_geometry(directory)
end
