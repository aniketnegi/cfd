# Cell-labelled surface feed on the periodic conical housing grid.
# This is not a resolved passage, a cut-cell mesh, or a supply boundary condition.
include(joinpath(@__DIR__, "grid.jl"))

"""
    surface_feed_grid(p, n_theta, n_axial; feed_diameter, feed_z)

Mark housing cells whose centres lie within the +Y intersection of a cylinder
of diameter `feed_diameter`, drilled along +Y at axial position `feed_z`.
All lengths are mm. `feed_cells[k]` corresponds to `quads[k]`; the cell-centre
mask approximates the curved mouth with stair-stepped cells, not a fitted rim.
"""
function surface_feed_grid(p::BearingParameters, n_theta::Int, n_axial::Int;
                           feed_diameter::Real, feed_z::Real)
    diameter, z0 = Float64(feed_diameter), Float64(feed_z)
    isfinite(diameter) && diameter > 0 ||
        throw(ArgumentError("feed_diameter must be positive and finite"))
    isfinite(z0) || throw(ArgumentError("feed_z must be finite"))
    radius = diameter / 2
    radius < z0 < p.bearing_length - radius ||
        throw(ArgumentError("feed footprint must lie strictly between axial ends"))
    smallest_local_radius = p.housing_midpoint_radius +
                            (p.half_bearing_length - z0 - radius) * p.slope
    radius < smallest_local_radius ||
        throw(ArgumentError("feed circle must fit on the +Y housing surface"))

    grid = housing_grid(p, n_theta, n_axial)
    feed_cells = Vector{Bool}(undef, length(grid.quads))
    for j in 1:n_axial, i in 1:n_theta
        theta = (i - 0.5) * (2pi / n_theta)
        z = (j - 0.5) * (p.bearing_length / n_axial)
        x, y, _ = bearing_point(p, theta, z, 1.0)
        # Cylinder cross-section is a circle in the XZ plane. Select only the
        # intended +Y opening, not its opposite-side intersection.
        feed_cells[i + n_theta * (j - 1)] = y > 0 && x^2 + (z - z0)^2 <= radius^2
    end
    any(feed_cells) || throw(ArgumentError("feed footprint has no cell centres; refine the grid"))
    return (; grid..., feed_cells)
end

function main()
    # Gangrade Table 1 geometry, illustrative eccentricity; the 5 mm drill
    # is this project's choice, not a dimension reported by Gangrade.
    p = BearingParameters(
        journal_midpoint_radius = 50.0,
        bearing_length = 100.0,
        radial_clearance = 0.050,
        semicone_angle_deg = 10.0,
        eccentricity_ratio = 0.6,
    )
    grid = surface_feed_grid(p, 512, 256; feed_diameter=5.0, feed_z=50.0)
    path = joinpath(@__DIR__, "..", "..", "..", "out", "conical-surface-feed.vtk")
    write_vtk(path, grid; feed_cells=grid.feed_cells)
    println("Wrote $(length(grid.quads)) quads to $(abspath(path))")
    println("Feed cells: $(count(grid.feed_cells)); selected analytic area: $(sum(grid.areas[grid.feed_cells])) mm²")
    println("Cell-centre radial gap: $(minimum(grid.radial_gaps)) to $(maximum(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
