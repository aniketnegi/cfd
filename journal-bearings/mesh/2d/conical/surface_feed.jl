# Quadrilateral mesh around the +Y inlet opening.
# Public lengths are mm and areas are mm². Caller loads bearing.jl and grid.jl first.
include(joinpath(@__DIR__, "inlet", "controlled.jl"))

"""
    surface_feed_grid(p, n_theta, n_axial; feed_diameter, feed_z)

Build a conical mesh around a circular +Y drilling at axial position feed_z.
Cells cover the solid surface around the opening; `feed_edges` identifies its rim.

Bearing dimensions come from p; feed_diameter and feed_z use mm.
n_theta and n_axial size the initial background directly. The inlet replaces
part of that background, so the returned cells do not form a rectangular array.
Node optimization runs at the requested resolution; there is no hidden coarse
grid or subsequent subdivision. Insufficient counts or a failed quality check
reject the mesh.

`nodes`, `radial_gaps` and `areas` use mm/mm/mm². `surface` retains the internal
surface coordinates in metres and the cell/face connectivity. Its coordinate
areas differ from physical conical areas. Eccentricity changes the gaps, not
housing vertices. Invalid inputs or unavailable transitions throw `ArgumentError`;
failed quality checks report the measurements and reject the mesh.
"""
function surface_feed_grid(
    p::BearingParameters,
    n_theta::Int,
    n_axial::Int;
    feed_diameter::Real,
    feed_z::Real,
    max_iterations::Int=ControlledSurfaceInlet.DEFAULT_MAX_ITERATIONS,
    movement_tolerance::Real=ControlledSurfaceInlet.DEFAULT_MOVEMENT_TOLERANCE,
    stability_steps::Int=ControlledSurfaceInlet.DEFAULT_STABILITY_STEPS,
)
    surface = ControlledSurfaceInlet.surface(p, n_theta, n_axial;
        feed_diameter, feed_z, max_iterations, movement_tolerance, stability_steps)

    # Internal coordinates are (R_ref * theta, z / cos(angle)), in metres.
    nodes = Vector{NTuple{3,Float64}}(undef, length(surface.nodes))
    for (node, (chart_x, chart_y)) in enumerate(surface.nodes)
        theta = chart_x / surface.R
        z_mm = ControlledSurfaceInlet.chart_axial_position_mm(surface.geometry_mm, surface.ca, chart_y)
        nodes[node] = bearing_point(p, theta, z_mm, 1.0)
    end

    radial_gaps = Vector{Float64}(undef, length(surface.quads))
    areas = Vector{Float64}(undef, length(surface.quads))
    for cell in eachindex(surface.quads)
        theta = surface.centre[cell, 1] / surface.R
        z_mm = ControlledSurfaceInlet.chart_axial_position_mm(
            surface.geometry_mm, surface.ca, surface.centre[cell, 2])
        geometry = bearing_radial_geometry(p, theta, z_mm)
        radial_gaps[cell] = geometry.radial_gap

        # The cone's area scale is linear in z. Evaluating it at the polygon
        # centroid exactly integrates the scale over that chart polygon.
        areas[cell] = surface.areas[cell] * geometry.housing_radius / (surface.R * 1e3) * 1e6
    end
    all(point -> all(isfinite, point), nodes) &&
        all(x -> isfinite(x) && x > 0, radial_gaps) &&
        all(x -> isfinite(x) && x > 0, areas) || error("invalid fitted geometry, gaps or areas")

    feed_faces = findall(==(:feed), surface.boundary)
    feed_edges = surface.edge_nodes[feed_faces]
    return (;
        nodes,
        quads=surface.quads,
        radial_gaps,
        areas,
        n_theta,
        n_axial,
        feed_edges,
        quality=surface.quality,
        acceptance=surface.acceptance,
        placement=surface.placement,
        surface,
        layout="fitted-land-direct-v1",
    )
end
