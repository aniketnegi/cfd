# Entry point for cylindrical/conical fitted inlet construction.
# Internal coordinates and lengths use metres; surface_feed_grid converts to mm.
module ControlledSurfaceInlet

using LinearAlgebra, SparseArrays
import ..BearingParameters, ..bearing_geometry, ..housing_radius
import ..background_grid, ..patch_axis

using ..SurfaceMeshQuality: cross2, cell_chart_area, cell_chart_centroid,
    chart_axial_position_mm, assemble_surface_mesh, passes_quality_limits,
    mesh_acceptance, local_edge_lengths, accumulate_motion!
const DEFAULT_MAX_ITERATIONS = 5000
const DEFAULT_MOVEMENT_TOLERANCE = 0.001
const DEFAULT_STABILITY_STEPS = 20

include("seed.jl")
include("placement.jl")

function surface(p::BearingParameters, n_theta::Int, n_axial::Int;
    feed_diameter::Real, feed_z::Real, max_iterations::Int=DEFAULT_MAX_ITERATIONS,
    movement_tolerance::Real=DEFAULT_MOVEMENT_TOLERANCE,
    stability_steps::Int=DEFAULT_STABILITY_STEPS)
    stability_steps > 0 || throw(ArgumentError("public meshes require a positive stability window"))
    initial_mesh = reference_seed(p, n_theta, n_axial; feed_diameter, feed_z)
    placement = controlled_layout(initial_mesh; max_iterations, movement_tolerance, stability_steps)
    mesh = placement.mesh

    quality = mesh.quality
    passes_quality_limits(quality) ||
        error("inlet mesh $n_theta×$n_axial does not meet the surface-quality limits: $quality")
    placement.acceptance.accepted ||
        error("inlet mesh $n_theta×$n_axial not accepted: $(placement.acceptance); quality=$quality")

    return (; mesh..., n_theta, n_axial, acceptance=placement.acceptance,
        placement=(; placement.accepted_steps, placement.termination_reason,
            placement.window_motion_ratio, placement.window_steps,
            placement.movement_tolerance, placement.stability_steps),
        patch_cells=initial_mesh.patch_cells,
        transition_layers=initial_mesh.transition_layers)
end

end
