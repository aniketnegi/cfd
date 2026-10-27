# Entry point for the conical inlet construction.
# Internal coordinates and lengths use metres; surface_feed_grid converts to mm.
module ControlledSurfaceInlet

using LinearAlgebra, SparseArrays
import ..BearingParameters, ..background_grid, ..patch_axis

include("metrics.jl")
include("seed.jl")
include("placement.jl")

function surface(p::BearingParameters, n_theta::Int, n_axial::Int;
    feed_diameter::Real, feed_z::Real)
    initial_mesh = reference_seed(p, n_theta, n_axial; feed_diameter, feed_z)
    placement = controlled_layout(initial_mesh)
    mesh = placement.mesh

    # These are the current surface-quality limits, not guarantees for 3D cells.
    quality = mesh.quality
    quality.minimum_scaled_corner > 0.9 &&
        quality.minimum_face_alignment > 0.9 &&
        quality.maximum_neighbour_area_ratio < 2 ||
        error("inlet mesh $n_theta×$n_axial does not meet the surface-quality limits: $quality")

    return (; mesh..., n_theta, n_axial,
        patch_cells=initial_mesh.patch_cells,
        transition_layers=initial_mesh.transition_layers)
end

end
