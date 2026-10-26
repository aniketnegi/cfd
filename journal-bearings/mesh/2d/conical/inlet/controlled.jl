"""Controlled surface inlet for the reviewed reference geometry only.

Internal chart lengths are metres. The public surface_feed_grid API converts
to mm. There are no dependencies on experimental scripts, Gmsh or fluid solvers.
"""
module ControlledSurfaceInlet
using LinearAlgebra, SparseArrays
include("metrics.jl")
include("seed.jl")
include("placement.jl")

function surface(level::Int)
    level in (1,2,4) || throw(ArgumentError("supported inlet refinement levels: 1,2,4"))
    placement=controlled_layout(reference_seed())
    mesh=subdivide_layout(placement.mesh,level)
    q=mesh.quality
    # Regression bounds from the accepted surface trial, NOT Hex8 OQ criteria.
    q.minimum_scaled_corner>.9 && q.minimum_face_alignment>.9 &&
        q.maximum_neighbour_area_ratio<2 || error("controlled inlet surface quality regression")
    return mesh
end
end
