# Count-based cylindrical API for signed full-film Reynolds.
include(joinpath(@__DIR__,"../../../bearing.jl"))
include(joinpath(@__DIR__,"../../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"full_film.jl"))
include(joinpath(@__DIR__,"pressure_vtk.jl"))

"""Signed full-film Reynolds on an unported uniform cylinder, with zero end
pressure. For supplied/adjusted meshes, use solve_full_film(grid; ...).
"""
function solve_cylindrical_full_film(p::BearingParameters, ntheta::Int, nz::Int;
        rpm::Real, viscosity_pa_s::Real)
    iszero(p.semicone_angle_deg) ||
        throw(ArgumentError("this entry point is cylindrical only"))
    return solve_full_film(housing_grid(p,ntheta,nz); rpm, viscosity_pa_s)
end
