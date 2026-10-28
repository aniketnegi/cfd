# Count-based cylindrical/conical entry points for signed full-film Reynolds.
# Surface-flux formulation and approximations: conical-formulation.md.
include(joinpath(@__DIR__, "cylindrical.jl"))

conical_radii(p::BearingParameters,z_m::Real) = film_radii(p,z_m)
conical_normal_gap(p::BearingParameters,theta::Real,z_m::Real) =
    normal_film_gap(p,theta,z_m)

function assemble_conical_full_film(p::BearingParameters, nt::Int, nz::Int; kwargs...)
    return assemble_full_film(housing_grid(p,nt,nz); kwargs...)
end

function solve_conical_full_film(p::BearingParameters, nt::Int, nz::Int; kwargs...)
    return solve_full_film(housing_grid(p,nt,nz); kwargs...)
end

"""Pressure force on the housing sidewall for a uniform periodic tensor grid,
in N. End faces and shear are excluded.
"""
function conical_housing_pressure_force(p::BearingParameters, pressure, nt, nz)
    nt >= 4 && nz >= 2 || throw(ArgumentError("grid must be at least 4×2"))
    theta_edges = [(i-1)*2pi/nt for i in 1:nt+1]
    z_edges = [(j-1)*p.bearing_length*1e-3/nz for j in 1:nz+1]
    return tensor_pressure_force(p,pressure,theta_edges,z_edges)
end
