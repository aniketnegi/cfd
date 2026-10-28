# Cylindrical entry point and fixed-input example for signed full-film Reynolds.
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__, "full_film.jl"))
include(joinpath(@__DIR__, "pressure_vtk.jl"))

const CYLINDER_CASE_DIR = joinpath(@__DIR__, "..", "..", "out", "cases",
    "cylindrical-unported-e0p6")
const CYLINDER_RUN_DIR = joinpath(CYLINDER_CASE_DIR, "runs", "2000rpm-vg32-40c")

"""Signed full-film Reynolds on an unported uniform cylinder, with zero end
pressure. For supplied/adjusted meshes, use solve_full_film(grid; ...).
"""
function solve_cylindrical_full_film(p::BearingParameters, ntheta::Int, nz::Int;
        rpm::Real, viscosity_pa_s::Real)
    iszero(p.semicone_angle_deg) ||
        throw(ArgumentError("this entry point is cylindrical only"))
    return solve_full_film(housing_grid(p,ntheta,nz); rpm, viscosity_pa_s)
end

function main()
    # Illustrative nominal VG32 viscosity and separately assumed density.
    p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
        radial_clearance=.05, semicone_angle_deg=0., eccentricity_ratio=.6)
    result = solve_cylindrical_full_film(p,512,256;
        rpm=2000., viscosity_pa_s=32e-6*860.)
    mesh_path = joinpath(CYLINDER_CASE_DIR,"mesh.vtk")
    path = joinpath(CYLINDER_RUN_DIR,"full-film.vtk")
    write_vtk(mesh_path,result.grid)
    write_pressure_vtk(path,result)
    println("Wrote signed full-film gauge pressure to $(abspath(path))")
    println("Pressure: $(minimum(result.pressure_pa)) to $(maximum(result.pressure_pa)) Pa")
    println("Conservation diagnostics: $(result.diagnostics)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
