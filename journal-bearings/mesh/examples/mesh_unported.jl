# Unported conical housing grid. This does not specify lubricant replenishment.
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/grid.jl"))
include(joinpath(@__DIR__, "../vtk.jl"))

function main(path=joinpath(@__DIR__, "../../out/cases/conical-unported-e0p6/mesh.vtk"))
    ispath(path) && throw(ArgumentError("output already exists: $path"))
    # Gangrade's geometry
    p = BearingParameters(
        journal_midpoint_radius = 50.0,
        bearing_length = 100.0,
        radial_clearance = 0.050,
        semicone_angle_deg = 10.0,
        eccentricity_ratio = 0.6,
    )

    grid = housing_grid(p, 512, 256)
    write_vtk(path, grid)

    println("Wrote $(length(grid.quads)) quads to $(abspath(path))")
    println("Cell-centre radial gap: $(minimum(grid.radial_gaps)) to $(maximum(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("usage: mesh_unported.jl [FRESH_OUTPUT.vtk]")
    isempty(ARGS) ? main() : main(only(ARGS))
end
