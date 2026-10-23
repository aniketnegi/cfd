# Unported conical housing grid. This does not specify lubricant replenishment.
include(joinpath(@__DIR__, "grid.jl"))

function main()
    # Gangrade Table 1 geometry; 0.6 is an example eccentricity within the
    # paper's study range, not the measured journal position in Fig. 8.
    p = BearingParameters(
        journal_midpoint_radius = 50.0,
        bearing_length = 100.0,
        radial_clearance = 0.050,
        semicone_angle_deg = 10.0,
        eccentricity_ratio = 0.6,
    )
    grid = housing_grid(p, 512, 256)
    path = joinpath(@__DIR__, "..", "..", "..", "out", "conical-unported.vtk")
    write_vtk(path, grid)
    println("Wrote $(length(grid.quads)) quads to $(abspath(path))")
    println("Cell-centre radial gap: $(minimum(grid.radial_gaps)) to $(maximum(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
