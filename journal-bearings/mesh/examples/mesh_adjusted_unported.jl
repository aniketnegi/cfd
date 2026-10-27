# Use the inlet's initial spacing on a complete cone, without making an opening.
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../vtk.jl"))

function main()
    p = BearingParameters(
        journal_midpoint_radius=50.0,
        bearing_length=100.0,
        radial_clearance=0.050,
        semicone_angle_deg=10.0,
        eccentricity_ratio=0.6,
    )
    n_theta = 128
    n_axial = 64
    patch_cells = 8
    patch_half_width = 7.5 # mm, matching the current 5 mm inlet construction

    theta_nodes = patch_axis(n_theta, n_theta ÷ 4, 2pi, pi/2,
        patch_half_width / p.housing_midpoint_radius, patch_cells)[1:end-1]
    z_nodes = patch_axis(n_axial, n_axial ÷ 2, p.bearing_length,
        p.half_bearing_length, patch_half_width, patch_cells)
    grid = housing_grid(p, n_theta, n_axial; theta_nodes, z_nodes)

    path = joinpath(@__DIR__, "../../out/cases/conical-unported-adjusted-e0p6/mesh.vtk")
    ispath(path) && throw(ArgumentError("output already exists: $path"))
    write_vtk(path, grid)
    println("Wrote $(length(grid.quads)) quads to $(abspath(path))")
    println("Cell-centre radial gap: $(extrema(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
