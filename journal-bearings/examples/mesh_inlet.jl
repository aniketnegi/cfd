# Generate the supported 5 mm inlet mesh. Run from journal-bearings/.
include(joinpath(@__DIR__, "../bearing.jl"))
include(joinpath(@__DIR__, "../mesh/2d/conical/surface_feed.jl"))
include(joinpath(@__DIR__, "../mesh/vtk.jl"))

function main()
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
        radial_clearance=.050,semicone_angle_deg=10.,eccentricity_ratio=.6)
    # A fresh name preserves the old stepped export. Named-case CLI writes
    # a reusable manifest as well; this standalone example is geometry only.
    path=joinpath(@__DIR__,"../out/cases/conical-surface-feed-fitted-e0p6-d5mm/mesh.vtk")
    ispath(path) && throw(ArgumentError("output already exists: $path"))
    grid=surface_feed_grid(p,512,256;feed_diameter=5.,feed_z=50.)
    write_vtk(path,grid)
    println("Wrote $(length(grid.quads)) land quads and $(length(grid.feed_edges)) inlet edges to $(abspath(path))")
    println("Surface quality (NOT Hex8 OQ): $(grid.quality)")
    println("Cell-centre radial gap: $(extrema(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE)==@__FILE__;main();end
