# Fitted, land-only surface inlet. Geometry only: no supply pressure or 3D passage.
include(joinpath(@__DIR__, "grid.jl"))
include(joinpath(@__DIR__, "inlet", "controlled.jl"))

"""
    surface_feed_grid(p, n_theta, n_axial; feed_diameter, feed_z)

Build the accepted controlled transition around the +Y circular-drill mouth.
Public lengths/areas are mm/mm². `quads` contains LAND cells only; `feed_edges`
contains pairs of node IDs on the inlet boundary. There is no `feed_cells` mask
and no rectangular i,j indexing of the returned cells. n_theta/n_axial describe
background resolution, not the actual number or indexing of fitted cells.

Currently supported: Rj=50 mm, L=100 mm, c=0.05 mm, semi-angle=10°, 5 mm mouth
at z=50 mm, and background pairs (128,64), (256,128), (512,256). Eccentricity may
vary within BearingParameters' limits: it changes gaps, not the housing mesh.
Other layouts throw ArgumentError rather than silently using the reference
geometry or falling back to the old staircase. `surface` exposes the SI cone
chart/face connectivity for a geometry-aware conservative operator; regular-grid
two-point stencils must not be applied to it. Surface scores are NOT Hex8 OQ.
"""
function surface_feed_grid(p::BearingParameters, n_theta::Int, n_axial::Int;
                           feed_diameter::Real, feed_z::Real)
    dimensions=(p.journal_midpoint_radius,p.bearing_length,p.radial_clearance,
                p.semicone_angle_deg,Float64(feed_diameter),Float64(feed_z))
    reference=(50.,100.,.05,10.,5.,50.)
    all(isfinite,dimensions) && dimensions==reference ||
        throw(ArgumentError("controlled fitted inlet currently supports Rj=50, L=100, c=0.05 mm, angle=10°, diameter=5 mm, z=50 mm only; other layouts need a checked construction"))
    (n_theta,n_axial) in ((128,64),(256,128),(512,256)) ||
        throw(ArgumentError("controlled fitted inlet supports background pairs 128×64, 256×128, 512×256; fitted cells are not an n_theta×n_axial array"))
    p.eccentricity < p.journal_midpoint_radius-p.half_bearing_length*p.slope ||
        throw(ArgumentError("journal eccentricity exceeds the small-end radius"))
    surface=ControlledSurfaceInlet.surface(n_theta÷128)
    nodes=[begin
        theta=x/surface.R;z=y*surface.ca
        rh=surface.R+(.05-z)*tand(10.)
        (1e3*rh*cos(theta),1e3*rh*sin(theta),1e3*z)
    end for (x,y) in surface.nodes]
    radial_gaps=Float64[];areas=Float64[]
    for k in eachindex(surface.quads)
        theta=surface.centre[k,1]/surface.R;z=surface.centre[k,2]*surface.ca*1e3
        rh=p.housing_midpoint_radius+(p.half_bearing_length-z)*p.slope
        xj,yj,_=bearing_point(p,theta,z,0.)
        push!(radial_gaps,rh-hypot(xj,yj))
        # Cone metric is linear in z: centroid integration is exact for the
        # polygonal chart cell, not an assertion of exact circular-mouth area.
        push!(areas,surface.areas[k]*rh/(surface.R*1e3)*1e6)
    end
    minimum(radial_gaps)>0 || error("journal intersects housing")
    feed_edges=surface.edge_nodes[findall(==(:feed),surface.boundary)]
    return (;nodes,quads=surface.quads,radial_gaps,areas,n_theta,n_axial,
            feed_edges,quality=surface.quality,surface,layout="controlled-fitted-land-v1")
end

function main()
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
        radial_clearance=.050,semicone_angle_deg=10.,eccentricity_ratio=.6)
    # A fresh name preserves the old stepped export. Named-case CLI writes
    # a reusable manifest as well; this standalone example is geometry only.
    path=joinpath(@__DIR__,"../../../out/cases/conical-surface-feed-fitted-e0p6-d5mm/mesh.vtk")
    ispath(path) && throw(ArgumentError("output already exists: $path"))
    grid=surface_feed_grid(p,512,256;feed_diameter=5.,feed_z=50.)
    write_vtk(path,grid)
    println("Wrote $(length(grid.quads)) land quads and $(length(grid.feed_edges)) inlet edges to $(abspath(path))")
    println("Surface quality (NOT Hex8 OQ): $(grid.quality)")
    println("Cell-centre radial gap: $(extrema(grid.radial_gaps)) mm")
end

if abspath(PROGRAM_FILE)==@__FILE__;main();end
