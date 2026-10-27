# Fixed-mesh force/flow references are not verified for the public mesher.
# Excluded from mesh/tests/runtests.jl; loads the experimental operator and Gmsh.
using Test

module CurrentMesh
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__, "../../mesh/2d/surface_feed.jl"))
end

module ReviewedTrial
include(joinpath(@__DIR__, "operator.jl"))
end

p=CurrentMesh.BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
    bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=.6)
coarse=CurrentMesh.surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=50.)

@testset "Inlet mesh with conservative fluid operator" begin
    reference=ReviewedTrial.BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    pose=ReviewedTrial.floating_bush_pose(reference,(-9.152501266813789e-7,1.7844427760798323e-7,0.))
    state=ReviewedTrial.solve_fitted(coarse.surface,pose)
    @test maximum(abs,state.force-[-18.24857748898862,-1002.1463387752476,321.13450985407525])<1e-6
    @test state.net_feed*60000≈.030887303265211042 rtol=1e-10
    @test minimum(state.pressure)>=0
    @test 0<=minimum(state.fill)<=maximum(state.fill)<=1
    @test state.imbalance<1e-10
    println("Coarse held-pose force=$(state.force), net feed=$(state.net_feed*60000) L/min; NOT equilibrium")
end
