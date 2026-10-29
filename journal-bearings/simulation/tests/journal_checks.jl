using Test,LinearAlgebra,SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../elasticity/src/housing.jl"))
include(joinpath(@__DIR__,"../elasticity/src/journal.jl"))

function check_journal()
@testset "Solid journal volume, pressure and outer interpolation" begin
    for angle in (0.,10.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=.6)
        g=housing_grid(p,16,8)
        a=assemble_journal(g;youngs_modulus_pa=210e9,poisson_ratio=.3,n_theta=16,n_axial=8,n_radial=2)
        r0,r1=journal_radius(p,0.)*1e-3,journal_radius(p,100.)*1e-3
        exact_volume=pi*.1*(r0^2+r0*r1+r1^2)/3
        @test a.minimum_jacobian>0
        @test a.volume_m3≈exact_volume*sin(2pi/16)/(2pi/16) rtol=1e-12
        load=a.pressure_operator*fill(1e6,length(g.quads))
        force=vec(sum(reshape(load,3,:),dims=2))
        @test force≈[0.,0.,-2pi*.05*.1*tand(angle)*1e6] atol=1e-8 rtol=1e-11
        @test sum(a.interface_areas_m2)≈2pi*.05*.1/cosd(angle) rtol=1e-12
        pressures=[1e6*(1+.4*cos((i-.5)*2pi/16)+.2*sin((i-.5)*2pi/16)+.3*(j-.5)/8) for j in 1:8 for i in 1:16]
        expected=zeros(3)
        for j in 1:8,i in 1:16
            lo,hi=(i-1)*2pi/16,i*2pi/16
            z0,z1=(j-1)*.1/8,j*.1/8
            area=(z1-z0)*journal_radius(p,(z0+z1)*500)*1e-3
            expected.-=pressures[i+16*(j-1)].*area.*[sin(hi)-sin(lo),cos(lo)-cos(hi),tand(angle)*(hi-lo)]
        end
        nonuniform=vec(sum(reshape(a.pressure_operator*pressures,3,:),dims=2))
        @test nonuniform≈expected rtol=1e-7 atol=1e-7
        state=solve_elasticity(a,load)
        @test state.relative_residual<1e-9
        @test norm(vec(sum(reshape(load+state.reaction_n,3,:),dims=2)))<1e-7
        points=[(theta,z) for theta in (0.,.1,pi,2pi) for z in (0.,.037,.1)]
        G=journal_gap_operator(a,points)
        translation=repeat([1e-7,2e-7,3e-7],size(a.nodes,1))
        expected=[cosd(angle)*(cos(t)*1e-7+sin(t)*2e-7)+sind(angle)*3e-7 for (t,z) in points]
        @test G*translation≈expected rtol=1e-12 atol=1e-20
        for (t,z) in points
            ids,N=journal_interface_weights(a.mesh,t,z)
            @test sum(N)≈1 rtol=1e-15
            @test all(x->x>=-1e-14,N)
            @test dot(collect(N),a.nodes[collect(ids),3])≈z atol=1e-15
        end
    end
end
@testset "Land plus reservoir trace on the complete journal" begin
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,semicone_angle_deg=10.,eccentricity_ratio=.6)
    g=surface_feed_grid(p,128,64;feed_diameter=7.,feed_z=50.)
    a=assemble_journal(g;youngs_modulus_pa=210e9,poisson_ratio=.3,n_theta=32,n_axial=16,n_radial=4)
    @test a.reservoir_area_m2>0
    @test sum(a.interface_areas_m2)+a.reservoir_area_m2≈2pi*.05*.1/cosd(10.) rtol=1e-12
    pressure=1e6
    load=pressure.*(a.pressure_operator*ones(length(g.quads))+a.reservoir_load_per_pa)
    force=vec(sum(reshape(load,3,:),dims=2))
    @test force≈[0.,0.,-2pi*.05*.1*tand(10.)*pressure] rtol=1e-10 atol=1e-7
    @test norm(a.reservoir_load_per_pa)>1e-8
end
@testset "Independent solid-cylinder compression reference" begin
    errors=Float64[]
    for nt in (16,32,64)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=0.)
        g=housing_grid(p,nt,4)
        mesh=journal_solid_mesh(p;n_theta=nt,n_axial=4,n_radial=nt÷8)
        raw=assemble_elasticity(mesh.nodes,mesh.cells;youngs_modulus_pa=210e9,poisson_ratio=.3)
        side=nt÷4;centre=1+(side÷2)+(side÷2)*(side+1)
        positive_x=1+side+(side÷2)*(side+1)
        fixed=[collect(3:3:size(raw.K,1));[3centre-2,3centre-1,3positive_x-1]]
        a=constrain_elasticity(raw,fixed)
        W=journal_pressure_operator(mesh,g)
        state=solve_elasticity(a,W.matrix*fill(1e6,length(g.quads)))
        coefficient=-(1+.3)*(1-2*.3)*1e6/210e9
        exact=vec(permutedims(hcat(coefficient.*a.nodes[:,1],coefficient.*a.nodes[:,2],zeros(size(a.nodes,1)))))
        push!(errors,norm(state.displacement_m-exact)/norm(exact))
    end
    println("Solid compressed-cylinder errors: $errors")
    @test all(errors[2:end] .< .4 .* errors[1:end-1])
    @test errors[end]<.002
end
end
if abspath(PROGRAM_FILE)==@__FILE__;check_journal();end
