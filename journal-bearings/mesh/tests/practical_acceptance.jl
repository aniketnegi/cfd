using Test, LinearAlgebra
include(joinpath(@__DIR__,"../examples/inlet_long_run.jl"))
const Q=SurfaceMeshQuality

@testset "Acceptance is geometry AND quality AND placement" begin
    good=(minimum_scaled_corner=.95,minimum_face_alignment=.96,maximum_neighbour_area_ratio=1.4)
    for valid in (false,true), stable in (false,true), quality in (good,merge(good,(minimum_scaled_corner=.9,)))
        a=Q.mesh_acceptance(quality;valid_geometry=valid,placement_stable=stable,termination_reason=:test)
        @test a.accepted==(valid && stable && quality.minimum_scaled_corner>.9)
    end
    @test !Q.passes_quality_limits(merge(good,(maximum_neighbour_area_ratio=2.,)))
    @test !Q.passes_quality_limits(merge(good,(minimum_face_alignment=NaN,)))
end

@testset "Physical distances, local edges and no cancellation" begin
    for angle in (0.,10.,25.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        mesh=ControlledSurfaceInlet.reference_seed(p,128,64;feed_diameter=5.,feed_z=50.)
        for (a,b) in (((.1,.04),(.101,.041)),((mesh.period-1e-4,.03),(1e-4,.031)),((.08,.05),(.08+1e-12,.05-2e-12)))
            # Independent Cartesian reference with actual BigFloat intermediates.
            reference=setprecision(BigFloat,256) do
                function point(v)
                    z=BigFloat(v[2])*BigFloat(mesh.ca)
                    r=(BigFloat(mesh.geometry_mm.journal_midpoint_radius)+
                       (BigFloat(mesh.geometry_mm.half_bearing_length)-1000z)*BigFloat(mesh.slope)+
                       BigFloat(mesh.geometry_mm.radial_clearance))/1000
                    theta=BigFloat(v[1])/BigFloat(mesh.R)
                    return [r*cos(theta),r*sin(theta),z]
                end
                Float64(norm(point(a)-point(b)))
            end
            @test Q.physical_distance(mesh,a,b) ≈ reference rtol=2e-11
        end
        h=Q.local_edge_lengths(mesh)
        @test all(i->isfinite(h[i])&&h[i]>0,ControlledSurfaceInlet.macro_data(mesh).free)
        initial=[mesh.R*pi/2,.05]; moved=initial+[1e-5,0.]
        path=zeros(1);scale=[.001]
        Q.accumulate_motion!(path,mesh,initial,moved,scale)
        Q.accumulate_motion!(path,mesh,moved,initial,scale)
        @test path[1] ≈ 2 * Q.physical_distance(mesh,initial,moved)/scale[1] rtol=1e-14
        @test path[1]>0 # Net displacement is zero; actual path is not.
    end
end

@testset "Unported analytic metrics equal general surface metrics" begin
    for angle in (0.,10.,25.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        for (theta,z) in ((collect(0:7).*2pi/8,[0.,25.,50.,75.,100.]),
            ([0.,1.2,2.8,4.6],[0.,25.,60.,100.]))
            grid=housing_grid(p,length(theta),length(z)-1;theta_nodes=theta,z_nodes=z)
            R=p.housing_midpoint_radius*1e-3;ca=cosd(angle)
            chart=[(R*t,axial*1e-3/ca) for axial in z for t in theta]
            full=Q.assemble_surface_mesh(chart,grid.quads,Set{Tuple{Int,Int}}();R,L=.1,ca,
                slope=p.slope,geometry_mm=bearing_geometry(p),feed_z_m=.05,inlet_radius_m=0.)
            for key in (:minimum_scaled_corner,:minimum_face_alignment,:maximum_neighbour_area_ratio,:maximum_edge_ratio)
                @test getproperty(grid.quality,key) ≈ getproperty(full.quality,key) rtol=3e-12
            end
            @test grid.acceptance.accepted && grid.acceptance.termination_reason==:direct_construction
            @test !haskey(grid,:placement) # No fake optimizer run for a tensor grid.
        end
        @test_throws ErrorException housing_grid(p,8,3;z_nodes=[0.,1.,50.,100.])
        @test_throws ArgumentError housing_grid(p,4,3;theta_nodes=[0.,.4,2.1,5.6])
    end
end

@testset "Wide inlet support cannot move the periodic seam" begin
    # The 4D movement halo reaches the seam; seam vertices must stay fixed.
    p=BearingParameters(journal_midpoint_radius=20.,bearing_length=70.,
        radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=.6)
    seed=ControlledSurfaceInlet.reference_seed(p,128,64;feed_diameter=8.,feed_z=35.)
    selection=ControlledSurfaceInlet.macro_data(seed)
    seam=findall(point->point[1]==0,seed.nodes)
    @test length(seam)==65
    @test 4*.008>seed.R*pi/2
    @test isempty(intersect(seam,selection.free))
    @test isempty(intersect(selection.boundary,selection.free))
    result=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=2)
    @test result.mesh.nodes[seam]==seed.nodes[seam]
    @test result.accepted_steps==2
end

@testset "Practical fitted stop and independent last-window measurement" begin
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
        semicone_angle_deg=10.,eccentricity_ratio=.6)
    seed=ControlledSurfaceInlet.reference_seed(p,128,64;feed_diameter=5.,feed_z=50.)
    selection=ControlledSurfaceInlet.macro_data(seed)
    recent=NamedTuple[]
    final_observation=Ref{Any}(nothing)
    function observe(s)
        s.phase==:final && (final_observation[]=s)
        push!(recent,(;iteration=s.iteration,variables=s.variables))
        length(recent)>21 && popfirst!(recent)
    end
    plain=ControlledSurfaceInlet.controlled_layout(seed)
    watched=ControlledSurfaceInlet.controlled_layout(seed;observer=observe,observe_every=1)
    @test plain==watched
    @test plain.acceptance.accepted
    @test plain.termination_reason==:placement_stable
    @test plain.window_steps==20
    @test final_observation[].acceptance==plain.acceptance
    @test final_observation[].window_motion_ratio==plain.window_motion_ratio
    @test final_observation[].placement_stable
    @test plain.window_motion_ratio<=.001
    @test !plain.gradient_tolerance_met # No compulsory old-gradient criterion.
    public_grid=surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=50.)
    @test public_grid.acceptance==plain.acceptance
    @test public_grid.surface.nodes==plain.mesh.nodes
    @test public_grid.placement.window_steps==20
    @test [r.iteration for r in recent]==collect(plain.accepted_steps-20:plain.accepted_steps)
    # Physical edge lengths directly from exported Cartesian bearing points.
    initial=physical_study_nodes(p,seed);h=fill(Inf,length(initial))
    for (a,b) in seed.edge_nodes
        d=norm(collect(initial[a]).-collect(initial[b]))*1e-3
        h[a]=min(h[a],d);h[b]=min(h[b],d)
    end
    path=zeros(length(selection.free))
    previous=free_cartesian(p,seed,selection,first(recent).variables)
    for state in recent[2:end]
        points=free_cartesian(p,seed,selection,state.variables)
        path .+= node_distances(previous,points).*1e-3./h[selection.free]
        previous=points
    end
    @test maximum(path) ≈ plain.window_motion_ratio rtol=2e-8
    @test ControlledSurfaceInlet.controlled_layout(seed;max_iterations=19).acceptance.accepted==false
    @test_throws ErrorException surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=50.,max_iterations=19)
    @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=50.,stability_steps=0)
    badseed=ControlledSurfaceInlet.reference_seed(p,64,32;feed_diameter=5.,feed_z=50.)
    bad=ControlledSurfaceInlet.controlled_layout(badseed;stability_steps=1,movement_tolerance=1e6)
    @test bad.placement_stable && !bad.acceptance.accepted
    @test bad.termination_reason==:quality_failure
    for kwargs in ((;stability_steps=-1),(;movement_tolerance=NaN),(;movement_tolerance=0.),(;max_iterations=0))
        @test_throws ArgumentError ControlledSurfaceInlet.controlled_layout(seed;kwargs...)
    end
end
