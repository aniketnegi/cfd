using Test
include(joinpath(@__DIR__,"../examples/inlet_long_run.jl"))

function telemetry_test_seed()
    p= BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
        semicone_angle_deg=10.,eccentricity_ratio=.6)
    return ControlledSurfaceInlet.reference_seed(p,64,32;feed_diameter=5.,feed_z=50.)
end

@testset "Observer cadence, chunks and numerical independence" begin
    seed=telemetry_test_seed()
    plain=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=7)
    events=NamedTuple[]
    watched=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=7,
        observer=state->push!(events,state),observe_every=2)
    @test plain==watched
    @test [e.phase for e in events]==[:initial,:update,:update,:update,:final]
    @test [e.iteration for e in events]==[0,2,4,6,7]
    @test reduce(vcat,[e.history_chunk for e in events])==plain.history
    @test last(events).energy===plain.final_energy
    @test last(events).gradient_norm===plain.final_gradient_norm
    @test last(events).evaluations==plain.evaluations
    @test last(events).termination_reason==:iteration_limit
    @test last(events).acceptance==plain.acceptance
    @test last(events).window_motion_ratio==plain.window_motion_ratio
    @test last(events).window_steps==plain.window_steps
    @test last(events).movement_tolerance==plain.movement_tolerance
    @test first(events).window_steps==0
    @test all(e->e.interval_motion_bound_m>=e.interval_max_step_m,events)
    @test all(e->length(e.position_changes)==length(e.gradient_changes),events)

    function alter_copies(state)
        fill!(state.variables,NaN)
        fill!(state.gradient,Inf)
        foreach(v->fill!(v,NaN),state.position_changes)
        foreach(v->fill!(v,NaN),state.gradient_changes)
        empty!(state.history_chunk)
        if !isnothing(state.diagnostics)
            fill!(state.diagnostics.base_variables,NaN)
            fill!(state.diagnostics.base_gradient,NaN)
            fill!(state.diagnostics.direction,NaN)
            empty!(state.diagnostics.trials)
        end
        return :ignored
    end
    altered=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=7,
        observer=alter_copies,observe_every=1)
    @test altered==plain
    @test_throws ArgumentError ControlledSurfaceInlet.controlled_layout(seed;observe_every=0)
    @test_throws ErrorException ControlledSurfaceInlet.controlled_layout(seed;observer=s->error("observer IO failure"))
end

@testset "Longer termination retains failed-search and original-anchor data" begin
    p= BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
        semicone_angle_deg=10.,eccentricity_ratio=.6)
    seed=ControlledSurfaceInlet.reference_seed(p,128,64;feed_diameter=5.,feed_z=50.)
    final=Ref{Any}(nothing)
    plain=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=720,stability_steps=0)
    watched=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=720,stability_steps=0,
        observer=s->(s.phase==:final && (final[]=s)),observe_every=100)
    @test watched==plain
    state=final[]
    @test state.energy===watched.final_energy
    @test state.gradient_norm===watched.final_gradient_norm
    @test state.accepted_steps==watched.accepted_steps
    @test state.termination_reason==watched.termination_reason
    selection=ControlledSurfaceInlet.macro_data(seed)
    energy,gradient=ControlledSurfaceInlet.layout_energy(seed,selection,state.variables)
    @test energy===watched.final_energy
    @test norm(gradient)===watched.final_gradient_norm
    if watched.termination_reason==:line_search_failure
        @test !state.diagnostics.accepted
        @test length(state.diagnostics.trials)==23
        @test state.diagnostics.base_variables==state.variables
    end
end

@testset "Stationary construction is observed without fictitious work" begin
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=10.,radial_clearance=.05,
        semicone_angle_deg=0.,eccentricity_ratio=0.)
    radius=p.housing_midpoint_radius*1e-3
    nodes=[(radius*theta,z) for z in (0.,.01) for theta in (0.,pi/2,pi,3pi/2)]
    quads=[(i,mod1(i+1,4),mod1(i+1,4)+4,i+4) for i in 1:4]
    seed=ControlledSurfaceInlet.assemble_surface_mesh(nodes,quads,Set{Tuple{Int,Int}}();
        R=radius,L=.01,ca=1.,slope=0.,geometry_mm=bearing_geometry(p),feed_z_m=.005,inlet_radius_m=.001)
    events=NamedTuple[]
    plain=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=4)
    watched=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=4,observer=s->push!(events,s))
    @test watched==plain
    @test [s.phase for s in events]==[:initial,:final]
    @test last(events).termination_reason==:stationary
    @test last(events).accepted_steps==0
    @test last(events).evaluations==1
    @test last(events).gradient_norm==0
    @test isnothing(last(events).diagnostics)
    @test reduce(vcat,[s.history_chunk for s in events])==plain.history
end

@testset "Exactly unchanged trials stop without accepting fake progress" begin
    # Large chart coordinates make the real first update round to zero despite
    # a nonzero gradient above tolerance. This small numerical fixture avoids
    # thousands of iterations and does not mock the objective or optimizer.
    p=BearingParameters(journal_midpoint_radius=1e7,bearing_length=4e7,
        radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=0.)
    radius=p.housing_midpoint_radius*1e-3
    nodes=[(radius*theta,z) for z in (0.,20000.,40000.) for theta in (0.,pi/2,pi,3pi/2)]
    nodes[2]=(nodes[2][1]+500.,nodes[2][2])
    nodes[7]=(nodes[7][1],nodes[7][2]+500.)
    quads=[(4*j+i,4*j+mod1(i+1,4),4*(j+1)+mod1(i+1,4),4*(j+1)+i)
        for j in 0:1 for i in 1:4]
    seed=ControlledSurfaceInlet.assemble_surface_mesh(nodes,quads,Set{Tuple{Int,Int}}();
        R=radius,L=40000.,ca=1.,slope=0.,geometry_mm=bearing_geometry(p),
        feed_z_m=20000.,inlet_radius_m=.0025)
    selection=ControlledSurfaceInlet.macro_data(seed)
    variables=vec(copy(selection.nodes[selection.free,:]))
    energy,gradient=ControlledSurfaceInlet.layout_energy(seed,selection,variables)
    @test selection.free==[6]
    @test norm(gradient)>1e-5
    @test variables-1e-7*gradient==variables

    events=NamedTuple[]
    plain=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=4,trace_iterations=(1,))
    watched=ControlledSurfaceInlet.controlled_layout(seed;max_iterations=4,
        trace_iterations=(1,),observer=s->push!(events,s),observe_every=1)
    @test watched==plain
    @test plain.termination_reason==:numerical_stagnation
    @test !plain.gradient_tolerance_met
    @test plain.accepted_steps==0
    @test length(plain.history)==1
    @test plain.evaluations==2 # Initial E/g plus the one recorded unchanged trial.
    @test plain.final_energy===energy
    @test plain.final_gradient_norm===norm(gradient)
    @test plain.mesh.nodes==seed.nodes
    @test plain.mesh.areas==seed.areas
    @test plain.mesh.quality==seed.quality
    @test ControlledSurfaceInlet.passes_quality_limits(plain.mesh.quality)
    @test [s.phase for s in events]==[:initial,:final]
    @test reduce(vcat,[s.history_chunk for s in events])==plain.history
    final=last(events)
    @test final.termination_reason==:numerical_stagnation
    @test final.variables==variables
    @test final.gradient==gradient
    @test isempty(final.position_changes) && isempty(final.gradient_changes)
    @test !final.diagnostics.accepted
    @test final.diagnostics.step_norm_m==0
    @test length(final.diagnostics.trials)==1
    @test only(final.diagnostics.trials).predicted_change==0
    @test only(only(plain.trace).trials).variables==variables
end

@testset "DEBUG output and persisted telemetry are observational" begin
    @test length(long_run_cases())==7
    # The Python freezer and Julia provenance must describe the same source
    # closure. In particular, neither may omit shared quality/metric includes.
    queue_source=read(joinpath(@__DIR__,"../examples/queue_inlet_long_runs.py"),String)
    declaration=match(r"SOURCE_FILES = \[(.*?)\]"s,queue_source).captures[1]
    frozen_files=[m.captures[1] for m in eachmatch(r"\"([^\"]+)\"",declaration)]
    @test Set(frozen_files)==Set(long_run_sources())
    root=normpath(joinpath(@__DIR__,"../.."))
    @test all(file->isfile(joinpath(root,file)),frozen_files)
    mktempdir() do folder
        paths=[joinpath(folder,"quiet"),joinpath(folder,"debug")]
        results=Any[]
        for (path,debug) in zip(paths,(false,true))
            open(joinpath(folder,"stderr-$debug.log"),"w") do io
                redirect_stderr(io) do
                    push!(results,run_long_case("reference",path;max_iterations=4,
                        debug,observe_every=1,quality_every=2))
                end
            end
        end
        @test results[1]==results[2]
        @test isempty(read(joinpath(folder,"stderr-false.log"),String))
        @test occursin("phase=final",read(joinpath(folder,"stderr-true.log"),String))
        for (path,result) in zip(paths,results)
            latest=TOML.parsefile(joinpath(path,"latest.toml"))
            @test latest["iteration"]==4
            @test latest["termination_reason"]=="iteration_limit"
            @test latest["energy"]===result.final_energy
            @test length(readlines(joinpath(path,"history.tsv")))==length(result.history)+1
            @test isfile(joinpath(path,"last-trials.tsv"))
            state=deserialize(joinpath(path,"latest-state.jls"))
            saved=deserialize(joinpath(path,"seed.jls"))
            energy,gradient=ControlledSurfaceInlet.layout_energy(saved.seed,saved.selection,state.variables)
            @test energy===result.final_energy
            @test norm(gradient)===result.final_gradient_norm
            @test saved.selection.initial_areas==saved.seed.areas
        end
        @test_throws ArgumentError run_long_case("reference",paths[1];max_iterations=4)
        bad=joinpath(folder,"invalid")
        @test_throws ArgumentError run_long_case("unknown",bad)
        @test !ispath(bad)
    end
end
