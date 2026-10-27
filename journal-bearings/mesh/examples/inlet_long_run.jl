# Fixed-budget L-BFGS study (stability_steps=0), not public acceptance.
using Serialization, Dates
include(joinpath(@__DIR__, "compare_inlet_iterations.jl"))

function long_run_cases()
    coarse, fine = iteration_cases("coarse"), iteration_cases("fine")
    reference = last(first(coarse))
    held = (radius=55.,length=110.,clearance=.07,angle=8.,eccentricity_ratio=.6,
        diameter=6.,feed_z=62.3,n_theta=129,n_axial=67)
    return (fine[1],coarse... ,"held-out"=>held,
        "reference-medium"=>(;reference...,n_theta=256,n_axial=128),fine[2])
end

function long_run_sources()
    return ("bearing.jl","mesh/2d/quality.jl","mesh/2d/grid.jl","mesh/2d/surface_feed.jl",
        "mesh/2d/inlet/controlled.jl","mesh/2d/metrics.jl",
        "mesh/2d/inlet/seed.jl","mesh/2d/inlet/placement.jl",
        "mesh/examples/compare_inlet_iterations.jl","mesh/examples/inlet_long_run.jl")
end

function atomic_toml(path, data)
    temporary = path * ".tmp"
    open(temporary,"w") do io
        TOML.print(io,data;sorted=true)
    end
    mv(temporary,path;force=true)
end

function atomic_state(path, state)
    temporary = path * ".tmp"
    serialize(temporary,state)
    mv(temporary,path;force=true)
end

function append_record(io, record)
    position(io)==0 && println(io,join(string.(keys(record)),'\t'))
    println(io,join(values(record),'\t'))
    flush(io)
end

study_layout(seed, selection, variables) =
    ControlledSurfaceInlet.assemble_layout(seed, selection, variables)

function free_cartesian(p, seed, selection, variables)
    positions = reshape(variables,:,2)
    return map(eachrow(positions)) do row
        z = ControlledSurfaceInlet.chart_axial_position_mm(seed.geometry_mm,seed.ca,row[2])
        bearing_point(p,row[1]/seed.R,z,1.)
    end
end

node_distances(a,b) = [sqrt(sum((x-y)^2 for (x,y) in zip(first,second)))
    for (first,second) in zip(a,b)]

function study_observer(directory,p,seed,selection;debug=false,quality_every=500)
    quality_every > 0 || throw(ArgumentError("quality_every must be positive"))
    history_io = open(joinpath(directory,"history.tsv"),"w")
    progress_io = open(joinpath(directory,"progress.tsv"),"w")
    quality_io = open(joinpath(directory,"quality.tsv"),"w")
    initial_points = physical_study_nodes(p,seed)
    local_scale = fill(Inf,length(seed.nodes))
    for (a,b) in seed.edge_nodes
        distance = sqrt(sum((x-y)^2 for (x,y) in zip(initial_points[a],initial_points[b])))
        local_scale[a] = min(local_scale[a],distance)
        local_scale[b] = min(local_scale[b],distance)
    end
    scale = local_scale[selection.free]
    previous_points = initial_points[selection.free]
    previous_energy = nothing
    started = time_ns()
    last_quality = Dict{String,Any}()
    snapshot_iteration = -1
    mkpath(joinpath(directory,"states"))

    function observe(state)
        for row in state.history_chunk
            append_record(history_io,row)
        end
        points = free_cartesian(p,seed,selection,state.variables)
        movement = node_distances(previous_points,points)
        from_initial = node_distances(initial_points[selection.free],points)
        diagnostics = state.diagnostics
        trial_count = isnothing(diagnostics) ? 0 : length(diagnostics.trials)
        inverse_scale = isnothing(diagnostics) ? NaN : diagnostics.inverse_curvature_scale
        step_norm = isnothing(diagnostics) ? 0. : diagnostics.step_norm_m
        maximum_step = isnothing(diagnostics) ? 0. : diagnostics.maximum_chart_step_m
        check_quality = state.phase != :update || state.iteration % quality_every == 0 ||
            state.iteration in (100,180,360,720)
        if check_quality
            layout = state.phase == :initial ? seed : study_layout(seed,selection,state.variables)
            q = layout.quality
            last_quality = Dict("quality_iteration"=>state.iteration,
                "quality_pass"=>ControlledSurfaceInlet.passes_quality_limits(q),
                "minimum_corner"=>q.minimum_scaled_corner,"minimum_alignment"=>q.minimum_face_alignment,
                "maximum_chart_area_ratio"=>q.maximum_neighbour_area_ratio,
                "maximum_edge_ratio"=>q.maximum_edge_ratio)
            append_record(quality_io,(;phase=string(state.phase),iteration=state.iteration,accepted_steps=state.accepted_steps,
                energy=state.energy,gradient_norm=state.gradient_norm,
                quality_pass=last_quality["quality_pass"],minimum_corner=q.minimum_scaled_corner,
                minimum_alignment=q.minimum_face_alignment,maximum_chart_area_ratio=q.maximum_neighbour_area_ratio,
                maximum_edge_ratio=q.maximum_edge_ratio,
                corner_violations=count(x->!(x>.9),layout.corner_scores),
                alignment_violations=count(x->!(x>.9),layout.face_alignment),
                area_ratio_violations=count(x->!(x<2),layout.neighbour_area_ratios)))
            checkpoint = joinpath(directory,"states","iteration-$(lpad(state.iteration,6,'0'))-$(state.phase).jls")
            atomic_state(checkpoint,state)
            cp(checkpoint,joinpath(directory,"latest-state.jls.tmp");force=true)
            mv(joinpath(directory,"latest-state.jls.tmp"),joinpath(directory,"latest-state.jls");force=true)
            snapshot_iteration = state.iteration
            if !isnothing(diagnostics)
                open(joinpath(directory,"last-trials.tsv.tmp"),"w") do io
                    for trial in diagnostics.trials
                        append_record(io,trial)
                    end
                end
                mv(joinpath(directory,"last-trials.tsv.tmp"),joinpath(directory,"last-trials.tsv");force=true)
            end
        end
        elapsed = (time_ns()-started)*1e-9
        record = (;phase=string(state.phase),iteration=state.iteration,accepted_steps=state.accepted_steps,
            evaluations=state.evaluations,elapsed_seconds=elapsed,energy=state.energy,
            gradient_norm=state.gradient_norm,scaled_gradient_norm=state.scaled_gradient_norm,
            energy_change=isnothing(previous_energy) ? NaN : state.energy-previous_energy,
            maximum_move_since_report_um=maximum(movement)*1e3,
            maximum_move_since_report_local_edge=maximum(movement./scale),
            maximum_move_from_seed_um=maximum(from_initial)*1e3,
            interval_chart_motion_bound_m=state.interval_motion_bound_m,
            interval_maximum_chart_step_m=state.interval_max_step_m,
            last_step_norm_m=step_norm,last_maximum_chart_step_m=maximum_step,
            last_trial_count=trial_count,inverse_curvature_scale=inverse_scale,
            memory_pairs=length(state.position_changes),snapshot_iteration,
            termination_reason=string(state.termination_reason))
        append_record(progress_io,record)
        status = merge(Dict(string(k)=>v for (k,v) in pairs(record)),last_quality,
            Dict("pid"=>getpid(),"updated_utc"=>string(now(UTC))))
        atomic_toml(joinpath(directory,"latest.toml"),status)
        if debug
            println(stderr,join(("$key=$value" for (key,value) in pairs(record))," "))
            flush(stderr)
        end
        previous_points,previous_energy = points,state.energy
        return nothing
    end
    close_observer() = foreach(close,(history_io,progress_io,quality_io))
    return observe,close_observer
end

function run_long_case(name,directory;max_iterations=20000,debug=false,observe_every=20,quality_every=500)
    ispath(directory) && throw(ArgumentError("output directory already exists: $directory"))
    max_iterations > 0 || throw(ArgumentError("iteration budget must be positive"))
    cases = long_run_cases()
    index = findfirst(pair->first(pair)==name,cases)
    isnothing(index) && throw(ArgumentError("unknown long-run case: $name"))
    g = last(cases[index])
    mkpath(directory)
    source_root = normpath(joinpath(@__DIR__,"../.."))
    metadata = Dict("name"=>name,"max_iterations"=>max_iterations,"stability_steps"=>0,"observe_every"=>observe_every,
        "quality_every"=>quality_every,"debug"=>debug,"pid"=>getpid(),"started_utc"=>string(now(UTC)),
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "geometry"=>Dict(string(k)=>v for (k,v) in pairs(g)),
        "source_sha256"=>Dict(path=>bytes2hex(sha256(read(joinpath(source_root,path)))) for path in long_run_sources()))
    atomic_toml(joinpath(directory,"inputs.toml"),metadata)
    atomic_toml(joinpath(directory,"latest.toml"),Dict("phase"=>"starting","pid"=>getpid(),"iteration"=>0))
    p = BearingParameters(journal_midpoint_radius=g.radius,bearing_length=g.length,
        radial_clearance=g.clearance,semicone_angle_deg=g.angle,eccentricity_ratio=g.eccentricity_ratio)
    seed = ControlledSurfaceInlet.reference_seed(p,g.n_theta,g.n_axial;feed_diameter=g.diameter,feed_z=g.feed_z)
    selection = ControlledSurfaceInlet.macro_data(seed)
    # Original mesh/selection preserve the original area anchors. State snapshots
    # are inspectable Julia-version-dependent data, NOT a terminal-mesh restart API.
    atomic_state(joinpath(directory,"seed.jls"),(;parameters=p,seed,selection))
    observer,finish = study_observer(directory,p,seed,selection;debug,quality_every)
    result = try
        ControlledSurfaceInlet.controlled_layout(seed;max_iterations,observer,observe_every,stability_steps=0)
    finally
        finish()
    end
    atomic_toml(joinpath(directory,"result.toml"),Dict("termination_reason"=>string(result.termination_reason),
        "accepted_steps"=>result.accepted_steps,"iterations_recorded"=>length(result.history),
        "evaluations"=>result.evaluations,"final_energy"=>result.final_energy,
        "final_gradient_norm"=>result.final_gradient_norm,
        "scaled_final_gradient_norm"=>result.scaled_final_gradient_norm,
        "gradient_tolerance_met"=>result.gradient_tolerance_met,
        "quality_pass"=>ControlledSurfaceInlet.passes_quality_limits(result.mesh.quality)))
    println("case=$name stop=$(result.termination_reason) accepted_steps=$(result.accepted_steps) " *
        "evaluations=$(result.evaluations) energy=$(result.final_energy) directory=$directory")
    return result
end

if abspath(PROGRAM_FILE)==@__FILE__
    args = filter(!=("--debug"),ARGS)
    2 <= length(args) <= 3 || error("usage: inlet_long_run.jl CASE FRESH_DIRECTORY [MAX_ITERATIONS] [--debug]")
    debug = "--debug" in ARGS || lowercase(get(ENV,"DEBUG","0")) in ("1","true","yes")
    budget = length(args)==3 ? parse(Int,args[3]) : 20000
    run_long_case(args[1],abspath(args[2]);max_iterations=budget,debug)
end
