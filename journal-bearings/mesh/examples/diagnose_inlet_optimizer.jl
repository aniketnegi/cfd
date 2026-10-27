# Bounded numerical diagnosis, not production progress logging or a new optimizer.
using LinearAlgebra, SparseArrays, Random, Serialization
include(joinpath(@__DIR__, "compare_inlet_iterations.jl"))

# Independent value-only oracle. All arithmetic, including coordinates and fixed
# inputs, is BigFloat. Frozen Float64 area anchors/caches are inputs to the SAME
# objective, not replaced by newly computed anchors. No production geometry or
# derivative routine is called here. Quads use diagonal-cross areas and corner
# cotangents use the cosine rule/Gram determinant instead of the production form.
function reference_layout_context(mesh, selection)
    return (; nodes=BigFloat.(selection.nodes), period=BigFloat(mesh.period),
        quads=mesh.quads, free=selection.free, active=selection.active,
        active_faces=selection.active_faces, owner=mesh.owner, neighbour=mesh.neighbour,
        initial_areas=BigFloat.(selection.initial_areas), anchors=BigFloat.(mesh.areas))
end

function reference_layout_energy(context, variables)
    points = copy(context.nodes)
    points[context.free, :] = reshape(variables, :, 2)
    areas = copy(context.initial_areas)
    energy = BigFloat(0)
    for cell in context.active
        quad = context.quads[cell]
        origin_x, origin_y = points[quad[1], :]
        vertices = ntuple(4) do corner
            dx = points[quad[corner], 1] - origin_x
            dx -= round(dx / context.period) * context.period
            (dx, points[quad[corner], 2] - origin_y)
        end
        a, b, c, d = vertices
        area = ((c[1]-a[1])*(d[2]-b[2]) - (c[2]-a[2])*(d[1]-b[1])) / 2
        area > 0 || return BigFloat(Inf)
        areas[cell] = area
        for corner in 1:4
            centre = vertices[corner]
            next = vertices[mod1(corner+1, 4)]
            previous = vertices[mod1(corner-1, 4)]
            u = (next[1]-centre[1], next[2]-centre[2])
            v = (previous[1]-centre[1], previous[2]-centre[2])
            u[1]*v[2] - u[2]*v[1] > 0 || return BigFloat(Inf)
            u_squared = u[1]^2 + u[2]^2
            v_squared = v[1]^2 + v[2]^2
            opposite_squared = (u[1]-v[1])^2 + (u[2]-v[2])^2
            numerator = (u_squared + v_squared - opposite_squared) / 2
            energy += numerator^2 / (u_squared*v_squared - numerator^2)
        end
    end
    for face in context.active_faces
        energy += log(areas[context.owner[face]] / areas[context.neighbour[face]])^2 / 4
    end
    for cell in context.active
        energy += BigFloat(0.01) * log(areas[cell] / context.anchors[cell])^2
    end
    return energy
end

function reference_directional_derivative(context, variables, direction, step)
    point, tangent = BigFloat.(variables), BigFloat.(direction)
    return (reference_layout_energy(context, point + step*tangent) -
        reference_layout_energy(context, point - step*tangent)) / (2step)
end

function diagnosis_directions(snapshot)
    gradient = snapshot.gradient
    random = rand(MersenneTwister(73), (-1., 1.), length(gradient))
    random ./= norm(random)
    coordinate = zeros(length(gradient))
    coordinate[argmax(abs.(gradient))] = 1
    return ("descent" => -gradient/norm(gradient), "random" => random,
        "coordinate" => coordinate, "search" => snapshot.direction/norm(snapshot.direction))
end

# Collect the production's scalar terms ONLY to isolate accumulation error.
# The left fold must reproduce its returned value exactly. This is a diagnostic
# decomposition, not an alternative objective used by the optimizer or oracle.
function diagnostic_layout_terms(mesh, selection, variables)
    positions = copy(selection.nodes)
    positions[selection.free,:] = reshape(variables,:,2)
    areas = copy(selection.initial_areas)
    terms = Float64[]
    for cell in selection.active
        coordinates = ControlledSurfaceInlet.cell_coordinates(mesh,positions,mesh.quads[cell])
        areas[cell] = ControlledSurfaceInlet.cell_chart_area(coordinates)
        areas[cell] > 0 || return [Inf]
        for corner in 1:4
            u = coordinates[mod1(corner+1,4),:] - coordinates[corner,:]
            v = coordinates[mod1(corner-1,4),:] - coordinates[corner,:]
            cross = u[1]*v[2]-u[2]*v[1]
            cross > 0 || return [Inf]
            push!(terms,(dot(u,v)/cross)^2)
        end
    end
    for face in selection.active_faces
        push!(terms,.25*log(areas[mesh.owner[face]]/areas[mesh.neighbour[face]])^2)
    end
    for cell in selection.active
        push!(terms,.01*log(areas[cell]/mesh.areas[cell])^2)
    end
    return terms
end

function diagnostic_accumulations(mesh, selection, variables, energy)
    terms = diagnostic_layout_terms(mesh,selection,variables)
    sequential = foldl(+,terms;init=0.)
    isequal(sequential,energy) || error("term decomposition disagrees with production")
    total, correction = 0., 0.
    for term in terms
        next = total+term
        correction += abs(total) >= abs(term) ? (total-next)+term : (term-next)+total
        total = next
    end
    return (; pairwise=sum(terms), compensated=total+correction, exact_term_sum=sum(BigFloat,terms))
end

function write_diagnostic_table(path, rows)
    isempty(rows) && error("empty diagnostic table: $path")
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), '\t'))
        for row in rows
            println(io, join(values(row), '\t'))
        end
    end
end

function derivative_diagnosis(mesh, selection, snapshots, context)
    rows = NamedTuple[]
    for snapshot in snapshots, (name, direction) in diagnosis_directions(snapshot)
        analytic = dot(snapshot.gradient, direction)
        reference = reference_directional_derivative(context, snapshot.variables, direction, big"1e-18")
        refined = reference_directional_derivative(context, snapshot.variables, direction, big"1e-22")
        for step in 10. .^ (-5:-1:-13)
            plus = first(ControlledSurfaceInlet.layout_energy(mesh, selection,
                snapshot.variables + step*direction))
            minus = first(ControlledSurfaceInlet.layout_energy(mesh, selection,
                snapshot.variables - step*direction))
            measured = (plus-minus) / (2step)
            push!(rows, (; iteration=snapshot.iteration, direction=name, step_m=step,
                analytic_per_m=analytic, reference_per_m=Float64(refined),
                analytic_error_per_m=Float64(BigFloat(analytic)-refined),
                reference_refinement_error_per_m=Float64(abs(reference-refined)),
                float_difference_per_m=measured, float_error_per_m=Float64(BigFloat(measured)-refined),
                valid_trials=isfinite(plus) && isfinite(minus)))
        end
    end
    return rows
end

function search_diagnosis(mesh, selection, snapshots, context)
    states, pairs, trials = NamedTuple[], NamedTuple[], NamedTuple[]
    for snapshot in snapshots
        reference_energy = reference_layout_energy(context, BigFloat.(snapshot.variables))
        accumulated = diagnostic_accumulations(mesh,selection,snapshot.variables,snapshot.energy)
        push!(states, (; iteration=snapshot.iteration, energy=snapshot.energy,
            gradient_norm_per_m=norm(snapshot.gradient),
            reference_energy=Float64(reference_energy),
            energy_error_ulps=Float64((BigFloat(snapshot.energy)-reference_energy)/eps(snapshot.energy)),
            term_rounding_error_ulps=Float64((accumulated.exact_term_sum-reference_energy)/eps(snapshot.energy)),
            pairwise_energy_error_ulps=Float64((BigFloat(accumulated.pairwise)-reference_energy)/eps(snapshot.energy)),
            compensated_energy_error_ulps=Float64((BigFloat(accumulated.compensated)-reference_energy)/eps(snapshot.energy)),
            inverse_curvature_scale_m2=snapshot.inverse_curvature_scale,
            gradient_fallback=snapshot.gradient_fallback,
            unclipped_node_motion_m=snapshot.maximum_unclipped_motion_m,
            step_limit_m=snapshot.step_limit_m,
            direction_norm_m=norm(snapshot.direction),
            gradient_dot_direction=dot(snapshot.gradient, snapshot.direction),
            memory_pairs=length(snapshot.position_changes)))
        for (index, (s,y)) in enumerate(zip(snapshot.position_changes, snapshot.gradient_changes))
            sy = dot(s,y)
            push!(pairs, (; iteration=snapshot.iteration, pair=index, position_norm_m=norm(s),
                gradient_change_norm_per_m=norm(y), curvature_product=sy,
                curvature_cosine=sy/(norm(s)*norm(y)), secant_curvature_per_m2=sy/dot(s,s)))
        end
        for trial in snapshot.trials
            change = trial.variables - snapshot.variables
            prediction = dot(snapshot.gradient, change)
            bound = snapshot.energy + 1e-4*prediction
            exact_change = BigFloat.(trial.variables) - BigFloat.(snapshot.variables)
            reference_trial = reference_layout_energy(context, BigFloat.(trial.variables))
            reference_prediction = dot(BigFloat.(snapshot.gradient), exact_change)
            accumulated_trial = diagnostic_accumulations(mesh,selection,trial.variables,trial.energy)
            push!(trials, (; iteration=snapshot.iteration, backtrack=trial.backtrack,
                step_norm_m=norm(change), maximum_node_motion_m=maximum(sqrt.(sum(reshape(change,:,2).^2,dims=2))),
                changed_coordinates=count(!iszero, change), predicted_change=prediction,
                required_change=1e-4*prediction, rounded_bound_change=bound-snapshot.energy,
                float_energy_change=trial.energy-snapshot.energy,
                reference_energy_change=Float64(reference_trial-reference_energy),
                pairwise_energy_change=accumulated_trial.pairwise-accumulated.pairwise,
                compensated_energy_change=accumulated_trial.compensated-accumulated.compensated,
                energy_ulp=eps(snapshot.energy),
                pairwise_armijo_pass=accumulated_trial.pairwise <= accumulated.pairwise+1e-4*prediction,
                compensated_armijo_pass=accumulated_trial.compensated <= accumulated.compensated+1e-4*prediction,
                float_armijo_pass=isfinite(trial.energy) && trial.energy <= bound,
                reference_armijo_pass=isfinite(reference_trial) &&
                    reference_trial <= reference_energy + BigFloat(1e-4)*reference_prediction))
        end
    end
    return (; states, pairs, trials)
end

# A smooth displacement vanishing towards fixed vertices, constructed from the
# graph Dirichlet Laplacian. It probes low spatial frequencies, NOT a proposed
# preconditioner. Curvature ratios below are directional samples, not a measured
# full Hessian condition number (nor proof that the Hessian is positive definite).
function smooth_diagnosis_direction(mesh, selection)
    free_index = Dict(node => index for (index,node) in enumerate(selection.free))
    n = length(selection.free)
    rows, columns, values = Int[], Int[], Float64[]
    for (a,b) in mesh.edge_nodes, (node, other) in ((a,b),(b,a))
        haskey(free_index,node) || continue
        i = free_index[node]
        push!(rows,i); push!(columns,i); push!(values,1.)
        if haskey(free_index,other)
            push!(rows,i); push!(columns,free_index[other]); push!(values,-1.)
        end
    end
    envelope = sparse(rows,columns,values,n,n) \ ones(n)
    direction = vcat(envelope,zeros(n))
    return direction/norm(direction)
end

function curvature_diagnosis(mesh, selection, snapshots)
    rows = NamedTuple[]
    smooth = smooth_diagnosis_direction(mesh, selection)
    for snapshot in snapshots
        directions = (diagnosis_directions(snapshot)..., "smooth-x" => smooth)
        for (name,direction) in directions, step in (1e-8,1e-9)
            _, plus = ControlledSurfaceInlet.layout_energy(mesh,selection,snapshot.variables+step*direction)
            _, minus = ControlledSurfaceInlet.layout_energy(mesh,selection,snapshot.variables-step*direction)
            action = (plus-minus)/(2step)
            push!(rows, (; iteration=snapshot.iteration, direction=name, step_m=step,
                directional_curvature_per_m2=dot(direction,action), action_norm_per_m2=norm(action)))
        end
    end
    return rows
end

function diagnose_optimizer(directory; name="reference", budget=720)
    ispath(directory) && throw(ArgumentError("output directory already exists: $directory"))
    budget isa Int && budget >= 2 || throw(ArgumentError("diagnosis budget must be an integer >= 2"))
    cases = (iteration_cases("coarse")..., iteration_cases("fine")...)
    index = findfirst(pair -> first(pair)==name, cases)
    isnothing(index) && throw(ArgumentError("unknown diagnosis case: $name"))
    g = last(cases[index])
    p = BearingParameters(journal_midpoint_radius=g.radius,bearing_length=g.length,
        radial_clearance=g.clearance,semicone_angle_deg=g.angle,eccentricity_ratio=g.eccentricity_ratio)
    mkpath(directory)
    sources = ("../../bearing.jl","../2d/conical/grid.jl","../2d/conical/surface_feed.jl",
        "../2d/conical/quality.jl","../2d/conical/inlet/seed.jl","../2d/conical/metrics.jl",
        "../2d/conical/inlet/placement.jl","../2d/conical/inlet/controlled.jl",
        "compare_inlet_iterations.jl","diagnose_inlet_optimizer.jl")
    metadata = Dict("name"=>name,"budget"=>budget,"stability_steps"=>0,"precision_bits"=>256,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "geometry"=>Dict(string(k)=>v for (k,v) in pairs(g)),
        "source_sha256"=>Dict(file=>bytes2hex(sha256(read(joinpath(@__DIR__,file)))) for file in sources))
    open(joinpath(directory,"inputs.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    initial = ControlledSurfaceInlet.reference_seed(p,g.n_theta,g.n_axial;
        feed_diameter=g.diameter,feed_z=g.feed_z)
    selection = ControlledSurfaceInlet.macro_data(initial)
    # Request all late steps of the original coarse failure, plus common budgets.
    result = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=budget,stability_steps=0,
        trace_iterations=(1,360,665,666,budget))
    serialize(joinpath(directory,"trace.bin"), result.trace)
    # If a different case stops early, capture its actual last attempted step by
    # replaying from the ORIGINAL seed, never by optimizing the terminal mesh.
    terminal_iteration = length(result.history)
    if isempty(result.trace) || last(result.trace).iteration != terminal_iteration
        result = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=budget,stability_steps=0,
            trace_iterations=(1,360,terminal_iteration))
        serialize(joinpath(directory,"trace.bin"), result.trace)
    end
    snapshots = unique(entry -> entry.iteration,
        [first(result.trace), result.trace[argmin(abs.([t.iteration-360 for t in result.trace]))],last(result.trace)])
    setprecision(BigFloat,256) do
        context = reference_layout_context(initial,selection)
        derivatives = derivative_diagnosis(initial,selection,snapshots,context)
        search = search_diagnosis(initial,selection,result.trace,context)
        write_diagnostic_table(joinpath(directory,"derivatives.tsv"), derivatives)
        write_diagnostic_table(joinpath(directory,"states.tsv"), search.states)
        write_diagnostic_table(joinpath(directory,"curvature-pairs.tsv"), search.pairs)
        write_diagnostic_table(joinpath(directory,"trials.tsv"), search.trials)
        write_diagnostic_table(joinpath(directory,"curvature-probes.tsv"),
            curvature_diagnosis(initial,selection,snapshots))
    end
    open(joinpath(directory,"result.toml"),"w") do io
        TOML.print(io,Dict("termination_reason"=>string(result.termination_reason),
            "accepted_steps"=>result.accepted_steps,"iterations_recorded"=>length(result.history),
            "final_energy"=>result.final_energy,"final_gradient_norm_per_m"=>result.final_gradient_norm,
            "quality_pass"=>ControlledSurfaceInlet.passes_quality_limits(result.mesh.quality));sorted=true)
    end
    println("$name: stop=$(result.termination_reason), steps=$(result.accepted_steps), " *
        "E=$(result.final_energy), |g|=$(result.final_gradient_norm); diagnostics=$directory")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 2 || error("usage: diagnose_inlet_optimizer.jl [FRESH_DIRECTORY] [CASE_NAME]")
    name = length(ARGS)==2 ? ARGS[2] : "reference"
    directory = isempty(ARGS) ? joinpath(@__DIR__,"../../out/cases/conical-inlet-diagnosis-$name") : ARGS[1]
    diagnose_optimizer(directory;name)
end
