# Diagnose iteration-budget sensitivity without changing production defaults.
using LinearAlgebra: norm, BLAS
using SHA, TOML
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))

function iteration_cases(group="coarse")
    reference = (radius=50.,length=100.,clearance=.05,angle=10.,eccentricity_ratio=.6,
        diameter=5.,feed_z=50.,n_theta=128,n_axial=64)
    if group == "coarse"
        return ("reference" => reference,
            "diameter-8" => (;reference...,diameter=8.),
            "combined-odd" => (radius=45.,length=90.,clearance=.06,angle=12.,
                eccentricity_ratio=.6,diameter=4.,feed_z=35.3,n_theta=127,n_axial=65))
    elseif group == "fine"
        return ("reference-fine" => (;reference...,n_theta=512,n_axial=256),
            "squarer-fine" => (;reference...,n_theta=640,n_axial=208))
    end
    throw(ArgumentError("study group must be coarse or fine"))
end

function physical_study_nodes(p, mesh)
    return map(mesh.nodes) do (x,y)
        theta = x / mesh.R
        z = ControlledSurfaceInlet.chart_axial_position_mm(mesh.geometry_mm,mesh.ca,y)
        bearing_point(p,theta,z,1.)
    end
end

function maximum_node_distance_mm(first, second)
    return maximum(sqrt(sum((a-b)^2 for (a,b) in zip(old,new)))
        for (old,new) in zip(first,second))
end

function study_terminal_state(initial, selection, result)
    positions = permutedims(hcat(collect.(result.mesh.nodes)...))
    variables = vec(copy(positions[selection.free,:]))
    # Preserve the original area anchors, not the rebuilt mesh's terminal areas.
    energy, gradient = ControlledSurfaceInlet.layout_energy(initial,selection,variables)
    gradient_norm = norm(gradient)
    energy === result.final_energy || error("terminal objective disagrees with optimizer")
    gradient_norm === result.final_gradient_norm || error("terminal gradient disagrees with optimizer")
    return (;energy,gradient_norm)
end

function compare_iterations(output_directory; cases=iteration_cases(), budgets=(180,360,720))
    ispath(output_directory) && throw(ArgumentError("output directory already exists: $output_directory"))
    !isempty(budgets) && all(b -> b isa Int && b > 0,budgets) &&
        all(>(0),diff(collect(budgets))) ||
        throw(ArgumentError("iteration budgets must be positive, strictly increasing integers"))
    isempty(cases) && throw(ArgumentError("at least one study case is required"))
    names = first.(cases)
    length(unique(names)) == length(names) || throw(ArgumentError("case names must be unique"))
    all(name -> !isempty(name) && all(c -> isletter(c) || isdigit(c) || c in ('-','_'),name),names) ||
        throw(ArgumentError("case names must be safe single path components"))
    mkpath(output_directory)

    source_files = ("../../bearing.jl","../2d/conical/quality.jl","../2d/conical/grid.jl","../2d/conical/surface_feed.jl",
        "../2d/conical/inlet/seed.jl","../2d/conical/metrics.jl",
        "../2d/conical/inlet/placement.jl","../2d/conical/inlet/controlled.jl",
        "compare_inlet_iterations.jl")
    metadata = Dict("julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(), "stability_steps"=>0, "budgets" => collect(budgets),
        "source_sha256" => Dict(file => bytes2hex(sha256(read(joinpath(@__DIR__,file)))) for file in source_files),
        "cases" => [merge(Dict("name"=>name),Dict(string(k)=>v for (k,v) in pairs(g))) for (name,g) in cases])
    open(joinpath(output_directory,"inputs.toml"),"w") do io
        TOML.print(io,metadata; sorted=true)
    end

    records = NamedTuple[]
    open(joinpath(output_directory,"comparison.tsv"),"w") do table
        for (name,g) in cases
            p = BearingParameters(journal_midpoint_radius=g.radius,bearing_length=g.length,
                radial_clearance=g.clearance,semicone_angle_deg=g.angle,
                eccentricity_ratio=g.eccentricity_ratio)
            seed_seconds = @elapsed initial = ControlledSurfaceInlet.reference_seed(p,g.n_theta,g.n_axial;
                feed_diameter=g.diameter,feed_z=g.feed_z)
            selection = ControlledSurfaceInlet.macro_data(initial)
            baseline_nodes = previous_nodes = nothing
            baseline_areas = previous_areas = nothing
            for budget in budgets
                # Every run restarts at the same seed; L-BFGS is not restarted
                # part-way through a longer run with its memory discarded.
                placement_seconds = @elapsed result = ControlledSurfaceInlet.controlled_layout(initial;
                    max_iterations=budget,stability_steps=0)
                evaluation_seconds = @elapsed terminal = study_terminal_state(initial,selection,result)
                nodes = physical_study_nodes(p,result.mesh)
                if isnothing(baseline_nodes)
                    baseline_nodes,baseline_areas = nodes,result.mesh.areas
                end
                move_previous_um = isnothing(previous_nodes) ? NaN :
                    1e3 * maximum_node_distance_mm(previous_nodes,nodes)
                area_previous = isnothing(previous_areas) ? NaN :
                    maximum(abs.(result.mesh.areas.-previous_areas)./previous_areas)
                quality = result.mesh.quality
                record = (;
                    name,n_theta=g.n_theta,n_axial=g.n_axial,diameter_mm=g.diameter,
                    budget,iterations_recorded=length(result.history),accepted_steps=result.accepted_steps,
                    termination_reason=string(result.termination_reason),cells=length(result.mesh.quads),
                    free_nodes=result.free_nodes,seed_seconds,placement_seconds,evaluation_seconds,
                    initial_energy=result.initial_energy,final_energy=terminal.energy,
                    gradient_norm_per_m=terminal.gradient_norm,
                    scaled_gradient_norm_per_m=result.scaled_final_gradient_norm,
                    gradient_tolerance_per_m=result.gradient_tolerance,
                    gradient_tolerance_met=result.gradient_tolerance_met,
                    minimum_corner=quality.minimum_scaled_corner,
                    minimum_alignment=quality.minimum_face_alignment,
                    maximum_chart_area_ratio=quality.maximum_neighbour_area_ratio,
                    maximum_edge_ratio=quality.maximum_edge_ratio,
                    quality_pass=ControlledSurfaceInlet.passes_quality_limits(quality),
                    move_previous_um,
                    move_baseline_um=1e3*maximum_node_distance_mm(baseline_nodes,nodes),
                    chart_area_change_previous=area_previous,
                    chart_area_change_baseline=maximum(abs.(result.mesh.areas.-baseline_areas)./baseline_areas),
                )
                isempty(records) && println(table,join(string.(keys(record)),'\t'))
                println(table,join(values(record),'\t'))
                flush(table)
                push!(records,record)

                prefix = joinpath(output_directory,"$name-$budget")
                open("$prefix-history.tsv","w") do io
                    println(io,"iteration\tenergy_before_step\tgradient_norm_before_step_per_m")
                    for entry in result.history
                        println(io,join((entry.iteration,entry.energy,entry.gradient_norm),'\t'))
                    end
                end
                # Diagnostic coordinates only, not an accepted VTK mesh. Keep
                # rejected terminal layouts visible in the table and this file.
                open("$prefix-chart-nodes.tsv","w") do io
                    println(io,"node\tchart_x_m\tchart_y_m")
                    for (node,(x,y)) in enumerate(result.mesh.nodes)
                        println(io,"$node\t$x\t$y")
                    end
                end
                println("$name $budget: E=$(terminal.energy), scaled |g|=$(result.scaled_final_gradient_norm), " *
                    "stop=$(result.termination_reason), quality=$(record.quality_pass), " *
                    "move from previous=$(move_previous_um) µm; $(placement_seconds) s")
                flush(stdout)
                previous_nodes,previous_areas = nodes,result.mesh.areas
            end
        end
    end
    return records
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 2 || error("usage: compare_inlet_iterations.jl [FRESH_DIRECTORY] [coarse|fine]")
    group = length(ARGS)==2 ? ARGS[2] : "coarse"
    directory = isempty(ARGS) ? joinpath(@__DIR__,"../../out/cases/conical-inlet-iterations-$group") : ARGS[1]
    compare_iterations(directory; cases=iteration_cases(group))
end
