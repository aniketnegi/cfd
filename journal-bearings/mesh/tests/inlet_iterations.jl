using Test
include(joinpath(@__DIR__, "../examples/compare_inlet_iterations.jl"))

@testset "Optimizer termination and terminal measurements" begin
    p = BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
        semicone_angle_deg=10.,eccentricity_ratio=.6)
    initial = ControlledSurfaceInlet.reference_seed(p,64,32;feed_diameter=5.,feed_z=50.)
    selection = ControlledSurfaceInlet.macro_data(initial)
    original_nodes = copy(initial.nodes)
    original_areas = copy(initial.areas)
    previous = nothing
    for budget in (1,2,4)
        result = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=budget,stability_steps=0)
        @test result.termination_reason == :iteration_limit
        @test result.accepted_steps == budget
        @test length(result.history) == budget
        @test result.final_energy <= last(result.history).energy
        @test all(diff([entry.energy for entry in result.history]) .<= 0)
        @test result.gradient_tolerance === 1e-5
        terminal = study_terminal_state(initial,selection,result)
        @test terminal.energy === result.final_energy
        @test terminal.gradient_norm === result.final_gradient_norm
        @test result.scaled_final_gradient_norm === result.final_gradient_norm*(initial.inlet_radius_m/.0025)
        @test result.gradient_tolerance_met == (result.scaled_final_gradient_norm < result.gradient_tolerance)
        @test initial.nodes == original_nodes
        @test initial.areas == original_areas
        if !isnothing(previous)
            @test result.history[1:length(previous.history)] == previous.history
            @test result.history[previous.accepted_steps+1].energy === previous.final_energy
            @test result.history[previous.accepted_steps+1].gradient_norm === previous.final_gradient_norm
        end
        previous = result
    end

    # Exercise an actual longer run without baking a platform-specific stopping
    # iteration into the test. Every recorded state has either an accepted step
    # or the one final gradient/line-search/stagnation stop.
    long_initial = ControlledSurfaceInlet.reference_seed(p,128,64;feed_diameter=5.,feed_z=50.)
    long_result = ControlledSurfaceInlet.controlled_layout(long_initial;max_iterations=720,stability_steps=0)
    @test long_result.termination_reason in (:iteration_limit,:gradient_tolerance,:line_search_failure,:numerical_stagnation)
    stopped_before_step = long_result.termination_reason != :iteration_limit
    @test long_result.accepted_steps == length(long_result.history)-Int(stopped_before_step)
    if stopped_before_step
        @test long_result.final_energy === last(long_result.history).energy
        @test long_result.final_gradient_norm === last(long_result.history).gradient_norm
    end
    if long_result.termination_reason == :gradient_tolerance
        @test long_result.gradient_tolerance_met
    end
    terminal = study_terminal_state(long_initial,ControlledSurfaceInlet.macro_data(long_initial),long_result)
    @test terminal.energy === long_result.final_energy

    # A mesh with no free vertices terminates before attempting any step.
    fixed_parameters = BearingParameters(journal_midpoint_radius=50.,bearing_length=10.,
        radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=0.)
    radius = fixed_parameters.housing_midpoint_radius*1e-3
    nodes = [(radius*theta,z) for z in (0.,.01) for theta in (0.,pi/2,pi,3pi/2)]
    quads = [(i,mod1(i+1,4),mod1(i+1,4)+4,i+4) for i in 1:4]
    fixed = ControlledSurfaceInlet.assemble_surface_mesh(nodes,quads,Set{Tuple{Int,Int}}();
        R=radius,L=.01,ca=1.,slope=0.,geometry_mm=bearing_geometry(fixed_parameters),
        feed_z_m=.005,inlet_radius_m=.001)
    result = ControlledSurfaceInlet.controlled_layout(fixed;max_iterations=4,stability_steps=0)
    @test result.termination_reason == :gradient_tolerance
    @test result.accepted_steps == 0
    @test length(result.history) == 1
    @test result.free_nodes == 0
    @test result.final_gradient_norm == 0
    @test result.gradient_tolerance_met
    @test result.mesh.nodes == nodes
end

@testset "One unchanged surface-quality predicate" begin
    quality = (minimum_scaled_corner=.95,minimum_face_alignment=.95,maximum_neighbour_area_ratio=1.5)
    @test ControlledSurfaceInlet.passes_quality_limits(quality)
    for changes in ((minimum_scaled_corner=.9,), (minimum_face_alignment=.9,),
        (maximum_neighbour_area_ratio=2.,), (minimum_scaled_corner=NaN,))
        @test !ControlledSurfaceInlet.passes_quality_limits(merge(quality,changes))
    end
    @test ControlledSurfaceInlet.passes_quality_limits((minimum_scaled_corner=nextfloat(.9),
        minimum_face_alignment=nextfloat(.9),maximum_neighbour_area_ratio=prevfloat(2.)))
end

@testset "Iteration study records independent runs without overwriting" begin
    @test length(iteration_cases("coarse")) == 3
    @test length(iteration_cases("fine")) == 2
    @test_throws ArgumentError iteration_cases("unknown")
    geometry = (;last(first(iteration_cases()))...,n_theta=64,n_axial=32)
    mktempdir() do directory
        output = joinpath(directory,"study")
        records = compare_iterations(output;cases=("smoke"=>geometry,),budgets=(1,2,4))
        @test length(records) == 3
        @test all(record -> record.termination_reason == "iteration_limit",records)
        @test all(record -> record.accepted_steps == record.budget,records)
        @test all(records) do record
            record.quality_pass == (record.minimum_corner > .9 && record.minimum_alignment > .9 &&
                record.maximum_chart_area_ratio < 2)
        end
        @test isnan(records[1].move_previous_um)
        @test records[1].move_baseline_um == 0
        @test records[end].move_baseline_um > 0
        @test records[end].final_energy <= records[1].final_energy
        @test length(readlines(joinpath(output,"comparison.tsv"))) == 4
        metadata = TOML.parsefile(joinpath(output,"inputs.toml"))
        @test metadata["budgets"] == [1,2,4]
        @test metadata["cases"][1]["name"] == "smoke"
        @test haskey(metadata["source_sha256"],"../2d/inlet/placement.jl")
        for record in records
            history = readlines(joinpath(output,"smoke-$(record.budget)-history.tsv"))
            nodes = readlines(joinpath(output,"smoke-$(record.budget)-chart-nodes.tsv"))
            @test length(history) == record.iterations_recorded+1
            @test length(nodes) > record.cells
        end
        table = joinpath(output,"comparison.tsv")
        before = sha256(read(table))
        @test_throws ArgumentError compare_iterations(output;cases=("smoke"=>geometry,),budgets=(1,2,4))
        @test sha256(read(table)) == before
        for budgets in ((),(0,),(2,1),(1,1),(1.,2.))
            invalid = joinpath(directory,"invalid")
            @test_throws ArgumentError compare_iterations(invalid;budgets)
            @test !ispath(invalid)
        end
    end
end
