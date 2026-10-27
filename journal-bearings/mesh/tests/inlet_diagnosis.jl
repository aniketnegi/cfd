using Test
include(joinpath(@__DIR__, "../examples/diagnose_inlet_optimizer.jl"))

@testset "Selected trace is observational and owns its data" begin
    p = BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
        radial_clearance=.05,semicone_angle_deg=10.,eccentricity_ratio=.6)
    initial = ControlledSurfaceInlet.reference_seed(p,64,32;feed_diameter=5.,feed_z=50.)
    selection = ControlledSurfaceInlet.macro_data(initial)
    plain = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=4)
    traced = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=4,trace_iterations=(1,3,4,99))
    @test isempty(plain.trace)
    @test [entry.iteration for entry in traced.trace] == [1,3,4]
    for key in keys(plain)
        key == :trace && continue
        @test isequal(getproperty(plain,key),getproperty(traced,key))
    end
    for snapshot in traced.trace
        energy, gradient = ControlledSurfaceInlet.layout_energy(initial,selection,snapshot.variables)
        @test snapshot.energy === energy
        @test snapshot.gradient == gradient
        @test energy === plain.history[snapshot.iteration].energy
        @test norm(gradient) === plain.history[snapshot.iteration].gradient_norm
        @test dot(snapshot.gradient,snapshot.direction) < 0
        @test 1 <= length(snapshot.trials) <= 23
        @test length(snapshot.position_changes) == length(snapshot.gradient_changes)
        for trial in snapshot.trials
            @test trial.energy === first(ControlledSurfaceInlet.layout_energy(initial,selection,trial.variables))
            @test trial.variables !== snapshot.variables
        end
    end
    # Perturb diagnostic copies, including pair vectors, after the run.
    traced.trace[1].variables[1] += 1
    traced.trace[1].gradient[1] += 1
    traced.trace[1].trials[1].variables[1] += 1
    traced.trace[end].position_changes[1][1] += 1
    @test traced.mesh == plain.mesh
    @test traced.history == plain.history
    @test initial.nodes == [Tuple(row) for row in eachrow(selection.nodes)]
end

@testset "Independent value oracle: shear, translation, seam and inversion" begin
    setprecision(BigFloat,256) do
        nodes = BigFloat[0 0; 1 0; 1 1; 0 1]
        context = (;nodes,period=BigFloat(10),quads=[(1,2,3,4)],free=[3,4],active=[1],
            active_faces=Int[],owner=Int[],neighbour=Int[],initial_areas=BigFloat[1],anchors=BigFloat[1])
        square = BigFloat[1,0,1,1]
        @test reference_layout_energy(context,square) == 0
        shear = BigFloat(1)/8
        variables = BigFloat[1+shear,shear,1,1]
        # Unit-area parallelogram: four cotangent squares, no area penalty.
        @test reference_layout_energy(context,variables) == 4shear^2
        @test reference_layout_energy(context,BigFloat[1,0,2,2]) == BigFloat(.01)*log(BigFloat(2))^2
        adjacent = (;context...,nodes=BigFloat[0 0; 1 0; 1 1; 0 1; 2 0; 2 1],
            quads=[(1,2,3,4),(2,5,6,3)],free=[1,4],active_faces=[1],owner=[1],neighbour=[2],
            initial_areas=BigFloat[1,1],anchors=BigFloat[1,1])
        @test isapprox(reference_layout_energy(adjacent,BigFloat[-1,-1,0,1]),
            (BigFloat(.25)+BigFloat(.01))*log(BigFloat(2))^2;rtol=big"1e-70")
        derivative = reference_directional_derivative(context,variables,BigFloat[1,1,0,0],BigFloat(10)^-20)
        @test isapprox(derivative,8shear;rtol=big"1e-40")
        shift = BigFloat[10^10 -10^9]
        translated = (;context...,nodes=context.nodes .+ shift)
        translated_variables = variables + BigFloat[10^10,10^10,-10^9,-10^9]
        @test reference_layout_energy(translated,translated_variables) == 4shear^2
        seam = (;context...,nodes=BigFloat[9.5 0; .5 0; .5 1; 9.5 1])
        @test reference_layout_energy(seam,BigFloat[.5,9.5,1,1]) == 0
        reversed = (;context...,quads=[(4,3,2,1)])
        @test isinf(reference_layout_energy(reversed,square))
        @test isinf(reference_layout_energy(context,BigFloat[.2,0,.2,1]))
    end
end

@testset "Gradient projections and term decomposition use original anchors" begin
    p = BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
        radial_clearance=.05,semicone_angle_deg=10.,eccentricity_ratio=.6)
    initial = ControlledSurfaceInlet.reference_seed(p,64,32;feed_diameter=5.,feed_z=50.)
    selection = ControlledSurfaceInlet.macro_data(initial)
    result = ControlledSurfaceInlet.controlled_layout(initial;max_iterations=4,trace_iterations=(1,4))
    setprecision(BigFloat,256) do
        context = reference_layout_context(initial,selection)
        @test eltype(context.nodes) == BigFloat
        for snapshot in result.trace
            exact = reference_layout_energy(context,BigFloat.(snapshot.variables))
            @test isapprox(exact,snapshot.energy;rtol=big"1e-12",atol=big"1e-12")
            terms = diagnostic_layout_terms(initial,selection,snapshot.variables)
            @test foldl(+,terms;init=0.) === snapshot.energy
            accumulated = diagnostic_accumulations(initial,selection,snapshot.variables,snapshot.energy)
            @test isapprox(accumulated.exact_term_sum,exact;rtol=big"1e-13",atol=big"1e-13")
            for (_,direction) in diagnosis_directions(snapshot)
                derivative = reference_directional_derivative(context,snapshot.variables,direction,BigFloat(10)^-20)
                @test isapprox(derivative,dot(snapshot.gradient,direction);rtol=big"1e-9",atol=big"1e-8")
            end
        end
        snapshot = last(result.trace)
        direction = last(first(diagnosis_directions(snapshot)))
        derivative256 = reference_directional_derivative(context,snapshot.variables,direction,BigFloat(10)^-20)
        setprecision(BigFloat,512) do
            higher = reference_layout_context(initial,selection)
            derivative512 = reference_directional_derivative(higher,snapshot.variables,direction,BigFloat(10)^-20)
            @test abs(derivative256-derivative512) < big"1e-40"
        end
    end
    smooth = smooth_diagnosis_direction(initial,selection)
    @test isapprox(norm(smooth),1.)
    @test all(iszero,smooth[length(selection.free)+1:end])
end

@testset "Diagnosis files and overwrite guard" begin
    mktempdir() do directory
        output = joinpath(directory,"diagnosis")
        result = diagnose_optimizer(output;budget=4)
        @test result.accepted_steps == 4
        @test TOML.parsefile(joinpath(output,"inputs.toml"))["precision_bits"] == 256
        @test TOML.parsefile(joinpath(output,"result.toml"))["termination_reason"] == "iteration_limit"
        @test deserialize(joinpath(output,"trace.bin")) == result.trace
        for name in ("derivatives","states","curvature-pairs","trials","curvature-probes")
            @test length(readlines(joinpath(output,"$name.tsv"))) > 1
        end
        before = sha256(read(joinpath(output,"result.toml")))
        @test_throws ArgumentError diagnose_optimizer(output;budget=4)
        @test sha256(read(joinpath(output,"result.toml"))) == before
        for budget in (0,1,2.)
            invalid = joinpath(directory,"invalid")
            @test_throws ArgumentError diagnose_optimizer(invalid;budget)
            @test !ispath(invalid)
        end
        invalid = joinpath(directory,"unknown")
        @test_throws ArgumentError diagnose_optimizer(invalid;name="unknown")
        @test !ispath(invalid)
    end
end
