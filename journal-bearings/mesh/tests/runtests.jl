# Each suite gets a fresh Julia process: source includes and test globals cannot
# leak between suites. Fluid experiments and snapshot comparisons are separate.
using Test

const suites = (
    "practical_acceptance", "inlet_telemetry", "meshing", "mesh_cli", "cylindrical",
    "inlet_counts", "inlet_geometry", "chart_geometry", "radial_geometry",
    "examples", "inlet_iterations", "inlet_diagnosis",
)

@testset "Surface mesher suites" begin
    for suite in suites
        @testset "$suite" begin
            project = normpath(joinpath(@__DIR__, "../.."))
            command = `$(Base.julia_cmd()) --project=$project --threads=1 --gcthreads=1 --check-bounds=yes $(joinpath(@__DIR__, "$suite.jl"))`
            command = addenv(command, "OPENBLAS_NUM_THREADS"=>"1", "OMP_NUM_THREADS"=>"1",
                "JULIA_NUM_THREADS"=>"1")
            @test success(run(ignorestatus(command)))
        end
    end
end
