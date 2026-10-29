include(joinpath(@__DIR__,"../../../../cli.jl"))
include(joinpath(@__DIR__,"../../../tests/fitted_jfo_checks.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;output=false,description="Check tensor and fitted-rim conical JFO transport.") do cfg
        check_tensor_jfo_contract()
        check_fitted_jfo()
    end)
end
