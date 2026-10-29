include(joinpath(@__DIR__,"../../../../tests/xi2024_jfo_checks.jl"))
include(joinpath(ROOT,"cli.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;output=false,
        description="Check the four-edge JFO boundary adapter, continuum forcing and independent nodal algebra.") do cfg
        check_xi2024()
    end)
end
