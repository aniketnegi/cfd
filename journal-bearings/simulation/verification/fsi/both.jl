# Same study interface and physics, with both approved elastic solids enabled.
include(joinpath(@__DIR__,"run.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,description="Run fixed-position cylinder/cone FSI with elastic journal and housing.") do cfg
        observer=stderr isa Base.TTY && !cfg.quiet && !cfg.verbose ? _->yield() : nothing
        run_fsi(cfg.output;full=cfg.full,observer,both_elastic=true)
    end)
end
