# Fixed-input cylindrical resolution sweep.
include(joinpath(@__DIR__, "../reynolds/src/convergence.jl"))

const CYLINDER_RUN_DIR = joinpath(@__DIR__, "../../out/cases",
    "cylindrical-unported-e0p6", "runs", "2000rpm-vg32-40c")

function main_convergence()
    p = BearingParameters(journal_midpoint_radius=50.0, bearing_length=100.0,
                          radial_clearance=0.050, semicone_angle_deg=0.0,
                          eccentricity_ratio=0.6)
    pairs = ((32, 16), (64, 32), (128, 64), (256, 128), (512, 256))
    path = joinpath(CYLINDER_RUN_DIR, "convergence.tsv")
    write_convergence(path, p, pairs; rpm=2000.0, viscosity_pa_s=32e-6*860.0)
    println("Wrote $(abspath(path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_convergence()
end
