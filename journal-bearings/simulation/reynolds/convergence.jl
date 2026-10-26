# Reproducible resolution sweep for any cylindrical unported full-film case;
# main_convergence() retains the original fixed illustrative example.
# Pressure is signed; force and flow diagnostics
# MUST NOT be interpreted as physical cavitating load or oil consumption.
using LinearAlgebra
include(joinpath(@__DIR__, "cylindrical.jl"))

function pressure_mismatch(coarse, fine, nt, nz)
    pc = reshape(coarse.pressure_pa, nt, nz)
    pf = reshape(fine.pressure_pa, 2nt, 2nz)
    fine_average = [(pf[2i-1, 2j-1] + pf[2i, 2j-1] +
                     pf[2i-1, 2j] + pf[2i, 2j]) / 4 for i in 1:nt, j in 1:nz]
    return norm(pc - fine_average) / norm(fine_average)
end

function signed_force_y(r, p, nt, nz)
    dtheta = 2pi / nt
    area = p.journal_midpoint_radius * 1e-3 * dtheta * (p.bearing_length * 1e-3 / nz)
    return -sum(r.pressure_pa[i + nt*(j-1)] * sin((i-0.5)*dtheta)
                for j in 1:nz for i in 1:nt) * area
end

function write_convergence(path, p, pairs; rpm, viscosity_pa_s)
    length(pairs) >= 2 || throw(ArgumentError("sweep requires at least two grid levels"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(("ntheta", "nz", "min_pressure_pa", "max_pressure_pa",
                          "signed_force_y_n", "absolute_axial_boundary_flux_m3_s",
                          "signed_net_axial_boundary_flux_m3_s", "relative_global_imbalance",
                          "relative_l1_cell_residual", "relative_max_cell_residual",
                          "relative_pressure_mismatch_from_previous"), '\t'))
        previous = nothing
        for (nt, nz) in pairs
            r = solve_cylindrical_full_film(p, nt, nz; rpm, viscosity_pa_s)
            d = r.diagnostics
            mismatch = isnothing(previous) ? NaN :
                pressure_mismatch(previous, r, nt ÷ 2, nz ÷ 2)
            println(io, join((nt, nz, minimum(r.pressure_pa), maximum(r.pressure_pa),
                              signed_force_y(r, p, nt, nz), d.absolute_side_flow_m3_s,
                              d.outward_side_flow_m3_s, d.relative_global_imbalance,
                              d.relative_l1_residual, d.relative_max_cell_residual,
                              mismatch), '\t'))
            previous = r
        end
    end
    return path
end

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
