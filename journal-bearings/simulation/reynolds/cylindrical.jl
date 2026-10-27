# Steady, isoviscous, incompressible Newtonian Reynolds equation on a
# *cylindrical* housing surface (semi-cone angle = 0). Signed full film only:
# negative gauge pressures here are numerical verification data, NOT cavitation.
# The θ-z finite-volume fluxes follow San Andrés, Notes 4, Eq. (4.1),
# https://rotorlab.tamu.edu/me626/Notes_pdf/Notes04%20Static%20Load%20JBs%2010.pdf
# SI units internally; bearing.jl and housing_grid use mm.
using LinearAlgebra, SparseArrays
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__, "../../mesh/vtk.jl"))

# Output locations for the illustrative main() run only. The solver function
# below is reusable and does not write files or impose these operating values.
const CYLINDER_CASE_DIR = joinpath(@__DIR__, "..", "..", "out", "cases",
                                   "cylindrical-unported-e0p6")
const CYLINDER_RUN_DIR = joinpath(CYLINDER_CASE_DIR, "runs", "2000rpm-vg32-40c")

"""
    solve_cylindrical_full_film(p, ntheta, nz; rpm, viscosity_pa_s)

Solve ∇·q=0, qθ=Uh/2-h³ ∂θp/(12μR), qz=-h³ ∂zp/(12μ),
U=ωRj, for prescribed journal position. θ is periodic; p=0 gauge at z=0,L.
The geometry is sampled from the same mapping as the exported 2D grid.
Face fluxes are shared by neighbouring cells (local conservation).
"""
function solve_cylindrical_full_film(p::BearingParameters, ntheta::Int, nz::Int;
                                     rpm::Real, viscosity_pa_s::Real)
    iszero(p.semicone_angle_deg) ||
        throw(ArgumentError("this verified formulation is cylindrical only"))
    ntheta >= 4 && nz >= 2 || throw(ArgumentError("grid must be at least 4×2"))
    speed, mu = Float64(rpm), Float64(viscosity_pa_s)
    isfinite(speed) || throw(ArgumentError("rpm must be finite"))
    isfinite(mu) && mu > 0 || throw(ArgumentError("viscosity must be positive and finite"))
    grid = housing_grid(p, ntheta, nz)
    dtheta = 2pi / ntheta
    dz = p.bearing_length * 1e-3 / nz
    R = p.housing_midpoint_radius * 1e-3
    U = (speed * 2pi / 60) * p.journal_midpoint_radius * 1e-3
    nc = ntheta * nz
    id(i, j) = i + ntheta * (j - 1)
    # Evaluate the clearance on faces rather than averaging h³ across cells.
    # Algebraic rearrangement of bearing.jl's radial intersection, avoiding
    # cancellation of two ~50 mm radii (and angle-dependent round-off at e=0).
    function h(theta)
        e, rj = p.eccentricity, p.journal_midpoint_radius
        t = e * sin(theta)
        return (p.radial_clearance - e*cos(theta) + t^2 / (rj + sqrt(rj^2 - t^2))) * 1e-3
    end
    rows, cols, vals = Int[], Int[], Float64[]
    b = zeros(nc)
    add(k, l, v) = (push!(rows, k); push!(cols, l); push!(vals, v))
    # Each circumferential face is visited once; west of i=1 is east of
    # i=ntheta. The signed advective flux is also inserted as an equal and
    # opposite pair, so the periodic seam cannot create or destroy mass.
    for j in 1:nz, i in 1:ntheta
        k, next = id(i, j), id(i == ntheta ? 1 : i + 1, j)
        hf = h(i * dtheta)
        conductance = hf^3 * dz / (12mu * R * dtheta)
        advective = U * hf * dz / 2
        add(k, k, conductance); add(next, next, conductance)
        add(k, next, -conductance); add(next, k, -conductance)
        b[k] -= advective
        b[next] += advective
    end
    # Axial faces have width R*dtheta. Boundary pressure is 0 gauge and
    # boundary-centre distance is dz/2, not dz.
    width = R * dtheta
    for j in 1:nz, i in 1:ntheta
        k = id(i, j)
        hf = h((i - 0.5) * dtheta)
        if j < nz
            next = id(i, j + 1)
            g = hf^3 * width / (12mu * dz)
            add(k, k, g); add(next, next, g)
            add(k, next, -g); add(next, k, -g)
        end
        if j == 1 || j == nz
            g = hf^3 * width / (12mu * (dz / 2))
            add(k, k, g)
        end
    end
    A = sparse(rows, cols, vals, nc, nc)
    pressure_pa = A \ b                # Julia SparseArrays/SuiteSparse, not a custom linear solver
    residual = A * pressure_pa - b
    # Signed outward side flux: signed full-film pressure can make positive
    # and negative contributions cancel. Report absolute flux as well.
    outward = 0.0
    absolute_outward = 0.0
    for j in (1, nz), i in 1:ntheta
        k = id(i, j)
        g = h((i - 0.5) * dtheta)^3 * width / (12mu * (dz / 2))
        f = g * pressure_pa[k]
        outward += f
        absolute_outward += abs(f)
    end
    # Dimensionless residual scale based on total *magnitude* of the
    # circumferential advective fluxes (not on signed net outflow).
    advective_scale = nz * sum(i -> abs(U * h(i * dtheta) * dz / 2), 1:ntheta)
    diagnostics = (; max_cell_residual_m3_s=norm(residual, Inf),
                   relative_max_cell_residual=norm(residual, Inf) /
                       max(advective_scale / nc, eps()),
                   l1_cell_residual_m3_s=norm(residual, 1),
                   outward_side_flow_m3_s=outward,
                   absolute_side_flow_m3_s=absolute_outward,
                   relative_global_imbalance=abs(outward) / max(absolute_outward, eps()),
                   relative_l1_residual=norm(residual, 1) / max(advective_scale, eps()))
    return (; grid, pressure_pa, diagnostics)
end

function write_pressure_vtk(path::AbstractString, result)
    # The shared mesh writer exports the identical geometry and cell order.
    write_vtk(path, result.grid)
    open(path, "a") do io
        println(io, "SCALARS pressure_gauge_pa double 1\nLOOKUP_TABLE default")
        for value in result.pressure_pa
            println(io, value)
        end
    end
end

function main()
    # Chosen reference: nominal ISO VG 32 at 40 °C is ν=32 mm²/s.
    # ρ=860 kg/m³ is a separate illustrative assumption from Gangrade,
    # *not* supplied by the ISO viscosity grade; μ=νρ=0.02752 Pa·s.
    p = BearingParameters(journal_midpoint_radius=50.0, bearing_length=100.0,
                          radial_clearance=0.050, semicone_angle_deg=0.0,
                          eccentricity_ratio=0.6)
    result = solve_cylindrical_full_film(p, 512, 256;
                                        rpm=2000.0, viscosity_pa_s=32e-6 * 860.0)
    mesh_path = joinpath(CYLINDER_CASE_DIR, "mesh.vtk")
    path = joinpath(CYLINDER_RUN_DIR, "full-film.vtk")
    write_vtk(mesh_path, result.grid)
    write_pressure_vtk(path, result)
    println("Wrote signed full-film gauge pressure to $(abspath(path))")
    println("Pressure: $(minimum(result.pressure_pa)) to $(maximum(result.pressure_pa)) Pa")
    println("Conservation diagnostics: $(result.diagnostics)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
