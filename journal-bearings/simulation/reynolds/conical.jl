# Rigid, steady, Newtonian, isoviscous conical Reynolds verification ONLY.
# Signed full film: no cavitation, central feed, elasticity or load equilibrium.
# Derivation and approximation choices: conical-formulation.md.
# Geometry and VTK dependencies are loaded through cylindrical.jl.
include(joinpath(@__DIR__, "cylindrical.jl"))

# All helpers below use SI metres; BearingParameters stores millimetres.
function conical_radii(p::BearingParameters, z_m::Real)
    rj = p.journal_midpoint_radius * 1e-3 +
         (p.half_bearing_length * 1e-3 - z_m) * p.slope
    return rj, rj + p.radial_clearance * 1e-3
end

function conical_normal_gap(p::BearingParameters, theta::Real, z_m::Real)
    rj, _ = conical_radii(p, z_m)
    e, c = p.eccentricity * 1e-3, p.radial_clearance * 1e-3
    t = e * sin(theta)
    # Stable radial-ray intersection, as in the cylinder, then normal projection.
    # This is a thin-film normal-gap approximation, not an exact normal ray cast.
    radial = c - e*cos(theta) + t^2 / (rj + sqrt(rj^2 - t^2))
    return radial * cosd(p.semicone_angle_deg)
end

"""
    assemble_conical_full_film(p, nt, nz; rpm, viscosity_pa_s,
                              end_pressure_pa=(0.0, 0.0))

Finite-volume operator for periodic theta and prescribed axial-end pressures.
End order is (large end z=0, small end z=L), pressures in Pa gauge.
The exposed assembly also permits independent manufactured-solution tests;
it does not accept physical feed sources or implement cavitation.
"""
function assemble_conical_full_film(p::BearingParameters, nt::Int, nz::Int;
                                   rpm::Real, viscosity_pa_s::Real,
                                   end_pressure_pa::Tuple{Real,Real}=(0.0, 0.0),
                                   eccentricity_angle_rad::Real=0.0)
    nt >= 4 && nz >= 2 || throw(ArgumentError("grid must be at least 4×2"))
    speed, mu = Float64(rpm), Float64(viscosity_pa_s)
    ends = Float64.(end_pressure_pa)
    isfinite(speed) || throw(ArgumentError("rpm must be finite"))
    isfinite(mu) && mu > 0 || throw(ArgumentError("viscosity must be finite and positive"))
    all(isfinite, ends) || throw(ArgumentError("end pressures must be finite"))
    angle = Float64(eccentricity_angle_rad)
    isfinite(angle) || throw(ArgumentError("eccentricity angle must be finite"))
    gap(theta,z) = conical_normal_gap(p,theta-angle,z)
    grid = housing_grid(p, nt, nz)  # also validates radial intersection and gaps
    # The housing grid stays fixed; only the journal offset direction changes.
    if !iszero(angle)
        gaps = [gap((i-0.5)*2pi/nt,(j-0.5)*p.bearing_length*1e-3/nz)*
                1e3/cosd(p.semicone_angle_deg) for i in 1:nt,j in 1:nz]
        grid = (;grid...,radial_gaps=vec(gaps))
    end
    dt = 2pi / nt
    dz = p.bearing_length * 1e-3 / nz
    ca = cosd(p.semicone_angle_deg)
    ds = dz / ca
    omega = speed * 2pi / 60
    id(i, j) = i + nt*(j-1)
    rows, cols, vals = Int[], Int[], Float64[]
    b = zeros(nt*nz)
    add(k, l, v) = (push!(rows, k); push!(cols, l); push!(vals, v))
    gt, adv = zeros(nt, nz), zeros(nt, nz)
    gz = zeros(nt, nz+1)
    hcell = zeros(nt, nz)
    for j in 1:nz, i in 1:nt
        z = (j-0.5)*dz
        rj, rh = conical_radii(p, z)
        hcell[i,j] = gap((i-0.5)*dt, z)
        h = gap(i*dt, z)
        # Circumferential faces: width ds, centre distance rh*dtheta.
        gt[i,j] = h^3 * ds / (12mu * rh * dt)
        adv[i,j] = omega*rj * h * ds / 2
        k, l = id(i,j), id(i == nt ? 1 : i+1,j)
        g = gt[i,j]
        add(k,k,g); add(l,l,g); add(k,l,-g); add(l,k,-g)
        b[k] -= adv[i,j]; b[l] += adv[i,j]
    end
    for j in 1:(nz+1), i in 1:nt
        z = (j-1)*dz
        _, rh = conical_radii(p, z)
        h = gap((i-0.5)*dt, z)
        distance = (j == 1 || j == nz+1) ? ds/2 : ds
        # Meridional faces: width rh*dtheta; no meridional Couette term.
        g = gz[i,j] = h^3 * rh * dt / (12mu * distance)
        if j == 1 || j == nz+1
            k = id(i, j == 1 ? 1 : nz)
            pend = ends[j == 1 ? 1 : 2]
            add(k,k,g); b[k] += g*pend
        else
            k, l = id(i,j-1), id(i,j)
            add(k,k,g); add(l,l,g); add(k,l,-g); add(l,k,-g)
        end
    end
    return (; A=sparse(rows,cols,vals,nt*nz,nt*nz), b, grid,
            normal_gap_m=vec(hcell), gt, adv, gz, end_pressure_pa=ends)
end

"""Pressure-only force ON THE HOUSING's conical wetted surface, in newtons.
Cell pressure is piecewise constant; oriented areas are integrated exactly.
This excludes end faces and shear, and is not a journal load-equilibrium solve.
"""
function conical_housing_pressure_force(p::BearingParameters, pressure, nt, nz)
    length(pressure) == nt*nz || throw(ArgumentError("pressure size does not match grid"))
    dt, dz = 2pi/nt, p.bearing_length*1e-3/nz
    force = zeros(3)
    for j in 1:nz, i in 1:nt
        _, rh = conical_radii(p, (j-0.5)*dz)
        lo, hi = (i-1)*dt, i*dt
        w = pressure[i+nt*(j-1)] * rh * dz
        force[1] += w * (sin(hi)-sin(lo))
        force[2] += w * (cos(lo)-cos(hi))
        force[3] += w * dt * p.slope
    end
    return Tuple(force)
end

function solve_conical_full_film(p::BearingParameters, nt::Int, nz::Int; kwargs...)
    a = assemble_conical_full_film(p, nt, nz; kwargs...)
    pressure_pa = a.A \ a.b
    pressure = reshape(pressure_pa, nt, nz)
    qt, qz = zeros(nt,nz), zeros(nt,nz+1)
    for j in 1:nz, i in 1:nt
        next = i == nt ? 1 : i+1
        qt[i,j] = a.adv[i,j] + a.gt[i,j]*(pressure[i,j]-pressure[next,j])
    end
    for j in 1:(nz+1), i in 1:nt
        left = j == 1 ? a.end_pressure_pa[1] : pressure[i,j-1]
        right = j == nz+1 ? a.end_pressure_pa[2] : pressure[i,j]
        qz[i,j] = a.gz[i,j]*(left-right) # positive toward +z (small end)
    end
    # Reconstruct continuity from face fluxes, not just the linear-system residual.
    balance = [qt[i,j]-qt[i == 1 ? nt : i-1,j]+qz[i,j+1]-qz[i,j]
               for i in 1:nt, j in 1:nz]
    outward = (-sum(qz[:,1]), sum(qz[:,end]))
    boundary_scale = sum(abs, qz[:,1]) + sum(abs, qz[:,end])
    flux_scale = sum(abs, qt) + sum(abs, qz)
    ratio(x, scale) = iszero(scale) ? (iszero(x) ? 0.0 : Inf) : x/scale
    diagnostics = (;
        algebraic_residual_inf_m3_s=norm(a.A*pressure_pa-a.b, Inf),
        max_cell_imbalance_m3_s=maximum(abs, balance),
        relative_l1_cell_imbalance=ratio(sum(abs,balance), flux_scale),
        relative_max_cell_imbalance=ratio(maximum(abs,balance), flux_scale/(nt*nz)),
        outward_end_flow_m3_s=outward,
        absolute_end_flow_m3_s=boundary_scale,
        relative_global_imbalance=ratio(abs(sum(outward)), boundary_scale),
        minimum_normal_gap_m=minimum(a.normal_gap_m))
    return (; grid=a.grid, pressure_pa, normal_gap_m=a.normal_gap_m,
            circumferential_face_flow_m3_s=qt, axial_face_flow_m3_s=qz,
            housing_pressure_force_n=conical_housing_pressure_force(p,pressure_pa,nt,nz),
            diagnostics)
end

function write_conical_pressure_vtk(path::AbstractString, result)
    write_pressure_vtk(path, result)
    open(path, "a") do io
        println(io, "SCALARS normal_gap_m double 1\nLOOKUP_TABLE default")
        for h in result.normal_gap_m
            println(io, h)
        end
    end
    return path
end
