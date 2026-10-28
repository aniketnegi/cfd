# Signed full-film Reynolds flow on a supplied periodic tensor grid.
# Load bearing.jl and mesh/2d/grid.jl before this file. All calculations use SI.
using LinearAlgebra, SparseArrays

function film_radii(p::BearingParameters, z_m::Real)
    z_mm = z_m * 1e3
    return journal_radius(p, z_mm) * 1e-3, housing_radius(p, z_mm) * 1e-3
end

function normal_film_gap(p::BearingParameters, theta::Real, z_m::Real)
    radius, _ = film_radii(p, z_m)
    eccentricity, clearance = p.eccentricity * 1e-3, p.radial_clearance * 1e-3
    transverse = eccentricity * sin(theta)
    radial_gap = clearance - eccentricity*cos(theta) +
        transverse^2 / (radius + sqrt(radius^2 - transverse^2))
    return radial_gap * cosd(p.semicone_angle_deg)
end

function full_film_axes(grid)
    hasproperty(grid, :feed_edges) &&
        throw(ArgumentError("full-film tensor solver does not support fitted inlet meshes"))
    all(name -> hasproperty(grid, name), (:parameters, :theta_nodes, :z_nodes)) ||
        throw(ArgumentError("provide an unported mesh from housing_grid with tensor coordinates"))
    p = grid.parameters
    nt, nz = grid.n_theta, grid.n_axial
    nt >= 4 && nz >= 2 || throw(ArgumentError("grid must be at least 4×2"))
    length(grid.theta_nodes) == nt && length(grid.z_nodes) == nz+1 ||
        throw(ArgumentError("tensor coordinate counts do not match the mesh"))
    background = background_grid(grid.theta_nodes, grid.z_nodes)
    first(grid.z_nodes) == 0 && last(grid.z_nodes) == p.bearing_length ||
        throw(ArgumentError("tensor coordinates must span the bearing length"))
    grid.acceptance.accepted || throw(ArgumentError("mesh has not passed acceptance"))
    grid.quads == background.quads &&
        grid.nodes == [bearing_point(p, theta, z, 1.) for (theta,z) in background.nodes] ||
        throw(ArgumentError("mesh geometry/connectivity disagrees with its tensor coordinates"))
    length(grid.areas) == length(grid.radial_gaps) == nt*nz &&
        all(x -> isfinite(x) && x > 0, grid.areas) &&
        all(x -> isfinite(x) && x > 0, grid.radial_gaps) ||
        throw(ArgumentError("mesh areas and gaps must be finite and positive"))
    theta_edges = [grid.theta_nodes; first(grid.theta_nodes)+2pi]
    theta_widths = diff(theta_edges)
    all(width -> 0 < width < pi, theta_widths) ||
        throw(ArgumentError("each angular cell must span less than pi"))
    z_edges = grid.z_nodes .* 1e-3
    theta_centres = (theta_edges[1:end-1] + theta_edges[2:end]) ./ 2
    z_centres = (z_edges[1:end-1] + z_edges[2:end]) ./ 2
    return (; theta_edges, theta_widths, theta_centres, z_edges, z_centres,
        z_widths=diff(z_edges))
end

"""
    assemble_full_film(grid; rpm, viscosity_pa_s,
                       end_pressure_pa=(0., 0.), eccentricity_angle_rad=0.)

Conservative cell-centred Reynolds operator on a supplied `housing_grid` mesh.
Uniform and adjusted tensor spacing use the same assembly; zero cone angle is
cylindrical. Pressures are Pa in one reference, ordered (z=0, z=L). No cavitation
or feed is imposed. Returns `A*p=b` and the face coefficients for flux checks.
"""
function assemble_full_film(grid; rpm::Real, viscosity_pa_s::Real,
        end_pressure_pa::Tuple{Real,Real}=(0.,0.), eccentricity_angle_rad::Real=0.)
    speed, mu, angle = Float64.((rpm, viscosity_pa_s, eccentricity_angle_rad))
    ends = Float64.(end_pressure_pa)
    isfinite(speed) || throw(ArgumentError("rpm must be finite"))
    isfinite(mu) && mu > 0 || throw(ArgumentError("viscosity must be finite and positive"))
    isfinite(angle) || throw(ArgumentError("eccentricity angle must be finite"))
    all(isfinite, ends) || throw(ArgumentError("end pressures must be finite"))
    axes = full_film_axes(grid)
    p = grid.parameters
    nt, nz = grid.n_theta, grid.n_axial
    ca = cosd(p.semicone_angle_deg)
    omega = speed * 2pi / 60
    gap(theta,z) = normal_film_gap(p, theta-angle, z)
    rows, columns, values = Int[], Int[], Float64[]
    add(row,column,value) = (push!(rows,row); push!(columns,column); push!(values,value))
    id(i,j) = i + nt*(j-1)
    b = zeros(nt*nz)
    gt, adv, hcell = zeros(nt,nz), zeros(nt,nz), zeros(nt,nz)
    gz = zeros(nt,nz+1)

    for j in 1:nz, i in 1:nt
        next = mod1(i+1, nt)
        z = axes.z_centres[j]
        rj, rh = film_radii(p,z)
        hcell[i,j] = gap(axes.theta_centres[i],z)
        h = gap(axes.theta_edges[i+1],z)
        width = axes.z_widths[j] / ca
        distance = rh * (axes.theta_widths[i]+axes.theta_widths[next]) / 2
        conductance = gt[i,j] = h^3 * width / (12mu * distance)
        adv[i,j] = omega*rj * h * width / 2
        owner, neighbour = id(i,j), id(next,j)
        add(owner,owner,conductance); add(neighbour,neighbour,conductance)
        add(owner,neighbour,-conductance); add(neighbour,owner,-conductance)
        b[owner] -= adv[i,j]; b[neighbour] += adv[i,j]
    end
    for j in 1:nz+1, i in 1:nt
        z = axes.z_edges[j]
        _, radius = film_radii(p,z)
        h = gap(axes.theta_centres[i],z)
        distance = if j == 1
            axes.z_widths[1] / (2ca)
        elseif j == nz+1
            axes.z_widths[end] / (2ca)
        else
            (axes.z_widths[j-1]+axes.z_widths[j]) / (2ca)
        end
        conductance = gz[i,j] = h^3 * radius * axes.theta_widths[i] / (12mu * distance)
        if j == 1 || j == nz+1
            owner = id(i, j == 1 ? 1 : nz)
            add(owner,owner,conductance)
            b[owner] += conductance*ends[j == 1 ? 1 : 2]
        else
            owner, neighbour = id(i,j-1), id(i,j)
            add(owner,owner,conductance); add(neighbour,neighbour,conductance)
            add(owner,neighbour,-conductance); add(neighbour,owner,-conductance)
        end
    end
    all(x -> isfinite(x) && x > 0, gt) &&
        all(x -> isfinite(x) && x > 0, gz) &&
        all(x -> isfinite(x) && x > 0, hcell) && all(isfinite, adv) && all(isfinite,b) ||
        throw(ArgumentError("nonfinite or nonpositive film coefficients"))
    normal_gap_m = vec(hcell)
    # Unannotated housing_grid gaps have their offset along +X.
    gap_angle = get(grid, :eccentricity_angle_rad, 0.)
    result_grid = angle == gap_angle ? grid :
        (; grid..., radial_gaps=normal_gap_m .* (1e3/ca), eccentricity_angle_rad=angle)
    return (; A=sparse(rows,columns,values,nt*nz,nt*nz), b, grid=result_grid, axes,
        normal_gap_m, gt, adv, gz, end_pressure_pa=ends,
        eccentricity_angle_rad=angle)
end

"""Shared oriented face flows. Theta faces point towards increasing theta;
axial faces point towards +z, including the two prescribed-pressure boundaries.
"""
function full_film_fluxes(a, pressure_pa)
    nt, nz = size(a.gt)
    length(pressure_pa) == nt*nz || throw(ArgumentError("pressure size does not match mesh"))
    pressure = reshape(pressure_pa,nt,nz)
    qt, qz = zeros(nt,nz), zeros(nt,nz+1)
    for j in 1:nz, i in 1:nt
        qt[i,j] = a.adv[i,j] + a.gt[i,j]*(pressure[i,j]-pressure[mod1(i+1,nt),j])
    end
    for j in 1:nz+1, i in 1:nt
        left = j == 1 ? a.end_pressure_pa[1] : pressure[i,j-1]
        right = j == nz+1 ? a.end_pressure_pa[2] : pressure[i,j]
        qz[i,j] = a.gz[i,j]*(left-right)
    end
    balance = [qt[i,j]-qt[mod1(i-1,nt),j]+qz[i,j+1]-qz[i,j] for i in 1:nt,j in 1:nz]
    return (; qt, qz, balance=vec(balance))
end

# Piecewise-constant pressure; exact oriented cone area for each tensor cell.
function tensor_pressure_force(p, pressure, theta_edges, z_edges)
    nt, nz = length(theta_edges)-1, length(z_edges)-1
    length(pressure) == nt*nz || throw(ArgumentError("pressure size does not match mesh"))
    all(isfinite, pressure) || throw(ArgumentError("pressure must be finite"))
    force = zeros(3)
    for j in 1:nz, i in 1:nt
        z = (z_edges[j]+z_edges[j+1])/2
        _, radius = film_radii(p,z)
        low, high = theta_edges[i], theta_edges[i+1]
        weight = pressure[i+nt*(j-1)] * radius * (z_edges[j+1]-z_edges[j])
        force[1] += weight * (sin(high)-sin(low))
        force[2] += weight * (cos(low)-cos(high))
        force[3] += weight * (high-low) * p.slope
    end
    return Tuple(force)
end

"""Pressure force on the wetted housing sidewall, N. Excludes end faces and shear;
this is not a journal-force or load-equilibrium calculation.
"""
function housing_pressure_force(grid, pressure_pa)
    axes = full_film_axes(grid)
    return tensor_pressure_force(grid.parameters,pressure_pa,axes.theta_edges,axes.z_edges)
end

"""Solve signed full-film Reynolds on a supplied unported cylindrical/conical
mesh. Negative pressures are retained. The supplied mesh is not replaced;
no file IO occurs.
"""
function solve_full_film(grid; kwargs...)
    a = assemble_full_film(grid; kwargs...)
    pressure_pa = a.A \ a.b
    all(isfinite, pressure_pa) || error("nonfinite full-film pressure solution")
    flux = full_film_fluxes(a,pressure_pa)
    outward = (-sum(flux.qz[:,1]), sum(flux.qz[:,end]))
    boundary_scale = sum(abs,flux.qz[:,1]) + sum(abs,flux.qz[:,end])
    flux_scale = sum(abs,flux.qt) + sum(abs,flux.qz)
    advective_scale = sum(abs,a.adv)
    ratio(value,scale) = iszero(scale) ? (iszero(value) ? 0. : Inf) : value/scale
    max_imbalance = maximum(abs,flux.balance)
    residual = a.A*pressure_pa-a.b
    cells = length(pressure_pa)
    diagnostics = (;
        algebraic_residual_inf_m3_s=norm(residual,Inf),
        max_cell_imbalance_m3_s=max_imbalance,
        relative_l1_cell_imbalance=ratio(sum(abs,flux.balance),flux_scale),
        relative_max_cell_imbalance=ratio(max_imbalance,flux_scale/cells),
        outward_end_flow_m3_s=outward, absolute_end_flow_m3_s=boundary_scale,
        relative_global_imbalance=ratio(abs(sum(outward)),boundary_scale),
        minimum_normal_gap_m=minimum(a.normal_gap_m),
        max_cell_residual_m3_s=norm(residual,Inf), l1_cell_residual_m3_s=norm(residual,1),
        relative_max_cell_residual=norm(residual,Inf)/max(advective_scale/cells,eps()),
        relative_l1_residual=norm(residual,1)/max(advective_scale,eps()),
        outward_side_flow_m3_s=sum(outward), absolute_side_flow_m3_s=boundary_scale)
    return (; grid=a.grid, pressure_pa, normal_gap_m=a.normal_gap_m,
        circumferential_face_flow_m3_s=flux.qt, axial_face_flow_m3_s=flux.qz,
        housing_pressure_force_n=tensor_pressure_force(grid.parameters,pressure_pa,
            a.axes.theta_edges,a.axes.z_edges), diagnostics)
end
