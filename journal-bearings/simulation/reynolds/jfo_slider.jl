# Steady 1D rigid/isoviscous/incompressible JFO verification solver.
# This does NOT yet extend conical.jl, prescribe a bearing supply or solve EHD.
# See jfo-slider.md for source, flux, boundary and complementarity conventions.
include(joinpath(@__DIR__,"jfo_active_set.jl"))

"""
    assemble_jfo_slider(gap, n; length_m, viscosity_pa_s, entrainment_m_s,
                        end_pressure_pa, cavitation_pressure_pa,
                        inflow_fill_fraction=1.0)

SI units. gap(x) is the normal film thickness on x in [-length/2,length/2].
Entrainment is the AVERAGE of wall speeds, not their relative speed.
Pressure values share one reference; the unknown pressure is p-p_cav.
The explicit incoming fill applies at the upwind end only. Outgoing content
is transported from the last cell, not forced back to one at the outflow.
Returns A*p_above_cavitation + B*void = b, void=1-fill, per-unit-width fluxes.
"""
function assemble_jfo_slider(gap, n::Int; length_m::Real, viscosity_pa_s::Real,
                            entrainment_m_s::Real, end_pressure_pa::Tuple{Real,Real},
                            cavitation_pressure_pa::Real, inflow_fill_fraction::Real=1.0)
    n >= 2 || throw(ArgumentError("slider needs at least two cells"))
    L,mu,u,pc,fin = Float64.((length_m,viscosity_pa_s,entrainment_m_s,
                             cavitation_pressure_pa,inflow_fill_fraction))
    all(isfinite,(L,mu,u,pc,fin,end_pressure_pa...)) ||
        throw(ArgumentError("all slider inputs must be finite"))
    L > 0 && mu > 0 || throw(ArgumentError("length and viscosity must be positive"))
    !iszero(u) || throw(ArgumentError("zero entrainment needs a separate liquid-inventory model"))
    0 <= fin <= 1 || throw(ArgumentError("inflow fill must be in [0,1]"))
    ends = Float64.(end_pressure_pa) .- pc
    minimum(ends) >= 0 || throw(ArgumentError("end pressure cannot be below cavitation pressure"))
    ends[u > 0 ? 1 : 2] > 0 && fin != 1 &&
        throw(ArgumentError("pressurized inflow must be full: fill=1"))
    dx = L/n
    faces = collect(range(-L/2,L/2;length=n+1))
    centres = [(i-0.5)*dx-L/2 for i in 1:n]
    hf,hc = Float64.(gap.(faces)),Float64.(gap.(centres))
    valid_gap(h) = isfinite(h) && h > 0
    all(valid_gap,hf) && all(valid_gap,hc) ||
        throw(ArgumentError("sampled gaps must be finite and positive"))
    g = hf.^3 ./ (12mu*dx)
    g[1] *= 2; g[end] *= 2 # prescribed pressures are half a cell from centres
    adv = u .* hf
    ar,ac,av,br,bc,bv = Int[],Int[],Float64[],Int[],Int[],Float64[]
    b = zeros(n)
    for face in 1:n+1
        left,right = face-1,face
        upstream = u > 0 ? left : right
        internal_upstream = 1 <= upstream <= n
        constant = adv[face]*(internal_upstream ? 1.0 : fin)
        left == 0 && (constant += g[face]*ends[1])
        right == n+1 && (constant -= g[face]*ends[2])
        # One face flux, inserted with opposite signs into neighbouring cells.
        for (row,sign) in ((left,1.0),(right,-1.0))
            1 <= row <= n || continue
            b[row] -= sign*constant
            if left > 0
                push!(ar,row);push!(ac,left);push!(av,sign*g[face])
            end
            if right <= n
                push!(ar,row);push!(ac,right);push!(av,-sign*g[face])
            end
            if internal_upstream
                push!(br,row);push!(bc,upstream);push!(bv,-sign*adv[face])
            end
        end
    end
    return (; A=sparse(ar,ac,av,n,n),B=sparse(br,bc,bv,n,n),b,
             x_m=centres,face_x_m=faces,gap_m=hc,face_gap_m=hf,g,adv,
             end_above_cavitation_pa=ends,cavitation_pressure_pa=pc,
             entrainment_m_s=u,inflow_fill_fraction=fin,
             pressure_scale_pa=max(maximum(ends),12mu*abs(u)*L/minimum(hf)^2,1.0))
end

"""Solve the mixed pressure/void system by active-set column selection.
No post-solve pressure or liquid-content clipping is performed. A failed
active set, inadmissible content or continuity residual throws an error.
This steady solver is only verified for the documented 1D cases.
"""
function solve_jfo_slider(gap,n::Int; max_iterations::Int=2000,
                          tolerance::Real=1e-10,kwargs...)
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    isfinite(tolerance) && 0 < tolerance < 1 || throw(ArgumentError("invalid tolerance"))
    a = assemble_jfo_slider(gap,n;kwargs...)
    state = solve_jfo_active_set(a.A,a.B,a.b,a.pressure_scale_pa;max_iterations)
    p,v,iterations = state.pressure,state.void,state.iterations
    fill = 1 .- v
    # Reconstruct flux separately from the matrix residual. Per-unit-width m²/s.
    flux = zeros(n+1)
    for face in 1:n+1
        pl = face == 1 ? a.end_above_cavitation_pa[1] : p[face-1]
        pr = face == n+1 ? a.end_above_cavitation_pa[2] : p[face]
        upstream = a.entrainment_m_s > 0 ? face-1 : face
        fraction = 1 <= upstream <= n ? fill[upstream] : a.inflow_fill_fraction
        flux[face] = a.g[face]*(pl-pr)+a.adv[face]*fraction
    end
    scale = maximum(abs,a.adv) # nominal fully liquid Couette flux, even for dry inflow
    local_residual = maximum(abs,diff(flux))/scale
    complementarity = maximum(abs,(p ./ a.pressure_scale_pa).*v)
    boundary_imbalance = abs(flux[end]-flux[1])/scale
    flux_spread = (maximum(flux)-minimum(flux))/scale
    minimum(p) >= -tolerance*a.pressure_scale_pa &&
        minimum(fill) >= -tolerance && maximum(fill) <= 1+tolerance &&
        complementarity <= tolerance && local_residual <= tolerance &&
        boundary_imbalance <= tolerance && flux_spread <= tolerance ||
        error("pressure/content solution failed admissibility or continuity checks")
    diagnostics = (; iterations,pressure_scale_pa=a.pressure_scale_pa,
        minimum_pressure_above_cavitation_pa=minimum(p),minimum_fill=minimum(fill),
        maximum_fill=maximum(fill),complementarity_residual=complementarity,
        max_cell_imbalance_m2_s=maximum(abs,diff(flux)),
        relative_max_cell_imbalance=local_residual,
        relative_boundary_imbalance=boundary_imbalance,
        relative_flux_spread=flux_spread)
    return (; x_m=a.x_m,face_x_m=a.face_x_m,gap_m=a.gap_m,
             pressure_above_cavitation_pa=p,pressure_pa=p .+ a.cavitation_pressure_pa,
             fill_fraction=fill,face_flux_m2_s=flux,diagnostics)
end
