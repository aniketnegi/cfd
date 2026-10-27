# Steady rigid conical JFO solver. See conical-jfo.md.
# No elasticity, thermal solve, passage resistance, or load equilibrium.
include(joinpath(@__DIR__,"conical.jl"))
include(joinpath(@__DIR__,"jfo_active_set.jl"))

"""Reuse the checked cone pressure mobility and add conservative upwind void
transport on its periodic theta faces. P=p-p_cav; v=1-fill; A*P+B*v=b.
No axial Couette transport, hence no independent axial advective fill datum.
This low-level operator alone does not define a liquid inventory or supply.
"""
function assemble_conical_jfo(p::BearingParameters, nt::Int, nz::Int;
                              rpm::Real, viscosity_pa_s::Real,
                              end_pressure_pa::Tuple{Real,Real},
                              cavitation_pressure_pa::Real,
                              eccentricity_angle_rad::Real=0.0)
    pc = Float64(cavitation_pressure_pa)
    isfinite(pc) || throw(ArgumentError("cavitation pressure must be finite"))
    ends = Float64.(end_pressure_pa) .- pc
    minimum(ends) >= 0 || throw(ArgumentError("end pressure below cavitation pressure"))
    a = assemble_conical_full_film(p,nt,nz;rpm,viscosity_pa_s,end_pressure_pa=ends,eccentricity_angle_rad)
    rows,cols,values = Int[],Int[],Float64[]
    for j in 1:nz,i in 1:nt
        k,l = i+nt*(j-1),(i == nt ? 1 : i+1)+nt*(j-1)
        c = a.adv[i,j]
        upstream = c >= 0 ? k : l
        # Flux Qtheta = g*(Pleft-Pright) + adv*(1-v_upstream).
        push!(rows,k);push!(cols,upstream);push!(values,-c)
        push!(rows,l);push!(cols,upstream);push!(values,c)
    end
    scale = max(maximum(ends),12Float64(viscosity_pa_s)*abs(Float64(rpm)*2pi/60)*
        (p.journal_midpoint_radius*1e-3)*(2pi*p.housing_midpoint_radius*1e-3)/
        minimum(a.normal_gap_m)^2,1.0)
    return (; a..., B=sparse(rows,cols,values,nt*nz,nt*nz),
             cavitation_pressure_pa=pc,pressure_scale_pa=scale)
end

"""Reconstruct oriented liquid-equivalent face flows (m³/s), independently of
matrix multiplication. Circumferential flow includes the nt→1 periodic face;
axial flow is positive towards z=L. Balance is outward flow from each cell.
"""
function conical_jfo_fluxes(a,pressure,fill)
    nt,nz = size(a.gt)
    length(pressure) == length(fill) == nt*nz || throw(ArgumentError("field size mismatch"))
    P,F = reshape(pressure,nt,nz),reshape(fill,nt,nz)
    qt,qz = zeros(nt,nz),zeros(nt,nz+1)
    for j in 1:nz,i in 1:nt
        next = i == nt ? 1 : i+1
        upstream = a.adv[i,j] >= 0 ? i : next
        qt[i,j] = a.gt[i,j]*(P[i,j]-P[next,j])+a.adv[i,j]*F[upstream,j]
    end
    for j in 1:nz+1,i in 1:nt
        left = j == 1 ? a.end_pressure_pa[1] : P[i,j-1]
        right = j == nz+1 ? a.end_pressure_pa[2] : P[i,j]
        qz[i,j] = a.gz[i,j]*(left-right)
    end
    balance = [qt[i,j]-qt[i == 1 ? nt : i-1,j]+qz[i,j+1]-qz[i,j]
               for i in 1:nt,j in 1:nz]
    return (; qt,qz,balance=vec(balance))
end

"""
    solve_conical_jfo(p,nt,nz; rpm,viscosity_pa_s,end_pressure_pa,
                     cavitation_pressure_pa,feed_cells=nothing,feed_pressure_pa=nothing)

Requires a pressure supply above p_cav at an end or an explicitly supplied
cell mask. A feed cell is a fully liquid, bidirectional prescribed-pressure
reservoir, NOT a resolved hole, check valve or prescribed flow. Its required
source/sink is recovered from the face balance. No mask is selected implicitly.
All pressures use the caller's common reference. No underfilled periodic
inventory is silently selected when all reservoirs equal p_cav.
"""
function solve_conical_jfo(p::BearingParameters,nt::Int,nz::Int;
                           feed_cells::Union{Nothing,AbstractVector{Bool}}=nothing,
                           feed_pressure_pa::Union{Nothing,Real}=nothing,
                           max_iterations::Int=2000,tolerance::Real=1e-10,kwargs...)
    isfinite(tolerance) && 0 < tolerance < 1 || throw(ArgumentError("invalid tolerance"))
    a = assemble_conical_jfo(p,nt,nz;kwargs...)
    n = nt*nz
    isnothing(feed_cells) == isnothing(feed_pressure_pa) ||
        throw(ArgumentError("provide both feed_cells and feed_pressure_pa, or neither"))
    fixed = isnothing(feed_cells) ? falses(n) : BitVector(feed_cells)
    length(fixed) == n || throw(ArgumentError("feed mask must have one value per cell"))
    fp = 0.0
    if !isnothing(feed_cells)
        any(fixed) || throw(ArgumentError("empty feed mask; resolve or redefine the numerical supply"))
        all(fixed) && throw(ArgumentError("feed mask leaves no free film cells"))
        fp = Float64(feed_pressure_pa)-a.cavitation_pressure_pa
        isfinite(fp) && fp > 0 || throw(ArgumentError("feed must be above cavitation pressure and fully liquid"))
    end
    any(a.end_pressure_pa .> 0) || any(fixed) || throw(ArgumentError(
        "no pressure supply above p_cav: periodic liquid inventory is unspecified; ambient ends alone are not replenishment"))
    free = findall(.!fixed)
    P,F = zeros(n),ones(n)
    P[fixed] .= fp
    rhs = a.b[free]-a.A[free,findall(fixed)]*P[fixed]
    scale = max(a.pressure_scale_pa,fp)
    state = solve_jfo_active_set(a.A[free,free],a.B[free,free],rhs,scale;max_iterations)
    P[free] = state.pressure
    F[free] = 1 .- state.void
    flux = conical_jfo_fluxes(a,P,F)
    source = zeros(n)
    source[fixed] = flux.balance[fixed] # positive injection, negative return to reservoir
    residual = flux.balance-source
    outward = (-sum(flux.qz[:,1]),sum(flux.qz[:,end]))
    supply = sum(source)
    # Explicit normalization, including a pressure-driven zero-speed scale.
    reservoir_scale = max(maximum(a.end_pressure_pa),fp)
    flow_scale = max(maximum(abs,a.adv),maximum(abs,flux.qt),maximum(abs,flux.qz),
        max(maximum(a.gt),maximum(a.gz))*reservoir_scale)
    global_scale = max(sum(abs,source)+sum(abs,flux.qz[:,1])+sum(abs,flux.qz[:,end]),flow_scale)
    local_error = maximum(abs,residual)/flow_scale
    global_error = abs(sum(outward)-supply)/global_scale
    complementarity = maximum(abs,(P ./ scale).*(1 .- F))
    minimum(P) >= -tolerance*scale && minimum(F) >= -tolerance &&
        maximum(F) <= 1+tolerance && complementarity <= tolerance &&
        local_error <= tolerance && global_error <= tolerance ||
        error("conical JFO failed bounds, complementarity or source-aware continuity checks")
    diagnostics = (; iterations=state.iterations,pressure_scale_pa=scale,
        flow_scale_m3_s=flow_scale,global_flow_scale_m3_s=global_scale,
        max_cell_imbalance_m3_s=maximum(abs,residual),relative_max_cell_imbalance=local_error,
        relative_global_imbalance=global_error,complementarity_residual=complementarity,
        outward_end_flow_m3_s=outward,net_feed_flow_m3_s=supply,
        gross_feed_inflow_m3_s=sum(max.(source,0.0)),gross_feed_return_m3_s=-sum(min.(source,0.0)),
        minimum_fill=minimum(F),cavitated_cells=count(F .< 1-1e-8))
    return (; grid=a.grid,normal_gap_m=a.normal_gap_m,pressure_pa=P .+ a.cavitation_pressure_pa,
        pressure_above_cavitation_pa=P,fill_fraction=F,feed_cells=fixed,
        feed_source_m3_s=source,circumferential_face_flow_m3_s=flux.qt,
        axial_face_flow_m3_s=flux.qz,diagnostics)
end

function write_conical_jfo_vtk(path::AbstractString,result)
    ispath(path) && error("refusing to overwrite $path")
    write_vtk(path,result.grid;feed_cells=result.feed_cells)
    open(path,"a") do io
        for (name,data) in (("pressure_pa",result.pressure_pa),
                            ("pressure_above_cavitation_pa",result.pressure_above_cavitation_pa),
                            ("fill_fraction",result.fill_fraction),
                            ("normal_gap_m",result.normal_gap_m),
                            ("feed_source_m3_s",result.feed_source_m3_s))
            println(io,"SCALARS $name double 1\nLOOKUP_TABLE default")
            foreach(x->println(io,x),data)
        end
    end
    return path
end
