# Pressure-only obstacle on the shared unported tensor pressure assembly.
# Load bearing.jl and mesh/2d/grid.jl first; see ../docs/swift-stieber.md.
using LinearAlgebra, SparseArrays
isdefined(@__MODULE__, :assemble_full_film) || include(joinpath(@__DIR__,"full_film.jl"))
isdefined(@__MODULE__, :solve_jfo_active_set) || include(joinpath(@__DIR__,"jfo_active_set.jl"))

function pressure_obstacle(A,b,pressure_scale,flow_scale;max_iterations=2000)
    isempty(b) && throw(ArgumentError("obstacle requires free pressure cells"))
    isfinite(flow_scale) && flow_scale > 0 || throw(ArgumentError("invalid flow scale"))
    # A*P - reaction = b; P>=0, reaction>=0, P*reaction=0.
    # Reuse the algebraic kernel, NOT its liquid-content interpretation.
    state = solve_jfo_active_set(A,-spdiagm(0=>fill(flow_scale,length(b))),b,
                                 pressure_scale;max_iterations)
    P = state.pressure
    reaction = A*P-b
    minimum(P) >= -1e-10*pressure_scale && minimum(reaction) >= -1e-10*flow_scale &&
        maximum(abs,(P/pressure_scale).*(reaction/flow_scale)) < 1e-10 ||
        error("pressure-only obstacle failed complementarity")
    return (;pressure=P,reaction,iterations=state.iterations)
end

"""
    solve_swift_stieber(grid; rpm, viscosity_pa_s, end_pressure_pa=(0.,0.),
                       cavitation_pressure_pa=0., eccentricity_angle_rad=0.,
                       pinned_cells=Int[], max_iterations=2000)

Pressure-only Reynolds obstacle on a supplied unported tensor mesh. Pressures
use one common reference; end pressures cannot be below cavitation. Optional
cell indices prescribe p=p_cav, not a finite inlet geometry or liquid inventory.
The result reports fully-liquid extrapolated fluxes and their cavity defects,
NOT a transported cavity fill or validated oil consumption. No file IO occurs.
"""
function solve_swift_stieber(grid;rpm::Real,viscosity_pa_s::Real,
        end_pressure_pa::Tuple{Real,Real}=(0.,0.),cavitation_pressure_pa::Real=0.,
        eccentricity_angle_rad::Real=0.,pinned_cells::AbstractVector{Int}=Int[],
        max_iterations::Int=2000)
    pc=Float64(cavitation_pressure_pa)
    isfinite(pc) || throw(ArgumentError("cavitation pressure must be finite"))
    ends=Float64.(end_pressure_pa).-pc
    all(isfinite,ends) && minimum(ends)>=0 ||
        throw(ArgumentError("end pressures must be finite and not below cavitation"))
    a=assemble_full_film(grid;rpm,viscosity_pa_s,end_pressure_pa=ends,eccentricity_angle_rad)
    n=length(a.b);p=grid.parameters
    all(i->1<=i<=n,pinned_cells) || throw(ArgumentError("invalid pinned cell index"))
    fixed=falses(n);fixed[pinned_cells].=true
    free=findall(.!fixed)
    isempty(free) && throw(ArgumentError("pinned cells leave no free pressure unknowns"))
    pressure_scale=max(maximum(ends),12Float64(viscosity_pa_s)*abs(Float64(rpm)*2pi/60)*
        (p.journal_midpoint_radius*1e-3)*(2pi*p.housing_midpoint_radius*1e-3)/minimum(a.normal_gap_m)^2,1.)
    flow_scale=max(maximum(abs,a.adv),maximum(a.gz)*maximum(ends),1e-30)
    state=pressure_obstacle(a.A[free,free],a.b[free],pressure_scale,flow_scale;max_iterations)
    P=zeros(n);P[free]=state.pressure
    reaction=a.A*P-a.b
    flux=full_film_fluxes(a,P)
    positive=P .> 1e-10*pressure_scale
    residual=any(positive) ? maximum(abs,reaction[positive])/flow_scale : 0.
    complementarity=maximum(abs,(P[free]/pressure_scale).*(reaction[free]/flow_scale))
    reconstruction=maximum(abs,flux.balance-reaction)/flow_scale
    outward=(-sum(flux.qz[:,1]),sum(flux.qz[:,end]))
    cavity=sum(reaction[.!positive .& .!fixed]);exchange=sum(reaction[fixed])
    global_scale=max(sum(abs,outward)+abs(cavity)+abs(exchange),flow_scale)
    accounting=abs(sum(outward)-cavity-exchange)/global_scale
    residual<1e-10 && reconstruction<1e-10 && accounting<1e-10 ||
        error("Swift-Stieber failed pressure/flux-defect diagnostics")
    pressure_pa=P .+ pc
    diagnostics=(;iterations=state.iterations,pressure_scale_pa=pressure_scale,
        flow_scale_m3_s=flow_scale,positive_pressure_relative_residual=residual,
        complementarity_relative=complementarity,flux_reconstruction_relative=reconstruction,
        defect_accounting_relative=accounting,cavity_full_film_flux_defect_m3_s=cavity,
        pinned_line_exchange_m3_s=exchange,extrapolated_outward_end_flow_m3_s=outward,
        extrapolated_end_outflow_m3_s=sum(outward))
    return (;grid=a.grid,pressure_pa,pressure_above_cavitation_pa=P,
        normal_gap_m=a.normal_gap_m,reaction_m3_s=reaction,pinned_cells=fixed,
        circumferential_face_flow_m3_s=flux.qt,axial_face_flow_m3_s=flux.qz,
        housing_pressure_force_n=housing_pressure_force(a.grid,pressure_pa),diagnostics)
end
