# Conservative liquid transport on the shared oriented pressure faces.
function jfo_content_operator(faces, cells)
    rows,columns,values = Int[],Int[],Float64[]
    for f in eachindex(faces.owner)
        k,l,c = faces.owner[f],faces.neighbour[f],faces.advective[f]
        upstream = c >= 0 ? k : l
        # An incoming boundary carries the fully liquid reservoir value.
        upstream == 0 && continue
        push!(rows,f);push!(columns,upstream);push!(values,c)
    end
    return sparse(rows,columns,values,length(faces.owner),cells)
end

"""Liquid-equivalent m³/s, oriented out of each face owner. Boundary inflow
carries fill=1; outgoing Couette transport carries the adjacent film content.
Pressure mobility includes the existing implicit nonorthogonal correction.
"""
function jfo_face_fluxes(a,pressure,fill)
    n=length(a.b)
    length(pressure)==length(fill)==n || throw(ArgumentError("field size mismatch"))
    all(isfinite,pressure) && all(isfinite,fill) || throw(ArgumentError("nonfinite JFO field"))
    faces=a.faces
    flow=faces.pressure_flux*pressure+faces.constant-faces.advective
    for f in eachindex(flow)
        k,l,c=faces.owner[f],faces.neighbour[f],faces.advective[f]
        upstream=c >= 0 ? k : l
        flow[f] += c*(upstream==0 ? 1. : fill[upstream])
    end
    return (;face_flow=flow,balance=faces.divergence*flow,
        boundary_outflow=boundary_outflows(faces,flow))
end

function assemble_fitted_jfo(grid;rpm,viscosity_pa_s,end_pressure_pa,
        cavitation_pressure_pa,eccentricity_angle_rad,feed_pressure_pa,gap_correction_m=nothing)
    pc=Float64(cavitation_pressure_pa)
    isfinite(pc) || throw(ArgumentError("cavitation pressure must be finite"))
    ends=Float64.(end_pressure_pa).-pc
    all(isfinite,ends) && minimum(ends)>=0 || throw(ArgumentError("invalid end pressures"))
    isnothing(feed_pressure_pa) && throw(ArgumentError("fitted JFO requires feed_pressure_pa"))
    feed=Float64(feed_pressure_pa)-pc
    isfinite(feed) && feed>0 || throw(ArgumentError("feed must be above cavitation pressure and fully liquid"))
    a=assemble_full_film(grid;rpm,viscosity_pa_s,end_pressure_pa=ends,
        eccentricity_angle_rad,feed_pressure_pa=feed,gap_correction_m)
    any(==(:feed),a.faces.boundary) || throw(ArgumentError("fitted JFO requires actual inlet rim faces"))
    C=jfo_content_operator(a.faces,length(a.b))
    p=a.grid.parameters
    scale=max(maximum(ends),feed,12Float64(viscosity_pa_s)*abs(Float64(rpm)*2pi/60)*
        (p.journal_midpoint_radius*1e-3)*(2pi*p.housing_midpoint_radius*1e-3)/minimum(a.normal_gap_m)^2,1.)
    return (;a...,B=-a.faces.divergence*C,content_flux=C,
        cavitation_pressure_pa=pc,pressure_scale_pa=scale)
end

function solve_fitted_jfo(grid;feed_pressure_pa,max_iterations,tolerance,kwargs...)
    a=assemble_conical_jfo(grid;feed_pressure_pa,kwargs...)
    state=solve_jfo_active_set(a.A,a.B,a.b,a.pressure_scale_pa;max_iterations)
    P,F=state.pressure,1 .-state.void
    flux=jfo_face_fluxes(a,P,F)
    faces=a.faces
    boundary=faces.neighbour.==0
    feed=faces.boundary.==:feed
    flow_scale=max(maximum(abs,faces.advective),maximum(abs,faces.constant),
        maximum(abs,flux.face_flow),1e-30)
    global_scale=max(sum(abs,flux.face_flow[boundary]),flow_scale)
    local_error=maximum(abs,flux.balance)/flow_scale
    global_error=abs(sum(values(flux.boundary_outflow)))/global_scale
    reconstruction=maximum(abs,flux.balance-(a.A*P+a.B*(1 .-F)-a.b))/flow_scale
    complementarity=maximum(abs,(P/a.pressure_scale_pa).*(1 .-F))
    minimum(P)>=-tolerance*a.pressure_scale_pa && minimum(F)>=-tolerance &&
        maximum(F)<=1+tolerance && complementarity<=tolerance &&
        local_error<=tolerance && global_error<=tolerance && reconstruction<=tolerance ||
        error("fitted JFO failed bounds, complementarity or face continuity checks")
    outward=(flux.boundary_outflow.large,flux.boundary_outflow.small)
    diagnostics=(;iterations=state.iterations,pressure_scale_pa=a.pressure_scale_pa,
        flow_scale_m3_s=flow_scale,global_flow_scale_m3_s=global_scale,
        max_cell_imbalance_m3_s=maximum(abs,flux.balance),relative_max_cell_imbalance=local_error,
        relative_global_imbalance=global_error,flux_reconstruction_relative=reconstruction,
        complementarity_residual=complementarity,outward_end_flow_m3_s=outward,
        outward_feed_flow_m3_s=flux.boundary_outflow.feed,
        net_feed_flow_m3_s=-flux.boundary_outflow.feed,
        gross_feed_inflow_m3_s=-sum(min.(flux.face_flow[feed],0.)),
        gross_feed_return_m3_s=sum(max.(flux.face_flow[feed],0.)),
        gross_boundary_inflow_m3_s=-sum(min.(flux.face_flow[boundary],0.)),
        gross_boundary_outflow_m3_s=sum(max.(flux.face_flow[boundary],0.)),
        minimum_fill=minimum(F),cavitated_cells=count(F.<1-1e-8))
    return (;grid=a.grid,normal_gap_m=a.normal_gap_m,pressure_pa=P.+a.cavitation_pressure_pa,
        pressure_above_cavitation_pa=P,fill_fraction=F,face_flow_m3_s=flux.face_flow,
        boundary_outflow_m3_s=flux.boundary_outflow,
        housing_pressure_force_n=housing_pressure_force(a.grid,P.+a.cavitation_pressure_pa),diagnostics)
end
