# Floating, nonrotating bush; steady JFO; full overlap on a longer cone.
# Rigid translations only. All displacement/load API inputs are SI.
include(joinpath(@__DIR__,"conical_jfo.jl"))

function floating_bush_pose(reference::BearingParameters, displacement_m)
    length(displacement_m)==3 && all(isfinite,displacement_m) ||
        throw(ArgumentError("expected three finite bush translations in metres"))
    bx,by,bz=Float64.(displacement_m).*1e3
    # Housing-local z corresponds to shaft z+bz. Its inner cone is unchanged.
    clearance=reference.radial_clearance+bz*reference.slope
    clearance>0 || throw(ArgumentError("axial seating closes the reference clearance"))
    e=hypot(bx,by)
    e<clearance || throw(ArgumentError("journal contacts/intersects the translated bush"))
    p=BearingParameters(journal_midpoint_radius=reference.journal_midpoint_radius-bz*reference.slope,
        bearing_length=reference.bearing_length,radial_clearance=clearance,
        semicone_angle_deg=reference.semicone_angle_deg,eccentricity_ratio=e/clearance)
    angle=iszero(e) ? 0.0 : atan(-by,-bx) # journal offset relative to the bush
    return (;parameters=p,eccentricity_angle_rad=angle,
        displacement_m=Tuple(Float64.(displacement_m)),radial_clearance_m=clearance*1e-3,
        minimum_normal_gap_m=(clearance-e)*cosd(p.semicone_angle_deg)*1e-3,
        required_large_end_extension_m=max(-bz,0.0)*1e-3,
        required_small_end_extension_m=max(bz,0.0)*1e-3)
end

"""Housing-fixed +Y drill footprint; geometry only, never an implicit supply BC."""
function floating_bush_feed(reference,nt,nz;diameter_m=0.005,z_m=0.05)
    isfinite(diameter_m) && diameter_m>0 && isfinite(z_m) || throw(ArgumentError("invalid feed geometry"))
    L=reference.bearing_length*1e-3
    diameter_m/2<z_m<L-diameter_m/2 || throw(ArgumentError("feed must lie between axial ends"))
    mask=vec([begin
        theta=(i-0.5)*2pi/nt; z=(j-0.5)*L/nz
        _,rh=conical_radii(reference,z)
        rh*sin(theta)>0 && (rh*cos(theta))^2+(z-z_m)^2<=(diameter_m/2)^2
    end for i in 1:nt,j in 1:nz])
    any(mask) || throw(ArgumentError("feed has no cells; refine grid"))
    return BitVector(mask)
end

function floating_bush_state(reference,displacement_m,nt,nz;feed_cells,
                             rpm=2000.0,viscosity_pa_s=0.0277,feed_pressure_pa=5e5,
                             end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0,
                             jfo_max_iterations=2000)
    pose=floating_bush_pose(reference,displacement_m)
    p=pose.parameters
    film=solve_conical_jfo(p,nt,nz;feed_cells,rpm,viscosity_pa_s,feed_pressure_pa,
        end_pressure_pa,cavitation_pressure_pa,eccentricity_angle_rad=pose.eccentricity_angle_rad,
        max_iterations=jfo_max_iterations)
    # No solid bush wall at the aperture. External NET loads must include any
    # unmodelled feed-passage/pipe, weight or fixture forces. No pressure is
    # silently applied to a fictitious wall across the inlet opening.
    nominal=conical_housing_pressure_force(p,film.pressure_pa,nt,nz)
    land_pressure=copy(film.pressure_pa); land_pressure[feed_cells].=0.0
    force=conical_housing_pressure_force(p,land_pressure,nt,nz)
    return (;pose,film,bush_force_n=force,nominal_surface_force_n=nominal,
        omitted_aperture_force_n=Tuple(collect(nominal).-collect(force)))
end

"""Damped finite-difference Newton solve for a rigid bush's steady force balance.
fixed_axial=true checks transverse equilibrium at the seed's axial seating.
Failure is reported, not converted into a claim of nonexistence or stability.
"""
function equilibrate_floating_bush(reference,nt,nz;feed_cells,
        external_force_n=(0.0,1000.0,-200.0),initial_displacement_m=(-2.5e-6,2.5e-6,0.0),
        fixed_axial=false,force_tolerance_n=0.05,max_iterations=40,fd_step=1e-5,kwargs...)
    isfinite(force_tolerance_n) && force_tolerance_n>0 && isfinite(fd_step) && fd_step>0 ||
        throw(ArgumentError("invalid equilibrium tolerances"))
    max_iterations>=1 || throw(ArgumentError("max_iterations must be positive"))
    length(external_force_n)==3 && all(isfinite,external_force_n) || throw(ArgumentError("invalid external load"))
    reference.slope>0 || throw(ArgumentError("axial pose scaling requires a nonzero cone angle"))
    length(initial_displacement_m)==3 && all(isfinite,initial_displacement_m) ||
        throw(ArgumentError("invalid initial bush displacement"))
    c=reference.radial_clearance*1e-3
    lengths=[c,c,c/reference.slope]
    q=Float64.(collect(initial_displacement_m))./lengths
    axes=fixed_axial ? [1,2] : [1,2,3]
    external=Float64.(collect(external_force_n))
    scales=[max(hypot(external[1],external[2]),1.0),max(hypot(external[1],external[2]),1.0),max(abs(external[3]),1.0)]
    history=NamedTuple[]; failures=String[]; evaluations=Ref(0)
    function evaluate(q)
        evaluations[]+=1
        state=floating_bush_state(reference,q.*lengths,nt,nz;feed_cells,kwargs...)
        residual=collect(state.bush_force_n)+external
        return (;state,residual,scaled=residual[axes]./scales[axes])
    end
    function attempt(q)
        try
            return evaluate(q)
        catch err
            msg=sprint(showerror,err)
            if err isa ArgumentError || occursin("active set failed",msg) || occursin("conical JFO failed",msg)
                push!(failures,msg);return nothing
            end
            rethrow()
        end
    end
    current=evaluate(q)
    reason="iteration limit"
    converged=false
    for iteration in 1:max_iterations
        push!(history,(;iteration,displacement_m=Tuple(q.*lengths),force_n=current.state.bush_force_n,
            residual_n=Tuple(current.residual),scaled_norm=norm(current.scaled)))
        if maximum(abs,current.residual[axes])<=force_tolerance_n
            converged=true;reason="force tolerance";break
        end
        J=zeros(length(axes),length(axes))
        valid=true
        for (column,k) in enumerate(axes)
            step=fd_step*max(1.0,abs(q[k])); trial=copy(q);trial[k]+=step
            perturbed=attempt(trial)
            if isnothing(perturbed)
                step=-step;trial=copy(q);trial[k]+=step;perturbed=attempt(trial)
            end
            if isnothing(perturbed)
                valid=false;break
            end
            J[:,column]=(perturbed.scaled-current.scaled)/step
        end
        if !valid
            reason="no admissible finite-difference perturbation";break
        end
        # Small (2x2/3x3) SVD guards nearly singular force-position maps.
        S=svd(J)
        if minimum(S.S)<1e-12*maximum(S.S) || iszero(maximum(S.S))
            reason="singular force-position Jacobian";break
        end
        step=-(S.V*((S.U'*current.scaled)./S.S))
        limit=max(1.0,norm(q[axes]))
        norm(step)>limit && (step .*= limit/norm(step))
        accepted=false
        for backtrack in 0:18
            fraction=2.0^(-backtrack)
            trial=copy(q);trial[axes].+=fraction.*step
            candidate=attempt(trial)
            if !isnothing(candidate) && norm(candidate.scaled)<(1-1e-4*fraction)*norm(current.scaled)
                q=trial;current=candidate;accepted=true;break
            end
        end
        if !accepted
            reason="line search stalled";break
        end
    end
    # The last accepted state may be the step taken on the final iteration.
    converged=maximum(abs,current.residual[axes])<=force_tolerance_n
    converged && (reason="force tolerance")
    return (;converged,reason,state=current.state,residual_n=Tuple(current.residual),history,
        evaluations=evaluations[],trial_failures=failures,fixed_axial,
        required_axial_support_n=fixed_axial ? -current.residual[3] : (converged ? 0.0 : NaN))
end
