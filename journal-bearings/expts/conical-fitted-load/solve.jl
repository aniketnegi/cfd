# Explicit fitted-case dispatch into the existing Newton equilibrium algorithm.
# No rectangular operator, feed mask, or fitted-cell reshaping is used here.
include(joinpath(@__DIR__,"../conical-fitted-inlet/run.jl"))
module PublicGeometry
include(joinpath(@__DIR__,"../../mesh/2d/conical/surface_feed.jl"))
end

struct FittedLoadProblem{M}
    reference::BearingParameters
    mesh::M
    # The unchanged Newton routine uses these two reference scaling quantities.
    radial_clearance::Float64
    slope::Float64
end

function fitted_load_problem(level)
    p=PublicGeometry.BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    grid=PublicGeometry.surface_feed_grid(p,128level,64level;feed_diameter=5.,feed_z=50.)
    ref=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    return FittedLoadProblem(ref,grid.surface,ref.radial_clearance,ref.slope)
end

function floating_bush_state(problem::FittedLoadProblem,d,nt,nz;
        feed_cells,jfo_max_iterations=150)
    feed_cells===nothing || throw(ArgumentError("fitted land problem has a rim, not a feed-cell mask"))
    (nt,nz)==(128problem.mesh.level,64problem.mesh.level) ||
        throw(ArgumentError("counts label the mesh background only; they do not index fitted cells"))
    pose=floating_bush_pose(problem.reference,d)
    # The housing mesh is fixed. Only the journal pose/gaps/speed radii change.
    film=solve_fitted(problem.mesh,pose;rpm=2000.,max_iterations=jfo_max_iterations)
    return (;pose,film,bush_force_n=Tuple(film.force))
end

function radial_equilibrium(problem,clearance_um,seed_xy;
        force_tolerance_n=.02,fd_step=1e-6,max_iterations=12,axial_load_n=200.)
    isfinite(axial_load_n) && axial_load_n>0 || throw(ArgumentError("positive axial load magnitude required"))
    isfinite(clearance_um) && clearance_um>0 || throw(ArgumentError("positive clearance required"))
    c=clearance_um*1e-6; bz=(c-problem.radial_clearance*1e-3)/problem.slope
    level=problem.mesh.level
    return equilibrate_floating_bush(problem,128level,64level;feed_cells=nothing,
        fixed_axial=true,initial_displacement_m=(seed_xy[1],seed_xy[2],bz),
        external_force_n=(0.,1000.,-Float64(axial_load_n)),force_tolerance_n,fd_step,max_iterations,
        jfo_max_iterations=150)
end

"""Bisection of a continuous bracketed scalar residual. The caller must ensure
that each function evaluation is a converged radial equilibrium. A sign bracket
is necessary here; a failed radial solve must never be returned as a residual.
"""
function bisect_axial(f,a,b;force_tolerance_n=.02,max_iterations=24)
    isfinite(a) && isfinite(b) && a<b || throw(ArgumentError("ordered finite interval required"))
    force_tolerance_n>0 && isfinite(force_tolerance_n) && max_iterations>0 ||
        throw(ArgumentError("invalid root tolerances"))
    evaluate(x)=begin y=f(x);isfinite(y) || error("nonfinite axial residual");y end
    fa=evaluate(a);fb=evaluate(b);history=NamedTuple[]
    abs(fa)<=force_tolerance_n && return (;converged=true,c=a,G=fa,history)
    abs(fb)<=force_tolerance_n && return (;converged=true,c=b,G=fb,history)
    signbit(fa)!=signbit(fb) || throw(ArgumentError("no sign bracket; nonexistence is not established"))
    c=a;fc=fa
    for iteration in 1:max_iterations
        c=(a+b)/2;fc=evaluate(c);push!(history,(;iteration,a,b,c,G=fc))
        abs(fc)<=force_tolerance_n && return (;converged=true,c,G=fc,history)
        if signbit(fc)==signbit(fa);a=c;fa=fc;else;b=c;fb=fc;end
    end
    return (;converged=false,c,G=fc,history)
end
