# Normal-gap, inertia, dissipation and wall-speed sensitivity diagnostics.
include(joinpath(@__DIR__,"..","..","simulation","reynolds","conical.jl"))
using TOML, SHA

function review_conical(out)
    p = BearingParameters(journal_midpoint_radius=50.0,bearing_length=100.0,
        radial_clearance=0.05,semicone_angle_deg=10.0,eccentricity_ratio=0.6)
    nt,nz,rpm,mu,rho = 256,128,2000.0,0.02752,860.0
    omega,ca,sa,e = rpm*2pi/60,cosd(10.0),sind(10.0),30e-6
    dt,dz = 2pi/nt,0.1/nz
    a = assemble_conical_full_film(p,nt,nz;rpm,viscosity_pa_s=mu)
    baseline = a.A \ a.b
    aspect,reduced_re,re_gap = Float64[],Float64[],Float64[]
    couette_dissipation = 0.0
    for j in 1:nz,i in 1:nt
        theta,z = (i-0.5)*dt,(j-0.5)*dz
        rj,rh = conical_radii(p,z)
        h = conical_normal_gap(p,theta,z)
        U = omega*rj
        push!(aspect,h/rh)
        push!(re_gap,rho*U*h/mu)
        push!(reduced_re,rho*U*h^2/(mu*rh))
        couette_dissipation += mu*U^2/h * rh*dt*dz/ca
    end
    # Independent exact normal-ray geometry, on extended cones at the ends.
    max_gap_discrepancy = 0.0
    for z in (0.0,0.05,0.1),theta in (0:255).*(2pi/256)
        rj,rh = conical_radii(p,z)
        lo,hi = 0.0,100e-6
        for _ in 1:60
            h = (lo+hi)/2
            r = rh-h*ca
            residual = hypot(r*cos(theta)-e,r*sin(theta))-(rj+h*sa*tand(10.0))
            if residual > 0
                lo = h
            else
                hi = h
            end
        end
        exact = (lo+hi)/2
        max_gap_discrepancy = max(max_gap_discrepancy,
            abs(conical_normal_gap(p,theta,z)/exact-1))
    end
    # Change only tangential wall-speed projections; hold pressure operator fixed.
    # At the approximate normal foot, omega ez×(x-e ex) projects onto:
    # Utheta=omega*(rfoot-e*cos(theta)), Us=omega*e*sin(alpha)*sin(theta).
    # This is an isolated sensitivity, NOT a complete higher-order Reynolds model.
    b = zeros(nt*nz)
    id(i,j) = i+nt*(j-1)
    for j in 1:nz,i in 1:nt
        theta,z = i*dt,(j-0.5)*dz
        _,rh = conical_radii(p,z)
        h = conical_normal_gap(p,theta,z)
        Utheta = omega*(rh-h*ca-e*cos(theta))
        flow = Utheta*h*dz/(2ca)
        k,l = id(i,j),id(i==nt ? 1 : i+1,j)
        b[k] -= flow; b[l] += flow
    end
    for j in 1:(nz+1),i in 1:nt
        theta,z = (i-0.5)*dt,(j-1)*dz
        _,rh = conical_radii(p,z)
        h = conical_normal_gap(p,theta,z)
        flow = omega*e*sa*sin(theta)*h*rh*dt/2
        j > 1 && (b[id(i,j-1)] -= flow)
        j <= nz && (b[id(i,j)] += flow)
    end
    perturbed = a.A \ b
    paths = [joinpath(@__DIR__,"..","..","simulation","reynolds",name)
             for name in ("conical.jl","cylindrical.jl","full_film.jl","pressure_vtk.jl",
                          "conical_checks.jl","conical-formulation.md")]
    result = Dict(
        "julia_version"=>string(VERSION),"grid"=>[nt,nz],"rpm"=>rpm,
        "viscosity_pa_s"=>mu,"assumed_density_kg_m3"=>rho,
        "normal_gap_over_housing_radius_range"=>collect(extrema(aspect)),
        "gap_reynolds_number_range"=>collect(extrema(re_gap)),
        "reduced_inertia_parameter_range"=>collect(extrema(reduced_re)),
        "couette_dissipation_component_w"=>couette_dissipation,
        "sampled_relative_normal_gap_discrepancy"=>max_gap_discrepancy,
        "wall_projection_pressure_relative_l2_change"=>norm(perturbed-baseline)/norm(baseline),
        "warning"=>"Scale estimates and isolated sensitivity, not physical validation or thermal prediction",
        "source_sha256"=>Dict(basename(path)=>bytes2hex(open(sha256,path)) for path in paths))
    ispath(out) && error("refusing overwrite")
    mkpath(out)
    open(joinpath(out,"scales-and-sensitivity.toml"),"w") do io
        TOML.print(io,result;sorted=true)
    end
    TOML.print(stdout,result;sorted=true)
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: geometry_review.jl FRESH_OUTPUT")
    review_conical(ARGS[1])
end
