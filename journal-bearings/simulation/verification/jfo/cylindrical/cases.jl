using LinearAlgebra, SparseArrays, Serialization, DelimitedFiles, Printf, Test, TOML, SHA
const ROOT=normpath(joinpath(@__DIR__,"../../../.."))
include(joinpath(ROOT,"bearing.jl"))
include(joinpath(ROOT,"mesh/2d/grid.jl"))
include(joinpath(ROOT,"simulation/reynolds/src/conical_jfo.jl"))

# One-radius/cosine-gap reference fixture, not the production circle geometry.
struct JFOCylinder
    journal_midpoint_radius::Float64
    housing_midpoint_radius::Float64
    bearing_length::Float64
    radial_clearance::Float64
    eccentricity_ratio::Float64
    semicone_angle_deg::Float64
    slope::Float64
end
film_radii(p::JFOCylinder,z::Real)=(p.journal_midpoint_radius*1e-3,p.journal_midpoint_radius*1e-3)
normal_film_gap(p::JFOCylinder,t::Real,z::Real)=p.radial_clearance*1e-3*(1-p.eccentricity_ratio*cos(t))
bearing_point(p::JFOCylinder,t,z,xi)=(p.journal_midpoint_radius*cos(t),p.journal_midpoint_radius*sin(t),Float64(z))

const MU=.0127
const RPM=48.1*60/(2pi)
const OMEGA=RPM*2pi/60
const R=.05
const CLEARANCE=.00291R
const PRESSURE_SCALE=MU*OMEGA*R^2/CLEARANCE^2

function cylinder_parameters(aspect,ecc;model="classical")
    model in ("classical","production") || throw(ArgumentError("unknown geometry model"))
    model=="classical" && return JFOCylinder(50.,50.,100aspect,1e3CLEARANCE,ecc,0.,0.)
    return BearingParameters(journal_midpoint_radius=50.,bearing_length=100aspect,
        radial_clearance=1e3CLEARANCE,eccentricity_ratio=ecc,semicone_angle_deg=0.)
end

function cylinder_grid(p::JFOCylinder,nt,nz;graded=false)
    theta=[2pi*(i-1)/nt for i in 1:nt]
    axial=collect(range(0.,p.bearing_length;length=nz+1))
    if graded
        theta .+= .2sin.(theta)
        axial .+= .1p.bearing_length.*sin.(2pi.*axial/p.bearing_length)./(2pi)
        axial[1],axial[end]=0.,p.bearing_length
    end
    bg=background_grid(theta,axial)
    nodes=[bearing_point(p,t,z,1.) for (t,z) in bg.nodes]
    areas=[p.journal_midpoint_radius*mod(theta[mod1(i+1,nt)]-theta[i],2pi)*(axial[j+1]-axial[j])
        for j in 1:nz for i in 1:nt]
    gaps=[1e3*normal_film_gap(p,theta[i]+mod(theta[mod1(i+1,nt)]-theta[i],2pi)/2,0.)
        for j in 1:nz for i in 1:nt]
    geometry=(;journal_midpoint_radius=p.journal_midpoint_radius,radial_clearance=0.,half_bearing_length=p.bearing_length/2,slope=0.)
    quality=SurfaceMeshQuality.tensor_grid_quality(theta,axial,geometry,1.)
    acceptance=SurfaceMeshQuality.mesh_acceptance(quality;valid_geometry=true,placement_stable=true,termination_reason=:direct_construction)
    return (;nodes,quads=bg.quads,radial_gaps=gaps,areas,n_theta=nt,n_axial=nz,quality,acceptance,
        parameters=p,theta_nodes=theta,z_nodes=axial)
end
function cylinder_grid(p::BearingParameters,nt,nz;graded=false)
    theta=[2pi*(i-1)/nt for i in 1:nt]
    axial=collect(range(0.,p.bearing_length;length=nz+1))
    if graded
        theta .+= .2sin.(theta)
        axial .+= .1p.bearing_length.*sin.(2pi.*axial/p.bearing_length)./(2pi)
        axial[1],axial[end]=0.,p.bearing_length
    end
    return housing_grid(p,nt,nz;theta_nodes=theta,z_nodes=axial)
end

function axial_count(aspect,nt)
    nt>=64 && nt%64==0 && ispow2(nt÷64) || throw(ArgumentError("counts must double from64"))
    return 2round(Int,64*max(aspect/2,.125)/2)*(nt÷64)
end
function cylinder_solution(aspect,ecc,pv0,nt;nz=axial_count(aspect,nt),model="classical")
    p=cylinder_parameters(aspect,ecc;model)
    return solve_conical_jfo(cylinder_grid(p,nt,nz);rpm=RPM,viscosity_pa_s=MU,
        end_pressure_pa=(0.,0.),cavitation_pressure_pa=pv0*PRESSURE_SCALE)
end

function plane_profile(grid,values;z_m=grid.parameters.bearing_length*.0005)
    axes=full_film_axes(grid);P=reshape(values,grid.n_theta,grid.n_axial)
    j=searchsortedlast(axes.z_centres,z_m)
    j>=1 && z_m<=last(axes.z_centres) || throw(ArgumentError("plane outside centre interpolation range"))
    z_m==axes.z_centres[j] && return copy(P[:,j])
    f=(z_m-axes.z_centres[j])/(axes.z_centres[j+1]-axes.z_centres[j])
    return (1-f).*P[:,j]+f.*P[:,j+1]
end
function periodic_profile(theta,values,query)
    return [begin
        q=mod(t-theta[1],2pi)+theta[1]
        i=searchsortedlast(theta,q);j=mod1(i+1,length(theta))
        right=j==1 ? theta[1]+2pi : theta[j]
        fraction=(q-theta[i])/(right-theta[i])
        (1-fraction)*values[i]+fraction*values[j]
    end for t in query]
end
function field_at(axes,p,theta,z)
    j=min(searchsortedlast(axes.z_centres,z),size(p,2)-1)
    1<=j && first(axes.z_centres)<=z<=last(axes.z_centres) || throw(ArgumentError("sample outside axial centres"))
    q=mod(theta-first(axes.theta_centres),2pi)+first(axes.theta_centres)
    i=searchsortedlast(axes.theta_centres,q);next=mod1(i+1,size(p,1))
    right=next==1 ? first(axes.theta_centres)+2pi : axes.theta_centres[next]
    t=(q-axes.theta_centres[i])/(right-axes.theta_centres[i])
    f=(z-axes.z_centres[j])/(axes.z_centres[j+1]-axes.z_centres[j])
    return (1-f)*((1-t)*p[i,j]+t*p[next,j])+f*((1-t)*p[i,j+1]+t*p[next,j+1])
end
function field_change(coarse,fine)
    axes=full_film_axes(coarse.grid);fine_axes=full_film_axes(fine.grid)
    fine_p=reshape(fine.pressure_pa,fine.grid.n_theta,fine.grid.n_axial)
    fine_f=reshape(fine.fill_fraction,fine.grid.n_theta,fine.grid.n_axial)
    pressure=Float64[];fill=Float64[]
    for z in axes.z_centres,t in axes.theta_centres
        push!(pressure,field_at(fine_axes,fine_p,t,z))
        push!(fill,field_at(fine_axes,fine_f,t,z))
    end
    areas=coarse.grid.areas
    return (;pressure_relative_L2=sqrt(sum(areas.*(pressure-coarse.pressure_pa).^2)/sum(areas.*pressure.^2)),
        fill_mean_absolute=sum(areas.*abs.(fill-coarse.fill_fraction))/sum(areas))
end

function observables(r)
    p=r.grid.parameters;L=p.bearing_length*1e-3
    force=housing_pressure_force(r.grid,r.pressure_pa)
    load=hypot(force[1],force[2])
    outward=[-r.axial_face_flow_m3_s[:,1];r.axial_face_flow_m3_s[:,end]]
    gross_out=sum(max.(outward,0.));gross_in=-sum(min.(outward,0.))
    mid=plane_profile(r.grid,r.pressure_pa)
    return (;P0=load*CLEARANCE^2/(MU*OMEGA*R^3*L),beta_deg=rad2deg(atan(abs(force[2]),abs(force[1]))),
        Q0=gross_out/(OMEGA*R^2*CLEARANCE),gross_out_m3_s=gross_out,gross_in_m3_s=gross_in,
        net_end_m3_s=sum(outward),midplane_peak_p0=maximum(mid)/PRESSURE_SCALE,
        minimum_fill=minimum(r.fill_fraction),cavitated_cells=r.diagnostics.cavitated_cells,
        iterations=r.diagnostics.iterations,local_balance=r.diagnostics.relative_max_cell_imbalance,
        global_balance=r.diagnostics.relative_global_imbalance,complementarity=r.diagnostics.complementarity_residual)
end
function row!(io,row;header=false)
    header && println(io,join(keys(row),'\t'))
    println(io,join(values(row),'\t'));flush(io)
end
function reference_rows()
    data,header=readdlm(joinpath(@__DIR__,"reference.csv"),',';header=true)
    return [(;L_D=Float64(r[1]),epsilon=Float64(r[2]),pv0=Float64(r[3]),P0=Float64(r[4]),
        beta_deg=Float64(r[5]),Q0=Float64(r[6])) for r in eachrow(data)]
end
function fresh_output(path)
    ispath(path) && error("refusing overwrite: $path")
    mkpath(path)
    paths=["cli.jl","bearing.jl","Project.toml","Manifest.toml","mesh/2d/grid.jl","mesh/2d/quality.jl",
        "mesh/2d/metrics.jl","mesh/vtk.jl","simulation/reynolds/src/conical_jfo.jl",
        "simulation/reynolds/src/fitted_jfo.jl",
        "simulation/reynolds/src/jfo_active_set.jl","simulation/reynolds/src/conical.jl",
        "simulation/reynolds/src/cylindrical.jl","simulation/reynolds/src/full_film.jl",
        "simulation/reynolds/src/face_fluxes.jl","simulation/reynolds/src/fitted_full_film.jl",
        "simulation/reynolds/src/pressure_vtk.jl","simulation/verification/command.py"]
    append!(paths,[relpath(joinpath(@__DIR__,name),ROOT) for name in readdir(@__DIR__)
        if endswith(name,".jl") || endswith(name,".py") || endswith(name,".csv") || endswith(name,".toml")])
    metadata=Dict("scope"=>"Steady finite cylindrical JFO numerical verification and literal historical comparison",
        "julia_version"=>string(VERSION),"physical_validation"=>false,
        "radius_m"=>R,"clearance_m"=>CLEARANCE,"mu_pa_s"=>MU,"omega_rad_s"=>OMEGA,"rpm"=>RPM,
        "pressure_scale_pa"=>PRESSURE_SCALE,"end_pressure_gauge_pa"=>[0.,0.],
        "supply"=>"Ungrooved submerged ends above subambient cavity; no hole or angular pin",
        "source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths))
    open(joinpath(path,"inputs.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    return path
end
