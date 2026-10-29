using LinearAlgebra, DelimitedFiles, Serialization, Printf, Test, TOML, SHA
const SOURCE = normpath(joinpath(@__DIR__, "../../.."))
include(joinpath(SOURCE, "bearing.jl"))
include(joinpath(SOURCE, "mesh/2d/grid.jl"))
include(joinpath(SOURCE, "simulation/reynolds/src/full_film.jl"))

struct DonaldsonCylinder
    radius_mm::Float64
    bearing_length::Float64
    radial_clearance::Float64
    eccentricity_ratio::Float64
    semicone_angle_deg::Float64
    slope::Float64
end
film_radii(p::DonaldsonCylinder, z::Real) = (p.radius_mm*1e-3, p.radius_mm*1e-3)
normal_film_gap(p::DonaldsonCylinder, theta::Real, z::Real) =
    p.radial_clearance*1e-3*(1+p.eccentricity_ratio*cos(theta))
function bearing_point(p::DonaldsonCylinder, theta, z, xi)
    xi == 1 || error("nominal housing mapping only")
    return (p.radius_mm*cos(theta), p.radius_mm*sin(theta), Float64(z))
end

function refinement_axial_count(aspect, nt)
    isfinite(aspect) && aspect > 0 || throw(ArgumentError("invalid L/D"))
    nt >= 128 && nt % 128 == 0 && ispow2(nt ÷ 128) ||
        throw(ArgumentError("refinement counts must be 128 times a power of two"))
    base_nz = 2round(Int, 128*clamp(aspect/2,.25,1.)/2)
    return base_nz*(nt ÷ 128)
end

function benchmark_grid(p, nt, nz)
    theta_nodes = [(i-1)*2pi/nt for i in 1:nt]
    z_nodes = collect(range(0., p.bearing_length; length=nz+1))
    background = background_grid(theta_nodes, z_nodes)
    nodes = [bearing_point(p,t,z,1.) for (t,z) in background.nodes]
    gaps = [1e3*normal_film_gap(p,(i-.5)*2pi/nt,0.) for j in 1:nz for i in 1:nt]
    areas = fill(p.radius_mm*(2pi/nt)*(p.bearing_length/nz), nt*nz)
    geometry = (;journal_midpoint_radius=p.radius_mm,radial_clearance=0.,
        half_bearing_length=p.bearing_length/2,slope=0.)
    quality = SurfaceMeshQuality.tensor_grid_quality(theta_nodes,z_nodes,geometry,1.)
    acceptance = SurfaceMeshQuality.mesh_acceptance(quality;valid_geometry=true,
        placement_stable=true,termination_reason=:direct_construction)
    return (;nodes,quads=background.quads,radial_gaps=gaps,areas,n_theta=nt,n_axial=nz,
        quality,acceptance,parameters=p,theta_nodes,z_nodes)
end

const RPM = 2000.
const MU = .02752
function infinite_load(p)
    r, length, c = p.radius_mm*1e-3, p.bearing_length*1e-3, p.radial_clearance*1e-3
    e = p.eccentricity_ratio
    return 12pi*MU*(RPM*2pi/60)*r^3*length*e/(c^2*(2+e^2)*sqrt(1-e^2))
end

# Periodic linear interpolation between cell centres; no pressure extrapolation.
function angular_value(values, theta)
    nt = length(values)
    coordinate = mod(theta,2pi)*nt/(2pi)+.5
    lower = floor(Int,coordinate)
    fraction = coordinate-lower
    return (1-fraction)*values[mod1(lower,nt)] + fraction*values[mod1(lower+1,nt)]
end

function observables(p, result)
    nt,nz = result.grid.n_theta,result.grid.n_axial
    pressure = reshape(result.pressure_pa,nt,nz)
    r,c,length = p.radius_mm*1e-3,p.radial_clearance*1e-3,p.bearing_length*1e-3
    scale = 6MU*(RPM*2pi/60)*r^2/c^2
    midline = (pressure[:,nz÷2]+pressure[:,nz÷2+1])/2/scale
    # Native FV end flux divided by face mobility gives the signed dp/dz.
    h = [normal_film_gap(p,(i-.5)*2pi/nt,0.) for i in 1:nt]
    mobility = h.^3 .* (r*2pi/nt)/(12MU)
    lower_gradient = -result.axial_face_flow_m3_s[:,1]./mobility .* (r/scale)
    upper_gradient = -result.axial_face_flow_m3_s[:,end]./mobility .* (r/scale)
    force = result.housing_pressure_force_n
    load = hypot(force[1],force[2])
    couette_torque = MU*(RPM*2pi/60)*r^3*length*(2pi/nt)*sum(1 ./ h)
    pressure_torque = .5*p.eccentricity_ratio*c*force[2]
    friction = (couette_torque+pressure_torque)/(c*load)
    index = argmax(midline)
    return (;load_ratio=load/infinite_load(p),
        pressure_160=angular_value(midline,deg2rad(160.)),
        lower_end_gradient_160=angular_value(lower_gradient,deg2rad(160.)),
        upper_end_gradient_160=angular_value(upper_gradient,deg2rad(160.)),
        midline_maximum=midline[index],maximum_angle_deg=(index-.5)*360/nt,
        shaft_friction_R_c=friction,shaft_torque_nm=couette_torque+pressure_torque,
        force_parallel_relative=abs(force[1])/load,
        odd_symmetry_relative=maximum(abs,pressure+reverse(pressure;dims=1))/maximum(abs,pressure),
        axial_symmetry_relative=maximum(abs,pressure-reverse(pressure;dims=2))/maximum(abs,pressure),
        relative_cell_imbalance=result.diagnostics.relative_l1_cell_imbalance,
        relative_global_imbalance=result.diagnostics.relative_global_imbalance)
end

function axial_profile(pressure, z_centres, target_z)
    size(pressure,2) == length(z_centres) || throw(DimensionMismatch("axial coordinates"))
    all(isfinite,z_centres) && all(diff(z_centres) .> 0) ||
        throw(ArgumentError("axial centres must be finite and strictly increasing"))
    isfinite(target_z) && first(z_centres) <= target_z <= last(z_centres) ||
        throw(ArgumentError("sampling outside cell centres requires explicit boundary data"))
    lower = searchsortedlast(z_centres,target_z)
    (lower == length(z_centres) || z_centres[lower] == target_z) &&
        return copy(pressure[:,lower])
    weight = (target_z-z_centres[lower])/(z_centres[lower+1]-z_centres[lower])
    return (1-weight).*pressure[:,lower] + weight.*pressure[:,lower+1]
end

function save_profiles(p,result,tag;directory=joinpath(@__DIR__,"profiles"))
    mkpath(directory)
    nt,nz=result.grid.n_theta,result.grid.n_axial
    r,c=p.radius_mm*1e-3,p.radial_clearance*1e-3
    scale=6MU*(RPM*2pi/60)*r^2/c^2
    pressure=reshape(result.pressure_pa,nt,nz)/scale
    z_centres=(result.grid.z_nodes[1:end-1]+result.grid.z_nodes[2:end])/2
    midplane_z,quarter_z=p.bearing_length/2,p.bearing_length/4
    midplane=axial_profile(pressure,z_centres,midplane_z)
    quarter=axial_profile(pressure,z_centres,quarter_z)
    open(joinpath(directory,"$(tag)-circumferential.tsv"),"w") do io
        println(io,"theta_deg\tP_midplane\tP_quarter\tP_first_axial_cell\tz_midplane_mm\tz_quarter_mm\tz_first_axial_cell_mm")
        for i in 1:nt
            println(io,join(((i-.5)*360/nt,midplane[i],quarter[i],pressure[i,1],
                midplane_z,quarter_z,z_centres[1]),'\t'))
        end
    end
    open(joinpath(directory,"$(tag)-axial.tsv"),"w") do io
        println(io,"w_from_midplane\tP_at_160_deg")
        println(io,join((-p.bearing_length/(2p.radius_mm),0.),'\t'))
        for j in 1:nz
            println(io,join((((j-.5)*p.bearing_length/nz-p.bearing_length/2)/p.radius_mm,
                angular_value(pressure[:,j],deg2rad(160.))),'\t'))
        end
        println(io,join((p.bearing_length/(2p.radius_mm),0.),'\t'))
    end
end

function production_infinite_load(p;samples=4096)
    rh=p.housing_midpoint_radius*1e-3
    h=[normal_film_gap(p,(i-.5)*2pi/samples,0.) for i in 1:samples]
    effective=sum(h.^(-2))/sum(h.^(-3))
    integral=2pi/samples*sum((h.-effective).*cos.([(i-.5)*2pi/samples for i in 1:samples])./h.^3)
    speed=RPM*2pi/60*p.journal_midpoint_radius*1e-3
    return abs(6MU*speed*rh^2*p.bearing_length*1e-3*integral)
end

function fresh_output(path; suite)
    ispath(path) && error("refusing to overwrite $path")
    mkpath(path)
    files=["cli.jl","simulation/verification/command.py","bearing.jl","mesh/2d/grid.jl","mesh/2d/quality.jl","mesh/2d/metrics.jl",
        "simulation/reynolds/src/full_film.jl","simulation/reynolds/src/face_fluxes.jl",
        "simulation/reynolds/src/fitted_full_film.jl","Project.toml","Manifest.toml"]
    append!(files,[relpath(joinpath(@__DIR__,name),SOURCE) for name in readdir(@__DIR__)
        if endswith(name,".jl") || endswith(name,".py") || endswith(name,".csv") || endswith(name,".toml")])
    metadata=Dict("suite"=>suite,"julia_version"=>string(VERSION),
        "rpm"=>RPM,"viscosity_pa_s"=>MU,"physical_validation"=>false,
        "source_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(SOURCE,name)))) for name in files))
    open(joinpath(path,"inputs.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    return path
end
