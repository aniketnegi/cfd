using LinearAlgebra, SparseArrays, Serialization, DelimitedFiles, Printf, Test, TOML, SHA
const SOURCE=normpath(joinpath(@__DIR__,"../../.."))
include(joinpath(SOURCE,"bearing.jl"))
include(joinpath(SOURCE,"mesh/2d/grid.jl"))
include(joinpath(SOURCE,"simulation/reynolds/src/swift_stieber.jl"))
const MU=.0277
cone(;angle=10.,ecc=.6)=BearingParameters(journal_midpoint_radius=50.,
    bearing_length=100.,radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=ecc)
function grid_for(p,nt,nz;graded=false)
    if graded
        t=(0:nt-1)./nt;s=(0:nz)./nz
        theta=.17*2pi/nt .+ 2pi .* t .+ 0.25 .* sin.(2pi .* t)
        axial=p.bearing_length .* (s .+ 0.15 .* sin.(2pi .* s)./(2pi))
        axial[1],axial[end]=0.,p.bearing_length
        return housing_grid(p,nt,nz;theta_nodes=theta,z_nodes=axial)
    end
    return housing_grid(p,nt,nz;theta_nodes=[(i-.5)*2pi/nt for i in 1:nt])
end
function assemble_case(p,nt,nz;mu=MU,rpm=2000.,pose=pi,graded=false)
    assemble_full_film(grid_for(p,nt,nz;graded);rpm,viscosity_pa_s=mu,
        end_pressure_pa=(0.,0.),eccentricity_angle_rad=pose)
end
function obstacle(a;mu=MU,rpm=2000.,groove=false,source=zeros(length(a.b)),computed=nothing)
    p=a.grid.parameters;nt,nz=size(a.gt)
    pressure_scale=max(6mu*abs(rpm*2pi/60)*(.05)^2/(50e-6)^2,1.)
    rhs=a.b+source
    flow_scale=max(maximum(abs,a.adv),maximum(abs,rhs),1e-30)
    fixed=falses(nt*nz);groove && (fixed[nt:nt:end].=true)
    free=findall(.!fixed)
    state=isnothing(computed) ? pressure_obstacle(a.A[free,free],rhs[free],pressure_scale,flow_scale) : computed
    pressure=zeros(nt*nz);pressure[free]=state.pressure
    reaction=a.A*pressure-rhs
    flux=full_film_fluxes(a,pressure)
    positive=pressure .> 1e-10*pressure_scale
    force=housing_pressure_force(a.grid,pressure)
    ends=(-sum(flux.qz[:,1]),sum(flux.qz[:,end]))
    cavity=sum(reaction[.!positive .& .!fixed]);pin=sum(reaction[fixed])
    accounting=sum(ends)-sum(source)-cavity-pin
    global_scale=max(abs(sum(ends))+sum(abs,source)+abs(cavity)+abs(pin),flow_scale)
    stats=(;iterations=state.iterations,
        positive_pressure_residual=any(positive) ? maximum(abs,reaction[positive])/flow_scale : 0.,
        complementarity=maximum(abs,(pressure[free]/pressure_scale).*(reaction[free]/flow_scale)),
        minimum_free_reaction_relative=minimum(reaction[free])/flow_scale,
        flux_reconstruction_relative=maximum(abs,flux.balance-source-reaction)/flow_scale,
        defect_accounting_relative=abs(accounting)/global_scale,
        cavity_full_liquid_flux_defect_m3_s=cavity,pinned_line_exchange_m3_s=pin,
        lower_outflow_m3_s=ends[1],upper_outflow_m3_s=ends[2],
        positive_cells=count(positive),force_x_N=force[1],force_y_N=force[2],force_z_N=force[3])
    return (;a,pressure,reaction,source,flux,force,ends,positive,pressure_scale,flow_scale,stats)
end
function solve_case(p,nt,nz;mu=MU,rpm=2000.,pose=pi,graded=false,groove=false)
    graded && groove && error("pinned-column test requires its explicitly aligned uniform grid")
    a=assemble_case(p,nt,nz;mu,rpm,pose,graded)
    pins=groove ? collect(nt:nt:nt*nz) : Int[]
    result=solve_swift_stieber(a.grid;rpm,viscosity_pa_s=mu,end_pressure_pa=(0.,0.),
        eccentricity_angle_rad=pose,pinned_cells=pins)
    computed=(;pressure=result.pressure_above_cavitation_pa[findall(.!result.pinned_cells)],iterations=result.diagnostics.iterations)
    obstacle(a;mu,rpm,groove,computed)
end
function midplane(a,pressure)
    z=a.axes.z_centres;target=a.grid.parameters.bearing_length*.0005
    matrix=reshape(pressure,size(a.gt))
    index=searchsortedlast(z,target)
    if z[index]==target
        return (;profile=matrix[:,index],target,lower=z[index],upper=z[index])
    end
    fraction=(target-z[index])/(z[index+1]-z[index])
    return (;profile=(1-fraction)*matrix[:,index]+fraction*matrix[:,index+1],target,
        lower=z[index],upper=z[index+1])
end
function save_field(OUTPUT,name,result;extra=(;))
    a=result.a;m=midplane(a,result.pressure)
    mkpath(joinpath(OUTPUT,"fields"))
    serialize(joinpath(OUTPUT,"fields",name*".jls"),(;
        pressure=result.pressure,reaction=result.reaction,source=result.source,
        theta=a.axes.theta_centres,z=a.axes.z_centres,gap=a.normal_gap_m,
        theta_edges=a.axes.theta_edges,z_edges=a.axes.z_edges,areas=a.grid.areas,
        force=result.force,ends=result.ends,midplane=m,extra...))
end
function row!(io,row;header=false)
    header && println(io,join(keys(row),'\t'))
    println(io,join(values(row),'\t'));flush(io)
end
function check_result(result)
    s=result.stats
    @test minimum(result.pressure)>=0
    @test s.positive_pressure_residual<1e-10
    @test s.complementarity<1e-10
    @test s.minimum_free_reaction_relative> -1e-10
    @test s.flux_reconstruction_relative<1e-10
    @test s.defect_accounting_relative<1e-10
end
function gauss(n)
    eigenproblem=eigen(SymTridiagonal(zeros(n),[i/sqrt(4i^2-1) for i in 1:n-1]))
    return eigenproblem.values,2 .* eigenproblem.vectors[1,:].^2
end

function fresh_output(path;suite)
    ispath(path) && error("refusing overwrite: $path")
    mkpath(path)
    files=["cli.jl","simulation/verification/command.py","bearing.jl","mesh/2d/grid.jl","mesh/2d/quality.jl","mesh/2d/metrics.jl",
        "simulation/reynolds/src/full_film.jl","simulation/reynolds/src/swift_stieber.jl",
        "simulation/reynolds/src/jfo_active_set.jl","simulation/reynolds/src/jfo_slider.jl",
        "Project.toml","Manifest.toml"]
    append!(files,[relpath(joinpath(@__DIR__,name),SOURCE) for name in readdir(@__DIR__)
        if endswith(name,".jl") || endswith(name,".py") || endswith(name,".csv") || endswith(name,".toml")])
    open(joinpath(path,"inputs.toml"),"w") do io
        TOML.print(io,Dict("suite"=>suite,"julia_version"=>string(VERSION),
            "physical_validation"=>false,"rpm"=>2000.,"mu_conical_pa_s"=>MU,"mu_cylindrical_pa_s"=>.02752,
            "boundary_scope"=>"Explicit numerical test conditions, not a recovered Gangrade inlet",
            "reference_journal_radius_mm"=>50.,"radial_clearance_mm"=>.05,
            "length_rule_mm"=>"Cylinder:100*L_D; cone:100",
            "end_pressure_pa"=>[0.,0.],"cavitation_pressure_pa"=>0.,
            "command_line"=>copy(ARGS),
            "source_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(SOURCE,name)))) for name in files));sorted=true)
    end
    path
end

struct ClassicalCylinder
    radius_mm::Float64
    bearing_length::Float64
    radial_clearance::Float64
    eccentricity_ratio::Float64
    semicone_angle_deg::Float64
    slope::Float64
end
film_radii(p::ClassicalCylinder,z::Real)=(p.radius_mm*1e-3,p.radius_mm*1e-3)
normal_film_gap(p::ClassicalCylinder,t::Real,z::Real)=p.radial_clearance*1e-3*(1+p.eccentricity_ratio*cos(t))
bearing_point(p::ClassicalCylinder,t,z,xi)=(p.radius_mm*cos(t),p.radius_mm*sin(t),Float64(z))

function cylinder_grid(p::ClassicalCylinder,nt,nz)
    theta=[(i-.5)*2pi/nt for i in 1:nt]
    axial=collect(range(0.,p.bearing_length;length=nz+1))
    bg=background_grid(theta,axial)
    nodes=[bearing_point(p,t,z,1.) for (t,z) in bg.nodes]
    radial_gaps=[1e3*normal_film_gap(p,i*2pi/nt,0.) for j in 1:nz for i in 1:nt]
    areas=fill(p.radius_mm*(2pi/nt)*(p.bearing_length/nz),nt*nz)
    geometry=(;journal_midpoint_radius=p.radius_mm,radial_clearance=0.,half_bearing_length=p.bearing_length/2,slope=0.)
    quality=SurfaceMeshQuality.tensor_grid_quality(theta,axial,geometry,1.)
    acceptance=SurfaceMeshQuality.mesh_acceptance(quality;valid_geometry=true,placement_stable=true,termination_reason=:direct_construction)
    return (;nodes,quads=bg.quads,radial_gaps,areas,n_theta=nt,n_axial=nz,quality,acceptance,parameters=p,theta_nodes=theta,z_nodes=axial)
end
function cylinder_grid(p::BearingParameters,nt,nz)
    return housing_grid(p,nt,nz;theta_nodes=[(i-.5)*2pi/nt for i in 1:nt])
end


function cylindrical_case(aspect,ecc,model,nt;groove=true)
    nz=round(Int,nt*max(aspect/2,.125))
    p=model=="classical" ? ClassicalCylinder(50.,100aspect,.05,ecc,0.,0.) :
        BearingParameters(journal_midpoint_radius=50.,bearing_length=100aspect,
            radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=ecc)
    angle=model=="classical" ? 0. : pi
    a=assemble_full_film(cylinder_grid(p,nt,nz);rpm=2000.,viscosity_pa_s=.02752,
        end_pressure_pa=(0.,0.),eccentricity_angle_rad=angle)
    fixed=falses(nt*nz)
    groove && (fixed[nt:nt:end].=true)
    free=findall(.!fixed)
    scale=6*.02752*(2000*2pi/60)*.05^2/(.05e-3)^2
    flow_scale=maximum(abs,a.adv)
    elapsed=@elapsed state=pressure_obstacle(a.A[free,free],a.b[free],scale,flow_scale)
    pressure=zeros(nt*nz);pressure[free]=state.pressure
    reaction=a.A*pressure-a.b
    flux=full_film_fluxes(a,pressure)
    @assert maximum(abs,flux.balance-reaction)/flow_scale<1e-10
    @assert minimum(pressure)>=0
    @assert minimum(reaction[free])/flow_scale>=-1e-10
    positive=pressure .> 1e-9*scale
    residual=maximum(abs,reaction[positive])/flow_scale
    complementarity=maximum(abs,(pressure[free]/scale).*(reaction[free]/flow_scale))
    @assert residual<1e-10 && complementarity<1e-10
    force=tensor_pressure_force(p,pressure,a.axes.theta_edges,a.axes.z_edges)
    load=hypot(force[1],force[2])
    outflow=-sum(flux.qz[:,1])+sum(flux.qz[:,end])
    cavity=sum(reaction[.!positive .& .!fixed])
    exchange=sum(reaction[fixed])
    @assert abs(outflow-cavity-exchange)/max(abs(outflow),flow_scale)<1e-10
    signed=a.A\a.b
    difference=norm(pressure-max.(signed,0.))/norm(pressure)
    matrix=reshape(pressure,nt,nz)
    last_full=findlast(>(1e-9*scale),matrix[:,max(1,nz÷2)])
    omega=2000*2pi/60;r=.05;c=.05e-3;length=100aspect*1e-3
    load0=load*c^2/(.02752*omega*r^3*length)
    flow0=outflow/(r*(omega*r)*c)
    row=(;model,L_D=aspect,epsilon=ecc,nt,nz,groove,load_n=load,load_P0=load0,
        attitude_deg=rad2deg(atan(abs(force[2]),abs(force[1]))),flow_Q0=flow0,
        max_pressure_pa=maximum(pressure),positive_cells=count(positive),iterations=state.iterations,
        positive_pressure_residual=residual,complementarity,
        clipping_relative_difference=difference,
        midplane_last_positive_theta_deg=last_full*360/nt,
        cavity_fully_liquid_flux_defect_m3_s=cavity,pinned_line_fully_liquid_exchange_m3_s=exchange,
        pressure_end_outflow_m3_s=outflow,seconds=elapsed)
    return row,(;pressure_pa=pressure,normal_gap_m=a.normal_gap_m,force_n=force,
        theta_centres=a.axes.theta_centres,z_centres=a.axes.z_centres,
        reaction_m3_s=reaction,signed_pressure_pa=signed)
end
