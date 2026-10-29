# Prescribed-pose cone/cylinder verification and the explicit 7mm fitted inlet.
using Test,LinearAlgebra,SparseArrays,TOML,SHA,Serialization
include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"../../../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../../../reynolds/src/pressure_sampling.jl"))
include(joinpath(@__DIR__,"../../swift-stieber/fitted_reference.jl"))
include(joinpath(@__DIR__,"transport_reference.jl"))
const TRANSPORT_ROOT=normpath(joinpath(@__DIR__,"../../../.."))
transport_parameters(angle;ecc=.6)=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
    radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=ecc)
function transport_mesh(angle,nt;fitted=false,graded=false,observer=nothing)
    p=transport_parameters(angle)
    fitted && return surface_feed_grid(p,nt,nt÷2;feed_diameter=7.,feed_z=50.,observer)
    theta=[2pi*t+(graded ? .2sin(2pi*t) : 0.) for t in (0:nt-1)./nt]
    z=[100*(t+(graded ? .1sin(2pi*t)/(2pi) : 0.)) for t in (0:nt÷2)./(nt÷2)]
    return housing_grid(p,nt,nt÷2;theta_nodes=theta,z_nodes=z)
end
function transport_pose(grid,ecc)
    p=transport_parameters(grid.parameters.semicone_angle_deg;ecc)
    return (;grid...,parameters=p)
end
function transport_solve(grid;ecc=.6,rpm=2000.)
    g=transport_pose(grid,ecc)
    fitted=hasproperty(g,:surface)
    return solve_conical_jfo(g;rpm,viscosity_pa_s=.0277,cavitation_pressure_pa=0.,
        end_pressure_pa=fitted ? (0.,0.) : (3e5,3e5),
        feed_pressure_pa=fitted ? 5e5 : nothing,eccentricity_angle_rad=3pi/2)
end

function transport_quad(grid,cell,axes=nothing)
    hasproperty(grid,:surface) && return grid.surface.localnodes[cell]
    a=isnothing(axes) ? full_film_axes(grid) : axes
    nt=grid.n_theta
    i,j=mod1(cell,nt),(cell-1)÷nt+1
    R=grid.parameters.housing_midpoint_radius*1e-3;ca=cosd(grid.parameters.semicone_angle_deg)
    x0,x1=a.theta_edges[i]*R,a.theta_edges[i+1]*R
    y0,y1=a.z_edges[j]/ca,a.z_edges[j+1]/ca
    return [x0 y0;x1 y0;x1 y1;x0 y1]
end
function transport_integral(grid,f;order=3)
    points,weights=fitted_gauss(order)
    result=zeros(length(grid.quads))
    axes=hasproperty(grid,:surface) ? nothing : full_film_axes(grid)
    for cell in eachindex(result)
        xy=transport_quad(grid,cell,axes)
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
            shape=[(1-u)*(1-v),(1+u)*(1-v),(1+u)*(1+v),(1-u)*(1+v)]/4
            du=[-(1-v),1-v,1+v,-(1+v)]/4;dv=[-(1-u),-(1+u),1+u,1-u]/4
            x,xu,xv=xy'*shape,xy'*du,xy'*dv
            jac=xu[1]*xv[2]-xu[2]*xv[1]
            result[cell]+=f(x...)*jac*wu*wv
        end
    end
    return result
end
function transport_centres(grid)
    if hasproperty(grid,:surface)
        s=grid.surface
        return s.centre[:,1]./s.R,s.centre[:,2].*s.ca
    end
    a=full_film_axes(grid)
    return vec([t for t in a.theta_centres,z in a.z_centres]),vec([z for t in a.theta_centres,z in a.z_centres])
end

function transport_exact_boundary(grid,angle,rpm;order=160)
    points,weights=fitted_gauss(order);ends=zeros(2)
    half=acos(-.25)
    for (j,z) in enumerate((0.,.1)),(u,w) in zip(points,weights)
        theta=pi/2+half*u
        ref=ConicalTransportReference.manufactured(theta,z,angle,rpm)
        g=ConicalTransportReference.geometry(theta,z,angle,.6;pose=.37)
        ends[j]+=(j==1 ? 1 : -1)*g.K*g.rh*g.ca*ref.pz*half*w
    end
    feed=0.
    if hasproperty(grid,:surface)
        s=grid.surface;nodes,ws=fitted_gauss(8)
        for face in findall(==(:feed),s.boundary)
            a,b=fitted_reference_edge(s,face);t=b-a
            for (u,w) in zip(nodes,ws)
                x=a+(u+1)/2*t
                ref=ConicalTransportReference.manufactured(x[1]/s.R,x[2]*s.ca,angle,rpm)
                feed+=w/2*(ref.qx*t[2]-ref.qy*t[1])
            end
        end
    end
    return [ends;feed]
end

function transport_exact_force(grid,angle;order=64)
    points,weights=fitted_gauss(order);force=zeros(3);k=tand(angle)
    for (v,wv) in zip(points,weights)
        z=(v+1)*.05;half=acos(-.25cos(2pi*z/.1))
        for (u,wu) in zip(points,weights)
            theta=pi/2+half*u
            ref=ConicalTransportReference.manufactured(theta,z,angle,0.)
            rh=.05005+(.05-z)*k
            force .+= ref.pressure*rh*half*wu*.05*wv.*[cos(theta),sin(theta),k]
        end
    end
    if hasproperty(grid,:surface)
        s=grid.surface;centre=[s.R*pi/2,.05/s.ca];points,weights=fitted_gauss(8)
        for face in findall(==(:feed),s.boundary)
            a,b=fitted_reference_edge(s,face);first,second=a-centre,b-centre
            jac=abs(first[1]*second[2]-first[2]*second[1])
            for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
                t,q=(u+1)/2,(v+1)/2;x=centre+t*first+(1-t)*q*second
                theta,z=x[1]/s.R,x[2]*s.ca;rh=.05005+(.05-z)*k
                ref=ConicalTransportReference.manufactured(theta,z,angle,0.)
                weight=jac*(1-t)*wu*wv/4*rh/s.R
                force .-= ref.pressure*weight.*[s.ca*cos(theta),s.ca*sin(theta),sind(angle)]
            end
        end
    end
    return force
end

function transport_mms(grid,rpm)
    angle=grid.parameters.semicone_angle_deg
    R=grid.parameters.housing_midpoint_radius*1e-3;ca=cosd(angle)
    ref(t,z)=ConicalTransportReference.manufactured(t,z,angle,rpm)
    if hasproperty(grid,:surface)
        s=grid.surface
        boundary_pressure(x,y,tag)=ref(x/R,y*ca).pressure
        faces=fitted_face_operator(s,grid.parameters;rpm,viscosity_pa_s=.0277,
            eccentricity_angle_rad=.37,boundary_pressure)
        @test all(ref(s.mid[f][1]/R,s.mid[f][2]*ca).fill==1. for f in findall(==(:feed),s.boundary))
        C=jfo_content_operator(faces,length(grid.quads))
        a=(;A=faces.divergence*faces.pressure_flux,B=-faces.divergence*C,
            b=-faces.divergence*faces.constant,faces,pressure_scale_pa=1e7)
    else
        a=assemble_conical_jfo(grid;rpm,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),
            cavitation_pressure_pa=0.,eccentricity_angle_rad=.37)
    end
    source=transport_integral(grid,(x,y)->ref(x/R,y*ca).source)
    t,z=transport_centres(grid);exact=[ref(ti,zi) for (ti,zi) in zip(t,z)]
    ep,ef=getproperty.(exact,:pressure),getproperty.(exact,:fill)
    state=solve_jfo_active_set(a.A,a.B,a.b+source,a.pressure_scale_pa)
    P,F=state.pressure,1 .-state.void
    flux=jfo_face_fluxes(a,P,F)
    weighted=grid.areas
    qe=transport_exact_boundary(grid,angle,rpm);fe=transport_exact_force(grid,angle)
    force=collect(housing_pressure_force(grid,P))
    return (;pressure_relative_L2=sqrt(sum((P-ep).^2 .*weighted)/sum(ep.^2 .*weighted)),
        fill_mean_absolute=sum(abs.(F-ef).*weighted)/sum(weighted),
        boundary_flow_relative_L2=norm(collect(values(flux.boundary_outflow))-qe)/norm(qe),
        force_relative_L2=norm(force-fe)/norm(fe),
        local_balance=norm(flux.balance-source,Inf)/max(maximum(abs,a.faces.advective),maximum(abs,source)),
        minimum_pressure=minimum(P),minimum_fill=minimum(F),maximum_fill=maximum(F),
        complementarity=maximum(abs,(P/a.pressure_scale_pa).*(1 .-F)),iterations=state.iterations)
end

# Fitted comparisons use actual geometric positions, not row/column restriction.
function transport_sample(result,theta,z;axes=nothing)
    grid=result.grid
    if !hasproperty(grid,:surface)
        a=isnothing(axes) ? full_film_axes(grid) : axes
        nt=grid.n_theta
        period=2pi;tt=mod(theta-a.theta_centres[1],period)+a.theta_centres[1]
        i=searchsortedlast(a.theta_centres,tt);next=mod1(i+1,nt)
        right=i==nt ? a.theta_centres[1]+period : a.theta_centres[next]
        w=(tt-a.theta_centres[i])/(right-a.theta_centres[i])
        zz=[0.;a.z_centres;.1];j=clamp(searchsortedlast(zz,z),1,length(zz)-1)
        u=(z-zz[j])/(zz[j+1]-zz[j])
        function field(values,boundary)
            M=reshape(values,nt,grid.n_axial)
            row(k)=k==1 ? boundary : k==length(zz) ? boundary : (1-w)*M[i,k-1]+w*M[next,k-1]
            return (1-u)*row(j)+u*row(j+1)
        end
        return field(result.pressure_pa,3e5),field(result.fill_fraction,1.)
    end
    error("use fitted plane sampling for fitted fields")
end
function transport_field_change(coarse,fine)
    theta,z=transport_centres(coarse.grid);axes=full_film_axes(fine.grid)
    values=[transport_sample(fine,t,y;axes) for (t,y) in zip(theta,z)]
    p,f=first.(values),last.(values);w=coarse.grid.areas
    return (;pressure_change=sqrt(sum((p-coarse.pressure_pa).^2 .*w)/sum(p.^2 .*w)),
        fill_change=sum(abs.(f-coarse.fill_fraction).*w)/sum(w))
end
