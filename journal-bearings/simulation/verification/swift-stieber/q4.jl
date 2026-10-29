# Verification-only nodal Q4 Galerkin assembly of the shared reduced cone PDE.
# The public housing_grid owns all nodes, coordinates and element connectivity.
using LinearAlgebra, SparseArrays

function q4_gauss(order)
    order>=2 || throw(ArgumentError("Q4 needs at least 2x2 quadrature"))
    rule=eigen(SymTridiagonal(zeros(order),[i/sqrt(4i^2-1) for i in 1:order-1]))
    return rule.values,2 .* rule.vectors[1,:].^2
end

function assemble_q4(grid;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),
        eccentricity_angle_rad=3pi/2,quadrature_order=2)
    axes=full_film_axes(grid)
    all(isfinite,(rpm,viscosity_pa_s,end_pressure_pa...,eccentricity_angle_rad)) && viscosity_pa_s>0 ||
        throw(ArgumentError("operating inputs must be finite with positive viscosity"))
    minimum(end_pressure_pa)>=0 || throw(ArgumentError("ends cannot be below the zero gauge pressure floor"))
    p=grid.parameters;nt,nz=grid.n_theta,grid.n_axial
    count=length(grid.nodes);ca=cosd(p.semicone_angle_deg)
    theta=[t for z in axes.z_edges for t in grid.theta_nodes]
    z=[z for z in axes.z_edges for t in grid.theta_nodes]
    hnode=[normal_film_gap(p,theta[i]-eccentricity_angle_rad,z[i]) for i in 1:count]
    fixed=falses(count);values=zeros(count)
    fixed[1:nt].=true;values[1:nt].=end_pressure_pa[1]
    fixed[end-nt+1:end].=true;values[end-nt+1:end].=end_pressure_pa[2]
    points,weights=q4_gauss(quadrature_order)
    rows,columns,entries=Int[],Int[],Float64[];b=zeros(count)
    for (cell,ids) in enumerate(grid.quads)
        i,j=mod1(cell,nt),(cell-1)÷nt+1
        dt,dz=axes.theta_widths[i],axes.z_widths[j]
        local_A,local_b=zeros(4,4),zeros(4)
        heights=hnode[collect(ids)]
        for (xi,wx) in zip(points,weights),(eta,wy) in zip(points,weights)
            N=[(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)]/4
            Nt=[-(1-eta),1-eta,1+eta,-(1+eta)]/(2dt)
            Nz=[-(1-xi),-(1+xi),1+xi,1-xi]/(2dz)
            zq=axes.z_edges[j]+(1+eta)*dz/2
            rj,rh=film_radii(p,zq)
            h=dot(N,heights);K=h^3/(12viscosity_pa_s)
            weight=wx*wy*dt*dz/4
            local_A .+= weight.*(K/(rh*ca).*(Nt*Nt')+K*rh*ca.*(Nz*Nz'))
            local_b .+= weight*(rpm*2pi/60*rj*h/(2ca)).*Nt
        end
        for a in 1:4
            b[ids[a]]+=local_b[a]
            for c in 1:4
                push!(rows,ids[a]);push!(columns,ids[c]);push!(entries,local_A[a,c])
            end
        end
    end
    pressure_scale=max(maximum(end_pressure_pa),12viscosity_pa_s*abs(rpm*2pi/60)*
        p.journal_midpoint_radius*1e-3*(2pi*p.housing_midpoint_radius*1e-3)/minimum(hnode)^2,1.)
    return (;A=sparse(rows,columns,entries,count,count),b,grid,axes,theta,z,hnode,fixed,values,pressure_scale)
end

function q4_plane(a,pressure;z_m=a.grid.parameters.bearing_length*.0005)
    nt,nz=a.grid.n_theta,a.grid.n_axial
    z=a.axes.z_edges
    first(z)<=z_m<=last(z) || throw(ArgumentError("plane outside bearing"))
    index=min(searchsortedlast(z,z_m),nz)
    fraction=(z_m-z[index])/(z[index+1]-z[index])
    P=reshape(pressure,nt,nz+1)
    return (1-fraction).*P[:,index]+fraction.*P[:,index+1]
end

function solve_q4(grid;closure=:swift,kwargs...)
    closure in (:swift,:signed) || throw(ArgumentError("unknown closure"))
    a=assemble_q4(grid;kwargs...)
    free=findall(.!a.fixed);fixed=findall(a.fixed)
    rhs=a.b[free]-a.A[free,fixed]*a.values[fixed]
    A=a.A[free,free];flow_scale=max(maximum(abs,rhs),maximum(diag(A))*max(maximum(a.values),1.),1e-30)
    P=copy(a.values);iterations=1
    if closure==:swift
        state=pressure_obstacle(A,rhs,a.pressure_scale,flow_scale)
        P[free]=state.pressure;iterations=state.iterations
    else
        P[free]=A\rhs
    end
    reaction=A*P[free]-rhs
    positive=P[free].>1e-10*a.pressure_scale
    residual=any(positive) ? maximum(abs,reaction[positive])/flow_scale : 0.
    complementarity=maximum(abs,(P[free]/a.pressure_scale).*(reaction/flow_scale))
    residual<1e-10 || error("Q4 pressure equations did not converge")
    if closure==:swift
        minimum(P)>=0 && minimum(reaction)>=-1e-10*flow_scale && complementarity<1e-10 ||
            error("Q4 pressure obstacle failed complementarity")
    end
    return (;a,pressure_pa=P,midplane_pressure_pa=q4_plane(a,P),iterations,
        positive_pressure_relative_residual=residual,complementarity)
end

function q4_periodic_sample(theta,pressure,query)
    length(theta)==length(pressure) || throw(ArgumentError("sample size mismatch"))
    dt=2pi/length(theta)
    return [begin
        coordinate=mod(t-theta[1],2pi)/dt
        lower=floor(Int,coordinate)+1;fraction=coordinate-floor(coordinate)
        (1-fraction)*pressure[lower]+fraction*pressure[mod1(lower+1,length(pressure))]
    end for t in query]
end
