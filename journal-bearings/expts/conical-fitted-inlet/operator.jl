# Experimental conservative multi-point surface FV. Production solver untouched.
include(joinpath(@__DIR__,"mesh.jl"))
include(joinpath(@__DIR__,"../../simulation/reynolds/conical_equilibrium.jl"))

function fitted_operator(m,pose;rpm=2000.,viscosity=.0277,unit_mobility=false,
        boundary_pressure=(x,y,tag)->tag==:feed ? 5e5 : 0.,
        incoming_fill=(x,y,tag)->1.)
    N=length(m.quads);F=length(m.owner)
    # x=Rref*theta, y=z/cos(alpha): the cone diffusion tensor is
    # K*diag(Rref/rh,rh/Rref), and entrainment vector is (U*h/2,0).
    wrap(x)=mod(x+m.period/2,m.period)-m.period/2
    adjacent=[Int[] for _ in 1:N]
    di=Int[];dj=Int[];dv=Float64[]
    for f in 1:F
        k,l=m.owner[f],m.neighbour[f];push!(adjacent[k],f)
        push!(di,k);push!(dj,f);push!(dv,1.)
        if l>0
            push!(adjacent[l],f);push!(di,l);push!(dj,f);push!(dv,-1.)
        end
    end
    divergence=sparse(di,dj,dv,N,F)
    gi=Int[];gj=Int[];gx=Float64[];gy=Float64[];g0=zeros(N,2)
    pb=[m.neighbour[f]==0 ? boundary_pressure(m.mid[f]...,m.boundary[f]) : 0. for f in 1:F]
    for k in 1:N
        vectors=Vector{Float64}[];others=Int[];boundaryvalues=Float64[]
        for f in adjacent[k]
            l=m.owner[f]==k ? m.neighbour[f] : m.owner[f]
            d=l>0 ? m.centre[l,:]-m.centre[k,:] : m.mid[f]-m.centre[k,:]
            d[1]=wrap(d[1]);push!(vectors,d);push!(others,l);push!(boundaryvalues,pb[f])
        end
        T=sum(d*d'/dot(d,d) for d in vectors)
        minimum(eigvals(Symmetric(T)))>1e-8 || error("singular gradient stencil")
        for (d,l,p) in zip(vectors,others,boundaryvalues)
            a=T\(d/dot(d,d))
            push!(gi,k);push!(gj,k);push!(gx,-a[1]);push!(gy,-a[2])
            if l>0
                push!(gi,k);push!(gj,l);push!(gx,a[1]);push!(gy,a[2])
            else
                g0[k,:].+=a*p
            end
        end
    end
    Gx=sparse(gi,gj,gx,N,N);Gy=sparse(gi,gj,gy,N,N)
    ti=Int[];tj=Int[];tv=Float64[];hi=Int[];hj=Int[];hv=Float64[]
    ei=Int[];ej=Int[];ev=Float64[];cx=zeros(F);cy=zeros(F);constant=zeros(F);adv=zeros(F)
    gap(x,y)=conical_normal_gap(pose.parameters,x/m.R-pose.eccentricity_angle_rad,y*m.ca)
    for f in 1:F
        k,l=m.owner[f],m.neighbour[f];x=m.mid[f];n=m.normals[f];length=m.lengths[f]
        rj,rh=conical_radii(pose.parameters,x[2]*m.ca)
        h=gap(x...);K=unit_mobility ? 1. : h^3/(12viscosity)
        D=[m.R/rh,rh/m.R];d=l>0 ? m.centre[l,:]-m.centre[k,:] : x-m.centre[k,:];d[1]=wrap(d[1])
        dot(n,d)>0 || error("nonpositive pressure connection")
        a=dot(n,D.*n)/dot(n,d)
        g=length*K*a
        push!(ti,f);push!(tj,k);push!(tv,g)
        if l>0
            push!(ti,f);push!(tj,l);push!(tv,-g)
        else
            constant[f]=-g*pb[f]
        end
        correction=-length*K*(D.*n-a*d);cx[f],cy[f]=correction
        push!(ei,f);push!(ej,k);push!(ev,l>0 ? .5 : 1.)
        if l>0;push!(ei,f);push!(ej,l);push!(ev,.5);end
        adv[f]=unit_mobility ? 0. : rpm*2pi/60*rj*h/2*n[1]*length
        upstream=adv[f]>=0 ? k : l
        if upstream>0
            push!(hi,f);push!(hj,upstream);push!(hv,adv[f])
        else
            adv[f]*=incoming_fill(x...,m.boundary[f])
        end
    end
    E=sparse(ei,ej,ev,F,N)
    pressure_flux=sparse(ti,tj,tv,F,N)+spdiagm(0=>cx)*E*Gx+spdiagm(0=>cy)*E*Gy
    constant.+=cx.*(E*g0[:,1])+cy.*(E*g0[:,2])
    transport=sparse(hi,hj,hv,F,N)
    A=divergence*pressure_flux;B=-divergence*transport;b=-divergence*(constant+adv)
    h=[gap(m.centre[k,:]...) for k in 1:N]
    return (;A,B,b,pressure_flux,transport,constant,adv,divergence,h,Gx,Gy,g0)
end

function fitted_force(m,pressure)
    force=zeros(3);area=0.
    for (k,xy) in enumerate(m.localnodes),u in (-1/sqrt(3),1/sqrt(3)),v in (-1/sqrt(3),1/sqrt(3))
        shape=[(1-u)*(1-v),(1+u)*(1-v),(1+u)*(1+v),(1-u)*(1+v)]/4
        du=[-(1-v),1-v,1+v,-(1+v)]/4;dv=[-(1-u),-(1+u),1+u,1-u]/4
        x=xy'*shape;Xu=xy'*du;Xv=xy'*dv;J=Xu[1]*Xv[2]-Xu[2]*Xv[1]
        J>0 || error("inverted bilinear quadrature mapping")
        theta=x[1]/m.R;z=x[2]*m.ca;rh=m.R+(.05-z)*tand(10.)
        area+=rh/m.R*J
        force.+=pressure[k]*rh/m.R*J.*[m.ca*cos(theta),m.ca*sin(theta),sind(10.)]
    end
    return (;force,area)
end

function solve_fitted(m,pose;rpm=2000.,max_iterations=150)
    a=fitted_operator(m,pose;rpm)
    scale=max(5e5,12*.0277*abs(rpm*2pi/60)*.05*m.period/minimum(a.h)^2)
    s=solve_jfo_active_set(a.A,a.B,a.b,scale;max_iterations)
    pressure=s.pressure;fill=1 .- s.void
    flux=a.pressure_flux*pressure+a.constant+a.adv-a.transport*s.void
    residual=a.divergence*flux
    feed=m.boundary.==:feed;large=m.boundary.==:large;small=m.boundary.==:small
    supply=-sum(flux[feed]);outward=[sum(flux[large]),sum(flux[small])]
    global_scale=sum(abs,flux[feed])+sum(abs,flux[large])+sum(abs,flux[small])
    local_error=maximum(abs,residual)/maximum(abs,flux)
    imbalance=abs(sum(outward)-supply)/global_scale
    minimum(pressure)>=-1e-10*scale && minimum(fill)>=-1e-10 && maximum(fill)<=1+1e-10 || error("invalid pressure/content bounds")
    local_error<1e-9 && imbalance<1e-9 || error("failed conservative flux gate")
    force=fitted_force(m,pressure).force
    return (;pressure,fill,force,iterations=s.iterations,net_feed=supply,
        injection=-sum(min.(flux[feed],0.)),return_flow=sum(max.(flux[feed],0.)),outward,
        local_error,imbalance,peak_pressure=maximum(pressure),minimum_fill=minimum(fill))
end
