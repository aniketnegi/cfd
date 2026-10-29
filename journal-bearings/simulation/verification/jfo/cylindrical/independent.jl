module CylindricalJFOReference
using LinearAlgebra, SparseArrays

# Nodal expanded finite differences, not production geometry or face coefficients.
# P is pressure ABOVE cavity, scaled by mu*|omega|*R^2/c^2; V=1-F.
function assemble(aspect,ecc,pv0,nt,nz;direction=1,source=(t,z)->0.)
    nt>=8 && nz>=4 || throw(ArgumentError("FD grid too small"))
    direction in (-1,1) || throw(ArgumentError("direction must be +/-1"))
    dt,dz=2pi/nt,2aspect/nz
    theta=[(i-1)*dt for i in 1:nt]
    zeta=[-aspect+j*dz for j in 0:nz]
    H=1 .-ecc.*cos.(theta);Ht=ecc.*sin.(theta)
    left=-1/dt^2 .+1.5.*Ht./H./dt
    right=-1/dt^2 .-1.5.*Ht./H./dt
    maximum(left)<=0 && maximum(right)<=0 || throw(ArgumentError("FD angular drift is underresolved"))
    diag=2/dt^2+2/dz^2;axial=-1/dz^2
    up=[mod1(i-direction,nt) for i in 1:nt]
    vd=-6 ./ H.^2 ./ dt
    vu=6 .* H[up] ./ H.^3 ./ dt
    forcing=[-6direction*Ht[i]/H[i]^3+12source(theta[i],zeta[j+1])/H[i]^3 for i in 1:nt,j in 1:nz-1]
    b=vec(copy(forcing));n=length(b)
    rows=Int[];cols=Int[];vals=Float64[];br=Int[];bc=Int[];bv=Float64[]
    for j in 1:nz-1,i in 1:nt
        k=i+nt*(j-1)
        for (column,value) in ((k,diag),(mod1(i-1,nt)+nt*(j-1),left[i]),(mod1(i+1,nt)+nt*(j-1),right[i]))
            push!(rows,k);push!(cols,column);push!(vals,value)
        end
        for nextj in (j-1,j+1)
            if 1<=nextj<=nz-1
                push!(rows,k);push!(cols,i+nt*(nextj-1));push!(vals,axial)
            else
                b[k]-=axial*(-pv0)
            end
        end
        push!(br,k);push!(bc,k);push!(bv,vd[i])
        push!(br,k);push!(bc,up[i]+nt*(j-1));push!(bv,vu[i])
    end
    return (;A=sparse(rows,cols,vals,n,n),B=sparse(br,bc,bv,n,n),b,forcing,
        nt,nz,aspect,ecc,pv0,theta,zeta,H,Ht,left,right,diag,axial,up,vd,vu,dt,dz,direction)
end

function balance(a,P,V)
    residual=0.
    @inbounds for j in 2:a.nz,i in 1:a.nt
        l,r=mod1(i-1,a.nt),mod1(i+1,a.nt)
        diff=(2P[i,j]-P[l,j]-P[r,j])/a.dt^2-1.5a.Ht[i]/a.H[i]*(P[r,j]-P[l,j])/a.dt+
            (2P[i,j]-P[i,j-1]-P[i,j+1])/a.dz^2
        transport=6*(a.H[a.up[i]]*V[a.up[i],j]-a.H[i]*V[i,j])/(a.H[i]^3*a.dt)
        residual=max(residual,abs(diff+transport-a.forcing[i,j-1])*a.H[i]^3)
    end
    return residual*a.dt/(6maximum(a.H))
end

# A scalar projected SOR in X=P-V; P=max(X,0), V=max(-X,0).
# No production active-set solver or production sparse factorization is used.
function solve(a;relaxation=1.,tolerance=1e-10,max_sweeps=200000)
    0<relaxation<2 || throw(ArgumentError("invalid relaxation"))
    P=fill(-a.pv0,a.nt,a.nz+1);V=zeros(size(P))
    order=a.direction>0 ? (1:a.nt) : (a.nt:-1:1)
    residual=Inf
    for sweep in 1:max_sweeps
        @inbounds for j in 2:a.nz,i in order
            l,r=mod1(i-1,a.nt),mod1(i+1,a.nt)
            rhs=a.forcing[i,j-1]-a.left[i]*P[l,j]-a.right[i]*P[r,j]-
                a.axial*(P[i,j-1]+P[i,j+1])-a.vu[i]*V[a.up[i],j]
            candidate=rhs>=0 ? rhs/a.diag : -rhs/a.vd[i]
            x=(1-relaxation)*(P[i,j]-V[i,j])+relaxation*candidate
            P[i,j]=max(x,0.);V[i,j]=max(-x,0.)
        end
        if sweep%20==0
            yield()
            residual=balance(a,P,V)
            if residual<tolerance
                return (;a,P,F=1 .-V,sweeps=sweep,residual,relaxation,tolerance)
            end
        end
    end
    error("Independent projected SOR did not converge: residual=$residual")
end

function observables(r)
    a=r.a
    pressure=r.P .+a.pv0
    fx=fy=0.
    for j in 2:a.nz,i in 1:a.nt
        fx+=pressure[i,j]*cos(a.theta[i])*a.dt*a.dz/(2a.aspect)
        fy+=pressure[i,j]*sin(a.theta[i])*a.dt*a.dz/(2a.aspect)
    end
    lower=(-3 .* r.P[:,1] .+ 4 .* r.P[:,2] .- r.P[:,3]) ./ (2a.dz)
    upper=(3 .* r.P[:,end] .- 4 .* r.P[:,end-1] .+ r.P[:,end-2]) ./ (2a.dz)
    out=[a.H.^3 .* lower .* a.dt ./ 12; -a.H.^3 .* upper .* a.dt ./ 12]
    gross=sum(max.(out,0.));incoming=-sum(min.(out,0.))
    return (;P0=hypot(fx,fy),beta_deg=rad2deg(atan(abs(fy),abs(fx))),Q0=gross,
        gross_in_Q0=incoming,net_Q0=sum(out),minimum_fill=minimum(r.F),sweeps=r.sweeps,residual=r.residual)
end
end
