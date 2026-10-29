# Verification mathematics only; no production assembly or solver includes.
module ConicalTransportReference
using LinearAlgebra,SparseArrays

function geometry(theta,z,angle,ecc;pose=0.,mu=.0277)
    ca,k=cosd(angle),tand(angle)
    rj=.05+(.05-z)*k;rh=rj+50e-6;e=ecc*50e-6
    st,ct=sincos(theta-pose);d=sqrt(rj^2-e^2*st^2)
    h=ca*(50e-6-e*ct+e^2*st^2/(rj+d))
    ht=ca*(e*st+e^2*st*ct/d);hz=ca*k*(rj/d-1)
    K=h^3/(12mu);Kt=h^2*ht/(4mu);Kz=h^2*hz/(4mu)
    return (;ca,k,rj,rh,h,ht,hz,K,Kt,Kz)
end

# Smooth complementary fields with a curved front; the +Y inlet is fully wet.
function manufactured(theta,z,angle,rpm;phase=pi/2)
    g=geometry(theta,z,angle,.6;pose=.37)
    t=cos(theta-phase)+.25cos(2pi*z/.1)
    tt,ttt=-sin(theta-phase),-cos(theta-phase)
    tz,tzz=-.25*2pi/.1*sin(2pi*z/.1),-.25*(2pi/.1)^2*cos(2pi*z/.1)
    f,df,ddf=z*(.1-z)/.1^2,(.1-2z)/.1^2,-2/.1^2
    P0=1e7
    P=pt=ptt=pz=pzz=Ft=0.;F=1.
    if t>0
        P=P0*f*t^4;pt=P0*f*4t^3*tt
        ptt=P0*f*(12t^2*tt^2+4t^3*ttt)
        pz=P0*(df*t^4+f*4t^3*tz)
        pzz=P0*(ddf*t^4+8df*t^3*tz+f*(12t^2*tz^2+4t^3*tzz))
    else
        F=1-.2*t^4;Ft=-.8*t^3*tt
    end
    omega=rpm*2pi/60
    diffusion=(g.Kt*pt+g.K*ptt)/(g.rh*g.ca)+g.ca*((g.Kz*g.rh-g.K*g.k)*pz+g.K*g.rh*pzz)
    transport=omega*g.rj*(g.ht*F+g.h*Ft)/(2g.ca)
    source=(transport-diffusion)*g.ca/.05005
    qx=omega*g.rj*g.h/2*F-g.K/g.rh*pt
    qy=-g.K*g.rh/.05005*g.ca*pz
    return (;pressure=P,fill=F,pt,pz,source,qx,qy)
end

# Expanded nodal theta-z finite differences. Pressure is scaled, content is not.
function assemble(angle,ecc,nt,nz;rpm=2000.,ends=(3e5,3e5),pose=3pi/2)
    nt>=8 && nz>=4 || throw(ArgumentError("FD grid too small"))
    rpm!=0 || throw(ArgumentError("rotating FD reference requires nonzero speed"))
    dt,dz=2pi/nt,.1/nz;direction=sign(rpm)
    scale=.0277*abs(rpm*2pi/60)*.05^2/(50e-6)^2
    theta=collect(0:nt-1).*dt;z=collect(0:nz).*dz
    n=nt*(nz-1)
    left,right,low,high,diagonal,vd,vu,forcing=(zeros(nt,nz-1) for _ in 1:8)
    H=[geometry(t,y,angle,ecc;pose).h for t in theta,y in z]
    up=[mod1(i-Int(direction),nt) for i in 1:nt]
    rows,cols,vals=Int[],Int[],Float64[];br,bc,bv=Int[],Int[],Float64[];b=zeros(n)
    add(k,l,v)=(push!(rows,k);push!(cols,l);push!(vals,v))
    for j in 1:nz-1,i in 1:nt
        g=geometry(theta[i],z[j+1],angle,ecc;pose)
        beta=(g.rh*g.ca)^2;drift=beta*(3g.hz/g.h-g.k/g.rh)
        left[i,j]=-1/dt^2+1.5g.ht/(g.h*dt)
        right[i,j]=-1/dt^2-1.5g.ht/(g.h*dt)
        low[i,j]=-beta/dz^2+drift/(2dz)
        high[i,j]=-beta/dz^2-drift/(2dz)
        diagonal[i,j]=2/dt^2+2beta/dz^2
        factor=6*.0277*abs(rpm*2pi/60)*g.rj*g.rh/(g.h^3*scale)
        vd[i,j]=-factor*g.h/dt;vu[i,j]=factor*H[up[i],j+1]/dt
        forcing[i,j]=-direction*factor*g.ht
        k=i+nt*(j-1);b[k]=forcing[i,j]
        add(k,k,diagonal[i,j]);add(k,mod1(i-1,nt)+nt*(j-1),left[i,j]);add(k,mod1(i+1,nt)+nt*(j-1),right[i,j])
        j==1 ? (b[k]-=low[i,j]*ends[1]/scale) : add(k,k-nt,low[i,j])
        j==nz-1 ? (b[k]-=high[i,j]*ends[2]/scale) : add(k,k+nt,high[i,j])
        push!(br,k);push!(bc,k);push!(bv,vd[i,j])
        push!(br,k);push!(bc,up[i]+nt*(j-1));push!(bv,vu[i,j])
    end
    maximum(left)<=0 && maximum(right)<=0 && maximum(low)<=0 && maximum(high)<=0 || error("FD drift underresolved")
    return (;A=sparse(rows,cols,vals,n,n),B=sparse(br,bc,bv,n,n),b,nt,nz,theta,z,dt,dz,scale,
        left,right,low,high,diagonal,vd,vu,forcing,up,ends,direction,angle,ecc,pose)
end

function nodal_residual(a,P,V)
    residual=0.
    @inbounds for j in 1:a.nz-1,i in 1:a.nt
        k=j+1;l,r=mod1(i-1,a.nt),mod1(i+1,a.nt)
        value=a.left[i,j]*(P[l,k]-P[i,k])+a.right[i,j]*(P[r,k]-P[i,k])+
            a.low[i,j]*(P[i,k-1]-P[i,k])+a.high[i,j]*(P[i,k+1]-P[i,k])+
            a.vd[i,j]*V[i,k]+a.vu[i,j]*V[a.up[i],k]-a.forcing[i,j]
        residual=max(residual,abs(value)/max(a.diagonal[i,j],-a.vd[i,j]))
    end
    return residual
end

function solve(a;relaxation=1.,tolerance=1e-11,max_sweeps=200000)
    P=ones(a.nt,a.nz+1);P[:,1].=a.ends[1]/a.scale;P[:,end].=a.ends[2]/a.scale
    V=zeros(size(P));residual=Inf
    order=a.direction>0 ? (1:a.nt) : (a.nt:-1:1)
    for sweep in 1:max_sweeps
        @inbounds for j in 1:a.nz-1,i in order
            k=j+1;l,r=mod1(i-1,a.nt),mod1(i+1,a.nt)
            rhs=a.forcing[i,j]-a.left[i,j]*P[l,k]-a.right[i,j]*P[r,k]-
                a.low[i,j]*P[i,k-1]-a.high[i,j]*P[i,k+1]-a.vu[i,j]*V[a.up[i],k]
            candidate=rhs>=0 ? rhs/a.diagonal[i,j] : -rhs/a.vd[i,j]
            x=(1-relaxation)*(P[i,k]-V[i,k])+relaxation*candidate
            P[i,k]=max(x,0.);V[i,k]=max(-x,0.)
        end
        if sweep%20==0
            yield();residual=nodal_residual(a,P,V)
            residual<tolerance && return (;a,pressure=P.*a.scale,fill=1 .-V,sweeps=sweep,residual)
        end
    end
    error("Independent cone FD failed: residual=$residual")
end

# Same-matrix scalar pressure/void iteration, independent of active-set policy.
function projected_gs(A,B,b,scale;relaxation=1.,tolerance=1e-11,max_sweeps=200000)
    C=A*scale;ar=sparse(transpose(C));br=sparse(transpose(B))
    ad,bd=diag(C),diag(B);all(>(0),ad) && all(<(0),bd) || error("invalid GS diagonals")
    p=zeros(length(b));v=zeros(length(b));residual=Inf
    for sweep in 1:max_sweeps
        for i in eachindex(b)
            rhs=b[i]
            for k in nzrange(ar,i)
                j=ar.rowval[k];j==i || (rhs-=ar.nzval[k]*p[j])
            end
            for k in nzrange(br,i)
                j=br.rowval[k];j==i || (rhs-=br.nzval[k]*v[j])
            end
            candidate=rhs>=0 ? rhs/ad[i] : -rhs/bd[i]
            x=(1-relaxation)*(p[i]-v[i])+relaxation*candidate
            p[i]=max(x,0.);v[i]=max(-x,0.)
        end
        if sweep%20==0
            yield();residual=norm((C*p+B*v-b)./max.(ad,-bd),Inf)
            residual<tolerance && return (;pressure=p.*scale,fill=1 .-v,sweeps=sweep,residual)
        end
    end
    error("Independent same-matrix GS failed: residual=$residual")
end
end
