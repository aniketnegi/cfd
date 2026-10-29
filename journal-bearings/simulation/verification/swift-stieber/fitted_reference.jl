# Verification-only analytical fields and a separate algebraic LCP iteration.
using LinearAlgebra, SparseArrays

function fitted_gauss(n)
    rule=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k^2-1) for k in 1:n-1]))
    return rule.values,2 .* rule.vectors[1,:].^2
end

function fitted_mms(theta,z,angle,rpm;front=.37)
    L,c,e,mu,P=.1,50e-6,30e-6,.0277,1e7
    ca,k=cosd(angle),tand(angle)
    rj=.05+(.05-z)*k;r=rj+c
    st,ct=sincos(theta-.37)
    root=sqrt(rj^2-e^2*st^2)
    h=ca*(c-e*ct+e^2*st^2/(rj+root))
    ht,hz=ca*(e*st+e^2*st*ct/root),ca*k*(rj/root-1)
    K,Kt,Kz=h^3/(12mu),h^2*ht/(4mu),h^2*hz/(4mu)
    t=mod(theta-front,2pi);wet=0<t<pi
    s,q=sincos(t)
    g=wet ? s^4 : 0.;gt=wet ? 4s^3*q : 0.;gtt=wet ? 12s^2*q^2-4s^4 : 0.
    f,df=z*(L-z)/L^2,(L-2z)/L^2
    pressure=P*g*f
    pt,ptt,pz,pzz=P*gt*f,P*gtt*f,P*g*df,-2P*g/L^2
    diffusion=(Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
    wedge=rpm*2pi/60*rj*ht/(2ca)
    reaction=wet ? 0. : .001*(1+.2cos(theta))
    return (;pressure,pt,pz,h,source=(wedge-diffusion)*ca/.05005-reaction,reaction)
end

function fitted_exact_ends(angle,front;n=160)
    points,weights=fitted_gauss(n)
    flow=zeros(2)
    for (j,z) in enumerate((0.,.1)),(x,w) in zip(points,weights)
        theta=front+(x+1)*pi/2
        ref=fitted_mms(theta,z,angle,0.;front)
        r=.05005+(.05-z)*tand(angle)
        flow[j]+=w*pi/2*ref.h^3/(12*.0277)*r*cosd(angle)*1e7*sin(theta-front)^4/.1
    end
    return flow
end

# Recover an owner's oriented chart edge, rather than reusing its cached normal.
function fitted_reference_edge(s,face)
    quad=s.quads[s.owner[face]];edge=s.edge_nodes[face]
    i=findfirst(i->minmax(quad[i],quad[mod1(i+1,4)])==edge,1:4)
    xy=s.localnodes[s.owner[face]]
    return xy[i,:],xy[mod1(i+1,4),:]
end

function fitted_exact_feed(s,angle,rpm,front;n=8)
    points,weights=fitted_gauss(n)
    flow=0.
    for face in findall(==(:feed),s.boundary)
        a,b=fitted_reference_edge(s,face);t=b-a
        for (u,w) in zip(points,weights)
            x=a+(u+1)/2*t
            theta,z=x[1]/s.R,x[2]*s.ca
            ref=fitted_mms(theta,z,angle,rpm;front)
            rj=.05+(.05-z)*tand(angle);rh=rj+50e-6
            K=ref.h^3/(12*.0277)
            qx=rpm*2pi/60*rj*ref.h/2-K/rh*ref.pt
            qy=-K*rh/s.R*s.ca*ref.pz
            flow+=w/2*(qx*t[2]-qy*t[1])
        end
    end
    return flow
end

# Analytic whole-cone force minus an independently integrated opening polygon.
function fitted_exact_force(s,angle,front;n=6)
    radial=1e7*.05005*.1/6
    force=radial.*[-16sin(front)/15,16cos(front)/15,tand(angle)*3pi/8]
    points,weights=fitted_gauss(n)
    centre=[s.R*pi/2,.05/s.ca]
    for face in findall(==(:feed),s.boundary)
        a,b=fitted_reference_edge(s,face)
        first,second=a-centre,b-centre
        jacobian=abs(first[1]*second[2]-first[2]*second[1])
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
            t,q=(u+1)/2,(v+1)/2
            x=centre+t*first+(1-t)*q*second
            theta,z=x[1]/s.R,x[2]*s.ca
            pressure=fitted_mms(theta,z,angle,0.;front).pressure
            radius=.05005+(.05-z)*tand(angle)
            area=jacobian*(1-t)*wu*wv/4*radius/s.R
            force .-= pressure*area.*[s.ca*cos(theta),s.ca*sin(theta),sind(angle)]
        end
    end
    return force
end

"""Projected row-wise SOR for cross-checking the same algebraic LCP.
No production assembly or active-set kernel is called. This checks the nonlinear
solve, not an independently discretized Reynolds equation.
"""
function fitted_projected_sor(A,b;relaxation=1.4,tolerance=1e-12,max_sweeps=100000)
    rows=sparse(transpose(A));diagonal=diag(A)
    all(>(0),diagonal) || error("projected SOR requires positive diagonal entries")
    p=zeros(length(b));residual=Inf
    for sweep in 1:max_sweeps
        for i in eachindex(p)
            rhs=b[i]
            for index in nzrange(rows,i)
                j=rows.rowval[index]
                j==i || (rhs-=rows.nzval[index]*p[j])
            end
            p[i]=max(0.,(1-relaxation)*p[i]+relaxation*rhs/diagonal[i])
        end
        if sweep%20==0
            slack=(A*p-b)./diagonal
            residual=norm(min.(p,slack),Inf)/max(norm(p,Inf),1e-30)
            residual<tolerance && return (;pressure=p,sweeps=sweep,residual)
        end
    end
    error("projected SOR failed after $max_sweeps sweeps; residual=$residual")
end
