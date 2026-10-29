# Continuous 3D Navier-equation manufactured solution; not K*exact forcing.
function elastic_manufactured(x,youngs_modulus_pa,poisson_ratio)
    E,nu=youngs_modulus_pa,poisson_ratio
    mu=E/(2*(1+nu));lambda=E*nu/((1+nu)*(1-2nu))
    sx,sy,sz=sinpi.(x);cx,cy,cz=cospi.(x);amplitude=1e-6
    u=[amplitude*sx*sy*sz,0.,0.]
    force=amplitude*pi^2 .* [(lambda+4mu)*sx*sy*sz,
        -(lambda+mu)*cx*cy*sz,-(lambda+mu)*cx*sy*cz]
    return (;displacement=u,force)
end
function elastic_manufactured_load(a)
    force=zeros(size(a.K,1))
    points,weights=(-sqrt(3/5),0.,sqrt(3/5)),(5/9,8/9,5/9)
    for cell in a.cells
        X=a.nodes[collect(cell),:];localforce=zeros(24)
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights),(w,ww) in zip(points,weights)
            N,dN=hex8_shape(u,v,w);x=X'*N
            f=elastic_manufactured(x,a.youngs_modulus_pa,a.poisson_ratio).force
            weight=det(X'*dN)*wu*wv*ww
            for i in 1:8,k in 1:3;localforce[3*(i-1)+k]+=N[i]*f[k]*weight;end
        end
        force[elastic_dofs(cell)].+=localforce
    end
    return force
end
function elastic_manufactured_error(a,u)
    numerator=denominator=0.
    points,weights=(-sqrt(3/5),0.,sqrt(3/5)),(5/9,8/9,5/9)
    for cell in a.cells
        X=a.nodes[collect(cell),:];U=reshape(u[elastic_dofs(cell)],3,8)
        for (x,wx) in zip(points,weights),(y,wy) in zip(points,weights),(z,wz) in zip(points,weights)
            N,dN=hex8_shape(x,y,z)
            exact=elastic_manufactured(X'*N,a.youngs_modulus_pa,a.poisson_ratio).displacement
            weight=det(X'*dN)*wx*wy*wz
            numerator+=sum(abs2,U*N-exact)*weight;denominator+=sum(abs2,exact)*weight
        end
    end
    return sqrt(numerator/denominator)
end
