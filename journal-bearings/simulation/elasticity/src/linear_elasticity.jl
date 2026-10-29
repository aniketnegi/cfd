# Small-strain, homogeneous isotropic 3D elasticity. All quantities use SI.
using LinearAlgebra,SparseArrays

function isotropic_stiffness(youngs_modulus_pa,poisson_ratio)
    E,nu=Float64(youngs_modulus_pa),Float64(poisson_ratio)
    isfinite(E) && E>0 || throw(ArgumentError("Young's modulus must be positive and finite"))
    isfinite(nu) && -1<nu<.5 || throw(ArgumentError("Poisson ratio must lie in (-1,.5)"))
    mu=E/(2*(1+nu));lambda=E*nu/((1+nu)*(1-2nu))
    C=zeros(6,6)
    C[1:3,1:3].=lambda
    for i in 1:3;C[i,i]+=2mu;C[i+3,i+3]=mu;end
    return C
end

# Standard Hex8 order, with the first four nodes on w=-1.
const HEX8_SIGNS=((-1.,-1.,-1.),(1.,-1.,-1.),(1.,1.,-1.),(-1.,1.,-1.),
                  (-1.,-1.,1.),(1.,-1.,1.),(1.,1.,1.),(-1.,1.,1.))
function hex8_shape(u,v,w)
    N=zeros(8);dN=zeros(8,3)
    for (i,(a,b,c)) in enumerate(HEX8_SIGNS)
        N[i]=(1+a*u)*(1+b*v)*(1+c*w)/8
        dN[i,1]=a*(1+b*v)*(1+c*w)/8
        dN[i,2]=b*(1+a*u)*(1+c*w)/8
        dN[i,3]=c*(1+a*u)*(1+b*v)/8
    end
    return N,dN
end
function elastic_strain_matrix(gradient)
    B=zeros(6,24)
    for i in 1:8
        x,y,z=gradient[i,:];k=3*(i-1)
        B[1,k+1]=x;B[2,k+2]=y;B[3,k+3]=z
        B[4,k+1]=y;B[4,k+2]=x
        B[5,k+2]=z;B[5,k+3]=y
        B[6,k+1]=z;B[6,k+3]=x
    end
    return B
end
elastic_dofs(nodes)=[3*(node-1)+component for node in nodes for component in 1:3]

function assemble_elasticity(nodes,cells;youngs_modulus_pa,poisson_ratio)
    size(nodes,2)==3 && all(isfinite,nodes) || throw(ArgumentError("solid nodes must be finite N×3 coordinates in metres"))
    C=isotropic_stiffness(youngs_modulus_pa,poisson_ratio)
    entries=length(cells)*24^2
    rows,cols,values=Vector{Int}(undef,entries),Vector{Int}(undef,entries),Vector{Float64}(undef,entries)
    cursor=0;minimum_jacobian=Inf;volume=0.
    points=(-1/sqrt(3),1/sqrt(3))
    quadrature=[hex8_shape(u,v,w) for u in points,v in points,w in points]
    for cell in cells
        length(cell)==8 && all(i->1<=i<=size(nodes,1),cell) || throw(ArgumentError("invalid Hex8 connectivity"))
        X=nodes[collect(cell),:];Ke=zeros(24,24)
        for (_,dN) in quadrature
            J=X'*dN;detJ=det(J)
            isfinite(detJ) && detJ>0 || throw(ArgumentError("nonpositive solid Jacobian"))
            minimum_jacobian=min(minimum_jacobian,detJ);volume+=detJ
            B=elastic_strain_matrix(dN/J)
            Ke .+= (B'*C*B).*detJ
        end
        dofs=elastic_dofs(cell)
        for j in 1:24,i in 1:24
            cursor+=1;rows[cursor]=dofs[i];cols[cursor]=dofs[j];values[cursor]=Ke[i,j]
        end
    end
    K=sparse(rows,cols,values,3size(nodes,1),3size(nodes,1))
    return (;K,C,nodes,cells,minimum_jacobian,volume_m3=volume,
        youngs_modulus_pa=Float64(youngs_modulus_pa),poisson_ratio=Float64(poisson_ratio))
end

function constrain_elasticity(a,fixed_dofs)
    fixed=sort!(unique!(collect(Int,fixed_dofs)))
    isempty(fixed) && throw(ArgumentError("elasticity requires displacement supports"))
    all(i->1<=i<=size(a.K,1),fixed) || throw(ArgumentError("invalid fixed displacement DOF"))
    free=setdiff(1:size(a.K,1),fixed)
    isempty(free) && throw(ArgumentError("no free solid degrees of freedom"))
    factor=cholesky(Symmetric(a.K[free,free]))
    return (;a...,fixed,free,factor)
end

function solve_elasticity(a,load;fixed_values=zeros(length(a.fixed)))
    length(load)==size(a.K,1) && all(isfinite,load) || throw(ArgumentError("invalid nodal load"))
    length(fixed_values)==length(a.fixed) && all(isfinite,fixed_values) || throw(ArgumentError("invalid support displacement"))
    u=zeros(length(load));u[a.fixed]=fixed_values
    rhs=load[a.free]-a.K[a.free,a.fixed]*fixed_values
    u[a.free]=a.factor\rhs
    all(isfinite,u) || error("nonfinite elastic displacement")
    residual=a.K*u-load
    scale=max(norm(load,Inf),norm(a.K*u,Inf),1e-30)
    relative=norm(residual[a.free],Inf)/scale
    relative<1e-9 || error("elastic equilibrium residual exceeds tolerance: $relative")
    reaction=zeros(length(load));reaction[a.fixed]=residual[a.fixed]
    strain_energy=dot(u,a.K*u)/2
    return (;displacement_m=u,reaction_n=reaction,relative_residual=relative,strain_energy_j=strain_energy)
end

function elastic_stress(a,displacement)
    length(displacement)==size(a.K,1) || throw(ArgumentError("displacement size mismatch"))
    points=(-1/sqrt(3),1/sqrt(3));maximum_von_mises=0.;maximum_strain=0.
    for cell in a.cells
        X=a.nodes[collect(cell),:];u=displacement[elastic_dofs(cell)]
        for x in points,y in points,z in points
            _,dN=hex8_shape(x,y,z)
            strain=elastic_strain_matrix(dN/(X'*dN))*u;stress=a.C*strain
            sx,sy,sz,txy,tyz,txz=stress
            vm=sqrt(((sx-sy)^2+(sy-sz)^2+(sz-sx)^2)/2+3*(txy^2+tyz^2+txz^2))
            maximum_von_mises=max(maximum_von_mises,vm)
            maximum_strain=max(maximum_strain,norm(strain,Inf))
        end
    end
    return (;maximum_von_mises_pa=maximum_von_mises,maximum_engineering_strain=maximum_strain)
end
