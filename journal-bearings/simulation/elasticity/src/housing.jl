# Reference housing volume with a matching inner fluid-land mesh.
# Caller loads bearing.jl and the public surface mesh constructors first.
isdefined(@__MODULE__,:assemble_elasticity) || include(joinpath(@__DIR__,"linear_elasticity.jl"))

function housing_solid_mesh(grid;wall_thickness_m::Real,n_radial::Int,
        feed_diameter_m=nothing,feed_z_m=nothing)
    t=Float64(wall_thickness_m)
    isfinite(t) && t>0 && n_radial>=1 || throw(ArgumentError("positive wall thickness and radial count required"))
    grid.acceptance.accepted || throw(ArgumentError("unaccepted interface mesh"))
    fitted=hasproperty(grid,:feed_edges)
    fitted && (isnothing(feed_diameter_m) || isnothing(feed_z_m)) &&
        throw(ArgumentError("fitted solid requires the physical drilling diameter and axial position"))
    if fitted
        isfinite(feed_diameter_m) && feed_diameter_m>0 && isfinite(feed_z_m) || throw(ArgumentError("invalid drilling geometry"))
    else
        isnothing(feed_diameter_m) && isnothing(feed_z_m) || throw(ArgumentError("unported solid has no drilling"))
    end
    n=length(grid.nodes);nodes=zeros((n_radial+1)*n,3)
    radius=fitted ? Float64(feed_diameter_m)/2 : 0.
    for layer in 0:n_radial,i in 1:n
        x,y,z=grid.nodes[i].*1e-3
        if layer==0
            nodes[i,:]=[x,y,z]
            continue
        end
        r=hypot(x,y);d=t*layer/n_radial
        X,Y=x*(1+d/r),y*(1+d/r)
        if fitted && y>0
            distance=hypot(x,z-feed_z_m)
            blend=clamp((distance-radius)/max(2radius,2t),0.,1.)
            weight=(1-blend)^2*(1+2blend)
            if weight>0
                # Rim columns retain x,z: a straight +Y polygonal drilling,
                # not a radially flaring port. The surrounding map is smooth.
                X=x*(1+d/r*(1-weight));Y=sqrt((r+d)^2-X^2)
            end
        end
        nodes[i+layer*n,:]=[X,Y,z]
    end
    cells=NTuple{8,Int}[]
    for layer in 0:n_radial-1,quad in grid.quads
        push!(cells,Tuple([collect(quad).+layer*n;collect(quad).+(layer+1)*n]))
    end
    if fitted
        rim=unique(vcat(first.(grid.feed_edges),last.(grid.feed_edges)))
        for node in rim
            x,_,z=nodes[node,:]
            abs(hypot(x,z-feed_z_m)-radius)<1e-9 || throw(ArgumentError("drilling geometry disagrees with interface rim"))
            for layer in 1:n_radial
                nodes[node+layer*n,1]==x && nodes[node+layer*n,3]==z || error("drill column is not parallel to +Y")
            end
        end
    end
    return (;nodes,cells,grid,n_radial,wall_thickness_m=t,base_nodes=n,
        outer_nodes=collect(n*n_radial+1:n*(n_radial+1)),
        feed_diameter_m,feed_z_m)
end

function housing_interface_chart(grid,cell)
    hasproperty(grid,:surface) && return grid.surface.localnodes[cell]
    nt=grid.n_theta;i,j=mod1(cell,nt),(cell-1)÷nt+1
    theta=grid.theta_nodes;z=grid.z_nodes.*1e-3
    R=grid.parameters.housing_midpoint_radius*1e-3;ca=cosd(grid.parameters.semicone_angle_deg)
    low=theta[i];high=i==nt ? theta[1]+2pi : theta[i+1]
    return [R*low z[j]/ca;R*high z[j]/ca;R*high z[j+1]/ca;R*low z[j+1]/ca]
end
function interface_quad_shape(u,v)
    N=[(1-u)*(1-v),(1+u)*(1-v),(1+u)*(1+v),(1-u)*(1+v)]/4
    du=[-(1-v),1-v,1+v,-(1+v)]/4
    dv=[-(1-u),-(1+u),1+u,1-u]/4
    return N,du,dv
end

# Consistent pressure work on the nominal curved bore. Partition of unity
# preserves the integrated fluid pressure force. Structural geometry is Hex8;
# its geometric torque approximation is checked separately under refinement.
function housing_pressure_operator(mesh)
    grid=mesh.grid;p=grid.parameters
    R=p.housing_midpoint_radius*1e-3;ca=cosd(p.semicone_angle_deg)
    rule=eigen(SymTridiagonal(zeros(4),[k/sqrt(4k^2-1) for k in 1:3]))
    points,weights=rule.values,2 .*rule.vectors[1,:].^2
    rows,cols,vals=Int[],Int[],Float64[];areas=zeros(length(grid.quads))
    for (cell,quad) in enumerate(grid.quads)
        xy=housing_interface_chart(grid,cell)
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
            N,du,dv=interface_quad_shape(u,v)
            x,xu,xv=xy'*N,xy'*du,xy'*dv
            jac=xu[1]*xv[2]-xu[2]*xv[1]
            jac>0 || error("inverted pressure transfer quad")
            theta,z=x[1]/R,x[2]*ca
            rh=housing_radius(p,z*1e3)*1e-3
            weight=rh/R*jac*wu*wv
            normal=[ca*cos(theta),ca*sin(theta),sind(p.semicone_angle_deg)]
            areas[cell]+=weight
            for i in 1:4,k in 1:3
                push!(rows,3*(quad[i]-1)+k);push!(cols,cell);push!(vals,N[i]*normal[k]*weight)
            end
        end
    end
    return (;matrix=sparse(rows,cols,vals,3size(mesh.nodes,1),length(grid.quads)),areas_m2=areas)
end

# Reservoir pressure acts on the polygonal drilling wall; passage losses and
# pipe mechanics are not modeled. Returns a load per pascal of gauge pressure.
function housing_passage_load(mesh)
    load=zeros(3size(mesh.nodes,1));grid=mesh.grid
    hasproperty(grid,:feed_edges) || return load
    s=grid.surface;n=mesh.base_nodes;points=(-1/sqrt(3),1/sqrt(3))
    for face in findall(==(:feed),s.boundary)
        quad=s.quads[s.owner[face]];edge=s.edge_nodes[face]
        index=findfirst(i->minmax(quad[i],quad[mod1(i+1,4)])==edge,1:4)
        a,b=quad[index],quad[mod1(index+1,4)]
        for layer in 0:mesh.n_radial-1
            ids=[a+layer*n,b+layer*n,b+(layer+1)*n,a+(layer+1)*n]
            X=mesh.nodes[ids,:]
            for u in points,v in points
                N,du,dv=interface_quad_shape(u,v)
                force=-cross(X'*du,X'*dv)
                for i in 1:4,k in 1:3;load[3*(ids[i]-1)+k]+=N[i]*force[k];end
            end
        end
    end
    return load
end

function assemble_housing(grid;youngs_modulus_pa,poisson_ratio,wall_thickness_m,
        n_radial,feed_diameter_m=nothing,feed_z_m=nothing)
    mesh=housing_solid_mesh(grid;wall_thickness_m,n_radial,feed_diameter_m,feed_z_m)
    a=assemble_elasticity(mesh.nodes,mesh.cells;youngs_modulus_pa,poisson_ratio)
    constrained=constrain_elasticity(a,elastic_dofs(mesh.outer_nodes))
    pressure=housing_pressure_operator(mesh)
    return (;constrained...,mesh,pressure_operator=pressure.matrix,
        interface_areas_m2=pressure.areas_m2,passage_load_per_pa=housing_passage_load(mesh))
end
function solve_housing(a,pressure_pa;passage_pressure_pa=0.)
    length(pressure_pa)==length(a.mesh.grid.quads) || throw(ArgumentError("pressure must have one value per fluid cell"))
    isfinite(passage_pressure_pa) || throw(ArgumentError("invalid passage pressure"))
    land_load=a.pressure_operator*pressure_pa
    load=land_load+a.passage_load_per_pa*passage_pressure_pa
    state=solve_elasticity(a,load)
    return (;state...,load_n=load,land_load_n=land_load,
        normal_displacement_cell_mean_m=(a.pressure_operator'*state.displacement_m)./a.interface_areas_m2)
end
