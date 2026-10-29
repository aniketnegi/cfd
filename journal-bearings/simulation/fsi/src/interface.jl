# Continuous Q1 housing displacement sampled at retained FV cell/face points.
# Positive displacement along the nominal housing normal opens the film.
struct HousingGapField
    index::Dict{Tuple{Float64,Float64},Int}
    values::Vector{Float64}
end
function (field::HousingGapField)(theta,z)
    key=(mod(Float64(theta),2pi),Float64(z))
    i=get(field.index,key,0)
    i>0 || throw(ArgumentError("gap query is not a retained point of this interface mesh"))
    return field.values[i]
end

function housing_quad_weights(xy,point)
    uv=zeros(2)
    for _ in 1:20
        N,du,dv=interface_quad_shape(uv...)
        residual=xy'*N-point
        correction=hcat(xy'*du,xy'*dv)\residual
        uv .-=correction
        if norm(correction,Inf)<1e-12
            maximum(abs,uv)<=1+1e-8 || throw(ArgumentError("interface query lies outside its owner"))
            return first(interface_quad_shape(uv...))
        end
    end
    error("interface coordinate inversion did not converge")
end
function housing_gap_operator(solid)
    grid=solid.mesh.grid;p=grid.parameters
    R=p.housing_midpoint_radius*1e-3;ca=cosd(p.semicone_angle_deg);period=2pi*R
    index=Dict{Tuple{Float64,Float64},Int}()
    rows,cols,vals=Int[],Int[],Float64[]
    points=Tuple{Float64,Float64}[]
    function add(cell,theta,z)
        key=(mod(Float64(theta),2pi),Float64(z))
        haskey(index,key) && return
        row=length(index)+1;index[key]=row;push!(points,key)
        xy=housing_interface_chart(grid,cell)
        centre=sum(xy[:,1])/4
        x=centre+mod(R*theta-centre+period/2,period)-period/2
        N=housing_quad_weights(xy,[x,z/ca])
        normal=[ca*cos(theta),ca*sin(theta),sind(p.semicone_angle_deg)]
        for (node,w) in zip(grid.quads[cell],N),k in 1:3
            push!(rows,row);push!(cols,3*(node-1)+k);push!(vals,w*normal[k])
        end
    end
    if hasproperty(grid,:surface)
        s=grid.surface
        for i in axes(s.centre,1);add(i,s.centre[i,1]/s.R,s.centre[i,2]*s.ca);end
        for f in eachindex(s.owner);add(s.owner[f],s.mid[f][1]/s.R,s.mid[f][2]*s.ca);end
    else
        a=full_film_axes(grid);nt,nz=grid.n_theta,grid.n_axial
        for j in 1:nz,i in 1:nt
            add(i+nt*(j-1),a.theta_centres[i],a.z_centres[j])
            add(i+nt*(j-1),a.theta_edges[i+1],a.z_centres[j])
        end
        for j in 1:nz+1,i in 1:nt
            add(i+nt*(j==1 ? 0 : j-2),a.theta_centres[i],a.z_edges[j])
        end
    end
    return (;matrix=sparse(rows,cols,vals,length(index),size(solid.K,1)),index,points)
end
housing_gap_field(transfer,u)=HousingGapField(transfer.index,transfer.matrix*u)
