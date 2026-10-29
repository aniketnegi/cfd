# Linear reconstruction at an actual axial plane of a fitted land mesh.
# Load full_film.jl first. Coordinates are SI, pressures in the boundary reference.
function fitted_pressure_slice(grid,pressure_pa;z_m,boundary_pressure)
    s=fitted_full_film_geometry(grid)
    length(pressure_pa)==length(s.quads) && all(isfinite,pressure_pa) ||
        throw(ArgumentError("pressure must have one finite value per land cell"))
    isfinite(z_m) && 0<z_m<s.L || throw(ArgumentError("sampling plane must lie inside the bearing"))
    target=Float64(z_m)/s.ca
    values=[s.neighbour[i]==0 ? Float64(boundary_pressure(s.mid[i]...,s.boundary[i])) : 0. for i in eachindex(s.owner)]
    all(isfinite,values) || throw(ArgumentError("boundary pressure must be finite"))
    gradient=fitted_pressure_gradients(s,values)
    gx=gradient.x*pressure_pa+gradient.constant[:,1]
    gy=gradient.y*pressure_pa+gradient.constant[:,2]
    cells=Int[];theta=NTuple{2,Float64}[];pressure=NTuple{2,Float64}[]
    wrap(dx)=mod(dx+s.period/2,s.period)-s.period/2
    maximum_pressure=-Inf;peak_theta=NaN
    function consider(x,p)
        if p>maximum_pressure
            maximum_pressure=p;peak_theta=mod(x/s.R,2pi)
        end
    end
    for (cell,xy) in enumerate(s.localnodes)
        # Half-open ownership prevents double counting a shared horizontal edge.
        minimum(xy[:,2])<=target<maximum(xy[:,2]) || continue
        crossings=Float64[]
        for i in 1:4
            j=mod1(i+1,4);y0,y1=xy[i,2],xy[j,2]
            if y0==y1
                y0==target && append!(crossings,(xy[i,1],xy[j,1]))
            elseif min(y0,y1)<=target<=max(y0,y1)
                push!(crossings,xy[i,1]+(target-y0)/(y1-y0)*(xy[j,1]-xy[i,1]))
            end
        end
        length(crossings)>=2 || continue
        low,high=extrema(crossings)
        high>low || continue
        sample(x)=pressure_pa[cell]+gx[cell]*wrap(x-s.centre[cell,1])+gy[cell]*(target-s.centre[cell,2])
        pl,ph=sample(low),sample(high)
        push!(cells,cell);push!(theta,(low/s.R,high/s.R));push!(pressure,(pl,ph))
        consider(low,pl);consider(high,ph)
    end
    # Include the prescribed rim trace where the plane intersects the opening.
    boundary_theta=Float64[];boundary_values=Float64[]
    for face in findall(==(:feed),s.boundary)
        a,b=s.nodes[s.edge_nodes[face][1]],s.nodes[s.edge_nodes[face][2]]
        x0,y0=a;x1,y1=b
        min(y0,y1)<=target<=max(y0,y1) || continue
        xs=y0==y1 ? (x0,x1) : (x0+(target-y0)/(y1-y0)*wrap(x1-x0),)
        for x in xs
            p=Float64(boundary_pressure(mod(x,s.period),target,:feed))
            isfinite(p) || throw(ArgumentError("nonfinite rim trace pressure"))
            push!(boundary_theta,mod(x/s.R,2pi));push!(boundary_values,p)
            consider(x,p)
        end
    end
    isempty(cells) && throw(ArgumentError("sampling plane has no land segments"))
    return (;z_m=Float64(z_m),cells,theta,pressure,boundary_theta,boundary_values,
        maximum_pressure_pa=maximum_pressure,peak_theta_rad=peak_theta)
end
