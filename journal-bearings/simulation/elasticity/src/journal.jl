# Solid nominal journal: square-core O-grid, no axial bore or overhangs.
isdefined(@__MODULE__,:housing_interface_chart) || include(joinpath(@__DIR__,"housing.jl"))

function journal_solid_mesh(p;n_theta::Int,n_axial::Int,n_radial::Int)
    n_theta>=8 && n_theta%4==0 && n_axial>=1 && n_radial>=1 ||
        throw(ArgumentError("journal needs angular count divisible by four and positive axial/radial counts"))
    side=n_theta÷4;origin=-pi/4
    id(i,j)=i+1+j*(side+1)
    cross_nodes=Tuple{Float64,Float64}[]
    for j in 0:side,i in 0:side
        push!(cross_nodes,(.25*(2i/side-1),.25*(2j/side-1)))
    end
    perimeter=[ [id(side,j) for j in 0:side-1];[id(i,side) for i in side:-1:1];
        [id(0,j) for j in side:-1:1];[id(i,0) for i in 0:side-1] ]
    cross_cells=NTuple{4,Int}[]
    for j in 0:side-1,i in 0:side-1
        push!(cross_cells,(id(i,j),id(i+1,j),id(i+1,j+1),id(i,j+1)))
    end
    previous=perimeter
    for layer in 1:n_radial
        ring=Int[];t=layer/n_radial
        for i in 1:n_theta
            theta=origin+2pi*(i-1)/n_theta;x,y=cross_nodes[perimeter[i]]
            push!(cross_nodes,((1-t)*x+t*cos(theta),(1-t)*y+t*sin(theta)))
            push!(ring,length(cross_nodes))
        end
        for i in 1:n_theta
            j=mod1(i+1,n_theta)
            push!(cross_cells,(previous[i],ring[i],ring[j],previous[j]))
        end
        previous=ring
    end
    nc=length(cross_nodes);L=p.bearing_length*1e-3
    nodes=zeros(nc*(n_axial+1),3)
    for j in 0:n_axial
        z=L*j/n_axial;r=journal_radius(p,z*1e3)*1e-3
        r>0 || throw(ArgumentError("journal radius must stay positive"))
        for (i,(x,y)) in enumerate(cross_nodes);nodes[i+j*nc,:]=[r*x,r*y,z];end
    end
    cells=NTuple{8,Int}[]
    for j in 0:n_axial-1,q in cross_cells
        push!(cells,Tuple([collect(q).+j*nc;collect(q).+(j+1)*nc]))
    end
    outer_nodes=[i+j*nc for j in 0:n_axial for i in previous]
    end_nodes=[collect(1:nc);collect(n_axial*nc+1:(n_axial+1)*nc)]
    return (;p,nodes,cells,n_theta,n_axial,n_radial,origin,length_m=L,
        cross_nodes=nc,outer_ring_nodes=previous,outer_nodes,end_nodes)
end

# Continuous outer-surface Q1 interpolation on the independent journal mesh.
function journal_interface_weights(mesh,theta,z)
    isfinite(theta) && isfinite(z) && -1e-12<=z<=mesh.length_m+1e-12 ||
        throw(ArgumentError("journal sample outside its axial span"))
    x=mod(theta-mesh.origin,2pi)*mesh.n_theta/(2pi)
    i=min(floor(Int,x),mesh.n_theta-1);u=x-i
    y=clamp(z,0.,mesh.length_m)*mesh.n_axial/mesh.length_m
    j=min(floor(Int,y),mesh.n_axial-1);v=y-j
    a,b=mesh.outer_ring_nodes[i+1],mesh.outer_ring_nodes[mod1(i+2,mesh.n_theta)]
    n=mesh.cross_nodes
    nodes=(a+j*n,b+j*n,b+(j+1)*n,a+(j+1)*n)
    weights=((1-u)*(1-v),u*(1-v),u*v,(1-u)*v)
    return nodes,weights
end

# Convex polygon clipping in the (theta,z) chart. Split at structural element
# edges before pressure quadrature; no nearest-cell pressure projection.
function journal_clip(poly,axis,bound,keep_above)
    isempty(poly) && return poly
    out=Vector{Vector{Float64}}();previous=poly[end]
    inside(p)=keep_above ? p[axis]>=bound : p[axis]<=bound
    for current in poly
        if inside(current)!=inside(previous)
            t=(bound-previous[axis])/(current[axis]-previous[axis])
            push!(out,previous+t*(current-previous))
        end
        inside(current) && push!(out,current)
        previous=current
    end
    return out
end
function journal_surface_quadrature(mesh,poly,visit)
    dt,dz=2pi/mesh.n_theta,mesh.length_m/mesh.n_axial
    imin=floor(Int,(minimum(p[1] for p in poly)-mesh.origin)/dt)
    imax=floor(Int,(maximum(p[1] for p in poly)-mesh.origin)/dt)
    jmin=clamp(floor(Int,minimum(p[2] for p in poly)/dz),0,mesh.n_axial-1)
    jmax=clamp(floor(Int,maximum(p[2] for p in poly)/dz),0,mesh.n_axial-1)
    rule=((1/3,1/3,1/3,.225),
        (.059715871789770,.470142064105115,.470142064105115,.132394152788506),
        (.470142064105115,.059715871789770,.470142064105115,.132394152788506),
        (.470142064105115,.470142064105115,.059715871789770,.132394152788506),
        (.797426985353087,.101286507323456,.101286507323456,.125939180544827),
        (.101286507323456,.797426985353087,.101286507323456,.125939180544827),
        (.101286507323456,.101286507323456,.797426985353087,.125939180544827))
    ca=cosd(mesh.p.semicone_angle_deg)
    for i in imin:imax,j in jmin:jmax
        q=poly
        for (axis,bound,above) in ((1,mesh.origin+i*dt,true),(1,mesh.origin+(i+1)*dt,false),
                (2,j*dz,true),(2,(j+1)*dz,false))
            q=journal_clip(q,axis,bound,above)
        end
        for k in 2:length(q)-1
            a,b,c=q[1],q[k],q[k+1]
            area=abs((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2
            area>0 || continue
            for (wa,wb,wc,w) in rule
                theta,z=wa*a+wb*b+wc*c
                radius=journal_radius(mesh.p,z*1e3)*1e-3
                visit(theta,z,w*area*radius/ca)
            end
        end
    end
end

function journal_pressure_operator(mesh,grid)
    rows,cols,vals=Int[],Int[],Float64[]
    areas=zeros(length(grid.quads));reservoir=zeros(3size(mesh.nodes,1));reservoir_area=0.
    R=grid.parameters.housing_midpoint_radius*1e-3;ca=cosd(mesh.p.semicone_angle_deg)
    function integrate(poly,column)
        journal_surface_quadrature(mesh,poly,(theta,z,weight)->begin
            ids,N=journal_interface_weights(mesh,theta,z)
            traction=-[ca*cos(theta),ca*sin(theta),sind(mesh.p.semicone_angle_deg)]
            if column>0
                areas[column]+=weight
            else
                reservoir_area+=weight
            end
            for i in 1:4,k in 1:3
                value=N[i]*traction[k]*weight;row=3*(ids[i]-1)+k
                if column==0
                    reservoir[row]+=value
                else
                    push!(rows,row);push!(cols,column);push!(vals,value)
                end
            end
        end)
    end
    for cell in eachindex(grid.quads)
        Q=housing_interface_chart(grid,cell)
        integrate([[Q[i,1]/R,Q[i,2]*ca] for i in 1:4],cell)
    end
    if hasproperty(grid,:feed_edges)
        s=grid.surface;rim=unique(vcat(first.(grid.feed_edges),last.(grid.feed_edges)))
        points=[[s.nodes[i][1]/R,s.nodes[i][2]*ca] for i in rim]
        origin=atan(sum(sin(p[1]) for p in points),sum(cos(p[1]) for p in points))
        points=[[origin+mod(p[1]-origin+pi,2pi)-pi,p[2]] for p in points]
        centre=reduce(+,points)/length(points)
        sort!(points;by=p->atan(p[2]-centre[2],R*(p[1]-centre[1])))
        for i in eachindex(points)
            integrate([centre,points[i],points[mod1(i+1,length(points))]],0)
        end
    end
    return (;matrix=sparse(rows,cols,vals,3size(mesh.nodes,1),length(grid.quads)),
        areas_m2=areas,reservoir_load_per_pa=reservoir,reservoir_area_m2=reservoir_area)
end
function assemble_journal(grid;youngs_modulus_pa,poisson_ratio,n_theta,n_axial,n_radial)
    mesh=journal_solid_mesh(grid.parameters;n_theta,n_axial,n_radial)
    a=assemble_elasticity(mesh.nodes,mesh.cells;youngs_modulus_pa,poisson_ratio)
    constrained=constrain_elasticity(a,elastic_dofs(mesh.end_nodes))
    pressure=journal_pressure_operator(mesh,grid)
    return (;constrained...,mesh,pressure_operator=pressure.matrix,
        interface_areas_m2=pressure.areas_m2,reservoir_load_per_pa=pressure.reservoir_load_per_pa,
        reservoir_area_m2=pressure.reservoir_area_m2)
end
function journal_gap_operator(a,points)
    rows,cols,vals=Int[],Int[],Float64[]
    ca=cosd(a.mesh.p.semicone_angle_deg)
    for (row,(theta,z)) in enumerate(points)
        ids,N=journal_interface_weights(a.mesh,theta,z)
        normal=[ca*cos(theta),ca*sin(theta),sind(a.mesh.p.semicone_angle_deg)]
        for i in 1:4,k in 1:3
            push!(rows,row);push!(cols,3*(ids[i]-1)+k);push!(vals,N[i]*normal[k])
        end
    end
    return sparse(rows,cols,vals,length(points),3size(a.nodes,1))
end
