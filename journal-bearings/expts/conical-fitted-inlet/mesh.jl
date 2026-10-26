# One case-specific surface candidate. Reuse the old layout, not its Hex8 sweep.
module RetainedInletSurface
include(joinpath(@__DIR__,"../surface-inlet-hex/mesh.jl"))
end
using LinearAlgebra, SparseArrays

function fitted_mesh(level=1;regular=false)
    level in (1,2,4) || throw(ArgumentError("supported comparison levels: 1,2,4"))
    nt,nz=128level,64level
    R,L,ca=.05005,.1,cosd(10.)
    period=2pi*R
    if regular
        nodes=[(R*(i-1)*2pi/nt,(j-1)*L/nz/ca) for j in 1:nz+1 for i in 1:nt]
        id(i,j)=mod(i-1,nt)+1+nt*(j-1)
        quads=[(id(i,j),id(i+1,j),id(i+1,j+1),id(i,j+1)) for j in 1:nz for i in 1:nt]
        rim=Set{Tuple{Int,Int}}()
    else
        p=RetainedInletSurface.BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
            bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
        uv,allquads,feed=RetainedInletSurface.bearing_surface(p;ntheta=nt,nz,
            inlet_radius=2.5,inlet_z=50.,q=8level,inner_layers=1,
            outer_layers=9level,patch_half_width=7.5,curved_transition=true,core_fraction=.9)
        nodes=[(R*t,z*1e-3/ca) for (t,z) in uv]
        quads=allquads[.!feed]
        # Prescribed-pressure reservoir interior has no unknowns. Its rim is
        # retained; this is NOT cutting a pipe/void through a 3D film volume.
        landedges=Set(Tuple(sort([q[i],q[mod1(i+1,4)]])) for q in quads for i in 1:4)
        rim=Set(Tuple(sort([q[i],q[mod1(i+1,4)]])) for q in allquads[feed] for i in 1:4)
        intersect!(rim,landedges)
    end
    return assemble_surface_mesh(nodes,quads,rim;R,L,ca,period,regular,level)
end

# Same geometry/quality calculation for a supplied conforming surface layout.
function assemble_surface_mesh(nodes,quads,rim;R=.05005,L=.1,ca=cosd(10.),
        period=2pi*R,regular=false,level=1)
    centre=zeros(length(quads),2);areas=zeros(length(quads));corners=Float64[];ratios=Float64[]
    localnodes=Matrix{Float64}[]
    for (k,q) in enumerate(quads)
        xy=hcat([collect(nodes[i]) for i in q]...)'
        xy=Matrix(xy)
        maximum(xy[:,1])-minimum(xy[:,1])>period/2 && (xy[xy[:,1].<period/2,1].+=period)
        push!(localnodes,xy)
        cross2(a,b)=a[1]*b[2]-a[2]*b[1]
        cross=[cross2(xy[i,:],xy[mod1(i+1,4),:]) for i in 1:4]
        area=sum(cross)/2
        area>0 || error("nonpositive cell $k")
        centre[k,:]=sum((xy[i,:]+xy[mod1(i+1,4),:])*cross[i] for i in 1:4)/(6area)
        centre[k,1]=mod(centre[k,1],period);areas[k]=area
        lengths=Float64[]
        for i in 1:4
            z=xy[i,2]*ca;rh=R+(.05-z)*tand(10.)
            S=[rh/R,1.]
            a=(xy[mod1(i+1,4),:]-xy[i,:]).*S
            b=(xy[mod1(i-1,4),:]-xy[i,:]).*S
            push!(corners,cross2(a,b)/(norm(a)*norm(b)));push!(lengths,norm(a))
        end
        push!(ratios,maximum(lengths)/minimum(lengths))
    end
    owner=Int[];neighbour=Int[];mid=Vector{Float64}[];normals=Vector{Float64}[];lengths=Float64[];boundary=Symbol[]
    edgeids=Tuple{Int,Int}[];lookup=Dict{Tuple{Int,Int},Int}()
    for (k,q) in enumerate(quads),i in 1:4
        pair=Tuple(sort([q[i],q[mod1(i+1,4)]]))
        if haskey(lookup,pair)
            f=lookup[pair];neighbour[f]==0 || error("nonmanifold edge");neighbour[f]=k
        else
            xy=localnodes[k];a=xy[i,:];b=xy[mod1(i+1,4),:];t=b-a;l=norm(t)
            m=(a+b)/2;m[1]=mod(m[1],period)
            push!(owner,k);push!(neighbour,0);push!(mid,m);push!(normals,[t[2],-t[1]]/l);push!(lengths,l)
            push!(edgeids,pair);lookup[pair]=length(owner);push!(boundary,:internal)
        end
    end
    wrapdx(x)=mod(x+period/2,period)-period/2
    cosines=Float64[];growth=Float64[]
    for f in eachindex(owner)
        k,l=owner[f],neighbour[f];x=mid[f]
        if l==0
            boundary[f]=edgeids[f] in rim ? :feed : abs(x[2])<1e-12 ? :large : abs(x[2]-L/ca)<1e-12 ? :small : error("unclassified boundary")
        end
        d=l==0 ? x-centre[k,:] : centre[l,:]-centre[k,:];d[1]=wrapdx(d[1])
        rh=R+(.05-x[2]*ca)*tand(10.);S=[rh/R,1.]
        normal=normals[f]./S
        push!(cosines,dot(d.*S,normal)/(norm(d.*S)*norm(normal)))
        l>0 && push!(growth,max(areas[k],areas[l])/min(areas[k],areas[l]))
    end
    minimum(corners)>0 && minimum(cosines)>0 || error("invalid corner or face orientation")
    quality=(;minimum_scaled_corner=minimum(corners),maximum_edge_ratio=maximum(ratios),
        minimum_face_alignment=minimum(cosines),maximum_nonorthogonality_deg=acosd(minimum(cosines)),
        maximum_neighbour_area_ratio=maximum(growth),cells=length(quads),rim_edges=count(==(:feed),boundary))
    return (;nodes,quads,centre,areas,owner,neighbour,mid,normals,lengths,boundary,localnodes,R,L,ca,period,quality,regular,level,
        corner_scores=vec(minimum(reshape(corners,4,:),dims=1)),face_alignment=cosines,neighbour_area_ratios=growth)
end
