# A fixed macro-layout experiment, not an implementation of an entire published
# mesher. Existing surface connectivity is retained; angles and cell-size
# transitions are optimized together ONCE, before geometry-preserving refinement.
# Include expts/conical-fitted-inlet/operator.jl or run.jl first.

function macro_data(m)
    nodes=hcat([collect(x) for x in m.nodes]...)'
    nodes=Matrix(nodes)
    boundary_edges=Dict{Tuple{Int,Int},Int}()
    valence=zeros(Int,length(m.nodes))
    for q in m.quads
        for i in 1:4
            valence[q[i]]+=1
            e=Tuple(sort([q[i],q[mod1(i+1,4)]]));boundary_edges[e]=get(boundary_edges,e,0)+1
        end
    end
    boundary=Set(v for (e,n) in boundary_edges if n==1 for v in e)
    rim=Set(e for (e,n) in boundary_edges if n==1 &&
        all(v->1e-10<nodes[v,2]<m.L/m.ca-1e-10,e))
    # Keep the exact mouth, axial ends, seam and remote background unchanged.
    # The movable halo includes the ring/background support interface.
    free=[i for i in eachindex(valence) if valence[i]>0 && !(i in boundary) &&
        abs(nodes[i,1]-m.R*pi/2)<.020 && abs(nodes[i,2]-.05/m.ca)<.020/m.ca]
    lookup=Dict(v=>k for (k,v) in enumerate(free))
    active=findall(q->any(v->haskey(lookup,v),q),m.quads)
    return (;nodes,free,active,rim,valence,boundary)
end

# Analytic gradient: angle energy cot(angle)^2, plus neighbouring log-area
# differences. Weights are explicit experiment choices, not literature constants.
function layout_energy(m,d,variables;area_weight=.25,anchor_weight=.01)
    X=copy(d.nodes);X[d.free,:]=reshape(variables,:,2)
    gradient=zeros(size(X));areas=zeros(length(m.quads));area_gradient=zeros(length(m.quads))
    energy=0.
    for (k,q) in enumerate(m.quads)
        # Cells at the periodic seam are fixed and excluded from angle updates.
        xy=Matrix(X[collect(q),:])
        maximum(xy[:,1])-minimum(xy[:,1])>m.period/2 && (xy[xy[:,1].<m.period/2,1].+=m.period)
        # Translate coordinates before shoelace area to avoid cancellation.
        xy .-= xy[1:1,:]
        A=sum(xy[i,1]*xy[mod1(i+1,4),2]-xy[i,2]*xy[mod1(i+1,4),1] for i in 1:4)/2
        A>0 || return (Inf,zeros(length(variables)))
        areas[k]=A
        k in d.active || continue
        for i in 1:4
            j,l=mod1(i+1,4),mod1(i-1,4)
            a=xy[j,:]-xy[i,:];b=xy[l,:]-xy[i,:]
            cross=a[1]*b[2]-a[2]*b[1];dotab=dot(a,b)
            cross>0 || return (Inf,zeros(length(variables)))
            energy+=(dotab/cross)^2
            ga=2dotab/cross^2*b-2dotab^2/cross^3*[b[2],-b[1]]
            gb=2dotab/cross^2*a-2dotab^2/cross^3*[-a[2],a[1]]
            gradient[q[j],:].+=ga;gradient[q[l],:].+=gb;gradient[q[i],:].-=ga+gb
        end
    end
    for f in eachindex(m.owner)
        k,l=m.owner[f],m.neighbour[f]
        l==0 && continue
        (k in d.active || l in d.active) || continue
        t=log(areas[k]/areas[l]);energy+=area_weight*t^2
        area_gradient[k]+=2area_weight*t/areas[k];area_gradient[l]-=2area_weight*t/areas[l]
    end
    for k in d.active
        # A small area anchor discourages collapse/drift without freezing sizes.
        t=log(areas[k]/m.areas[k]);energy+=anchor_weight*t^2
        area_gradient[k]+=2anchor_weight*t/areas[k]
    end
    for (k,q) in enumerate(m.quads)
        area_gradient[k]==0 && continue
        for i in 1:4
            j,l=mod1(i+1,4),mod1(i-1,4)
            dx=X[q[l],1]-X[q[j],1];dx=mod(dx+m.period/2,m.period)-m.period/2
            gradient[q[i],1]+=area_gradient[k]*(X[q[j],2]-X[q[l],2])/2
            gradient[q[i],2]+=area_gradient[k]*dx/2
        end
    end
    return energy,vec(gradient[d.free,:])
end

function controlled_layout(m;max_iterations=180)
    d=macro_data(m);x=vec(copy(d.nodes[d.free,:]));f,g=layout_energy(m,d,x)
    ss=Vector{Float64}[];ys=Vector{Float64}[];history=NamedTuple[]
    for iteration in 1:max_iterations
        push!(history,(;iteration,energy=f,gradient_norm=norm(g)))
        norm(g)<1e-5 && break
        q=copy(g);alpha=zeros(length(ss))
        for i in reverse(eachindex(ss))
            alpha[i]=dot(ss[i],q)/dot(ys[i],ss[i]);q.-=alpha[i].*ys[i]
        end
        gamma=isempty(ss) ? 1e-7 : dot(ss[end],ys[end])/dot(ys[end],ys[end])
        direction=gamma.*q
        for i in eachindex(ss)
            beta=dot(ys[i],direction)/dot(ys[i],ss[i]);direction.+=ss[i].*(alpha[i]-beta)
        end
        direction .*= -1
        dot(direction,g)<0 || (direction=-g*1e-7)
        maxmove=maximum(sqrt.(sum(reshape(direction,:,2).^2,dims=2)))
        maxmove>.0002 && (direction.*=.0002/maxmove) # bounded 0.2 mm step
        accepted=false
        for backtrack in 0:22
            trial=x+2.0^(-backtrack)*direction;ft,gt=layout_energy(m,d,trial)
            if isfinite(ft) && ft<=f+1e-4*dot(g,trial-x)
                s=trial-x;y=gt-g
                if dot(s,y)>1e-12*norm(s)*norm(y)
                    push!(ss,s);push!(ys,y)
                    length(ss)>8 && (popfirst!(ss);popfirst!(ys))
                end
                x,f,g=trial,ft,gt;accepted=true;break
            end
        end
        accepted || break
    end
    nodes=copy(d.nodes);nodes[d.free,:]=reshape(x,:,2)
    mesh=assemble_surface_mesh([Tuple(nodes[i,:]) for i in axes(nodes,1)],m.quads,d.rim)
    return (;mesh,history,initial_energy=first(history).energy,final_energy=f,
        free_nodes=length(d.free),maximum_move_m=maximum(sqrt.(sum((nodes-d.nodes).^2,dims=2))))
end

# Refine the SAME macro map, not a newly solved harmonic layout. A Coons-type
# correction follows the exact drilling rim; straight macro edges are shared.
function subdivide_layout(m,level)
    level in (1,2,4) || throw(ArgumentError("comparison levels are 1,2,4"))
    d=macro_data(m);nodes=NTuple{2,Float64}[];quads=NTuple{4,Int}[]
    ids=Dict{Any,Int}();rim=Set{Tuple{Int,Int}}()
    function add(key,xy)
        if haskey(ids,key)
            old=nodes[ids[key]];dx=mod(xy[1]-old[1]+m.period/2,m.period)-m.period/2
            hypot(dx,xy[2]-old[2])<1e-12 || error("nonmatching shared macro edge")
            return ids[key]
        end
        push!(nodes,(mod(xy[1],m.period),xy[2]));ids[key]=length(nodes)
        return length(nodes)
    end
    for (k,q) in enumerate(m.quads)
        xy=m.localnodes[k];boundaryside=0
        for side in 1:4
            Tuple(sort([q[side],q[mod1(side+1,4)]])) in d.rim && (boundaryside=side)
        end
        angle=Float64[]
        if boundaryside>0
            for i in (boundaryside,mod1(boundaryside+1,4))
                x,y=xy[i,:];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                push!(angle,atan(z-.05,rh*cos(x/m.R)))
            end
        end
        grid=zeros(Int,level+1,level+1)
        for j in 0:level,i in 0:level
            u,v=i/level,j/level
            position=vec(xy'*[(1-u)*(1-v),u*(1-v),u*v,(1-u)*v])
            if boundaryside>0
                t,w=boundaryside==1 ? (u,1-v) : boundaryside==2 ? (v,u) : boundaryside==3 ? (1-u,v) : (1-v,1-u)
                delta=mod(angle[2]-angle[1]+pi,2pi)-pi;theta=angle[1]+t*delta
                X=.0025*cos(theta);z=.05+.0025*sin(theta);rh=m.R+(.05-z)*tand(10.)
                exact=[m.R*acos(X/rh),z/m.ca]
                chord=(1-t)*xy[boundaryside,:]+t*xy[mod1(boundaryside+1,4),:]
                position.+=w.*(exact-chord)
            end
            if i in (0,level) && j in (0,level)
                vertex=j==0 ? (i==0 ? q[1] : q[2]) : (i==0 ? q[4] : q[3])
                key=(:vertex,vertex)
            elseif j==0 || i==level || j==level || i==0
                side,t=j==0 ? (1,i) : i==level ? (2,j) : j==level ? (3,level-i) : (4,level-j)
                a,b=q[side],q[mod1(side+1,4)]
                key=(:edge,min(a,b),max(a,b),a<b ? t : level-t)
            else
                key=(:cell,k,i,j)
            end
            grid[i+1,j+1]=add(key,position)
        end
        for j in 1:level,i in 1:level
            child=(grid[i,j],grid[i+1,j],grid[i+1,j+1],grid[i,j+1]);push!(quads,child)
            boundaryside==1 && j==1 && push!(rim,Tuple(sort([child[1],child[2]])))
            boundaryside==2 && i==level && push!(rim,Tuple(sort([child[2],child[3]])))
            boundaryside==3 && j==level && push!(rim,Tuple(sort([child[3],child[4]])))
            boundaryside==4 && i==1 && push!(rim,Tuple(sort([child[4],child[1]])))
        end
    end
    return assemble_surface_mesh(nodes,quads,rim;level)
end
