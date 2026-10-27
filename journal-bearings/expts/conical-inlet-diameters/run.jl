# One fixed-budget diameter-relative recipe. No per-diameter parameter search.
using Test,TOML,SHA,Printf
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../../mesh/vtk.jl"))
module FieldReference
include(joinpath(@__DIR__,"../conical-fitted-inlet/checks.jl"))
end
const C=ControlledSurfaceInlet
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const DIAMETERS=4:10
quality_pass(m)=m.quality.minimum_scaled_corner>.9 &&
    m.quality.minimum_face_alignment>.9 && m.quality.maximum_neighbour_area_ratio<2

function geometry_export(m)
    nodes=[begin
        x,y=point;z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
        (1e3*rh*cos(x/m.R),1e3*rh*sin(x/m.R),1e3*z)
    end for point in m.nodes]
    areas=m.areas.*[(m.R+(.05-m.centre[k,2]*m.ca)*tand(10.))/m.R*1e6 for k in eachindex(m.quads)]
    # Geometry-only centred reference: no physical operating-point solve.
    return (;nodes,quads=m.quads,areas,radial_gaps=fill(.05,length(m.quads)),
        feed_edges=m.edge_nodes[findall(==(:feed),m.boundary)],surface=m)
end

function save_coarse(out,d,m)
    open(joinpath(out,"coarse-$d-nodes.tsv"),"w") do io
        foreach(x->println(io,join(x,'\t')),m.nodes)
    end
    open(joinpath(out,"coarse-$d-quads.tsv"),"w") do io
        foreach(q->println(io,join(q,'\t')),m.quads)
    end
end

function load_coarse(out,d)
    nodes=[Tuple(parse.(Float64,split(line,'\t'))) for line in readlines(joinpath(out,"coarse-$d-nodes.tsv"))]
    quads=[Tuple(parse.(Int,split(line,'\t'))) for line in readlines(joinpath(out,"coarse-$d-quads.tsv"))]
    edges=Dict{Tuple{Int,Int},Int}()
    for q in quads,i in 1:4
        e=Tuple(sort([q[i],q[mod1(i+1,4)]]));edges[e]=get(edges,e,0)+1
    end
    ca=cosd(10.)
    rim=Set(e for (e,n) in edges if n==1 && all(v->1e-10<nodes[v][2]<.1/ca-1e-10,e))
    return C.assemble_surface_mesh(nodes,quads,rim;inlet_radius_m=d*.0005)
end

function patch_export(path,m)
    g=geometry_export(m)
    face_score=ones(length(m.quads))
    for f in eachindex(m.owner)
        k,l=m.owner[f],m.neighbour[f]
        face_score[k]=min(face_score[k],m.face_alignment[f])
        l>0 && (face_score[l]=min(face_score[l],m.face_alignment[f]))
    end
    open(path,"w") do io
        println(io,"cell\tx1_mm\tz1_mm\tx2_mm\tz2_mm\tx3_mm\tz3_mm\tx4_mm\tz4_mm\tcorner\talignment")
        for (k,q) in enumerate(m.quads)
            sum(g.nodes[v][2] for v in q)>0 || continue
            p=[(g.nodes[v][1],g.nodes[v][3]-50) for v in q]
            println(io,join((k,Iterators.flatten(p)...,m.corner_scores[k],face_score[k]),'\t'))
        end
    end
    f=argmin(m.face_alignment);x,y=m.mid[f];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
    return (1e3*rh*cos(x/m.R),1e3*(z-.05))
end

function true_mouth_area(d)
    r=d*.0005;n=20000;dz=2r/n
    sum(1:n) do j
        z=.05-r+(j-.5)*dz;rh=.05005+(.05-z)*tand(10.)
        rh*2asin(sqrt(r^2-(z-.05)^2)/rh)*dz/cosd(10.)
    end
end

function diameter_geometry(out)
    ispath(out) && error("refusing to overwrite diameter evidence: $out")
    mkpath(out)
    firstrow=true
    open(joinpath(out,"quality.tsv"),"w") do io
        @testset "Seven diameters: topology, exact rim and independent area/traction" begin
            for d in DIAMETERS
                timed=@timed C.controlled_layout(C.reference_seed(;inlet_diameter_mm=d))
                p=timed.value;coarse=p.mesh;save_coarse(out,d,coarse)
                @test load_coarse(out,d).nodes==coarse.nodes
                @test load_coarse(out,d).quads==coarse.quads
                @test p.final_energy<p.initial_energy
                open(joinpath(out,"placement-$d.tsv"),"w") do f
                    println(f,"iteration\tenergy\tgradient_norm")
                    foreach(h->println(f,join(values(h),'\t')),p.history)
                end
                exact_area=true_mouth_area(d);area_errors=Float64[]
                for level in (1,2,4)
                    mt=@timed C.subdivide_layout(coarse,level);m=mt.value;g=geometry_export(m)
                    q=m.quality;info=C.macro_data(m)
                    @test length(m.quads)==8416level^2
                    @test q.rim_edges==32level
                    @test count(==(:large),m.boundary)==128level
                    @test count(==(:small),m.boundary)==128level
                    @test length(m.nodes)-length(m.owner)+length(m.quads)==-1
                    @test all(v->v in (2,4,5),info.valence)
                    @test count(==(5),info.valence)==4
                    @test all(v->info.valence[v]==2,info.boundary)
                    @test all(g.feed_edges) do edge
                        all(edge) do v
                            x,y,z=g.nodes[v]
                            y>0 && abs(hypot(x,z-50)-d/2)<1e-10
                        end
                    end
                    @test minimum(m.corner_scores)>0 && minimum(m.face_alignment)>0
                    uniform=FieldReference.fitted_force(m,fill(5e5,length(m.quads)))
                    area=2pi*m.R*m.L/m.ca-uniform.area
                    ae=abs(area/exact_area-1);push!(area_errors,ae)
                    @test ae<.008/level^2
                    @test abs(uniform.force[1])<.002
                    @test abs(uniform.force[2]/(-5e5*pi*(d*.0005)^2)-1)<.008/level^2
                    @test abs(uniform.force[3]-5e5*sind(10.)*(2pi*m.R*m.L/m.ca-exact_area))<5e5*sind(10.)*exact_area*.008/level^2+.0001
                    write_vtk(joinpath(out,"diameter-$d-level-$level.vtk"),g)
                    wx,wz=patch_export(joinpath(out,"patch-$d-$level.tsv"),m)
                    row=(;diameter_mm=d,level,cells=length(m.quads),rim_edges=q.rim_edges,
                        minimum_corner=q.minimum_scaled_corner,minimum_alignment=q.minimum_face_alignment,
                        nonorthogonality_deg=q.maximum_nonorthogonality_deg,max_area_ratio=q.maximum_neighbour_area_ratio,
                        max_edge_ratio=q.maximum_edge_ratio,quality_pass=quality_pass(m),
                        worst_face_x_mm=wx,worst_face_z_from_mid_mm=wz,relative_mouth_area_error=ae,
                        placement_wall_s=timed.time,refinement_wall_s=mt.time,free_nodes=p.free_nodes,
                        placement_iterations=length(p.history),initial_energy=p.initial_energy,final_energy=p.final_energy,
                        maximum_node_move_mm=p.maximum_move_m*1e3)
                    firstrow && println(io,join(keys(row),'\t'));firstrow=false
                    println(io,join(values(row),'\t'));flush(io)
                    @printf("D=%d mm level=%d corner=%.5f alignment=%.5f area ratio=%.4f accepted=%s\n",d,level,row.minimum_corner,row.minimum_alignment,row.max_area_ratio,string(row.quality_pass));flush(stdout)
                end
                @test area_errors[3]<area_errors[1]/10
            end
        end
    end
    sources=["bearing.jl","Project.toml","mesh/2d/grid.jl","mesh/2d/surface_feed.jl"]
    append!(sources,["mesh/2d/"*s for s in ("quality.jl","metrics.jl","inlet/controlled.jl","inlet/placement.jl","inlet/seed.jl")])
    append!(sources,["expts/conical-inlet-diameters/"*s for s in ("run.jl","operator_checks.jl","plot.py")])
    append!(sources,["expts/conical-fitted-inlet/"*s for s in ("mesh.jl","operator.jl","checks.jl")])
    append!(sources,["simulation/reynolds/"*s for s in ("conical_equilibrium.jl","conical_jfo.jl","conical.jl","jfo_active_set.jl","cylindrical.jl")])
    push!(sources,"expts/surface-inlet-hex/mesh.jl")
    meta=Dict("diameters_mm"=>collect(DIAMETERS),"levels"=>[1,2,4],
        "status"=>"Geometry and manufactured-field study; no physical diameter/load tuning or equilibrium",
        "recipe"=>"q=8, outer layers=9, background128x64; patch half-width1.5D, halo4D, max placement step0.04D; inverse-Hessian and gradient scales follow D/5",
        "area_weight"=>.25,"area_anchor_weight"=>.01,"placement_iteration_cap"=>180,
        "minimum_corner_bound"=>.9,"minimum_alignment_bound"=>.9,"maximum_area_ratio_bound"=>2.,
        "geometry_export_eccentricity_ratio"=>0.,"radial_clearance_m"=>50e-6,"minimum_normal_gap_m"=>50e-6*cosd(10.),
        "julia_version"=>string(VERSION),"blas_threads"=>FieldReference.BLAS.get_num_threads(),
        "julia_threads"=>Threads.nthreads(),"source_sha256"=>Dict(s=>bytes2hex(sha256(read(joinpath(ROOT,s)))) for s in sources),
        "coarse_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(out,f)))) for f in readdir(out) if startswith(f,"coarse-")))
    open(joinpath(out,"geometry.toml"),"w") do io;TOML.print(io,meta;sorted=true);end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide a fresh evidence directory")
    diameter_geometry(ARGS[1])
end
