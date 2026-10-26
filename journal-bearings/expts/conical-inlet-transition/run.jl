using TOML,SHA,Printf
include(joinpath(@__DIR__,"../conical-fitted-inlet/run.jl"))
include(joinpath(@__DIR__,"layout.jl"))

function export_patch(path,m)
    d=macro_data(m)
    open(path,"w") do io
        println(io,"cell\tx1_mm\tz1_mm\tx2_mm\tz2_mm\tx3_mm\tz3_mm\tx4_mm\tz4_mm\tcorner")
        for (k,xy) in enumerate(m.localnodes)
            sin(sum(xy[:,1])/4/m.R)>0 || continue
            p=[begin
                x,y=xy[i,:];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                (1e3*rh*cos(x/m.R),1e3*(z-.05))
            end for i in 1:4]
            all(v->abs(v[1])<24 && abs(v[2])<24,p) || continue
            println(io,join((k,Iterators.flatten(p)...,m.corner_scores[k]),'\t'))
        end
    end
    open(replace(path,".tsv"=>"-junctions.tsv"),"w") do io
        println(io,"x_mm\tz_mm\tvalence")
        for i in findall(==(5),d.valence)
            x,y=m.nodes[i];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
            println(io,join((1e3*rh*cos(x/m.R),1e3*(z-.05),5),'\t'))
        end
    end
end

function transition_run(out)
    ispath(out) && error("refusing to overwrite transition evidence: $out")
    mkpath(out)
    root=FITTED_ROOT
    prior=joinpath(root,"out/cases/conical-floating-bush/runs/fitted-inlet-diagnostic-reviewed")
    oldrecord=TOML.parsefile(joinpath(prior,"comparison.toml"))
    for (file,expected) in oldrecord["source_sha256"]
        source=file=="expts/conical-fitted-inlet/mesh.jl" ?
            joinpath(root,"out/reviews/fitted-inlet-before-transition/mesh.jl") :
            file=="mesh/2d/conical/grid.jl" ?
            joinpath(root,"out/reviews/pre-surface-inlet-promotion/mesh/2d/conical/grid.jl") : joinpath(root,file)
        bytes2hex(sha256(read(source)))==expected || error("historical source mismatch: $file")
    end
    oldlines=split.(readlines(joinpath(prior,"comparison.tsv")),'\t');oldheader=first(oldlines)
    base=fitted_mesh();placement=@timed controlled_layout(base);result=placement.value
    open(joinpath(out,"placement.tsv"),"w") do io
        println(io,"iteration\tenergy\tgradient_norm")
        foreach(h->println(io,join(values(h),'\t')),result.history)
    end
    open(joinpath(out,"node-movement.tsv"),"w") do io
        println(io,"id\tx_old_mm\tz_old_mm\tx_new_mm\tz_new_mm")
        for i in macro_data(base).free
            values=Float64[]
            for m in (base,result.mesh)
                x,y=m.nodes[i];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                append!(values,[1e3*rh*cos(x/m.R),1e3*(z-.05)])
            end
            println(io,join((i,values...),'\t'))
        end
    end
    reference=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    displacement=(-9.152501266813789e-7,1.7844427760798323e-7,0.)
    pose=floating_bush_pose(reference,displacement)
    rows=NamedTuple[]
    open(joinpath(out,"comparison.tsv"),"w") do io
        for level in (1,2,4),kind in ("original","frozen","controlled")
            m=kind=="original" ? fitted_mesh(level) : subdivide_layout(kind=="frozen" ? base : result.mesh,level)
            if kind=="original"
                old=only(filter(r->r[1]=="fitted" && parse(Int,r[2])==level,oldlines[2:end]))
                value(key)=parse(Float64,old[findfirst(==(key),oldheader)])
                s=(;force=[value("Fx_N"),value("Fy_N"),value("Fz_N")],net_feed=value("net_feed_m3_s"),
                    injection=value("gross_rim_out_m3_s"),return_flow=value("gross_rim_in_m3_s"),
                    outward=[value("large_end_m3_s"),value("small_end_m3_s")],peak_pressure=1e6*value("land_peak_MPa"),
                    minimum_fill=value("minimum_fill"),iterations=Int(value("iterations")),imbalance=value("imbalance"))
                isapprox(m.quality.minimum_scaled_corner,value("minimum_corner");atol=1e-13) || error("historical geometry changed")
                wall=NaN
            else
                elapsed=@timed solve_fitted(m,pose);s=elapsed.value;wall=elapsed.time
                write_fitted_vtk(joinpath(out,"$kind-$level.vtk"),m,s)
            end
            q=m.quality
            row=(;kind,level,cells=length(m.quads),corner=q.minimum_scaled_corner,alignment=q.minimum_face_alignment,
                nonorthogonality_deg=q.maximum_nonorthogonality_deg,area_ratio=q.maximum_neighbour_area_ratio,
                edge_ratio=q.maximum_edge_ratio,Fx_N=s.force[1],Fy_N=s.force[2],Fz_N=s.force[3],
                radial_residual_N=hypot(s.force[1],s.force[2]+1000),axial_residual_N=s.force[3]-200,
                net_feed_L_min=s.net_feed*60000,rim_out_L_min=s.injection*60000,rim_in_L_min=s.return_flow*60000,
                large_end_L_min=s.outward[1]*60000,small_end_L_min=s.outward[2]*60000,
                land_peak_MPa=s.peak_pressure/1e6,minimum_fill=s.minimum_fill,iterations=s.iterations,
                flow_imbalance=s.imbalance,solve_wall_s=wall)
            isempty(rows) && println(io,join(keys(row),'\t'));push!(rows,row)
            println(io,join(values(row),'\t'));flush(io)
            export_patch(joinpath(out,"patch-$kind-$level.tsv"),m)
            @printf("%s level=%d cells=%d corner=%.6f alignment=%.6f area ratio=%.4f Fz=%.5f N flow=%.8f L/min\n",kind,level,row.cells,row.corner,row.alignment,row.area_ratio,row.Fz_N,row.net_feed_L_min)
            flush(stdout)
        end
    end
    files=["expts/conical-inlet-transition/"*s for s in ("layout.jl","run.jl","checks.jl","plot.py")]
    append!(files,["expts/conical-fitted-inlet/"*s for s in ("mesh.jl","operator.jl","checks.jl","run.jl")])
    append!(files,["expts/surface-inlet-hex/mesh.jl","bearing.jl","mesh/2d/conical/grid.jl"])
    append!(files,["simulation/reynolds/"*s for s in ("conical_equilibrium.jl","conical_jfo.jl","conical.jl","jfo_active_set.jl","cylindrical.jl")])
    metadata=Dict("status"=>"Surface transition experiment; no equilibrium re-solve, 3D extrusion or production promotion",
        "method"=>"Same coarse connectivity; coupled cotangent-angle and log-area placement, then fixed-map subdivision with exact-rim Coons correction. Not a new transition-template topology or a reproduction of a full published elliptic mesher.",
        "area_weight"=>.25,"area_anchor_weight"=>.01,"iteration_cap"=>180,
        "placement_wall_s"=>placement.time,"initial_energy"=>result.initial_energy,"final_energy"=>result.final_energy,
        "movable_nodes"=>result.free_nodes,"maximum_node_movement_m"=>result.maximum_move_m,
        "movable_halo_halfwidth_m"=>.020,"displacement_m"=>collect(displacement),
        "eccentricity_ratio"=>pose.parameters.eccentricity_ratio,"minimum_normal_gap_m"=>pose.minimum_normal_gap_m,
        "radial_clearance_m"=>pose.radial_clearance_m,"feed_diameter_m"=>.005,"feed_pressure_pa_gauge"=>5e5,
        "end_pressure_pa_gauge"=>[0.,0.],"cavitation_pressure_pa_gauge"=>0.,"rpm"=>2000.,"viscosity_pa_s"=>.0277,
        "flow_note"=>"Gross rim exchange is NOT gross pump-line injection/return; net supply is comparable.",
        "julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
        "historical_comparison_sha256"=>bytes2hex(sha256(read(joinpath(prior,"comparison.tsv")))),
        "historical_manifest_sha256"=>bytes2hex(sha256(read(joinpath(prior,"comparison.toml")))),
        "source_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(root,f)))) for f in files))
    open(joinpath(out,"comparison.toml"),"w") do io;TOML.print(io,metadata;sorted=true);end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide a fresh output directory")
    transition_run(ARGS[1])
end
