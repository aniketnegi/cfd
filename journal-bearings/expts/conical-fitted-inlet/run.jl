using TOML,SHA,Printf
include(joinpath(@__DIR__,"operator.jl"))
const FITTED_ROOT=normpath(joinpath(@__DIR__,"../.."))

function write_fitted_vtk(path,m,s)
    ispath(path) && error("refusing to overwrite $path")
    open(path,"w") do io
        println(io,"# vtk DataFile Version 3.0\nFitted surface diagnostic, not a 3D mesh\nASCII\nDATASET UNSTRUCTURED_GRID")
        println(io,"POINTS $(length(m.nodes)) double")
        for (x,y) in m.nodes
            theta=x/m.R;z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
            println(io,"$(1e3*rh*cos(theta)) $(1e3*rh*sin(theta)) $(1e3*z)")
        end
        println(io,"CELLS $(length(m.quads)) $(5length(m.quads))")
        for q in m.quads;println(io,"4 ",join(q.-1,' '));end
        println(io,"CELL_TYPES $(length(m.quads))");for q in m.quads;println(io,9);end
        println(io,"CELL_DATA $(length(m.quads))")
        for (name,field) in (("pressure_pa",s.pressure),("fill_fraction",s.fill),("scaled_corner",m.corner_scores))
            println(io,"SCALARS $name double 1\nLOOKUP_TABLE default");foreach(p->println(io,p),field)
        end
    end
end

# These are exchanges across the film/reservoir RIM, not pump-line gross flow.
function regular_rim_exchange(film,nt,nz)
    mask=reshape(film.feed_cells,nt,nz);exchange=Float64[]
    for j in 1:nz,i in 1:nt
        next=mod1(i+1,nt)
        mask[i,j]!=mask[next,j] && push!(exchange,(mask[i,j] ? 1. : -1.)*film.circumferential_face_flow_m3_s[i,j])
        j<nz && mask[i,j]!=mask[i,j+1] && push!(exchange,(mask[i,j] ? 1. : -1.)*film.axial_face_flow_m3_s[i,j+1])
    end
    return (sum(max.(exchange,0.)),-sum(min.(exchange,0.)))
end

function run_fitted(out)
    ispath(out) && error("refusing to overwrite fitted-inlet evidence: $out")
    mkpath(out)
    reference=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    displacement=(-9.152501266813789e-7,1.7844427760798323e-7,0.)
    pose=floating_bush_pose(reference,displacement)
    open(joinpath(out,"comparison.tsv"),"w") do io
        println(io,"mesh\tlevel\tcells\tFx_N\tFy_N\tFz_N\tnet_feed_m3_s\tgross_rim_out_m3_s\tgross_rim_in_m3_s\tsource_injection_m3_s\tsource_return_m3_s\tlarge_end_m3_s\tsmall_end_m3_s\tland_peak_MPa\tminimum_fill\titerations\timbalance\tsolve_wall_s\tminimum_corner\tminimum_face_alignment\tmaximum_neighbour_area_ratio")
        for level in (1,2,4)
            nt,nz=128level,64level
            mask=floating_bush_feed(reference,nt,nz)
            t=@timed floating_bush_state(reference,displacement,nt,nz;feed_cells=mask,jfo_max_iterations=150)
            s=t.value;d=s.film.diagnostics
            rimflow=regular_rim_exchange(s.film,nt,nz)
            isapprox(rimflow[1]-rimflow[2],d.net_feed_flow_m3_s;rtol=1e-10) || error("rim/source accounting mismatch")
            println(io,join(("regular",level,nt*nz,s.bush_force_n...,d.net_feed_flow_m3_s,rimflow...,d.gross_feed_inflow_m3_s,
                d.gross_feed_return_m3_s,d.outward_end_flow_m3_s...,maximum(s.film.pressure_pa[.!mask])/1e6,
                d.minimum_fill,d.iterations,d.relative_global_imbalance,t.time,NaN,NaN,NaN),'\t'));flush(io)
            println("regular level=$level force=$(s.bush_force_n) supply=$(d.net_feed_flow_m3_s)");flush(stdout)
            m=fitted_mesh(level)
            t=@timed solve_fitted(m,pose)
            s=t.value;q=m.quality
            println(io,join(("fitted",level,length(m.quads),s.force...,s.net_feed,s.injection,s.return_flow,NaN,NaN,s.outward...,
                s.peak_pressure/1e6,s.minimum_fill,s.iterations,s.imbalance,t.time,q.minimum_scaled_corner,
                q.minimum_face_alignment,q.maximum_neighbour_area_ratio),'\t'));flush(io)
            println("fitted level=$level force=$(s.force) supply=$(s.net_feed) quality=$q");flush(stdout)
            write_fitted_vtk(joinpath(out,"fitted-$nt.vtk"),m,s)
            open(joinpath(out,"patch-$nt.tsv"),"w") do patch
                println(patch,"cell\tx1_mm\tz1_mm\tx2_mm\tz2_mm\tx3_mm\tz3_mm\tx4_mm\tz4_mm\tcorner")
                for (k,xy) in enumerate(m.localnodes)
                    # Exclude the opposite -Y face: its XZ projection overlaps.
                    sin(sum(xy[:,1])/4/m.R)>0 || continue
                    physical=[begin
                        x,y=xy[i,:];z=y*m.ca;rh=m.R+(.05-z)*tand(10.)
                        (1e3*rh*cos(x/m.R),1e3*(z-.05))
                    end for i in 1:4]
                    all(v->abs(v[1])<17 && abs(v[2])<17,physical) || continue
                    println(patch,join((k,Iterators.flatten(physical)...,m.corner_scores[k]),'\t'))
                end
            end
        end
    end
    files=["expts/conical-fitted-inlet/"*s for s in ("mesh.jl","operator.jl","checks.jl","run.jl","plot.py")]
    append!(files,["expts/surface-inlet-hex/mesh.jl","bearing.jl","mesh/2d/conical/grid.jl"])
    append!(files,["simulation/reynolds/"*s for s in ("conical_equilibrium.jl","conical_jfo.jl","conical.jl","jfo_active_set.jl","cylindrical.jl")])
    record=Dict("status"=>"Diagnostic mesh/operator comparison only; NOT an accepted production or 3D Hex8 mesh",
        "displacement_m"=>collect(displacement),"eccentricity_ratio"=>pose.parameters.eccentricity_ratio,
        "radial_clearance_m"=>pose.radial_clearance_m,"minimum_normal_gap_m"=>pose.minimum_normal_gap_m,
        "rpm"=>2000.,"viscosity_pa_s"=>.0277,"feed_diameter_m"=>.005,"feed_pressure_pa_gauge"=>5e5,
        "end_pressure_pa_gauge"=>[0.,0.],"cavitation_pressure_pa_gauge"=>0.,
        "mesh_recipe"=>"retained curved/harmonic layout; q=8*level, outer layers=9*level, fixed half-width=7.5mm; reservoir interior eliminated from pressure unknowns",
        "flow_accounting"=>"Matched gross flows cross the film/reservoir rim; they are NOT gross pump-line injection/return. Original per-cell reservoir source sums reported separately; unavailable after eliminating the fitted reservoir interior.",
        "julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(FITTED_ROOT,f)))) for f in files))
    open(joinpath(out,"comparison.toml"),"w") do io;TOML.print(io,record;sorted=true);end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide a fresh output directory")
    run_fitted(ARGS[1])
end
