include(joinpath(@__DIR__,"solve.jl"))
const ROOT=normpath(joinpath(@__DIR__,"../.."))

function run_continuation(out)
    ispath(out) && error("refusing to overwrite fitted load evidence: $out")
    BLAS.set_num_threads(1)
    config=TOML.parsefile(joinpath(@__DIR__,"case.toml"));n=config["numerics"]
    # This case-specific adapter intentionally exposes no physical-parameter tuning.
    @assert config["physical"]["external_force_n"]==[0.,1000.,-200.]
    @assert config["physical"]["rpm"]==2000. && config["physical"]["viscosity_pa_s"]==.0277
    @assert Tuple(config["physical"][k] for k in ("journal_midpoint_radius_mm","bearing_length_mm",
        "reference_radial_clearance_um","semicone_angle_deg","feed_diameter_mm","feed_z_mm",
        "feed_pressure_pa_gauge","exit_pressure_pa_gauge","cavity_pressure_pa_gauge")) ==
        (50.,100.,50.,10.,5.,50.,5e5,0.,0.)
    @assert config["physical"]["feed_location"]=="+Y; housing fixed"
    @assert n["jfo_iteration_cap"]==150
    mkpath(out)
    sources=String["bearing.jl","Project.toml"]
    for directory in ("mesh/2d/conical","simulation/reynolds","expts/conical-fitted-inlet",
                      "expts/surface-inlet-hex","expts/conical-fitted-load")
        for (dir,_,files) in walkdir(joinpath(ROOT,directory)),file in files
            endswith(file,".jl") || endswith(file,".toml") || endswith(file,".py") || continue
            push!(sources,relpath(joinpath(dir,file),ROOT))
        end
    end
    hashes=Dict(f=>bytes2hex(sha256(read(joinpath(ROOT,f)))) for f in sources)
    for file in sources
        dest=joinpath(out,"source-snapshot",file);mkpath(dirname(dest));cp(joinpath(ROOT,file),dest)
    end
    meta=Dict{String,Any}("case"=>config,"source_sha256"=>hashes,"status"=>"running",
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "coordinates"=>"VTK in fixed housing-local frame, millimetres; displacements in states.tsv are relative to located shaft",
        "flow_definitions"=>"rim_out/rim_in are film-reservoir perimeter exchanges, NOT pump injection/return; eliminated reservoir source decomposition unavailable",
        "method"=>"Existing damped transverse Newton with explicit FittedLoadProblem dispatch; bracketed bisection only after converged radial states change axial residual sign",
        "reference_clearance_um"=>50.,"steady_root_is_not_dynamic_stability"=>true)
    save_meta()=open(io->TOML.print(io,meta;sorted=true),joinpath(out,"run.toml"),"w")
    save_meta();rows=NamedTuple[];solutions=Dict{Tuple{Int,Float64},Any}();calls=Ref(0)
    start=time();outcome="failed";problems=Dict();brackets=Dict{Int,Vector{Tuple{Float64,Float64}}}()
    try
        open(joinpath(out,"iterations.tsv"),"w") do history
            println(history,"label\titeration\tbx_m\tby_m\tbz_m\tFx_N\tFy_N\tFz_N\tRx_N\tRy_N\tRz_N\tscaled_norm")
            function solve_at(level,c,seed;kind="branch",fd=n["finite_difference_step"])
                calls[]+=1;calls[]<=n["radial_call_cap"] || error("radial-call budget exhausted")
                time()-start<n["wall_cap_seconds"] || error("wall-time budget exhausted")
                problem=problems[level]
                label="$(kind)-L$(level)-c$(replace(string(c),'.'=>'p'))-$(calls[])"
                println("START $label seed=$seed");flush(stdout)
                t=@timed radial_equilibrium(problem,c,seed;fd_step=fd,
                    force_tolerance_n=n["force_component_tolerance_n"],max_iterations=n["radial_iteration_cap"])
                result=t.value;state=result.state;p=state.pose;s=state.film;m=problem.mesh
                physical_areas=m.areas.*((m.R .+ (.05 .-m.centre[:,2].*m.ca).*tand(10.))./m.R)
                e=hypot(p.displacement_m[1:2]...)
                row=(;label,kind,level,cells=length(m.quads),radial_converged=result.converged,
                    joint_balanced=result.converged && abs(result.residual_n[3])<=n["force_component_tolerance_n"],reason=result.reason,
                    bx_um=p.displacement_m[1]*1e6,by_um=p.displacement_m[2]*1e6,bz_mm=p.displacement_m[3]*1e3,
                    clearance_um=p.radial_clearance_m*1e6,eccentricity_um=e*1e6,eccentricity_ratio=p.parameters.eccentricity_ratio,
                    min_normal_gap_um=p.minimum_normal_gap_m*1e6,max_normal_gap_um=(p.radial_clearance_m+e)*m.ca*1e6,
                    Fx_N=s.force[1],Fy_N=s.force[2],Fz_N=s.force[3],Rx_N=result.residual_n[1],Ry_N=result.residual_n[2],G_N=result.residual_n[3],
                    radial_residual_N=hypot(result.residual_n[1:2]...),support_if_axially_held_N=-result.residual_n[3],
                    land_peak_MPa=s.peak_pressure/1e6,land_mean_MPa=sum(s.pressure.*physical_areas)/sum(physical_areas)/1e6,
                    land_min_Pa=minimum(s.pressure),net_feed_L_min=s.net_feed*60000,rim_out_L_min=s.injection*60000,rim_in_L_min=s.return_flow*60000,
                    large_exit_L_min=s.outward[1]*60000,small_exit_L_min=s.outward[2]*60000,
                    minimum_fill=minimum(s.fill),maximum_fill=maximum(s.fill),
                    underfull_area_fraction=sum(physical_areas[s.fill.<1-1e-8])/sum(physical_areas),
                    local_flux_error=s.local_error,global_imbalance=s.imbalance,
                    large_extension_mm=p.required_large_end_extension_m*1e3,small_extension_mm=p.required_small_end_extension_m*1e3,
                    force_concentration=tand(10.)*hypot(s.force[1:2]...)/s.force[3],
                    field_evaluations=result.evaluations,final_active_set_iterations=s.iterations,wall_s=t.time,fd_step=fd)
                open(joinpath(out,"states.tsv"),isempty(rows) ? "w" : "a") do io
                    isempty(rows) && println(io,join(keys(row),'\t'));println(io,join(values(row),'\t'))
                end
                push!(rows,row)
                for h in result.history
                    println(history,join((label,h.iteration,h.displacement_m...,h.force_n...,h.residual_n...,h.scaled_norm),'\t'))
                end
                flush(history)
                open(joinpath(out,label*"-failures.txt"),"w") do io
                    foreach(x->println(io,x),result.trial_failures)
                end
                if result.converged
                    solutions[(level,Float64(c))]=result
                    # Actual pressure/content on land cells; aperture interior is absent.
                    path=joinpath(out,label*".vtk");write_fitted_vtk(path,m,s)
                    gaps=[conical_normal_gap(p.parameters,m.centre[k,1]/m.R-p.eccentricity_angle_rad,m.centre[k,2]*m.ca) for k in eachindex(m.quads)]
                    open(path,"a") do io
                        println(io,"SCALARS normal_gap_m double 1\nLOOKUP_TABLE default");foreach(h->println(io,h),gaps)
                    end
                    if level==4 && c==50.
                        open(joinpath(out,"fixed-seat-fields.tsv"),"w") do io
                            println(io,"theta_deg\tz_mm\tpressure_MPa\tfill\tnormal_gap_um\tarea_m2")
                            for k in eachindex(m.quads)
                                println(io,join((m.centre[k,1]/m.R*180/pi,m.centre[k,2]*m.ca*1000,s.pressure[k]/1e6,s.fill[k],gaps[k]*1e6,physical_areas[k]),'\t'))
                            end
                        end
                    end
                end
                @printf("END %s converged=%s c=%.6g um e/c=%.8g Fz=%.6f G=%.6f N radial residual=%.6g N time=%.2fs evals=%d\n",
                    label,string(result.converged),c,row.eccentricity_ratio,row.Fz_N,row.G_N,row.radial_residual_N,t.time,result.evaluations);flush(stdout)
                return result
            end
            function seed_for(level,c)
                # Prefer the same-clearance solution on the preceding mesh.
                if haskey(solutions,(level÷2,Float64(c)))
                    return solutions[(level÷2,Float64(c))].state.pose.displacement_m[1:2]
                end
                keys_here=filter(k->k[1]==level,collect(keys(solutions)))
                isempty(keys_here) && return (-9.152501266813789e-7,1.7844427760798323e-7).*(c/50)^3
                k=keys_here[argmin([abs(log(c/key[2])) for key in keys_here])]
                return solutions[k].state.pose.displacement_m[1:2].*(c/k[2])^3
            end
            for level in (1,2,4)
                problems[level]=fitted_load_problem(level) # fixed housing built ONCE per resolution
                clearances=level==4 ? sort(unique(vcat(n["fine_anchor_clearances_um"],
                    [x for pair in get(brackets,2,Tuple{Float64,Float64}[]) for x in pair]))) : n["branch_clearances_um"]
                order=vcat([50.],sort(filter(<(50.),clearances);rev=true),sort(filter(>(50.),clearances)))
                for c in order
                    solve_at(level,c,seed_for(level,c))
                end
                brackets[level]=Tuple{Float64,Float64}[]
                for i in 1:length(clearances)-1
                    a,b=Float64.(clearances[i:i+1])
                    haskey(solutions,(level,a)) && haskey(solutions,(level,b)) || continue
                    fa,fb=(solutions[(level,x)].residual_n[3] for x in (a,b))
                    signbit(fa)!=signbit(fb) && push!(brackets[level],(a,b))
                end
                for (a,b) in brackets[level]
                    function residual(c)
                        if !haskey(solutions,(level,Float64(c)))
                            result=solve_at(level,c,seed_for(level,c);kind="axial-root")
                            result.converged || error("radial solve failed inside axial bracket; root not certified")
                        end
                        return solutions[(level,Float64(c))].residual_n[3]
                    end
                    root=bisect_axial(residual,a,b;force_tolerance_n=n["force_component_tolerance_n"],max_iterations=n["bisection_iteration_cap"])
                    open(joinpath(out,"root-L$level-$(a)-$(b).toml"),"w") do io
                        TOML.print(io,Dict("converged"=>root.converged,"clearance_um"=>root.c,"axial_residual_N"=>root.G,
                            "history"=>[Dict(string(k)=>v for (k,v) in pairs(h)) for h in root.history]);sorted=true)
                    end
                end
            end
            # Check FD sensitivity and return continuation at an already solved point.
            if haskey(solutions,(2,12.5)) && haskey(solutions,(2,50.))
                low=solutions[(2,12.5)].state.pose.displacement_m[1:2].*(50/12.5)^3
                solve_at(2,50.,low;kind="return-check",fd=n["finite_difference_step"]/2)
            end
        end
        outcome="completed"
    catch err
        meta["error"]=sprint(showerror,err);rethrow()
    finally
        meta["status"]=outcome;meta["elapsed_wall_s"]=time()-start;meta["radial_calls"]=calls[]
        meta["completed_states"]=length(rows)
        meta["all_saved_radial_states_converged"]=all(r->r.radial_converged,rows)
        meta["joint_balanced_states"]=count(r->r.joint_balanced,rows)
        meta["sign_brackets"]=Dict(string(k)=>[collect(p) for p in pairs] for (k,pairs) in brackets)
        meta["conclusion_scope"]="Sampled branch and tested clearances only; no sign bracket does not prove global nonexistence or uniqueness. No dynamic stability or EHD result."
        save_meta()
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide a fresh output directory")
    run_continuation(ARGS[1])
end
