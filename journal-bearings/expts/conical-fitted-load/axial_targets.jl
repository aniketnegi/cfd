include(joinpath(@__DIR__,"solve.jl"))
const TARGET_ROOT=normpath(joinpath(@__DIR__,"../.."))

function run_axial_targets(out)
    ispath(out) && error("refusing to overwrite axial target evidence: $out")
    cfg=TOML.parsefile(joinpath(@__DIR__,"axial-targets.toml"))
    @assert (cfg["radial_load_n"],cfg["feed_pressure_pa_gauge"],cfg["feed_diameter_mm"],cfg["rpm"],cfg["viscosity_pa_s"],cfg["exit_and_cavity_pressure_pa_gauge"])==(1000.,5e5,5.,2000.,.0277,0.)
    @assert cfg["jfo_iteration_cap"]==150
    BLAS.set_num_threads(1);mkpath(out)
    prior=joinpath(TARGET_ROOT,"out/cases/conical-floating-bush/runs/fitted-load-continuation/states.tsv")
    lines=split.(readlines(prior),'\t');header=first(lines)
    warm=Dict((parse(Int,r[findfirst(==("level"),header)]),round(parse(Float64,r[findfirst(==("clearance_um"),header)]);digits=8))=>
        (parse(Float64,r[findfirst(==("bx_um"),header)])*1e-6,parse(Float64,r[findfirst(==("by_um"),header)])*1e-6)
        for r in lines[2:end] if r[findfirst(==("radial_converged"),header)]=="true")
    cp(prior,joinpath(out,"prior-states.tsv"))
    sources=["bearing.jl","Project.toml"]
    for directory in ("mesh/2d","simulation/reynolds","expts/conical-fitted-inlet","expts/surface-inlet-hex","expts/conical-fitted-load")
        for (dir,_,files) in walkdir(joinpath(TARGET_ROOT,directory)),file in files
            endswith(file,".jl") || endswith(file,".toml") || endswith(file,".py") || continue
            push!(sources,relpath(joinpath(dir,file),TARGET_ROOT))
        end
    end
    hashes=Dict(f=>bytes2hex(sha256(read(joinpath(TARGET_ROOT,f)))) for f in sources)
    for f in sources
        dest=joinpath(out,"source-snapshot",f);mkpath(dirname(dest));cp(joinpath(TARGET_ROOT,f),dest)
    end
    meta=Dict{String,Any}("case"=>cfg,"status"=>"running","source_sha256"=>hashes,
        "prior_states_sha256"=>bytes2hex(sha256(read(prior))),"julia_version"=>string(VERSION),
        "julia_threads"=>Threads.nthreads(),"startup_blas_interface_threads"=>BLAS.get_num_threads(),
        "thread_note"=>"Interface count is not proof of serial native execution; previous run used additional native workers",
        "root_selection"=>"First sign bracket per target and mesh; no uniqueness claim",
        "scope"=>"Formal steady model only. Large-gap thin-film/inertia validity and dynamic stability NOT established",
        "flow_definitions"=>"Gross rim exchange is not pump-line injection/return; source decomposition unavailable")
    savemeta()=open(io->TOML.print(io,meta;sorted=true),joinpath(out,"run.toml"),"w")
    savemeta();start=time();calls=Ref(0);rows=NamedTuple[];cache=Dict();problems=Dict();targets=NamedTuple[]
    brackets=Dict{Tuple{Int,Float64},Vector{Tuple{Float64,Float64}}}();outcome="failed"
    function seed_for(level,c)
        haskey(warm,(level,c)) && return warm[(level,c)]
        level==2 && haskey(cache,(1,c)) && return cache[(1,c)].state.pose.displacement_m[1:2]
        candidates=[k for k in keys(cache) if k[1]==level]
        if isempty(candidates)
            candidates=[k for k in keys(cache) if k[1]==1]
        end
        isempty(candidates) && error("no admissible continuation seed")
        lower=filter(k->k[2]<=c,candidates);upper=filter(k->k[2]>=c,candidates)
        if !isempty(lower) && !isempty(upper)
            a=lower[argmax(last.(lower))];b=upper[argmin(last.(upper))]
            a[2]==b[2] && return cache[a].state.pose.displacement_m[1:2]
            t=(c-a[2])/(b[2]-a[2])
            return (1-t).*cache[a].state.pose.displacement_m[1:2].+t.*cache[b].state.pose.displacement_m[1:2]
        end
        # Preserve e/c during outward extrapolation; c^3 is unsafe at high e/c.
        key=candidates[argmin([abs(log(c/k[2])) for k in candidates])]
        return cache[key].state.pose.displacement_m[1:2].*(c/key[2])
    end
    try
        open(joinpath(out,"iterations.tsv"),"w") do history
            println(history,"label\titeration\tbx_m\tby_m\tbz_m\tFx_N\tFy_N\tFz_N\tRx_N\tRy_N\tRz_N")
            function solve_at(level,c;target=500.,kind="branch")
                key=(level,Float64(c));haskey(cache,key) && return cache[key]
                calls[]+=1;calls[]<=cfg["radial_call_cap"] || error("radial call budget exhausted")
                time()-start<cfg["wall_cap_seconds"] || error("wall budget exhausted")
                label="$(kind)-L$(level)-c$(replace(string(c),'.'=>'p'))-$(calls[])"
                seed=seed_for(level,c);println("START $label axial residual reference=$target N");flush(stdout)
                result=nothing;elapsed=0.
                try
                    elapsed=@elapsed result=radial_equilibrium(problems[level],c,seed;axial_load_n=target,
                        force_tolerance_n=cfg["force_component_tolerance_n"],max_iterations=cfg["radial_iteration_cap"])
                catch err
                    msg=sprint(showerror,err)
                    if err isa ArgumentError || occursin("active set failed",msg) || occursin("conservative flux gate",msg) || occursin("pressure/content bounds",msg)
                        open(io->println(io,msg),joinpath(out,label*"-exception.txt"),"w")
                        println("FAILED $label: $msg");flush(stdout);return nothing
                    end
                    rethrow()
                end
                p=result.state.pose;s=result.state.film;m=problems[level].mesh;e=hypot(p.displacement_m[1:2]...)
                areas=m.areas.*((m.R .+ (.05 .-m.centre[:,2].*m.ca).*tand(10.))./m.R)
                row=(;label,kind,level,cells=length(m.quads),radial_converged=result.converged,reason=result.reason,
                    axial_residual_reference_N=target,joint_balanced=result.converged && abs(s.force[3]-target)<=.02,
                    clearance_um=p.radial_clearance_m*1e6,bx_um=p.displacement_m[1]*1e6,by_um=p.displacement_m[2]*1e6,bz_mm=p.displacement_m[3]*1e3,
                    eccentricity_um=e*1e6,eccentricity_ratio=p.parameters.eccentricity_ratio,
                    min_normal_gap_um=p.minimum_normal_gap_m*1e6,max_normal_gap_um=(p.radial_clearance_m+e)*m.ca*1e6,
                    Fx_N=s.force[1],Fy_N=s.force[2],Fz_N=s.force[3],Rx_N=result.residual_n[1],Ry_N=result.residual_n[2],G_N=s.force[3]-target,
                    radial_residual_N=hypot(result.residual_n[1:2]...),support_if_axially_held_N=target-s.force[3],
                    land_min_Pa=minimum(s.pressure),land_peak_MPa=s.peak_pressure/1e6,land_mean_MPa=sum(s.pressure.*areas)/sum(areas)/1e6,
                    net_feed_L_min=s.net_feed*60000,rim_out_L_min=s.injection*60000,rim_in_L_min=s.return_flow*60000,
                    large_exit_L_min=s.outward[1]*60000,small_exit_L_min=s.outward[2]*60000,
                    minimum_fill=minimum(s.fill),maximum_fill=maximum(s.fill),underfull_area_fraction=sum(areas[s.fill.<1-1e-8])/sum(areas),
                    local_flux_error=s.local_error,global_imbalance=s.imbalance,
                    large_extension_mm=p.required_large_end_extension_m*1e3,small_extension_mm=p.required_small_end_extension_m*1e3,
                    max_normal_gap_over_mouth_diameter=(p.radial_clearance_m+e)*m.ca/.005,
                    field_evaluations=result.evaluations,final_active_set_iterations=s.iterations,wall_s=elapsed)
                open(joinpath(out,"states.tsv"),isempty(rows) ? "w" : "a") do io
                    isempty(rows) && println(io,join(keys(row),'\t'));println(io,join(values(row),'\t'))
                end
                push!(rows,row)
                for h in result.history;println(history,join((label,h.iteration,h.displacement_m...,h.force_n...,h.residual_n...),'\t'));end
                flush(history)
                open(joinpath(out,label*"-failures.txt"),"w") do io;foreach(s->println(io,s),result.trial_failures);end
                println("END $label radial=$(result.converged) e/c=$(row.eccentricity_ratio) Fz=$(row.Fz_N) N residual=$(row.radial_residual_N) N time=$elapsed s");flush(stdout)
                if result.converged
                    cache[key]=result
                    # Export every solved state, including bisection samples.
                    path=joinpath(out,label*".vtk");write_fitted_vtk(path,m,s)
                    open(path,"a") do io
                        println(io,"SCALARS normal_gap_m double 1\nLOOKUP_TABLE default")
                        for k in eachindex(m.quads);println(io,conical_normal_gap(p.parameters,m.centre[k,1]/m.R-p.eccentricity_angle_rad,m.centre[k,2]*m.ca));end
                    end
                    return result
                end
                return nothing
            end
            coarse_sampled=Float64[]
            for level in cfg["levels"]
                problems[level]=fitted_load_problem(level)
                cs=if level==1
                    cfg["coarse_clearances_um"]
                else
                    # Independent mesh check at every bracket end, plus range ends.
                    sort(unique(vcat([100.,maximum(coarse_sampled)],
                        [x for target in cfg["axial_targets_n"] for pair in get(brackets,(1,target),Tuple{Float64,Float64}[]) for x in pair])))
                end
                completed=Float64[]
                for c in cs
                    s=solve_at(level,c)
                    isnothing(s) && break
                    push!(completed,c);level==1 && push!(coarse_sampled,c)
                end
                isempty(completed) && error("no converged radial state on level $level")
                for target in cfg["axial_targets_n"]
                    found=Tuple{Float64,Float64}[]
                    for i in 1:length(completed)-1
                        a,b=completed[i:i+1];fa=cache[(level,a)].state.bush_force_n[3]-target;fb=cache[(level,b)].state.bush_force_n[3]-target
                        signbit(fa)!=signbit(fb) && push!(found,(a,b))
                    end
                    brackets[(level,target)]=found
                    root_c=NaN;root_Fz=NaN;success=false;reason="no sign bracket in completed samples"
                    for (a,b) in found
                        f=c->begin
                            s=solve_at(level,c;target,kind="root-$(Int(target))")
                            isnothing(s) && error("radial failure inside axial bracket; root not certified")
                            s.state.bush_force_n[3]-target
                        end
                        root=bisect_axial(f,a,b;force_tolerance_n=cfg["force_component_tolerance_n"],max_iterations=cfg["bisection_iteration_cap"])
                        root_c=root.c;root_Fz=root.G+target;success=root.converged
                        reason=success ? "formal steady root; physical validity and stability unverified" : "bisection iteration limit"
                        open(joinpath(out,"root-L$level-target$(Int(target))-$(a)-$(b).toml"),"w") do io
                            TOML.print(io,Dict("converged"=>success,"clearance_um"=>root_c,"Fz_N"=>root_Fz,
                                "history"=>[Dict(string(k)=>v for (k,v) in pairs(h)) for h in root.history]);sorted=true)
                        end
                        break # First detected bracket only; no uniqueness claim.
                    end
                    row=(;level,target_axial_N=target,converged=success,reason,clearance_um=root_c,Fz_N=root_Fz,
                        residual_N=root_Fz-target,lowest_sample_clearance_um=minimum(completed),highest_sample_clearance_um=maximum(completed),
                        minimum_sample_Fz_N=minimum(cache[(level,c)].state.bush_force_n[3] for c in completed),
                        maximum_sample_Fz_N=maximum(cache[(level,c)].state.bush_force_n[3] for c in completed))
                    open(joinpath(out,"targets.tsv"),isempty(targets) ? "w" : "a") do io
                        isempty(targets) && println(io,join(keys(row),'\t'));println(io,join(values(row),'\t'))
                    end
                    push!(targets,row);println("TARGET $row");flush(stdout)
                end
            end
        end
        outcome="completed"
    catch err
        meta["error"]=sprint(showerror,err);rethrow()
    finally
        meta["status"]=outcome;meta["elapsed_wall_s"]=time()-start;meta["radial_calls"]=calls[];meta["saved_states"]=length(rows)
        meta["converged_targets"]=count(r->r.converged,targets)
        meta["source_scope"]="No passage, elasticity, inertia correction or hidden support; same pressure operator. Branch reuse is valid because axial target does not enter the fixed-seat radial equations."
        savemeta()
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide fresh output directory")
    run_axial_targets(ARGS[1])
end
