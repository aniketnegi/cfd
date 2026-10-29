include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"q4.jl"))
include(joinpath(SOURCE,"simulation/reynolds/src/pressure_vtk.jl"))

function fv_midplane(result)
    grid=result.grid;nt,nz=grid.n_theta,grid.n_axial
    z=(grid.z_nodes[1:end-1]+grid.z_nodes[2:end]).*.0005
    P=reshape(result.pressure_pa,nt,nz);target=grid.parameters.bearing_length*.0005
    lower=searchsortedlast(z,target)
    if z[lower]==target
        profile=P[:,lower]
    else
        weight=(target-z[lower])/(z[lower+1]-z[lower])
        profile=(1-weight).*P[:,lower]+weight.*P[:,lower+1]
    end
    theta=(grid.theta_nodes+[grid.theta_nodes[2:end];grid.theta_nodes[1]+2pi])/2
    return (;theta,pressure=profile,nz,alpha0_max=maximum(abs,(P[1,:]+P[end,:])/2))
end

function record_agreement(output,ecc,fe,fv)
    sample=q4_periodic_sample(fv.theta,fv.pressure,fe.a.grid.theta_nodes)
    relative_L2=norm(fe.midplane_pressure_pa-sample)/norm(sample)
    relative_Linf=maximum(abs,fe.midplane_pressure_pa-sample)/maximum(abs,sample)
    relative_peak_difference=maximum(fe.midplane_pressure_pa)/maximum(fv.pressure)-1
    path=joinpath(output,"unported-agreement.tsv");header=!ispath(path)
    open(path,"a") do io
        row!(io,(;epsilon=ecc,Q4_nt=fe.a.grid.n_theta,Q4_nz=fe.a.grid.n_axial,
            FV_nt=length(fv.theta),FV_nz=fv.nz,relative_L2,relative_Linf,relative_peak_difference);header)
    end
    @test relative_L2<.001 && relative_Linf<.002 && abs(relative_peak_difference)<.001
end

function export_q4_coarse(output,grid,fe,fv,ecc)
    nt,nz=grid.n_theta,grid.n_axial
    paper_node(k)=((k-1)%nt)*(nz+1)+(nz-(k-1)÷nt+1)
    open(joinpath(output,"q4-nodes-e$(ecc).tsv"),"w") do io
        println(io,"public_node\tpaper_node\ttheta_deg\tz_mm\tx_mm\ty_mm\tpressure_pa\tnormal_gap_m")
        for k in eachindex(grid.nodes)
            println(io,join((k,paper_node(k),rad2deg(fe.a.theta[k]),fe.a.z[k]*1e3,
                grid.nodes[k][1],grid.nodes[k][2],fe.pressure_pa[k],fe.a.hnode[k]),'\t'))
        end
    end
    open(joinpath(output,"q4-elements-e$(ecc).tsv"),"w") do io
        println(io,"public_element\tpaper_element\tpublic_n1\tpublic_n2\tpublic_n3\tpublic_n4\tpaper_n1\tpaper_n2\tpaper_n3\tpaper_n4")
        for (cell,ids) in enumerate(grid.quads)
            i,j=mod1(cell,nt),(cell-1)÷nt+1
            println(io,join((cell,(i-1)*nz+nz-j+1,ids...,map(paper_node,ids)...),'\t'))
        end
    end
    path=joinpath(output,"q4-nodal-e$(ecc).vtk")
    write_vtk(path,fv.grid)
    open(path,"a") do io
        println(io,"POINT_DATA $(length(grid.nodes))\nSCALARS nodal_pressure_pa double 1\nLOOKUP_TABLE default")
        foreach(p->println(io,p),fe.pressure_pa)
    end
    write_conical_pressure_vtk(joinpath(output,"fv-cell-e$(ecc).vtk"),fv)
end

function unported_grid_study(output;full=false)
    fresh_output(output;suite="unported-Q4-FV-grid-comparison")
    reference=TOML.parsefile(joinpath(SOURCE,"expts/gangrade-figure6/reference.toml"))
    targets=reference["points"]
    full || (targets=filter(p->p["eccentricity_ratio"] in (.6,.8),targets))
    scale=reference["inputs"]["pressure_scale_pa"]
    y0=reference["reading"]["y_pbar_0_px"];height=y0-reference["reading"]["y_pbar_100_px"]
    levels=(12,24,48,96,192,384)
    inputs=TOML.parsefile(joinpath(output,"inputs.toml"))
    inputs["reference_sha256"]=bytes2hex(sha256(read(joinpath(SOURCE,"expts/gangrade-figure6/reference.toml"))))
    inputs["unported_protocol"]=Dict("end_pressure_gauge_pa"=>[0.,0.],"cavitation_pressure_gauge_pa"=>0.,
        "eccentricity_angle_rad"=>3pi/2,"angular_start_rad"=>0.,"angular_boundary"=>"periodic; no additional pinned line",
        "inlet"=>"none; published0.5MPa used only to normalize reported pressure",
        "Q4_gap"=>"bilinear nodal interpolation of the unchanged production gap",
        "quadrature"=>"2x2 Gauss, chosen reconstruction; authors' exact quadrature unknown",
        "coordinates"=>"public uniform housing_grid theta/z; original coordinate file unavailable",
        "refinement_counts"=>["$(n)x$(n÷3)" for n in levels],"additional_FV_counts"=>"512x256 at epsilon.6/.8")
    inputs["source_sha256"]["simulation/reynolds/src/pressure_vtk.jl"]=bytes2hex(sha256(read(joinpath(SOURCE,"simulation/reynolds/src/pressure_vtk.jl"))))
    inputs["source_sha256"]["mesh/vtk.jl"]=bytes2hex(sha256(read(joinpath(SOURCE,"mesh/vtk.jl"))))
    open(joinpath(output,"inputs.toml"),"w") do io
        TOML.print(io,inputs;sorted=true)
    end
    mkpath(joinpath(output,"profiles"))
    @testset "Unported public mesh Q4/FV coarse-to-fine comparison" begin
        open(joinpath(output,"unported-grid.tsv"),"w") do io
            first=true
            open(joinpath(output,"unported-sampling.tsv"),"w") do sampling
                sample_first=true
                for point in targets
                    ecc=point["eccentricity_ratio"];p=cone(;ecc)
                    fea=100*(y0-point["fea_y_px"])/height;cfd=100*(y0-point["cfd_y_px"])/height
                    fe_fine=nothing;fv_fine=nothing;coarse=nothing
                    for nt in levels
                        nz=nt÷3;grid=housing_grid(p,nt,nz)
                        fv=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=MU,eccentricity_angle_rad=3pi/2)
                        fe=solve_q4(grid;rpm=2000.,viscosity_pa_s=MU,eccentricity_angle_rad=3pi/2)
                        mid=fv_midplane(fv)
                        for (method,profile,theta,dofs,free,iterations,residual,kkt,zero_line) in (
                            ("FV",mid.pressure,mid.theta,length(grid.quads),length(grid.quads),fv.diagnostics.iterations,
                             fv.diagnostics.positive_pressure_relative_residual,fv.diagnostics.complementarity_relative,mid.alpha0_max),
                            ("Q4",fe.midplane_pressure_pa,grid.theta_nodes,length(grid.nodes),count(.!fe.a.fixed),fe.iterations,
                             fe.positive_pressure_relative_residual,fe.complementarity,maximum(abs,fe.pressure_pa[1:nt:end])))
                            peak=maximum(profile);angle=theta[argmax(profile)]
                            row!(io,(;method,epsilon=ecc,nt,nz,elements=length(grid.quads),pressure_dofs=dofs,free_pressure_dofs=free,
                                midplane_maximum_pa=peak,midplane_pbar=peak/scale,peak_theta_deg=rad2deg(angle),
                                fea_pbar=fea,cfd_pbar=cfd,alpha0_max_pa=zero_line,
                                iterations,positive_pressure_residual=residual,complementarity=kkt);header=first)
                            first=false
                            @test residual<1e-10 && kkt<1e-10
                            open(joinpath(output,"profiles","$(method)-e$(ecc)-$(nt)x$(nz).tsv"),"w") do profile_io
                                println(profile_io,"theta_rad\tpressure_pa")
                                for (t,value) in zip(theta,profile)
                                    println(profile_io,"$t\t$value")
                                end
                            end
                            @printf("%s eps=%g mesh=%dx%d dofs=%d midplane=%.9g MPa alpha0=%.3g Pa\n",
                                method,ecc,nt,nz,dofs,peak/1e6,zero_line);flush(stdout)
                        end
                        if nt==12
                            @test length(grid.quads)==48 && length(grid.nodes)==60
                            @test fe.a.grid===grid
                            coarse=(;theta=copy(grid.theta_nodes),peak=maximum(fe.midplane_pressure_pa))
                            ecc in (.6,.8) && export_q4_coarse(output,grid,fe,fv,ecc)
                        end
                        fe_fine=fe;fv_fine=mid
                        GC.gc();yield()
                    end
                    for (method,theta,profile) in (("Q4",fe_fine.a.grid.theta_nodes,fe_fine.midplane_pressure_pa),
                            ("FV",fv_fine.theta,fv_fine.pressure))
                        samples=q4_periodic_sample(theta,profile,coarse.theta)
                        row!(sampling,(;epsilon=ecc,method,coarse_solved_Q4_peak_pa=coarse.peak,
                            fine_peak_pa=maximum(profile),fine_sampled_at_12_nodes_peak_pa=maximum(samples));header=sample_first)
                        sample_first=false
                    end
                    record_agreement(output,ecc,fe_fine,fv_fine)
                    if ecc in (.6,.8)
                        fv=solve_swift_stieber(housing_grid(p,512,256);rpm=2000.,viscosity_pa_s=MU,eccentricity_angle_rad=3pi/2)
                        mid=fv_midplane(fv)
                        row!(io,(;method="FV",epsilon=ecc,nt=512,nz=256,elements=512*256,
                            pressure_dofs=512*256,free_pressure_dofs=512*256,midplane_maximum_pa=maximum(mid.pressure),
                            midplane_pbar=maximum(mid.pressure)/scale,peak_theta_deg=rad2deg(mid.theta[argmax(mid.pressure)]),
                            fea_pbar=fea,cfd_pbar=cfd,alpha0_max_pa=mid.alpha0_max,
                            iterations=fv.diagnostics.iterations,
                            positive_pressure_residual=fv.diagnostics.positive_pressure_relative_residual,
                            complementarity=fv.diagnostics.complementarity_relative))
                        record_agreement(output,ecc,fe_fine,mid)
                    end
                end
            end
        end
    end
end
include(joinpath(@__DIR__,"../../../cli.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,description="Compare unported Q4 and finite-volume pressure on public coarse/refined meshes.") do cfg
        unported_grid_study(cfg.output;full=cfg.full)
    end)
end
