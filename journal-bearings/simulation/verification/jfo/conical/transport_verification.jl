# Tensor cone verification followed by resolved-rim cylinder/cone transport.
include(joinpath(@__DIR__,"transport_cases.jl"))
include(joinpath(@__DIR__,"transport_sampling.jl"))
include(joinpath(@__DIR__,"../../../../cli.jl"))

function transport_write(io,row;header=false)
    header && println(io,join(string.(keys(row)),'\t'))
    println(io,join(values(row),'\t'));flush(io)
end
function transport_flow_vector(r)
    if hasproperty(r.grid,:surface)
        tags=r.grid.surface.boundary;q=r.face_flow_m3_s
    else
        nt=r.grid.n_theta
        tags=[fill(:large,nt);fill(:small,nt)]
        q=[-r.axial_face_flow_m3_s[:,1];r.axial_face_flow_m3_s[:,end]]
    end
    out=[sum(max.(q[tags.==tag],0.)) for tag in (:large,:small,:feed)]
    return [out;[-sum(min.(q[tags.==tag],0.)) for tag in (:large,:small,:feed)]]
end
function transport_save_field(out,name,r)
    theta,z=transport_centres(r.grid)
    open(joinpath(out,name*".tsv"),"w") do io
        println(io,"theta_rad\tz_m\tpressure_pa\tfill_fraction\tnormal_gap_m")
        for i in eachindex(theta)
            println(io,join((theta[i],z[i],r.pressure_pa[i],r.fill_fraction[i],r.normal_gap_m[i]),'\t'))
        end
    end
    serialize(joinpath(out,name*".jls"),r)
    write_conical_jfo_vtk(joinpath(out,name*".vtk"),r)
end
function transport_fd_comparison(angle,nt;rpm=2000.)
    a=ConicalTransportReference.assemble(angle,.6,nt,nt÷2;rpm)
    r=ConicalTransportReference.solve(a)
    fv=transport_solve(transport_mesh(angle,nt);rpm)
    axes=full_film_axes(fv.grid)
    samples=[transport_sample(fv,t,z;axes) for t in a.theta,z in a.z[2:end-1]]
    P,F=first.(samples),last.(samples)
    ep,ef=r.pressure[:,2:end-1],r.fill[:,2:end-1]
    weight=[(.05005+(.05-z)*tand(angle))/cosd(angle) for t in a.theta,z in a.z[2:end-1]]
    boundary=Float64[]
    for (j,z) in enumerate((0.,.1)),i in 1:nt
        derivative=j==1 ? (-3r.pressure[i,1]+4r.pressure[i,2]-r.pressure[i,3])/(2a.dz) :
            (3r.pressure[i,end]-4r.pressure[i,end-1]+r.pressure[i,end-2])/(2a.dz)
        g=ConicalTransportReference.geometry(a.theta[i],z,angle,.6;pose=3pi/2)
        push!(boundary,(j==1 ? 1 : -1)*g.K*g.rh*g.ca*derivative*a.dt)
    end
    return (;angle_deg=angle,nt,nz=nt÷2,rpm,
        pressure_relative_L2=sqrt(sum((P-ep).^2 .*weight)/sum(ep.^2 .*weight)),
        pressure_relative_Linf=maximum(abs,P-ep)/maximum(abs,ep),
        fill_mean_absolute=sum(abs.(F-ef).*weight)/sum(weight),minimum_fill=minimum(r.fill),
        FD_signed_net_m3_s=sum(boundary),FD_relative_net=abs(sum(boundary))/sum(abs,boundary),
        sweeps=r.sweeps,residual=r.residual)
end
function transport_metadata(out,full)
    ispath(out) && error("refusing to overwrite $out")
    mkpath(out)
    paths=["bearing.jl","cli.jl","Project.toml","Manifest.toml","simulation/tests/fitted_jfo_checks.jl",
        "simulation/verification/swift-stieber/fitted_reference.jl"]
    for base in ("simulation/reynolds/src","mesh","simulation/verification/jfo/conical")
        for (dir,_,files) in walkdir(joinpath(TRANSPORT_ROOT,base)),file in files
            endswith(file,".jl") && push!(paths,relpath(joinpath(dir,file),TRANSPORT_ROOT))
        end
    end
    inputs=Dict("scope"=>"Prescribed-eccentricity conical JFO verification and resolved7mm inlet",
        "julia_version"=>string(VERSION),"full"=>full,"physical_validation"=>false,
        "journal_midpoint_radius_m"=>.05,"axial_length_m"=>.1,"radial_clearance_m"=>50e-6,
        "mu_pa_s"=>.0277,"eccentricity_angle_rad"=>3pi/2,"rpm_magnitude"=>2000.,
        "fitted_feed_diameter_m"=>.007,"fitted_feed_axial_m"=>.05,"fitted_feed_angle_rad"=>pi/2,
        "fitted_feed_pressure_gauge_pa"=>5e5,"fitted_end_pressure_gauge_pa"=>[0.,0.],
        "cavitation_pressure_gauge_pa"=>0.,"tensor_end_pressure_gauge_pa"=>[3e5,3e5],
        "supply_status"=>"7mm maximum-gap rim and0.5MPa are explicit design choices, not recovered rig facts; tensor end supply is a separate numerical fixture",
        "source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(TRANSPORT_ROOT,p)))) for p in unique(paths)))
    open(joinpath(out,"inputs.toml"),"w") do io;TOML.print(io,inputs;sorted=true);end
end
function verify_conical_transport(out;full=false,observer=nothing)
    transport_metadata(out,full)
    tensor_grids=full ? (64,128,256,512) : (64,128,256)
    fitted_grids=full ? (128,256,512) : (128,256)
    meshes=Dict{Tuple{String,Float64,Int},Any}()
    @testset "Conical and fitted JFO transport" begin
        for angle in (0.,10.,20.),nt in tensor_grids
            meshes[("tensor",angle,nt)]=transport_mesh(angle,nt)
        end
        for angle in (0.,10.),nt in fitted_grids
            println("mesh fitted angle=$angle grid=$nt × $(nt÷2)")
            meshes[("fitted",angle,nt)]=transport_mesh(angle,nt;fitted=true,observer)
            yield()
        end
        open(joinpath(out,"physical.tsv"),"w") do io
            first=true
            for layout in ("tensor","fitted"),angle in (layout=="tensor" ? (0.,10.,20.) : (0.,10.)),
                ecc in (layout=="tensor" ? (.6,.8) : (.2,.6,.8)),rpm in (layout=="tensor" ? (2000.,) : (-2000.,2000.))
                previous=nothing;previous_flow=nothing;previous_force=nothing
                grids=layout=="tensor" ? tensor_grids : fitted_grids
                for nt in grids
                    g=meshes[(layout,angle,nt)];elapsed=@elapsed r=transport_solve(g;ecc,rpm)
                    flow=transport_flow_vector(r);force=collect(housing_pressure_force(r.grid,r.pressure_pa))
                    e=isnothing(previous) ? (;pressure_change=NaN,fill_change=NaN,common_land_samples=0) :
                        layout=="tensor" ? (;transport_field_change(previous,r)...,common_land_samples=length(previous.pressure_pa)) : transport_profile_change(previous,r)
                    fc=isnothing(previous) ? NaN : norm(flow-previous_flow)/norm(flow)
                    wc=isnothing(previous) ? NaN : norm(force-previous_force)/norm(force)
                    d=r.diagnostics
                    row=(;layout,angle_deg=angle,epsilon=ecc,rpm,nt,nz=nt÷2,cells=length(g.quads),
                        peak_pa=maximum(r.pressure_pa),minimum_fill=minimum(r.fill_fraction),
                        force_x_n=force[1],force_y_n=force[2],force_z_n=force[3],
                        gross_out_large_m3_s=flow[1],gross_out_small_m3_s=flow[2],gross_out_feed_m3_s=flow[3],
                        gross_in_large_m3_s=flow[4],gross_in_small_m3_s=flow[5],gross_in_feed_m3_s=flow[6],
                        net_feed_m3_s=d.net_feed_flow_m3_s,pressure_change=e.pressure_change,
                        fill_change=e.fill_change,common_land_samples=e.common_land_samples,
                        flow_change=fc,force_change=wc,local_balance=d.relative_max_cell_imbalance,
                        global_balance=d.relative_global_imbalance,complementarity=d.complementarity_residual,
                        iterations=d.iterations,elapsed_s=elapsed)
                    transport_write(io,row;header=first);first=false
                    @test minimum(r.pressure_above_cavitation_pa)>=0
                    @test 0<=minimum(r.fill_fraction)<=maximum(r.fill_fraction)<=1
                    @test d.relative_max_cell_imbalance<1e-10 && d.relative_global_imbalance<1e-10
                    @test d.complementarity_residual<1e-10
                    if nt==last(grids)
                        @test e.pressure_change<.01 && e.fill_change<.01
                        @test fc<.02 && wc<.01
                        if layout=="fitted"
                            open(joinpath(out,"profiles-$angle-$ecc-$rpm.tsv"),"w") do stream
                                for (i,p) in enumerate(transport_profiles(r));transport_write(stream,p;header=i==1);end
                            end
                            if rpm>0 && ecc>=.6
                                transport_save_field(out,"fitted-$angle-$ecc",r)
                            end
                        end
                    end
                    println("physical $layout angle=$angle epsilon=$ecc rpm=$rpm nt=$nt peak=$(maximum(r.pressure_pa)) fill=$(minimum(r.fill_fraction))")
                    previous,previous_flow,previous_force=r,flow,force
                    yield()
                end
            end
        end
        open(joinpath(out,"manufactured.tsv"),"w") do io
            first=true
            for layout in ("tensor","fitted"),angle in (layout=="tensor" ? (0.,10.,20.) : (0.,10.)),rpm in (-2000.,2000.)
                previous=nothing
                for nt in (layout=="tensor" ? (64,128,256) : fitted_grids)
                    g=meshes[(layout,angle,nt)];m=transport_mms(g,rpm)
                    transport_write(io,(;layout,angle_deg=angle,rpm,nt,nz=nt÷2,m...);header=first);first=false
                    @test m.minimum_pressure>=0 && 0<=m.minimum_fill<=m.maximum_fill<=1
                    @test m.complementarity==0 && m.local_balance<1e-10
                    if !isnothing(previous)
                        @test m.pressure_relative_L2<.8previous.pressure_relative_L2
                        @test m.fill_mean_absolute<.7previous.fill_mean_absolute
                    end
                    previous=m
                    println("MMS $layout angle=$angle rpm=$rpm nt=$nt pressure=$(m.pressure_relative_L2) fill=$(m.fill_mean_absolute)")
                    yield()
                end
                @test previous.pressure_relative_L2<.003 && previous.fill_mean_absolute<.003
                @test previous.boundary_flow_relative_L2<.01 && previous.force_relative_L2<.003
            end
            # Additional graded tensor coordinates, both rotation signs.
            for rpm in (-2000.,2000.),nt in (64,128,256)
                m=transport_mms(transport_mesh(10.,nt;graded=true),rpm)
                transport_write(io,(;layout="graded",angle_deg=10.,rpm,nt,nz=nt÷2,m...))
                @test m.local_balance<1e-10 && m.complementarity==0
            end
        end
        open(joinpath(out,"independent.tsv"),"w") do io
            first=true
            for angle in (10.,20.)
                previous=nothing
                for nt in (64,128,256)
                    row=transport_fd_comparison(angle,nt)
                    transport_write(io,row;header=first);first=false
                    @test 0<=row.minimum_fill<=1 && row.residual<1e-11
                    isnothing(previous) || @test row.pressure_relative_L2<.6previous.pressure_relative_L2
                    previous=row
                    println("independent angle=$angle nt=$nt pressure=$(row.pressure_relative_L2) fill=$(row.fill_mean_absolute)")
                end
                @test previous.pressure_relative_L2<.003 && previous.fill_mean_absolute<.006
            end
        end
        g=meshes[("fitted",10.,128)]
        a=assemble_conical_jfo(g;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),
            cavitation_pressure_pa=0.,feed_pressure_pa=5e5,eccentricity_angle_rad=3pi/2)
        policy=solve_jfo_active_set(a.A,a.B,a.b,a.pressure_scale_pa)
        open(joinpath(out,"algebraic.tsv"),"w") do io
            for (i,(relaxation,tolerance)) in enumerate(((1.,1e-12),(.8,3e-13)))
                r=ConicalTransportReference.projected_gs(a.A,a.B,a.b,a.pressure_scale_pa;relaxation,tolerance)
                pe=norm(r.pressure-policy.pressure)/norm(policy.pressure)
                fe=norm(r.fill-(1 .-policy.void),Inf)
                @test pe<1e-6 && fe<1e-5
                transport_write(io,(;relaxation,tolerance,pressure_relative_L2=pe,fill_Linf=fe,sweeps=r.sweeps,residual=r.residual);header=i==1)
            end
        end
        # Continuous cylinder limit without a new feed/inventory assumption.
        cylinder=transport_solve(meshes[("tensor",0.,64)])
        open(joinpath(out,"cylindrical-limit.tsv"),"w") do io
            previous=Inf
            for (i,angle) in enumerate((.1,.01,.001,.0001,.00001))
                r=transport_solve(transport_mesh(angle,64))
                pe=norm(r.pressure_pa-cylinder.pressure_pa)/norm(cylinder.pressure_pa)
                fe=sum(abs,r.fill_fraction-cylinder.fill_fraction)/length(r.fill_fraction)
                @test pe<previous;previous=pe
                transport_write(io,(;angle_deg=angle,pressure_relative_L2=pe,fill_mean_absolute=fe);header=i==1)
            end
            @test previous<1e-6
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,description="Verify shared conical and fitted-rim JFO liquid transport.") do cfg
        observer=stderr isa Base.TTY && !cfg.quiet && !cfg.verbose ? _->yield() : nothing
        verify_conical_transport(cfg.output;full=cfg.full,observer)
    end)
end
