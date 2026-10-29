include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(ROOT,"cli.jl"))

function fresh_xi_output(path)
    ispath(path) && error("refusing overwrite: $path")
    mkpath(path)
    paths=["bearing.jl","cli.jl","Project.toml","Manifest.toml","mesh/2d/grid.jl","mesh/2d/quality.jl",
        "mesh/2d/metrics.jl","mesh/vtk.jl","simulation/verification/command.py",
        "simulation/tests/xi2024_jfo_checks.jl"]
    append!(paths,["simulation/reynolds/src/"*name for name in readdir(joinpath(ROOT,"simulation/reynolds/src")) if endswith(name,".jl")])
    append!(paths,["simulation/verification/jfo/cylindrical/"*name for name in ("cases.jl","independent.jl")])
    append!(paths,[relpath(joinpath(@__DIR__,name),ROOT) for name in readdir(@__DIR__)
        if any(ext->endswith(name,ext),(".jl",".py",".csv",".toml"))])
    reference=TOML.parsefile(joinpath(@__DIR__,"reference.toml"))
    metadata=Dict("scope"=>"Xi2024 section3.3 finite-cylinder JFO pressure comparison",
        "julia_version"=>string(VERSION),"physical_validation"=>false,
        "geometry_convention_fully_recovered"=>false,"gap_origin_status"=>reference["fixture_assumptions"]["source_limit"],
        "inputs"=>Dict(string(k)=>v for (k,v) in pairs(XI)),"rpm"=>XI_RPM,
        "n_theta"=>[64,128,256,512,1024],"n_axial"=>[16,32,64,128,256],
        "source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in paths))
    open(joinpath(path,"inputs.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
end

function xi_fd_comparison(fv,fd)
    a=fd.a;axes=full_film_axes(fv.grid)
    locations=(;theta_centres=a.theta,z_centres=(a.zeta.+a.aspect).*XI.radius_m)
    pressure=Float64[];fill=Float64[]
    for z in axes.z_centres,t in axes.theta_centres
        push!(pressure,field_at(locations,fd.pressure_pa,t,z))
        push!(fill,field_at(locations,fd.fill_fraction,t,z))
    end
    areas=fv.grid.areas
    return (;pressure_relative_L2=sqrt(sum(areas.*(pressure-fv.pressure_pa).^2)/sum(areas.*pressure.^2)),
        pressure_relative_Linf=maximum(abs,pressure-fv.pressure_pa)/maximum(abs,pressure),
        fill_mean_absolute=sum(areas.*abs.(fill-fv.fill_fraction))/sum(areas))
end

function save_xi_field(out,model,r)
    a=r.cut.periodic;nt,nz=size(a.gt)
    open(joinpath(out,"fields-$model.tsv"),"w") do io
        println(io,"x_m\ty_m\tpressure_absolute_pa\tfill_fraction")
        writedlm(io,hcat(repeat(a.axes.theta_centres.*XI.radius_m,nz),
            repeat(a.axes.z_centres;inner=nt),r.pressure_pa,r.fill_fraction),'\t')
    end
    path=joinpath(out,"$model.vtk")
    write_vtk(path,r.grid)
    open(path,"a") do io
        for (name,field) in (("pressure_absolute_pa",r.pressure_pa),("fill_fraction",r.fill_fraction))
            println(io,"SCALARS $name double 1\nLOOKUP_TABLE default")
            foreach(value->println(io,value),field)
        end
    end
end

function run_xi(path)
    fresh_xi_output(path)
    controls=Dict{Int,Any}()
    open(joinpath(path,"results.tsv"),"w") do results
        open(joinpath(path,"profiles.tsv"),"w") do profiles
            firstrow=firstprofile=true
            for model in ("classical","production")
                previous=nothing
                for n in (64,128,256,512,1024)
                    started=time_ns();r=xi_solve(n,n÷4;model)
                    @test minimum(r.pressure_pa)>=0
                    @test 0<=minimum(r.fill_fraction)<=maximum(r.fill_fraction)<=1
                    @test r.stats.local_balance<1e-10 && r.stats.local_face_balance<1e-10
                    @test r.stats.global_balance<1e-10 && r.stats.complementarity<1e-10
                    change=isnothing(previous) ? (;pressure_relative_L2=NaN,fill_mean_absolute=NaN) : field_change(previous,r)
                    row=(;model,nt=n,nz=n÷4,r.stats...,pressure_mismatch_previous=change.pressure_relative_L2,
                        fill_mismatch_previous=change.fill_mean_absolute,elapsed_s=(time_ns()-started)/1e9)
                    row!(results,row;header=firstrow);firstrow=false
                    a=r.cut.periodic;x=[0.;a.axes.theta_centres.*XI.radius_m;XI.circumference_m]
                    pressure=xi_profile(r,x);fill=xi_profile(r,x;field=:fill_fraction)
                    for i in eachindex(x)
                        row!(profiles,(;model,nt=n,x_m=x[i],y_m=XI.length_m/2,pressure_absolute_pa=pressure[i],fill_fraction=fill[i]);header=firstprofile)
                        firstprofile=false
                    end
                    model=="classical" && n<=512 && (controls[n]=r)
                    n==1024 && save_xi_field(path,model,r)
                    println("$model $n × $(n÷4): $(r.stats)");flush(stdout);yield()
                    previous=r
                end
            end
        end
    end
    open(joinpath(path,"manufactured.tsv"),"w") do io
        firstrow=true
        for graded in (false,true),rpm in (-XI_RPM,XI_RPM)
            rows=[xi_manufactured(n,n÷4;graded,rpm) for n in (64,128,256)]
            for r in rows
                @test r.minimum_pressure>0 && r.fill_error==0
                @test r.local_balance<1e-10
                row!(io,r;header=firstrow);firstrow=false
            end
            @test rows[end].pressure_relative_L2<.002
            @test rows[end].pressure_relative_L2<rows[1].pressure_relative_L2/2
            println("MMS graded=$graded rpm=$rpm: $rows");flush(stdout);yield()
        end
    end
    open(joinpath(path,"independent.tsv"),"w") do io
        firstrow=true
        for n in (64,128,256,512)
            fd=xi_fd(n,n÷4);fv=controls[n]
            comparison=xi_fd_comparison(fv,fd)
            @test minimum(fd.pressure_pa)>=0 && 0<=minimum(fd.fill_fraction)<=maximum(fd.fill_fraction)<=1
            @test fd.residual<1e-10
            if n==512
                @test comparison.pressure_relative_L2<.01
                @test comparison.fill_mean_absolute<.01
            end
            row!(io,(;nt=n,nz=n÷4,comparison...,sweeps=fd.sweeps,residual=fd.residual);header=firstrow);firstrow=false
            serialize(joinpath(path,"fd-$n.jls"),(;theta=fd.a.theta,z_m=(fd.a.zeta.+fd.a.aspect).*XI.radius_m,
                pressure_pa=fd.pressure_pa,fill_fraction=fd.fill_fraction))
            println("FD $n × $(n÷4): $comparison");flush(stdout);yield()
        end
    end
    open(joinpath(path,"algebraic.tsv"),"w") do io
        baseline=xi_fd(64,16)
        firstrow=true
        for (label,r) in (("relaxation-0.8",xi_fd(64,16;relaxation=.8)),("stopping-1e-11",xi_fd(64,16;tolerance=1e-11)))
            pe=norm(r.pressure_pa-baseline.pressure_pa)/norm(r.pressure_pa)
            fe=maximum(abs,r.fill_fraction-baseline.fill_fraction)
            @test pe<1e-7 && fe<1e-7
            row!(io,(;control=label,pressure_relative_L2=pe,fill_Linf=fe);header=firstrow);firstrow=false
        end
    end
    return path
end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;description="Compare Xi2024 Fig.7 with a declared cosine-gap convention and four pressure edges.") do cfg
        @testset "Xi2024 finite-cylinder JFO" begin
            run_xi(cfg.output)
        end
    end)
end
