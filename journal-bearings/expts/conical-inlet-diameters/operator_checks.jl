# Independent manufactured fields on each saved diameter-specific coarse map.
include(joinpath(@__DIR__,"run.jl"))

function diameter_operator_checks(out)
    path=joinpath(out,"operator.tsv")
    ispath(path) && error("refusing to overwrite diameter operator evidence: $path")
    metadata=TOML.parsefile(joinpath(out,"geometry.toml"))
    for (file,digest) in metadata["source_sha256"]
        bytes2hex(sha256(read(joinpath(ROOT,file))))==digest || error("study source changed: $file")
    end
    for (file,digest) in metadata["coarse_sha256"]
        bytes2hex(sha256(read(joinpath(out,file))))==digest || error("coarse map changed: $file")
    end
    passed=false
    try
        open(path,"w") do io
            println(io,"diameter_mm\tkind\trpm\tlevel\tpressure_error\tfill_error\tlocal_residual\twall_s")
            @testset "Manufactured fields across seven inlet diameters" begin
                for d in DIAMETERS
                    @testset "D=$d mm" begin
                        coarse=load_coarse(out,d)
                        @testset "smooth pressure" begin
                            errors=Float64[]
                            for level in (1,2,4)
                                t=@elapsed begin
                                    m=C.subdivide_layout(coarse,level);kx=2pi/m.period;ky=pi/(m.L/m.ca)
                                    exact(x,y)=1+.2*cos(kx*x)*sin(ky*y)
                                    a=FieldReference.fitted_operator(m,FieldReference.pose;unit_mobility=true,
                                        boundary_pressure=(x,y,tag)->exact(x,y))
                                    P=[exact(m.centre[k,:]...) for k in eachindex(m.quads)]
                                    source=[begin
                                        x,y=m.centre[k,:];rh=m.R+(.05-y*m.ca)*tand(10.)
                                        .2*cos(kx*x)*((m.R/rh*kx^2+rh/m.R*ky^2)*sin(ky*y)+sind(10.)/m.R*ky*cos(ky*y))
                                    end for k in eachindex(m.quads)].*m.areas
                                    result=a.A\(a.b+source)
                                    err=sqrt(sum((result-P).^2 .*m.areas)/sum((P.-1).^2 .*m.areas));push!(errors,err)
                                    flux=a.pressure_flux*result+a.constant
                                    residual=maximum(abs,a.divergence*flux-source)/maximum(abs,flux)
                                    @test minimum(result)>0
                                    @test residual<1e-9
                                end
                                println(io,join((d,"pressure",0,level,err,NaN,residual,t),'\t'));flush(io)
                            end
                            @test errors[1]/errors[2]>2.5
                            @test errors[2]/errors[3]>2.5
                            @test errors[3]<.002
                        end
                        manufactured_pose=FieldReference.floating_bush_pose(FieldReference.reference,(-30e-6,0.,0.))
                        for rpm in (2000.,-2000.)
                            @testset "curved cavity rpm=$rpm" begin
                                pe=Float64[];fe=Float64[]
                                for level in (1,2)
                                    t=@elapsed begin
                                        m=C.subdivide_layout(coarse,level)
                                        a=FieldReference.fitted_operator(m,manufactured_pose;rpm,
                                            boundary_pressure=(x,y,tag)->FieldReference.fitted_mms_point(x,y,m;rpm).P,
                                            incoming_fill=(x,y,tag)->FieldReference.fitted_mms_point(x,y,m;rpm).F)
                                        exact=[FieldReference.fitted_mms_point(m.centre[k,:]...,m;rpm) for k in eachindex(m.quads)]
                                        source=[f.source for f in exact].*m.areas
                                        state=FieldReference.solve_jfo_active_set(a.A,a.B,a.b+source,1e8;max_iterations=250)
                                        P=[f.P for f in exact];F=[f.F for f in exact]
                                        ep=sqrt(sum((state.pressure-P).^2 .*m.areas)/sum(P.^2 .*m.areas))
                                        ef=sum(abs.(1 .-state.void-F).*m.areas)/sum(m.areas)
                                        push!(pe,ep);push!(fe,ef)
                                        flux=a.pressure_flux*state.pressure+a.constant+a.adv-a.transport*state.void
                                        residual=maximum(abs,a.divergence*flux-source)/maximum(abs,a.adv)
                                        @test minimum(state.pressure)>=-1e-5
                                        @test minimum(1 .-state.void)>=0 && maximum(1 .-state.void)<=1+1e-10
                                        @test residual<1e-8
                                    end
                                    println(io,join((d,"JFO",rpm,level,ep,ef,residual,t),'\t'));flush(io)
                                end
                                @test pe[2]<pe[1]/1.5
                                @test fe[2]<fe[1]/1.5
                                @test pe[2]<.01 && fe[2]<.005
                            end
                        end
                    end
                    println("Finished manufactured-field checks for D=$d mm");flush(stdout)
                end
            end
        end
        passed=true
    finally
        open(joinpath(out,"operator.toml"),"w") do io
            TOML.print(io,Dict("checks_passed"=>passed,"meaning"=>"Smooth manufactured fields and conservation, NOT physical load predictions or an independent discontinuous-front benchmark",
                "source_geometry_manifest_sha256"=>bytes2hex(sha256(read(joinpath(out,"geometry.toml")))),
                "operator_table_sha256"=>bytes2hex(sha256(read(path))));sorted=true)
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("provide the completed geometry evidence directory")
    diameter_operator_checks(ARGS[1])
end
