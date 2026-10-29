include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(SOURCE,"mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"fitted_reference.jl"))

function fitted_obstacle_checks(A,b,P,pscale,qscale)
    reaction=A*P-b
    @test minimum(P)>=0
    @test minimum(reaction)/qscale > -1e-10
    @test maximum(abs,(P/pscale).*(reaction/qscale)) < 1e-10
    return reaction
end

function fitted_verification_mesh(angle,nt;observer=nothing)
    @printf("Fitted mesh angle=%g counts=%dx%d\n",angle,nt,nt÷2);flush(stdout)
    return surface_feed_grid(cone(;angle),nt,nt÷2;feed_diameter=5.,feed_z=50.,observer)
end

function fitted_verification(output;full=false,observer=nothing)
    fresh_output(output;suite="fitted-pressure-obstacle")
    inputs=TOML.parsefile(joinpath(output,"inputs.toml"))
    merge!(inputs,Dict("fixture_feed_pressure_pa"=>3e5,"fixture_viscosity_pa_s"=>MU,
        "eccentricity_angle_rad"=>.37,
        "feed_diameter_mm"=>5.,"feed_z_mm"=>50.,"mms_front_rad"=>[.37,pi/2],
        "grids"=>full ? [128,256,512] : [128,256]))
    open(joinpath(output,"inputs.toml"),"w") do io
        TOML.print(io,inputs;sorted=true)
    end
    grids=full ? (128,256,512) : (128,256)
    meshes=Dict((angle,nt)=>fitted_verification_mesh(angle,nt;observer) for angle in (0.,10.),nt in grids)
    @testset "Fitted manufactured pressure obstacle, force and boundary flux" begin
        open(joinpath(output,"fitted-mms.tsv"),"w") do io
            first=true
            for angle in (0.,10.),front in (.37,pi/2),rpm in (-2000.,2000.)
                errors,flow_errors,force_errors=Float64[],Float64[],Float64[]
                qe=fitted_exact_ends(angle,front)
                @test norm(qe-fitted_exact_ends(angle,front;n=80))/norm(qe)<1e-10
                for nt in grids
                    grid=meshes[(angle,nt)];s=grid.surface;p=grid.parameters
                    boundary_pressure(x,y,tag)=fitted_mms(x/s.R,y*s.ca,angle,rpm;front).pressure
                    faces=fitted_face_operator(s,p;rpm,viscosity_pa_s=MU,
                        eccentricity_angle_rad=.37,boundary_pressure)
                    refs=[fitted_mms(s.centre[i,1]/s.R,s.centre[i,2]*s.ca,angle,rpm;front) for i in axes(s.centre,1)]
                    exact=[ref.pressure for ref in refs]
                    source=[refs[i].source*s.areas[i] for i in eachindex(refs)]
                    A=faces.divergence*faces.pressure_flux
                    b=-(faces.divergence*faces.constant)+source
                    pscale=1e7;qscale=max(maximum(abs,faces.advective),maximum(abs,b))
                    state=pressure_obstacle(A,b,pscale,qscale)
                    P=state.pressure;reaction=fitted_obstacle_checks(A,b,P,pscale,qscale)
                    flow=faces.pressure_flux*P+faces.constant
                    boundary=boundary_outflows(faces,flow)
                    @test norm(faces.divergence*flow-source-reaction,Inf)/qscale<1e-10
                    @test abs(sum(values(boundary))-sum(source)-sum(reaction))/
                        max(sum(abs,source)+sum(abs,reaction),qscale)<1e-10
                    qfeed=fitted_exact_feed(s,angle,rpm,front)
                    qexact=[qe;qfeed]
                    @test abs(qfeed-fitted_exact_feed(s,angle,rpm,front;n=12))/norm(qexact)<1e-8
                    fexact=fitted_exact_force(s,angle,front)
                    @test norm(fexact-fitted_exact_force(s,angle,front;n=10))/norm(fexact)<1e-10
                    error=sqrt(sum((P-exact).^2 .* grid.areas)/sum(exact.^2 .* grid.areas))
                    flow_error=norm(collect(values(boundary))-qexact)/norm(qexact)
                    force_error=norm(collect(housing_pressure_force(grid,P))-fexact)/norm(fexact)
                    push!(errors,error);push!(flow_errors,flow_error);push!(force_errors,force_error)
                    row!(io,(;angle_deg=angle,front_rad=front,rpm,nt,nz=nt÷2,
                        pressure_relative_L2=error,boundary_flow_relative_L2=flow_error,
                        force_relative_L2=force_error,lower_flow=boundary.large,upper_flow=boundary.small,
                        feed_flow=boundary.feed,exact_lower_flow=qe[1],exact_upper_flow=qe[2],exact_feed_flow=qfeed,
                        iterations=state.iterations,complementarity=maximum(abs,(P/pscale).*(reaction/qscale)));
                        header=first)
                    first=false
                    @printf("Fitted MMS angle=%g front=%g rpm=%g nt=%d P=%.6g Q=%.6g F=%.6g\n",
                        angle,front,rpm,nt,error,flow_error,force_error);flush(stdout)
                    if nt==last(grids) && angle==10. && front==pi/2 && rpm>0
                        serialize(joinpath(output,"manufactured-field.jls"),(;theta=s.centre[:,1]/s.R,
                            z=s.centre[:,2]*s.ca,P,exact,source,reaction,normal_gap_m=[ref.h for ref in refs],flow))
                    end
                    GC.gc();yield()
                end
                @test all(errors[2:end] .< .65 .* errors[1:end-1])
                @test all(flow_errors[2:end] .< .75 .* flow_errors[1:end-1])
                @test all(force_errors[2:end] .< .65 .* force_errors[1:end-1])
                @test errors[end]<.003
                @test flow_errors[end]<.003
                @test force_errors[end]<.003
            end
        end
    end
    @testset "Fitted rotating physical-equation refinement" begin
        open(joinpath(output,"fitted-refinement.tsv"),"w") do io
            first=true
            for angle in (0.,10.),ecc in (.2,.6,.8),rpm in (-2000.,2000.)
                previous=nothing
                for nt in grids
                    grid=(;meshes[(angle,nt)]...,parameters=cone(;angle,ecc))
                    result=solve_swift_stieber(grid;rpm,viscosity_pa_s=MU,feed_pressure_pa=3e5,
                        eccentricity_angle_rad=.37)
                    s=result.diagnostics
                    @test s.complementarity_relative<1e-10
                    @test s.positive_pressure_relative_residual<1e-10
                    @test s.defect_accounting_relative<1e-10
                    @test any(iszero,result.pressure_above_cavitation_pa)
                    force=collect(result.housing_pressure_force_n)
                    flow=collect(values(result.boundary_outflow_m3_s))
                    peak=maximum(result.pressure_pa)
                    dp=isnothing(previous) ? NaN : abs(peak/previous.peak-1)
                    df=isnothing(previous) ? NaN : norm(force-previous.force)/norm(force)
                    dq=isnothing(previous) ? NaN : norm(flow-previous.flow)/norm(flow)
                    row!(io,(;angle_deg=angle,epsilon=ecc,rpm,nt,nz=nt÷2,
                        maximum_pressure_pa=peak,pressure_peak_relative_change=dp,
                        force_relative_change=df,boundary_flow_relative_change=dq,
                        force_x_N=force[1],force_y_N=force[2],force_z_N=force[3],
                        lower_outflow_m3_s=flow[1],upper_outflow_m3_s=flow[2],feed_outflow_m3_s=flow[3],
                        s...);header=first)
                    first=false
                    if nt==last(grids)
                        @test max(dp,df,dq)<.02
                    end
                    previous=(;peak,force,flow)
                    @printf("Fitted refinement angle=%g ecc=%g rpm=%g nt=%d peak=%.8g Pa\n",
                        angle,ecc,rpm,nt,peak);flush(stdout)
                    GC.gc();yield()
                end
            end
        end
    end
    @testset "Independent algebraic projected-SOR comparison" begin
        grid=meshes[(10.,128)]
        a=assemble_full_film(grid;rpm=2000.,viscosity_pa_s=MU,feed_pressure_pa=3e5,
            eccentricity_angle_rad=.37)
        result=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=MU,feed_pressure_pa=3e5,
            eccentricity_angle_rad=.37)
        scale=result.diagnostics.pressure_scale_pa
        reference=fitted_projected_sor(a.A,a.b/scale)
        difference=norm(reference.pressure-result.pressure_pa/scale)/norm(reference.pressure)
        @test difference<1e-7
        @test reference.residual<1e-12
        control=fitted_projected_sor(a.A,a.b/scale;relaxation=1.,tolerance=1e-11)
        @test control.residual<1e-11
        @test norm(reference.pressure-control.pressure)/norm(reference.pressure)<1e-7
        row=(;pressure_relative_L2=difference,sweeps=reference.sweeps,natural_residual=reference.residual,
            control_sweeps=control.sweeps,control_natural_residual=control.residual)
        open(joinpath(output,"fitted-algebraic.tsv"),"w") do io
            row!(io,row;header=true)
        end
        println("Projected-SOR comparison: ",row);flush(stdout)
    end
end
include(joinpath(@__DIR__,"../../../cli.jl"))
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,description="Verify fitted pressure-only cavitation and its boundary fluxes.") do cfg
        observer=stderr isa Base.TTY && !cfg.quiet && !cfg.verbose ? _->yield() : nothing
        fitted_verification(cfg.output;full=cfg.full,observer)
    end)
end
