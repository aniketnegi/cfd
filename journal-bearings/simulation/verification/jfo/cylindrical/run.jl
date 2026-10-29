include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"manufactured.jl"))
include(joinpath(@__DIR__,"independent.jl"))
include(joinpath(ROOT,"cli.jl"))

function check_state(r)
    s=observables(r)
    @test minimum(r.pressure_above_cavitation_pa)>=0
    @test 0<=minimum(r.fill_fraction)<=maximum(r.fill_fraction)<=1
    @test s.local_balance<1e-10
    @test s.global_balance<1e-10
    @test s.complementarity<1e-10
    @test abs(s.net_end_m3_s)/max(s.gross_out_m3_s+s.gross_in_m3_s,1e-30)<1e-10
end

function fd_comparison(fv,fd)
    axes=full_film_axes(fv.grid);a=fd.a
    locations=(;theta_centres=a.theta,z_centres=(a.zeta.+a.aspect).*R)
    pressure=Float64[];fill=Float64[]
    for z in axes.z_centres,t in axes.theta_centres
        push!(pressure,PRESSURE_SCALE*(field_at(locations,fd.P,t,z)+a.pv0))
        push!(fill,field_at(locations,fd.F,t,z))
    end
    areas=fv.grid.areas
    return (;pressure_relative_L2=sqrt(sum(areas.*(pressure-fv.pressure_pa).^2)/sum(areas.*pressure.^2)),
        pressure_relative_Linf=maximum(abs,pressure-fv.pressure_pa)/maximum(abs,pressure),
        fill_mean_absolute=sum(areas.*abs.(fill-fv.fill_fraction))/sum(areas))
end

function independent_study(out,fv_controls;full=false)
    cases=[(.6,-1.),(.6,-.1)]
    full && append!(cases,[(.4,-.1),(.8,-.1)])
    open(joinpath(out,"independent.tsv"),"w") do io
        firstrow=true
        for (ecc,pv) in cases,nt in (64,128,256)
            a=CylindricalJFOReference.assemble(1.,ecc,pv,nt,nt÷2)
            fd=CylindricalJFOReference.solve(a)
            fs=CylindricalJFOReference.observables(fd)
            fv=fv_controls[(ecc,pv,nt)];vs=observables(fv)
            comparison=fd_comparison(fv,fd)
            @test minimum(fd.P)>=0
            @test 0<=minimum(fd.F)<=maximum(fd.F)<=1
            @test fd.residual<1e-10
            row=(;epsilon=ecc,pv0=pv,nt,nz=nt÷2,fd_P0=fs.P0,fd_beta_deg=fs.beta_deg,fd_Q0=fs.Q0,
                fd_net_Q0=fs.net_Q0,fv_P0=vs.P0,fv_beta_deg=vs.beta_deg,fv_Q0=vs.Q0,
                comparison...,sweeps=fd.sweeps,residual=fd.residual,relaxation=fd.relaxation,tolerance=fd.tolerance)
            row!(io,row;header=firstrow);firstrow=false
            serialize(joinpath(out,"fd-$(ecc)-$(pv)-$(nt).jls"),(;theta=a.theta,zeta=a.zeta,P=fd.P,F=fd.F,stats=fs))
            println("Independent FD epsilon=$ecc pv0=$pv $nt × $(nt÷2): $fs");yield()
        end
    end
    open(joinpath(out,"algebraic.tsv"),"w") do io
        firstrow=true
        for (ecc,pv) in cases
            a=CylindricalJFOReference.assemble(1.,ecc,pv,64,32)
            baseline=CylindricalJFOReference.solve(a)
            for (name,computed) in (("sor-relaxation-0.8",CylindricalJFOReference.solve(a;relaxation=.8)),
                    ("sor-tighter",CylindricalJFOReference.solve(a;tolerance=1e-11)))
                p=norm(computed.P-baseline.P)/norm(computed.P)
                f=maximum(abs,computed.F-baseline.F)
                @test p<1e-7 && f<1e-7
                row!(io,(;epsilon=ecc,pv0=pv,control=name,pressure_relative_L2=p,fill_Linf=f);header=firstrow);firstrow=false
            end
            # Same independent FD matrices, solved by the production policy kernel:
            # this tests algebra only, not a second spatial method.
            state=solve_jfo_active_set(a.A,a.B,a.b,1.)
            p=norm(state.pressure-vec(baseline.P[:,2:end-1]))/norm(state.pressure)
            f=maximum(abs,(1 .-state.void)-vec(baseline.F[:,2:end-1]))
            @test p<1e-7 && f<1e-7
            row!(io,(;epsilon=ecc,pv0=pv,control="policy-same-fd-matrices",pressure_relative_L2=p,fill_Linf=f))
        end
    end
end

function manufactured_study(out)
    open(joinpath(out,"manufactured.tsv"),"w") do io
        firstrow=true
        for graded in (false,true),rpm in (-RPM,RPM)
            rows=[manufactured_case(n,n÷2;rpm,graded) for n in (64,128,256)]
            for r in rows
                @test r.minimum_pressure>=0
                @test 0<=r.minimum_fill<=r.maximum_fill<=1
                @test r.local_balance<1e-10 && r.complementarity<1e-10
                row!(io,r;header=firstrow);firstrow=false
                println("Manufactured $r");yield()
            end
            @test rows[end].pressure_relative_L2<.01
            @test rows[end].fill_mean_absolute<.005
            @test rows[end].pressure_relative_L2<rows[1].pressure_relative_L2/2
            @test rows[end].fill_mean_absolute<rows[1].fill_mean_absolute/2
        end
    end
end

function run_cylindrical_jfo(path;full=false)
    out=fresh_output(path)
    metadata=TOML.parsefile(joinpath(out,"inputs.toml"))
    metadata["full"]=full
    metadata["n_theta"]=[64,128,256,512]
    metadata["independent_n_theta"]=[64,128,256]
    open(joinpath(out,"inputs.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    points=reference_rows()
    full || (points=filter(r->r.L_D==1. && r.epsilon==.6,points))
    fv_controls=Dict{Tuple{Float64,Float64,Int},Any}()
    open(joinpath(out,"results.tsv"),"w") do io
        open(joinpath(out,"profiles.tsv"),"w") do profiles
            firstrow=firstprofile=true
            for ref in points,model in ("classical","production")
                previous=nothing
                for nt in (64,128,256,512)
                    started=time_ns()
                    r=cylinder_solution(ref.L_D,ref.epsilon,ref.pv0,nt;model)
                    check_state(r);s=observables(r)
                    change=isnothing(previous) ? (;pressure_relative_L2=NaN,fill_mean_absolute=NaN) : field_change(previous,r)
                    row=(;model,L_D=ref.L_D,epsilon=ref.epsilon,pv0=ref.pv0,nt,nz=r.grid.n_axial,s...,
                        pressure_mismatch_previous=change.pressure_relative_L2,fill_mismatch_previous=change.fill_mean_absolute,
                        reference_P0=ref.P0,reference_beta_deg=ref.beta_deg,reference_Q0=ref.Q0,
                        elapsed_s=(time_ns()-started)/1e9)
                    row!(io,row;header=firstrow);firstrow=false
                    if ref.L_D==1. && model=="classical" && nt<=256 &&
                            (ref.epsilon==.6 || (ref.pv0==-.1 && ref.epsilon in (.4,.8)))
                        fv_controls[(ref.epsilon,ref.pv0,nt)]=r
                    end
                    println("$model LD=$(ref.L_D) epsilon=$(ref.epsilon) pv0=$(ref.pv0) $nt × $(r.grid.n_axial): P0=$(s.P0) beta=$(s.beta_deg) Q0=$(s.Q0)")
                    if nt==512
                        axes=full_film_axes(r.grid)
                        pressure=plane_profile(r.grid,r.pressure_pa);fill=plane_profile(r.grid,r.fill_fraction)
                        phi=mod.(axes.theta_centres.+pi,2pi).-pi
                        for i in sortperm(phi)
                            profile=(;model,L_D=ref.L_D,epsilon=ref.epsilon,pv0=ref.pv0,nt,phi_rad=phi[i],
                                z_m=ref.L_D*R,p0=pressure[i]/PRESSURE_SCALE,fill_fraction=fill[i])
                            row!(profiles,profile;header=firstprofile);firstprofile=false
                        end
                        if ref.L_D==1. && ref.epsilon==.6
                            name="$(model)-$(ref.pv0)"
                            serialize(joinpath(out,"$name.jls"),(;pressure_pa=r.pressure_pa,fill_fraction=r.fill_fraction,
                                theta=axes.theta_centres,z_m=axes.z_centres,stats=s))
                            write_conical_jfo_vtk(joinpath(out,"$name.vtk"),r)
                            if model=="classical"
                                open(joinpath(out,"field-$(ref.pv0).tsv"),"w") do field
                                    println(field,"theta_rad\tz_m\tp0\tfill_fraction")
                                    writedlm(field,hcat(repeat(axes.theta_centres,r.grid.n_axial),
                                        repeat(axes.z_centres;inner=nt),r.pressure_pa/PRESSURE_SCALE,r.fill_fraction),'\t')
                                end
                            end
                        end
                    end
                    previous=r;yield()
                end
            end
        end
    end
    manufactured_study(out)
    independent_study(out,fv_controls;full)
    return out
end

if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,
        description="Verify steady finite-cylinder JFO; --full compares all 24 historical rows.") do cfg
        @testset "Finite-cylinder JFO verification" begin
            run_cylindrical_jfo(cfg.output;full=cfg.full)
        end
    end)
end
