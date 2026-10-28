include(joinpath(@__DIR__,"cases.jl"))
const PREF=1e7
const FRONT=.37
function manufactured(theta,z,angle,rpm)
    L,c,e,mu=.1,50e-6,30e-6,MU
    k,ca=tand(angle),cosd(angle)
    rj=.05+(L/2-z)*k;r=rj+c
    st,ct=sincos(theta-pi)
    D=sqrt(rj^2-e^2*st^2)
    h=ca*(c-e*ct+e^2*st^2/(rj+D))
    ht=ca*(e*st+e^2*st*ct/D);hz=ca*k*(rj/D-1)
    K,Kt,Kz=h^3/(12mu),h^2*ht/(4mu),h^2*hz/(4mu)
    t=mod(theta-FRONT,2pi)
    wet=0<t<pi
    s,q=sincos(t)
    g=wet ? s^4 : 0.;gt=wet ? 4s^3*q : 0.;gtt=wet ? 12s^2*q^2-4s^4 : 0.
    f,df=z*(L-z)/L^2,(L-2z)/L^2
    pressure=PREF*g*f
    pt,ptt,pz,pzz=PREF*gt*f,PREF*gtt*f,PREF*g*df,-2PREF*g/L^2
    diffusion=(Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
    wedge=rpm*2pi/60*rj*ht/(2ca)
    reaction_density=wet ? 0. : PREF*(c*ca)^3/(12mu*.05005*ca)*(1+.25cos(theta))
    return (;pressure,source=wedge-diffusion-reaction_density,h,reaction_density)
end
function exact_ends(angle,n)
    x,w=gauss(n)
    theta=FRONT .+ (x.+1).*pi/2
    result=zeros(2)
    for (j,z) in enumerate((0.,.1))
        r=.05005+(.05-z)*tand(angle)
        for i in eachindex(x)
            ref=manufactured(theta[i],z,angle,2000.)
            g=sin(theta[i]-FRONT)^4
            result[j]+=w[i]*pi/2*ref.h^3/(12MU)*r*cosd(angle)*PREF*g/.1
        end
    end
    return result
end
function run_mms(output)
    OUTPUT=fresh_output(output; suite="conical-manufactured")
    ispath(joinpath(OUTPUT,"mms.tsv")) && error("refusing overwrite")
    families=vec([(angle,rpm,false) for angle in (0.,10.,20.),rpm in (-2000.,2000.)])
    append!(families,[(10.,rpm,true) for rpm in (-2000.,2000.)])
    @testset "Conical pressure-obstacle manufactured fields and end flows" begin
        open(joinpath(OUTPUT,"mms.tsv"),"w") do io
            first=true
            for (angle,rpm,graded) in families
                errors=Float64[];flow_errors=Float64[];force_errors=Float64[]
                qe=exact_ends(angle,160)
                @test norm(qe-exact_ends(angle,80))/norm(qe)<1e-12
                for nt in (32,64,128,256)
                    nz=nt÷2
                    a=assemble_case(cone(;angle),nt,nz;rpm,graded)
                    exact,source=zeros(nt,nz),zeros(nt,nz)
                    for j in 1:nz,i in 1:nt
                        ref=manufactured(a.axes.theta_centres[i],a.axes.z_centres[j],angle,rpm)
                        exact[i,j]=ref.pressure
                        source[i,j]=ref.source*a.axes.theta_widths[i]*a.axes.z_widths[j]
                    end
                    result=obstacle(a;rpm,source=vec(source))
                    check_result(result)
                    weights=a.grid.areas
                    error=sqrt(sum((result.pressure-vec(exact)).^2 .* weights)/sum(vec(exact).^2 .* weights))
                    flow_error=norm(collect(result.ends)-qe)/norm(qe)
                    radial=PREF*.05005*.1/6
                    fe=radial.*[-16sin(FRONT)/15,16cos(FRONT)/15,tand(angle)*3pi/8]
                    force_error=norm(collect(result.force)-fe)/norm(fe)
                    push!(errors,error);push!(flow_errors,flow_error);push!(force_errors,force_error)
                    row!(io,(;angle_deg=angle,rpm,graded,nt,nz,pressure_relative_L2=error,
                        end_flow_relative_L2=flow_error,force_relative_L2=force_error,
                        exact_lower_flow=qe[1],exact_upper_flow=qe[2],result.stats...);header=first)
                    first=false
                    @printf("MMS angle=%g rpm=%g graded=%s nt=%d pressure=%.4g flow=%.4g force=%.4g\n",
                        angle,rpm,graded,nt,error,flow_error,force_error);flush(stdout)
                    if nt==256 && angle==10. && rpm==2000.
                        save_field(OUTPUT,"mms-graded$(graded)",result;extra=(;exact=vec(exact)))
                    end
                end
                @test all(errors[2:end] .< 0.35 .* errors[1:end-1])
                @test all(flow_errors[2:end] .< 0.4 .* flow_errors[1:end-1])
                @test all(force_errors[2:end] .< 0.4 .* force_errors[1:end-1])
                @test errors[end]<.002
                @test flow_errors[end]<.002
                @test force_errors[end]<.002
            end
        end
    end
end
include(joinpath(@__DIR__, "../../../cli.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    exit(BearingCLI.study_main(; full=false, description="Run conical-manufactured verification.") do cfg
        run_mms(cfg.output)
    end)
end
