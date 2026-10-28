include(joinpath(@__DIR__,"cases.jl"))
function independent_gap(p,theta,z;pose=pi)
    setprecision(256) do
        angle=BigFloat(p.semicone_angle_deg)*big(pi)/180
        radius=BigFloat(p.journal_midpoint_radius)/1000+
            (BigFloat(p.bearing_length)/2000-BigFloat(z))*tan(angle)
        c=BigFloat(p.radial_clearance)/1000;e=BigFloat(p.eccentricity)/1000
        t=BigFloat(theta)-BigFloat(pose)
        Float64((radius+c-e*cos(t)-sqrt(radius^2-e^2*sin(t)^2))*cos(angle))
    end
end
function quadrature_force(p,pressure,te,ze)
    x,w=gauss(8);v,u=gauss(2)
    nt,nz=length(te)-1,length(ze)-1
    force=zeros(3);k=tan(deg2rad(p.semicone_angle_deg))
    for j in 1:nz,i in 1:nt
        value=pressure[i+nt*(j-1)]
        for q in eachindex(v),s in eachindex(x)
            theta=(te[i]+te[i+1])/2+(te[i+1]-te[i])*x[s]/2
            z=(ze[j]+ze[j+1])/2+(ze[j+1]-ze[j])*v[q]/2
            radius=(p.journal_midpoint_radius+p.radial_clearance)*1e-3+
                (p.bearing_length*.0005-z)*k
            tangent_theta=[-radius*sin(theta),radius*cos(theta),0.]
            tangent_z=[-k*cos(theta),-k*sin(theta),1.]
            force .+= value*w[s]*u[q]*(te[i+1]-te[i])*(ze[j+1]-ze[j])/4 .* cross(tangent_theta,tangent_z)
        end
    end
    force
end
function run_controls(output)
    OUTPUT=fresh_output(output; suite="conical-controls")
    @testset "Independent cone metric, gap and face-coefficient checks" begin
        for angle in (0.,5.,10.,20.)
            p=cone(;angle,ecc=.8);nt,nz=32,16
            a=assemble_case(p,nt,nz;graded=true)
            ca,k=cosd(angle),tan(deg2rad(angle))
            rh(z)=.05005+(.05-z)*k
            rj(z)=.05+(.05-z)*k
            gt,adv,gz=zeros(nt,nz),zeros(nt,nz),zeros(nt,nz+1)
            gap=zeros(nt,nz)
            for j in 1:nz,i in 1:nt
                t,z=a.axes.theta_centres[i],a.axes.z_centres[j]
                gap[i,j]=independent_gap(p,t,z)
                face=a.axes.theta_edges[i+1]
                h=independent_gap(p,face,z)
                width=a.axes.z_widths[j]/ca
                distance=rh(z)*(a.axes.theta_widths[i]+a.axes.theta_widths[mod1(i+1,nt)])/2
                gt[i,j]=h^3*width/(12MU*distance)
                adv[i,j]=(2000*2pi/60)*rj(z)*h*width/2
            end
            for j in 1:nz+1,i in 1:nt
                z=a.axes.z_edges[j]
                h=independent_gap(p,a.axes.theta_centres[i],z)
                distance=j==1 ? a.axes.z_widths[1]/(2ca) : j==nz+1 ?
                    a.axes.z_widths[end]/(2ca) : (a.axes.z_widths[j-1]+a.axes.z_widths[j])/(2ca)
                gz[i,j]=h^3*rh(z)*a.axes.theta_widths[i]/(12MU*distance)
            end
            @test vec(gap)≈a.normal_gap_m rtol=1e-12
            @test gt≈a.gt rtol=3e-12
            @test gz≈a.gz rtol=3e-12
            @test adv≈a.adv rtol=1e-12
            @test sum(a.grid.areas)*1e-6≈2pi*.05005*.1/ca rtol=2e-14
            @test issymmetric(a.A)
            @test isposdef(Symmetric(Matrix(a.A)))
            @test all(x->x<=0,nonzeros(a.A-spdiagm(0=>diag(a.A))))
            values=[1e6*(1+.4cos(t)-.3sin(t))*z*(.1-z)/.1^2
                for t in a.axes.theta_centres,z in a.axes.z_centres]
            flux=full_film_fluxes(a,vec(values))
            qt=[adv[i,j]+gt[i,j]*(values[i,j]-values[mod1(i+1,nt),j]) for i in 1:nt,j in 1:nz]
            qz=[gz[i,j]*((j==1 ? 0. : values[i,j-1])-(j==nz+1 ? 0. : values[i,j])) for i in 1:nt,j in 1:nz+1]
            @test flux.qt≈qt rtol=3e-12
            @test flux.qz≈qz rtol=3e-12
            fq=quadrature_force(p,vec(values),a.axes.theta_edges,a.axes.z_edges)
            f=collect(housing_pressure_force(a.grid,vec(values)))
            @test norm(f-fq)/norm(fq)<1e-12
        end
    end
    @testset "Analytical continuum housing-pressure-force convergence" begin
        open(joinpath(OUTPUT,"force-analytic.tsv"),"w") do io
            first=true
            for angle in (0.,10.,20.)
                errors=Float64[]
                for nt in (32,64,128,256)
                    nz=nt÷2;p=cone(;angle)
                    g=grid_for(p,nt,nz;graded=true);axes=full_film_axes(g)
                    pressure=[1e6*(1+.4cos(t)-.3sin(t))*z*(.1-z)/.1^2 for t in axes.theta_centres,z in axes.z_centres]
                    f=collect(housing_pressure_force(g,vec(pressure)))
                    exact=1e6*.05005*.1/6 .* [.4pi,-.3pi,2pi*tand(angle)]
                    error=norm(f-exact)/norm(exact);push!(errors,error)
                    row!(io,(;angle_deg=angle,nt,nz,force_relative_L2=error);header=first);first=false
                end
                @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
                @test errors[end]<1e-4
            end
        end
    end
    @testset "Swift-Stieber cylinder reproduction and small-angle limit" begin
        fixture=TOML.parsefile(joinpath(@__DIR__,"reference.toml"))["production_cylinder_regression"]
        zero=solve_case(cone(;angle=0.,ecc=.8),128,64;mu=.02752,groove=true)
        check_result(zero)
        @test hypot(zero.force[1],zero.force[2])≈fixture["load_N"] rtol=1e-11
        @test maximum(zero.pressure)≈fixture["maximum_pressure_Pa"] rtol=1e-11
        @test sum(zero.ends)/((2000*2pi/60)*.05^2*50e-6)≈fixture["Q0"] rtol=1e-11
        @test zero.force[3]==0.
        open(joinpath(OUTPUT,"zero-angle.tsv"),"w") do io
            first=true
            for groove in (false,true)
                reference=solve_case(cone(;angle=0.),128,64;groove)
                errors=Float64[]
                for angle in (.1,.01,.001,.0001,.00001)
                    result=solve_case(cone(;angle),128,64;groove)
                    check_result(result)
                    error=norm(result.pressure-reference.pressure)/norm(reference.pressure)
                    push!(errors,error)
                    row!(io,(;groove,angle_deg=angle,pressure_relative_L2=error,force_z_N=result.force[3]);header=first);first=false
                end
                @test all(errors[2:end] .< 0.11 .* errors[1:end-1])
                @test errors[end]<1e-6
                @test reference.force[3]==0.
            end
        end
    end
    @testset "Obstacle symmetries, scaling and no-clipping controls" begin
        p=cone();nt,nz=64,32
        for groove in (false,true)
            base=solve_case(p,nt,nz;groove);check_result(base)
            twice_mu=solve_case(p,nt,nz;groove,mu=2MU)
            @test norm(twice_mu.pressure-2base.pressure)/norm(base.pressure)<1e-10
            @test norm(twice_mu.flux.qz-base.flux.qz)/norm(base.flux.qz)<1e-10
            @test norm(collect(twice_mu.force)-2collect(base.force))/norm(collect(base.force))<1e-10
            twice_rpm=solve_case(p,nt,nz;groove,rpm=4000.)
            @test norm(twice_rpm.pressure-2base.pressure)/norm(base.pressure)<1e-10
            @test norm(twice_rpm.flux.qz-2base.flux.qz)/norm(base.flux.qz)<1e-10
            negative=solve_case(p,nt,nz;groove,rpm=-2000.)
            reflection=[mod1(nt-i,nt) for i in 1:nt]
            expected=reshape(base.pressure,nt,nz)[reflection,:]
            @test norm(negative.pressure-vec(expected))/norm(base.pressure)<1e-10
            @test negative.force[1]≈base.force[1] rtol=1e-10
            @test negative.force[2]≈-base.force[2] rtol=1e-10
            @test negative.force[3]≈base.force[3] rtol=1e-10
            free=groove ? findall(i->i%nt!=0,1:nt*nz) : collect(1:nt*nz)
            signed=zeros(nt*nz);signed[free]=base.a.A[free,free]\base.a.b[free]
            @test norm(base.pressure-max.(signed,0.))/norm(base.pressure)>.01
        end
        base=solve_case(p,nt,nz)
        rotated=solve_case(p,nt,nz;pose=pi+pi/2)
        @test norm(rotated.pressure-vec(circshift(reshape(base.pressure,nt,nz),(nt÷4,0))))/norm(base.pressure)<1e-10
        @test rotated.force[1]≈-base.force[2] rtol=1e-10
        @test rotated.force[2]≈base.force[1] rtol=1e-10
        @test rotated.force[3]≈base.force[3] rtol=1e-10
        @test maximum(abs,solve_case(cone(;ecc=0.),nt,nz).pressure)==0.
        @test maximum(abs,solve_case(p,nt,nz;rpm=0.).pressure)==0.
    end
end
include(joinpath(@__DIR__, "../../../cli.jl"))
if abspath(PROGRAM_FILE) == @__FILE__
    exit(BearingCLI.study_main(; full=false, description="Run conical-controls verification.") do cfg
        run_controls(cfg.output)
    end)
end
