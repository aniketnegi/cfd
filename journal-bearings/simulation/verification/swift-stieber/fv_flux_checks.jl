include(joinpath(@__DIR__,"cases.jl"))
length(ARGS)==1 || error("usage: fv_flux_checks.jl FRESH_OUTPUT")
const OUTPUT=fresh_output(ARGS[1];suite="analytical-FV-end-flux")

@testset "FV end flux: independent linear and quadratic solutions" begin
    mu,c,Pstar,beta=.02752,50e-6,1e6,1.
    open(joinpath(OUTPUT,"fv-flux-checks.tsv"),"w") do io
        println(io,"R_m\tepsilon\tnt\tnz\tlinear_lower_relative_error\tlinear_upper_relative_error\tquadratic_lower_relative_error\tquadratic_upper_relative_error\tanalytic_sample_flux_ratio\texpected_sample_flux_ratio")
        for r in (.037,.05),ecc in (0.,.8),nt in (32,64,128)
            nz=nt÷2;length=2r*beta
            p=ClassicalCylinder(1000r,1000length,1000c,ecc,0.,0.)
            grid=cylinder_grid(p,nt,nz)
            integral_h3=2pi*c^3*(1+1.5ecc^2)
            a=assemble_full_film(grid;rpm=0.,viscosity_pa_s=mu,end_pressure_pa=(2Pstar,Pstar))
            solved=a.A\a.b
            expected_pressure=[2Pstar-Pstar*z/length for z in a.axes.z_centres for i in 1:nt]
            @test norm(solved-expected_pressure)/norm(expected_pressure)<1e-11
            flux=full_film_fluxes(a,solved)
            expected_linear=r*integral_h3/(12mu)*Pstar/length
            low,high=-sum(flux.qz[:,1]),sum(flux.qz[:,end])
            linear_errors=(abs(low/(-expected_linear)-1),abs(high/expected_linear-1))
            @test maximum(linear_errors)<1e-10
            @test abs(low+high)/expected_linear<1e-10

            a=assemble_full_film(grid;rpm=0.,viscosity_pa_s=mu,end_pressure_pa=(0.,0.))
            # Analytic manufactured injection = -K*p_zz, integrated by cell area.
            rhs=[(c*(1+ecc*cos(t)))^3/(12mu)*2Pstar/(r*beta)^2 *
                (r*2pi/nt)*(length/nz) for z in a.axes.z_centres for t in a.axes.theta_centres]
            solved=a.A\rhs
            sampled=[Pstar*(1-((z-length/2)/(r*beta))^2) for z in a.axes.z_centres for i in 1:nt]
            @test norm(solved-sampled.-Pstar/nz^2)/norm(sampled)<1e-11
            flux=full_film_fluxes(a,solved)
            expected_quadratic=integral_h3*Pstar/(6mu*beta)
            low,high=-sum(flux.qz[:,1]),sum(flux.qz[:,end])
            quadratic_errors=(abs(low/expected_quadratic-1),abs(high/expected_quadratic-1))
            @test maximum(quadratic_errors)<1e-10
            @test abs(low-high)/expected_quadratic<1e-10
            @test abs(sum(rhs)-2expected_quadratic)/expected_quadratic<1e-12
            sampled_flux=full_film_fluxes(a,sampled)
            sampled_ratio=(-sum(sampled_flux.qz[:,1])+sum(sampled_flux.qz[:,end]))/(2expected_quadratic)
            # A raw analytic cell-centre quadratic has only a first-order boundary secant.
            # The solved FV pressure offsets it correctly; do not confuse the two tests.
            @test sampled_ratio ≈ 1-1/(2nz) rtol=1e-12
            println(io,join((r,ecc,nt,nz,linear_errors...,quadratic_errors...,
                sampled_ratio,1-1/(2nz)),'\t'))
        end
    end
end
