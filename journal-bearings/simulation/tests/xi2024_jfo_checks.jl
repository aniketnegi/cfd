using Test
include(joinpath(@__DIR__,"../verification/jfo/cylindrical/xi2024/cases.jl"))

function check_xi2024()
    @testset "Four-edge JFO boundary adapter" begin
        for graded in (false,true),rpm in (-XI_RPM,0.,XI_RPM)
            cut=xi_assembly(16,8;graded,rpm);a=cut.periodic;nt,nz=size(a.gt)
            P=[2e5*(1+sin(k)) for k in 1:nt*nz];F=[.7+.2cos(k) for k in 1:nt*nz]
            flux=cut_fluxes(cut,P,F)
            @test flux.balance≈cut.A*P+cut.B*(1 .-F)-cut.b rtol=1e-12
            @test sum(flux.balance)≈sum(flux.outward) atol=1e-18
            unequal=cut_angular_seam(a,(2e5,3e5))
            @test cut_fluxes(unequal,P,F).balance≈unequal.A*P+unequal.B*(1 .-F)-unequal.b rtol=1e-12
            @test norm(Matrix(cut.A-cut.A'))==0
            for j in (1,8)
                l,k=1+nt*(j-1),nt*j
                @test cut.A[l,k]==cut.A[k,l]==0
                @test cut.B[l,k]==cut.B[k,l]==0
                @test cut.A[l,l]≈a.A[l,l]+cut.low[j]-a.gt[end,j] rtol=1e-13
                @test cut.A[k,k]≈a.A[k,k]+cut.high[j]-a.gt[end,j] rtol=1e-13
                if !graded
                    @test cut.low[j]≈2a.gt[end,j]
                    @test cut.high[j]≈2a.gt[end,j]
                end
                @test flux.qt[1,j]≈cut.low[j]*(1e5-P[l])+a.adv[end,j]*(rpm>=0 ? 1 : F[l])
                @test flux.qt[end,j]≈cut.high[j]*(P[k]-1e5)+a.adv[end,j]*(rpm>=0 ? F[k] : 1)
            end
            @test Matrix(cut.A)[2:15,2:15]==Matrix(a.A)[2:15,2:15]
            @test_throws ArgumentError cut_angular_seam(a,(0.,1e5))
        end
        stationary=xi_solve(32,8;rpm=0.)
        @test stationary.pressure_pa≈fill(XI.boundary_absolute_pa,32*8) rtol=1e-12
        @test stationary.fill_fraction==ones(32*8)
        @test maximum(abs,stationary.flux.qt)<1e-20
        @test maximum(abs,stationary.flux.qz)<1e-20
        cut=xi_assembly(16,8)
        @test maximum(abs,cut.A*ones(128)-cut.periodic.A*ones(128))>0
        a=assemble_conical_jfo(cut.periodic.grid;rpm=XI_RPM,viscosity_pa_s=XI.mu_pa_s,
            end_pressure_pa=(1e8,1e8),cavitation_pressure_pa=0.,eccentricity_angle_rad=pi)
        flooded=cut_angular_seam(a,(1e8,1e8))
        state=solve_jfo_active_set(flooded.A,flooded.B,flooded.b,a.pressure_scale_pa)
        @test state.void==zeros(length(state.void))
        @test state.pressure≈flooded.A\flooded.b rtol=1e-12
    end
    @testset "Xi cylinder: scaling, reflection and actual-plane sampling" begin
        r=xi_solve(64,16);reverse=xi_solve(64,16;rpm=-XI_RPM)
        @test reshape(reverse.pressure_pa,64,16)≈Base.reverse(reshape(r.pressure_pa,64,16);dims=1) rtol=1e-10 atol=1e-7
        @test reshape(reverse.fill_fraction,64,16)≈Base.reverse(reshape(r.fill_fraction,64,16);dims=1) rtol=1e-10
        @test r.stats.local_balance<1e-10 && r.stats.global_balance<1e-10
        @test r.stats.complementarity<1e-10
        @test minimum(r.pressure_pa)>=0 && 0<minimum(r.fill_fraction)<1
        @test maximum(r.fill_fraction)<=1
        @test r.stats.gross_in_m3_s≈r.stats.gross_out_m3_s rtol=1e-10
        @test xi_profile(r,[0.,XI.circumference_m])≈[1e5,1e5]
        @test xi_profile(r,[0.,XI.circumference_m];field=:fill_fraction)==[1.,1.]
        a=r.cut.periodic
        @test minimum(a.normal_gap_m)>0
        @test normal_film_gap(a.grid.parameters,0.,0.)≈XI.clearance_m*(1-.8)
        @test a.eccentricity_angle_rad==Float64(pi)
        @test isapprox(XI_RPM,150.;rtol=1e-15)
        @test_throws ArgumentError xi_profile(r,[-.01])
        for graded in (false,true),rpm in (-XI_RPM,XI_RPM)
            coarse=xi_manufactured(32,8;graded,rpm)
            fine=xi_manufactured(64,16;graded,rpm)
            @test fine.pressure_relative_L2<coarse.pressure_relative_L2/2
            @test fine.minimum_pressure>0
            @test fine.fill_error==0
            @test fine.local_balance<1e-10
        end
    end
    @testset "Independent nodal cut-domain algebra and convergence" begin
        fd=xi_fd(32,8;tolerance=1e-11);a=fd.a
        fixed=collect(1:a.nt:length(a.b));free=setdiff(collect(eachindex(a.b)),fixed)
        rhs=a.b[free]-a.A[free,fixed]*fill(-a.pv0,length(fixed))
        state=solve_jfo_active_set(a.A[free,free],a.B[free,free],rhs,1.)
        @test norm(state.pressure-vec(fd.pressure_pa[:,2:end-1])[free]/XI_PS)/norm(state.pressure)<1e-7
        @test maximum(abs,state.void-vec(1 .-fd.fill_fraction[:,2:end-1])[free])<1e-7
        @test fd.pressure_pa[1,:]≈fill(XI.boundary_absolute_pa,size(fd.pressure_pa,2)) rtol=1e-14
        @test all(fd.fill_fraction[1,:].==1)
        @test 0<=minimum(fd.fill_fraction)<=maximum(fd.fill_fraction)<=1
        @test fd.residual<1e-11
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    check_xi2024()
end
