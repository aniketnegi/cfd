using Test
include(joinpath(@__DIR__,"../reynolds/src/conical_equilibrium.jl"))
reference=BearingParameters(journal_midpoint_radius=50.0,radial_clearance=0.05,
    bearing_length=100.0,semicone_angle_deg=10.0,eccentricity_ratio=0.0)

@testset "Floating-bush geometry and frame" begin
    for d in ((0.,0.,0.),(12e-6,-8e-6,1e-4),(-5e-6,10e-6,-1e-4))
        pose=floating_bush_pose(reference,d); p=pose.parameters
        @test p.housing_midpoint_radius ≈ reference.housing_midpoint_radius atol=1e-13
        @test pose.radial_clearance_m ≈ 50e-6+d[3]*tand(10)
        @test pose.required_large_end_extension_m==max(-d[3],0.)
        @test pose.required_small_end_extension_m==max(d[3],0.)
        for theta in (0.,0.9,2.1,4.2),z in (0.01,0.08)
            rh=0.05005+(0.05-z)*tand(10)
            rj=0.05+(0.05-z-d[3])*tand(10)
            parallel=-(d[1]*cos(theta)+d[2]*sin(theta))
            perpendicular=d[1]*sin(theta)-d[2]*cos(theta)
            intersection=parallel+sqrt(rj^2-perpendicular^2)
            expected=(rh-intersection)*cosd(10)
            actual=conical_normal_gap(p,theta-pose.eccentricity_angle_rad,z)
            @test actual ≈ expected rtol=1e-10
            @test hypot(d[1]+intersection*cos(theta),d[2]+intersection*sin(theta)) ≈ rj rtol=1e-13
        end
    end
    @test_throws ArgumentError floating_bush_pose(reference,(0.,0.,-1e-3))
    @test_throws ArgumentError floating_bush_pose(reference,(50e-6,0.,0.))
    @test_throws ArgumentError floating_bush_pose(reference,(NaN,0.,0.))
end

@testset "Fixed housing mesh, rotated wedge and pressure force" begin
    nt,nz=32,16
    p=floating_bush_pose(reference,(-30e-6,0.,0.)).parameters
    old=assemble_conical_full_film(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277)
    zero=assemble_conical_full_film(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277,eccentricity_angle_rad=0.)
    @test old.A==zero.A && old.b==zero.b
    shift=5; angle=shift*2pi/nt
    a=assemble_conical_full_film(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277,eccentricity_angle_rad=angle)
    @test reshape(a.normal_gap_m,nt,nz) ≈ circshift(reshape(old.normal_gap_m,nt,nz),(shift,0)) rtol=1e-13
    @test a.grid.nodes==old.grid.nodes
    @test a.grid.radial_gaps ≈ 1e3.*a.normal_gap_m/cosd(10)
    base=solve_conical_jfo(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(5e5,0.),cavitation_pressure_pa=0.)
    rotated=solve_conical_jfo(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(5e5,0.),cavitation_pressure_pa=0.,eccentricity_angle_rad=angle)
    @test reshape(rotated.pressure_pa,nt,nz) ≈ circshift(reshape(base.pressure_pa,nt,nz),(shift,0)) rtol=1e-10
    @test reshape(rotated.fill_fraction,nt,nz) ≈ circshift(reshape(base.fill_fraction,nt,nz),(shift,0)) atol=1e-10
    f=conical_housing_pressure_force(p,base.pressure_pa,nt,nz)
    fr=conical_housing_pressure_force(p,rotated.pressure_pa,nt,nz)
    @test collect(fr) ≈ [cos(angle)*f[1]-sin(angle)*f[2],sin(angle)*f[1]+cos(angle)*f[2],f[3]] rtol=1e-10
    uniform=conical_housing_pressure_force(reference,fill(5e5,nt*nz),nt,nz)
    @test hypot(uniform[1],uniform[2])<1e-10
    @test uniform[3] ≈ 5e5*pi*((.05005+.05*tand(10))^2-(.05005-.05*tand(10))^2) rtol=1e-13
    @test_throws ArgumentError assemble_conical_full_film(p,nt,nz;rpm=2000.,viscosity_pa_s=.0277,eccentricity_angle_rad=NaN)
end

@testset "Aperture traction and equilibrium iteration" begin
    nt,nz=128,64; mask=floating_bush_feed(reference,nt,nz)
    state=floating_bush_state(reference,(-1e-6,0.3e-6,0.),nt,nz;feed_cells=mask)
    aperture=zeros(nt*nz);aperture[mask].=5e5
    omitted=conical_housing_pressure_force(reference,aperture,nt,nz)
    @test collect(state.omitted_aperture_force_n) ≈ collect(omitted) atol=1e-10
    @test abs(omitted[1])<1e-10 && omitted[2]>0 && omitted[3]>0
    @test state.bush_force_n[3]>=tand(10)*hypot(state.bush_force_n[1:2]...)
    @test state.film.diagnostics.relative_global_imbalance<1e-10
    fixed=equilibrate_floating_bush(reference,nt,nz;feed_cells=mask,fixed_axial=true,fd_step=1e-5)
    @test fixed.converged
    @test maximum(abs,collect(fixed.residual_n)[1:2])<0.05
    @test fixed.state.pose.displacement_m[3]==0
    @test fixed.required_axial_support_n ≈ 200-fixed.state.bush_force_n[3]
    # An inverse solve checks ONLY root finding, not independent physical validity.
    known=(-50e-6,25e-6,0.0008)
    forward=floating_bush_state(reference,known,nt,nz;feed_cells=mask)
    inverse=equilibrate_floating_bush(reference,nt,nz;feed_cells=mask,
        external_force_n=Tuple(-collect(forward.bush_force_n)),initial_displacement_m=collect(known).*[1.03,.98,1.02],fd_step=1e-5)
    @test inverse.converged
    @test norm(collect(inverse.state.pose.displacement_m).-collect(known))<1e-7
    failed=equilibrate_floating_bush(reference,nt,nz;feed_cells=mask,max_iterations=1)
    @test !failed.converged
    @test failed.reason=="iteration limit"
    @test_throws ArgumentError floating_bush_feed(reference,4,2)
    @test_throws ArgumentError equilibrate_floating_bush(reference,nt,nz;feed_cells=mask,fd_step=0.)
    println("Fixed seating force=$(fixed.state.bush_force_n); synthetic pose recovery=$(inverse.state.pose.displacement_m)")
end
