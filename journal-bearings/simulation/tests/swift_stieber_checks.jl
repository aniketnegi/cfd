using Test, LinearAlgebra
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../reynolds/src/swift_stieber.jl"))

ss_parameters(angle)=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
    radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=.6)
@testset "Supplied tensor pressure obstacle" begin
    @test !isdefined(@__MODULE__,:write_vtk)
    for angle in (0.,10.),graded in (false,true)
        p=ss_parameters(angle)
        nt,nz=32,16
        theta=graded ? [2pi*t+.25sin(2pi*t) for t in (0:nt-1)./nt] : nothing
        axial=graded ? [100*(t+.15sin(2pi*t)/(2pi)) for t in (0:nz)./nz] : nothing
        !isnothing(axial) && (axial[1]=0.;axial[end]=100.)
        grid=housing_grid(p,nt,nz;theta_nodes=theta,z_nodes=axial)
        before=deepcopy(grid)
        result=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277)
        @test grid==before
        @test result.grid===grid
        @test minimum(result.pressure_pa)>=0
        @test any(iszero,result.pressure_pa)
        @test result.diagnostics.positive_pressure_relative_residual<1e-10
        @test result.diagnostics.complementarity_relative<1e-10
        @test result.diagnostics.defect_accounting_relative<1e-10
        signed=solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277)
        @test norm(result.pressure_pa-max.(signed.pressure_pa,0.))/norm(result.pressure_pa)>.01
        flooded=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(3e7,3e7))
        full=solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(3e7,3e7))
        @test flooded.pressure_pa≈full.pressure_pa rtol=1e-11
        doubled=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0554)
        @test doubled.pressure_pa≈2result.pressure_pa rtol=1e-11
        @test doubled.axial_face_flow_m3_s≈result.axial_face_flow_m3_s rtol=1e-11
        rotated=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,eccentricity_angle_rad=.7)
        @test rotated.grid.nodes===grid.nodes
        @test rotated.grid.radial_gaps≈rotated.normal_gap_m.*(1e3/cosd(angle)) rtol=1e-13
    end
end
@testset "Pressure reference, pinned lines and quiescent limits" begin
    grid=housing_grid(ss_parameters(10.),32,16)
    base=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,pinned_cells=1:32:512)
    shifted=solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,
        cavitation_pressure_pa=101325.,end_pressure_pa=(101325.,101325.),pinned_cells=1:32:512)
    @test shifted.pressure_above_cavitation_pa≈base.pressure_pa rtol=1e-12
    @test shifted.axial_face_flow_m3_s≈base.axial_face_flow_m3_s rtol=1e-12
    @test all(shifted.pressure_pa[1:32:end].==101325.)
    @test shifted.diagnostics.pinned_line_exchange_m3_s≈base.diagnostics.pinned_line_exchange_m3_s rtol=1e-12
    idle=solve_swift_stieber(grid;rpm=0.,viscosity_pa_s=.0277)
    @test all(iszero,idle.pressure_pa)
    @test all(iszero,idle.axial_face_flow_m3_s)
    pressure=solve_swift_stieber(grid;rpm=0.,viscosity_pa_s=.0277,end_pressure_pa=(2e5,2e5))
    @test pressure.pressure_pa≈fill(2e5,512) rtol=1e-11
end
@testset "Invalid obstacle contracts" begin
    grid=housing_grid(ss_parameters(10.),16,8)
    for kwargs in ((;pinned_cells=[0]),(;pinned_cells=[129]),(;pinned_cells=collect(1:128)),
        (;cavitation_pressure_pa=NaN),(;end_pressure_pa=(-1.,0.)),
        (;cavitation_pressure_pa=1.),(;viscosity_pa_s=0.),(;rpm=Inf))
        defaults=(;rpm=2000.,viscosity_pa_s=.0277)
        @test_throws ArgumentError solve_swift_stieber(grid;merge(defaults,kwargs)...)
    end
    @test_throws ArgumentError solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,max_iterations=0)
    @test_throws ErrorException solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277,max_iterations=1)
end
