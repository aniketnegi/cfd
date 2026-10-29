using Test,LinearAlgebra,SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../fsi/src/steady.jl"))

function check_both_fsi()
@testset "Two-solid fixed-position FSI" begin
    for angle in (0.,10.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=.6)
        g=housing_grid(p,24,12)
        base=prepare_housing_fsi(g;youngs_modulus_pa=97e9,poisson_ratio=.3,wall_thickness_m=.005,n_radial=2)
        model=add_elastic_journal(base;youngs_modulus_pa=210e9,poisson_ratio=.3,n_theta=16,n_axial=8,n_radial=2)
        options=(;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(3e5,3e5),cavitation_pressure_pa=0.,eccentricity_angle_rad=3pi/2)
        one=solve_one_way_fsi(model;fluid_options=options)
        house=solve_one_way_fsi(base;fluid_options=options)
        @test one.fluid.pressure_pa==house.fluid.pressure_pa
        @test one.solids.housing.displacement_m==house.solid.displacement_m
        @test one.solids.journal.strain_energy_j>0
        u=(housing=one.solids.housing.displacement_m,journal=one.solids.journal.displacement_m)
        @test fsi_gap_values(model,u)≈base.gap.matrix*u.housing-model.bodies.journal.gap_matrix*u.journal rtol=1e-15
        @test one.conservative_gap_lower_bound_m>0
        common=(housing=repeat([1e-7,-2e-7,3e-7],size(base.solid.nodes,1)),
            journal=repeat([1e-7,-2e-7,3e-7],size(model.bodies.journal.assembly.nodes,1)))
        @test norm(fsi_gap_values(model,common),Inf)<1e-20
        @test all(iszero,one.solids.journal.displacement_m[model.bodies.journal.assembly.fixed])
        two=solve_two_way_fsi(model;fluid_options=options)
        @test two.history[end].solid_relative_residual<=1e-5
        @test all(s->s.relative_residual<=1e-5,values(two.solids))
        @test all(s->s.maximum_engineering_strain<.01,values(two.stresses))
        @test two.fluid.diagnostics.relative_global_imbalance<1e-10
        @test minimum(two.fluid.normal_gap_m)>0
        u=fsi_map(s->s.displacement_m,two.solids)
        a=assemble_conical_jfo(g;options...,gap_correction_m=fsi_gap_field(model,u))
        @test a.normal_gap_m==two.fluid.normal_gap_m
        shifted=(;options...,end_pressure_pa=(401325.,401325.),cavitation_pressure_pa=101325.)
        datum=solve_two_way_fsi(model;fluid_options=shifted,external_pressure_pa=101325.)
        @test datum.fluid.pressure_above_cavitation_pa≈two.fluid.pressure_above_cavitation_pa rtol=1e-10
        @test datum.solids.journal.displacement_m≈two.solids.journal.displacement_m rtol=1e-10
        @test datum.solids.housing.displacement_m≈two.solids.housing.displacement_m rtol=1e-10
        for state in values(two.solids)
            force=vec(sum(reshape(state.load_n+state.reaction_n,3,:),dims=2))
            scale=sum(norm(view(state.load_n,3i-2:3i)) for i in 1:length(state.load_n)÷3)
            @test norm(force)<1e-4*scale
        end
        stiff=add_elastic_journal(base;youngs_modulus_pa=210e13,poisson_ratio=.3,n_theta=16,n_axial=8,n_radial=2)
        limit=solve_two_way_fsi(stiff;fluid_options=options)
        housing=solve_two_way_fsi(base;fluid_options=options)
        @test norm(limit.fluid.pressure_pa-housing.fluid.pressure_pa)/norm(housing.fluid.pressure_pa)<1e-4
        tighter=solve_two_way_fsi(model;fluid_options=options,relaxation=.8,gap_tolerance_m=3e-11,pressure_tolerance=3e-6,solid_tolerance=3e-6)
        @test norm(two.fluid.pressure_pa-tighter.fluid.pressure_pa)/norm(tighter.fluid.pressure_pa)<1e-4
        @test_throws ArgumentError add_elastic_journal(model;youngs_modulus_pa=210e9,poisson_ratio=.3,n_theta=16,n_axial=8,n_radial=2)
        println("Both-elastic angle=$angle iterations=$(length(two.history)) peak=$(maximum(two.fluid.pressure_pa))")
    end
end
end
if abspath(PROGRAM_FILE)==@__FILE__;check_both_fsi();end
