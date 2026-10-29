using Test,LinearAlgebra,SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../fsi/src/steady.jl"))
include(joinpath(@__DIR__,"../elasticity/src/vtk.jl"))

function check_fsi()
@testset "Normal-gap correction and unchanged rigid operator" begin
    for angle in (0.,10.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        g=housing_grid(p,16,8)
        kw=(;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(3e5,3e5),cavitation_pressure_pa=0.,eccentricity_angle_rad=3pi/2)
        old=assemble_conical_jfo(g;kw...);zero=assemble_conical_jfo(g;kw...,gap_correction_m=(t,z)->0.)
        @test old.A==zero.A && old.B==zero.B && old.b==zero.b
        @test old.normal_gap_m==zero.normal_gap_m
        @test_throws ArgumentError assemble_conical_jfo(g;kw...,gap_correction_m=(t,z)->-1e-3)
        concentric=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=0.)
        cg=housing_grid(concentric,16,8);h0=50e-6*cosd(angle);factor=1.2
        a=assemble_conical_jfo(cg;kw...)
        b=assemble_conical_jfo(cg;kw...,gap_correction_m=(t,z)->(factor-1)*h0)
        @test b.gt≈factor^3*a.gt rtol=1e-14
        @test b.gz≈factor^3*a.gz rtol=1e-14
        @test b.adv≈factor*a.adv rtol=1e-14
    end
end
@testset "Fixed-position housing one-way and two-way FSI" begin
    for angle in (0.,10.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        g=housing_grid(p,24,12)
        model=prepare_housing_fsi(g;youngs_modulus_pa=97e9,poisson_ratio=.3,wall_thickness_m=.005,n_radial=2)
        translation=repeat([1e-7,-2e-7,3e-7],size(model.solid.nodes,1))
        field=housing_gap_field(model.gap,translation)
        for (t,z) in model.gap.points
            expected=cosd(angle)*(cos(t)*1e-7-sin(t)*2e-7)+sind(angle)*3e-7
            @test field(t,z)≈expected atol=1e-20 rtol=1e-12
        end
        axial=repeat([0.,0.,1e-7],size(model.solid.nodes,1))
        @test model.gap.matrix*axial≈fill(sind(angle)*1e-7,length(model.gap.points)) atol=1e-20
        options=(;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(3e5,3e5),cavitation_pressure_pa=0.,eccentricity_angle_rad=3pi/2)
        one=solve_one_way_fsi(model;fluid_options=options)
        @test !one.pressure_recomputed_on_deformation
        @test one.minimum_predicted_gap_m>0
        @test one.conservative_gap_lower_bound_m>0
        controlled=solve_one_way_fsi(model;fluid_options=(;options...,max_iterations=2000,tolerance=1e-10))
        @test controlled.fluid.pressure_pa==one.fluid.pressure_pa
        @test norm(one.gap_correction_at_assembly_points_m,Inf)>1e-9
        @test one.solid.relative_residual<1e-9
        shifted=(;options...,end_pressure_pa=(401325.,401325.),cavitation_pressure_pa=101325.)
        reference=solve_one_way_fsi(model;fluid_options=shifted,external_pressure_pa=101325.)
        @test reference.solid.displacement_m≈one.solid.displacement_m rtol=1e-11
        @test_throws ArgumentError solve_one_way_fsi(model;fluid_options=(;options...,gap_correction_m=(t,z)->0.))
        two=solve_two_way_fsi(model;fluid_options=options)
        last=two.history[end]
        @test two.pressure_recomputed_on_deformation
        @test length(two.history)>1
        @test last.gap_defect_m<=1e-10 && last.pressure_relative_change<=1e-5
        @test last.solid_relative_residual<=1e-5
        @test minimum(two.fluid.normal_gap_m)>0
        @test norm(two.fluid.pressure_pa-one.fluid.pressure_pa)/norm(one.fluid.pressure_pa)>1e-5
        @test minimum(two.fluid.pressure_above_cavitation_pa)>=0
        @test 0<=minimum(two.fluid.fill_fraction)<=maximum(two.fluid.fill_fraction)<=1
        @test two.fluid.diagnostics.relative_global_imbalance<1e-10
        gap=housing_gap_field(model.gap,two.solid.displacement_m)
        a=assemble_conical_jfo(g;options...,gap_correction_m=gap)
        @test a.normal_gap_m==two.fluid.normal_gap_m
        @test two.conservative_gap_lower_bound_m>0
        shifted_two=solve_two_way_fsi(model;fluid_options=shifted,external_pressure_pa=101325.)
        @test shifted_two.solid.displacement_m≈two.solid.displacement_m rtol=1e-10
        @test shifted_two.fluid.pressure_above_cavitation_pa≈two.fluid.pressure_above_cavitation_pa rtol=1e-10
        stiff=prepare_housing_fsi(g;youngs_modulus_pa=97e13,poisson_ratio=.3,wall_thickness_m=.005,n_radial=2)
        rigid_limit=solve_two_way_fsi(stiff;fluid_options=options)
        @test norm(rigid_limit.fluid.pressure_pa-one.fluid.pressure_pa)/norm(one.fluid.pressure_pa)<5e-5
        tighter=solve_two_way_fsi(model;fluid_options=options,relaxation=.8,gap_tolerance_m=3e-11,
            pressure_tolerance=3e-6,solid_tolerance=3e-6)
        @test norm(two.fluid.pressure_pa-tighter.fluid.pressure_pa)/norm(tighter.fluid.pressure_pa)<1e-4
        @test norm(model.gap.matrix*(two.solid.displacement_m-tighter.solid.displacement_m),Inf)<1e-9
        @test_throws ArgumentError solve_two_way_fsi(model;fluid_options=options,relaxation=0.)
        @test_throws ErrorException solve_two_way_fsi(model;fluid_options=options,max_iterations=1)
        println("FSI angle=$angle iterations=$(length(two.history)) peak=$(maximum(two.fluid.pressure_pa)) maxnormal=$(maximum(two.gap_correction_at_assembly_points_m))")
    end
end
@testset "Matching 7mm solid drill and fitted gap feedback" begin
    p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
        semicone_angle_deg=10.,eccentricity_ratio=.6)
    g=surface_feed_grid(p,128,64;feed_diameter=7.,feed_z=50.)
    m=prepare_housing_fsi(g;youngs_modulus_pa=97e9,poisson_ratio=.3,wall_thickness_m=.005,
        n_radial=1,feed_diameter_m=.007,feed_z_m=.05)
    @test m.solid.minimum_jacobian>0
    @test m.solid.mesh.nodes[1:length(g.nodes),:] == reduce(vcat,permutedims.(collect.(g.nodes))).*1e-3
    @test norm(reshape(m.solid.passage_load_per_pa,3,:)[2,:],Inf)<1e-14
    s=g.surface;n=length(g.nodes)
    for f in findall(==(:feed),s.boundary)
        edge=s.edge_nodes[f];quad=s.quads[s.owner[f]]
        k=findfirst(i->minmax(quad[i],quad[mod1(i+1,4)])==edge,1:4)
        a,b=quad[k],quad[mod1(k+1,4)]
        ids=[a,b,b+n,a+n];X=m.solid.nodes[ids,:]
        N,du,dv=interface_quad_shape(0.,0.);point=X'*N
        load=-cross(X'*du,X'*dv)
        @test dot(load,[point[1],0.,point[3]-.05])>0
    end
    options=(;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),cavitation_pressure_pa=0.,
        feed_pressure_pa=5e5,eccentricity_angle_rad=3pi/2)
    one=solve_one_way_fsi(m;fluid_options=options)
    corrected=housing_gap_field(m.gap,one.solid.displacement_m)
    a=assemble_conical_jfo(g;options...);b=assemble_conical_jfo(g;options...,gap_correction_m=corrected)
    @test a.A!=b.A && a.B!=b.B
    @test a.faces.advective!=b.faces.advective
    @test minimum(b.normal_gap_m)>0
    zero=assemble_conical_jfo(g;options...,gap_correction_m=(t,z)->0.)
    @test a.A==zero.A && a.B==zero.B && a.b==zero.b
    @test a.normal_gap_m==zero.normal_gap_m
    @test one.solid.relative_residual<1e-9
    mktempdir() do dir
        path=write_elasticity_vtk(joinpath(dir,"solid.vtk"),m.solid,one.solid)
        lines=readlines(path)
        @test "POINT_DATA $(size(m.solid.nodes,1))" in lines
        @test "CELL_TYPES $(length(m.solid.cells))" in lines
        @test "VECTORS displacement_m double" in lines
        @test_throws ArgumentError write_elasticity_vtk(path,m.solid,one.solid)
    end
    @test_throws ArgumentError assemble_conical_jfo(g;options...,gap_correction_m=(t,z)->-1e-3)
end
end
if abspath(PROGRAM_FILE)==@__FILE__;check_fsi();end
