# Fixed-position elastic-housing pilot; remote-sized sweeps require --full.
using Test,TOML,SHA,Serialization
include(joinpath(@__DIR__,"../../../bearing.jl"))
include(joinpath(@__DIR__,"../../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../../fsi/src/steady.jl"))
include(joinpath(@__DIR__,"../../elasticity/src/vtk.jl"))
include(joinpath(@__DIR__,"../../../cli.jl"))
const FSI_ROOT=normpath(joinpath(@__DIR__,"../../.."))

function fsi_write(io,row;header=false)
    header && println(io,join(string.(keys(row)),'\t'))
    println(io,join(values(row),'\t'));flush(io)
end
function fsi_position(model,ecc)
    p=model.grid.parameters
    parameters=BearingParameters(journal_midpoint_radius=p.journal_midpoint_radius,
        bearing_length=p.bearing_length,radial_clearance=p.radial_clearance,
        semicone_angle_deg=p.semicone_angle_deg,eccentricity_ratio=ecc)
    return (;model...,grid=(;model.grid...,parameters))
end
function fsi_metadata(out,full,both_elastic)
    ispath(out) && error("refusing to overwrite $out");mkpath(out)
    paths=["Project.toml","Manifest.toml","bearing.jl","cli.jl"]
    for base in ("simulation/reynolds/src","simulation/elasticity/src","simulation/fsi/src","simulation/verification/fsi","mesh")
        for (dir,_,files) in walkdir(joinpath(FSI_ROOT,base)),file in files
            endswith(file,".jl") && push!(paths,relpath(joinpath(dir,file),FSI_ROOT))
        end
    end
    metadata=Dict("scope"=>(both_elastic ? "Fixed-position elastic journal and housing; linear small-strain normal-gap coupling" : "Fixed-position, rigid journal and elastic housing; linear small-strain normal-gap coupling"),
        "full"=>full,"journal_midpoint_radius_m"=>.05,"length_m"=>.1,"radial_clearance_m"=>50e-6,
        "youngs_modulus_pa"=>97e9,"poisson_ratio"=>.3,"radial_wall_thickness_m"=>.005,
        "supports"=>"outer surface fixed, axial end faces free","material_status"=>"user-adopted homogeneous reference alloy, not recovered rig material",
        "feed_diameter_m"=>.007,"feed_z_m"=>.05,"feed_pressure_gauge_pa"=>5e5,
        "passage_wall_pressure"=>"uniform reservoir gauge pressure; no passage loss model",
        "end_pressure_gauge_pa"=>[0.,0.],"cavity_pressure_gauge_pa"=>0.,
        "rpm"=>2000.,"mu_pa_s"=>.0277,"offset_angle_rad"=>3pi/2,
        "elastic_journal"=>both_elastic,"journal_youngs_modulus_pa"=>210e9,"journal_poisson_ratio"=>.3,
        "journal_supports"=>"all displacement components fixed at both axial end faces",
        "journal_geometry"=>"solid nominal 100mm cylinder/frustum, no bore or overhangs; rigid pose in gap law",
        "journal_opening_trace"=>"uniform reservoir pressure beneath the actual housing opening",
        "load_control"=>false,"physical_validation"=>false,
        "source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(FSI_ROOT,p)))) for p in unique(paths)))
    open(joinpath(out,"inputs.toml"),"w") do io;TOML.print(io,metadata;sorted=true);end
end
function run_fsi(out;full=false,observer=nothing,both_elastic=false)
    fsi_metadata(out,full,both_elastic)
    first=true
    @testset "Fixed-position bearing FSI" begin
    open(joinpath(out,"results.tsv"),"w") do io
        for angle in (0.,10.),nt in (full ? (128,256) : (128,))
            p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
                semicone_angle_deg=angle,eccentricity_ratio=.6)
            println("FSI mesh angle=$angle background=$nt × $(nt÷2)")
            mesh_seconds=@elapsed g=surface_feed_grid(p,nt,nt÷2;feed_diameter=7.,feed_z=50.,observer)
            nr=nt==128 ? 2 : 4
            prepare_seconds=@elapsed begin
                base=prepare_housing_fsi(g;youngs_modulus_pa=97e9,poisson_ratio=.3,
                    wall_thickness_m=.005,n_radial=nr,feed_diameter_m=.007,feed_z_m=.05)
                if both_elastic
                    jn=nt÷4
                    base=add_elastic_journal(base;youngs_modulus_pa=210e9,poisson_ratio=.3,
                        n_theta=jn,n_axial=jn÷2,n_radial=jn÷8)
                end
            end
            println("Housing nodes=$(size(base.solid.nodes,1)) DOFs=$(size(base.solid.K,1)) free=$(length(base.solid.free)) prepared in $prepare_seconds seconds")
            if both_elastic
                ja=base.bodies.journal.assembly
                println("Journal nodes=$(size(ja.nodes,1)) DOFs=$(size(ja.K,1)) free=$(length(ja.free))")
            end
            for ecc in (full ? (.2,.6,.8) : (.6,))
                model=fsi_position(base,ecc)
                options=(;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),
                    cavitation_pressure_pa=0.,feed_pressure_pa=5e5,eccentricity_angle_rad=3pi/2)
                one_seconds=@elapsed one=solve_one_way_fsi(model;fluid_options=options)
                two_seconds=@elapsed two=solve_two_way_fsi(model;fluid_options=options,
                    observer=row->begin
                        println("coupling angle=$angle epsilon=$ecc iteration=$(row.iteration) gapdefect=$(row.gap_defect_m) pressurechange=$(row.pressure_relative_change)")
                        yield()
                    end)
                for (mode,result,seconds) in (("one-way",one,one_seconds),("two-way",two,two_seconds))
                    minimum_gap=mode=="one-way" ? one.minimum_predicted_gap_m : minimum(result.fluid.normal_gap_m)
                    row=(;angle_deg=angle,epsilon=ecc,nt,nz=nt÷2,n_radial=nr,
                        fluid_cells=length(g.quads),solid_nodes=size(base.solid.nodes,1),solid_dofs=size(base.solid.K,1),
                        journal_nodes=both_elastic ? size(base.bodies.journal.assembly.nodes,1) : 0,
                        journal_dofs=both_elastic ? size(base.bodies.journal.assembly.K,1) : 0,
                        journal_nt=both_elastic ? nt÷4 : 0,journal_nz=both_elastic ? nt÷8 : 0,
                        journal_nr=both_elastic ? nt÷32 : 0,
                        mode,mesh_seconds,prepare_seconds,solve_seconds=seconds,
                        peak_pressure_pa=maximum(result.fluid.pressure_pa),minimum_gap_m=minimum_gap,
                        conservative_gap_lower_bound_m=result.conservative_gap_lower_bound_m,
                        minimum_fill=minimum(result.fluid.fill_fraction),net_feed_m3_s=result.fluid.diagnostics.net_feed_flow_m3_s,
                        max_gap_correction_m=maximum(result.gap_correction_at_assembly_points_m),
                        min_gap_correction_m=minimum(result.gap_correction_at_assembly_points_m),
                        max_housing_normal_m=maximum(base.gap.matrix*result.solids.housing.displacement_m),
                        min_housing_normal_m=minimum(base.gap.matrix*result.solids.housing.displacement_m),
                        max_journal_outward_normal_m=both_elastic ? maximum(base.bodies.journal.gap_matrix*result.solids.journal.displacement_m) : 0.,
                        min_journal_outward_normal_m=both_elastic ? minimum(base.bodies.journal.gap_matrix*result.solids.journal.displacement_m) : 0.,
                        elastic_residual=maximum(s.relative_residual for s in values(result.solids)),
                        strain_energy_j=sum(s.strain_energy_j for s in values(result.solids)),
                        housing_energy_j=result.solids.housing.strain_energy_j,
                        journal_energy_j=both_elastic ? result.solids.journal.strain_energy_j : 0.,
                        maximum_von_mises_pa=maximum(s.maximum_von_mises_pa for s in values(result.stresses)),
                        maximum_strain=maximum(s.maximum_engineering_strain for s in values(result.stresses)),
                        liquid_global_balance=result.fluid.diagnostics.relative_global_imbalance,
                        coupling_iterations=mode=="one-way" ? 1 : length(two.history),peak_process_rss_bytes=Sys.maxrss())
                    fsi_write(io,row;header=first);first=false
                    @test minimum_gap>0 && row.maximum_strain<.01
                    @test result.conservative_gap_lower_bound_m>0
                    @test result.fluid.diagnostics.relative_global_imbalance<1e-10
                    @test row.elastic_residual<=(mode=="one-way" ? 1e-9 : 1e-5)
                    stem="$mode-$angle-$ecc-$nt"
                    serialize(joinpath(out,stem*".jls"),result)
                    write_elasticity_vtk(joinpath(out,stem*"-solid.vtk"),base.solid,result.solid)
                    write_conical_jfo_vtk(joinpath(out,stem*"-fluid.vtk"),result.fluid)
                    if both_elastic
                        ja=base.bodies.journal.assembly;js=result.solids.journal
                        write_elasticity_vtk(joinpath(out,stem*"-journal.vtk"),ja,js)
                        open(joinpath(out,stem*"-journal.tsv"),"w") do stream
                            println(stream,"x_m\ty_m\tz_m\tux_m\tuy_m\tuz_m")
                            for k in axes(ja.nodes,1)
                                println(stream,join([ja.nodes[k,:];js.displacement_m[3k-2:3k]],'\t'))
                            end
                        end
                    end
                    open(joinpath(out,stem*"-fluid.tsv"),"w") do stream
                        println(stream,"theta_rad\tz_m\tpressure_pa\tfill_fraction\tnormal_gap_m")
                        s=g.surface
                        for k in eachindex(result.fluid.pressure_pa)
                            println(stream,join((s.centre[k,1]/s.R,s.centre[k,2]*s.ca,
                                result.fluid.pressure_pa[k],result.fluid.fill_fraction[k],result.fluid.normal_gap_m[k]),'\t'))
                        end
                    end
                    open(joinpath(out,stem*"-solid.tsv"),"w") do stream
                        println(stream,"x_m\ty_m\tz_m\tux_m\tuy_m\tuz_m")
                        for k in axes(base.solid.nodes,1)
                            println(stream,join([base.solid.nodes[k,:];result.solid.displacement_m[3k-2:3k]],'\t'))
                        end
                    end
                end
                open(joinpath(out,"coupling-$angle-$ecc-$nt.tsv"),"w") do stream
                    for (i,row) in enumerate(two.history);fsi_write(stream,row;header=i==1);end
                end
                println("Finished angle=$angle epsilon=$ecc: one-way=$one_seconds seconds, two-way=$two_seconds seconds")
            end
            base=nothing;model=nothing;one=nothing;two=nothing;ja=nothing;js=nothing
            GC.gc()
        end
    end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    exit(BearingCLI.study_main(;full=true,description="Run fixed-position cylinder/cone elastic-housing one-/two-way FSI.") do cfg
        observer=stderr isa Base.TTY && !cfg.quiet && !cfg.verbose ? _->yield() : nothing
        run_fsi(cfg.output;full=cfg.full,observer)
    end)
end
