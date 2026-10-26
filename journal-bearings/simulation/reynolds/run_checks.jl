using Test
include(joinpath(@__DIR__, "run.jl"))

function make_case(root, name; angle=0.0, fed=false, ntheta=16, nz=8)
    p = BearingParameters(journal_midpoint_radius=50.0, bearing_length=100.0,
                          radial_clearance=0.05, semicone_angle_deg=angle,
                          eccentricity_ratio=0.6)
    dir = case_directory(name; cases_root=root)
    mkdir(dir)
    grid = housing_grid(p, ntheta, nz)
    # The rejection test needs a valid manifest and mesh checksum; feed
    # footprint geometry itself is exercised by mesh/2d/conical/checks.jl.
    write_vtk(joinpath(dir, "mesh.vtk"), grid;
              feed_cells=fed ? falses(length(grid.quads)) : nothing)
    save_experiment_manifest(dir, name, fed ? "surface-feed" : "unported", p, ntheta, nz;
                             feed_diameter=fed ? 5.0 : nothing,
                             feed_z=fed ? 50.0 : nothing)
    return dir
end

@testset "Named Reynolds run CLI" begin
    @test occursin("SYNOPSIS", RUN_USAGE)
    @test isnothing(main_run(["--help"]))
    mktempdir() do root
        dir = make_case(root, "test-cylinder")
        args = ["--experiment", "test-cylinder", "--run", "baseline",
                "--rpm", "2000", "--viscosity-pa-s=0.02752"]
        cfg = parse_run_args(args; cases_root=root)
        @test (cfg.exp.ntheta, cfg.exp.nz, cfg.rpm, cfg.mu) == (16, 8, 2000, 0.02752)
        main_run(args; cases_root=root)
        rundir = joinpath(dir, "runs", "baseline")
        @test isfile(joinpath(rundir, "run.toml"))
        @test isfile(joinpath(rundir, "full-film.vtk"))
        @test occursin("SCALARS pressure_gauge_pa double 1", read(joinpath(rundir, "full-film.vtk"), String))
        diagnostic = TOML.parsefile(joinpath(rundir, "diagnostics.toml"))
        @test diagnostic["relative_global_imbalance"] < 1e-11
        @test diagnostic["min_pressure_gauge_pa"] < 0
        saved_hash = mesh_hash(joinpath(rundir, "run.toml"))
        @test_throws ArgumentError main_run(args; cases_root=root) # no overwrite
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--rpm", "2100"]; cases_root=root)
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--viscosity-pa-s", "0.04"]; cases_root=root)
        main_run(["--experiment", "test-cylinder", "--run", "baseline", "--replace"];
                 cases_root=root)  # read both operating inputs from saved run.toml
        @test mesh_hash(joinpath(rundir, "run.toml")) == saved_hash
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "new"]; cases_root=root)
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "new", "--rpm", "2000", "--viscosity-pa-s", "0"]; cases_root=root)
        @test_throws ArgumentError parse_run_args(["--experiment", "../test-cylinder",
            "--run", "other"]; cases_root=root)
        main_run(["--experiment", "test-cylinder", "--run", "second",
                  "--rpm", "1000", "--viscosity-pa-s", "0.02752", "--sweep"];
                 cases_root=root) # parent runs/ already exists
        @test isfile(joinpath(dir, "runs", "second", "full-film.vtk"))
        sweep_path = joinpath(dir, "runs", "second", "convergence.tsv")
        @test length(readlines(sweep_path)) == 4 # header + 4×2, 8×4, 16×8
        @test_throws ArgumentError main_run(["--experiment", "test-cylinder",
            "--run", "second", "--sweep"]; cases_root=root)
        main_run(["--experiment", "test-cylinder", "--run", "second",
                  "--replace", "--sweep"]; cases_root=root)
        @test length(readlines(sweep_path)) == 4
        open(joinpath(dir, "mesh.vtk"), "a") do io
            print(io, "corrupted")
        end
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--replace"]; cases_root=root)

        make_case(root, "test-cone"; angle=10.0)
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cone",
            "--run", "baseline", "--rpm", "2000", "--viscosity-pa-s", "0.02752"];
            cases_root=root)
        make_case(root, "test-fed"; fed=true, ntheta=128, nz=64)
        @test_throws ArgumentError parse_run_args(["--experiment", "test-fed",
            "--run", "baseline", "--rpm", "2000", "--viscosity-pa-s", "0.02752"];
            cases_root=root)
    end
end
