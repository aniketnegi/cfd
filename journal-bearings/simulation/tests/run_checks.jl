using Test
include(joinpath(@__DIR__, "../reynolds/run.jl"))

function make_case(root, name; angle=0.0, fed=false, ntheta=16, nz=8)
    p = BearingParameters(journal_midpoint_radius=50.0, bearing_length=100.0,
                          radial_clearance=0.05, semicone_angle_deg=angle,
                          eccentricity_ratio=0.6)
    dir = case_directory(name; cases_root=root)
    mkdir(dir)
    grid = housing_grid(p, ntheta, nz)
    # The rejection test needs a valid manifest and mesh checksum; feed
    # footprint geometry itself is exercised by mesh/tests/meshing.jl.
    write_vtk(joinpath(dir, "mesh.vtk"), grid;
              feed_cells=fed ? falses(length(grid.quads)) : nothing)
    save_experiment_manifest(dir, name, fed ? "surface-feed" : "unported", p, ntheta, nz;
                             feed_diameter=fed ? 5.0 : nothing,
                             feed_z=fed ? 50.0 : nothing)
    return dir
end

function capture_run(f)
    mktemp() do outpath, out
        mktemp() do errpath, err
            value = redirect_stdout(out) do
                redirect_stderr(err) do
                    f()
                end
            end
            flush(out); flush(err)
            return (;value, stdout=read(outpath, String), stderr=read(errpath, String))
        end
    end
end

@testset "Named Reynolds run CLI" begin
    help = capture_run(() -> main_run(["--help"]))
    @test isnothing(help.value)
    @test occursin("usage:", help.stdout)
    @test occursin("--viscosity-pa-s", help.stdout)
    @test isempty(help.stderr)
    mktempdir() do root
        dir = make_case(root, "test-cylinder")
        args = ["--experiment", "test-cylinder", "--run", "baseline",
                "--rpm", "2000", "--viscosity-pa-s=0.02752"]
        cfg = parse_run_args(args; cases_root=root)
        @test (cfg.exp.ntheta, cfg.exp.nz, cfg.rpm, cfg.mu) == (16, 8, 2000, 0.02752)
        for invalid in (
            String[], ["--unknown"], ["--rpm"],
            [args; "--rpm"; "1000"], [args; "--experiment=test-cylinder"],
            [args; "--run=other"], [args; "--viscosity-pa-s=.03"],
            [args; "--sweep"; "--sweep"], [args; "--replace"; "--replace"],
            [args; "-v"; "-q"], [args; "--quiet=true"],
            [args; "unexpected-positional"])
            @test_throws ArgumentError parse_run_args(invalid; cases_root=root)
        end
        for speed in ("NaN", "Inf", "not-a-number")
            @test_throws ArgumentError parse_run_args(["--experiment=test-cylinder",
                "--run=invalid", "--rpm=$speed", "--viscosity-pa-s=.02752"]; cases_root=root)
        end
        for viscosity in ("NaN", "Inf", "-0.1", "0", "not-a-number")
            @test_throws ArgumentError parse_run_args(["--experiment=test-cylinder",
                "--run=invalid", "--rpm=-2000", "--viscosity-pa-s=$viscosity"]; cases_root=root)
        end
        negative = parse_run_args(["--experiment=test-cylinder", "--run=reverse",
            "--rpm", "-2000", "--viscosity-pa-s", "0.02752"]; cases_root=root)
        @test negative.rpm == -2000
        normal = capture_run(() -> main_run(args; cases_root=root))
        @test isempty(normal.stdout)
        @test occursin("run test-cylinder/baseline grid=16x8", normal.stderr)
        @test occursin("pressure min=", normal.stderr)
        @test occursin("output ", normal.stderr)
        @test !occursin("algebraic_residual_inf_m3_s=", normal.stderr)
        rundir = joinpath(dir, "runs", "baseline")
        @test isfile(joinpath(rundir, "run.toml"))
        @test isfile(joinpath(rundir, "full-film.vtk"))
        @test occursin("SCALARS pressure_gauge_pa double 1", read(joinpath(rundir, "full-film.vtk"), String))
        diagnostic = TOML.parsefile(joinpath(rundir, "diagnostics.toml"))
        @test diagnostic["relative_global_imbalance"] < 1e-11
        @test diagnostic["min_pressure_gauge_pa"] < 0
        @test length(diagnostic["outward_end_flow_m3_s"]) == 2
        saved_hash = mesh_hash(joinpath(rundir, "run.toml"))
        @test_throws ArgumentError main_run(args; cases_root=root) # no overwrite
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--rpm", "2100"]; cases_root=root)
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--viscosity-pa-s", "0.04"]; cases_root=root)
        main_run(["--experiment", "test-cylinder", "--run", "baseline", "--replace"];
                 cases_root=root)  # read both operating inputs from saved run.toml
        @test mesh_hash(joinpath(rundir, "run.toml")) == saved_hash
        artifacts = ("run.toml", "full-film.vtk", "diagnostics.toml")
        hashes = [mesh_hash(joinpath(rundir, name)) for name in artifacts]
        rerun = ["--experiment=test-cylinder", "--run=baseline", "--replace"]
        quiet = capture_run(() -> main_run([rerun; "-q"]; cases_root=root))
        @test isempty(quiet.stdout) && isempty(quiet.stderr)
        verbose = capture_run(() -> main_run([rerun; "-v"]; cases_root=root))
        @test isempty(verbose.stdout)
        @test occursin("solve ", verbose.stderr)
        @test occursin("algebraic_residual_inf_m3_s=", verbose.stderr)
        @test hashes == [mesh_hash(joinpath(rundir, name)) for name in artifacts]
        rejected = capture_run(() -> run_cli([args; "-q"]; cases_root=root))
        @test rejected.value == 2
        @test isempty(rejected.stdout)
        @test occursin("result already exists", rejected.stderr)
        @test hashes == [mesh_hash(joinpath(rundir, name)) for name in artifacts]
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
        sweep_hash = mesh_hash(sweep_path)
        swept = capture_run(() -> main_run(["--experiment=test-cylinder", "--run=second",
            "--replace", "--sweep", "-v"]; cases_root=root))
        @test occursin("sweep 4x2, 8x4, 16x8", swept.stderr)
        @test mesh_hash(sweep_path) == sweep_hash
        open(joinpath(dir, "mesh.vtk"), "a") do io
            print(io, "corrupted")
        end
        @test_throws ArgumentError parse_run_args(["--experiment", "test-cylinder",
            "--run", "baseline", "--replace"]; cases_root=root)

        blocked = make_case(root, "blocked")
        write(joinpath(blocked, "runs"), "not a directory")
        failure = capture_run(() -> run_cli(["--experiment=blocked", "--run=trial",
            "--rpm=2000", "--viscosity-pa-s=.02752", "-q"]; cases_root=root))
        @test failure.value == 1
        @test isempty(failure.stdout)
        @test occursin("run.jl:", failure.stderr)
        @test !occursin("Stacktrace", failure.stderr)
        @test read(joinpath(blocked, "runs"), String) == "not a directory"

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

@testset "Run entry-point exit codes and streams" begin
    project = normpath(joinpath(@__DIR__, "../.."))
    script = joinpath(project, "simulation/reynolds/run.jl")
    for (args, expected) in ((["--help"], 0), (String[], 2), (["--unknown"], 2))
        result = capture_run() do
            command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script $args`
            process = run(pipeline(ignorestatus(command); stdout=stdout, stderr=stderr))
            process.exitcode
        end
        @test result.value == expected
        if expected == 0
            @test occursin("usage:", result.stdout)
            @test isempty(result.stderr)
        else
            @test isempty(result.stdout)
            @test occursin("run.jl:", result.stderr)
            @test !occursin("Stacktrace", result.stderr)
        end
    end
end
