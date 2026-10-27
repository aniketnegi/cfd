using Test
include(joinpath(@__DIR__, "../2d/generate.jl"))

const BASE_ARGS = ["--geometry", "unported", "--journal-radius-mm", "50",
                   "--length-mm", "100", "--clearance-mm", "0.05",
                   "--cone-angle-deg", "10", "--eccentricity-ratio", "0.6"]

@testset "2D mesh CLI" begin
    @test !isdefined(@__MODULE__, :solve_conical_jfo)
    @test all(id.name != "Gmsh" for id in keys(Base.loaded_modules))
    @test occursin("NAME\n", USAGE)
    @test occursin("SYNOPSIS\n", USAGE)
    @test occursin("-h, --help", USAGE)
    @test isnothing(main(["-h"]))
    mktempdir() do dir
        path = joinpath(dir, "unported.vtk")
        args = vcat(BASE_ARGS, ["--output", path, "--n-theta", "16", "--n-axial", "4"])
        cfg = parse_grid_args(args)
        @test (cfg.n_theta, cfg.n_axial) == (16, 4)
        equals_path = joinpath(dir, "equals.vtk")
        equals = parse_grid_args(vcat(BASE_ARGS, ["--output=$equals_path",
                                                    "--n-theta=16", "--n-axial=4"]))
        @test (equals.n_theta, equals.n_axial) == (16, 4)
        main(args)
        @test occursin("POLYGONS 64 320", read(path, String))
        @test_throws ArgumentError main(args) # no accidental overwrite

        fed_path = joinpath(dir, "fed.vtk")
        fed = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
                   ["--output", fed_path, "--n-theta", "128", "--n-axial", "64",
                    "--feed-diameter-mm", "5", "--feed-z-mm", "50"])
        main(fed)
        @test occursin("SCALARS feed_rim int 1", read(fed_path, String))
        @test occursin("LINES 32 96", read(fed_path, String))
        @test !occursin("SCALARS feed_cell int 1", read(fed_path, String))

        alternate_path = joinpath(dir, "alternate-counts.vtk")
        alternate = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
            ["--output", alternate_path, "--n-theta", "160", "--n-axial", "52",
             "--feed-diameter-mm", "5", "--feed-z-mm", "50"])
        main(alternate)
        @test occursin("POLYGONS 8544 42720", read(alternate_path, String))
        @test_throws ArgumentError main(alternate)

        varied_path = joinpath(dir, "varied-geometry.vtk")
        varied = ["--geometry", "surface-feed", "--journal-radius-mm", "60",
            "--length-mm", "120", "--clearance-mm", "0.08", "--cone-angle-deg", "15",
            "--eccentricity-ratio", "0.6", "--n-theta", "160", "--n-axial", "64",
            "--feed-diameter-mm", "6", "--feed-z-mm", "70.3", "--output", varied_path]
        varied_cfg = parse_grid_args(varied)
        @test (varied_cfg.p.journal_midpoint_radius, varied_cfg.p.bearing_length,
            varied_cfg.p.radial_clearance, varied_cfg.p.semicone_angle_deg) == (60.,120.,.08,15.)
        main(varied)
        lines = readlines(varied_path)
        points = findfirst(line -> startswith(line,"POINTS "), lines)
        count = parse(Int,split(lines[points])[2])
        vertices = [Tuple(parse.(Float64,split(line))) for line in lines[points+1:points+count]]
        @test vertices[1][1] ≈ 60.08 + 60tan(deg2rad(15.)) atol=1e-11
        @test vertices[1][2:3] == (0.,0.)
        rim = findfirst(line -> startswith(line,"LINES "), lines)
        edges = parse(Int,split(lines[rim])[2])
        @test all(lines[rim+1:rim+edges]) do line
            all(parse.(Int,split(line)[2:end]).+1) do node
                x,y,z = vertices[node]
                y > 0 && isapprox(x^2+(z-70.3)^2,9.; atol=1e-10, rtol=0)
            end
        end
        unavailable_path = joinpath(dir,"unavailable-transition.vtk")
        unavailable = copy(varied)
        unavailable[findfirst(==("--output"),unavailable)+1] = unavailable_path
        unavailable[findfirst(==("--feed-z-mm"),unavailable)+1] = "5"
        @test_throws r"physical support crosses domain boundary" main(unavailable)
        @test !ispath(unavailable_path)

        rejected_path = joinpath(dir, "rejected-quality.vtk")
        rejected = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
            ["--output", rejected_path, "--n-theta", "12", "--n-axial", "6",
             "--feed-diameter-mm", "5", "--feed-z-mm", "50"])
        @test_throws r"does not meet the surface-quality limits" main(rejected)
        @test !ispath(rejected_path)
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "bad.txt"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--n-theta", "abc"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--n-axial", "0"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--bogus", "1"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--length-mm=1"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--feed-z-mm", "50"]))
        fed_missing = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
                           ["--output", "missing.vtk", "--feed-diameter-mm", "5"])
        @test_throws ArgumentError parse_grid_args(fed_missing)

        cases = joinpath(dir, "cases")
        mkpath(cases)
        named = vcat(BASE_ARGS, ["--experiment", "trial-a", "--n-theta", "16", "--n-axial", "4"])
        main(named; cases_root=cases)
        exp = load_experiment("trial-a"; cases_root=cases)
        @test (exp.ntheta, exp.nz, exp.variant) == (16, 4, "unported")
        @test exp.p.semicone_angle_deg == 10.0
        @test isfile(joinpath(exp.dir, "experiment.toml"))
        @test_throws ArgumentError main(named; cases_root=cases) # immutable named case
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS,
            ["--experiment", "../unsafe"]); cases_root=cases)
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS,
            ["--experiment", "trial-b", "--output", "both.vtk"]); cases_root=cases)
        @test_throws ArgumentError parse_grid_args(BASE_ARGS; cases_root=cases)
        open(exp.mesh, "a") do io
            println(io, "altered mesh")
        end
        @test_throws ArgumentError load_experiment("trial-a"; cases_root=cases)

        # The case loader does not support fitted-land-direct-v1.
        # Reject the request before creating an output directory.
        named_feed = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
                         ["--experiment", "trial-feed", "--n-theta", "128", "--n-axial", "64",
                          "--feed-diameter-mm", "5", "--feed-z-mm", "50"])
        @test_throws ArgumentError main(named_feed; cases_root=cases)
        @test !ispath(joinpath(cases, "trial-feed"))

        # Schema-1 manifests without a layout identify cell-centre masks.
        legacy_dir = joinpath(cases, "legacy-mask"); mkdir(legacy_dir)
        legacy_grid = housing_grid(cfg.p,128,64)
        write_vtk(joinpath(legacy_dir,"mesh.vtk"),legacy_grid;feed_cells=falses(length(legacy_grid.quads)))
        save_experiment_manifest(legacy_dir,"legacy-mask","surface-feed",cfg.p,128,64;
            feed_diameter=5.,feed_z=50.)
        @test load_experiment("legacy-mask";cases_root=cases).layout == "cell-centre-mask"

        # A named fitted request must not leave a partially written case.
        unsupported = vcat(["--geometry","surface-feed"],BASE_ARGS[3:end],
            ["--experiment","unsupported","--n-theta","256","--n-axial","64",
             "--feed-diameter-mm","5","--feed-z-mm","50"])
        @test_throws ArgumentError main(unsupported;cases_root=cases)
        @test !ispath(joinpath(cases,"unsupported"))
    end
end
