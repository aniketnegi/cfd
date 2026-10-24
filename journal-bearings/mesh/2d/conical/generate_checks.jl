using Test
include(joinpath(@__DIR__, "generate.jl"))

const BASE_ARGS = ["--geometry", "unported", "--journal-radius-mm", "50",
                   "--length-mm", "100", "--clearance-mm", "0.05",
                   "--cone-angle-deg", "10", "--eccentricity-ratio", "0.6"]

@testset "2D mesh CLI" begin
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
        @test occursin("SCALARS feed_cell int 1", read(fed_path, String))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "bad.txt"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--n-theta", "abc"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--n-axial", "0"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--bogus", "1"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--length-mm=1"]))
        @test_throws ArgumentError parse_grid_args(vcat(BASE_ARGS, ["--output", "x.vtk", "--feed-z-mm", "50"]))
        fed_missing = vcat(["--geometry", "surface-feed"], BASE_ARGS[3:end],
                           ["--output", "missing.vtk", "--feed-diameter-mm", "5"])
        @test_throws ArgumentError parse_grid_args(fed_missing)
    end
end
