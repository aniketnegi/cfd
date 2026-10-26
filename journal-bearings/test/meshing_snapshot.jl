# Record or compare every returned mesh field and the exact VTK bytes.
# Run the same script at the source checkpoint and after a structural refactor.
using SHA, Serialization, TOML, Test
include(joinpath(@__DIR__, "../bearing.jl"))
include(joinpath(@__DIR__, "../mesh/2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../mesh/2d/conical/surface_feed.jl"))
include(joinpath(@__DIR__, "../mesh/vtk.jl"))

function field_hash(value)
    io = IOBuffer()
    serialize(io, value)
    return bytes2hex(sha256(take!(io)))
end

function mesh_snapshot()
    records = Dict{String,Any}()
    mktempdir() do dir
        for angle in (0.0, 10.0), level in (1, 2, 4)
            p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
                radial_clearance=.05, semicone_angle_deg=angle, eccentricity_ratio=.6)
            name = "regular-angle$(Int(angle))-level$level"
            grid = housing_grid(p, 128level, 64level)
            path = joinpath(dir, "$name.vtk")
            write_vtk(path, grid)
            records[name] = Dict(
                "fields" => Dict(string(k) => field_hash(v) for (k, v) in pairs(grid)),
                "vtk_sha256" => bytes2hex(sha256(read(path))),
                "cells" => length(grid.quads), "nodes" => length(grid.nodes))
        end
        p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
            radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
        for level in (1, 2, 4)
            name = "inlet-level$level"
            grid = surface_feed_grid(p, 128level, 64level; feed_diameter=5., feed_z=50.)
            path = joinpath(dir, "$name.vtk")
            write_vtk(path, grid)
            record = Dict{String,Any}(
                "fields" => Dict(string(k) => field_hash(v) for (k, v) in pairs(grid)),
                "surface_fields" => Dict(string(k) => field_hash(v) for (k, v) in pairs(grid.surface)),
                "vtk_sha256" => bytes2hex(sha256(read(path))),
                "cells" => length(grid.quads), "nodes" => length(grid.nodes),
                "quality" => Dict(string(k) => v for (k, v) in pairs(grid.quality)))
            records[name] = record
            println("$name: $(record["cells"]) cells, $(record["nodes"]) nodes; $(grid.quality)")
        end
    end
    return records
end

function main(args)
    length(args) == 2 && args[1] in ("record", "compare") ||
        error("usage: julia --project=. test/meshing_snapshot.jl record|compare SNAPSHOT.toml")
    mode, path = args
    mode == "record" && ispath(path) && error("refusing to overwrite baseline: $path")
    mode == "compare" && !isfile(path) && error("missing baseline: $path")
    actual = mesh_snapshot()
    if mode == "record"
        mkpath(dirname(path))
        open(path, "w") do io
            TOML.print(io, actual; sorted=true)
        end
        println("Recorded all mesh fields and exports in $path")
    else
        expected = TOML.parsefile(path)
        @testset "Unchanged mesh fields and VTK output" begin
            @test Set(keys(actual)) == Set(keys(expected))
            for name in sort!(collect(keys(expected)))
                @testset "$name" begin
                    @test Set(keys(actual[name])) == Set(keys(expected[name]))
                    for field in keys(expected[name])
                        @test isequal(actual[name][field], expected[name][field])
                    end
                end
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
