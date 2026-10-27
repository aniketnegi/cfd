using Test, SHA

@testset "Single-mesh examples use fresh outputs" begin
    project = abspath(joinpath(@__DIR__, "../.."))
    examples = (("mesh_unported",131072,131584,0),
        ("mesh_adjusted_unported",8192,8320,0),
        ("mesh_inlet",134252,134823,120))
    mktempdir() do directory
        for (name,cells,vertices,rim_edges) in examples
            @testset "$name" begin
                script = joinpath(@__DIR__,"../examples/$name.jl")
                path = joinpath(directory,"$name.vtk")
                command = `$(Base.julia_cmd()) --project=$project --check-bounds=yes $script $path`
                output = IOBuffer()
                process = run(pipeline(ignorestatus(command),stdout=output,stderr=output))
                message = String(take!(output))
                success(process) || println(message)
                @test success(process)
                @test occursin("Wrote ",message)
                @test isfile(path)

                lines = readlines(path)
                @test "POINTS $vertices double" in lines
                @test "POLYGONS $cells $(5cells)" in lines
                @test "CELL_DATA $(cells+rim_edges)" in lines
                @test "SCALARS radial_gap_mm double 1" in lines
                @test "SCALARS analytic_area_mm2 double 1" in lines
                points = findfirst(line -> startswith(line,"POINTS "),lines)
                x,y,z = parse.(Float64,split(lines[points+1]))
                @test x ≈ 50.05+50tan(deg2rad(10.)) atol=1e-11
                @test (y,z) == (0.,0.)
                @test !any(line -> startswith(line,"SCALARS feed_cell "),lines)
                if rim_edges > 0
                    @test "LINES $rim_edges $(3rim_edges)" in lines
                    @test "SCALARS feed_rim int 1" in lines
                else
                    @test !any(line -> startswith(line,"LINES "),lines)
                end

                original_hash = sha256(read(path))
                process = run(pipeline(ignorestatus(command),stdout=output,stderr=output))
                @test !success(process)
                @test occursin("output already exists",String(take!(output)))
                @test sha256(read(path)) == original_hash
            end
        end
    end
end
