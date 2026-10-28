# Pressure exporters; construction and solution do not load this file.
include(joinpath(@__DIR__,"../../../mesh/vtk.jl"))

function write_pressure_vtk(path::AbstractString, result)
    write_vtk(path, result.grid)
    open(path, "a") do io
        println(io, "SCALARS pressure_gauge_pa double 1\nLOOKUP_TABLE default")
        for value in result.pressure_pa
            println(io, value)
        end
    end
    return path
end

function write_conical_pressure_vtk(path::AbstractString, result)
    write_pressure_vtk(path, result)
    open(path, "a") do io
        println(io, "SCALARS normal_gap_m double 1\nLOOKUP_TABLE default")
        for gap in result.normal_gap_m
            println(io, gap)
        end
    end
    return path
end
