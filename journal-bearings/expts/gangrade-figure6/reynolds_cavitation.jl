# Count-based adapter for case-specific Swift-Stieber comparisons.
# No Gangrade inlet or start-line alignment is selected implicitly.
isdefined(@__MODULE__, :BearingParameters) || include(joinpath(@__DIR__,"../../bearing.jl"))
isdefined(@__MODULE__, :housing_grid) || include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
isdefined(@__MODULE__, :solve_swift_stieber) || include(joinpath(@__DIR__,"../../simulation/reynolds/src/swift_stieber.jl"))

function fig6_reynolds(p,nt,nz;rpm=2000.0,viscosity_pa_s=0.0277,
        end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0,
        pinned_column::Union{Nothing,Int}=nothing,max_iterations=2000)
    isnothing(pinned_column) || 1<=pinned_column<=nt || throw(ArgumentError("invalid pinned column"))
    pins=isnothing(pinned_column) ? Int[] : collect(pinned_column:nt:nt*nz)
    result=solve_swift_stieber(housing_grid(p,nt,nz);rpm,viscosity_pa_s,
        end_pressure_pa,cavitation_pressure_pa,pinned_cells=pins,max_iterations)
    d=result.diagnostics
    return (;pressure_pa=result.pressure_pa,
        pressure_above_cavitation_pa=result.pressure_above_cavitation_pa,
        reaction_m3_s=result.reaction_m3_s,iterations=d.iterations,
        positive_pressure_relative_residual=d.positive_pressure_relative_residual,
        complementarity_relative=d.complementarity_relative,
        cavity_full_film_flux_defect_m3_s=d.cavity_full_film_flux_defect_m3_s,
        pinned_line_exchange_m3_s=d.pinned_line_exchange_m3_s,
        extrapolated_end_outflow_m3_s=d.extrapolated_end_outflow_m3_s)
end
