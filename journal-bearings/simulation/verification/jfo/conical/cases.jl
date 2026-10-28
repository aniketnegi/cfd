# Artificial-reservoir fixture and refinement comparison for JFO checks.
include(joinpath(@__DIR__, "../../../reynolds/src/conical_jfo.jl"))
include(joinpath(@__DIR__, "conical_jfo_mms.jl"))

function conical_jfo_numerical_case(nt,nz)
    p = BearingParameters(journal_midpoint_radius=50.0,bearing_length=100.0,
        radial_clearance=0.05,semicone_angle_deg=10.0,eccentricity_ratio=0.6)
    mask = vec([pi/4<(i-0.5)*2pi/nt<3pi/4 && 0.375<(j-0.5)/nz<0.625
                for i in 1:nt,j in 1:nz])
    return solve_conical_jfo(p,nt,nz;rpm=2000.0,viscosity_pa_s=0.02752,
        end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0,
        feed_cells=mask,feed_pressure_pa=5e5)
end

function conical_jfo_nested_error(coarse,fine,nt,nz)
    weights = reshape(fine.grid.areas,2nt,2nz)
    function average(field)
        f = reshape(field,2nt,2nz)
        [sum(f[2i-1:2i,2j-1:2j].*weights[2i-1:2i,2j-1:2j]) /
         sum(weights[2i-1:2i,2j-1:2j]) for i in 1:nt,j in 1:nz]
    end
    p,f = vec(average(fine.pressure_pa)),vec(average(fine.fill_fraction))
    return (; pressure_relative_l2=norm(coarse.pressure_pa-p)/norm(p),
        fill_area_mean_absolute=sum(abs.(coarse.fill_fraction-f).*coarse.grid.areas)/sum(coarse.grid.areas))
end
