include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"sampling.jl"))
length(ARGS) in (1,2) && (length(ARGS)==1 || ARGS[2]=="--full") || error("usage: run_conical.jl FRESH_OUTPUT [--full]")
const FULL="--full" in ARGS
const OUTPUT=fresh_output(ARGS[1];suite=FULL ? "conical-refinement-full" : "conical-refinement-representative")
function main()
    ispath(joinpath(OUTPUT,"results.tsv")) && error("refusing overwrite")
    families=[(10.,ecc,groove) for ecc in (.2,.6,.8) for groove in (false,true)]
    push!(families,(20.,.6,false))
    FULL || filter!(item->item[1]==10. && item[2]==.6,families)
    grids=FULL ? (64,128,256,512) : (64,128,256)
    @testset "Conical Swift-Stieber physical-equation refinement" begin
        open(joinpath(OUTPUT,"results.tsv"),"w") do io
            first=true
            for (angle,ecc,groove) in families
                previous=nothing;field_changes=Float64[];force_changes=Float64[];flow_changes=Float64[]
                for nt in grids
                    nz=nt÷2;p=cone(;angle,ecc)
                    elapsed=@elapsed result=solve_case(p,nt,nz;groove)
                    check_result(result)
                    m=midplane(result.a,result.pressure)
                    @test m.lower<=m.target<=m.upper
                    field_change=force_change=flow_change=NaN
                    if !isnothing(previous)
                        field_change=coordinate_mismatch(previous,result)
                        force_change=norm(collect(result.force)-collect(previous.force))/norm(collect(result.force))
                        flow_change=abs(sum(result.ends)-sum(previous.ends))/abs(sum(result.ends))
                        push!(field_changes,field_change);push!(force_changes,force_change);push!(flow_changes,flow_change)
                    end
                    peak_index=argmax(m.profile)
                    row=(;angle_deg=angle,epsilon=ecc,groove,nt,nz,
                        maximum_pressure_pa=maximum(result.pressure),midplane_maximum_pa=maximum(m.profile),
                        midplane_sample_z_m=m.target,midplane_lower_z_m=m.lower,midplane_upper_z_m=m.upper,
                        midplane_peak_theta_rad=result.a.axes.theta_centres[peak_index],
                        pressure_common_coordinate_mismatch=field_change,force_relative_change=force_change,
                        total_end_flow_relative_change=flow_change,seconds=elapsed,result.stats...)
                    row!(io,row;header=first);first=false
                    @printf("angle=%g eps=%g groove=%s %dx%d mid=%.7fMPa global=%.7fMPa residual=%.3g\n",
                        angle,ecc,groove,nt,nz,maximum(m.profile)/1e6,maximum(result.pressure)/1e6,result.stats.positive_pressure_residual)
                    flush(stdout)
                    if nt==last(grids)
                        save_field("angle$(angle)-e$(ecc)-groove$(groove)",result)
                        mkpath(joinpath(OUTPUT,"profiles"))
                        open(joinpath(OUTPUT,"profiles","angle$(angle)-e$(ecc)-groove$(groove).tsv"),"w") do profile
                            println(profile,"theta_rad\tmidplane_pressure_pa\tlower_outward_flux_m3_s\tupper_outward_flux_m3_s")
                            for i in 1:nt
                                println(profile,join((result.a.axes.theta_centres[i],m.profile[i],
                                    -result.flux.qz[i,1],result.flux.qz[i,end]),'\t'))
                            end
                        end
                    end
                    previous=result
                    GC.gc()
                end
                @test all(field_changes[2:end] .< 0.4 .* field_changes[1:end-1])
                @test all(force_changes[2:end] .< 0.4 .* force_changes[1:end-1])
                @test field_changes[end]<.001
                @test force_changes[end]<.001
                @test flow_changes[end]<.001
            end
        end
    end
end
main()
