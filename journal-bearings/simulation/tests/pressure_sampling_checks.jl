using Test, LinearAlgebra
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../reynolds/src/full_film.jl"))
include(joinpath(@__DIR__,"../reynolds/src/pressure_sampling.jl"))

@testset "Actual-plane fitted pressure reconstruction" begin
    for angle in (0.,10.)
        errors=Float64[]
        for nt in (128,256)
            p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
                semicone_angle_deg=angle,eccentricity_ratio=.6)
            grid=surface_feed_grid(p,nt,nt÷2;feed_diameter=5.,feed_z=50.)
            s=grid.surface
            constant=fitted_pressure_slice(grid,fill(2e5,length(grid.quads));z_m=.05,
                boundary_pressure=(x,y,tag)->2e5)
            @test all(value->all(p->abs(p/2e5-1)<1e-10,value),constant.pressure)
            @test constant.maximum_pressure_pa≈2e5 rtol=1e-10
            coverage=sum(pair[2]-pair[1] for pair in constant.theta)
            @test coverage≈2pi-2asin(2.5/50.05) rtol=1e-10
            affine=1e5 .+ 3e6.*s.centre[:,2]
            for z in (.025,.05,.075)
                slice=fitted_pressure_slice(grid,affine;z_m=z,boundary_pressure=(x,y,tag)->1e5+3e6*y)
                target=1e5+3e6*z/s.ca
                @test all(pair->all(p->abs(p/target-1)<1e-10,pair),slice.pressure)
            end
            exact(theta,z)=1e6*(2+cos(theta))*4z*(.1-z)/.1^2
            values=[exact(s.centre[i,1]/s.R,s.centre[i,2]*s.ca) for i in axes(s.centre,1)]
            slice=fitted_pressure_slice(grid,values;z_m=.05,
                boundary_pressure=(x,y,tag)->exact(x/s.R,y*s.ca))
            error=maximum(abs(slice.pressure[i][j]-exact(slice.theta[i][j],.05)) for i in eachindex(slice.cells),j in 1:2)/3e6
            push!(errors,error)
            @test abs(slice.maximum_pressure_pa/3e6-1)<.002
            @test abs(atan(sin(slice.peak_theta_rad),cos(slice.peak_theta_rad)))<2pi/nt
            @test length(slice.boundary_values)>=2
            @test_throws ArgumentError fitted_pressure_slice(grid,values;z_m=0.,boundary_pressure=(x,y,tag)->0.)
            @test_throws ArgumentError fitted_pressure_slice(grid,[1.];z_m=.05,boundary_pressure=(x,y,tag)->0.)
            println("Fitted slice angle=$angle nt=$nt endpoint_relative_Linf=$error")
        end
        @test errors[2]<.4errors[1]
        @test errors[2]<.0005
    end
end
