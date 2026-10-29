using Test, LinearAlgebra, SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../reynolds/src/swift_stieber.jl"))
include(joinpath(@__DIR__,"../verification/swift-stieber/fitted_reference.jl"))

function enumerate_obstacle(A,b)
    solutions=Vector{Float64}[];n=length(b)
    for mask in 0:2^n-1
        free=findall(i->!iszero(mask & (1<<(i-1))),1:n)
        p=zeros(n)
        isempty(free) || (p[free]=A[free,free]\b[free])
        reaction=A*p-b
        if minimum(p)>=-1e-12 && minimum(reaction)>=-1e-12 && maximum(abs,p.*reaction)<1e-12
            push!(solutions,p)
        end
    end
    isempty(solutions) && error("enumeration found no obstacle solution")
    @test all(p->norm(p-solutions[1],Inf)<1e-12,solutions)
    return first(solutions)
end

@testset "Nonsymmetric pressure obstacle against exhaustive algebraic reference" begin
    A=[4. 1. 0. 0.;-1. 3. .5 0.;.2 -.7 2.5 .3;0. .1 -.6 3.]
    @test !issymmetric(A)
    @test isposdef(Symmetric((A+A')/2))
    @test any(>(0),nonzeros(sparse(A-Diagonal(diag(A)))))
    for b in ([1.,-2.,.5,-.4],[-1.,2.,-.5,.4],ones(4),-ones(4))
        exact=enumerate_obstacle(A,b)
        for pscale in (.3,10.)
            result=pressure_obstacle(sparse(A),b,pscale,.7)
            @test result.pressure≈exact atol=1e-12
            @test minimum(result.reaction)>-1e-12
        end
        reference=fitted_projected_sor(sparse(A),b)
        @test reference.pressure≈exact atol=1e-11
    end
end

parameters(angle)=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
    radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=.6)

@testset "Fitted pressure-only cavitation contract and fully-liquid recovery" begin
    @test !isdefined(@__MODULE__,:write_vtk)
    for angle in (0.,10.)
        grid=surface_feed_grid(parameters(angle),128,64;feed_diameter=5.,feed_z=50.)
        original=deepcopy(grid)
        kw=(;rpm=2000.,viscosity_pa_s=.0277,feed_pressure_pa=3e5,eccentricity_angle_rad=.37)
        result=solve_swift_stieber(grid;kw...)
        @test grid==original
        @test result.grid.nodes===grid.nodes
        @test result.grid.surface===grid.surface
        @test minimum(result.pressure_pa)>=0
        @test any(iszero,result.pressure_pa)
        @test maximum(result.pressure_pa)>3e5
        @test !hasproperty(result,:axial_face_flow_m3_s)
        @test result.grid.radial_gaps≈result.normal_gap_m.*(1e3/cosd(angle)) rtol=1e-13
        a=assemble_full_film(grid;kw...)
        @test !issymmetric(a.A)
        @test isposdef(Symmetric((a.A+a.A')/2))
        @test any(>(0),nonzeros(a.A-spdiagm(0=>diag(a.A))))
        reaction=a.A*result.pressure_pa-a.b
        @test reaction==result.reaction_m3_s
        @test minimum(reaction)/result.diagnostics.flow_scale_m3_s > -1e-10
        @test norm(a.faces.divergence*result.face_flow_m3_s-reaction,Inf)/
            result.diagnostics.flow_scale_m3_s<1e-10
        @test result.diagnostics.positive_pressure_relative_residual<1e-10
        @test result.diagnostics.complementarity_relative<1e-10
        @test result.diagnostics.flux_reconstruction_relative<1e-10
        @test result.diagnostics.defect_accounting_relative<1e-10
        @test result.diagnostics.extrapolated_outward_feed_flow_m3_s==result.boundary_outflow_m3_s.feed
        @test collect(result.diagnostics.extrapolated_outward_end_flow_m3_s)==
            [result.boundary_outflow_m3_s.large,result.boundary_outflow_m3_s.small]
        @test abs(sum(values(result.boundary_outflow_m3_s))-sum(reaction))<1e-10*
            max(sum(abs,values(result.boundary_outflow_m3_s)),result.diagnostics.flow_scale_m3_s)
        signed=solve_full_film(grid;kw...)
        @test norm(result.pressure_pa-max.(signed.pressure_pa,0.))/norm(result.pressure_pa)>.01
        shifted=solve_swift_stieber(grid;kw...,cavitation_pressure_pa=101325.,
            end_pressure_pa=(101325.,101325.),feed_pressure_pa=401325.)
        @test shifted.pressure_above_cavitation_pa==result.pressure_pa
        @test shifted.face_flow_m3_s==result.face_flow_m3_s
        @test shifted.pressure_pa==result.pressure_pa .+ 101325.
        doubled=solve_swift_stieber(grid;kw...,viscosity_pa_s=.0554,feed_pressure_pa=6e5)
        @test doubled.pressure_pa≈2result.pressure_pa rtol=1e-11
        @test doubled.face_flow_m3_s≈result.face_flow_m3_s rtol=1e-11
        flooded=solve_swift_stieber(grid;kw...,end_pressure_pa=(3e7,3e7),feed_pressure_pa=3e7)
        full=solve_full_film(grid;kw...,end_pressure_pa=(3e7,3e7),feed_pressure_pa=3e7)
        @test minimum(full.pressure_pa)>0
        @test flooded.pressure_pa≈full.pressure_pa rtol=1e-11
        @test flooded.face_flow_m3_s≈full.face_flow_m3_s rtol=1e-11
        idle=solve_swift_stieber(grid;kw...,rpm=0.,feed_pressure_pa=0.)
        @test all(iszero,idle.pressure_pa)
        @test all(iszero,idle.face_flow_m3_s)
        pins=[1,5,19]
        pinned=solve_swift_stieber(grid;kw...,pinned_cells=pins)
        @test all(iszero,pinned.pressure_pa[pins])
        @test pinned.diagnostics.defect_accounting_relative<1e-10
        constant=solve_swift_stieber(grid;kw...,rpm=0.,end_pressure_pa=(2e5,2e5),feed_pressure_pa=2e5)
        @test constant.pressure_pa≈fill(2e5,length(grid.quads)) rtol=1e-10
        @test maximum(abs,constant.face_flow_m3_s)<1e-16
        reused=solve_swift_stieber(result.grid;kw...)
        @test reused.pressure_pa==result.pressure_pa
        @test reused.grid===result.grid
        @test_throws ArgumentError solve_swift_stieber(grid;rpm=2000.,viscosity_pa_s=.0277)
        @test_throws ArgumentError solve_swift_stieber(grid;kw...,feed_pressure_pa=-1.)
        @test_throws ArgumentError solve_swift_stieber(grid;kw...,feed_pressure_pa=NaN)
        @test_throws ArgumentError solve_swift_stieber(grid;kw...,cavitation_pressure_pa=4e5,end_pressure_pa=(4e5,4e5))
        @test_throws ArgumentError solve_swift_stieber(grid;kw...,pinned_cells=[length(grid.quads)+1])
        @test_throws ArgumentError solve_swift_stieber(grid;kw...,max_iterations=0)
        @test_throws ErrorException solve_swift_stieber(grid;kw...,max_iterations=1)
    end
end
