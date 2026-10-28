using Test, LinearAlgebra
include(joinpath(@__DIR__,"../verification/jfo/conical/cases.jl"))

jc(; angle=10.0,ecc=0.6) = BearingParameters(journal_midpoint_radius=50.0,
    bearing_length=100.0,radial_clearance=0.05,semicone_angle_deg=angle,eccentricity_ratio=ecc)
# An intentionally broad rectangular numerical reservoir, NOT the rig's port.
patch(nt,nz; centre=pi/2) = vec([abs(mod((i-0.5)*2pi/nt-centre+pi,2pi)-pi)<pi/4 &&
    0.375<(j-0.5)/nz<0.625 for i in 1:nt,j in 1:nz])
js(p=jc(),nt=32,nz=16; rpm=2000.0,mu=0.02752,ends=(0.0,0.0),pc=0.0,
   mask=patch(nt,nz),feed=5e5,kwargs...) = solve_conical_jfo(p,nt,nz;
    rpm,viscosity_pa_s=mu,end_pressure_pa=ends,cavitation_pressure_pa=pc,
    feed_cells=mask,feed_pressure_pa=feed,kwargs...)

@testset "Shared conservative 2D operator and seam" begin
    for rpm in (-2000.0,2000.0)
        nt,nz = 12,6
        a = assemble_conical_jfo(jc(),nt,nz;rpm,viscosity_pa_s=0.02752,
            end_pressure_pa=(2e5,1e5),cavitation_pressure_pa=0.0)
        P = [1e6*(1+sin(k)) for k in 1:nt*nz]
        F = [0.5+0.2*cos(k) for k in 1:nt*nz]
        q = conical_jfo_fluxes(a,P,F)
        @test q.balance ≈ a.A*P+a.B*(1 .- F)-a.b rtol=1e-12
        @test maximum(abs,vec(sum(a.B;dims=1))) == 0
        @test sum(q.balance) ≈ sum(q.qz[:,end])-sum(q.qz[:,1]) atol=1e-19
        for j in (1,nz)
            upstream = rpm > 0 ? nt : 1
            seam = a.gt[nt,j]*(P[nt*j]-P[1+nt*(j-1)])+
                a.adv[nt,j]*F[upstream+nt*(j-1)]
            @test q.qt[nt,j] ≈ seam rtol=1e-14
        end
        # Void transport is purely circumferential and preserves each row's inventory.
        for j in 1:nz
            @test abs(sum(reshape(a.B*(1 .- F),nt,nz)[:,j])) < 1e-20
        end
    end
end

@testset "Full-film limit and stationary cone" begin
    for angle in (0.0,10.0,30.0)
        p = jc(;angle)
        for rpm in (0.0,2000.0)
            r = js(p;rpm,ends=(3e7,3e7),mask=nothing,feed=nothing)
            old = solve_conical_full_film(p,32,16;rpm,viscosity_pa_s=0.02752,
                end_pressure_pa=(3e7,3e7))
            @test r.pressure_pa ≈ old.pressure_pa rtol=1e-11
            @test r.fill_fraction == ones(32*16)
            @test r.circumferential_face_flow_m3_s ≈ old.circumferential_face_flow_m3_s atol=1e-16
            @test r.axial_face_flow_m3_s ≈ old.axial_face_flow_m3_s atol=1e-16
        end
    end
    # Independent log-radius stationary centred-cone solution.
    p = jc(ecc=0.0)
    r = js(p,16,128;rpm=0.0,ends=(5e5,1e5),mask=nothing,feed=nothing)
    k = tand(10.0); r0 = 0.05005+0.05k; rL = 0.05005-0.05k
    exact = [5e5-4e5*log((r0-k*(j-0.5)*0.1/128)/r0)/log(rL/r0)
             for i in 1:16,j in 1:128]
    @test norm(r.pressure_pa-vec(exact))/norm(exact) < 1e-5
    @test r.fill_fraction == ones(16*128)
end

@testset "Curved 2D cavitation manufactured solution, both rotation signs" begin
    for rpm in (-2000.0,2000.0)
        results = [conical_jfo_mms(n,n÷2;rpm) for n in (32,64,128,256)]
        pe,fe = [r.pressure_error for r in results],[r.fill_error for r in results]
        @test all(pe[2:end] .< 0.3 .* pe[1:end-1])
        @test pe[end] < 1e-3
        @test all(fe[2:end] .< 0.6 .* fe[1:end-1])
        @test fe[end] < 1e-3
        for r in results
            @test r.minimum_pressure_pa >= 0
            @test 0 <= r.minimum_fill < r.maximum_fill <= 1
            @test r.complementarity == 0
            @test r.source_balance_error < 1e-10
        end
        println("2D MMS rpm=$rpm pressure errors=$pe; fill errors=$fe")
    end
end

@testset "Pressure reservoir source accounting, including seam and return flow" begin
    for centre in (0.0,pi/2)
        nt,nz = 64,32
        mask = patch(nt,nz;centre)
        r = js(jc(),nt,nz;mask)
        @test all(r.pressure_pa[mask] .== 5e5)
        @test all(r.fill_fraction[mask] .== 1)
        @test minimum(r.pressure_pa) == 0
        @test 0 < minimum(r.fill_fraction) < 1
        @test all(r.feed_source_m3_s[.!mask] .== 0)
        @test r.diagnostics.relative_max_cell_imbalance < 1e-10
        @test r.diagnostics.relative_global_imbalance < 1e-10
        @test sum(r.diagnostics.outward_end_flow_m3_s) ≈ r.diagnostics.net_feed_flow_m3_s rtol=1e-10
        @test minimum(-r.axial_face_flow_m3_s[:,1]) >= 0
        @test minimum(r.axial_face_flow_m3_s[:,end]) >= 0
        # Independently sum just the reservoir/film interface, not reservoir cells.
        m = reshape(mask,nt,nz)
        contour = 0.0
        for j in 1:nz,i in 1:nt
            next = i == nt ? 1 : i+1
            contour += (m[i,j]-m[next,j])*r.circumferential_face_flow_m3_s[i,j]
        end
        for j in 2:nz,i in 1:nt
            contour += (m[i,j-1]-m[i,j])*r.axial_face_flow_m3_s[i,j]
        end
        @test contour ≈ sum(r.feed_source_m3_s) rtol=1e-12
        if centre == 0
            @test all(m[1,13:20]) && all(m[end,13:20])
            @test r.diagnostics.gross_feed_return_m3_s > 0 # reservoir is not a check valve
            @test r.diagnostics.gross_feed_inflow_m3_s > r.diagnostics.gross_feed_return_m3_s
        end
    end
end

@testset "Reversal, pressure reference, scaling and cylindrical limit" begin
    nt,nz = 32,16
    r = js()
    shifted = js(;pc=101325.0,ends=(101325.0,101325.0),feed=601325.0)
    @test shifted.pressure_above_cavitation_pa ≈ r.pressure_pa rtol=1e-12
    @test shifted.fill_fraction ≈ r.fill_fraction rtol=1e-12
    @test shifted.feed_source_m3_s ≈ r.feed_source_m3_s rtol=1e-12
    reflected_mask = vec(reverse(reshape(patch(nt,nz),nt,nz);dims=1))
    reversed = js(;rpm=-2000.0,mask=reflected_mask)
    @test reshape(reversed.pressure_pa,nt,nz) ≈ reverse(reshape(r.pressure_pa,nt,nz);dims=1) rtol=1e-11
    @test reshape(reversed.fill_fraction,nt,nz) ≈ reverse(reshape(r.fill_fraction,nt,nz);dims=1) rtol=1e-11
    @test reversed.circumferential_face_flow_m3_s ≈ -circshift(reverse(r.circumferential_face_flow_m3_s;dims=1),(-1,0)) rtol=1e-11
    @test reversed.axial_face_flow_m3_s ≈ reverse(r.axial_face_flow_m3_s;dims=1) rtol=1e-11
    faster = js(;rpm=4000.0,feed=1e6)
    @test faster.pressure_pa ≈ 2r.pressure_pa rtol=1e-11
    @test faster.fill_fraction ≈ r.fill_fraction rtol=1e-11
    @test faster.feed_source_m3_s ≈ 2r.feed_source_m3_s rtol=1e-11
    viscous = js(;mu=0.05504,feed=1e6)
    @test viscous.pressure_pa ≈ 2r.pressure_pa rtol=1e-11
    @test viscous.fill_fraction ≈ r.fill_fraction rtol=1e-11
    @test viscous.feed_source_m3_s ≈ r.feed_source_m3_s rtol=1e-11
    cylinder,near = js(jc(angle=0.0)),js(jc(angle=1e-5))
    @test norm(cylinder.pressure_pa-near.pressure_pa)/norm(cylinder.pressure_pa) < 1e-6
    for ends in ((5e5,0.0),(0.0,5e5))
        endfed = js(;ends,mask=nothing,feed=nothing)
        @test endfed.diagnostics.cavitated_cells > 0
        @test endfed.diagnostics.relative_global_imbalance < 1e-10
        @test sum(endfed.diagnostics.outward_end_flow_m3_s) ≈ 0 atol=1e-15
    end
end

@testset "Unforced patch-fed cone refinement" begin
    grids = (32,64,128,256)
    runs = [conical_jfo_numerical_case(n,n÷2) for n in grids]
    errors = [conical_jfo_nested_error(runs[k-1],runs[k],grids[k-1],grids[k-1]÷2)
              for k in 2:4]
    pe = [e.pressure_relative_l2 for e in errors]
    fe = [e.fill_area_mean_absolute for e in errors]
    @test all(pe[2:end] .< 0.8 .* pe[1:end-1])
    @test pe[end] < 0.02
    @test fe[end] < 0.02
    @test fe[end] < fe[1]
    flows = [r.diagnostics.net_feed_flow_m3_s for r in runs]
    @test abs(flows[end]/flows[end-1]-1) < 0.03
    for r in runs
        @test minimum(r.pressure_pa) >= 0
        @test 0 <= minimum(r.fill_fraction) < 1
        @test r.diagnostics.relative_global_imbalance < 1e-10
    end
    println("Unforced patch refinement: pressure mismatches=$pe; fill mismatches=$fe; feed flows=$flows")
end

@testset "No silent inventory, invalid supply, and export guards" begin
    @test_throws ArgumentError js(;mask=nothing,feed=nothing) # all reservoirs at pcav
    @test_throws ArgumentError js(;rpm=0.0,mask=nothing,feed=nothing)
    @test_throws ArgumentError js(;ends=(-1.0,0.0))
    @test_throws ArgumentError js(;pc=NaN)
    @test_throws ArgumentError js(;mask=nothing)
    @test_throws ArgumentError js(;feed=nothing)
    @test_throws ArgumentError js(;mask=falses(512))
    @test_throws ArgumentError js(;mask=trues(512))
    @test_throws ArgumentError js(;mask=trues(3))
    @test_throws ArgumentError js(;feed=0.0)
    @test_throws ArgumentError js(;feed=Inf)
    @test_throws ArgumentError js(;mu=0.0)
    @test_throws ArgumentError js(;rpm=NaN)
    @test_throws ArgumentError js(;tolerance=0.0)
    @test_throws ArgumentError js(;max_iterations=0)
    @test_throws ErrorException js(;max_iterations=1)
    mktempdir() do dir
        r = js()
        path = write_conical_jfo_vtk(joinpath(dir,"cavitation.vtk"),r)
        data = read(path,String)
        for name in ("fill_fraction","pressure_above_cavitation_pa","feed_source_m3_s","feed_cell")
            @test occursin("SCALARS $name",data)
        end
        @test_throws ErrorException write_conical_jfo_vtk(path,r)
    end
end
