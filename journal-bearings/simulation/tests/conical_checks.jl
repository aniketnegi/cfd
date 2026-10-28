using Test, LinearAlgebra
include(joinpath(@__DIR__,"../reynolds/src/conical.jl"))

cone(; angle=10.0, ecc=0.6) = BearingParameters(
    journal_midpoint_radius=50.0, bearing_length=100.0,
    radial_clearance=0.05, semicone_angle_deg=angle, eccentricity_ratio=ecc)
const CONICAL_MU = 32e-6*860.0 # illustrative VG32 reference, not measured rig oil
cs(p, nt=32, nz=16; rpm=2000.0, mu=CONICAL_MU, ends=(0.0,0.0)) =
    solve_conical_full_film(p, nt, nz; rpm, viscosity_pa_s=mu, end_pressure_pa=ends)

# Compare at nested cell centres using the mean of four fine pressures.
function conical_pressure_mismatch(coarse, fine, nt, nz)
    pc, pf = reshape(coarse.pressure_pa,nt,nz), reshape(fine.pressure_pa,2nt,2nz)
    # Area-weight the fine cells: housing radius changes along z.
    af = reshape(fine.grid.areas,2nt,2nz)
    avg = [sum(pf[2i-1:2i,2j-1:2j] .* af[2i-1:2i,2j-1:2j]) /
           sum(af[2i-1:2i,2j-1:2j]) for i in 1:nt, j in 1:nz]
    return norm(pc-avg)/norm(avg)
end

@testset "Conical geometry and zero-angle regression" begin
    p = cone()
    for z in (0.0,0.05,0.1), theta in (0.0,0.7,pi,4.0)
        xh,yh,_ = bearing_point(p,theta,z*1e3,1.0)
        xj,yj,_ = bearing_point(p,theta,z*1e3,0.0)
        @test conical_normal_gap(p,theta,z) ≈
              (hypot(xh,yh)-hypot(xj,yj))*1e-3*cosd(10.0) rtol=1e-11
        # Independent normal-ray intersection of the analytic (extended) cones.
        # Quantify the projection approximation, rather than calling it exact.
        rj = 0.05+(0.05-z)*tand(10.0)
        rh, e = rj+50e-6,30e-6
        lo,hi = 0.0,100e-6
        for _ in 1:60
            h = (lo+hi)/2
            radius = rh-h*cosd(10.0)
            distance = hypot(radius*cos(theta)-e,radius*sin(theta)) -
                       (rj+h*sind(10.0)*tand(10.0))
            if distance > 0
                lo = h
            else
                hi = h
            end
        end
        @test conical_normal_gap(p,theta,z) ≈ (lo+hi)/2 rtol=1e-7
    end
    for z in (0.0,0.05,0.1), theta in (0.0,1.0,3.0)
        @test conical_normal_gap(cone(ecc=0.0),theta,z) ≈ 50e-6*cosd(10.0) rtol=1e-14
    end
    p0 = cone(angle=0.0)
    old = solve_cylindrical_full_film(p0,32,16;rpm=2000.0,viscosity_pa_s=CONICAL_MU)
    new = cs(p0)
    @test new.pressure_pa ≈ old.pressure_pa rtol=2e-12
    @test new.grid.nodes == old.grid.nodes
    @test new.grid.quads == old.grid.quads
    @test new.normal_gap_m ≈ new.grid.radial_gaps .* 1e-3 rtol=1e-11
    @test abs(new.housing_pressure_force_n[3]) == 0
    near = cs(cone(angle=1e-5))
    @test norm(near.pressure_pa-new.pressure_pa)/norm(new.pressure_pa) < 1e-6
end

@testset "Conical conservation, linearity, pressure reference and export" begin
    p = cone()
    r = cs(p)
    @test minimum(r.pressure_pa) < 0 < maximum(r.pressure_pa)
    @test r.diagnostics.relative_l1_cell_imbalance < 1e-11
    @test r.diagnostics.relative_max_cell_imbalance < 1e-10
    @test r.diagnostics.relative_global_imbalance < 1e-10
    @test cs(p;rpm=-2000.0).pressure_pa ≈ -r.pressure_pa rtol=1e-11
    @test cs(p;mu=2CONICAL_MU).pressure_pa ≈ 2r.pressure_pa rtol=1e-11
    @test maximum(abs,cs(cone(ecc=0.0)).pressure_pa) < 1e-8
    @test maximum(abs,cs(p;rpm=0.0).pressure_pa) == 0
    @test cs(p;rpm=0.0,ends=(123456.0,123456.0)).pressure_pa ≈
          fill(123456.0,32*16) rtol=2e-12
    shifted = cs(p;ends=(2e5,2e5))
    @test shifted.pressure_pa ≈ r.pressure_pa .+ 2e5 rtol=1e-11
    @test shifted.axial_face_flow_m3_s ≈ r.axial_face_flow_m3_s rtol=1e-10
    @test shifted.circumferential_face_flow_m3_s ≈ r.circumferential_face_flow_m3_s rtol=1e-11
    # A cone is not axially symmetric: reversing z does not reproduce the field.
    field = reshape(r.pressure_pa,32,16)
    @test norm(field-reverse(field;dims=2))/norm(field) > 0.01
    @test_throws ArgumentError cs(p,3,2)
    @test_throws ArgumentError cs(p,8,1)
    @test_throws ArgumentError cs(p;mu=0.0)
    @test_throws ArgumentError cs(p;mu=NaN)
    @test_throws ArgumentError cs(p;rpm=Inf)
    @test_throws ArgumentError cs(p;ends=(NaN,0.0))
    mktempdir() do dir
        path = write_conical_pressure_vtk(joinpath(dir,"full-film.vtk"),r)
        output = read(path,String)
        @test occursin("SCALARS pressure_gauge_pa double 1",output)
        values = split(output,"SCALARS normal_gap_m double 1\nLOOKUP_TABLE default\n")
        @test length(values) == 2
        @test parse.(Float64,split(strip(values[2]))) == r.normal_gap_m
    end
end

@testset "Independent centred-cone pressure-driven analytical solution" begin
    mu, p0, pL = CONICAL_MU, 5e5, 0.0
    for angle in (0.0,10.0)
        p = cone(;angle,ecc=0.0)
        L, r0 = 0.1, (50.05+50*tand(angle))*1e-3
        k, ca = tand(angle), cosd(angle)
        integral(z) = iszero(k) ? z/r0 : -log1p(-k*z/r0)/k
        h = 50e-6*ca
        exact_flow = 2pi*h^3*ca/(12mu)*(p0-pL)/integral(L)
        errors, flow_errors = Float64[],Float64[]
        for nz in (8,16,32,64,128)
            r = cs(p,8,nz;rpm=0.0,ends=(p0,pL))
            exact = [p0+(pL-p0)*integral((j-0.5)*L/nz)/integral(L)
                     for i in 1:8,j in 1:nz]
            err = norm(r.pressure_pa-vec(exact))/norm(exact)
            q = vec(sum(r.axial_face_flow_m3_s;dims=1))
            push!(errors,err)
            push!(flow_errors,maximum(abs.(q.-exact_flow))/abs(exact_flow))
            @test maximum(q)-minimum(q) < 1e-10*exact_flow
            @test r.diagnostics.relative_global_imbalance < 1e-10
            @test r.pressure_pa ≈ vec(repeat(reshape(r.pressure_pa,8,nz)[1:1,:],8,1)) rtol=1e-12
            @test r.diagnostics.outward_end_flow_m3_s[1] < 0 < r.diagnostics.outward_end_flow_m3_s[2]
        end
        if iszero(angle)
            @test maximum(errors) < 1e-12
            @test maximum(flow_errors) < 1e-11
        else
            @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
            @test all(flow_errors[2:end] .< 0.3 .* flow_errors[1:end-1])
            @test errors[end] < 1e-5
            @test flow_errors[end] < 1e-5
        end
        println("Axisymmetric angle=$angle: pressure errors=$errors; flow errors=$flow_errors")
    end
end

@testset "Two-dimensional manufactured conical pressure" begin
    # Constant normal gap, but variable radius and both metric coefficients.
    # Exact p=P*sin(theta)*z*(L-z)/L², with independently differentiated source.
    # S = -∂theta(K/(r*cosα)*∂theta p)-∂z(K*r*cosα*∂z p), K=h³/(12μ).
    p = cone(ecc=0.0)
    L, r0, k, ca, P = 0.1,(50.05+50*tand(10.0))*1e-3,tand(10.0),cosd(10.0),1e5
    K = (50e-6*ca)^3/(12CONICAL_MU)
    errors = Float64[]
    for (nt,nz) in ((16,8),(32,16),(64,32),(128,64))
        a = assemble_conical_full_film(p,nt,nz;rpm=0.0,viscosity_pa_s=CONICAL_MU)
        dt,dz = 2pi/nt,L/nz
        source = [begin
            z = (j-0.5)*dz
            r = r0-k*z
            f, df = z*(L-z)/L^2,(L-2z)/L^2
            # Exact angular integral; midpoint integration along z.
            K*P*(cos((i-1)*dt)-cos(i*dt))*dz *
                (f/(r*ca)+ca*(2r/L^2+k*df))
        end for i in 1:nt,j in 1:nz]
        exact = [P*sin((i-0.5)*dt)*((j-0.5)*dz)*(L-(j-0.5)*dz)/L^2
                 for i in 1:nt,j in 1:nz]
        computed = a.A \ (a.b+vec(source))
        push!(errors,norm(computed-vec(exact))/norm(exact))
    end
    @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
    @test errors[end] < 5e-4
    println("Manufactured pressure relative errors: $errors")
end

@testset "Rotating eccentric cone manufactured solution" begin
    # The centred MMS above has no wedge source and constant h.
    # Exercise nonzero Couette forcing AND variable h^3 in both coordinates.
    # All derivatives below are analytic, independent of the assembly stencil.
    p = cone()
    L, c, e, ca, k, P = 0.1,50e-6,30e-6,cosd(10.0),tand(10.0),1e7
    for rpm in (-2000.0,2000.0)
        omega = rpm*2pi/60
        errors = Float64[]
        for (nt,nz) in ((32,16),(64,32),(128,64),(256,128))
            a = assemble_conical_full_film(p,nt,nz;rpm,viscosity_pa_s=CONICAL_MU)
            dt,dz = 2pi/nt,L/nz
            exact,source = zeros(nt,nz),zeros(nt,nz)
            for j in 1:nz,i in 1:nt
                theta,z = (i-0.5)*dt,(j-0.5)*dz
                st,ct = sincos(theta)
                rj = 0.05+(L/2-z)*k
                r = rj+c
                D = sqrt(rj^2-e^2*st^2)
                h = ca*(c-e*ct+e^2*st^2/(rj+D))
                ht = ca*(e*st+e^2*st*ct/D)
                hz = ca*k*(rj/D-1)
                K,Kt,Kz = h^3/(12CONICAL_MU),h^2*ht/(4CONICAL_MU),h^2*hz/(4CONICAL_MU)
                f,df = z*(L-z)/L^2,(L-2z)/L^2
                pressure = exact[i,j] = P*st*f
                pt,ptt,pz,pzz = P*ct*f,-pressure,P*st*df,-2P*st/L^2
                diffusion = (Kt*pt+K*ptt)/(r*ca) + ca*((Kz*r-K*k)*pz+K*r*pzz)
                wedge = omega*rj*ht/(2ca)
                # A approximates -diffusion, b approximates -wedge.
                # Artificial injection = wedge-div(K grad p); midpoint integration.
                source[i,j] = (wedge-diffusion)*dt*dz
            end
            computed = a.A \ (a.b+vec(source))
            push!(errors,norm(computed-vec(exact))/norm(exact))
        end
        @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
        @test errors[end] < 1e-3
        println("Rotating eccentric MMS rpm=$rpm errors=$errors")
    end
end

@testset "Diffusion operator and end-boundary invariants" begin
    p = cone()
    a = assemble_conical_full_film(p,16,8;rpm=0.0,viscosity_pa_s=CONICAL_MU)
    @test issymmetric(a.A)
    @test isposdef(Symmetric(Matrix(a.A)))
    @test all(x -> x <= 0, nonzeros(a.A-spdiagm(0=>diag(a.A))))
    pressure = cs(p;rpm=0.0,ends=(2e5,5e5))
    @test minimum(pressure.pressure_pa) >= 2e5
    @test maximum(pressure.pressure_pa) <= 5e5
    @test pressure.diagnostics.outward_end_flow_m3_s[1] > 0
    @test pressure.diagnostics.outward_end_flow_m3_s[2] < 0
    reverse = cs(p;rpm=0.0,ends=(5e5,2e5))
    @test reverse.pressure_pa ≈ 7e5 .- pressure.pressure_pa rtol=1e-11
    @test reverse.axial_face_flow_m3_s ≈ -pressure.axial_face_flow_m3_s rtol=1e-11
end

@testset "Oriented pressure-force integration" begin
    p = cone()
    nt,nz,P = 128,32,2e5
    f = conical_housing_pressure_force(p,fill(P,nt*nz),nt,nz)
    large,small = (50.05+50*tand(10.0))*1e-3,(50.05-50*tand(10.0))*1e-3
    @test abs(f[1]) < 1e-9
    @test abs(f[2]) < 1e-9
    @test f[3] ≈ P*pi*(large^2-small^2) rtol=1e-12
    # Independent continuum integrals for p=P*cos(theta) and P*sin(theta).
    for (component,fn) in ((1,cos),(2,sin))
        pressure = [P*fn((i-0.5)*2pi/nt) for j in 1:nz for i in 1:nt]
        f = conical_housing_pressure_force(p,pressure,nt,nz)
        @test f[component] ≈ P*pi*0.05005*0.1 rtol=2e-4
        @test abs(f[3]) < 1e-9
        @test abs(f[3-component]) < 1e-9
    end
    @test_throws ArgumentError conical_housing_pressure_force(p,[1.0],nt,nz)
end

@testset "Eccentric rotating cone refinement" begin
    p = cone()
    runs = [cs(p,nt,nz) for (nt,nz) in ((32,16),(64,32),(128,64))]
    errors = [conical_pressure_mismatch(runs[j],runs[j+1],32*2^(j-1),16*2^(j-1)) for j in 1:2]
    @test 0 < errors[2] < 0.3errors[1]
    @test errors[2] < 0.005
    force_changes = [norm(collect(runs[j+1].housing_pressure_force_n).-collect(runs[j].housing_pressure_force_n)) /
                     norm(collect(runs[j+1].housing_pressure_force_n)) for j in 1:2]
    @test force_changes[2] < 0.3force_changes[1]
    @test force_changes[2] < 0.005
    println("Eccentric pressure mismatches: $errors; force changes: $force_changes")
end
