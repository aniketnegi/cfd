using Test, LinearAlgebra, SparseArrays
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../reynolds/src/full_film.jl"))

parameters(; angle=10., ecc=.6, radius=50., length=100.) = BearingParameters(
    journal_midpoint_radius=radius, bearing_length=length, radial_clearance=.05,
    semicone_angle_deg=angle, eccentricity_ratio=ecc)

# Smooth coordinate maps give related graded meshes under refinement.
function graded_grid(p, nt, nz)
    theta_nodes = [.17*2pi/nt+2pi*t+.25sin(2pi*t) for t in (0:nt-1)./nt]
    z_nodes = [p.bearing_length*(t+.15sin(2pi*t)/(2pi)) for t in (0:nz)./nz]
    z_nodes[1], z_nodes[end] = 0., p.bearing_length
    return housing_grid(p,nt,nz; theta_nodes, z_nodes)
end

@testset "Full-film dependencies and supplied tensor geometry" begin
    @test !isdefined(@__MODULE__, :write_vtk)
    @test !isdefined(@__MODULE__, :solve_jfo_active_set)
    @test all(id.name != "Gmsh" for id in keys(Base.loaded_modules))
    for angle in (0.,10.)
        p = parameters(;angle)
        grid = graded_grid(p,32,16)
        original = deepcopy(grid)
        a = assemble_full_film(grid;rpm=2000.,viscosity_pa_s=.0277)
        r = solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277)
        @test grid == original
        @test r.grid === grid
        @test grid.parameters === p
        @test a.axes.theta_edges[1:end-1] == grid.theta_nodes
        @test a.axes.z_edges == grid.z_nodes.*1e-3
        @test length(r.pressure_pa) == length(grid.quads)
        uniform = solve_full_film(housing_grid(p,32,16);rpm=2000.,viscosity_pa_s=.0277)
        @test norm(r.pressure_pa-uniform.pressure_pa)/norm(uniform.pressure_pa) > .05
        @test r.diagnostics.relative_l1_cell_imbalance < 1e-11
        @test r.diagnostics.relative_global_imbalance < 1e-10
        @test solve_full_film(grid;rpm=-2000.,viscosity_pa_s=.0277).pressure_pa ≈ -r.pressure_pa rtol=1e-11
        @test solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0554).pressure_pa ≈ 2r.pressure_pa rtol=1e-11
        shifted = solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(2e5,2e5))
        @test shifted.pressure_pa ≈ r.pressure_pa .+ 2e5 rtol=1e-11
        @test shifted.axial_face_flow_m3_s ≈ r.axial_face_flow_m3_s rtol=1e-10
        rotated = assemble_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,eccentricity_angle_rad=.7)
        @test rotated.grid.nodes === grid.nodes
        @test rotated.grid.radial_gaps ≈ rotated.normal_gap_m .* (1e3/cosd(angle)) rtol=1e-14
        @test grid == original
    end
    p = parameters()
    theta, z = [0.,1.2,2.8,4.6], [0.,25.,60.,100.]
    grid = housing_grid(p,4,3;theta_nodes=theta,z_nodes=z)
    theta[1], z[1] = .1, 1.
    @test grid.theta_nodes == [0.,1.2,2.8,4.6]
    @test grid.z_nodes == [0.,25.,60.,100.]
end

@testset "Positive conductances, energy identity and shared face balance" begin
    for angle in (0.,10.), rpm in (-2000.,2000.)
        grid = graded_grid(parameters(;angle),12,7)
        a = assemble_full_film(grid;rpm,viscosity_pa_s=.0277,end_pressure_pa=(2e5,1e5))
        nt,nz = grid.n_theta,grid.n_axial
        @test issymmetric(a.A)
        @test isposdef(Symmetric(Matrix(a.A)))
        @test all(>(0),a.gt) && all(>(0),a.gz)
        @test all(<=(0),nonzeros(a.A-spdiagm(0=>diag(a.A))))
        x = reshape([sin(k)+.1cos(3k) for k in 1:nt*nz],nt,nz)
        energy = sum(a.gt[i,j]*(x[i,j]-x[mod1(i+1,nt),j])^2 for i in 1:nt,j in 1:nz)
        energy += sum(a.gz[i,j]*(x[i,j-1]-x[i,j])^2 for i in 1:nt,j in 2:nz)
        energy += sum(a.gz[i,1]*x[i,1]^2+a.gz[i,end]*x[i,end]^2 for i in 1:nt)
        @test dot(vec(x),a.A*vec(x)) ≈ energy rtol=2e-14
        pressure = vec(x).*1e6
        flux = full_film_fluxes(a,pressure)
        @test flux.balance ≈ a.A*pressure-a.b rtol=1e-12
        @test sum(flux.balance) ≈ sum(flux.qz[:,end])-sum(flux.qz[:,1]) atol=1e-18
        @test flux.qt[end,1] ≈ a.adv[end,1]+a.gt[end,1]*(x[end,1]-x[1,1])*1e6 rtol=1e-14
    end
end

@testset "Exact stationary cylinder and analytical cone on adjusted grids" begin
    for angle in (0.,10.)
        p = parameters(;angle,ecc=0.)
        L, r0 = .1, .05005+.05tand(angle)
        slope, ca = tand(angle),cosd(angle)
        integral(z) = iszero(slope) ? z/r0 : -log1p(-slope*z/r0)/slope
        exact_flow = 2pi*(50e-6*ca)^3*ca/(12*.0277)*4e5/integral(L)
        errors, flow_errors = Float64[], Float64[]
        for nz in (8,16,32,64)
            grid = graded_grid(p,12,nz)
            r = solve_full_film(grid;rpm=0.,viscosity_pa_s=.0277,end_pressure_pa=(5e5,1e5))
            z = (grid.z_nodes[1:end-1]+grid.z_nodes[2:end]).*.0005
            exact = [5e5-4e5*integral(value)/integral(L) for i in 1:12,value in z]
            push!(errors,norm(r.pressure_pa-vec(exact))/norm(exact))
            push!(flow_errors,maximum(abs,vec(sum(r.axial_face_flow_m3_s;dims=1)).-exact_flow)/abs(exact_flow))
            @test r.diagnostics.relative_global_imbalance < 1e-10
            constant = solve_full_film(grid;rpm=0.,viscosity_pa_s=.0277,end_pressure_pa=(2e5,2e5))
            @test maximum(abs,constant.pressure_pa.-2e5)/2e5 < 1e-11
            @test constant.diagnostics.max_cell_imbalance_m3_s < 1e-18
        end
        if angle==0
            @test maximum(errors) < 1e-11
            @test maximum(flow_errors) < 1e-10
        else
            @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
            @test all(flow_errors[2:end] .< 0.3 .* flow_errors[1:end-1])
            @test errors[end] < 3e-5
        end
        println("Adjusted stationary angle=$angle pressure=$errors flow=$flow_errors")
    end
end

# Independently differentiated continuous theta-z equation, not A*exact.
function rotating_reference(theta,z,angle,rpm)
    L,c,e,mu,P = .1,50e-6,30e-6,.0277,1e7
    ca,k = cosd(angle),tand(angle)
    st,ct = sincos(theta)
    rj = .05+(L/2-z)*k
    r = rj+c
    D = sqrt(rj^2-e^2*st^2)
    h = ca*(c-e*ct+e^2*st^2/(rj+D))
    ht,hz = ca*(e*st+e^2*st*ct/D),ca*k*(rj/D-1)
    K,Kt,Kz = h^3/(12mu),h^2*ht/(4mu),h^2*hz/(4mu)
    f,df = z*(L-z)/L^2,(L-2z)/L^2
    pressure = P*st*f
    pt,ptt,pz,pzz = P*ct*f,-pressure,P*st*df,-2P*st/L^2
    diffusion = (Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
    wedge = rpm*2pi/60*rj*ht/(2ca)
    return (;pressure,source=wedge-diffusion)
end

@testset "Rotating eccentric cylinder/cone MMS on graded meshes" begin
    for angle in (0.,10.), rpm in (-2000.,2000.)
        errors = Float64[]
        for nt in (32,64,128)
            nz = nt÷2
            grid = graded_grid(parameters(;angle),nt,nz)
            a = assemble_full_film(grid;rpm,viscosity_pa_s=.0277)
            exact, source = zeros(nt,nz),zeros(nt,nz)
            for j in 1:nz,i in 1:nt
                value = rotating_reference(a.axes.theta_centres[i],a.axes.z_centres[j],angle,rpm)
                exact[i,j] = value.pressure
                source[i,j] = value.source*a.axes.theta_widths[i]*a.axes.z_widths[j]
            end
            computed = a.A \ (a.b+vec(source))
            weights = vec([a.axes.theta_widths[i]*a.axes.z_widths[j] for i in 1:nt,j in 1:nz])
            push!(errors,sqrt(sum((computed-vec(exact)).^2 .* weights)/sum(vec(exact).^2 .* weights)))
        end
        @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
        @test errors[end] < .002
        println("Adjusted rotating MMS angle=$angle rpm=$rpm errors=$errors")
    end
end

@testset "Housing force includes the actual cone angle and cell extents" begin
    for angle in (0.,10.,15.)
        p = parameters(;angle,radius=35.,length=70.)
        grid = graded_grid(p,32,16)
        pressure = 2e5
        force = housing_pressure_force(grid,fill(pressure,length(grid.quads)))
        large = (.03505+.035tand(angle))
        small = (.03505-.035tand(angle))
        @test abs(force[1]) < 1e-9
        @test abs(force[2]) < 1e-9
        if angle==0
            @test force[3] == 0
        else
            @test force[3] ≈ pressure*pi*(large^2-small^2) rtol=2e-13
        end
    end
end

@testset "Invalid meshes and operating inputs reject explicitly" begin
    grid = graded_grid(parameters(),16,8)
    for kwargs in ((;rpm=Inf,viscosity_pa_s=.0277),(;rpm=0.,viscosity_pa_s=0.),
            (;rpm=0.,viscosity_pa_s=NaN),(;rpm=0.,viscosity_pa_s=.0277,end_pressure_pa=(NaN,0.)),
            (;rpm=0.,viscosity_pa_s=.0277,eccentricity_angle_rad=Inf))
        @test_throws ArgumentError solve_full_film(grid;kwargs...)
    end
    bad_theta = copy(grid.theta_nodes);bad_theta[2]=bad_theta[1]
    bad_nodes = copy(grid.nodes);bad_nodes[1]=(0.,0.,0.)
    bad_quads = copy(grid.quads);bad_quads[1]=reverse(bad_quads[1])
    bad_z = copy(grid.z_nodes);bad_z[end]=99.
    for bad in ((;grid...,theta_nodes=bad_theta),(;grid...,nodes=bad_nodes),
            (;grid...,quads=bad_quads),(;grid...,z_nodes=bad_z),
            (;grid...,areas=fill(NaN,length(grid.quads))),
            (;grid...,acceptance=(;grid.acceptance...,accepted=false)),
            (;grid...,theta_nodes=grid.theta_nodes[1:end-1]),
            housing_grid(parameters(),3,2),housing_grid(parameters(),4,1))
        @test_throws ArgumentError solve_full_film(bad;rpm=0.,viscosity_pa_s=.0277)
    end
    @test_throws ArgumentError housing_pressure_force(grid,[1.])
    @test_throws ArgumentError housing_pressure_force(grid,fill(NaN,length(grid.quads)))
    @test_throws ArgumentError full_film_fluxes(assemble_full_film(grid;rpm=0.,viscosity_pa_s=.0277),[1.])
end

include(joinpath(@__DIR__,"../reynolds/src/pressure_vtk.jl"))
@testset "Reused meshes retain current-pose gaps and unchanged housing geometry" begin
    for angle in (0.,10.), make_grid in (housing_grid,graded_grid)
        p = parameters(;angle)
        grid = make_grid(p,32,16)
        original = deepcopy(grid)
        previous = grid
        for pose_angle in (.7,.7,0.,0.,-.4,0.)
            before = deepcopy(previous)
            result = solve_full_film(previous;rpm=2000.,viscosity_pa_s=.0277,
                eccentricity_angle_rad=pose_angle)
            fresh = solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,
                eccentricity_angle_rad=pose_angle)
            @test get(result.grid,:eccentricity_angle_rad,0.) == pose_angle
            if get(previous,:eccentricity_angle_rad,0.) == pose_angle
                @test result.grid === previous
            end
            @test result.pressure_pa == fresh.pressure_pa
            @test result.normal_gap_m == fresh.normal_gap_m
            @test result.grid.radial_gaps ≈ result.normal_gap_m .* (1e3/cosd(angle)) rtol=1e-11
            for name in (:nodes,:quads,:areas,:theta_nodes,:z_nodes,:parameters,:quality,:acceptance)
                @test getproperty(result.grid,name) === getproperty(grid,name)
            end
            @test previous == before
            @test grid == original
            mktempdir() do directory
                path = write_conical_pressure_vtk(joinpath(directory,"pressure.vtk"),result)
                lines = readlines(path)
                radial_start = findfirst(==("SCALARS radial_gap_mm double 1"),lines)+2
                normal_start = findfirst(==("SCALARS normal_gap_m double 1"),lines)+2
                count = length(grid.quads)
                radial = parse.(Float64,lines[radial_start:radial_start+count-1])
                normal = parse.(Float64,lines[normal_start:normal_start+count-1])
                @test radial ≈ normal .* (1e3/cosd(angle)) rtol=1e-11
            end
            previous = result.grid
        end
    end
end

include(joinpath(@__DIR__, "../../mesh/2d/surface_feed.jl"))
@testset "Fitted meshes require an explicit rim pressure" begin
    fitted = surface_feed_grid(parameters(),128,64;feed_diameter=5.,feed_z=50.)
    @test fitted.acceptance.accepted
    @test_throws ArgumentError full_film_axes(fitted)
    @test_throws ArgumentError solve_full_film(fitted;rpm=0.,viscosity_pa_s=.0277)
end
