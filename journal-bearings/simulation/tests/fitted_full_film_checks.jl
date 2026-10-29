using Test, LinearAlgebra, SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/surface_feed.jl"))
include(joinpath(@__DIR__,"../reynolds/src/full_film.jl"))

parameters(angle;ecc=.6) = BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,
    radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=ecc)

@testset "Tensor coefficients expose the common oriented face operator" begin
    @test !isdefined(@__MODULE__,:write_vtk)
    @test all(id.name != "Gmsh" for id in keys(Base.loaded_modules))
    for angle in (0.,10.)
        p = parameters(angle)
        grid = housing_grid(p,24,12;
            theta_nodes=[2pi*t+.2sin(2pi*t) for t in (0:23)./24],
            z_nodes=[100*(t+.1sin(2pi*t)/(2pi)) for t in (0:12)./12])
        a = assemble_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(2e5,1e5))
        @test a.A ≈ a.faces.divergence*a.faces.pressure_flux rtol=2e-15
        @test a.b ≈ -a.faces.divergence*a.faces.constant rtol=2e-14
        trial = [1e6*sin(i) for i in eachindex(a.b)]
        old = full_film_fluxes(a,trial)
        common = full_film_face_fluxes(a,trial)
        @test common.face_flow ≈ old.face_flow rtol=2e-14
        @test common.balance ≈ old.balance rtol=2e-14
        @test collect(values(common.boundary_outflow)) ≈ collect(values(old.boundary_outflow)) rtol=2e-14
        @test sum(common.balance) ≈ sum(values(common.boundary_outflow)) atol=1e-18
        @test all(sum(a.faces.divergence[:,f]) == 0 for f in findall(>(0),a.faces.neighbour))
        @test_throws ArgumentError assemble_full_film(grid;rpm=0.,viscosity_pa_s=.0277,feed_pressure_pa=0.)
        @test_throws ArgumentError full_film_face_fluxes(a,[1.])
    end
end

# Cartesian surface tangents provide an independent pressure-area vector integral.
function cartesian_force(grid, pressure)
    s,p = grid.surface,grid.parameters
    rule = eigen(SymTridiagonal(zeros(6),[k/sqrt(4k^2-1) for k in 1:5]))
    points,weights = rule.values,2 .* rule.vectors[1,:].^2
    result = zeros(3)
    ca,k = cosd(p.semicone_angle_deg),tand(p.semicone_angle_deg)
    for (cell,xy) in enumerate(s.localnodes), (u,wu) in zip(points,weights), (v,wv) in zip(points,weights)
        shape = [(1-u)*(1-v),(1+u)*(1-v),(1+u)*(1+v),(1-u)*(1+v)]/4
        du,dv = [-(1-v),1-v,1+v,-(1+v)]/4,[-(1-u),-(1+u),1+u,1-u]/4
        x,xu,xv = xy'*shape,xy'*du,xy'*dv
        theta,z = x[1]/s.R,x[2]*ca
        radius = (p.journal_midpoint_radius+p.radial_clearance)*1e-3+(p.bearing_length*.0005-z)*k
        tx = [-radius/s.R*sin(theta),radius/s.R*cos(theta),0.]
        ty = [-k*ca*cos(theta),-k*ca*sin(theta),ca]
        result .+= pressure[cell]*wu*wv.*cross(tx*xu[1]+ty*xu[2],tx*xv[1]+ty*xv[2])
    end
    return result
end

# Independently differentiated theta-z equation, transformed to the metre chart.
function reference(theta,z,angle,rpm;phi=.37)
    L,c,e,mu,P = .1,50e-6,30e-6,.0277,1e7
    ca,k = cosd(angle),tand(angle)
    st,ct = sincos(theta-phi)
    rj = .05+(L/2-z)*k
    r = rj+c
    root = sqrt(rj^2-e^2*st^2)
    h = ca*(c-e*ct+e^2*st^2/(rj+root))
    ht,hz = ca*(e*st+e^2*st*ct/root),ca*k*(rj/root-1)
    K,Kt,Kz = h^3/(12mu),h^2*ht/(4mu),h^2*hz/(4mu)
    f,df = z*(L-z)/L^2,(L-2z)/L^2
    pressure = P*sin(theta)*f
    pt,ptt,pz,pzz = P*cos(theta)*f,-pressure,P*sin(theta)*df,-2P*sin(theta)/L^2
    diffusion = (Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
    wedge = rpm*2pi/60*rj*ht/(2ca)
    return (;pressure,source=(wedge-diffusion)*ca/.05005)
end

@testset "Fitted signed full film: boundaries, fluxes, scaling and pose reuse" begin
    for angle in (0.,10.)
        p = parameters(angle)
        grid = surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=50.)
        original = deepcopy(grid)
        @test grid.parameters === p
        @test length(grid.quads) != grid.n_theta*grid.n_axial
        @test_throws ArgumentError full_film_axes(grid)
        @test_throws ArgumentError solve_full_film(grid;rpm=0.,viscosity_pa_s=.0277)
        for feed in (Inf,NaN)
            @test_throws ArgumentError solve_full_film(grid;rpm=0.,viscosity_pa_s=.0277,feed_pressure_pa=feed)
        end
        constant = solve_full_film(grid;rpm=0.,viscosity_pa_s=.0277,
            end_pressure_pa=(2e5,2e5),feed_pressure_pa=2e5)
        @test maximum(abs,constant.pressure_pa.-2e5)/2e5 < 1e-10
        @test maximum(abs,constant.face_flow_m3_s) < 1e-16
        @test !hasproperty(constant,:axial_face_flow_m3_s)
        pressure_force = cartesian_force(grid,constant.pressure_pa)
        @test collect(constant.housing_pressure_force_n) ≈ pressure_force rtol=1e-10 atol=1e-9
        angle == 0 && @test constant.housing_pressure_force_n[3] == 0.
        kw = (;viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),feed_pressure_pa=3e5,
            eccentricity_angle_rad=.37)
        still = solve_full_film(grid;kw...,rpm=0.)
        moving = solve_full_film(grid;kw...,rpm=2000.)
        reverse = solve_full_film(grid;kw...,rpm=-2000.)
        @test minimum(moving.pressure_pa) < 0
        @test (moving.pressure_pa+reverse.pressure_pa)/2 ≈ still.pressure_pa rtol=1e-9
        @test still.boundary_outflow_m3_s.feed < 0
        @test still.boundary_outflow_m3_s.large > 0
        @test still.boundary_outflow_m3_s.small > 0
        @test abs(sum(values(moving.boundary_outflow_m3_s))) < 1e-10*sum(abs,values(moving.boundary_outflow_m3_s))
        @test moving.diagnostics.relative_l1_cell_imbalance < 1e-10
        @test moving.diagnostics.relative_global_imbalance < 1e-10
        shifted = solve_full_film(grid;rpm=2000.,viscosity_pa_s=.0277,
            eccentricity_angle_rad=.37,end_pressure_pa=(123456.,123456.),feed_pressure_pa=423456.)
        @test shifted.pressure_pa ≈ moving.pressure_pa .+ 123456 rtol=1e-10
        @test shifted.face_flow_m3_s ≈ moving.face_flow_m3_s rtol=1e-10
        viscous = solve_full_film(grid;kw...,rpm=2000.,viscosity_pa_s=.0554)
        @test viscous.pressure_pa ≈ 2moving.pressure_pa-still.pressure_pa rtol=1e-10
        @test moving.grid.nodes === grid.nodes
        @test moving.grid.surface === grid.surface
        @test moving.grid.radial_gaps ≈ moving.normal_gap_m .* (1e3/cosd(angle)) rtol=1e-13
        reused = solve_full_film(moving.grid;kw...,rpm=2000.)
        @test reused.pressure_pa == moving.pressure_pa
        @test reused.grid === moving.grid
        reset = solve_full_film(moving.grid;kw...,rpm=2000.,eccentricity_angle_rad=0.)
        fresh = solve_full_film(grid;kw...,rpm=2000.,eccentricity_angle_rad=0.)
        @test reset.pressure_pa == fresh.pressure_pa
        @test reset.normal_gap_m == fresh.normal_gap_m
        @test grid == original
        # Affine axial pressure has an exact face flux even on skew fitted cells.
        # On a cone its independently integrated volume source is constant.
        s = grid.surface
        concentric = parameters(angle;ecc=0.)
        pressure_gradient = 3e6
        boundary_pressure(x,y,tag) = 1e5+pressure_gradient*y
        affine_faces = fitted_face_operator(s,concentric;rpm=0.,viscosity_pa_s=.0277,
            eccentricity_angle_rad=0.,boundary_pressure)
        affine = 1e5 .+ pressure_gradient.*s.centre[:,2]
        K = (50e-6*cosd(angle))^3/(12*.0277)
        expected_flux = [-K*(.05005+(.05-x[2]*s.ca)*tand(angle))/s.R*
            normal[2]*width*pressure_gradient for (x,normal,width) in zip(s.mid,s.normals,s.lengths)]
        @test affine_faces.pressure_flux*affine+affine_faces.constant ≈ expected_flux rtol=2e-12 atol=1e-20
        source = K*pressure_gradient*tand(angle)*s.ca/s.R .* s.areas
        affine_A = affine_faces.divergence*affine_faces.pressure_flux
        affine_b = -(affine_faces.divergence*affine_faces.constant)
        @test affine_A \ (affine_b+source) ≈ affine rtol=1e-10
        a = assemble_full_film(grid;kw...,rpm=2000.)
        @test a.A == a.faces.divergence*a.faces.pressure_flux
        @test a.b == -(a.faces.divergence*a.faces.constant)
        flux = full_film_fluxes(a,moving.pressure_pa)
        @test norm(flux.balance-(a.A*moving.pressure_pa-a.b),Inf) < 1e-18
        @test all(a.faces.neighbour[a.faces.boundary .== :feed] .== 0)
        @test_throws ArgumentError full_film_fluxes(a,zeros(length(a.b)-1))
        @test_throws ArgumentError housing_pressure_force(grid,[1.])
        @test_throws ArgumentError solve_full_film((;grid...,acceptance=(;accepted=false));kw...,rpm=0.)
        @test_throws ArgumentError solve_full_film((;grid...,parameters=parameters(5.));kw...,rpm=0.)
    end
end

@testset "Fitted rotating eccentric cylinder/cone manufactured fields" begin
    for angle in (0.,10.)
        errors = Dict(rpm=>Float64[] for rpm in (-2000.,2000.))
        for nt in (128,256)
            p = parameters(angle)
            grid = surface_feed_grid(p,nt,nt÷2;feed_diameter=5.,feed_z=50.)
            s = grid.surface
            for rpm in (-2000.,2000.)
                boundary_pressure(x,y,tag) = reference(x/s.R,y*s.ca,angle,rpm).pressure
                faces = fitted_face_operator(s,p;rpm,viscosity_pa_s=.0277,eccentricity_angle_rad=.37,boundary_pressure)
                exact = [reference(s.centre[i,1]/s.R,s.centre[i,2]*s.ca,angle,rpm).pressure for i in axes(s.centre,1)]
                source = [reference(s.centre[i,1]/s.R,s.centre[i,2]*s.ca,angle,rpm).source*s.areas[i] for i in axes(s.centre,1)]
                A = faces.divergence*faces.pressure_flux
                b = -(faces.divergence*faces.constant)
                computed = A \ (b+source)
                error = sqrt(sum((computed-exact).^2 .* grid.areas)/sum(exact.^2 .* grid.areas))
                push!(errors[rpm],error)
                @test norm(A*computed-b-source,Inf)/maximum(abs,source) < 1e-9
                println("Fitted MMS angle=$angle rpm=$rpm nt=$nt relative_L2=$error")
            end
        end
        for history in values(errors)
            @test history[2] < .65history[1]
            @test history[2] < .003
        end
    end
end

include(joinpath(@__DIR__,"../reynolds/src/pressure_vtk.jl"))
@testset "Fitted force and export use supplied dimensions and cone angle" begin
    p = BearingParameters(journal_midpoint_radius=40.,bearing_length=80.,radial_clearance=.08,
        semicone_angle_deg=7.,eccentricity_ratio=.4)
    grid = surface_feed_grid(p,160,52;feed_diameter=4.,feed_z=35.)
    s = grid.surface
    pressure = [2e5*(1+.1cos(s.centre[i,1]/s.R)+.2s.centre[i,2]*s.ca/.08) for i in axes(s.centre,1)]
    @test collect(housing_pressure_force(grid,pressure)) ≈ cartesian_force(grid,pressure) rtol=1e-11
    result = solve_full_film(grid;rpm=0.,viscosity_pa_s=.04,
        feed_pressure_pa=2e5,end_pressure_pa=(2e5,2e5))
    @test result.pressure_pa ≈ fill(2e5,length(grid.quads)) rtol=1e-10
    mktempdir() do directory
        path = write_conical_pressure_vtk(joinpath(directory,"pressure.vtk"),result)
        lines = readlines(path)
        rim = length(grid.feed_edges)
        cells = rim+length(grid.quads)
        @test "CELL_DATA $cells" in lines
        scalars = findall(startswith("SCALARS "),lines)
        for (index,start) in enumerate(scalars)
            stop = index == length(scalars) ? length(lines) : scalars[index+1]-1
            @test length(lines[start+2:stop]) == cells
        end
        for (name,expected) in (("pressure_gauge_pa",result.pressure_pa),("normal_gap_m",result.normal_gap_m))
            start = findfirst(==("SCALARS $name double 1"),lines)+2
            actual = parse.(Float64,lines[start:start+cells-1])
            @test all(isnan,actual[1:rim])
            @test actual[rim+1:end] == expected
        end
        @test any(startswith("LINES "),lines)
    end
end

include(joinpath(@__DIR__,"../reynolds/src/swift_stieber.jl"))
@testset "Fitted cavitation requires an explicit rim pressure" begin
    # Neither a missing supply value nor an incomplete mesh is inferred.
    @test_throws ArgumentError solve_swift_stieber((;feed_edges=Tuple{Int,Int}[]);rpm=0.,viscosity_pa_s=.0277)
end
