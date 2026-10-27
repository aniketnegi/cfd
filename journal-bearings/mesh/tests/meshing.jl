using Test
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))
include(joinpath(@__DIR__, "../vtk.jl"))

@testset "Eccentricity checked when bearing parameters are created" begin
    for angle in (0., 30.)
        small_end_radius = 1. - tan(deg2rad(angle))
        parameters(clearance) = BearingParameters(journal_midpoint_radius=1.,
            bearing_length=2., radial_clearance=clearance,
            semicone_angle_deg=angle, eccentricity_ratio=.5)

        # Exercise the strict boundary, not just clearly invalid dimensions.
        valid = parameters(2 * prevfloat(small_end_radius))
        @test valid.eccentricity < small_end_radius
        @test_throws ArgumentError parameters(2 * small_end_radius)
        @test_throws ArgumentError parameters(2 * nextfloat(small_end_radius))
        @test all(isfinite, bearing_point(valid, pi/2, 2., 0.))
        @test_throws ArgumentError bearing_point(valid, 0., -1., 1.)
        @test_throws ArgumentError bearing_point(valid, 0., 3., 1.)
        @test_throws ArgumentError bearing_point(valid, NaN, 0., 1.)
        @test_throws ArgumentError bearing_point(valid, 0., 0., 1.1)
    end
end

p=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
                    bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=.6)
hole_errors=Float64[]
coarse=nothing
# Independent integration of the true drilling footprint on the cone.
r,z0,n=2.5,50.,10_000
dz=2r/n
mouth_area=sum(1:n) do j
    z=z0-r+(j-.5)*dz;rh=p.housing_midpoint_radius+(p.half_bearing_length-z)*p.slope
    rh*2asin(sqrt(r^2-(z-z0)^2)/rh)*dz/cosd(p.semicone_angle_deg)
end
whole_area=2pi*p.housing_midpoint_radius*p.bearing_length/cosd(p.semicone_angle_deg)

@testset "Fitted inlet geometry and connectivity" begin
    for (n_theta, n_axial) in ((128,64), (256,128), (512,256))
        g=surface_feed_grid(p,n_theta,n_axial;feed_diameter=5.,feed_z=50.)
        n_theta==128 && (global coarse=g)
        @test !hasproperty(g,:feed_cells)
        @test g.layout=="fitted-land-direct-v1"
        @test (g.surface.n_theta, g.surface.n_axial) == (n_theta, n_axial)
        patch_cells = g.surface.patch_cells
        layers = g.surface.transition_layers
        @test length(g.quads) == g.n_theta*g.n_axial - patch_cells^2 + 4patch_cells*layers
        @test length(g.feed_edges) == 4patch_cells
        @test length(g.nodes) == length(g.quads) + g.n_theta + 2patch_cells - 1
        @test length(g.nodes)==length(unique(g.nodes))
        @test g.quality.minimum_scaled_corner>.9
        @test g.quality.minimum_face_alignment>.9
        @test g.quality.maximum_neighbour_area_ratio<2
        @test all(g.nodes) do (x,y,z)
            rh=p.housing_midpoint_radius+(p.half_bearing_length-z)*p.slope
            abs(hypot(x,y)-rh)<1e-11 && -1e-12<=z<=p.bearing_length+1e-12
        end
        @test all(g.feed_edges) do (a,b)
            all((a,b)) do v
                x,y,z=g.nodes[v]
                y>0 && abs(hypot(x,z-z0)-r)<1e-11
            end
        end
        edgecounts=Dict{Tuple{Int,Int},Int}();valence=zeros(Int,length(g.nodes))
        for q in g.quads
            for i in 1:4
                valence[q[i]]+=1
                e=Tuple(sort([q[i],q[mod1(i+1,4)]]));edgecounts[e]=get(edgecounts,e,0)+1
            end
        end
        @test all(v->v in (1,2),values(edgecounts))
        @test all(v->v in (2,4,5),valence)
        @test count(==(5),valence)==4
        @test length(g.nodes)-length(edgecounts)+length(g.quads)==-1
        @test count(==(1),values(edgecounts))==2g.n_theta+length(g.feed_edges)
        @test all(e->edgecounts[e]==1,g.feed_edges)
        @test minimum(g.radial_gaps)>0
        @test all(p.minimum_gap-1e-11 .<= g.radial_gaps .<= p.maximum_gap+1e-11)
        @test isapprox(sum(g.areas),sum(g.surface.areas.*
            ((g.surface.R .+ (g.surface.L/2 .- g.surface.centre[:,2].*g.surface.ca).*p.slope)./g.surface.R))*1e6;rtol=1e-13)
        hole=whole_area-sum(g.areas);push!(hole_errors,abs(hole/mouth_area-1))
        @test 0<hole<1.001mouth_area
        @test hole_errors[end] < .007*(8/patch_cells)^2
        println("Inlet $n_theta×$n_axial $(g.quality); mouth area error=$(hole_errors[end])")
    end
    @test hole_errors[3]<hole_errors[1]/10
end

@testset "VTK export, units and inlet boundary" begin
    g=coarse
    mktempdir() do dir
        file=joinpath(dir,"fitted.vtk");write_vtk(file,g);lines=readlines(file)
        ip=findfirst(s->startswith(s,"POINTS "),lines)
        @test parse(Int,split(lines[ip])[2])==length(g.nodes)
        @test [Tuple(parse.(Float64,split(s))) for s in lines[ip+1:ip+length(g.nodes)]]==g.nodes
        il=findfirst(s->startswith(s,"LINES "),lines)
        @test parse(Int,split(lines[il])[2])==length(g.feed_edges)
        @test [Tuple(parse.(Int,split(s)[2:end]).+1) for s in lines[il+1:il+length(g.feed_edges)]]==g.feed_edges
        iq=findfirst(s->startswith(s,"POLYGONS "),lines)
        @test [Tuple(parse.(Int,split(s)[2:end]).+1) for s in lines[iq+1:iq+length(g.quads)]]==g.quads
        nc=length(g.quads)+length(g.feed_edges)
        @test "CELL_DATA $nc" in lines
        for name in ("radial_gap_mm","analytic_area_mm2","feed_rim","surface_corner_score")
            i=findfirst(s->startswith(s,"SCALARS $name "),lines)
            a=parse.(Float64,lines[i+2:i+1+nc])
            if name=="feed_rim"
                @test all(==(1),a[1:length(g.feed_edges)])
                @test all(==(0),a[length(g.feed_edges)+1:end])
            elseif name=="analytic_area_mm2"
                @test a[length(g.feed_edges)+1:end]==g.areas
            elseif name=="radial_gap_mm"
                @test a[length(g.feed_edges)+1:end]==g.radial_gaps
            else
                @test a[length(g.feed_edges)+1:end]==g.surface.corner_scores
            end
        end
        @test !any(s->startswith(s,"SCALARS feed_cell "),lines)
        @test_throws ArgumentError write_vtk(file,g)
        @test_throws ArgumentError write_vtk(joinpath(dir,"mask.vtk"),g;feed_cells=falses(length(g.quads)))
        plain=housing_grid(p,16,4);file2=joinpath(dir,"plain.vtk");write_vtk(file2,plain)
        text=read(file2,String)
        @test occursin("POLYGONS 64 320",text)
        @test occursin("CELL_DATA 64",text)
        @test !occursin("LINES ",text)
        @test sum(plain.areas)≈whole_area rtol=1e-12
    end
    centred=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    g0=surface_feed_grid(centred,128,64;feed_diameter=5.,feed_z=50.)
    @test g0.nodes==g.nodes && g0.quads==g.quads
    @test all(isapprox.(g0.radial_gaps,.05;atol=1e-12))
    @test g0.radial_gaps!=g.radial_gaps
end

@testset "Invalid inputs and unavailable inlet transitions" begin
    for diameter in (0.,-1.,Inf,NaN),z in (50.,)
        @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=diameter,feed_z=z)
    end
    for z in (NaN,Inf,-1.,0.,2.5,97.5,100.,102.5)
        @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=z)
    end
    for (nt,nz) in ((2,64),(128,0),(-1,64),(8,4))
        @test_throws ArgumentError surface_feed_grid(p,nt,nz;feed_diameter=5.,feed_z=50.)
    end
    # These openings fit axially, but the surrounding transition does not.
    for z in (5.,95.)
        @test_throws r"physical support crosses domain boundary" surface_feed_grid(p,128,64;
            feed_diameter=5.,feed_z=z)
    end
    @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=100.,feed_z=50.)
end

# The inlet deliberately shares the background builder, without loading IO or solvers.
module InletDependencies
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))
end

@testset "Mesh dependency boundaries" begin
    @test isdefined(InletDependencies, :surface_feed_grid)
    @test isdefined(InletDependencies, :housing_grid)
    @test InletDependencies.ControlledSurfaceInlet.BearingParameters === InletDependencies.BearingParameters
    @test InletDependencies.ControlledSurfaceInlet.background_grid === InletDependencies.background_grid
    @test InletDependencies.ControlledSurfaceInlet.patch_axis === InletDependencies.patch_axis
    @test !isdefined(InletDependencies, :write_vtk)
    @test !isdefined(@__MODULE__, :solve_conical_jfo)
    @test !isdefined(@__MODULE__, :equilibrate_floating_bush)
    @test all(id.name != "Gmsh" for id in keys(Base.loaded_modules))
end

@testset "Regular cylinder and cone mesh" begin
    for angle in (0., 10.)
        parameters = BearingParameters(journal_midpoint_radius=50., radial_clearance=.05,
            bearing_length=100., semicone_angle_deg=angle, eccentricity_ratio=.6)
        grid = housing_grid(parameters, 16, 4)
        @test length(grid.nodes) == 16 * 5
        @test length(grid.quads) == 16 * 4
        @test length(grid.nodes) == length(unique(grid.nodes))
        @test grid.quads[16] == (16, 1, 17, 32)
        @test minimum(grid.radial_gaps) > 0
        @test all(parameters.minimum_gap-1e-11 .<= grid.radial_gaps .<= parameters.maximum_gap+1e-11)
        @test sum(grid.areas) ≈ 2pi * parameters.housing_midpoint_radius *
            parameters.bearing_length / cosd(angle) rtol=1e-12
        @test_throws ArgumentError housing_grid(parameters, 2, 4)
        @test_throws ArgumentError housing_grid(parameters, 16, 0)
        mktempdir() do dir
            file = joinpath(dir, "regular.vtk")
            write_vtk(file, grid)
            original = read(file, String)
            # Regular exports intentionally accept an existing temporary file:
            # the named fluid runner uses that to replace results atomically.
            write_vtk(file, grid)
            @test read(file, String) == original
            @test_throws ArgumentError write_vtk(joinpath(dir, "bad.vtk"), grid;
                feed_cells=[true])
        end
    end
end

@testset "Uniform axial endpoints remain on the bearing" begin
    for length_mm in (0.1, 100.), n_axial in (11, 97)
        parameters = BearingParameters(journal_midpoint_radius=50., radial_clearance=.05,
            bearing_length=length_mm, semicone_angle_deg=10., eccentricity_ratio=.6)
        grid = housing_grid(parameters, 8, n_axial)
        @test first(grid.nodes)[3] == 0
        @test last(grid.nodes)[3] == length_mm
        @test all(node -> 0 <= node[3] <= length_mm, grid.nodes)
        @test sum(grid.areas) ≈ 2pi * parameters.housing_midpoint_radius *
            length_mm / cosd(parameters.semicone_angle_deg) rtol=1e-12
    end
end

@testset "Shared periodic background construction" begin
    theta_nodes = [(i - 1) * 2pi/8 for i in 1:8]
    z_nodes = [0., 50., 100.]
    background = background_grid(theta_nodes, z_nodes)
    @test (background.n_theta, background.n_axial) == (8, 2)
    @test length(background.nodes) == 24
    @test length(background.quads) == 16
    @test background.nodes[1] == (0., 0.)
    @test background.nodes[9] == (0., 50.)
    @test background.nodes[24] == (theta_nodes[end], 100.)
    @test background.quads[1] == (1, 2, 10, 9)
    @test background.quads[8] == (8, 1, 9, 16)
    @test all(quad -> length(unique(quad)) == 4, background.quads)
    @test Set(Iterators.flatten(background.quads)) == Set(eachindex(background.nodes))

    # Constructor output must not retain aliases to the supplied position vectors.
    original_first = first(background.nodes)
    theta_nodes[1] = .1
    z_nodes[1] = 1.
    @test first(background.nodes) == original_first
    @test eltype(background_grid([0, 1, 3, 5], [0, 25, 100]).nodes) == NTuple{2,Float64}

    for angles in ([0., 1.], [0., 1., 1.], [0., 2., 1.], [-.1, 1., 2.],
                   [0., 1., 2pi], [0., NaN, 2.], [0., 1., Inf])
        @test_throws ArgumentError background_grid(angles, [0., 100.])
    end
    for axial_positions in ([0.], [0., 0.], [1., 0.], [0., Inf], [0., NaN])
        @test_throws ArgumentError background_grid([0., 1., 2.], axial_positions)
    end
end

@testset "Housing geometry with supplied spacing" begin
    custom_theta = [0., .4, 2.1, 5.6]
    custom_z = [0., 7., 42., 100.]
    for angle in (0., 10.)
        parameters = BearingParameters(journal_midpoint_radius=50., radial_clearance=.05,
            bearing_length=100., semicone_angle_deg=angle, eccentricity_ratio=.6)
        for spacing in ((; theta_nodes=custom_theta, z_nodes=custom_z),
                        (; theta_nodes=custom_theta), (; z_nodes=custom_z))
            grid = housing_grid(parameters, 4, 3; spacing...)
            @test length(grid.nodes) == 16
            @test length(grid.quads) == 12
            @test minimum(grid.areas) > 0
            @test sum(grid.areas) ≈ 2pi * parameters.housing_midpoint_radius *
                parameters.bearing_length / cosd(angle) rtol=1e-12

            angles = haskey(spacing, :theta_nodes) ? custom_theta : [(i-1)*2pi/4 for i in 1:4]
            axial_positions = haskey(spacing, :z_nodes) ? custom_z : [(j-1)*100/3 for j in 1:4]
            for j in 1:3, i in 1:4
                theta_low = angles[i]
                theta_high = i == 4 ? angles[1] + 2pi : angles[i+1]
                theta_mid = (theta_low + theta_high) / 2
                z_mid = (axial_positions[j] + axial_positions[j+1]) / 2
                journal_radius = parameters.journal_midpoint_radius +
                    (parameters.half_bearing_length-z_mid) * parameters.slope
                housing_radius = journal_radius + parameters.radial_clearance
                journal_ray = parameters.eccentricity * cos(theta_mid) +
                    sqrt(journal_radius^2 - (parameters.eccentricity*sin(theta_mid))^2)
                cell = i + 4*(j-1)
                @test grid.radial_gaps[cell] ≈ housing_radius-journal_ray atol=1e-11
                @test grid.areas[cell] ≈ housing_radius * (theta_high-theta_low) *
                    (axial_positions[j+1]-axial_positions[j]) / cosd(angle) rtol=1e-12
            end
        end
        @test_throws DimensionMismatch housing_grid(parameters, 4, 3; theta_nodes=[0., 1., 2.])
        @test_throws DimensionMismatch housing_grid(parameters, 4, 3; z_nodes=[0., 100.])
        @test_throws ArgumentError housing_grid(parameters, 4, 3; z_nodes=[1., 7., 42., 100.])
        @test_throws ArgumentError housing_grid(parameters, 4, 3; z_nodes=[0., 7., 42., 99.])

        # The cut may start at a nonzero angle; its last face still closes a full turn.
        shifted = housing_grid(parameters, 4, 3; theta_nodes=custom_theta .+ .1, z_nodes=custom_z)
        @test sum(shifted.areas) ≈ sum(housing_grid(parameters, 4, 3).areas) rtol=1e-12
    end
end

@testset "Adjusted spacing is shared with the unported mesh" begin
    parameters = BearingParameters(journal_midpoint_radius=50., radial_clearance=.05,
        bearing_length=100., semicone_angle_deg=10., eccentricity_ratio=.6)
    theta_nodes = patch_axis(128, 32, 2pi, pi/2, 7.5/parameters.housing_midpoint_radius, 8)[1:end-1]
    z_nodes = patch_axis(64, 32, parameters.bearing_length, 50., 7.5, 8)
    original_theta = copy(theta_nodes)
    original_z = copy(z_nodes)
    adjusted = housing_grid(parameters, 128, 64; theta_nodes, z_nodes)
    uniform = housing_grid(parameters, 128, 64)
    @test adjusted.quads == uniform.quads
    @test adjusted.nodes != uniform.nodes
    @test adjusted.areas != uniform.areas
    @test sum(adjusted.areas) ≈ sum(uniform.areas) rtol=1e-12
    @test all(diff(z_nodes[25:41]) .≈ 1.875)
    @test theta_nodes == original_theta
    @test z_nodes == original_z
    @test_throws ArgumentError patch_axis(16, 8, 100., 50., 7.5, 8)
    @test_throws ArgumentError patch_axis(64, 32, 100., 5., 7.5, 8)
end
