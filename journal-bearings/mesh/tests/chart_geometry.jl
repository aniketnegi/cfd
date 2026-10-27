using Test
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))

const Inlet = ControlledSurfaceInlet

# Use the untranslated shoelace formula at high precision as an independent
# reference for the exact Float64 vertices, not for unrepresentable ideal inputs.
function reference_polygon(coordinates)
    return setprecision(BigFloat, 256) do
        points = BigFloat.(coordinates)
        corners = size(points, 1)
        crosses = [points[i,1]*points[mod1(i+1,corners),2] -
            points[mod1(i+1,corners),1]*points[i,2] for i in 1:corners]
        area = sum(crosses)/2
        centroid = [sum((points[i,axis]+points[mod1(i+1,corners),axis])*crosses[i]
            for i in 1:corners)/(6area) for axis in 1:2]
        return area, centroid
    end
end

function check_polygon(coordinates)
    original = copy(coordinates)
    reference_area, reference_centroid = reference_polygon(coordinates)
    area = Inlet.cell_chart_area(coordinates)
    centroid = Inlet.cell_chart_centroid(coordinates, area)
    @test abs(BigFloat(area)-reference_area) <= 8eps(abs(Float64(reference_area)))
    for axis in 1:2
        span = maximum(coordinates[:,axis])-minimum(coordinates[:,axis])
        scale = max(span,abs(Float64(reference_centroid[axis])))
        @test abs(BigFloat(centroid[axis])-reference_centroid[axis]) <= 8eps(scale)
    end
    @test Inlet.cell_chart_centroid(coordinates) == centroid
    @test coordinates == original
end

@testset "Stable signed polygon areas and centroids" begin
    shapes = (
        [0. 0.; 4. 0.; 4. 2.; 0. 2.],
        [0. 0.; 4. 0.; 3. 2.; 0. 1.],
        [0. 0.; 2. 1.; 3. 3.; 1. 2.],
        [0. 0.; 3. 0.; 0. 2.],
    )
    for shape in shapes, scale in (1.,2.0^-10,2.0^-20),
        offset in ([0.,0.],[.31,.09],[2.0^20,-2.0^20])
        coordinates = scale .* shape .+ transpose(offset)
        check_polygon(coordinates)
        check_polygon(reverse(coordinates; dims=1))
    end

    # Dyadic inputs remain exactly representable under these translations.
    # Otherwise translating vertices can change their represented geometry.
    shape = 2.0^-20 .* shapes[2]
    area = Inlet.cell_chart_area(shape)
    for offset in ([1.,-1.],[2.0^20,-2.0^20])
        translated = shape .+ transpose(offset)
        @test Inlet.cell_chart_area(translated) === area
        @test Inlet.cell_chart_area(reverse(translated; dims=1)) === -area
        for shift in 0:3
            check_polygon(circshift(translated,(shift,0)))
        end
    end
    collapsed = [0. 0.; 1. 0.; 2. 0.; 3. 0.]
    @test iszero(Inlet.cell_chart_area(collapsed))
    @test_throws ArgumentError Inlet.cell_chart_centroid(collapsed)
end

@testset "Unwrapped seam geometry and orientation" begin
    nodes = [(x,y) for y in (0.,.01) for x in (.02,.14,.26,.38)]
    quads = [(i,mod1(i+1,4),mod1(i+1,4)+4,i+4) for i in 1:4]
    parameters = BearingParameters(journal_midpoint_radius=.4/(2pi)*1e3-.05,
        bearing_length=10.,radial_clearance=.05,semicone_angle_deg=0.,eccentricity_ratio=0.)
    geometry = (R=parameters.housing_midpoint_radius*1e-3,L=.01,ca=1.,slope=0.,
        geometry_mm=bearing_geometry(parameters),feed_z_m=.005,inlet_radius_m=.001)
    mesh = Inlet.assemble_surface_mesh(nodes,quads,Set{Tuple{Int,Int}}(); geometry...)
    positions = permutedims(hcat(collect.(nodes)...))
    @test count(==(:large),mesh.boundary) == 4
    @test count(==(:small),mesh.boundary) == 4
    @test sum(mesh.areas) ≈ mesh.period*mesh.L rtol=1e-14
    for cell in eachindex(quads)
        coordinates = mesh.localnodes[cell]
        reference_area, reference_centroid = reference_polygon(coordinates)
        relative = Inlet.cell_coordinates(mesh,positions,quads[cell])
        @test mesh.areas[cell] === Inlet.cell_chart_area(relative)
        @test abs(BigFloat(mesh.areas[cell])-reference_area) <= 8eps(Float64(reference_area))
        expected_x = Float64(mod(reference_centroid[1],BigFloat(mesh.period)))
        distance = mod(mesh.centre[cell,1]-expected_x+mesh.period/2,mesh.period)-mesh.period/2
        @test abs(distance) <= 8eps(mesh.period)
        @test mesh.centre[cell,2] ≈ Float64(reference_centroid[2]) atol=8eps(mesh.L) rtol=0
        check_polygon(coordinates)
    end
    @test maximum(mesh.localnodes[end][:,1]) > mesh.period
    reversed = copy(quads)
    reversed[1] = reverse(reversed[1])
    @test_throws r"nonpositive cell 1" Inlet.assemble_surface_mesh(nodes,reversed,
        Set{Tuple{Int,Int}}(); geometry...)
end

@testset "Assembly and optimization share exactly the same areas" begin
    cases = (
        (50.,100.,.05,10.,5.,50.,128,64),
        (50.,100.,.05,10.,5.,50.,512,256),
        (45.,90.,.06,12.,4.,35.3,127,65),
    )
    for (radius,length_mm,clearance,angle,diameter,feed_z,n_theta,n_axial) in cases
        @testset "$n_theta×$n_axial, angle=$angle" begin
            p = BearingParameters(journal_midpoint_radius=radius,bearing_length=length_mm,
                radial_clearance=clearance,semicone_angle_deg=angle,eccentricity_ratio=.6)
            mesh = Inlet.reference_seed(p,n_theta,n_axial;feed_diameter=diameter,feed_z)
            selection = Inlet.macro_data(mesh)
            @test mesh.areas == selection.initial_areas
            variables = vec(copy(selection.nodes[selection.free,:]))
            energy, gradient = Inlet.layout_energy(mesh,selection,variables)
            unanchored_energy, unanchored_gradient = Inlet.layout_energy(mesh,selection,variables;
                anchor_weight=0.)
            # The area-change anchor must be exactly zero at its reference mesh.
            @test energy === unanchored_energy
            @test gradient == unanchored_gradient
            samples = unique(vcat(round.(Int,range(1,length(mesh.quads);length=64)),
                n_theta==512 ? [130333] : Int[]))
            for cell in samples
                coordinates = mesh.localnodes[cell]
                reference_area, reference_centroid = reference_polygon(coordinates)
                @test abs(BigFloat(mesh.areas[cell])-reference_area) <= 8eps(Float64(reference_area))
                expected_x = Float64(mod(reference_centroid[1],BigFloat(mesh.period)))
                distance = mod(mesh.centre[cell,1]-expected_x+mesh.period/2,mesh.period)-mesh.period/2
                @test abs(distance) <= 8eps(mesh.period)
                @test abs(BigFloat(mesh.centre[cell,2])-reference_centroid[2]) <= 8eps(mesh.L/mesh.ca)
            end
        end
    end
end
