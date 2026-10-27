using Test
include(joinpath(@__DIR__, "../../bearing.jl"))
include(joinpath(@__DIR__, "../2d/conical/grid.jl"))
include(joinpath(@__DIR__, "../2d/conical/surface_feed.jl"))

# Independent circle/ray calculation using the same represented dimensions,
# but high-precision trigonometry and arithmetic rather than Cartesian hypot.
function reference_radii(p, theta, z)
    return setprecision(BigFloat,256) do
        radius = BigFloat(p.journal_midpoint_radius) +
            (BigFloat(p.half_bearing_length)-BigFloat(z))*BigFloat(p.slope)
        housing = radius + BigFloat(p.radial_clearance)
        sine,cosine = sincos(BigFloat(theta))
        eccentricity = BigFloat(p.eccentricity)
        ray = eccentricity*cosine + sqrt(radius^2-(eccentricity*sine)^2)
        return radius,housing,ray,housing-ray
    end
end

@testset "One radial geometry for point mapping and measurement" begin
    for angle in (0.,10.,15.), length_mm in (80.,123.4), eccentricity in (0.,.6,.99)
        p = BearingParameters(journal_midpoint_radius=50.,bearing_length=length_mm,
            radial_clearance=.05,semicone_angle_deg=angle,eccentricity_ratio=eccentricity)
        dimensions = bearing_geometry(p)
        @test !hasproperty(dimensions,:eccentricity)
        @test dimensions.slope === p.slope
        @test dimensions.half_bearing_length === p.half_bearing_length
        for z in (0.,length_mm/3,length_mm/2,length_mm), theta in (0.,.37,pi/2,pi,5.2,2pi)
            theta = Float64(theta)
            radial = bearing_radial_geometry(p,theta,z)
            radius,housing,ray,gap = reference_radii(p,theta,z)
            tolerance = 8eps(Float64(housing))
            @test abs(BigFloat(journal_radius(p,z))-radius) <= tolerance
            @test abs(BigFloat(radial.housing_radius)-housing) <= tolerance
            @test abs(BigFloat(radial.journal_ray_radius)-ray) <= tolerance
            @test abs(BigFloat(radial.radial_gap)-gap) <= tolerance
            @test housing_radius(p,z) === housing_radius(dimensions,z)
            @test journal_radius(p,z) === journal_radius(dimensions,z)
            @test radial.radial_gap === radial.housing_radius-radial.journal_ray_radius
            @test radial.radial_gap > 0
            for fraction in (0.,.25,.7,1.)
                point = bearing_point(p,theta,z,fraction)
                expected_radius = fraction==1 ? radial.housing_radius :
                    radial.journal_ray_radius+fraction*radial.radial_gap
                @test point == (expected_radius*radial.cos_theta,expected_radius*radial.sin_theta,z)
            end
            x,y,_ = bearing_point(p,theta,z,0.)
            @test hypot(x-p.eccentricity,y) ≈ journal_radius(p,z) atol=tolerance rtol=0
        end
        for theta in (NaN,Inf), z in (0.,length_mm)
            @test_throws ArgumentError bearing_radial_geometry(p,theta,z)
            @test_throws ArgumentError bearing_point(p,theta,z,1.)
        end
        for z in (NaN,Inf,-1.,nextfloat(length_mm))
            @test_throws ArgumentError bearing_radial_geometry(p,0.,z)
            @test_throws ArgumentError bearing_point(p,0.,z,1.)
        end
    end
end

@testset "Housing points and dimension records are pose independent" begin
    dimensions = (journal_midpoint_radius=40.,bearing_length=123.4,
        radial_clearance=.1,semicone_angle_deg=15.)
    centred = BearingParameters(;dimensions...,eccentricity_ratio=0.)
    offset = BearingParameters(;dimensions...,eccentricity_ratio=.99)
    @test bearing_geometry(centred) == bearing_geometry(offset)
    for z in (0.,61.7,123.4), theta in range(0.,2pi;length=17)
        @test bearing_point(centred,theta,z,1.) == bearing_point(offset,theta,z,1.)
    end
end

@testset "Uniform and adjusted meshes use the shared radii exactly" begin
    p = BearingParameters(journal_midpoint_radius=40.,bearing_length=123.4,
        radial_clearance=.1,semicone_angle_deg=15.,eccentricity_ratio=.6)
    theta_nodes = [0.,1.2,2.8,4.6]
    z_nodes = [0.,30.,70.,123.4]
    for spacing in ((;),(;theta_nodes,z_nodes))
        grid = housing_grid(p,4,3;spacing...)
        @test all(node -> 0 <= node[3] <= p.bearing_length,grid.nodes)
        for j in 1:3, i in 1:4
            if isempty(spacing)
                theta = (i-.5)*(2pi/4)
                z = (j-.5)*(p.bearing_length/3)
                dtheta = 2pi/4
                ds = (p.bearing_length/3)/cosd(p.semicone_angle_deg)
            else
                high_theta = i==4 ? theta_nodes[1]+2pi : theta_nodes[i+1]
                theta = (theta_nodes[i]+high_theta)/2
                z = (z_nodes[j]+z_nodes[j+1])/2
                dtheta = high_theta-theta_nodes[i]
                ds = (z_nodes[j+1]-z_nodes[j])/cosd(p.semicone_angle_deg)
            end
            geometry = bearing_radial_geometry(p,theta,z)
            cell = i+4*(j-1)
            @test grid.radial_gaps[cell] === geometry.radial_gap
            @test grid.areas[cell] === geometry.housing_radius*dtheta*ds
        end
    end
end

@testset "Fitted point, gap and metric consumers use the same dimensions" begin
    p = BearingParameters(journal_midpoint_radius=45.,bearing_length=123.4,
        radial_clearance=.06,semicone_angle_deg=12.,eccentricity_ratio=.6)
    grid = surface_feed_grid(p,127,65;feed_diameter=4.,feed_z=61.7)
    surface = grid.surface
    dimensions = bearing_geometry(p)
    @test surface.geometry_mm == dimensions
    @test all(eachindex(grid.nodes)) do node
        chart_x,chart_y = surface.nodes[node]
        theta = chart_x/surface.R
        z = ControlledSurfaceInlet.chart_axial_position_mm(dimensions,surface.ca,chart_y)
        grid.nodes[node] == bearing_point(p,theta,z,1.)
    end
    @test all(eachindex(grid.quads)) do cell
        theta = surface.centre[cell,1]/surface.R
        z = ControlledSurfaceInlet.chart_axial_position_mm(dimensions,surface.ca,surface.centre[cell,2])
        radial = bearing_radial_geometry(p,theta,z)
        grid.radial_gaps[cell] === radial.radial_gap &&
            grid.areas[cell] === surface.areas[cell]*radial.housing_radius/(surface.R*1e3)*1e6
    end
    @test ControlledSurfaceInlet.chart_axial_position_mm(dimensions,surface.ca,0.) === 0.
    @test ControlledSurfaceInlet.chart_axial_position_mm(dimensions,surface.ca,
        surface.L/surface.ca) === p.bearing_length
    @test ControlledSurfaceInlet.housing_radius === housing_radius
    @test ControlledSurfaceInlet.bearing_geometry === bearing_geometry
end
