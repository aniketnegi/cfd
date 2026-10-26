using Test
include(joinpath(@__DIR__, "surface_feed.jl"))

# Historical implementation is a comparison fixture only, never a production
# dependency. It includes its own BearingParameters in a separate module.
module ReviewedTrial
include(joinpath(@__DIR__, "../../../expts/conical-fitted-inlet/operator.jl"))
include(joinpath(@__DIR__, "../../../expts/conical-inlet-transition/layout.jl"))
end

p=BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
                    bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=.6)
prior=ReviewedTrial.controlled_layout(ReviewedTrial.fitted_mesh()).mesh
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

@testset "Promoted fitted inlet: exact match to reviewed geometry" begin
    for level in (1,2,4)
        g=surface_feed_grid(p,128level,64level;feed_diameter=5.,feed_z=50.)
        level==1 && (global coarse=g)
        reference=ReviewedTrial.subdivide_layout(prior,level)
        # Exact equality, not merely matching a minimum-quality scalar.
        for field in (:nodes,:quads,:centre,:areas,:owner,:neighbour,:mid,:normals,
                      :lengths,:boundary,:localnodes,:corner_scores,:face_alignment)
            @test getproperty(g.surface,field)==getproperty(reference,field)
        end
        @test !hasproperty(g,:feed_cells)
        @test g.layout=="controlled-fitted-land-v1"
        @test length(g.quads)==8416level^2
        @test length(g.feed_edges)==32level
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
            ((g.surface.R .+ (.05 .- g.surface.centre[:,2].*g.surface.ca).*tand(10.))./g.surface.R))*1e6;rtol=1e-13)
        hole=whole_area-sum(g.areas);push!(hole_errors,abs(hole/mouth_area-1))
        @test 0<hole<1.001mouth_area
        @test hole_errors[end]<.007/level^2
        println("Promoted level=$level $(g.quality); mouth area error=$(hole_errors[end])")
    end
    @test hole_errors[3]<hole_errors[1]/10
end

@testset "Public fitted export, units and non-stepped boundary" begin
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

@testset "Promoted chart with the checked conservative field operator" begin
    reference=ReviewedTrial.BearingParameters(journal_midpoint_radius=50.,radial_clearance=.05,
        bearing_length=100.,semicone_angle_deg=10.,eccentricity_ratio=0.)
    pose=ReviewedTrial.floating_bush_pose(reference,(-9.152501266813789e-7,1.7844427760798323e-7,0.))
    state=ReviewedTrial.solve_fitted(coarse.surface,pose)
    @test maximum(abs,state.force-[-18.24857748898862,-1002.1463387752476,321.13450985407525])<1e-6
    @test state.net_feed*60000≈.030887303265211042 rtol=1e-10
    @test minimum(state.pressure)>=0
    @test 0<=minimum(state.fill)<=maximum(state.fill)<=1
    @test state.imbalance<1e-10
    println("Promoted coarse held-pose force=$(state.force), net feed=$(state.net_feed*60000) L/min; NOT equilibrium")
end

@testset "No silent fallback or unsupported geometry substitution" begin
    for diameter in (0.,-1.,Inf,NaN,4.,6.),z in (50.,)
        @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=diameter,feed_z=z)
    end
    for z in (NaN,-1.,2.5,49.,97.5,102.5)
        @test_throws ArgumentError surface_feed_grid(p,128,64;feed_diameter=5.,feed_z=z)
    end
    for (nt,nz) in ((2,64),(128,0),(256,64),(1024,512),(8,4))
        @test_throws ArgumentError surface_feed_grid(p,nt,nz;feed_diameter=5.,feed_z=50.)
    end
    for (R,L,c,angle) in ((50.,100.,.05,0.),(51.,100.,.05,10.),(50.,90.,.05,10.),(50.,100.,.06,10.))
        other=BearingParameters(journal_midpoint_radius=R,bearing_length=L,radial_clearance=c,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        @test_throws ArgumentError surface_feed_grid(other,128,64;feed_diameter=5.,feed_z=50.)
    end
end
