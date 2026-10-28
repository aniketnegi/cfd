include(joinpath(@__DIR__,"benchmark.jl"))

@testset "Exact factor-two grid families" begin
    for aspect in (.1,.5,.6,.7,.8,.9,1.,1.25,10.)
        counts=refinement_axial_count.(aspect,(128,256,512))
        @test counts[2]==2counts[1]
        @test counts[3]==2counts[2]
    end
    @test refinement_axial_count.(.7,(128,256,512))==(44,88,176)
    @test_throws ArgumentError refinement_axial_count(.7,192)
    @test_throws ArgumentError refinement_axial_count(-1.,128)
end

@testset "Axial interpolation and saved-coordinate contract" begin
    for nz in (2,7,8,9,256)
        nodes=collect(range(0.,100.;length=nz+1))
        centres=(nodes[1:end-1]+nodes[2:end])/2
        pressure=[i+3z for i in 1:4,z in centres]
        @test axial_profile(pressure,centres,25.) ≈ collect(1:4).+75. rtol=1e-14
        @test axial_profile(pressure,centres,50.) ≈ collect(1:4).+150. rtol=1e-14
        @test axial_profile(pressure,centres,first(centres))==pressure[:,1]
        @test axial_profile(pressure,centres,last(centres))==pressure[:,end]
        p=(;radius_mm=50.,bearing_length=100.,radial_clearance=.05)
        scale=6MU*(RPM*2pi/60)*.05^2/(.05e-3)^2
        result=(;pressure_pa=vec(pressure*scale),grid=(;n_theta=4,n_axial=nz,z_nodes=nodes))
        mktempdir() do directory
            save_profiles(p,result,"linear";directory)
            data,_=readdlm(joinpath(directory,"linear-circumferential.tsv"),'\t',Float64;header=true)
            @test data[:,3] ≈ collect(1:4).+75. rtol=1e-14
            @test data[:,2] ≈ collect(1:4).+150. rtol=1e-14
            @test all(data[:,6].==25.)
            @test all(data[:,5].==50.)
        end
    end
    centres=[1.,3.,8.,10.]
    pressure=[2+5z for i in 1:3,z in centres]
    @test axial_profile(pressure,centres,2.5) ≈ fill(14.5,3)
    @test_throws ArgumentError axial_profile(pressure,centres,0.)
    @test_throws ArgumentError axial_profile(pressure,reverse(centres),2.5)
    @test_throws DimensionMismatch axial_profile(pressure,centres[1:3],2.5)
end
