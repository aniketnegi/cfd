include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(@__DIR__,"sampling.jl"))
@testset "Coordinate sampling on shifted and graded periodic tensors" begin
    for graded in (false,true)
        errors=Float64[]
        for nt in (32,64,128,256)
            nz=nt÷2
            cg=grid_for(cone(),nt,nz;graded);fg=grid_for(cone(),2nt,2nz;graded)
            c,f=full_film_axes(cg),full_film_axes(fg)
            values=[sin(t)+.7cos(2t)+2z+.2*(z/.1)^2 for t in f.theta_centres,z in f.z_centres]
            sampled=[sample_pressure(f.theta_centres,f.z_centres,vec(values),t,z)
                for t in c.theta_centres,z in c.z_centres]
            exact=[sin(t)+.7cos(2t)+2z+.2*(z/.1)^2 for t in c.theta_centres,z in c.z_centres]
            error=norm(sampled-exact)/norm(exact);push!(errors,error)
            linear=[3+2z for t in f.theta_centres,z in f.z_centres]
            restored=[sample_pressure(f.theta_centres,f.z_centres,vec(linear),t,z)
                for t in c.theta_centres,z in c.z_centres]
            @test maximum(abs,restored-[3+2z for t in c.theta_centres,z in c.z_centres])<3e-15
            @test sample_pressure(f.theta_centres,f.z_centres,vec(values),0.,.05)≈
                sample_pressure(f.theta_centres,f.z_centres,vec(values),2pi,.05) atol=1e-14
        end
        @test errors[end]<1e-4
        @test all(errors[2:end] .< 0.3 .* errors[1:end-1])
    end
    nt,nz=64,32
    c=full_film_axes(grid_for(cone(),nt,nz));f=full_film_axes(grid_for(cone(),2nt,2nz))
    fine=[sin(t)*(1+z) for t in f.theta_centres,z in f.z_centres]
    exact=[sin(t)*(1+z) for t in c.theta_centres,z in c.z_centres]
    sampled=[sample_pressure(f.theta_centres,f.z_centres,vec(fine),t,z)
        for t in c.theta_centres,z in c.z_centres]
    wrong=[sum(fine[2i-1:2i,2j-1:2j])/4 for i in 1:nt,j in 1:nz]
    @test norm(sampled-exact)/norm(exact)<1e-13
    @test norm(wrong-exact)/norm(exact)>.01
end
