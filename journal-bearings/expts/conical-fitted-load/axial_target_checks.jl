include(joinpath(@__DIR__,"checks.jl"))
@testset "Changing axial target leaves fixed-seat radial equations unchanged" begin
    seed=(-.9174126793137724e-6,.16737680363128737e-6)
    a=radial_equilibrium(problem,50.,seed)
    b=radial_equilibrium(problem,50.,seed;axial_load_n=500.)
    @test a.converged && b.converged
    @test a.state.pose.displacement_m==b.state.pose.displacement_m
    @test a.state.film.pressure==b.state.film.pressure
    @test a.state.film.fill==b.state.film.fill
    @test b.residual_n[3]≈a.residual_n[3]-300 atol=1e-10
    @test b.required_axial_support_n≈500-b.state.bush_force_n[3]
    @test_throws ArgumentError radial_equilibrium(problem,50.,seed;axial_load_n=0.)
    @test_throws ArgumentError radial_equilibrium(problem,50.,seed;axial_load_n=NaN)
end
