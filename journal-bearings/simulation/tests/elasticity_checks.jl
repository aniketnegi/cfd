using Test,LinearAlgebra,SparseArrays
include(joinpath(@__DIR__,"../../bearing.jl"))
include(joinpath(@__DIR__,"../../mesh/2d/grid.jl"))
include(joinpath(@__DIR__,"../reynolds/src/full_film.jl"))
include(joinpath(@__DIR__,"../elasticity/src/housing.jl"))
include(joinpath(@__DIR__,"../verification/fsi/elastic_reference.jl"))

function elastic_cube(n)
    id(i,j,k)=i+1+(n+1)*j+(n+1)^2*k
    nodes=zeros((n+1)^3,3)
    for k in 0:n,j in 0:n,i in 0:n;nodes[id(i,j,k),:]=[i/n,j/n,k/n];end
    cells=NTuple{8,Int}[]
    for k in 0:n-1,j in 0:n-1,i in 0:n-1
        push!(cells,(id(i,j,k),id(i+1,j,k),id(i+1,j+1,k),id(i,j+1,k),
            id(i,j,k+1),id(i+1,j,k+1),id(i+1,j+1,k+1),id(i,j+1,k+1)))
    end
    boundary=[i for i in axes(nodes,1) if any(x->x==0 || x==1,nodes[i,:])]
    return nodes,cells,boundary
end

function check_elasticity()
@testset "Hex8 elasticity consistency and affine patch" begin
    for x in (-.7,.3),y in (-.4,.8),z in (-.8,.6)
        N,dN=hex8_shape(x,y,z)
        @test sum(N)≈1. atol=1e-15
        @test norm(vec(sum(dN;dims=1)),Inf)<1e-15
    end
    X,cells,boundary=elastic_cube(2)
    a=assemble_elasticity(X,cells;youngs_modulus_pa=97e9,poisson_ratio=.3)
    @test a.volume_m3≈1. rtol=1e-14
    @test norm(a.K-a.K',Inf)/norm(a.K,Inf)<1e-14
    translate=repeat([.1,-.2,.3],size(X,1));rotate=vec(hcat([-x[2] for x in eachrow(X)],[x[1] for x in eachrow(X)],zeros(size(X,1)))')
    @test norm(a.K*translate,Inf)/norm(a.K,Inf)<1e-15
    @test norm(a.K*rotate,Inf)/norm(a.K,Inf)<1e-15
    affine=[1e-4 2e-5 -3e-5;4e-5 -2e-4 1e-5;3e-5 1e-5 5e-5]
    exact=vec((X*affine')')
    constrained=constrain_elasticity(a,elastic_dofs(boundary))
    r=solve_elasticity(constrained,zeros(size(a.K,1));fixed_values=exact[constrained.fixed])
    @test r.displacement_m≈exact rtol=1e-12
    @test r.relative_residual<1e-12
    strain=[affine[1,1],affine[2,2],affine[3,3],affine[1,2]+affine[2,1],affine[2,3]+affine[3,2],affine[1,3]+affine[3,1]]
    @test r.strain_energy_j≈dot(strain,a.C*strain)/2 rtol=1e-12
    @test elastic_stress(a,r.displacement_m).maximum_engineering_strain≈norm(strain,Inf) rtol=1e-12
    @test_throws ArgumentError isotropic_stiffness(0.,.3)
    @test_throws ArgumentError isotropic_stiffness(97e9,.5)
    @test_throws ArgumentError constrain_elasticity(a,Int[])
    inverted=[(1,4,3,2,5,8,7,6)]
    @test_throws ArgumentError assemble_elasticity(X[collect(cells[1]),:],inverted;youngs_modulus_pa=1.,poisson_ratio=.3)
end
@testset "Continuous 3D elastic manufactured solution" begin
    errors=Float64[]
    for n in (2,4,8,16)
        X,cells,boundary=elastic_cube(n)
        a=assemble_elasticity(X,cells;youngs_modulus_pa=97e9,poisson_ratio=.3)
        a=constrain_elasticity(a,elastic_dofs(boundary))
        r=solve_elasticity(a,elastic_manufactured_load(a))
        push!(errors,elastic_manufactured_error(a,r.displacement_m))
    end
    println("3D elastic MMS errors: $errors")
    @test all(errors[2:end] .< .4 .* errors[1:end-1])
    @test errors[end]<.03
end
@testset "Curved housing pressure force, work and restrained cylinder reference" begin
    for angle in (0.,10.)
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,
            semicone_angle_deg=angle,eccentricity_ratio=.6)
        grid=housing_grid(p,32,8)
        a=assemble_housing(grid;youngs_modulus_pa=97e9,poisson_ratio=.3,wall_thickness_m=.005,n_radial=2)
        pressure=[1e6*(1+.2sin(i)) for i in eachindex(grid.quads)]
        f=a.pressure_operator*pressure
        force=vec(sum(reshape(f,3,:);dims=2))
        @test force≈collect(housing_pressure_force(grid,pressure)) rtol=1e-12 atol=1e-8
        @test a.interface_areas_m2≈grid.areas.*1e-6 rtol=1e-12
        r=solve_housing(a,pressure)
        @test all(iszero,r.displacement_m[a.fixed])
        @test r.strain_energy_j>0
        @test r.relative_residual<1e-10
        @test vec(sum(reshape(r.load_n+r.reaction_n,3,:);dims=2))≈zeros(3) atol=1e-7
        work=sum(pressure.*a.interface_areas_m2.*r.normal_displacement_cell_mean_m)
        @test work≈dot(r.displacement_m,f) rtol=1e-13
        twice=solve_housing(a,2pressure)
        @test twice.displacement_m≈2r.displacement_m rtol=1e-12
        @test a.minimum_jacobian>0
    end
    errors=Float64[]
    for (nt,nr) in ((16,2),(32,4),(64,8))
        p=BearingParameters(journal_midpoint_radius=50.,bearing_length=100.,radial_clearance=.05,eccentricity_ratio=0.)
        g=housing_grid(p,nt,2);mesh=housing_solid_mesh(g;wall_thickness_m=.005,n_radial=nr)
        a=assemble_elasticity(mesh.nodes,mesh.cells;youngs_modulus_pa=97e9,poisson_ratio=.3)
        # This analytic control fixes all uz (plane strain), not the free-ended design.
        fixed=unique([elastic_dofs(mesh.outer_nodes);collect(3:3:size(a.K,1))])
        a=constrain_elasticity(a,fixed)
        transfer=housing_pressure_operator(mesh);r=solve_elasticity(a,transfer.matrix*fill(1e6,length(g.quads)))
        E,nu=97e9,.3;mu=E/(2*(1+nu));lambda=E*nu/((1+nu)*(1-2nu))
        inner=.05005;outer=inner+.005
        A=-1e6/(2*(lambda+mu)+2mu*outer^2/inner^2);B=-A*outer^2
        exact=Float64[];computed=Float64[]
        for node in 1:mesh.base_nodes
            x,y,_=mesh.nodes[node,:];radius=hypot(x,y)
            push!(exact,A*radius+B/radius)
            u=r.displacement_m[3node-2:3node]
            push!(computed,(x*u[1]+y*u[2])/radius)
        end
        push!(errors,norm(computed-exact)/norm(exact))
    end
    println("Restrained-cylinder displacement errors: $errors")
    @test all(errors[2:end] .< .4 .* errors[1:end-1])
    @test errors[end]<.01
end
end
if abspath(PROGRAM_FILE)==@__FILE__;check_elasticity();end
