using Test,LinearAlgebra,SparseArrays
include(joinpath(@__DIR__,"../verification/jfo/conical/transport_cases.jl"))
include(joinpath(@__DIR__,"../verification/jfo/conical/transport_sampling.jl"))

function check_tensor_jfo_contract()
@testset "Oriented JFO transport preserves the tensor contract" begin
    for angle in (0.,10.),rpm in (-2000.,2000.)
        g=transport_mesh(angle,24;graded=true)
        a=assemble_conical_jfo(g;rpm,viscosity_pa_s=.0277,end_pressure_pa=(2e5,1e5),cavitation_pressure_pa=0.)
        C=jfo_content_operator(a.faces,length(a.b))
        @test a.B≈-a.faces.divergence*C rtol=2e-15
        P=[1e6*(1+sin(i)) for i in eachindex(a.b)];F=[.6+.3cos(i) for i in eachindex(a.b)]
        old=conical_jfo_fluxes(a,P,F);new=jfo_face_fluxes(a,P,F)
        @test new.balance≈old.balance rtol=2e-13
        @test new.face_flow[1:length(a.b)]≈vec(old.qt) rtol=2e-14
        @test collect(values(new.boundary_outflow))[1:2]≈[-sum(old.qz[:,1]),sum(old.qz[:,end])] rtol=1e-14
    end
end
end

function check_fitted_jfo()
@testset "Fitted JFO inlet, flux algebra, scaling and export" begin
    for angle in (0.,10.)
        g=transport_mesh(angle,128;fitted=true)
        @test length(g.quads)!=g.n_theta*g.n_axial
        original=deepcopy(g)
        kw=(;viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),cavitation_pressure_pa=0.,
            feed_pressure_pa=5e5,eccentricity_angle_rad=3pi/2)
        for rpm in (-2000.,2000.)
            a=assemble_conical_jfo(g;kw...,rpm)
            signed=assemble_full_film(g;rpm,viscosity_pa_s=.0277,end_pressure_pa=(0.,0.),
                feed_pressure_pa=5e5,eccentricity_angle_rad=3pi/2)
            @test a.A==signed.A && a.b==signed.b
            @test a.B==-a.faces.divergence*a.content_flux
            P=[1e6*(1+.2sin(i)) for i in eachindex(a.b)];F=[.6+.2cos(i) for i in eachindex(a.b)]
            q=jfo_face_fluxes(a,P,F)
            @test q.balance≈a.A*P+a.B*(1 .-F)-a.b rtol=2e-12
            @test sum(q.balance)≈sum(values(q.boundary_outflow)) atol=1e-18
            @test jfo_face_fluxes(a,P,ones(length(P))).face_flow≈full_film_face_fluxes(a,P).face_flow rtol=2e-14
            incoming=findall(f->a.faces.neighbour[f]==0 && a.faces.advective[f]<0,eachindex(a.faces.owner))
            outgoing=findall(f->a.faces.neighbour[f]==0 && a.faces.advective[f]>0,eachindex(a.faces.owner))
            @test !isempty(incoming) && !isempty(outgoing)
            @test all(nnz(a.content_flux[f,:])==0 for f in incoming)
            @test all(a.content_flux[f,a.faces.owner[f]]==a.faces.advective[f] for f in outgoing)
            @test all(a.faces.boundary[f]==:feed for f in incoming)
            r=solve_conical_jfo(g;kw...,rpm)
            @test minimum(r.pressure_above_cavitation_pa)>=0
            @test 0<minimum(r.fill_fraction)<1
            @test maximum(r.fill_fraction)<=1
            @test r.diagnostics.complementarity_residual==0
            @test r.diagnostics.relative_max_cell_imbalance<1e-10
            @test r.diagnostics.relative_global_imbalance<1e-10
            @test sum(r.diagnostics.outward_end_flow_m3_s)≈r.diagnostics.net_feed_flow_m3_s rtol=1e-10
            @test !hasproperty(r,:axial_face_flow_m3_s) && !hasproperty(r,:feed_cells)
            @test r.grid.surface===g.surface && r.grid.nodes===g.nodes
        end
        r=solve_conical_jfo(g;kw...,rpm=2000.)
        shifted=solve_conical_jfo(g;kw...,rpm=2000.,cavitation_pressure_pa=101325.,
            end_pressure_pa=(101325.,101325.),feed_pressure_pa=601325.)
        @test shifted.pressure_above_cavitation_pa≈r.pressure_pa rtol=1e-11
        @test shifted.fill_fraction≈r.fill_fraction rtol=1e-11
        @test shifted.face_flow_m3_s≈r.face_flow_m3_s rtol=1e-11
        faster=solve_conical_jfo(g;kw...,rpm=4000.,feed_pressure_pa=1e6)
        viscous=solve_conical_jfo(g;kw...,rpm=2000.,viscosity_pa_s=.0554,feed_pressure_pa=1e6)
        @test faster.pressure_pa≈2r.pressure_pa rtol=1e-10
        @test faster.fill_fraction≈r.fill_fraction rtol=1e-10
        @test faster.face_flow_m3_s≈2r.face_flow_m3_s rtol=1e-10
        @test viscous.pressure_pa≈2r.pressure_pa rtol=1e-10
        @test viscous.fill_fraction≈r.fill_fraction rtol=1e-10
        @test viscous.face_flow_m3_s≈r.face_flow_m3_s rtol=1e-10
        reused=solve_conical_jfo(r.grid;kw...,rpm=2000.)
        @test reused.pressure_pa==r.pressure_pa && reused.fill_fraction==r.fill_fraction
        wet=solve_conical_jfo(g;kw...,rpm=2000.,end_pressure_pa=(1e8,1e8),feed_pressure_pa=1e8)
        full=solve_full_film(g;rpm=2000.,viscosity_pa_s=.0277,end_pressure_pa=(1e8,1e8),
            feed_pressure_pa=1e8,eccentricity_angle_rad=3pi/2)
        @test wet.fill_fraction==ones(length(g.quads))
        @test wet.pressure_pa≈full.pressure_pa rtol=1e-11
        @test wet.face_flow_m3_s≈full.face_flow_m3_s rtol=1e-10 atol=1e-16
        still=solve_conical_jfo(g;kw...,rpm=0.)
        @test still.fill_fraction==ones(length(g.quads))
        @test still.diagnostics.gross_feed_return_m3_s==0
        returning=solve_conical_jfo(g;kw...,rpm=2000.,eccentricity_angle_rad=pi/2)
        @test returning.diagnostics.gross_feed_return_m3_s>0
        @test returning.diagnostics.gross_feed_inflow_m3_s>returning.diagnostics.gross_feed_return_m3_s
        @test g==original
        # Exhaust all16 masks of a genuine mixed fitted subproblem.
        a=assemble_conical_jfo(g;kw...,rpm=2000.)
        transitions=[f for f in eachindex(a.faces.owner) if a.faces.neighbour[f]>0 &&
            (r.pressure_pa[a.faces.owner[f]]>0)!=(r.pressure_pa[a.faces.neighbour[f]]>0)]
        ids=unique([a.faces.owner[transitions[1:4]];a.faces.neighbour[transitions[1:4]]])[1:4]
        rest=setdiff(eachindex(a.b),ids);v=1 .-r.fill_fraction
        rhs=a.b[ids]-a.A[ids,rest]*r.pressure_pa[rest]-a.B[ids,rest]*v[rest]
        A,B=Matrix(a.A[ids,ids]),Matrix(a.B[ids,ids]);feasible=0
        for bits in 0:15
            full=[(bits>>(i-1))&1==1 for i in 1:4]
            M=A*Diagonal(a.pressure_scale_pa.*full)+B*Diagonal(Float64.(.!full))
            iszero(det(M)) && continue
            y=M\rhs;p=a.pressure_scale_pa .* y .* full;void=y .* .!full
            if minimum(p)>=-1e-9 && minimum(void)>=-1e-12
                feasible+=1
                @test p≈r.pressure_pa[ids] rtol=1e-9 atol=1e-6
                @test void≈v[ids] atol=1e-9
            end
        end
        @test feasible>=1
        angle==0 && @test r.housing_pressure_force_n[3]==0.
        mktempdir() do dir
            path=write_conical_jfo_vtk(joinpath(dir,"fitted.vtk"),r)
            lines=readlines(path);rim=length(g.feed_edges);cells=rim+length(g.quads)
            @test "CELL_DATA $cells" in lines
            for name in ("pressure_pa","pressure_above_cavitation_pa","fill_fraction","normal_gap_m")
                start=findfirst(==("SCALARS $name double 1"),lines)+2
                data=parse.(Float64,lines[start:start+cells-1])
                @test all(isnan,data[1:rim])
                @test data[rim+1:end]==getproperty(r,Symbol(name))
            end
            @test_throws ErrorException write_conical_jfo_vtk(path,r)
        end
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,feed_pressure_pa=nothing)
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,feed_pressure_pa=0.)
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,feed_pressure_pa=Inf)
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,feed_cells=falses(length(g.quads)))
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,end_pressure_pa=(-1.,0.))
        @test_throws ArgumentError solve_conical_jfo(g;kw...,rpm=2000.,max_iterations=0)
        @test_throws ErrorException solve_conical_jfo(g;kw...,rpm=2000.,max_iterations=1)
        @test_throws ArgumentError jfo_face_fluxes(assemble_conical_jfo(g;kw...,rpm=2000.),[0.],[1.])
    end
end
@testset "Continuous conical manufactured field and independent FD algebra" begin
    for angle in (0.,10.,20.),rpm in (-2000.,2000.),theta in (.4,1.8,3.9,5.5),z in (.023,.069)
        ref(t,y)=ConicalTransportReference.manufactured(t,y,angle,rpm)
        dt,dz=1e-4,1e-6
        divergence=(ref(theta+dt,z).qx-ref(theta-dt,z).qx)/(2dt*.05005)+
            cosd(angle)*(ref(theta,z+dz).qy-ref(theta,z-dz).qy)/(2dz)
        @test divergence≈ref(theta,z).source atol=1e-8 rtol=2e-6
    end
    errors=Float64[]
    for nt in (32,64,128)
        a=ConicalTransportReference.assemble(20.,.6,nt,nt÷2;ends=(0.,0.))
        rhs=zeros(nt,a.nz-1);exact=zeros(size(rhs))
        for j in 1:a.nz-1,i in 1:nt
            t,z=a.theta[i],a.z[j+1];g=ConicalTransportReference.geometry(t,z,20.,.6;pose=3pi/2)
            f,df,ddf=z*(.1-z)/.1^2,(.1-2z)/.1^2,-2/.1^2
            p=cos(t)*f;pt=-sin(t)*f;ptt=-p;pz=cos(t)*df;pzz=cos(t)*ddf
            beta=(g.rh*g.ca)^2;drift=beta*(3g.hz/g.h-g.k/g.rh)
            rhs[i,j]=-ptt-3g.ht/g.h*pt-beta*pzz-drift*pz;exact[i,j]=p
        end
        computed=a.A\vec(rhs);push!(errors,norm(computed-vec(exact))/norm(exact))
    end
    @test 1.9<log2(errors[1]/errors[2])<2.1
    @test 1.9<log2(errors[2]/errors[3])<2.1
    for angle in (0.,10.,20.),rpm in (-2000.,2000.)
        coarse=transport_mms(transport_mesh(angle,32),rpm)
        fine=transport_mms(transport_mesh(angle,64),rpm)
        @test fine.pressure_relative_L2<.5coarse.pressure_relative_L2
        @test fine.fill_mean_absolute<.7coarse.fill_mean_absolute
        @test fine.complementarity==0 && fine.local_balance<1e-10
    end
    a=ConicalTransportReference.assemble(10.,.6,32,16)
    fd=ConicalTransportReference.solve(a)
    p=vec(fd.pressure[:,2:end-1])./a.scale;v=vec(1 .-fd.fill[:,2:end-1])
    @test norm(a.A*p+a.B*v-a.b,Inf)/max(norm(a.b,Inf),1.)<1e-8
    policy=solve_jfo_active_set(a.A,a.B,a.b,1.)
    @test policy.pressure≈p rtol=1e-7
    @test policy.void≈v atol=1e-7
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    check_tensor_jfo_contract()
    check_fitted_jfo()
end
