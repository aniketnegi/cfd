include(joinpath(@__DIR__,"benchmark.jl"))
include(joinpath(@__DIR__,"spectral.jl"))

function cubic_angular_value(values,theta)
    nt=length(values);x=mod(theta,2pi)*nt/(2pi)+.5
    start=floor(Int,x)-1
    return sum(values[mod1(k,nt)]*prod((x-j)/(k-j) for j in start:start+3 if j!=k)
        for k in start:start+3)
end
function reference_field(basis,theta,w,beta)
    angular=sin.(theta .* transpose(1:size(basis.vectors,1)))*basis.vectors
    axial=[basis.amplitudes[n]*(1-cosh_ratio(basis.alpha[n],z,beta))
        for n in eachindex(basis.alpha),z in w]
    return angular*axial
end
function main(output;full=false)
    fresh_output(output;suite="full-film-fields")
    cases=vec([(aspect,ecc,nt) for aspect in (.1,1.,10.),ecc in (.2,.8,.95,.99),nt in (128,256,512)])
    append!(cases,[(.1,.99,1024),(.1,.99,2048),(1.,.8,1024),(10.,.99,1024)])
    full || filter!(case->case[1]==1. && case[2]==.8 && case[3]<=512,cases)
    open(joinpath(output,"fields.tsv"),"w") do io
        first=true
        for (aspect,ecc,nt) in cases
            nz=refinement_axial_count(aspect,nt)
            p=DonaldsonCylinder(50.,100aspect,.05,ecc,0.,0.)
            elapsed=@elapsed result=solve_full_film(benchmark_grid(p,nt,nz);
                rpm=RPM,viscosity_pa_s=MU,end_pressure_pa=(0.,0.))
            obs=observables(p,result)
            basis=spectral_basis(ecc,320)
            theta=[(i-.5)*2pi/nt for i in 1:nt]
            w=[((j-.5)/nz*2-1)*aspect for j in 1:nz]
            reference=reference_field(basis,theta,w,aspect)
            r,c=p.radius_mm*1e-3,p.radial_clearance*1e-3
            scale=6MU*(RPM*2pi/60)*r^2/c^2
            pressure=reshape(result.pressure_pa,nt,nz)/scale
            error=pressure-reference
            midline=(pressure[:,nz÷2]+pressure[:,nz÷2+1])/2
            mobility=[normal_film_gap(p,t,0.)^3*r*2pi/nt/(12MU) for t in theta]
            gradient=-result.axial_face_flow_m3_s[:,1]./mobility .* (r/scale)
            data=(;L_D=aspect,epsilon=ecc,nt,nz,obs...,
                spectral_load=spectral_load(basis,aspect),
                spectral_pressure_160=spectral_pressure(basis,deg2rad(160.),0.,aspect),
                spectral_gradient_160=spectral_end_gradient(basis,deg2rad(160.),aspect),
                pressure_160_cubic=cubic_angular_value(midline,deg2rad(160.)),
                gradient_160_cubic=cubic_angular_value(gradient,deg2rad(160.)),
                whole_field_relative_L2=norm(error)/norm(reference),
                whole_field_relative_Linf=maximum(abs,error)/maximum(abs,reference),seconds=elapsed)
            first && println(io,join(keys(data),'\t'));first=false
            println(io,join(values(data),'\t'));flush(io)
            @printf("LD=%.1f eps=%.2f %dx%d load=%.9f wholefieldL2=%.5g P160=%.8f gradient160=%.8f\n",
                aspect,ecc,nt,nz,obs.load_ratio,data.whole_field_relative_L2,
                data.pressure_160_cubic,data.gradient_160_cubic)
            flush(stdout);GC.gc()
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (1,2) && (length(ARGS)==1 || ARGS[2]=="--full") || error("usage: field_checks.jl FRESH_OUTPUT [--full]")
    main(ARGS[1];full="--full" in ARGS)
end
