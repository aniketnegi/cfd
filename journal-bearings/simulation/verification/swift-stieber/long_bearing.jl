include(joinpath(@__DIR__,"cases.jl"))
include(joinpath(SOURCE,"simulation/reynolds/src/jfo_slider.jl"))
length(ARGS)==1 || error("usage: long_bearing.jl FRESH_OUTPUT")
const OUTPUT=fresh_output(ARGS[1];suite="long-bearing-ODE")

function quadrature(order)
    pair=eigen(SymTridiagonal(zeros(order),[k/sqrt(4k^2-1) for k in 1:order-1]))
    nodes,weights=pair.values,2 .* pair.vectors[1,:].^2
    return (f,a,b)->(b-a)/2*sum(weights.*f.((a+b)/2 .+ (b-a)/2 .* nodes))
end
function continuum(ecc;order=192)
    integrate=quadrature(order)
    h(t)=1+ecc*cos(t)
    equation(r)=integrate(t->(h(t)-h(r))/h(t)^3,0.,r)
    lo,hi=pi,2pi
    @assert equation(lo)>0>equation(hi)
    for _ in 1:60
        r=(lo+hi)/2
        equation(r)>0 ? (lo=r) : (hi=r)
    end
    rupture=(lo+hi)/2
    gradient(t)=6*(h(t)-h(rupture))/h(t)^3
    pressure(t)=t<rupture ? integrate(gradient,0.,t) : 0.
    fx=-integrate(t->gradient(t)*sin(t),0.,rupture)
    fy=integrate(t->gradient(t)*cos(t),0.,rupture)
    return (;rupture,pressure,load=hypot(fx,fy),attitude=rad2deg(atan(abs(fy),abs(fx))))
end
function main()
    ispath(joinpath(OUTPUT,"long-bearing.tsv")) && error("refusing overwrite")
    # Printed finite-bearing study p111, infinite-width row.
    targets=[(.2130,2.73,60.3),(.3204,4.02,52.3),(.4289,5.44,44.4),
             (.5383,7.09,37.1),(.6479,9.26,29.9),(.8660,21.4,15.7)]
    r,c,mu,omega=.05,.05e-3,.02752,2000*2pi/60
    length=2pi*r;scale=mu*omega*r^2/c^2
    open(joinpath(OUTPUT,"long-bearing.tsv"),"w") do io
        println(io,"epsilon\tn\tpublished_P0\texact_P0\tcomputed_P0\tpublished_beta_deg\texact_beta_deg\tcomputed_beta_deg\texact_rupture_deg\tlast_positive_face_deg\tpressure_relative_L2\tquadrature_load_relative_change\titerations")
        for (ecc,published,beta) in targets
            exact=continuum(ecc);check=continuum(ecc;order=96)
            qerror=abs(check.load/exact.load-1)
            @assert qerror<1e-10
            gap(x)=c*(1+ecc*cos((x+length/2)/r))
            errors=Float64[]
            for n in (256,512,1024,2048)
                a=assemble_jfo_slider(gap,n;length_m=length,viscosity_pa_s=mu,
                    entrainment_m_s=omega*r/2,end_pressure_pa=(0.,0.),cavitation_pressure_pa=0.)
                state=pressure_obstacle(a.A,a.b,a.pressure_scale_pa,maximum(abs,a.adv))
                pressure=state.pressure/scale
                theta=(a.x_m .+ length/2)./r
                exact_values=exact.pressure.(theta)
                error=norm(pressure-exact_values)/norm(exact_values)
                push!(errors,error)
                edges=(a.face_x_m .+ length/2)./r
                fx=sum(pressure.*(sin.(edges[2:end])-sin.(edges[1:end-1])))
                fy=sum(pressure.*(cos.(edges[1:end-1])-cos.(edges[2:end])))
                last=findlast(>(1e-8),state.pressure)
                data=(ecc,n,published,exact.load,hypot(fx,fy),beta,exact.attitude,
                    rad2deg(atan(abs(fy),abs(fx))),rad2deg(exact.rupture),last*360/n,error,qerror,state.iterations)
                println(io,join(data,'\t'));flush(io)
            end
            @assert errors[end]<1e-5
            @printf("epsilon=%.4f exactP0=%.8f printed=%.3f beta=%.6f pressureL2=%.3e\n",ecc,exact.load,published,exact.attitude,errors[end])
        end
    end
end
main()
