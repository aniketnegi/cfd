# Independently differentiated continuous pressure/content fields; no A*exact forcing.
function cylinder_exact(theta,zeta,aspect,ecc;direction=1.)
    H=1-ecc*cos(theta);Ht=ecc*sin(theta)
    t=cos(theta)+.25cos(pi*zeta/aspect)
    tt,ttt=-sin(theta),-cos(theta)
    tz=-.25pi/aspect*sin(pi*zeta/aspect)
    tzz=-.25(pi/aspect)^2*cos(pi*zeta/aspect)
    f,df,ddf=1-(zeta/aspect)^2,-2zeta/aspect^2,-2/aspect^2
    P=Pt=Ptt=Pzz=0.;F=1.;Ft=0.;amplitude=.5
    if t>0
        P=amplitude*f*t^4
        Pt=amplitude*f*4t^3*tt
        Ptt=amplitude*f*(12t^2*tt^2+4t^3*ttt)
        Pzz=amplitude*(ddf*t^4+8df*t^3*tz+f*(12t^2*tz^2+4t^3*tzz))
    else
        F=1-.2*(-t)^4
        Ft=.8*(-t)^3*tt
    end
    source=direction*(Ht*F+H*Ft)/2-(3H^2*Ht*Pt+H^3*(Ptt+Pzz))/12
    return (;P,F,source)
end

function manufactured_case(nt,nz;rpm=RPM,graded=false)
    aspect,ecc=1.,.6
    grid=cylinder_grid(cylinder_parameters(aspect,ecc),nt,nz;graded)
    a=assemble_conical_jfo(grid;rpm,viscosity_pa_s=MU,end_pressure_pa=(0.,0.),cavitation_pressure_pa=0.)
    scale=MU*abs(rpm*2pi/60)*R^2/CLEARANCE^2
    source=zeros(nt*nz);exactP=similar(source);exactF=similar(source)
    points=(-sqrt(3/5),0.,sqrt(3/5));weights=(5/9,8/9,5/9)
    for j in 1:nz,i in 1:nt
        k=i+nt*(j-1)
        theta,z=a.axes.theta_centres[i],a.axes.z_centres[j]
        exact=cylinder_exact(theta,(z-.05)/R,aspect,ecc;direction=sign(rpm))
        exactP[k]=scale*exact.P;exactF[k]=exact.F
        dt,dz=a.axes.theta_widths[i],a.axes.z_widths[j]
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
            t=a.axes.theta_edges[i]+(u+1)*dt/2
            axial=a.axes.z_edges[j]+(v+1)*dz/2
            S=cylinder_exact(t,(axial-.05)/R,aspect,ecc;direction=sign(rpm)).source
            source[k]+=wu*wv*dt*dz/4*abs(rpm*2pi/60)*R*CLEARANCE*S
        end
    end
    state=solve_jfo_active_set(a.A,a.B,a.b+source,a.pressure_scale_pa)
    pressure,fill=state.pressure,1 .-state.void
    flux=conical_jfo_fluxes(a,pressure,fill)
    flow_scale=max(maximum(abs,a.adv),maximum(abs,flux.qt),maximum(abs,flux.qz))
    areas=grid.areas
    return (;nt,nz,rpm,graded,pressure_relative_L2=sqrt(sum(areas.*(pressure-exactP).^2)/sum(areas.*exactP.^2)),
        fill_mean_absolute=sum(areas.*abs.(fill-exactF))/sum(areas),
        local_balance=maximum(abs,flux.balance-source)/flow_scale,
        minimum_pressure=minimum(pressure),minimum_fill=minimum(fill),maximum_fill=maximum(fill),
        complementarity=maximum(abs,(pressure/a.pressure_scale_pa).*(1 .-fill)),iterations=state.iterations)
end
