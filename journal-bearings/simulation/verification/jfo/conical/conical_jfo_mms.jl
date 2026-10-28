# Independent continuum forcing for a 2D conical pressure/content manufactured
# solution. Artificial distributed injection/removal is TEST DATA, not a feed
# option in solve_conical_jfo. Include conical_jfo.jl before this file.
function conical_jfo_mms(nt,nz; rpm=2000.0)
    p = BearingParameters(journal_midpoint_radius=50.0,bearing_length=100.0,
        radial_clearance=0.05,semicone_angle_deg=10.0,eccentricity_ratio=0.6)
    mu,L,c,e,ca,k,P0 = 0.02752,0.1,50e-6,30e-6,cosd(10.0),tand(10.0),1e7
    a = assemble_conical_jfo(p,nt,nz;rpm,viscosity_pa_s=mu,
        end_pressure_pa=(0.0,0.0),cavitation_pressure_pa=0.0)
    dt,dz,omega = 2pi/nt,L/nz,rpm*2pi/60
    exactP,exactF,source = zeros(nt,nz),ones(nt,nz),zeros(nt,nz)
    for j in 1:nz,i in 1:nt
        theta,z = (i-0.5)*dt,(j-0.5)*dz
        st,ct = sincos(theta)
        rj = 0.05+(L/2-z)*k
        r = rj+c
        D = sqrt(rj^2-e^2*st^2)
        h = ca*(c-e*ct+e^2*st^2/(rj+D))
        ht = ca*(e*st+e^2*st*ct/D)
        hz = ca*k*(rj/D-1)
        K,Kt,Kz = h^3/(12mu),h^2*ht/(4mu),h^2*hz/(4mu)
        # Curved interfaces t(theta,z)=0; P and void have disjoint supports.
        t = ct+0.25*cos(2pi*z/L)
        tt,ttt = -st,-ct
        tz,tzz = -0.25*(2pi/L)*sin(2pi*z/L),-0.25*(2pi/L)^2*cos(2pi*z/L)
        f,df,ddf = z*(L-z)/L^2,(L-2z)/L^2,-2/L^2
        pt,ptt,pz,pzz,Ft = 0.0,0.0,0.0,0.0,0.0
        if t > 0
            exactP[i,j] = P0*f*t^4
            pt = P0*f*4t^3*tt
            ptt = P0*f*(12t^2*tt^2+4t^3*ttt)
            pz = P0*(df*t^4+f*4t^3*tz)
            pzz = P0*(ddf*t^4+8df*t^3*tz+f*(12t^2*tz^2+4t^3*tzz))
        else
            exactF[i,j] = 1-0.2*(-t)^4
            Ft = 0.8*(-t)^3*tt
        end
        diffusion = (Kt*pt+K*ptt)/(r*ca)+ca*((Kz*r-K*k)*pz+K*r*pzz)
        transport = omega*rj*(ht*exactF[i,j]+h*Ft)/(2ca)
        source[i,j] = (transport-diffusion)*dt*dz
    end
    state = solve_jfo_active_set(a.A,a.B,a.b+vec(source),a.pressure_scale_pa)
    P,F = state.pressure,1 .- state.void
    flux = conical_jfo_fluxes(a,P,F)
    return (; nt,nz,rpm,pressure_error=norm(P-vec(exactP))/norm(exactP),
        fill_error=sum(abs,F-vec(exactF))/(nt*nz),minimum_fill=minimum(F),
        minimum_pressure_pa=minimum(P),maximum_fill=maximum(F),
        complementarity=maximum(abs,(P ./ a.pressure_scale_pa).*(1 .- F)),
        source_balance_error=maximum(abs,flux.balance-vec(source))/maximum(abs,a.adv),
        iterations=state.iterations)
end
