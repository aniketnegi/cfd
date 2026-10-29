include(joinpath(@__DIR__,"../cases.jl"))
include(joinpath(@__DIR__,"../independent.jl"))

const XI=(;circumference_m=.08,length_m=.02,radius_m=.08/(2pi),clearance_m=17.5e-6,
    eccentricity_ratio=.8,mu_pa_s=.014,surface_speed_m_s=.2,boundary_absolute_pa=1e5,cavity_absolute_pa=0.)
const XI_RPM=XI.surface_speed_m_s/XI.radius_m*60/(2pi)
const XI_PS=XI.mu_pa_s*XI.surface_speed_m_s*XI.radius_m/XI.clearance_m^2

# Four actual pressure boundaries, not a finite-width reservoir-cell column.
# Internal production mobility/transport coefficients are retained unchanged.
function cut_angular_seam(a,pressure_above_cavity::Tuple{Real,Real})
    all(p->isfinite(p)&&p>0,pressure_above_cavity) || throw(ArgumentError("angular reservoirs must be above cavity pressure"))
    nt,nz=size(a.gt)
    A,B,b=copy(a.A),copy(a.B),copy(a.b)
    low,high=zeros(nz),zeros(nz)
    for j in 1:nz
        l,k=1+nt*(j-1),nt*j
        g,c=a.gt[nt,j],a.adv[nt,j]
        widths=a.axes.theta_widths
        low[j]=g*(widths[end]+widths[1])/widths[1]
        high[j]=g*(widths[end]+widths[1])/widths[end]
        A[k,l]+=g;A[l,k]+=g
        A[l,l]+=low[j]-g;A[k,k]+=high[j]-g
        b[l]+=low[j]*pressure_above_cavity[1]
        b[k]+=high[j]*pressure_above_cavity[2]
        # Incoming content is1; outgoing content uses the adjacent film cell.
        # The full-liquid RHS/diagonal stay unchanged, only the wrap coupling goes.
        B[c>=0 ? l : k,c>=0 ? k : l]=0.
    end
    dropzeros!(A);dropzeros!(B)
    return (;periodic=a,A,B,b,angular_pressure_pa=Float64.(pressure_above_cavity),low,high)
end

function cut_fluxes(cut,pressure,fill)
    a=cut.periodic;nt,nz=size(a.gt)
    native=conical_jfo_fluxes(a,pressure,fill)
    P,F=reshape(pressure,nt,nz),reshape(fill,nt,nz)
    qt=zeros(nt+1,nz);qt[2:nt,:]=native.qt[1:nt-1,:]
    for j in 1:nz
        c=a.adv[nt,j]
        qt[1,j]=cut.low[j]*(cut.angular_pressure_pa[1]-P[1,j])+c*(c>=0 ? 1. : F[1,j])
        qt[end,j]=cut.high[j]*(P[end,j]-cut.angular_pressure_pa[2])+c*(c>=0 ? F[end,j] : 1.)
    end
    qz=native.qz
    balance=vec([qt[i+1,j]-qt[i,j]+qz[i,j+1]-qz[i,j] for i in 1:nt,j in 1:nz])
    outward=[-qt[1,:];qt[end,:];-qz[:,1];qz[:,end]]
    return (;qt,qz,balance,outward)
end

function xi_assembly(nt,nz;model="classical",rpm=XI_RPM,graded=false)
    radius,clearance=1e3XI.radius_m,1e3XI.clearance_m
    p=if model=="classical"
        JFOCylinder(radius,radius,1e3XI.length_m,clearance,XI.eccentricity_ratio,0.,0.)
    elseif model=="production"
        BearingParameters(journal_midpoint_radius=radius,bearing_length=1e3XI.length_m,
            radial_clearance=clearance,eccentricity_ratio=XI.eccentricity_ratio,semicone_angle_deg=0.)
    else
        throw(ArgumentError("unknown geometry model"))
    end
    a=assemble_conical_jfo(cylinder_grid(p,nt,nz;graded);rpm,viscosity_pa_s=XI.mu_pa_s,
        end_pressure_pa=(XI.boundary_absolute_pa,XI.boundary_absolute_pa),cavitation_pressure_pa=XI.cavity_absolute_pa,
        eccentricity_angle_rad=pi)
    return cut_angular_seam(a,(XI.boundary_absolute_pa,XI.boundary_absolute_pa))
end
function xi_solve(nt,nz;kwargs...)
    cut=xi_assembly(nt,nz;kwargs...);a=cut.periodic
    state=solve_jfo_active_set(cut.A,cut.B,cut.b,a.pressure_scale_pa)
    pressure,fill=state.pressure,1 .-state.void
    flux=cut_fluxes(cut,pressure,fill)
    actual_scale=max(maximum(abs,a.adv),maximum(abs,flux.qt),maximum(abs,flux.qz))
    reservoir_scale=max(maximum(cut.low),maximum(cut.high),maximum(a.gz))*XI.boundary_absolute_pa
    scale=max(actual_scale,reservoir_scale)
    gross_out=sum(max.(flux.outward,0.));gross_in=-sum(min.(flux.outward,0.))
    stats=(;peak_pa=maximum(pressure),minimum_fill=minimum(fill),
        local_balance=maximum(abs,flux.balance)/scale,
        local_face_balance=maximum(abs,flux.balance)/max(actual_scale,1e-30),
        max_cell_imbalance_m3_s=maximum(abs,flux.balance),flow_scale_m3_s=scale,
        global_balance=abs(sum(flux.outward))/max(gross_out+gross_in,scale),
        complementarity=maximum(abs,(pressure/a.pressure_scale_pa).*(1 .-fill)),
        gross_out_m3_s=gross_out,gross_in_m3_s=gross_in,net_boundary_m3_s=sum(flux.outward),
        left_out_m3_s=-sum(flux.qt[1,:]),right_out_m3_s=sum(flux.qt[end,:]),
        lower_out_m3_s=-sum(flux.qz[:,1]),upper_out_m3_s=sum(flux.qz[:,end]),iterations=state.iterations)
    return (;grid=a.grid,pressure_pa=pressure,fill_fraction=fill,cut,flux,stats)
end

# The nodal FD seam is an exact Dirichlet node, not a FV cell reservoir.
# Negative e here is only a coordinate shift: H=1+.8cos(theta).
function xi_fd(nt,nz;relaxation=1.,tolerance=1e-10,max_sweeps=200000)
    aspect=XI.length_m/(2XI.radius_m)
    a=CylindricalJFOReference.assemble(aspect,-XI.eccentricity_ratio,-XI.boundary_absolute_pa/XI_PS,nt,nz)
    P=fill(-a.pv0,nt,nz+1);V=zeros(size(P));residual=Inf
    for sweep in 1:max_sweeps
        @inbounds for j in 2:nz,i in 2:nt
            l,r=i-1,mod1(i+1,nt)
            rhs=a.forcing[i,j-1]-a.left[i]*P[l,j]-a.right[i]*P[r,j]-
                a.axial*(P[i,j-1]+P[i,j+1])-a.vu[i]*V[a.up[i],j]
            candidate=rhs>=0 ? rhs/a.diag : -rhs/a.vd[i]
            x=(1-relaxation)*(P[i,j]-V[i,j])+relaxation*candidate
            P[i,j]=max(x,0.);V[i,j]=max(-x,0.)
        end
        if sweep%20==0
            yield();residual=0.
            @inbounds for j in 2:nz,i in 2:nt
                l,r=i-1,mod1(i+1,nt)
                pressure=(2P[i,j]-P[l,j]-P[r,j])/a.dt^2-
                    1.5a.Ht[i]/a.H[i]*(P[r,j]-P[l,j])/a.dt+(2P[i,j]-P[i,j-1]-P[i,j+1])/a.dz^2
                transport=6*(a.H[a.up[i]]*V[a.up[i],j]-a.H[i]*V[i,j])/(a.H[i]^3*a.dt)
                residual=max(residual,abs(pressure+transport-a.forcing[i,j-1])*a.H[i]^3)
            end
            residual*=a.dt/(6maximum(a.H))
            residual<tolerance && return (;a,pressure_pa=XI_PS.*P,fill_fraction=1 .-V,sweeps=sweep,residual)
        end
    end
    error("cut-domain FD did not converge: $residual")
end

# Continuous full-liquid manufactured field with the same four constant traces.
# Coupled cavity transport is checked by arbitrary-field balances and physical FD.
function xi_manufactured(nt,nz;graded=false,rpm=XI_RPM)
    cut=xi_assembly(nt,nz;graded,rpm);a=cut.periodic
    source=zeros(nt*nz);exact=similar(source)
    points=(-sqrt(3/5),0.,sqrt(3/5));weights=(5/9,8/9,5/9)
    amplitude=2e5;L=XI.length_m;R=XI.radius_m;mu=XI.mu_pa_s
    U=rpm*2pi/60*R
    for j in 1:nz,i in 1:nt
        k=i+nt*(j-1);t,z=a.axes.theta_centres[i],a.axes.z_centres[j]
        exact[k]=XI.boundary_absolute_pa+amplitude*sin(t/2)*sin(pi*z/L)
        dt,dz=a.axes.theta_widths[i],a.axes.z_widths[j]
        for (u,wu) in zip(points,weights),(v,wv) in zip(points,weights)
            theta=a.axes.theta_edges[i]+(u+1)*dt/2
            axial=a.axes.z_edges[j]+(v+1)*dz/2
            h=XI.clearance_m*(1+XI.eccentricity_ratio*cos(theta))
            ht=-XI.clearance_m*XI.eccentricity_ratio*sin(theta)
            pt=amplitude/2*cos(theta/2)*sin(pi*axial/L)
            ptt=-amplitude/4*sin(theta/2)*sin(pi*axial/L)
            pzz=-amplitude*(pi/L)^2*sin(theta/2)*sin(pi*axial/L)
            density=U*ht/2-(3h^2*ht*pt+h^3*ptt)/(12mu*R)-h^3*R*pzz/(12mu)
            source[k]+=wu*wv*dt*dz/4*density
        end
    end
    state=solve_jfo_active_set(cut.A,cut.B,cut.b+source,a.pressure_scale_pa)
    P,F=state.pressure,1 .-state.void;flux=cut_fluxes(cut,P,F)
    scale=max(maximum(abs,flux.qt),maximum(abs,flux.qz),maximum(abs,a.adv))
    return (;nt,nz,graded,rpm,pressure_relative_L2=norm(P-exact)/norm(exact),
        local_balance=maximum(abs,flux.balance-source)/scale,minimum_pressure=minimum(P),
        fill_error=maximum(abs,F.-1),iterations=state.iterations)
end

function xi_profile(r,x_m;field=:pressure_pa)
    a=r.cut.periodic;values=plane_profile(r.grid,getproperty(r,field);z_m=XI.length_m/2)
    x=a.axes.theta_centres.*XI.radius_m
    boundary=field==:pressure_pa ? XI.boundary_absolute_pa : 1.
    nodes=[0.;x;XI.circumference_m];samples=[boundary;values;boundary]
    return [begin
        0<=q<=XI.circumference_m || throw(ArgumentError("profile coordinate outside domain"))
        i=min(searchsortedlast(nodes,q),length(nodes)-1)
        f=(q-nodes[i])/(nodes[i+1]-nodes[i]);(1-f)*samples[i]+f*samples[i+1]
    end for q in x_m]
end
