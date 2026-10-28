using LinearAlgebra, Printf, Serialization, DelimitedFiles

# Independent nodal FD model: Pi=p*c^2/(mu*omega*R^2), w=z/R.
# -Pi_tt - 3(H_t/H)Pi_t - Pi_ww = -6H_t/H^3, H=1+eps*cos(t).
# The four physical rectangle edges have Pi=0. No production code is included.
function fd_coefficients(nt,nz;ecc=.8,beta=1.,forcing=nothing)
    dt,dw=2pi/nt,2beta/nz
    theta=collect(range(0.,2pi;length=nt+1))
    axial=collect(range(-beta,beta;length=nz+1))
    h=1 .+ ecc.*cos.(theta)
    drift=-3ecc.*sin.(theta)./h
    east=1/dt^2 .+ drift./(2dt)
    west=1/dt^2 .- drift./(2dt)
    minimum(east)>0 && minimum(west)>0 || error("FD stencil not monotone: refine theta")
    rhs=isnothing(forcing) ? [6ecc*sin(t)/(1+ecc*cos(t))^3 for t in theta,w in axial] :
        [forcing(t,w) for t in theta,w in axial]
    return (;nt,nz,dt,dw,theta,axial,h,east,west,rhs,az=1/dw^2,diag=2/dt^2+2/dw^2)
end

function doubled_initial(old)
    nx,ny=size(old);fine=zeros(2nx-1,2ny-1)
    fine[1:2:end,1:2:end]=old
    fine[2:2:end-1,1:2:end]=(old[1:end-1,:]+old[2:end,:])/2
    fine[1:2:end,2:2:end-1]=(old[:,1:end-1]+old[:,2:end])/2
    fine[2:2:end-1,2:2:end-1]=(old[1:end-1,1:end-1]+old[2:end,1:end-1]+
        old[1:end-1,2:end]+old[2:end,2:end])/4
    return fine
end

function fd_residual(c,p;obstacle=true)
    pscale=max(maximum(abs,p),1.);bscale=max(maximum(abs,c.rhs),1.)
    kkt,full_residual=0.,0.
    @inbounds for j in 2:c.nz,i in 2:c.nt
        residual=c.east[i]*(p[i,j]-p[i+1,j])+c.west[i]*(p[i,j]-p[i-1,j])+
            c.az*((p[i,j]-p[i,j+1])+(p[i,j]-p[i,j-1]))-c.rhs[i,j]
        violation=obstacle ? abs(min(p[i,j]/pscale,residual/bscale)) : abs(residual/bscale)
        kkt=max(kkt,violation)
        if !obstacle || p[i,j]>1e-8*pscale
            full_residual=max(full_residual,abs(residual)/bscale)
        end
    end
    return (;kkt,full_residual)
end

function solve_fd(nt,nz;ecc=.8,beta=1.,forcing=nothing,obstacle=true,
                  initial=nothing,relaxation=nothing,tolerance=1e-10,max_sweeps=10000,report_every=0)
    c=fd_coefficients(nt,nz;ecc,beta,forcing)
    p=isnothing(initial) ? zeros(nt+1,nz+1) : copy(initial)
    size(p)==(nt+1,nz+1) || error("initial pressure dimensions")
    p[1,:].=0;p[end,:].=0;p[:,1].=0;p[:,end].=0
    rho=(cos(pi/nt)/c.dt^2+cos(pi/nz)/c.dw^2)/(1/c.dt^2+1/c.dw^2)
    omega=isnothing(relaxation) ? 2/(1+sqrt(1-rho^2)) : relaxation
    0<omega<2 || error("SOR relaxation must lie between0and2")
    for sweep in 1:max_sweeps
        @inbounds for j in 2:nz,i in 2:nt
            residual=c.east[i]*(p[i,j]-p[i+1,j])+c.west[i]*(p[i,j]-p[i-1,j])+
                c.az*((p[i,j]-p[i,j+1])+(p[i,j]-p[i,j-1]))-c.rhs[i,j]
            updated=p[i,j]-omega*residual/c.diag
            p[i,j]=obstacle ? max(0.,updated) : updated
        end
        if sweep%20==0
            residual=fd_residual(c,p;obstacle)
            if report_every>0 && sweep%report_every==0
                println("sweep=$sweep residual=$residual");flush(stdout)
            end
            if residual.kkt<=tolerance && residual.full_residual<=tolerance
                return (;pressure=p,coefficients=c,sweeps=sweep,relaxation=omega,residual)
            end
        end
    end
    error("Independent projected SOR failed after $max_sweeps sweeps: $(fd_residual(c,p;obstacle))")
end

function end_gradients(p,dw,order)
    if order==1
        low=(p[:,2]-p[:,1])/dw
        high=(p[:,end]-p[:,end-1])/dw
    elseif order==2
        low=(-3p[:,1]+4p[:,2]-p[:,3])/(2dw)
        high=(3p[:,end]-4p[:,end-1]+p[:,end-2])/(2dw)
    elseif order==3
        low=(-11p[:,1]+18p[:,2]-9p[:,3]+2p[:,4])/(6dw)
        high=(11p[:,end]-18p[:,end-1]+9p[:,end-2]-2p[:,end-3])/(6dw)
    else
        error("end derivative order must be1,2or3")
    end
    return low,high
end
trapezoid(values,spacing)=spacing*(sum(values)-(first(values)+last(values))/2)
function fd_observables(result;beta=1.)
    p,c=result.pressure,result.coefficients
    flows=Float64[]
    low2,high2=end_gradients(p,c.dw,2)
    low_flux2,high_flux2=c.h.^3 .* low2/12,-c.h.^3 .* high2/12
    for order in 1:3
        low,high=end_gradients(p,c.dw,order)
        push!(flows,trapezoid(c.h.^3 .*(low-high)/12,c.dt))
    end
    integrated=[trapezoid(p[i,:],c.dw) for i in 1:c.nt+1]
    fx=trapezoid(integrated.*cos.(c.theta),c.dt)/(2beta)
    fy=trapezoid(integrated.*sin.(c.theta),c.dt)/(2beta)
    mid=p[:,c.nz÷2+1]
    last_full=findlast(>(1e-8*maximum(p)),mid)
    return (;Q0_first=flows[1],Q0_second=flows[2],Q0_third=flows[3],
        Q0_lower=trapezoid(low_flux2,c.dt),Q0_upper=trapezoid(high_flux2,c.dt),
        load_P0=hypot(fx,fy),attitude_deg=rad2deg(atan(abs(fy),abs(fx))),
        peak_midplane=maximum(mid),last_positive_midplane_deg=c.theta[last_full]*180/pi,
        minimum_lower_flow_density=minimum(low_flux2),
        negative_lower_flow_integral=trapezoid(min.(low_flux2,0.),c.dt))
end
