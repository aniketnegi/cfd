# Actual-plane piecewise-linear samples; never reinterpret fitted cells as a tensor.
function transport_slice_sampler(slice)
    order=sortperm(first.(slice.theta));left=first.(slice.theta)[order]
    function value(t)
        k=searchsortedlast(left,t)
        k==0 && return NaN
        segment=order[k];lo,hi=slice.theta[segment]
        lo<=t<=hi || return NaN
        p,q=slice.pressure[segment]
        return p+(q-p)*(t-lo)/(hi-lo)
    end
    return function(theta)
        for t in (theta,theta-2pi,theta+2pi)
            p=value(t);isfinite(p) && return p
        end
        return NaN
    end
end
function transport_profiles(result;n=1024)
    rows=NamedTuple[]
    for z in (.025,.05,.075)
        p=fitted_pressure_slice(result.grid,result.pressure_pa;z_m=z,
            boundary_pressure=(x,y,tag)->tag==:feed ? 5e5 : 0.)
        f=fitted_pressure_slice(result.grid,result.fill_fraction;z_m=z,
            boundary_pressure=(x,y,tag)->1.)
        psample,fsample=transport_slice_sampler(p),transport_slice_sampler(f)
        for i in 1:n
            theta=(i-.5)*2pi/n
            push!(rows,(;z_m=z,theta_rad=theta,pressure_pa=psample(theta),
                fill_fraction=fsample(theta)))
        end
    end
    return rows
end
function transport_profile_change(coarse,fine)
    a,b=transport_profiles(coarse),transport_profiles(fine)
    good=[i for i in eachindex(a) if isfinite(a[i].pressure_pa) && isfinite(b[i].pressure_pa)]
    p=[b[i].pressure_pa for i in good];f=[b[i].fill_fraction for i in good]
    return (;pressure_change=norm(p-[a[i].pressure_pa for i in good])/norm(p),
        fill_change=sum(abs,f-[a[i].fill_fraction for i in good])/length(good),
        common_land_samples=length(good))
end
