include(joinpath(@__DIR__,"spectral.jl"))

# Donaldson uses converged INFINITE-operator eigenpairs, not eigenpairs of an
# NxN truncation (section5.3.1,p134). Retain their first N Fourier coefficients
# and first N eigenfunctions, then enforce N boundary coefficients (Eq6.2-9).
# The independent symmetric Galerkin operator supplies converged eigenpairs;
# its weighted eigenvector normalization cancels in the boundary solve.
function donaldson_basis(ecc,modes)
    full=spectral_basis(ecc,max(320,4modes))
    vectors=full.vectors[1:modes,1:modes]
    root=sqrt(1-ecc^2)
    infinite_first=ecc/((1+ecc^2/2)*root)
    coefficients=[infinite_first*(-ecc)^(m-1)*(m+root)/(1+root)^m for m in 1:modes]
    return (;alpha=full.alpha[1:modes],vectors,
        amplitudes=vectors\coefficients,infinite_first)
end
function series_observables(ecc,basis,beta;theta=deg2rad(160.))
    a=basis.alpha
    first=basis.vectors[1,:].*basis.amplitudes
    ratio=1-sum(first.*tanh.(a.*beta)./(a.*beta))/basis.infinite_first
    correction=basis.vectors*(basis.amplitudes./cosh.(a.*beta))
    infinite_pressure=ecc/(1+ecc^2/2)*(1+ecc*cos(theta)/2)*sin(theta)/(1+ecc*cos(theta))^2
    pressure=infinite_pressure-sum(correction[m]*sin(m*theta) for m in eachindex(correction))
    end_correction=basis.vectors*basis.amplitudes
    end_pressure=infinite_pressure-sum(end_correction[m]*sin(m*theta) for m in eachindex(end_correction))
    derivative=basis.vectors*(basis.amplitudes.*a.*tanh.(a.*beta))
    gradient=sum(derivative[m]*sin(m*theta) for m in eachindex(derivative))
    friction=(2+ecc^2)/(6ecc*ratio)+ecc/2
    return (;ratio,pressure,end_pressure,gradient,friction)
end
