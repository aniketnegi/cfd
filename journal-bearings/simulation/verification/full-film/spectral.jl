using LinearAlgebra, DelimitedFiles, Printf

# Independent sine-Galerkin reference for (H^3 P_theta)_theta+H^3 P_ww=H_theta.
# H^3 has only four cosine coefficients. Axial dependence is solved analytically.
function spectral_basis(ecc, modes)
    coefficients=[1+1.5ecc^2,3ecc+.75ecc^3,1.5ecc^2,.25ecc^3]
    coefficient(k)=k<=3 ? coefficients[k+1]/(k==0 ? 1 : 2) : 0.
    mass=[coefficient(abs(m-n))-coefficient(m+n) for m in 1:modes,n in 1:modes]
    stiffness=[m*n*(coefficient(abs(m-n))+coefficient(m+n)) for m in 1:modes,n in 1:modes]
    eig=eigen(Symmetric(stiffness),Symmetric(mass))
    minimum(eig.values)>0 || error("positive spectrum")
    alpha=sqrt.(eig.values)
    amplitudes=ecc .* eig.vectors[1,:] ./ eig.values
    exact_infinite_first=ecc/((1+ecc^2/2)*sqrt(1-ecc^2))
    return (;vectors=eig.vectors,alpha,amplitudes,exact_infinite_first)
end
function spectral_load(basis,beta)
    factor=1 .- tanh.(basis.alpha.*beta)./(basis.alpha.*beta)
    return sum(basis.vectors[1,:].*basis.amplitudes.*factor)/basis.exact_infinite_first
end
# Stable cosh(alpha*w)/cosh(alpha*beta), including large alpha*beta.
cosh_ratio(alpha,w,beta)=exp(alpha*(abs(w)-beta))*(1+exp(-2alpha*abs(w)))/(1+exp(-2alpha*beta))
function spectral_pressure(basis,theta,w,beta)
    shape=[1-cosh_ratio(a,w,beta) for a in basis.alpha]
    coefficients=basis.vectors*(basis.amplitudes.*shape)
    return sum(coefficients[m]*sin(m*theta) for m in eachindex(coefficients))
end
function spectral_end_gradient(basis,theta,beta)
    coefficients=basis.vectors*(basis.amplitudes.*basis.alpha.*tanh.(basis.alpha.*beta))
    return sum(coefficients[m]*sin(m*theta) for m in eachindex(coefficients))
end
