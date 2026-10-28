# Independent continuum reference, not a call to the finite-volume solver.
# Lengiewicz, Wichrowski & Stupkiewicz (2014), Section 4.1 and author supplement
# https://bluebox.ippt.pan.pl/~sstupkie/files/bearing.pdf (version 1.3).
# Re-derive their constant-flux solution in SI using quadrature and bisection;
# do not copy the Mathematica symbolic antiderivative or use graph-read targets.
using LinearAlgebra

slider_gap(x) = 20e-6-5e-6*cos(2pi*x/0.125)

function slider_reference(; quadrature_order::Int=64)
    quadrature_order >= 8 || throw(ArgumentError("quadrature needs at least 8 points"))
    # Gauss-Legendre nodes/weights by the symmetric Jacobi matrix.
    beta = [k/sqrt(4k^2-1) for k in 1:quadrature_order-1]
    eig = eigen(SymTridiagonal(zeros(quadrature_order),beta))
    nodes,weights = eig.values,2 .* eig.vectors[1,:].^2
    integrate(f,a,b) = (b-a)/2 * sum(weights .* f.((a+b)/2 .+ (b-a)/2 .* nodes))
    function root(f,a,b)
        fa,fb = f(a),f(b)
        fa*fb < 0 || error("continuum reference root is not bracketed")
        for _ in 1:60
            m = (a+b)/2
            fm = f(m)
            iszero(fm) && return m
            if signbit(fm) == signbit(fa)
                a,fa = m,fm
            else
                b = m
            end
        end
        return (a+b)/2
    end
    L,mu,u,pb = 0.125,0.015,2.0,1e6
    # At rupture dp/dx=0, hence q=u*h(xr). Integrate from the left reservoir.
    rupture = root(0.0,L/2) do xr
        hs = slider_gap(xr)
        pb+12mu*u*integrate(x -> (slider_gap(x)-hs)/slider_gap(x)^3,-L/2,xr)
    end
    hs = slider_gap(rupture)
    gradient(x) = 12mu*u*(slider_gap(x)-hs)/slider_gap(x)^3
    reformation = root(x -> integrate(gradient,x,L/2)-pb,rupture,L/2)
    pressure(x) = x <= rupture ? pb+integrate(gradient,-L/2,x) :
                  x < reformation ? 0.0 : pb-integrate(gradient,x,L/2)
    fill(x) = rupture < x < reformation ? hs/slider_gap(x) : 1.0
    return (; pressure,fill,rupture_m=rupture,reformation_m=reformation,
             constant_flux_m2_s=u*hs,hstar_m=hs,peak_pressure_pa=pressure(-rupture),
             pressure_at_zero_pa=pressure(0.0),gradient)
end
