# Shared steady pressure/void complementarity kernel. Assemblies and physical
# supply contracts remain in their respective slider/conical solvers.
using LinearAlgebra, SparseArrays

function solve_jfo_active_set(A, B, b, pressure_scale;
                              max_iterations::Int=2000)
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    n = length(b)
    size(A) == size(B) == (n,n) || throw(ArgumentError("mixed-system size mismatch"))
    isfinite(pressure_scale) && pressure_scale > 0 || throw(ArgumentError("invalid pressure scale"))
    full = trues(n)
    for iteration in 1:max_iterations
        # P unknown, v=0 in full cells; P=0, v unknown in cavitated cells.
        M = A*spdiagm(0=>pressure_scale.*full) + B*spdiagm(0=>Float64.(.!full))
        y = M \ b
        all(isfinite,y) || error("non-finite pressure/content active-set solve")
        pressure = pressure_scale .* y .* full
        void = y .* .!full
        next_full = (pressure ./ pressure_scale .- void) .> 0
        next_full == full && return (; pressure,void,iterations=iteration)
        full = next_full
    end
    error("pressure/content active set failed after $max_iterations iterations; check supply/inventory and resolution")
end
