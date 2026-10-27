# Solver and mesh support

| Entry point | Supported domain | Model |
|---|---|---|
| `simulation/reynolds/run.jl` | Named unported cylinders | Signed full film |
| `solve_cylindrical_full_film` | Uniform periodic cylindrical grid | Signed full film |
| `solve_conical_full_film` | Uniform periodic cylindrical/conical grid | Signed full film |
| `solve_conical_jfo` | Uniform periodic cylindrical/conical grid; explicit reservoir cells or pressurized ends | Pressure and liquid content |
| `equilibrate_floating_bush` | Uniform conical grid with a cell-mask feed | Rigid pressure-load balance |
| `expts/gangrade-figure6/reynolds_cavitation.jl` | Uniform periodic grid; optional pinned pressure columns | Pressure-only Swift–Stieber obstacle |

The bearing solvers construct their own uniform grids from counts. They do not
consume supplied fitted meshes or the mesher's optional adjusted spacing.
Fitted cell counts are not `n_theta*n_axial`.

The mesh CLI supports standalone fitted VTK exports, but rejects named fitted
cases. `experiment.jl` accepts `controlled-fitted-land-v1` manifests with a
fixed geometry/count contract; it does not load the public mesher's
`fitted-land-direct-v1` layout. The named Reynolds runner rejects all fed and
nonzero-angle cases.

## Fitted simulation limitations

`expts/conical-fitted-load/solve.jl` is not usable with the public mesher:

- Its private module includes `surface_feed.jl` without the prerequisite
  `grid.jl`, leaving `SurfaceMeshQuality` and background helpers undefined.
- It accesses `mesh.level`, which the public direct-count mesh does not expose.
- The force integral in `expts/conical-fitted-inlet/operator.jl` and the
  experiment's VTK mapping hard-code a 10° cone and 50 mm axial midpoint.

The experimental fitted pressure operator uses least-squares gradients and
multipoint nonorthogonal corrections. Its matrix is nonsymmetric and can have
positive off-diagonals. The uniform-grid SPD/maximum-principle argument does not
apply to it. Surface-mesh acceptance is not verification of fitted fluid flow.

`expts/conical-fitted-inlet/legacy_meshing_regression.jl` contains fixed-mesh
force/flow references that are not verified for the public mesh. It is excluded
from `mesh/tests/runtests.jl`.

## Numerical guarantees and limits

With positive face conductances and Dirichlet anchoring, the uniform full-film
pressure matrix satisfies

```text
x' A x = sum_internal G_f (x_owner-x_neighbour)^2
       + sum_Dirichlet G_f x_owner^2 > 0  for x != 0.
```

It is symmetric positive definite. Shared face fluxes enter adjacent cells
with opposite signs, giving conservative internal cancellation. These algebraic
properties do not establish continuum accuracy; analytical and manufactured
solutions test consistency and refinement.

The pressure-only obstacle minimizes `0.5 P' A P - b' P` subject to `P >= 0`
for this SPD matrix. Its conditions are `P >= 0`, `A*P-b >= 0` and
`P*(A*P-b)=0`. It does not transport liquid content through cavities.

JFO checks pressure/content bounds, complementarity and reconstructed liquid
flux. Its active-set iteration has a finite cap and reports failure; no global
convergence or uniqueness guarantee covers arbitrary inputs.

Signed negative pressures are numerical diagnostics, not physical cavitation.
JFO needs an explicit supply/inventory specification. Load balance alone is not
mesh convergence or dynamic stability. Thermal effects and elasticity are not
implemented in these solvers.

See the [solver documentation](../simulation/reynolds/README.md) for checks.
Numerical verification is distinct from physical validation against matched
measurements; the [Gangrade comparison](gangrade-reference.md) is not a matched
validation case.
