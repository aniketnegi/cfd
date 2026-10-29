# Solver and mesh support

| Entry point                                     | Supported domain                                                                        | Model                                |
|-------------------------------------------------|-----------------------------------------------------------------------------------------|--------------------------------------|
| `simulation/reynolds/run.jl`                    | Named unported cylinders                                                                | Signed full film                     |
| `solve_full_film(grid; ...)` | Supplied cylindrical/conical unported tensor or fitted land grid; fitted rim pressure explicit | Signed full film |
| `solve_cylindrical_full_film`                   | Uniform periodic cylindrical grid                                                       | Signed full film                     |
| `solve_conical_full_film`                       | Uniform periodic cylindrical/conical grid                                               | Signed full film                     |
| `solve_conical_jfo(grid; ...)` | Supplied cylindrical/conical tensor or fitted land grid; tensor counts/masks or ends, fitted rim pressure explicit | Pressure and liquid content |
| `prepare_housing_fsi` / `solve_one_way_fsi` / `solve_two_way_fsi` | Supplied cylindrical/conical tensor or fitted land; fixed nominal journal position, radial housing wall and fixed outer support | Steady JFO and 3D linear-elastic housing; optional independently meshed solid journal via `add_elastic_journal` |
| `equilibrate_floating_bush`                     | Uniform conical grid with a cell-mask feed                                              | Rigid pressure-load balance          |
| `solve_swift_stieber(grid; ...)` | Supplied cylindrical/conical tensor or fitted land grid; fitted rim pressure explicit; optional floor-pressure cells | Pressure-only Swift–Stieber obstacle |

Full film, Swift–Stieber and JFO consume the parameters and tensor coordinates
retained by `housing_grid`. They also consume supplied `surface_feed_grid`
connectivity, with explicit `feed_pressure_pa` on rim faces.
Count-based full-film entry points construct uniform grids and delegate to the
same assembly. JFO accepts supplied uniform/adjusted tensors and fitted land connectivity;
equilibrium remains uniform-grid/count-based. Fitted
cell counts are not `n_theta*n_axial`; see [fitted full film](../simulation/reynolds/docs/fitted-full-film.md).

The mesh CLI supports standalone fitted VTK exports, but rejects named fitted
cases. `experiment.jl` accepts `controlled-fitted-land-v1` manifests with a
fixed geometry/count contract; it does not load the public mesher's
`fitted-land-direct-v1` layout. The named Reynolds runner rejects all fed and
nonzero-angle cases.

## Fitted simulation limitations

Maintained fitted support includes signed full film, pressure-only
Swift–Stieber and JFO pressure/content transport with an explicit above-cavity
rim pressure. Fitted load balance and named cases remain unsupported.
The corrected pressure matrix is not guaranteed symmetric or monotone; fitted
complementarity is not a tensor-style convex minimization.

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

With positive face conductances and Dirichlet anchoring, the tensor full-film
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
mesh convergence or dynamic stability. The standalone fluid solvers have no
structural equations; [housing FSI](../simulation/fsi/docs/steady.md) supplies a
separate linear-elastic model and relative two-body normal-gap feedback.
The supplied journal model is an end-clamped solid cylinder/frustum. It does
not establish coupled spatial convergence. Thermal effects remain
unsupported.

See the [solver documentation](../simulation/reynolds/README.md) for checks.
Numerical verification is distinct from physical validation against matched
measurements; the [Gangrade comparison](gangrade-reference.md) is not a matched
validation case.
