# Reynolds solvers

Rigid, steady, isothermal, incompressible Newtonian lubrication models.
Geometry inputs use millimetres; solver quantities use SI units.

| Model                                        | Documentation                                         | Implementation           |
|----------------------------------------------|-------------------------------------------------------|--------------------------|
| Cylindrical/conical signed full film         | This page; [formulation](docs/conical-formulation.md) | `src/full_film.jl`           |
| Cylindrical/conical pressure-only cavitation | [Swift–Stieber](docs/swift-stieber.md)                | `src/swift_stieber.jl`       |
| 1D mass-conserving cavitation                | [Slider benchmark](docs/jfo-slider.md)                | `src/jfo_slider.jl`          |
| Periodic conical mass-conserving cavitation  | [JFO formulation](docs/conical-jfo.md)                | `src/conical_jfo.jl`         |
| Rigid floating-bush force balance            | [Equilibrium](docs/conical-equilibrium.md)            | `src/conical_equilibrium.jl` |

The full-film and Swift–Stieber solvers accept supplied unported tensor meshes,
including adjusted spacing, and share one pressure assembly. Count-based
full-film functions delegate to that assembly. JFO and equilibrium interfaces
construct uniform grids from counts; fitted meshes are unsupported.
See [solver support](../../docs/solver-support.md) for integration limits.
Thermal effects, elastic deformation and fluid–structure coupling are absent.

Implementation lives in `src/`; `run.jl` is the named-case entry point.
`src/full_film.jl` owns shared pressure assembly, fluxes and forces.
`src/cylindrical.jl` and `src/conical.jl` provide count-based convenience APIs.

## Commands

Run tests from `journal-bearings/`:

```sh
julia --project=. simulation/tests/runtests.jl
```

See [test selection](../tests/README.md) for individual suites.

The named `run.jl` CLI supports only unported, zero-angle cylinders. See the
[root README](../../README.org) for mesh/run commands. It verifies the saved
mesh checksum, requires rpm and dynamic viscosity for a new run, and saves
`run.toml`, `full-film.vtk` and `diagnostics.toml`. `--sweep` also writes a
convergence TSV. Saved inputs cannot change under the same run name; reruns
require `--replace`. A standalone mesh `--output FILE` has no experiment
manifest and cannot be loaded as a named simulation case. Use `-v` for timings
and full diagnostics, or `-q` for errors only.

Fixed-input cylindrical examples:

```sh
julia --project=. simulation/examples/cylindrical_convergence.jl
julia --project=. simulation/examples/cylindrical_full_film.jl
```

See [examples](../examples/README.md) for inputs and outputs, or
[verification](../verification/README.md) for reference comparisons.

## Supplied-mesh full-film API

```julia
include("bearing.jl")
include("mesh/2d/grid.jl")
include("simulation/reynolds/src/full_film.jl")

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
grid = housing_grid(p,128,64) # optional theta_nodes / z_nodes adjust spacing
result = solve_full_film(grid; rpm=2000., viscosity_pa_s=.0277,
    end_pressure_pa=(0.,0.))
```

These are numerical-example inputs, not an adopted physical operating point.
Zero cone angle uses the same solver. `assemble_full_film` exposes `A*p=b` and
face coefficients; `full_film_fluxes` reconstructs continuity from those faces.
Results retain the supplied grid and return pressure, normal gap, face flows,
housing pressure force and residual diagnostics. `eccentricity_angle_rad`
rotates the journal offset relative to the fixed mesh. Returned grids annotate
updated radial gaps with `eccentricity_angle_rad`; reusing a result grid with
a different angle refreshes those gaps without moving housing vertices. The
solve keyword still defaults to zero, regardless of the input grid's annotation.

For export, separately include `simulation/reynolds/src/pressure_vtk.jl` and call
`write_conical_pressure_vtk(path,result)` for pressure and normal gap.
The solver core loads no exporter, cavitation code or experimental mesher.

Angular coordinates are radians and axial coordinates are mm in the mesh;
the solver converts lengths to metres. Tensor metadata must match the supplied
vertices/connectivity. Fitted meshes reject rather than being reshaped or
silently replaced. The named CLI remains cylindrical and uniform-grid only.

## Cylindrical equations

The steady Reynolds fluxes follow
[San Andrés, *Static Load Performance of Plain Journal Bearings*, Notes 4,
Eqs. 4.1 and 4.3–4.4](https://rotorlab.tamu.edu/me626/Notes_pdf/Notes04%20Static%20Load%20JBs%2010.pdf):

```text
qθ = U h/2 - h³/(12 μ R) ∂p/∂θ,
qz = -h³/(12 μ) ∂p/∂z,
U = ω Rj,   ∂θ(qθ)/R + ∂z(qz) = 0.
```

The fluxes have units m²/s per unit transverse width. Circumferential faces
are periodic. Axial end pressures are zero gauge, imposed at the physical
boundary using half-cell distances. One shared flux enters neighbouring
cells with opposite signs. Julia SparseArrays/SuiteSparse solves the linear
system. No cavitation floor is imposed: negative pressure is an intentional
full-film diagnostic, not a physically admissible cavitating prediction.

The fixed example uses Rj=50 mm, L=100 mm, radial clearance 0.050 mm,
e/c=0.6 along +X and +2000 rpm. Position is prescribed, not solved from load.
Its illustrative nominal ISO VG 32 viscosity at 40 °C is ν=32 mm²/s;
separately assuming ρ=860 kg/m³ gives μ=0.02752 Pa·s. ISO grade does not
specify density or an actual oil sample's operating temperature. Density is
only used for this conversion, not as an incompressible solver input.

## Numerical verification

`full_film_checks.jl` checks supplied/adjusted meshes, the discrete energy
identity, paired face balance, stationary exact solutions, rotating eccentric
manufactured solutions for both geometries, and force integration. Smoothly
graded refinement checks pressure accuracy; abrupt spacing changes do not carry
a blanket second-order accuracy guarantee.

`checks.jl` covers concentric zero pressure, speed/viscosity scaling, axial
symmetry, signed pressure, linear residuals, boundary balance and the
San Andrés short-bearing asymptote at L/D=0.01. That asymptote is not an exact
finite-width solution.

Independent finite-width classical references, printed Donaldson table inputs,
whole-field comparisons and refinement-analysis checks are available in
[full-film verification](../verification/full-film/README.md).
The benchmark explicitly distinguishes classical cosine geometry from the
production circle-intersection model. Its generated tables and profiles belong
to each run directory, not the solver API.

Self-convergence is not an absolute error bound. Positive and negative signed
end flows can cancel; their absolute sum is not oil consumption. Regression
thresholds are numerical guards, not physical accuracy specifications, and the
API minimum of 4×2 cells is not a verified operating resolution.

No independent finite-width experimental pressure field is matched by this
cylindrical example. It is not a [Gangrade reproduction](../../docs/gangrade-reference.md).
Clipping signed pressure is not JFO. Cavitating bearing predictions require
explicit replenishment, incoming content and cavitation pressure; force-balance
predictions additionally require a loaded-body and restraint convention.
