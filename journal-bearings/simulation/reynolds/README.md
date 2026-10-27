# Reynolds solvers

Rigid, steady, isothermal, incompressible Newtonian lubrication models.
Geometry inputs use millimetres; solver quantities use SI units.

| Model | Documentation | Implementation |
|---|---|---|
| Cylindrical signed full film | This page | `cylindrical.jl` |
| Conical signed full film | [Formulation](conical-formulation.md) | `conical.jl` |
| 1D mass-conserving cavitation | [Slider benchmark](jfo-slider.md) | `jfo_slider.jl` |
| Periodic conical mass-conserving cavitation | [JFO formulation](conical-jfo.md) | `conical_jfo.jl` |
| Rigid floating-bush force balance | [Equilibrium](conical-equilibrium.md) | `conical_equilibrium.jl` |

The bearing solvers construct uniform grids from counts. They do not accept
supplied fitted or adjusted-spacing meshes. The pressure-only Swift–Stieber
implementation is in `expts/gangrade-figure6/reynolds_cavitation.jl`.
See [solver support](../../docs/solver-support.md) for integration limits.
Thermal effects, elastic deformation and fluid–structure coupling are absent.

## Commands

From `journal-bearings/`, run each test in its own process:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/conical_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/jfo_slider_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/conical_jfo_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/conical_equilibrium_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/run_checks.jl
julia --project=. --check-bounds=yes expts/gangrade-figure6/reynolds_checks.jl
```

The named `run.jl` CLI supports only unported, zero-angle cylinders. See the
[root README](../../README.org) for mesh/run commands. It verifies the saved
mesh checksum, requires rpm and dynamic viscosity for a new run, and saves
`run.toml`, `full-film.vtk` and `diagnostics.toml`. `--sweep` also writes a
convergence TSV. Saved inputs cannot change under the same run name; reruns
require `--replace`. A standalone mesh `--output FILE` has no experiment
manifest and cannot be loaded as a named simulation case.

Fixed-input cylindrical examples:

```sh
julia --project=. simulation/reynolds/convergence.jl
julia --project=. simulation/reynolds/cylindrical.jl
```

These write under `out/cases/cylindrical-unported-e0p6/`. The result VTK includes
its geometry, gap and signed gauge pressure and can be opened independently.
Conical and slider verification scripts are documented on their model pages.

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

`checks.jl` covers concentric zero pressure, speed/viscosity scaling, axial
symmetry, signed pressure, linear residuals, boundary balance and the
San Andrés short-bearing asymptote at L/D=0.01. That asymptote is not an exact
finite-width solution.

The five-grid L/D=1 example compares each coarse pressure with the average of
its four refined cells. `E_p` is the relative L2 mismatch. Force is the signed
transverse pressure integral on a nominal cylindrical journal, not physical
load capacity.

| θ×z cells | max p (MPa gauge) | signed Fy (N) | E_p from coarser grid |
|---|---:|---:|---:|
| 32×16 | 16.597653 | 112518.01 | — |
| 64×32 | 16.668570 | 111806.08 | 0.013143 |
| 128×64 | 16.654109 | 111627.17 | 0.003276 |
| 256×128 | 16.657060 | 111582.39 | 0.0008183 |
| 512×256 | 16.655752 | 111571.19 | 0.0002045 |

Approximately fourfold reduction indicates second-order self-convergence,
not an absolute error bound. At 512×256, min/max signed pressure is ±16.655752
MPa and relative global imbalance is about 2.2e-14. Positive and negative end
flows cancel; their absolute sum is a residual normalizer, not oil consumption.

Test guards require the second refinement mismatch to be <0.4 times the first
and <0.005, L1/global residuals <1e-11, and max-cell normalized residual <1e-10.
The short-bearing check uses 24×12, 48×24 and 96×48; its finest relative L2 error
is about 0.0004083 against a 0.001 threshold. These are numerical regression
guards, not physical accuracy specifications. The API minimum of 4×2 cells is
not a verified operating resolution.

No independent finite-width experimental pressure field is matched by this
cylindrical example. It is not a [Gangrade reproduction](../../docs/gangrade-reference.md).
Clipping signed pressure is not JFO. Cavitating bearing predictions require
explicit replenishment, incoming content and cavitation pressure; force-balance
predictions additionally require a loaded-body and restraint convention.
