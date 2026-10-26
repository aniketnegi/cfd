# Reynolds numerical-verification checkpoints

A separate **conical signed full-film** solver is now implemented and tested.
Read [conical-formulation.md](conical-formulation.md) for the dimensional
surface-flux derivation, approximation choices, pressure-force signs, tests
and saved 10° refinement results. Run its checks with:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_checks.jl
```

It is not yet a cavitating, port-fed, load-equilibrium or EHD solver. The
existing named `run.jl` interface remains cylindrical-only; the conical page
provides a case-specific verification script without expanding that CLI.

## Published cavitation checkpoint (1D, separate from bearing solvers)

A conservative pressure/liquid-content solver now passes the Lengiewicz et al.
sinusoidal-slider analytical benchmark, including rupture and reformation.
Read [jfo-slider.md](jfo-slider.md) for the primary-source audit (including the
125 mm length correction), formulation, 149 checks, six-grid results and plot.

```sh
julia --project=. --check-bounds=yes simulation/reynolds/jfo_slider_checks.jl
```

This verifies a 1D method, not the periodic conical assembly or physical VG32
supply. No bearing cavitation pressure or feed assumption was chosen by it.

## 2D periodic conical cavitation checkpoint

[conical-jfo.md](conical-jfo.md) documents the pressure/content extension,
171 checks, a curved 2D manufactured solution, and numerical reservoir-fed
refinement. The active-set kernel is shared with the published slider test.

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_jfo_checks.jl
```

Supply/source accounting is explicit, including reservoir return flow. The
saved broad patch is not the rig's port; physical feed/cavitation parameters,
load equilibrium and elasticity remain separate decisions. The named CLI is
unchanged, and previous runs were preserved.

## Cylindrical checkpoint (retained unchanged)

From `journal-bearings/`, verify the solver and the named-case interface:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/run_checks.jl
julia --project=. --check-bounds=yes test/mesh_cli.jl
```

See [`README.org`](../../README.org) for mesh and run CLI examples. The
historical fixed-input example commands remain available:

```sh
julia --project=. simulation/reynolds/convergence.jl
julia --project=. simulation/reynolds/cylindrical.jl
```

The named mesh CLI writes `experiment.toml` with geometry, grid, feed variant,
and a SHA-256 of `mesh.vtk`. The simulation CLI reads and verifies it; a new
run needs an explicit rpm and **dynamic** viscosity (Pa·s). It saves
`run.toml` (including fixed zero-gauge axial pressure, periodic θ, no
cavitation), a result VTK, and `diagnostics.toml`. Use `--sweep` to also write
a convergence TSV. For the **same** run name, omitted operating flags are
reloaded from `run.toml`; differing values are rejected. Reruns require
`--replace`. A new geometry requires a new experiment name. The runner
currently rejects fed and nonzero-angle cases rather than silently applying
the cylindrical unported model. The old `--output FILE` mesh mode produces
only a standalone VTK and cannot serve as a saved simulation experiment.

All maintained outputs for the illustrative case live under
`out/cases/cylindrical-unported-e0p6/`: `mesh.vtk` contains geometry and gap
fields; `runs/2000rpm-vg32-40c/full-film.vtk` contains the same mesh plus
signed gauge pressure; and `runs/2000rpm-vg32-40c/convergence.tsv` contains
the grid sweep. A future plot can live alongside its run. The full-film VTK
intentionally contains its own geometry so it can be opened independently.
Both VTK files use the periodic cell layout from `mesh/2d/conical/grid.jl` at
**zero** semi-cone angle. The two conical mesh-only examples have their own
`out/cases/conical-*/mesh.vtk`. Older OpenFOAM cases and the large
mesh-research backup retain their previous paths; these were not migrated.
This is a cylindrical numerical-verification case, **not** a conical result,
a JFO solver, or a physically admissible bearing prediction.

## Model and sources

The control-volume form of the steady, isothermal, incompressible Newtonian
Reynolds equation and its periodic-θ, ambient-pressure axial-end conditions
are taken from [San Andrés, *Static Load Performance of Plain Journal
Bearings*, Notes 4, Eq. 4.1 and 4.3–4.4](https://rotorlab.tamu.edu/me626/Notes_pdf/Notes04%20Static%20Load%20JBs%2010.pdf).
Face mass-volume fluxes for a cylinder are

```text
qθ = U h/2 - h³/(12 μ R) ∂p/∂θ,   qz = -h³/(12 μ) ∂p/∂z,
U = ω Rj,   ∂θ(qθ)/R + ∂z(qz) = 0.
```

`qθ`, `qz` are volume flux per unit transverse length (m²/s). One shared flux
per internal face enforces cellwise conservation; circumferential faces wrap
at the seam. The axial boundary flux uses the **half-cell** distance to
p=0 gauge. `SparseArrays`/Julia's SuiteSparse solve the linear system; no
third-party Reynolds code was copied or adapted because the screened codes
in [`probe/reynolds-candidates.md`](../../../../probe/reynolds-candidates.md)
do not provide the matching periodic model. The bearing-specific flux
assembly is kept here for inspection. No cavitation pressure floor is imposed:
negative gauge pressures are an expected *diagnostic* of this full-film model.

The run prescribes `Rj=50 mm`, `L=100 mm`, clearance `0.050 mm`, journal
position `e/C=0.6` along +X and `2000 rpm` in increasing θ. This does **not**
solve for the journal position from an applied load. Its oil is an illustrative
**nominal ISO VG 32 at 40 °C**: [ISO 3448](https://cdn.standards.iteh.ai/samples/8774/3288791cc12a452ea9d1a8cf94dacf56/ISO-3448-1992.pdf)
specifies ν=32 mm²/s nominal (range 28.8–35.2 at 40 °C). We separately
assume ρ=860 kg/m³ from Gangrade's reference properties, giving
μ=νρ=0.02752 Pa·s. Neither the density nor the temperature of a particular
oil sample is fixed by its ISO VG grade. Density only converts ν to μ here;
it is **not** an independently required incompressible Reynolds input.

## Evidence, not validation

`checks.jl` tests zero-eccentricity, speed/viscosity scaling, axial symmetry,
signed pressure, local linear-system residual, net boundary volume flux and
convergence to San Andrés' **short-bearing asymptote** at `L/D=0.01`.
The asymptote is not an exact finite-L/D result, and no comparison to an
independent finite-width published pressure field or experiment is claimed.

The reproducible sweep (`convergence.jl`, Julia 1.13.2) uses five successively
doubled grids for `L/D=1`. `E_p` compares the pressure in each coarse cell
to the **average of its four refined cells**, then takes an L2 norm relative
to the norm of the refined-cell averages. Force shown is the integrated
*transverse component of signed pressure on a nominal cylindrical journal*,
not a physical load capacity (negative pressures make it inadmissible):

| θ×z cells | max p (MPa gauge) | signed Fy (N) | sum abs axial boundary flux (m³/s) | E_p from previous grid |
| --- | ---: | ---: | ---: | ---: |
| 32×16 | 16.597653 | 112518.01 | 4.794668e-5 | — |
| 64×32 | 16.668570 | 111806.08 | 4.802074e-5 | 0.013143 |
| 128×64 | 16.654109 | 111627.17 | 4.803941e-5 | 0.003276 |
| 256×128 | 16.657060 | 111582.39 | 4.804408e-5 | 0.0008183 |
| 512×256 | 16.655752 | 111571.19 | 4.804525e-5 | 0.0002045 |

Differences fall by about 4× per refinement. This is **observed
self-convergence**, not a finite-width analytical error bound. At 512×256:
min/max signed gauge pressure = ±16.655752 MPa; absolute net signed axial
outflow = `1.0644e-18 m³/s`; sum absolute axial boundary contributions =
`4.804525e-5 m³/s`; relative global imbalance = `2.2154e-14`.
The max absolute cell residual is `1.1563e-20 m³/s` and its ratio to the
*mean* circumferential advective-face flux magnitude is `1.1305e-13`.
Summed absolute cell residual is `8.6092e-17 m³/s`, divided by the sum of
advective-face flux magnitudes = `6.4222e-15`. These check discrete algebraic
conservation, **not** physical supply or pressure accuracy. Signed positive
and negative boundary flows cancel: the sum of absolute contributions is a
normalization, **not** physical oil consumption.

The tests start the finite-width grid sequence at **32×16**, check 32×16 →
64×32 and 64×32 → 128×64 against `E_p(second) < 0.4 E_p(first)` and
`E_p(second) < 0.005`. They check L1 residual/global imbalance < `1e-11`
and max-cell normalized residual < `1e-10` on 32×16. These are test guards,
**not** physical accuracy criteria. There is no pressure cutoff in this
full-film run. The solver interface permits ≥4×2 cells, but that minimum is
not a verified resolution. For the analytical short-bearing limit, tests use
`L=1 mm` (`L/D=0.01`) and 24×12, 48×24, 96×48; relative L2 errors against
the leading-order parabola are approximately 0.006889, 0.001679, 0.0004083;
the final test threshold is 0.001. This is a limiting formula, not a
finite-L/D pressure oracle; leading-order gap and discrete resolution leave
a nonzero floor. No independent finite-width published pressure field or
experiment has been matched.

**Not a Gangrade reproduction:** the reference run borrows the reported
50 mm midpoint journal radius, 100 mm length, 50 µm clearance, the Fig. 8
2000 rpm speed, and the table's 860 kg/m³ density. It changes semi-cone angle
**10° → 0°**, prescribes illustrative offset e=30 µm (ε=0.6), omits the
paper's listed 0.5 MPa supply and feed feature, uses nominal ISO VG 32
μ=0.02752 rather than the table's 0.0277 Pa·s, does not apply the Fig. 8
1000 N radial/200 N axial loads, and deliberately leaves negative pressures
unclipped. A match to Gangrade's oil temperature, cavitation and pressure
references has **not** been demonstrated.

## Current limits and next checkpoints

The [conical source review](../../docs/conical-solver-review.md) and
[independent 1D cavitation benchmark](jfo-slider.md) are complete within their
stated scopes. The [2D conical JFO extension](conical-jfo.md) now also checks
pressure/content transport and supply accounting. The subsequent
[Fig.6 comparison](../../docs/gangrade-figure6-comparison.md) and
[Reynolds-versus-JFO control](../../docs/gangrade-cavitation-comparison.md)
remain unmatched at high eccentricity. A subsequent
[60-node Q4 reconstruction](../../docs/gangrade-q4-grid-comparison.md) demonstrates
a large coarse-mesh effect; refined FE/FV pressures agree within 0.024%.
Resolve the remaining modelling/implementation
mismatch before using the paper as validation; load equilibrium and the physical
baseline remain unfinished. Clipping negative pressure from either
full-film output is **not** JFO. The signed unported
run's ambient ends do not establish a physical replenishment mechanism after
cavitation is added. Specify supply, incoming content and cavitation pressure
together before physical comparisons, then review load equilibrium and EHD.

The separate [floating-bush JFO driver](conical-equilibrium.md) now implements
translated geometry and pressure-load balance for an explicitly chosen local
feed/axial-freedom hypothesis. Fixed-seating radial equilibrium succeeds;
simultaneous 1000 N radial / 200 N axial equilibrium is not established.
87 new assertions and all 478 earlier solver assertions pass. See that report
for the unresolved axial search, feed refinement and force-body conventions.
