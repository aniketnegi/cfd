# Signed full-film verification

Check signed Reynolds pressure and load against Donaldson's published tables,
an independent spectral solution, and grid refinement. These cases have no
cavitation; negative pressures are retained.

Run from `journal-bearings/` with one native thread for repeatable comparisons:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. --check-bounds=yes simulation/verification/full-film/profile_checks.jl
python simulation/verification/full-film/postprocessing_checks.py
julia --project=. simulation/verification/full-film/run.jl out/verification/full-film-smoke
python simulation/verification/full-film/analyse.py out/verification/full-film-smoke
python simulation/verification/full-film/plot.py out/verification/full-film-smoke
```

Use a fresh output directory. Runs save inputs, source hashes, comparison tables
and profiles there. Plotting requires NumPy and Matplotlib.

The representative run covers `(L/D,epsilon)=(.7,.2),(1,.8),(.1,.99)` on both
models and three grids. Add `--full` after the output argument for all 320
Table 6.5 points, both models and three grids (1920 FV solves). Full-run plots
show the same three representative cases; `summary.json` includes every family.

## Models and independent reference

`reference.toml` identifies Donaldson's MIT1965 thesis and table location;
`table-6-5.csv` contains the printed load ratios. The classical benchmark-only
`DonaldsonCylinder` dispatch uses `h=c*(1+epsilon*cos(theta))` and one shaft
radius. The production comparison retains the circle-intersection gap and separate
journal/housing radii. Each finite load uses its own model's infinite-width
normalization. Neither case clips signed pressure or models cavitation.

`spectral.jl` independently assembles the classical sine-Galerkin operator and
solves its axial dependence analytically. It shares no FV matrix or flux
reconstruction. Its full-field accuracy claims apply to the classical model,
not the production circle-intersection geometry.

`donaldson_series.jl` reconstructs the source's converged-eigenpair then
finite-series procedure. In particular, a 15-term source approximation is not
an exact continuum target at high eccentricity. Comparing an independent
reference and a printed table tests different things; matching every printed
digit is not a solver acceptance condition.

## Field and geometry checks

```sh
julia --project=. simulation/verification/full-film/field_checks.jl out/verification/full-film-fields
julia --project=. simulation/verification/full-film/reference_profiles.jl out/verification/full-film-reference
julia --project=. simulation/verification/full-film/geometry_review.jl out/verification/full-film-geometry
python simulation/verification/full-film/mutation_checks.py out/verification/full-film-mutations
```

The field runner compares the `L/D=1,epsilon=.8` whole field, pressure, end
gradient and friction at 128/256/512 grids with the independent reference.
`--full` adds selected short/long high-eccentricity fields and finer local
runs. The reference-profile runner produces classical curves without FV.
The geometry review reports normal-gap and wall-speed approximations, inertia
scales and an isolated wall-speed sensitivity; it is not a thermal model.
Mutation checks run an unmodified control and deliberate defects in temporary
source copies, never in the working solver.

## Named conical case

`conical_verification.jl` runs signed full-film refinement for the saved
`conical-unported-e0p6` case. It requires an existing mesh and experiment manifest;
see [case setup and outputs](../../reynolds/docs/conical-formulation.md#commands).

## Refinement and interpretation

The axial base count is selected once at 128 angular cells, then both counts
are doubled. Physical profiles are interpolated at the actual requested
planes. Tests reject rounded non-self-similar families and out-of-span samples.
`postprocessing.py` rejects unsupported Richardson sequences instead of
silently emitting an order. Even an admissible estimate is not a rigorous
error bound or proof that a family is asymptotic.

The pressure force is on the housing sidewall. Integrated loads can converge
faster than local pressures/gradients; small residuals or matching source
hashes alone do not verify the physical equation. These are reduced-model
numerical checks, not experimental validation.
