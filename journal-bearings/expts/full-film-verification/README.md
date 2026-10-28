# Signed full-film verification

Run from `journal-bearings/` with one native thread for repeatable comparisons:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. --check-bounds=yes expts/full-film-verification/profile_checks.jl
python expts/full-film-verification/postprocessing_checks.py
julia --project=. expts/full-film-verification/run.jl out/verification/full-film-smoke
python expts/full-film-verification/analyse.py out/verification/full-film-smoke
python expts/full-film-verification/plot.py out/verification/full-film-smoke
```

Every solve reads the checked-out shared FV core. The output directory must
not exist. It receives input/source hashes, raw tables and requested physical
sampling coordinates; generated output stays outside the source directory.
No PDF, source snapshot, old review directory or saved binary field is needed.
Julia dependencies are in the root project; plotting uses NumPy/Matplotlib.

The representative run covers `(L/D,epsilon)=(.7,.2),(1,.8),(.1,.99)` on both
models and three grids. Add `--full` after the output argument for all320
Table6.5 points, both models and three grids (1920FV solves). Full-run plots
show the same three representative cases; summary.json includes every family.
The full sweep is an opt-in benchmark, not an ordinary unit test.

## Models and independent reference

`reference.toml` identifies Donaldson's MIT1965 thesis and table location;
`table-6-5.csv` contains the printed load ratios. The classical benchmark-only
`DonaldsonCylinder` dispatch uses `h=c*(1+epsilon*cos(theta))` and one shaft
radius. It does not modify production `BearingParameters` methods or defaults.
The production comparison retains the circle-intersection gap and separate
journal/housing radii. Each finite load uses its OWN model's infinite-width
normalization. Neither case clips signed pressure or models cavitation.

`spectral.jl` independently assembles the classical sine-Galerkin operator and
solves its axial dependence analytically. It shares no FV matrix or flux
reconstruction. Its full-field accuracy claims apply to the classical model,
not the production circle-intersection geometry.

`donaldson_series.jl` reconstructs the source's converged-eigenpair then
finite-series procedure. In particular, a15-term source approximation is not
an exact continuum target at high eccentricity. Comparing an independent
reference and a printed table tests different things; matching every printed
digit is not a solver acceptance condition.

## Field and geometry checks

```sh
julia --project=. expts/full-film-verification/field_checks.jl out/verification/full-film-fields
julia --project=. expts/full-film-verification/reference_profiles.jl out/verification/full-film-reference
julia --project=. expts/full-film-verification/geometry_review.jl out/verification/full-film-geometry
python expts/full-film-verification/mutation_checks.py out/verification/full-film-mutations
```

The field runner compares the `L/D=1,epsilon=.8` whole field, pressure, end
gradient and friction at128/256/512 grids with the independent reference.
`--full` adds selected short/long high-eccentricity fields and finer local
runs. The reference-profile runner produces classical curves without FV.
The geometry review reports normal-gap and wall-speed approximations, inertia
scales and an isolated wall-speed sensitivity; it is not a thermal model.
Mutation checks run an unmodified control and deliberate defects in temporary
source copies, never in the working solver.

## Refinement and interpretation

The axial base count is selected ONCE at128 angular cells, then both counts
are doubled. Physical profiles are interpolated at the actual requested
planes. Tests reject rounded non-self-similar families and out-of-span samples.
`postprocessing.py` rejects unsupported Richardson sequences instead of
silently emitting an order. Even an admissible estimate is not a rigorous
error bound or proof that a family is asymptotic.

The pressure force is on the housing sidewall. Integrated loads can converge
faster than local pressures/gradients; small residuals or matching source
hashes alone do not verify the physical equation. These are reduced-model
numerical checks, not experimental validation.
