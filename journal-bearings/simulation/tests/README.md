# Simulation tests

<!--  TODO: Too vague -->
Checks for Reynolds solvers, numerical accuracy and the named-case interface.
<!-- TODO: very inconvenient to be running from another directory always isn't it? -->
Run all suites from `journal-bearings/`:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. simulation/tests/runtests.jl
```

For a focused run, pass filenames:

<!-- TODO: this is a very bad interface! -->
```sh
julia --project=. simulation/tests/runtests.jl full_film_checks.jl run_checks.jl
```

Use `-v` for detailed output, `-q` for failures only, or `--help` for usage.

## Suites

<!-- TODO: bit more detail here would be very nice i feel -->
- `full_film_checks.jl`, `swift_stieber_checks.jl`: supplied-mesh solver contracts.
- `fitted_full_film_checks.jl`: fitted pressure/flux accuracy, inlet boundaries,
  curved-surface forces and the common face operator.
- `fitted_swift_stieber_checks.jl`: fitted obstacle/reference contracts,
  nonlinear algebra, pressure-reference shifts and fully-liquid recovery.
- `pressure_sampling_checks.jl`: actual-plane reconstruction on fitted land,
  including periodic seams, inlet cuts and analytical pressure profiles.
- `checks.jl`, `conical_checks.jl`: cylindrical and conical full-film checks.
- `jfo_slider_checks.jl`, `conical_jfo_checks.jl`: liquid-content transport and references.
- `cylindrical_jfo_checks.jl`: supplied tensor JFO, submerged end exchange,
  analytical limits, reference/scaling invariance and actual-coordinate sampling.
- `elasticity_checks.jl`: 3D solid patch, energy, continuum manufactured forcing
  and analytical restrained-cylinder displacement.
- `journal_checks.jl`: solid O-grid volume, independent journal interpolation,
  land/opening pressure transfer and analytical solid-cylinder compression.
- `both_fsi_checks.jl`: relative two-body gap and fixed-position elastic journal/
  housing coupling.
- `fsi_checks.jl`: normal-gap projection, one-/two-way fixed-position coupling,
  actual drilled housing geometry and liquid/elastic residuals.
- `fitted_jfo_checks.jl`: tensor/face transport equivalence, actual inlet-rim
  liquid transport, fitted mixed-algebra oracles and continuous conical fields.
- `xi2024_jfo_checks.jl`: four-edge supply boundary adapter, manufactured
  pressure, reflection and independent cut-domain algebra.
- `conical_equilibrium_checks.jl`: rigid position/load-balance checks.
- `run_checks.jl`: named-case interface, output protection and resolution sweeps.

For reproducible comparison tables and fields, see
[verification](../verification/README.md).
