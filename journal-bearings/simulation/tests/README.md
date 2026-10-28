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
- `checks.jl`, `conical_checks.jl`: cylindrical and conical full-film checks.
- `jfo_slider_checks.jl`, `conical_jfo_checks.jl`: liquid-content transport and references.
- `conical_equilibrium_checks.jl`: rigid position/load-balance checks.
- `run_checks.jl`: named-case interface, output protection and resolution sweeps.

For reproducible comparison tables and fields, see
[verification](../verification/README.md).
