# JFO verification

Check mass-conserving pressure and liquid-content calculations. Run from
`journal-bearings/`, using a fresh output directory for each study:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. simulation/verification/jfo/cylindrical/run.jl out/verification/jfo-cylinder
julia --project=. simulation/verification/jfo/slider/jfo_slider_verification.jl out/verification/jfo-slider
julia --project=. simulation/verification/jfo/conical/conical_jfo_verification.jl out/verification/jfo-cone
```

Use `--help` for options, `-v` for refinement diagnostics or `-q` for errors only.
The cylinder requires an explicit directory; slider/cone runners also offer
case-specific defaults. Pass a fresh directory for repeated studies.

## Cylinder

[`cylindrical/`](cylindrical/) checks finite, ungrooved, submerged cylinders:
load, attitude and gross outward exchange at two subambient cavity pressures.
It combines continuous manufactured fields, independent nodal pressure/content
finite differences and comparisons with 24 printed Jakobsson–Floberg1957 rows.
Use `--full` for all rows. The source is distinct from the atmospheric grooved
Swift–Stieber reference. Its README defines models, normalization and outputs.

## Slider

The Lengiewicz et al. sinusoidal slider compares pressure, liquid content,
rupture/reformation positions and flux against an independent constant-flux
reference. It writes `convergence.tsv`, `profiles.tsv`, `face-flux.tsv` and
`verification.toml`.

See the [benchmark formulation](../../reynolds/docs/jfo-slider.md) for published
inputs and normalization. `slider/jfo_slider_reference.jl` contains the reference;
`slider/jfo_slider_verification.jl` runs the refinement study.

## Cone

[Conical transport](conical/README.md) checks supplied tensors and actual7mm
fitted inlets, independent conical finite differences, curved manufactured
pressure/content, and boundary flow/force refinement. Run
`conical/transport_verification.jl DIRECTORY`, then `conical/analyse.py` and
`conical/plot.py`. Use `--full` to include512×256 fitted meshes.

The separate artificial-reservoir study checks an artificial reservoir case and a manufactured
pressure/content field at both rotation signs. It writes `refinement.tsv`,
`manufactured.tsv`, `fields.tsv`, `cavitation.vtk` and `verification.toml`.
The reservoir is a numerical test condition, not a reconstructed physical inlet.

See the [conical formulation](../../reynolds/docs/conical-jfo.md) for equations,
boundaries and error measures. `conical/cases.jl` defines the reservoir case;
`conical/conical_jfo_mms.jl` defines the manufactured solution.

Run the associated checks with:

```sh
julia --project=. simulation/tests/runtests.jl cylindrical_jfo_checks.jl jfo_slider_checks.jl conical_jfo_checks.jl fitted_jfo_checks.jl
```
