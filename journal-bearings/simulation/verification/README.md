# Numerical verification

Comparisons against analytical solutions, independent reference solvers and
grid refinement:

- [Full film](full-film/README.md): Donaldson load tables, spectral pressure
  fields and conical geometry checks.
- [Swift–Stieber](swift-stieber/README.md): independent finite differences,
  end flow and manufactured conical solutions.
- [JFO](jfo/README.md): a published slider benchmark and conical
  pressure/liquid-content refinement.

Run commands from `journal-bearings/` after the [project setup](../../README.org).
Each suite lists its inputs, commands and comparison limits. Use fresh output
directories under `out/verification/`.

## Start here

- Full film: `full-film/run.jl DIRECTORY`, then its `analyse.py` and `plot.py`.
- Swift–Stieber: `swift-stieber/run_cylindrical.jl DIRECTORY` or
  `run_conical.jl DIRECTORY`, then that suite's `analyse.py` and `plot.py`.
- JFO: the [slider and cone runners](jfo/README.md) save their own verification
  tables and manifests; the other suites' postprocessors do not consume them.

Julia study runners accept `--help`, `-v` and `-q`. Where supported, `--full`
selects the larger study. Python postprocessors accept the same reporting flags.

## Independent comparisons and specialist studies

Use Swift–Stieber's `compare_fd.jl` for the independent finite-width comparison.
The [full-film specialist studies](full-film/README.md#specialist-studies) cover
field/reference profiles, geometry sensitivity, mutation checks and a saved
conical case. The [Swift–Stieber specialist studies](swift-stieber/README.md#specialist-studies)
cover manufactured fields, end-flow reconstruction and analytical controls.
These are separate investigations, not extra steps required after every run.

For solver usage, see [Reynolds models](../reynolds/README.md).
For automated solver checks, see [tests](../tests/README.md). Check these command
interfaces with:

```sh
python simulation/verification/interface_checks.py
```
