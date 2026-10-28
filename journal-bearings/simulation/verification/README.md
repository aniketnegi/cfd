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

For solver usage, see [Reynolds models](../reynolds/README.md).
For automated checks, see [tests](../tests/README.md).
