# Simulation

Steady Reynolds lubrication models for rigid journal bearings. Start with
[Reynolds solvers](reynolds/README.md) for supported models, inputs, units and
Julia APIs. The [named-case example](../README.org#named-cylindrical-case)
shows the mesh-to-simulation workflow; that command-line entry point supports
unported cylinders only.

To assess numerical accuracy, use the [verification suites](verification/README.md):
published comparisons, independent reference solvers and refinement studies.
[Fixed-input examples](examples/README.md) demonstrate pressure export and resolution sweeps.
[Automated tests](tests/README.md) cover solver contracts, numerical checks and
the named-case interface.

Run commands from `journal-bearings/` with `julia --project=.`. Generated results
belong under `out/`, not alongside source. See [solver support](../docs/solver-support.md)
for mesh compatibility and physical-model limits.
