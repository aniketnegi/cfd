# Documentation

- [Surface meshing](../mesh/2d/README.md): API, units, acceptance, CLI and tests.
- [Reynolds solvers](../simulation/reynolds/README.md): models, commands and verification.
- [Swift–Stieber](../simulation/reynolds/swift-stieber.md): supplied-mesh
  pressure-only cavitation, boundary contracts and diagnostic limits.
- [Solver support](solver-support.md): mesh compatibility and integration limits.
- [Full-film verification](../expts/full-film-verification/README.md): independent
  classical reference, geometry comparisons and reproducible benchmark commands.
- [Gangrade reference](gangrade-reference.md): published inputs, boundary ambiguity
  and the scope of the comparison experiments.

Experimental mesh generators and diagnostics live in `expts/`; their local
READMEs describe their own inputs and limits. They are not public meshing APIs.
