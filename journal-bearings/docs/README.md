# Documentation

- [Surface meshing](../mesh/README.md): generation, examples and tests.
- [Reynolds solvers](../simulation/reynolds/README.md): models, commands and verification.
- [Swift–Stieber](../simulation/reynolds/docs/swift-stieber.md): supplied-mesh
  pressure-only cavitation, boundary contracts and diagnostic limits.
- [Solver support](solver-support.md): mesh compatibility and integration limits.
- [Full-film verification](../simulation/verification/full-film/README.md): independent
  classical reference, geometry comparisons and reproducible benchmark commands.
- [Swift–Stieber verification](../simulation/verification/swift-stieber/README.md):
  independent finite-width pressure/flux checks and conical manufactured cases.
- [JFO verification](../simulation/verification/jfo/README.md): cylindrical source
  tables, slider reference and conical pressure/liquid-content studies.
- [Gangrade reference](gangrade-reference.md): published inputs, boundary ambiguity
  and the scope of the comparison experiments.
- [Feed-hole sizing](feed-hole-sizing.md): manufacturer/literature dimensions,
  the standard-informed 7 mm example and its pressure/liquid-supply specification.

Experimental mesh generators and diagnostics live in `expts/`; their local
READMEs describe their own inputs and limits. They are not public meshing APIs.
