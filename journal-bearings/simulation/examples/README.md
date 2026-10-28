# Simulation examples

Fixed-input signed full-film calculations for an unported cylinder. Run from
`journal-bearings/`:

```sh
julia --project=. simulation/examples/cylindrical_full_film.jl
julia --project=. simulation/examples/cylindrical_convergence.jl
```

Both use a 50 mm journal radius, 100 mm length, 50 micrometre clearance,
eccentricity ratio 0.6, 2000 rpm and viscosity 0.02752 Pa s.
The first exports a 512 by 256 pressure field; the second sweeps five grids
from 32 by 16 to 512 by 256.

Results go to `out/cases/cylindrical-unported-e0p6/`, with pressure and convergence
files under `runs/2000rpm-vg32-40c/`. Repeating these examples overwrites their
outputs. Negative pressure is retained; these examples do not model cavitation.

For configurable named runs, see the [project example](../../README.org).
For direct Julia calls, see the [solver API](../reynolds/README.md).
