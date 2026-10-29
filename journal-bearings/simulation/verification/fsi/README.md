# Fixed-position bearing FSI pilots

Run one-way and two-way coupling on cylinders and 10° cones with an elastic
housing and either a rigid or elastic journal.
The actual 7 mm +Y midlength inlet is at 0.5 MPa gauge; ends/cavity are at zero gauge.
Nominal journal radius 50 mm, length 100 mm, clearance 50 μm, offset −Y,
2000 rpm and viscosity 0.0277 Pa·s. Housing wall 5 mm (radial), E = 97 GPa,
ν = 0.30, outer face fixed, axial ends free. These are declared project inputs, not an experimental FSI case.

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
julia --project=. simulation/tests/runtests.jl elasticity_checks.jl fsi_checks.jl
julia --project=. simulation/verification/fsi/run.jl out/verification/housing-fsi
python simulation/verification/fsi/analyse.py out/verification/housing-fsi
python simulation/verification/fsi/plot.py out/verification/housing-fsi

# Both bodies elastic; same CLI and a separate fresh output directory:
julia --project=. simulation/verification/fsi/both.jl out/verification/both-fsi
python simulation/verification/fsi/analyse.py out/verification/both-fsi
python simulation/verification/fsi/plot.py out/verification/both-fsi
```

The default is a local pilot: 128×64 background, two radial housing layers, ε = 0.6,
both geometries and both coupling modes. The pressure mesh is not a tensor
around the opening. Solid nodes/DOFs and actual pressure cells are reported.
The optional solid journal uses a 32×16 angular/axial mesh with four O-grid
radial layers and a square core. It has no bore/overhangs, E = 210 GPa,
ν = 0.30, and both end faces fully fixed relative to the nominal pose.
The opening's reservoir pressure loads the journal beneath the actual aperture.

`--full` adds 256×128 background with four housing layers and ε = 0.2/0.6/0.8.
The elastic journal then uses 64×32 angular/axial counts and eight radial layers.
This is
memory-heavier and is not an instruction to start an untested remote batch.
Benchmark a representative case on the target host first. Matching geometry
and material reuse the structural factorization across eccentricities; fluid
coefficients are recomputed on each deformed gap.

## Outputs

- `inputs.toml`: declared physical choices and numerical source fingerprints.
- `results.tsv`: mode, pressures, predicted/coupled gaps, normal displacement,
  energy, stress, liquid balance, coupling iterations, preparation/solve time
  and peak process RSS. Housing/journal counts, energies and displacements are
  separate; maximum stress/strain and elastic residual cover all active bodies.
  RSS is cumulative for the process, not per-case allocation.
- `coupling-*.tsv`: gap fixed-point defect, pressure change and actual structural
  residual on every iteration. The first pressure-change value is infinite
  because no previous pressure exists.
- Per-mode fluid/solid TSV and VTK fields and serialized results. Solid geometry
  uses SI coordinates with displacement vectors; fluid VTK geometry uses mm
  according to the existing exporter. Do not overlay them without unit conversion.
- `analysis.json`, coupling-history and field PNGs.

One-way pressure is the rigid pressure: its minimum-gap column is the predicted
elastic gap, not a second fluid solve. Two-way pressure and displacement belong
to the same converged gap iterate. Both are prescribed-position calculations;
forces are reactions, not solved journal-load equilibrium.

See [coupling](../../fsi/docs/steady.md) and
[elasticity](../../elasticity/docs/linear-elasticity.md) for equations and scope.
The independent solid patch/Navier/Lamé checks and established JFO verification
are separate from coupled spatial accuracy. Journal mesh refinement can be
studied independently of fluid/housing refinement. A successful pilot does
**not** establish mesh-independent FSI predictions. Contact, plasticity,
thermal effects, rotational inertia and centrifugal prestress are not included.
