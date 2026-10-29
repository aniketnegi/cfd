# Conical JFO transport

Verify prescribed-position pressure and liquid content on tensor cones and
land meshes around an actual circular inlet. The pressure mobility, gap and
wall-speed laws are shared with signed full film; the pressure/void kernel is
shared with cylinder JFO. This is not load equilibrium or FSI.

## Run

From `journal-bearings/`, use a fresh output directory:

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
julia --project=. simulation/verification/jfo/conical/transport_checks.jl
julia --project=. simulation/verification/jfo/conical/transport_verification.jl out/verification/jfo-conical
python simulation/verification/jfo/conical/analyse.py out/verification/jfo-conical
python simulation/verification/jfo/conical/plot.py out/verification/jfo-conical
```

`--full` adds512×256 physical grids and fitted manufactured fields. The default
uses64×32→128×64→256×128 tensors and128×64→256×128 fitted meshes. Mesh quality
criteria and optimizer settings are unchanged. The study is self-contained;
no saved binary meshes or historical review directories are needed.

## Inputs and separate supply problems

Both use journal midpoint radius50mm, axial length100mm, radial clearance50μm,
viscosity.0277Pa·s and prescribed eccentricity directed−Y. No physical input is
fitted to a historical pressure curve.

- **Unported tensors:** angles0/10/20°, eccentricity.6/.8,2000rpm, periodic
  circumference and both ends at.3MPa above the cavity. This is an explicit
  numerical end-supply fixture, not the drilled-inlet design.
- **Fitted land:** angles0/10°, eccentricity.2/.6/.8,±2000rpm,7mm midlength+Y
  inlet at maximum clearance. Its rim is.5MPa gauge, ends/cavity0gauge. These
  are chosen design inputs, not recovered Gangrade feed facts or a hole optimum.

The opening is excluded from the pressure domain. Incoming rim Couette flow
carries content1; outgoing flow uses the adjacent film state. Pressure exchange
is bidirectional and includes gross return, not a check valve or imposed flow.
Straight axial ends have no normal Couette velocity. Pressure-driven exchange
is still counted, including negative face flows. See the
[JFO formulation](../../../reynolds/docs/conical-jfo.md).

## Independent numerical checks

`transport_reference.jl` has no production includes. It contains:

- An expanded nodal theta-z finite-difference cone operator with independently
  evaluated radius/gap derivatives and scalar projected Gauss–Seidel. End-supplied
  10/20° cones refine64×32→128×64→256×128. The finite-difference end-gradient
  postprocessing is not exactly conservative at finite spacing; its signed
  imbalance is reported, not subjected to the FV roundoff conservation target.
- Continuous complementary pressure/content with curved interfaces in both
  coordinates. Artificial volume forcing is differentiated analytically and
  integrated with3×3Gauss quadrature on each actual chart quad. It is not
  `A*exact` and does not prescribe the active mask. The inlet lies in the
  fully-liquid part of this manufactured field. Known boundary pressure varies
  along the rim; this mathematical fixture is separate from the constant-pressure
  physical-equation study.
- A scalar projected pressure/void iteration for the fitted A/B matrices.
  This cross-checks nonlinear algebra, not a second fitted spatial method.

For the manufactured field, `L=.1m`, `t=cos(theta-pi/2)+.25cos(2pi*z/L)`,
`P=1e7*z*(L-z)/L^2*max(t,0)^4` and `fill=1-.2max(-t,0)^4`.
Its gap uses eccentricity.6 and offset angle.37rad, independently of the
physical-equation study's−Y offset.

Manufactured pressure/content errors use all cell centres with physical area
weights. Boundary-flow targets integrate the continuous flux along independent
oriented rim edges; force targets integrate the full curved cone minus the
opening polygon. Both rotation signs and graded tensor coordinates are checked.
The cylindrical limit reduces the cone angle without changing supply conditions.
Content is first-order upwind; measured pressure/content rates are reported.

## Outputs

`inputs.toml` fingerprints numerical source and inputs. Postprocessors separately
record their own code fingerprints. `analyse.py` uses explicit checks, including
when Python assertions are disabled.

- `physical.tsv`: peaks, physical cell counts, housing pressure-force components,
  gross in/out flows for each boundary, net feed, bounds/continuity diagnostics
  and self-changes. Housing sidewall force excludes end faces and shear; it is
  not a solved journal load or stability result.
- `manufactured.tsv`: continuous pressure/content, boundary-flow and force errors.
- `independent.tsv`: independent unported FD/FV field differences and FD end
  imbalance. `algebraic.tsv`: fitted same-matrix iteration controls.
- `cylindrical-limit.tsv`: pressure/content changes as angle approaches zero.
- `profiles-*.tsv`: reconstructed fields at actual25/50/75mm planes. Holes are
  NaN, not fabricated land cells. Fitted physical self-changes use common land
  samples on these planes, **not whole-field norms**. Tensor comparisons sample
  the fine solution at actual coarse cell centres. Linear reconstructions are
  unclipped and can overshoot near sharp content jumps; bounds diagnostics
  apply to the solved cell fields, not the reconstruction.
- `fitted-*.tsv/.jls/.vtk`: representative finest fields on actual geometry.
  VTK rim lines are annotations, with NaN pressure/content/gap values.
- `analysis.json` and refinement/manufactured/independent/profile/field PNGs.

Local/global balance and complementarity test discrete convergence. Continuous
manufactured fields and independent finite differences test numerical accuracy;
self-changes are not continuum error bounds. No independent published full-field
conical JFO experiment has been matched by this procedure.

## Artificial reservoir fixture

`conical_jfo_verification.jl` retains the separate rectangular pressure-reservoir
and original manufactured checks. Its cell-mask perimeter is a numerical source,
not a resolved circular inlet. Its tables/manifests are different from those
consumed by the transport postprocessors above. See the
[original fixture definition](../../../reynolds/docs/conical-jfo.md#unforced-numerical-reservoir-test).
