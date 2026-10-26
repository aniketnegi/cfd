# Floating-bush JFO load-equilibrium checkpoint

## Status

For the complete operating-point interpretation, new feed/rotation controls,
measured CPU/RAM costs and proposed next steps, see the
[operating-point report](../../docs/conical-load-operating-point-report.md).
A subsequent [fitted-inlet diagnostic](../../docs/conical-fitted-inlet-comparison.md)
compares the same held pose: forces/flow settle faster, but transition quality
worsens, so that candidate is not adopted and no equilibrium is re-solved.
A later [controlled-transition trial](../../docs/conical-controlled-transition.md)
improves the surface quality substantially, but still gives 322.87 N axially
at the held pose. The user has since approved promoting it to the
[public surface-inlet generator](../../mesh/2d/conical/README.md), replacing that
generator's stepped mask. The equilibrium driver below still uses the regular
grid; fitted field integration must be explicit. A subsequent case-specific
[fitted load continuation](../../docs/conical-fitted-load-balance.md) now obtains
18 radial equilibria/checks: finest axial forces are 319.35–331.07 N over the
sampled 12.5–100 µm clearance range, with no 200 N crossing found. This is not
a simultaneous free-axial equilibrium, global nonexistence proof or 3D Hex8
acceptance. The historical regular-grid checkpoint below is retained.

**Radial equilibrium is implemented and checked. Simultaneous radial/axial
equilibrium for the selected case has not been established.** Do not use a
failed three-coordinate iterate as an operating point or start elasticity from it.

The user approved these explicit hypotheses, not a recovered Fig.8 fixture:

- Rotating shaft axis located; nonrotating, untilted bush with three translations.
- Net external force on the bush `(0,+1000,-200)` N. +Y is the feed side;
  +Z points from the large towards the small cone end.
- Longer journal retaining full overlap with the 100 mm bush. Report the cone
  extension required by translation; no real stop/travel capacity is asserted.
- Reference: coincident axial midplanes, journal mean radius 50 mm, radial
  clearance 50 µm, semi-cone angle 10°, speed 2000 rpm.
- Housing-fixed 5 mm +Y feed at mid-length, ideal fully liquid bidirectional
  reservoir at 0.5 MPa gauge **at the film mouth**. No pipe/restrictor/check valve.
- JFO, ambient axial exits and provisional `p_cav=0` gauge. Lubricant-specific
  characterization is explicitly deferred in `../.tasks/PROJECT.org`.
- Constant `mu=0.0277 Pa s` is the paper reference, not a measured operating oil
  state. No temperature is inferred. Report volume flow; density is unnecessary
  for this incompressible pressure/volume-flow calculation.

## Geometry and force contract

Let the bush translation relative to the located shaft be `b=(bx,by,bz)`, and
let z be measured from the large end in housing coordinates. Then

```text
journal centre relative to housing = (-bx,-by),
rj(z;bz) = Rj_ref + (L/2-z-bz)*tan(alpha),
c(bz)    = c_ref + bz*tan(alpha),
rh(z)    = Rj_ref+c_ref+(L/2-z)*tan(alpha).
```

Thus the bush geometry and feed remain fixed in housing coordinates, while
axial movement changes the journal radius sampled beneath it. The housing mesh
does not rotate with the eccentricity vector. The shared cone assemblers now
accept `eccentricity_angle_rad` (default zero, preserving earlier +X cases).
Only the gap is evaluated at the shifted angle. Grid radial-gap metadata follows
that shift; housing coordinates and the feed mask do not.

The projected normal-gap model is unchanged. Contact/intersection and invalid
cone radii are rejected. Axial displacement is a **steady position**, not an
axial sliding velocity: no axial Couette/squeeze term is invented.

Pressure force is integrated on the **solid bush lands**, excluding feed-mask
cells: there is no solid wall spanning the ideal opening. The nominal full-
surface resultant and omitted aperture traction are also reported. At 512×256
the omitted nominal traction is about `(0,9.594,1.692)` N. The mask is still a
stair-step approximation, not a fitted hole boundary.

The specified load is the **net of all non-modelled forces** on this body,
including any weight, passage/pipe and fixture contributions. These individual
contributions have not been recovered or calculated. This is not a claim that
the paper's actuator readout equals that net. Bush shear traction, tilt and
elasticity are absent; orientation restraints may supply moments. This is a
leading-order pressure-load model, not a full rig free-body reconstruction.

Equilibrium requires `F_land(b)+F_external=0`. At fixed axial seating only the
two transverse equations are solved; the missing axial support reaction is
explicitly `-(F_land_z-200)` N. Without that reaction the state is not a
three-force equilibrium.

For nonnegative gauge pressure on the cone lands, geometry supplies the
necessary bound `Fz >= tan(alpha)*hypot(Fx,Fy)`. At 1000 N radial load its lower
bound is 176.3 N. Therefore **200 N is not ruled out just by the cone angle**;
the pressure distribution and supply determine whether it can be attained.

## Numerical implementation and verification

`conical_equilibrium.jl` reuses JFO and the checked housing force integral.
A scaled, damped finite-difference Newton iteration solves two or three force
equations, rejecting inadmissible poses and recording inner-solver failures.
It reports convergence/failure, not equilibrium existence or dynamical stability.
The three-equation solve can be poorly conditioned when axial position has
little influence on the radially balanced force distribution.

**87 new assertions pass**, covering:

- Independent world-coordinate circle/ray geometry for translated cones.
- Unchanged housing geometry, axial gap changes and required end extensions.
- Rotational covariance of pressure, liquid content and force, and fixed-grid
  gap metadata with arbitrary eccentricity orientation.
- Analytical uniform-pressure axial force; aperture traction exclusion.
- Radial equilibrium, explicit support reaction, invalid inputs and failure status.
- Recovery of a known synthetic three-coordinate pose from its calculated load.
  This last test checks root finding only; it is not independent physics validation.

All 478 prior cone-JFO, slider, signed-cone and cylinder assertions also pass:
**565 assertions across five suites**. Earlier zero-angle production source
versions were preserved in `out/reviews/pre-equilibrium-sources/`; historical
comparison outputs were not modified. Their old source hashes naturally predate
the new optional orientation argument.

## Results obtained

Fixed nominal 50 µm radial clearance, after transverse equilibrium:

| Grid | Eccentricity ratio | Axial film force | Additional axial support needed |
| --- | ---: | ---: | ---: |
| 128×64 | 0.019084 | 278.80 N | −78.80 N |
| 256×128 | 0.018665 | 305.27 N | −105.27 N |
| 512×256 | 0.018650 | 313.54 N | −113.54 N |

At 512×256, bush transverse displacement is approximately
`(-0.91525,+0.17844)` µm and the analytic minimum normal gap is **48.322 µm**.
The transverse force residual is below 0.0001 N. Net feed is approximately
**0.03014 L/min**, with no reservoir return in this state. The maximum pressure
is 0.5 MPa, imposed by the feed; it is not an independently predicted peak.
Normalized global flow imbalance is approximately 2.2e-13.

**Not grid independent:** axial force and feed still change by about 2.7–2.8%
on the final refinement. The small pressure patch's finite-cell representation
and upwind transport need further refinement before reporting converged flow or
axial force. These changes are substantially smaller than the observed ~114 N
axial mismatch at fixed reference seating, but are not an error bound.

The first three-coordinate attempt at 256×128 exceeded a 600-second command
limit. That partial run is preserved, not called a solution. A fresh bounded
run caps each JFO active-set solve at 100 iterations and attempts joint balance
at 128×64 for 12 outer iterations, using two finite-difference steps (1e-5 and
1e-6 in scaled coordinates). **Neither closes force balance.** Their final
radial residuals remain roughly 86–139 N and axial residuals roughly 64 N.
Changing the finite-difference step did not produce a satisfactory solution;
this is not a proof of nonexistence.

To interpret the failure, separate radial equilibria were obtained at several
axial positions on 128×64:

| Operating radial clearance | Axial film force after radial balance |
| ---: | ---: |
| 12.5 µm | 275.24 N |
| 25 µm | 275.98 N |
| 50 µm | 278.80 N |
| 100 µm | 288.65 N |
| 200 µm | 321.56 N |

All exceed the target 200 N. This is **a sampled coarse-grid radial-equilibrium
branch**, not a global search or a guarantee of uniqueness. The different gaps
are consequences of imposed diagnostic axial positions, not viscosity/geometry
fitting to the paper. Larger-gap states also require reviewing the lubrication
model's validity before treating them as physical predictions.

The result suggests that the assumed mouth pressure, applied net loads and
absence of axial support may be incompatible on this branch. It does not yet
identify which physical assumption is wrong. Do not silently restore a support,
reduce feed pressure or declare a joint equilibrium from the failed iterates.

## Evidence and next step

Case and commands:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_equilibrium_checks.jl
julia --project=. --check-bounds=yes expts/conical-load-equilibrium/run.jl FRESH_DIRECTORY
python expts/conical-load-equilibrium/plot.py FRESH_DIRECTORY
```

Inputs: `expts/conical-load-equilibrium/case.toml`.
Completed bounded evidence:
`out/cases/conical-floating-bush/runs/jfo-load-checkpoint-bounded/`:
`states.tsv`, `iterations.tsv`, fixed-seating VTK fields, `run.toml`, run/test
logs, trial-failure logs and `load-balance.png`.

In `states.tsv`, **converged refers only to the equations selected for that
row**: `fixed-*` and `radial-branch-*` rows enforce transverse balance only.
The `free-axial-*` rows attempt all three equations and are not converged.

The bounded run recomputed fixed-seat fields from initial guesses preserved in
the interrupted run's iteration table; its hash is recorded. Plotting needed a
TSV-delimiter correction after the numerical run; the original plotting source
is preserved in `plot-before-delimiter-fix.py`, and `plot.json` records the actual
plot source/input/output hashes separately. Numerical evidence was not rewritten.

Next strengthen the axial-equilibrium search with safeguarded continuation of
radial equilibria and refine the feed treatment. Establish whether a target
crossing exists on a numerically resolved admissible branch before changing
physical hypotheses. No Fig.8 experimental validation, stable operating-point
claim or elasticity coupling has been completed.
