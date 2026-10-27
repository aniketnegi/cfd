# Floating-bush JFO force balance

`conical_equilibrium.jl` solves two or three steady pressure-force equations
for an untilted, nonrotating bush translating relative to a located shaft.
It uses the uniform-grid JFO solver with an explicit feed-cell mask, not the
[fitted surface mesh](../../mesh/2d/README.md).

Radial equilibrium is tested. Simultaneous 1000 N radial / 200 N axial balance
is not established for the configured case. A failed iterate is not an
operating point; steady force balance is not dynamic stability.

## Configured physical assumptions

`expts/conical-load-equilibrium/case.toml` specifies:

- Journal reference midpoint radius 50 mm, bush length 100 mm, radial clearance
  50 µm at coincident axial midplanes, semi-cone angle 10° and speed 2000 rpm.
- Net external bush force `(0,+1000,-200)` N. +Y is the feed side; +Z points
  from the large toward the small end.
- A sufficiently long journal to retain full overlap during translation.
- Housing-fixed 5 mm +Y midpoint opening, represented by cell centres inside
  its footprint; ideal bidirectional full-liquid reservoir at 0.5 MPa gauge.
- Ambient axial exits and zero-gauge cavitation pressure.
- Constant dynamic viscosity 0.0277 Pa·s, not a measured operating oil state.

The feed diameter, net-load interpretation, full overlap and axial freedom are
project hypotheses, not a recovered experimental fixture. There is no passage,
restrictor, check valve, shear load, tilt, thermal model or elastic deformation.
Density is unnecessary for incompressible pressure and volume flow.

## Geometry and force contract

Displacement/load API inputs are SI. For bush translation `b=(bx,by,bz)` and
housing-local z measured from the large end:

```text
journal centre relative to housing = (-bx,-by)
rj(z;bz) = Rj_ref + (L/2-z-bz)*tan(alpha)
c(bz)    = c_ref + bz*tan(alpha)
rh(z)    = Rj_ref+c_ref+(L/2-z)*tan(alpha)
```

The housing and feed stay fixed. `eccentricity_angle_rad` rotates only the gap
and its metadata, not the mesh. The projected normal-gap approximation follows
[the conical formulation](conical-formulation.md). Invalid radii or contact are
rejected. Axial displacement is a steady position, not an axial sliding speed;
there is no axial Couette or squeeze term.

Pressure force acts on **solid bush lands**, excluding feed-mask cells.
Full-surface force and omitted aperture traction are also reported. The mask
is a staircase approximation, not a fitted hole boundary. External force is
net of all unmodelled weight, passage/pipe and fixture loads; those components
are not individually calculated.

Equilibrium requires `F_land(b)+F_external=0`. With `fixed_axial=true`, only
transverse balance is enforced at the initial axial seating. For the configured
200 N axial load, the required support reaction is `-(F_land_z-200)` N.
Without that reaction, radial balance alone is not full equilibrium.

Nonnegative gauge pressure on cone lands implies
`Fz >= tan(alpha)*hypot(Fx,Fy)`. At 10° and 1000 N radial force the lower bound
is 176.3 N: geometry permits, but does not guarantee, a 200 N axial solution.

## Solver and checks

A scaled, damped finite-difference Newton method solves the selected force
components. A small SVD guards nearly singular force-position maps. Invalid
poses and inner-solver failures are recorded; iteration exhaustion returns
failure, not a nonexistence claim.

Tests cover independently translated circle/ray geometry, fixed housing,
axial clearance/overlap changes, rotational covariance, uniform-pressure force,
aperture exclusion, radial balance, reaction forces and failure status.
Recovering a synthetic pose from its calculated force tests root finding,
not independent physical validity.

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_equilibrium_checks.jl
julia --project=. --check-bounds=yes expts/conical-load-equilibrium/run.jl FRESH_DIRECTORY
python expts/conical-load-equilibrium/plot.py FRESH_DIRECTORY
```

Run from `journal-bearings/`. The runner refuses an existing output directory
and writes `states.tsv`, `iterations.tsv`, `run.toml`, field VTKs and trial
failures. In `states.tsv`, `converged` applies only to that row's selected
equations: `fixed-*` and `radial-branch-*` enforce transverse balance only.

## Numerical limits of the configured case

At fixed 50 µm clearance, uniform-grid radial solutions give:

| Grid | e/c | Axial force | Additional axial support required |
|---|---:|---:|---:|
| 128×64 | 0.019084 | 278.80 N | −78.80 N |
| 256×128 | 0.018665 | 305.27 N | −105.27 N |
| 512×256 | 0.018650 | 313.54 N | −113.54 N |

At 512×256, the transverse offset is about 0.93248 µm, minimum normal gap
48.322 µm and net feed 0.03014 L/min. Axial force and flow change by about
2.7–2.8% on the final refinement; they are not grid independent. The 0.5 MPa
peak is prescribed by the reservoir, not independently predicted.

No matched Fig. 8 validation or stable free-axial operating point is established.
The fitted-load experiment is incompatible with the public mesher; see
[solver support](../../docs/solver-support.md).
