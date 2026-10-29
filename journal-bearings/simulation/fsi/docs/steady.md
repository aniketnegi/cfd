# Fixed-position housing FSI

One-way and two-way coupling between steady incompressible JFO and 3D
[elastic solids](../../elasticity/docs/linear-elasticity.md). Use an elastic
housing with either a rigid journal or an independently meshed elastic solid
journal. Nominal eccentricity and attitude are prescribed. Supports carry
reactions; no external-load equilibrium is imposed. This is not transient or
ALE coupling.

## One-way

Solve the rigid-film JFO problem, transfer its pressure to the housing, and
solve elasticity. The resulting deformed gap is checked, but pressure is **not**
recomputed. Its pressure remains the rigid result. Reported predicted gap and
undeformed fluid gap must not be confused.

## Two-way

On the fixed nominal housing coordinates:

```text
h(theta,z) = h_rigid(theta,z)
             + (u_housing(theta,z) - u_journal(theta,z)) dot n_housing(theta,z)
n_housing  = (cos(alpha)*cos(theta), cos(alpha)*sin(theta), sin(alpha))
```

Set `u_journal = 0` for the rigid-journal model. Inward journal displacement
opens the film; it must not be added with the housing sign.
The axial displacement term matters for a cone. The full vector is projected;
neither displacement magnitude nor radial displacement alone is the correction.
For a cylinder, purely axial displacement makes no normal-gap contribution.
The stiffness, reference surface metric and nominal wall speed remain fixed;
pressure mobility h³ and Couette transport h are recomputed from the new gap.
This is a small-displacement normal-gap EHL approximation, not a new deformed
solid/fluid geometry at every iteration.

The coupled fixed point solves fluid pressure at the current displacement,
then computes the elastic target under that pressure. Underrelaxation updates
displacement. Convergence requires:

- the normal-gap fixed-point defect below its absolute tolerance;
- relative change in pressure **above cavity**, independent of pressure datum;
- the actual free-DOF elastic residual of **each** body under the current
  fluid traction.

Returned fluid and displacement use the **same** gap iterate. A last elastic
target is not silently paired with pressure computed on another gap. Invalid
gaps or failure to converge cause an error; no contact replacement or gap
clipping is performed. A returned conservative lower bound subtracts the
largest interface-node displacement magnitude of each body from the minimum
rigid gap. A
positive bound certifies positivity throughout the Q1 interface, not only at
FV sampling points; a nonpositive bound alone does not prove contact.
Iterations represent equilibrium updates, not physical
time steps or a transported transient liquid inventory.

## Interface

The inner structural quads match the public FV land topology. Q1 displacement
is continuously reconstructed on those quads at the retained fluid cell/face
points. Coordinates are inverted on their known owner; no nearest-row lookup,
fake tensor reshape or extrapolation across the opening is used.

The structural pressure operator preserves integrated force and pressure work;
see the elastic-model page for the separate torque/geometric limitation.
The prescribed supply also loads the straight drilling wall. All pressures
must share a datum; `external_pressure_pa` defines the ambient reference for
solid tractions.

The solid journal has a square-core O-grid with no collapsed axis elements.
Its surface mesh can be coarser than the fluid mesh. Pressure transfer clips
FV land polygons at structural element boundaries before triangle quadrature;
there is no nearest-pressure lookup. Beneath the actual housing opening,
journal traction uses the uniform prescribed reservoir pressure. This closes
the ideal-reservoir model; it is not a recovered experimental passage field.
Journal loads use its own nominal radius/area and inward pressure normal,
not simply the negated housing force.

Both volumes are nominal reference meshes: rigid eccentric pose remains in
`h_rigid`, not in a geometric distance between faceted FE surfaces. Mesh chord
error must not be mistaken for film thickness. Rotational inertia, centrifugal
prestress, material-frame cycling/fatigue and thermal effects are omitted.
Structural tractions are pressure-only; viscous shear and shaft driving torque
are not transferred to the solids.

## API

Load `bearing.jl`, `mesh/2d/grid.jl`, and the fitted mesher when applicable,
then `simulation/fsi/src/steady.jl`.

```julia
model = prepare_housing_fsi(grid; youngs_modulus_pa=97e9, poisson_ratio=.3,
    wall_thickness_m=.005, n_radial=2, feed_diameter_m=.007, feed_z_m=.05)
fluid = (;rpm=2000., viscosity_pa_s=.0277, end_pressure_pa=(0.,0.),
    cavitation_pressure_pa=0., feed_pressure_pa=5e5,
    eccentricity_angle_rad=3pi/2)
one = solve_one_way_fsi(model; fluid_options=fluid)
two = solve_two_way_fsi(model; fluid_options=fluid, relaxation=.5)

both = add_elastic_journal(model; youngs_modulus_pa=210e9, poisson_ratio=.3,
    n_theta=32, n_axial=16, n_radial=4)
coupled = solve_two_way_fsi(both; fluid_options=fluid)
# coupled.solids.housing / coupled.solids.journal and coupled.stresses
```

For unported numerical end-supply cases, omit the drilling parameters and rim
pressure. Reservoir-cell masks are not a structural drilling model and are
rejected by this coupling. User-supplied `gap_correction_m` is also rejected:
FSI owns that feedback. The lower-level shared fluid assemblers accept such a
callback for mathematical tests; it must return a normal gap correction in
metres at absolute angle and axial position in metres.

The selected elastic journal is a homogeneous 100 mm solid cylinder/frustum,
E = 210 GPa, ν = 0.30, with both end faces fixed in all displacement components
relative to the prescribed nominal pose. It has no bore or overhangs. These are
explicit working assumptions, not recovered rig dimensions or supports.
The housing retains its separate 97 GPa, 5 mm wall, outer-fixed/end-free model.

Returned `solids` and `stresses` name every active body. For compatibility with
housing-only consumers, `solid` and `stress` are aliases for the housing, not
aggregates. The combined gap correction and convergence use all active bodies.
Each state's `normal_displacement_cell_mean_m` is pressure-work conjugate: it
is outward-positive for the housing but inward-positive for the journal.
The TSV journal outward-normal columns instead use the common geometric normal
and are negative under compression.
