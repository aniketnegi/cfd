# Swift–Stieber pressure-only cavitation

`swift_stieber.jl` applies a nonnegative pressure obstacle to the shared rigid,
steady, isoviscous Reynolds assembly. It supports supplied cylindrical/conical
tensor meshes, including adjusted spacing, and fitted land meshes with an
explicit prescribed rim pressure. Geometry, normal-gap and
wall-speed approximations are those of the [full-film formulation](conical-formulation.md).
This closure uses the full-film face operator for the supplied mesh; it does
not assemble a separate pressure discretization.

## Equations and unknowns

Let `P=p-p_cav`. After imposing prescribed pressure cells, the free equations
satisfy

```text
P >= 0,  r=A*P-b >= 0,  P*r=0.
```

In positive-pressure cells the Reynolds equation holds (`r=0`). In the cavity,
pressure is fixed to `p_cav`; its edge is unknown. For a smooth rupture front,
the classical condition is zero pressure above cavitation and zero streamwise
pressure slope. The active-set iteration changes both the pressure and the
pressure-supporting region; it is not a one-pass clamp of a signed solution.

For the shared SPD tensor matrix this is a convex pressure obstacle. The
implementation reuses `solve_jfo_active_set` with a diagonal reaction block as
an algebraic kernel. That auxiliary reaction is NOT a transported void field.
The iteration cap reports failure; no unconverged field is silently returned.

The [fitted operator](fitted-full-film.md) can be nonsymmetric with positive
off-diagonals. Its complementarity equations are the same, but not the gradient
conditions of the tensor convex energy. The same algebraic active-set kernel
uses general sparse solves, followed by pressure/reaction and flux checks.
No global convergence or uniqueness guarantee covers arbitrary fitted meshes.
Manufactured fields and a separate projected-SOR algebraic comparison check
covered configurations independently of a small residual.

## Supplied-mesh API

From `journal-bearings/`:

```julia
include("bearing.jl")
include("mesh/2d/grid.jl")
include("simulation/reynolds/src/swift_stieber.jl")

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
grid = housing_grid(p,128,64)
result = solve_swift_stieber(grid; rpm=2000., viscosity_pa_s=.0277,
    end_pressure_pa=(0.,0.), cavitation_pressure_pa=0.)
```

These are explicit numerical-example conditions, not adopted Gangrade inlet
conditions. Lengths in mesh inputs are mm; computed pressures are Pa, viscosity
Pa*s, normal gaps metres, face flows m3/s and housing forces N.

For an opening, additionally load `mesh/2d/surface_feed.jl`, construct
`surface_feed_grid`, and supply `feed_pressure_pa` to `solve_swift_stieber`.
It acts on rim faces, not on a ring of fixed land cells. The pressure must be
finite and at least `cavitation_pressure_pa`; no supply pressure is selected
implicitly. Unported grids reject the feed-pressure keyword.

Keywords:

- `rpm`, `viscosity_pa_s`: prescribed speed and positive dynamic viscosity.
- `end_pressure_pa=(p_large,p_small)`: pressures at physical axial faces,
  ordered `z=0,z=L`; each must be at least `cavitation_pressure_pa`.
- `cavitation_pressure_pa=0.`: constant pressure floor in the SAME reference.
- `feed_pressure_pa=nothing`: required rim pressure on fitted meshes, in that
  reference. It may equal the pressure floor; this is not a liquid-content or
  replenishment specification.
- `eccentricity_angle_rad=0.`: offset orientation relative to the fixed mesh;
  result-mesh pose/gap handling is the same as the full-film API.
- `pinned_cells=Int[]`: optional cell indices fixed to `p_cav`; at least one
  free pressure cell must remain. A column approximates a prescribed line
  only when its pressure-centre alignment is explicitly chosen.
- `max_iterations=2000`: active-set cap, not an accuracy claim.

A pinned cell/line is not a circular inlet, supply-pressure patch or incoming
liquid-content condition. No start line or groove is selected implicitly.
Axial reservoirs and periodic angular boundaries are numerical conditions;
they do not reconstruct an experimental oil passage.

## Results and interpretation

The result contains `grid`, `pressure_pa`, `pressure_above_cavitation_pa`,
`normal_gap_m`, `reaction_m3_s`, the Boolean `pinned_cells` mask,
`face_flow_m3_s`, `boundary_outflow_m3_s`, `housing_pressure_force_n` and
`diagnostics`. Pressure order follows `grid.quads`; for tensors this is
`i+n_theta*(j-1)`. Tensor results also retain theta/axial flow arrays. Fitted
results use their supplied owner/neighbour faces, with positive outward flow
and `large`, `small`, `feed` boundary totals. The input housing is not moved.

Fluxes extrapolate the fully-liquid Couette–Poiseuille law through all cells.
Positive-pressure residuals and complementarity are checked separately from
flux reconstruction. Diagnostics report:

- pressure/flow normalization scales and iteration count;
- positive-pressure residual and complementarity;
- reconstruction error between face balance and `A*P-b`;
- fully-liquid cavity flux defect and pinned-line exchange;
- outward axial end flows and their defect-accounting residual;
- for fitted meshes, outward feed flow, included in global boundary accounting.

The sum of boundary outflows equals the cavity defect plus pinned-cell
exchange, up to algebraic roundoff. The axial end total alone need not equal
feed injection in a pressure-only cavity model.

Cavity reaction accounting is NOT JFO mass-conserving transport or a predicted
liquid inventory. End fluxes are pressure-only model outputs, not experimentally
validated oil consumption. The method does not establish a physical reformation
supply contract when no supply geometry has been specified.

The reported force acts on the housing sidewall, excluding shear and end faces.
It is not journal force, load equilibrium or stability. A common pressure
reference shift preserves `P` and flux but changes the open cone's axial
sidewall force; complete structural force balance also requires other surfaces
and external pressure.

For pressure/normal-gap VTK, separately include `pressure_vtk.jl` and call
`write_conical_pressure_vtk(path,result)`. The named `run.jl` CLI remains a
signed-full-film cylindrical runner; it does not select this closure.

## Checks and limits

```sh
julia --project=. simulation/tests/runtests.jl swift_stieber_checks.jl fitted_swift_stieber_checks.jl
```

Checks cover supplied/adjusted meshes, pressure references, pinning, scaling,
quiescent limits, full-film recovery and distinction from signed clipping.
Residual thresholds establish algebraic convergence, not continuum accuracy.
[Swift–Stieber verification](../../verification/swift-stieber/README.md)
provides independently formulated finite-width references, analytical end-flux
checks, conical manufactured fields and coordinate-aware refinement runners.

Thermal effects, elasticity and load equilibrium are not implemented here. [JFO](conical-jfo.md) is a separate closure with liquid-content
transport and an explicit supply/inventory requirement. [Gangrade's reported
inputs](../../../docs/gangrade-reference.md) do not fully define an inlet model.
