# 2D periodic conical mass-conserving cavitation

## What is implemented

`conical_jfo.jl` solves steady rigid, isothermal, incompressible Newtonian
pressure **and liquid content** around the circumference and along the cone.
It uses the [conical surface formulation](conical-formulation.md);
it is not a stack of independent slider calculations.

**Numerical verification only.** The saved example uses an intentionally broad
artificial reservoir patch, not a recovered rig inlet or the project's 5 mm
opening. No physical VG32 cavitation pressure, experimental supply condition,
load equilibrium, thermal model or elasticity has been selected by this test.
No new published Gangrade evidence is claimed.

The [floating-bush driver](conical-equilibrium.md) uses this solver for rigid
force balance. `eccentricity_angle_rad` rotates the gap relative to the fixed
housing mesh/feed; its default is zero (+X).

`jfo_active_set.jl` supplies the pressure/void active-set kernel shared with
the [1D slider](jfo-slider.md). The two solvers have separate flux assemblies.

## Equations and supply contract

Coordinates, cone radii, normal-gap approximation and nominal wall speed follow [conical-formulation.md](conical-formulation.md):
`z=0` large end, positive theta about +Z, `ds=dz/cos(alpha)`,
`Utheta=omega*rj(z)`, `Us=0`. Geometry/VTK lengths are mm; the solver uses SI.

Let `P=p-p_cav`, `phi=liquid content`, `v=1-phi` and `K=h³/(12 mu)`:

```text
P >= 0,  0 <= phi <= 1,  P(1-phi) = 0.

Qtheta = gt (P_left-P_right) + adv phi_upwind,
gt = K ds/(r dtheta),  adv = Utheta h ds/2,

Qz = gz (P_lower_z-P_higher_z),
gz = K r dtheta/ds.
```

These are integrated **liquid-equivalent** face volume flows in m³/s; multiply
by constant liquid density for mass flow. Axial-end conductance uses the
half-cell distance. The theta faces are periodic, including the last-to-first
cell face, and upstream content changes with rotation sign. There is no axial
Couette flow, so no additional axial advective content boundary is prescribed.

The checked full-film assembler supplies `A`, pressure conductances and fully
liquid transport. A separately assembled conservative upwind matrix `B` gives
`A P + B v = b`. Each internal face is shared with opposite divergence signs.
Pressure mobility is not multiplied again by content: pressure is constant
inside a cavity, and content is one inside full film.

**Supply is mandatory for this steady public solver:**

- Axial reservoir pressure may be at or above p_cav. At least one end or
  an explicit internal reservoir must be strictly above p_cav.
- With no feed and both ends at p_cav, the solver **rejects the case** rather
  than inventing a circulating underfilled inventory or calling it flooded.
  A separate inventory/transient formulation would be required to select such
  a state. This is a deliberate scope restriction, not proof that no solution
  exists.
- A supplied Boolean `feed_cells` vector (theta-fast ordering) and
  `feed_pressure_pa > p_cav` prescribe fully liquid reservoir cells. Empty or
  all-domain masks, missing pressure/mask pairs and inconsistent inputs fail.
- Feed pressure is imposed at **cell centres**, with the ordinary cell-to-cell
  pressure conductance. This is not a fitted hole boundary or a half-cell
  pressure condition on the mask perimeter. Cell size affects that approximation.
- Feed is **bidirectional**: the pressure reservoir may inject or accept oil.
  There is no check valve, restrictor, pipe resistance or fixed flow rate. Both
  gross injection and gross return are reported; their difference is net feed.
- A feed mask is passed explicitly; the mesh's geometric opening marker does
  not silently become a pressure boundary. All pressure inputs use one caller-
  specified reference, and p_cav is required explicitly.

Fixed feed pressures/zero void are eliminated before the shared active-set
solve. No pressure or content clipping follows it. Nonconvergence, singular
linear systems or failed bounds/continuity checks cause failure; global
convergence for arbitrary cases is not asserted.

### Feed-flow accounting

After reconstructing face flows independently of matrix multiplication:

```text
balance_cell = Qtheta_east-Qtheta_west + Qz_high-Qz_low,
source_cell  = balance_cell on reservoir cells; zero elsewhere,
residual     = balance_cell-source_cell,
net_outflow  = -sum(Qz_large_end)+sum(Qz_small_end) = sum(source_cell).
```

A positive source is injection; negative is return to the reservoir. The
feed-cell residual is not called a conservation error: it is the flow the
prescribed-pressure constraint must supply. Tests independently sum the
reservoir/film **perimeter flux**, including a reservoir crossing the periodic
seam, and compare it with this recovered source.

Local balance is normalized by the largest nominal Couette, actual face or
reservoir-pressure conductance flow scale. Global balance uses gross feed/end
exchange, with that local scale as a zero-flow fallback. Values and normalizers
are saved, so roundoff balances cannot be mistaken for discretization errors.

## Verification, not a published 2D benchmark match

Checks cover:

- Arbitrary pressure/content fields: matrix residual versus independently
  reconstructed face divergence; periodic seam and row-wise transport conservation.
- Full-film reduction at 0°, 10° and 30°, with stationary and rotating walls;
  agreement with the signed solver when its pressure stays positive.
- Independent centred stationary cone log-radius pressure solution.
- A genuinely 2D manufactured cavitating field, for **both** rotation signs.
- Pressure-reference shifts, speed/viscosity scaling, reflection under rotation
  reversal, and the near-zero-angle cylindrical limit.
- Reservoir pressures/full content, gross and net feed, outward end flow,
  seam-crossing reservoirs with return flow, and both end-feed orientations.
- Unforced patch-fed refinement, invalid/underspecified supply rejection,
  nonconvergence, and VTK overwrite protection.

### Independent manufactured solution

`conical_jfo_mms.jl` analytically differentiates gap, pressure and content:

```text
t(theta,z) = cos(theta) + 0.25 cos(2 pi z/L),
f(z)       = z(L-z)/L²,
P          = 1e7 f(z) max(t,0)^4,
phi        = 1 - 0.2 max(-t,0)^4.
```

The cavity boundary `t=0` curves in both coordinates. Exact source density is
computed from the continuous conical divergence of Couette/content transport
and pressure diffusion, then midpoint-integrated over each cell. It is not
obtained by applying the discrete matrix to the exact field. This exercises
the variable cone metric, eccentric gap, rotating content transport and axial
pressure diffusion together.

Artificial distributed injection/removal exists **only in this low-level
verification calculation**, not as a physical supply option in the public
solver. The manufactured interface is smooth; it does not independently test
a discontinuous physical 2D reformation front. The published 1D benchmark still
checks rupture/reformation and a content jump through the shared kernel.

| Grid theta × z | Pressure relative L2 error | Content mean absolute error |
| --- | ---: | ---: |
| 32 × 16 | 3.17270e-2 | 7.52404e-3 |
| 64 × 32 | 7.18314e-3 | 3.72843e-3 |
| 128 × 64 | 1.73931e-3 | 1.85979e-3 |
| 256 × 128 | 4.30024e-4 | 9.29272e-4 |

The two rotation signs give the same errors to roundoff by reflection.
Pressure is approximately second order; upwind content is first order.
The finest pressure error is **0.043%**, not a physical-bearing accuracy claim.
No independent published full-field 2D conical cavitation solution or experiment
has yet been matched.

### Unforced numerical reservoir test

This separately tests the public solver without manufactured forcing:

- 50 mm mean journal radius, 100 mm axial length, 50 µm radial clearance,
  10° semi-cone, prescribed e/c=0.6 along +X.
- 2000 rpm, mu=0.02752 Pa·s, an illustrative nominal-property reference,
  not a measured operating-temperature oil.
- Numerical p_cav=0 and both axial ends at zero; pressure reservoir at
  0.5 MPa over **45° < theta < 135°, 37.5 < z < 62.5 mm**.
- This intentionally coarse-friendly rectangle is **not** the 5 mm port,
  not Gangrade's documented inlet, and not an adopted physical boundary model.

| Grid | Pressure mismatch to previous grid | Content area-mean mismatch | Net reservoir flow (m³/s) |
| --- | ---: | ---: | ---: |
| 32 × 16 | — | — | 1.177904e-5 |
| 64 × 32 | 2.23092e-2 | 1.56690e-2 | 1.290400e-5 |
| 128 × 64 | 8.83940e-3 | 7.72337e-3 | 1.346784e-5 |
| 256 × 128 | 4.21051e-3 | 3.92426e-3 | 1.375932e-5 |

Fine fields are area-averaged onto coarse cells before comparison. This is
**self-convergence**, not error against a continuum oracle. It is slower than
the smooth manufactured pressure convergence: content fronts and prescribed
pressure cells affect it. Net feed still changes **2.16%** on the last
refinement; oil consumption is not grid-independent.

At 256×128: net end outflow and net feed agree with normalized imbalance
4.16e-14; maximum free-cell residual is 1.01e-20 m³/s. There are 11,350 cells
with content below 1−1e-8; minimum content is 0.255606. The numerical peak is
18.9693 MPa, **not an experimental load-equilibrium prediction**. Fifteen
active-set iterations were required.

The plot shows pressure/content for this artificial patch case and, separately,
manufactured-solution errors. The third panel is not an error bound for the
first two panels. The JFO result does not automatically report a journal force;
loaded-body/pressure-reference conventions and the role of a physical opening
still need review before equilibrium.

## Files and reproduction

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_jfo_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/jfo_slider_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/conical_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/conical_jfo_verification.jl
python expts/conical-jfo-verification/plot.py
```

Saved under
`out/cases/conical-jfo-numerical/runs/transport-verification/`:
`verification.toml` (inputs, normalizers and source hashes), `refinement.tsv`,
`manufactured.tsv`, `fields.tsv`, `cavitation.vtk`, `verification.png`, and
four check/regression logs. Julia 1.13.2. Existing runs are refused; pass a
fresh output directory to the verification script and then the plot script to
repeat. NumPy/Matplotlib are needed only for the optional plot.

## Physical scope

The artificial reservoir does not specify a rig port, experimental axial
restraint or lubricant cavitation pressure. These must be explicit physical
inputs for a bearing prediction. See [Gangrade reference limits](../../docs/gangrade-reference.md).
