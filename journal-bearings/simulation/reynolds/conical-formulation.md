# Conical full-film formulation and numerical verification

`conical.jl` implements **rigid, steady, incompressible, Newtonian,
isoviscous, signed full-film Reynolds flow**. It does not implement cavitation,
a localized feed, elastic deformation or load equilibrium. Negative pressures
are retained only for numerical verification. These results must not be used
as physical load capacity, oil consumption or a Gangrade reproduction.

The returned force acts on the housing. A journal-load convention requires a
separate journal-force calculation.

The derivation below is our dimensional surface-flux formulation, not a claim
to have recovered Gangrade's Fortran or Fluent implementation. His paper's
pressure-boundary inconsistency is not used to define this numerical test.
See [the reference inputs and limits](../../docs/gangrade-reference.md).

## Coordinates, gap and approximation choices

Lengths in `BearingParameters` and the VTK geometry are millimetres. The
solver converts to **metres**, Pa, Pa·s, radians/s, m³/s and newtons.

- `z=0` is the large end; `z=L` is the small end.
- `theta=0` is +X; `theta=pi/2` is +Y (`theta=atan2(y,x)`).
- The journal is displaced by `e` along +X, without tilt or axial movement.
- Positive rpm means rotation about +Z in increasing theta.
- `alpha` is the semi-cone angle; `c` is radial clearance.

The existing geometry supplies nominal journal and housing radii:

```text
rj(z) = Rj_mid + (L/2 - z) tan(alpha)
r(z)  = rj(z) + c
```

The actual eccentric journal intersects a radial ray at

```text
rhoj(theta,z) = e cos(theta) + sqrt(rj(z)^2 - e^2 sin(theta)^2).
```

We evaluate the radial gap without subtracting nearly equal radii:

```text
g = c - e cos(theta)
      + e^2 sin(theta)^2 / [rj + sqrt(rj^2 - e^2 sin(theta)^2)]
h = g cos(alpha).
```

**`h` is a thin-film normal-gap approximation**, not an exact normal ray
intersection of the offset cones. It is exact for concentric parallel cones
and for the zero-angle cylinder. Its leading term is
`(c-e cos(theta)) cos(alpha)`, consistent with the radial-clearance convention
of Gangrade Eq. (8). Keeping the small radial-intersection correction retains
compatibility with `bearing.jl` and the existing cylindrical verification.
The tests independently intersect a housing-normal ray with the journal and
bound this approximation's relative discrepancy by 1e-7 at selected points
of the 10°, e/c=0.6 reference geometry. This is not a global error bound for
all clearances, angles or offsets. Near axial edges the comparison uses the
extended analytic cone, not a finite-end-face intersection.

The wall-speed approximation is `Utheta=omega*rj(z)`, `Us=0`: nominal local
journal radius, no meridional Couette transport. We neglect offset corrections
to tangential wall velocity, through-film curvature and higher-order thin-film
terms. This is a leading-order lubrication model, not exact moving-wall 3D
kinematics. In particular, retaining a geometric correction in `g` does not
make the entire model higher-order accurate. No validity claim for inertia,
thermal effects or the experimental operating regime follows from these tests.

## Surface metric and conservative equation

Meridional distance satisfies `ds=dz/cos(alpha)`. Surface distances and area are

```text
circumferential: r dtheta
meridional:     dz/cos(alpha)
dA:            r/cos(alpha) dtheta dz.
```

For local parallel-film Couette–Poiseuille flow, flux per unit transverse
width (m²/s) is

```text
K = h^3/(12 mu)
qtheta = Utheta h/2 - K/r * dp/dtheta
qs     = -K cos(alpha) * dp/dz.
```

Surface continuity is

```text
(1/r) d(qtheta)/dtheta + (1/r) d(r qs)/ds = 0.
```

Equivalently, in theta–z coordinates,

```text
d/dtheta [K/(r cos(alpha)) * dp/dtheta]
+ d/dz [K r cos(alpha) * dp/dz]
= d/dtheta [Utheta h/(2 cos(alpha))].
```

Each term of this integrated-coordinate equation has units m²/s; multiplying
by `dtheta dz` gives a cell volume balance in m³/s. The essential cone factors
are the varying `r(z)` and the **opposite** occurrences of `cos(alpha)` in
circumferential and axial conductances. Merely changing the radius in the
cylindrical solver would miss them.

At zero angle this reduces to the cylindrical formulation already documented
in [README.md](README.md), including its nominal journal speed and housing
surface metric. The common Couette–Poiseuille law is also the basis of
[San Andrés, Notes 4, Eq. 4.1](https://rotorlab.tamu.edu/me626/Notes_pdf/Notes04%20Static%20Load%20JBs%2010.pdf).
The cone metric, the analytical solution below and the manufactured source
are derived here; they are not independent published conical validation data.

## Discretization and boundary conditions

`conical.jl` uses a uniform periodic theta–z housing grid with cell-centred
pressure. Its assembly is separate from `solve_cylindrical_full_film`;
`cylindrical.jl` supplies geometry and VTK dependencies.

For each circumferential face, evaluated at the face centre:

```text
Gtheta = K * (dz/cos(alpha)) / (r dtheta)
Qtheta = Utheta h (dz/cos(alpha))/2 + Gtheta*(p_left-p_right).
```

For each meridional face:

```text
Gz = K * r dtheta / (dz/cos(alpha))
Qz = Gz*(p_large-side-p_small-side).
```

The same face flux belongs to both adjacent cells. Circumferential neighbours
wrap across the seam. End pressures are prescribed **at the physical faces**;
the distance to the first/last pressure cell is half a cell, doubling the
boundary conductance. End-pressure order is `(z=0, z=L)`, in Pa gauge.
The default is `(0,0)`. Nonzero ends prescribe end reservoirs; they do not
represent a localized feed port.

Assembly produces `A*p=b`, solved with Julia SparseArrays/SuiteSparse.
Diagnostics independently reconstruct cell continuity from returned face
flows. `Qz` is positive toward +z; outward end totals are `(-sum(Qz0),sum(QzL))`.
The L1 continuity ratio divides summed absolute cell imbalance by summed
absolute theta and z face flows; the max-cell ratio uses the mean of that
scale per cell. The global ratio divides absolute net outward flow by the
sum of absolute individual end-face flows. Ratios become uninformative near
zero flow; inspect absolute residuals for quiescent constant-pressure tests.
Signed inflow/outflow can cancel: the absolute-end-flow diagnostic is **not**
physical oil consumption in an uncavitated signed solution.

## Pressure-force convention

The normal from fluid into housing is

```text
n = cos(alpha)*(cos(theta),sin(theta),0) + sin(alpha)*(0,0,1).
```

We report `F_housing = integral(p*n*dA)` over the conical wetted surface.
For piecewise constant cell pressure, the oriented area is integrated exactly:

```text
Ax = r(z_mid) dz [sin(theta_hi)-sin(theta_lo)]
Ay = r(z_mid) dz [cos(theta_lo)-cos(theta_hi)]
Az = r(z_mid) dz dtheta tan(alpha).
```

This is **force on the housing**, not force on the journal. End faces and
viscous shear are excluded. A journal force calculation needs its own surface
projection and sign convention; do not silently reuse the housing vector as
an exact journal force or claim equilibrium with the paper's applied loads.
A uniform pressure gives zero transverse force and
`Fz = p*pi*(r_large^2-r_small^2)`. Tests check this and independent transverse
integrals for sinusoidal pressures. A constant gauge-pressure offset changes
this open-surface axial force; complete assembly force balance needs external
pressure and other loaded surfaces as well.

## Independent verification problems

### Centred, pressure-driven cone

With `e=0`, `h=c cos(alpha)` is constant and no rotation is necessary.
The equation reduces to `d/dz(r dp/dz)=0`. With `r=r0-kz`, `k=tan(alpha)`,

```text
I(z) = integral_0^z dz/r(z) = -log(1-kz/r0)/k
p(z) = p0 + (pL-p0)*I(z)/I(L)
Q(+z) = 2pi h^3 cos(alpha)/(12 mu) * (p0-pL)/I(L).
```

Use `I(z)=z/r0` when alpha=0 and `log1p` for numerical stability. This checks
varying-radius axial conductance, normal gap, end half-cells and flow signs.
The 0.5 MPa end pressure here is an **analytical test input**.

### Two-dimensional manufactured solution

For a centred cone take `p=P sin(theta) z(L-z)/L²`, zero end pressures and
zero rpm. Independently differentiate

```text
S = -d/dtheta[K/(r cos(alpha))*dp/dtheta]
    -d/dz[K r cos(alpha)*dp/dz].
```

Integrate `S` exactly in theta and by midpoint quadrature in z, and solve the
assembled operator with that artificial cell source. The source exists only
in the test, not in the physical solver interface. This tests both cone
metric terms together against a known two-dimensional pressure field.

### Rotating eccentric manufactured solution

For `p=P sin(theta) z(L-z)/L²`, e/c=0.6, alpha=10° and rpm=±2000,
the tests analytically differentiate the eccentric gap, diffusion coefficients
and Couette source independently of the discrete stencil. Relative pressure
errors decrease approximately fourfold per grid doubling; on 256×128 they
are 1.76e-4 and 1.29e-4 for positive and negative rotation.

### Other checks

- Zero-angle agreement with the cylinder and small-angle continuity.
- Stationary constant pressure and reversed end-pressure flow.
- Pressure-matrix symmetry, positive energy and maximum principle.
- Rotation reversal, viscosity scaling and pressure-offset invariance of flow.
- Housing-force integrals and pressure/force self-convergence.
- Invalid inputs, scalar field sizes and VTK normal-gap output.

## Commands and numerical results

From `journal-bearings/`:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/run_checks.jl
```

A case-specific script uses the existing checked mesh and geometry manifest:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/conical_verification.jl
```

It requires `out/cases/conical-unported-e0p6/experiment.toml` and matching
`mesh.vtk`, checks the fixed 50 mm / 100 mm / 50 µm / 10° / e/c=0.6 geometry
and 512×256 saved grid, and writes only a **new** run directory:

```text
out/cases/conical-unported-e0p6/runs/full-film-verification-vg32-40c/
  verification.toml
  convergence.tsv
  full-film.vtk
```

It refuses to overwrite that directory. To repeat, move the existing run first. `verification.toml` records inputs, source/mesh hashes,
completion status and diagnostics. It is deliberately not a general-run
`run.toml`; the existing `run.jl` CLI remains cylindrical-only. VTK coordinates
and radial gaps remain in mm; `normal_gap_m` is in metres and pressure in Pa.

The demonstration uses assumed 2000 rpm and nominal ISO VG 32 at 40 °C:
`nu=32e-6 m²/s`, separately assumed `rho=860 kg/m³`, hence `mu=0.02752 Pa·s`.
It does not establish the rig's oil temperature or density.

Numerical results:

- Centred 10° cone, 128 axial cells: relative L2 pressure error **4.78e-6**;
  maximum relative axial-flow error **1.33e-6**, against the analytical result.
  Both errors decrease by about 4× per refinement.
- Manufactured 2D field, 16×8 through 128×64: relative pressure errors
  **0.0138161, 0.00350343, 0.000879018, 0.000219953**.
- Eccentric cone, five grids below: area-weighted fine-to-coarse comparison
  gives observed second-order convergence, not an absolute error bound.

| Grid | max signed p (MPa) | housing Fy (N) | relative pressure mismatch from previous |
| --- | ---: | ---: | ---: |
| 32×16 | 17.400746 | -117878.22 | — |
| 64×32 | 17.436855 | -117272.79 | 0.0131067 |
| 128×64 | 17.414731 | -117120.10 | 0.00326819 |
| 256×128 | 17.414783 | -117081.85 | 0.000816569 |
| 512×256 | 17.413556 | -117072.28 | 0.000204113 |

At 512×256, relative L1 cell imbalance is **6.21e-15** and relative global
imbalance **1.98e-14**. The minimum pressure is approximately **-17.41 MPa
GAUGE**: the negative branch and the signed force above are deliberately
nonphysical verification outputs. The approximately zero axial force in this
antisymmetric signed solution does not show that a cavitating cone cannot
carry axial load. The prescribed offset is not the equilibrium for 1000 N.
These signed fields are not a comparison with Gangrade's measured pressure.

For cavitation with an explicit supply contract, see [conical JFO](conical-jfo.md).
Ambient ends alone do not guarantee replenishment when cavitation pressure is
also ambient and axial Couette transport is absent.
