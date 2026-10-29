# Signed full film on fitted land meshes

`solve_full_film` accepts the supplied `surface_feed_grid` mesh. The circular
opening is excluded from the film domain; its rim carries prescribed pressure.
The solver does not model flow through the drilling or determine its mouth
pressure from an upstream reservoir.

```julia
include("bearing.jl")
include("mesh/2d/grid.jl")
include("mesh/2d/surface_feed.jl")
include("simulation/reynolds/src/full_film.jl")

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
grid = surface_feed_grid(p,128,64; feed_diameter=5., feed_z=50.)
result = solve_full_film(grid; rpm=2000., viscosity_pa_s=.0277,
    end_pressure_pa=(0.,0.), feed_pressure_pa=3e5)
```

These are numerical-example inputs. No diameter or boundary interpretation is
inferred from Gangrade. The feed pressure must be explicit and finite; unported
meshes reject it. All pressures use one reference. Signed negative values are
retained, not clipped or interpreted as physical cavity pressures.

## Surface flux

Use metre chart coordinates `x=R_ref*theta`, `y=z/cos(alpha)`. The same normal
gap, local nominal journal wall speed and cone geometry approximations as the
[tensor formulation](conical-formulation.md) give

```text
K = h^3/(12*mu)
M = diag(R_ref/r_h, r_h/R_ref)
j = (omega*r_j*h/2, 0) - K*M*grad(p)
```

A face has an owner, a neighbour (zero at a boundary), chart length `ell` and
outward chart unit normal `n`. Its centre connection is `d`. Define

```text
a = (n' M n)/(n' d)
q_pressure = ell*K*a*(p_owner-p_neighbour)
             - ell*K*(M*n-a*d)'*grad(p)_face
```

At boundaries the face pressure replaces the neighbour pressure. Weighted
least-squares cell gradients include those prescribed face values. Interior
face gradients average the two cell gradients; a boundary uses its owner
reconstruction. The correction accounts for nonorthogonal connections and the
cone metric. Simply using two-point distances on fitted cells omits it.
Mobility, metric and Couette flux use face-midpoint values.

The assembled shared face contract is `q=F*p+c`, `A=D*F`, `b=-D*c`. Each internal
face enters its owner with +1 and neighbour with -1 in `D`; boundary flow is
positive outward from the film. The signed sparse linear solve is the same
`A\b` entry point as for tensor grids. Its matrix need not be symmetric or
monotone; the tensor SPD/convex-obstacle argument does not apply here.

## Results

Pressure and normal gaps follow `grid.quads`, not an `n_theta*n_axial` array.
`face_flow_m3_s` follows `grid.surface.owner` and uses its outward-owner normals.
`boundary_outflow_m3_s` contains `large` (z=0), `small` (z=L) and `feed` totals.
Negative feed flow means net injection into the land. Global balance includes
all three boundaries. End-flow diagnostics contain the two axial ends only.

Fitted results have no rectangular `circumferential_face_flow_m3_s` or
`axial_face_flow_m3_s` arrays. Tensor results retain both arrays, their existing
sign conventions and diagnostics, and additionally expose the common face flow.
For tensor faces the ordering is flattened circumferential faces followed by
flattened axial faces; lower-end faces reverse the axial-array sign to point
outward from their owner.

`assemble_full_film` exposes the face operator as `faces`; use
`full_film_face_fluxes(assembly,pressure)` to reconstruct `q` and `D*q` for either
layout. `full_film_fluxes` also retains the tensor `qt`/`qz` view.

Housing force integrates piecewise-constant pressure over the curved cone
sidewall using 3x3 quadrature in each bilinear chart quad. It uses the supplied
radii and angle, excludes the opening, and does not include end-face pressure,
shear or passage forces. Housing vertices and connectivity stay fixed when
`eccentricity_angle_rad` changes; returned gap annotations follow that pose.
The [pressure exporter](../src/pressure_vtk.jl) works on both layouts. Its
pressure/gap scalars use NaN on rim-line annotation cells, followed by the land
values in quad order.

## Sampling an axial plane

Separately load `src/pressure_sampling.jl` and call
`fitted_pressure_slice(grid,pressure; z_m, boundary_pressure)`. The boundary
callback takes metre-chart `(x,y,tag)` and returns pressure in the same reference.
The function intersects actual land quads with the requested plane and evaluates
a least-squares linear reconstruction, including the prescribed rim trace.
It returns angular segments, their endpoint pressures and the plane maximum.
Shared horizontal edges use half-open ownership; no tensor row is assumed.
Reconstructed values are not clipped or additional solved unknowns. Check
refinement of the sampled maximum separately from the global native-cell peak.

## Checks and scope

`fitted_full_film_checks.jl` covers the common tensor face representation,
constant/affine pressure and exact face flux, rotating manufactured fields,
pressure/viscosity/speed scaling, inlet/end balance, Cartesian force quadrature
and pose reuse on fitted cylinders and cones. Manufactured forcing is derived
from the continuous equation, not from applying the discrete matrix to a target.

[Fitted Swift–Stieber](swift-stieber.md) applies the pressure-only complementarity
conditions to this same operator, without importing the tensor SPD guarantee.
[Fitted JFO](conical-jfo.md#fitted-face-transport) transports liquid content on
these faces with an explicit above-cavity rim pressure. Fitted load balance and
named simulation cases remain unsupported; the named CLI accepts unported
cylinders only.
