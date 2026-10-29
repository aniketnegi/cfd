# Linear-elastic bearing solids

Static, small-strain, homogeneous isotropic elasticity in three Cartesian
coordinates, with SI lengths, forces and stresses:

```text
epsilon = (grad(u) + grad(u)')/2
sigma   = lambda*tr(epsilon)*I + 2*mu*epsilon
div(sigma) = 0
mu = E/(2*(1+nu)), lambda = E*nu/((1+nu)*(1-2*nu))
```

`src/linear_elasticity.jl` assembles isoparametric Hex8 elements with 2×2×2 Gauss
quadrature. The strain vector uses engineering shear components. Positive
Jacobians are required. Displacement constraints eliminate fixed DOFs; sparse
Cholesky solves the symmetric elastic stiffness. The factorization is reusable
for changing pressure loads in geometrically linear statics.

`src/housing.jl` constructs radial layers from the supplied fluid land quads.
The bore vertices match the fluid mesh. At a fitted +Y drilling, rim columns
retain their physical x,z coordinates through the wall; the opening is a
straight polygonal approximation of the drilling, not a radially flaring port.
A smooth surrounding map joins those columns to ordinary radial layers. Every
volume element is checked at all stiffness quadrature points.

## Selected housing model

The FSI example uses a constant 5 mm radial wall, homogeneous E = 97 GPa and ν = 0.30,
all outer-surface displacement components fixed, and traction-free axial ends.
These are explicit project inputs. The material/thickness values are informed
by [Liu et al. 2019, Table 3](https://www.mechanics-industry.org/articles/meca/full_html/2019/02/mi180284/T3.html),
not a reconstruction of its layered housing or the Gangrade rig.

Pressure acts on the lubricated bore; the drilling wall carries the prescribed
reservoir pressure, with no passage-loss or pipe-mechanics model. Loads are
pressure excess relative to the reference ambient state.

The bore load matrix integrates the Q1 virtual displacement against pressure
on the nominal curved surface. Partition of unity preserves its integrated
force. The transpose gives the work-conjugate mean normal displacement:

```text
f_land = W*p
mean_normal_displacement = W'*u / interface_area
u'*f_land = sum(p*interface_area*mean_normal_displacement)
```

The volume geometry is Hex8, whereas the load uses the nominal curved surface.
Geometric torque accuracy therefore needs refinement; exact moment preservation
against the curved continuum is not asserted. This distinction must not be
hidden by force or energy checks alone.

`solve_housing` returns displacement, reactions, energy and free-DOF residual.
`elastic_stress` evaluates Gauss-point stress and engineering strain. There is
no plasticity, yield criterion, thermal expansion or dynamics.
A small-strain numerical guard does not certify a particular alloy's strength.

## Solid journal

`src/journal.jl` constructs a solid cylinder/conical frustum with a square
core and surrounding O-grid quads, extruded along its axis. There is no
collapsed central Hex8 row. The journal uses its own angular, axial and radial
mesh counts, independently of the fluid and housing meshes.

The selected model spans the 100 mm bearing with no bore or overhangs,
E = 210 GPa and ν = 0.30; all displacement components of both end faces are
fixed relative to the prescribed nominal pose. These are user-adopted test
inputs. They do not reproduce a measured shaft/support arrangement.

The load operator transfers piecewise-constant FV land pressure after clipping
each land polygon at journal FE edges. The actual opening polygon is filled
by the specified reservoir-pressure trace. Pressure acts inward on the
journal's own nominal curved radius and area. It is not the negated housing
load. Surface interpolation is continuous Q1; volume geometry remains Hex8,
so curved-continuum moment accuracy still requires spatial verification.

`src/vtk.jl` exports the reference volume in metres, with separate displacement,
applied-load and reaction vectors. Use a displacement warp for visualization.

## Numerical checks

`simulation/tests/elasticity_checks.jl` checks partition of unity, rigid motions,
an affine displacement patch, energy, pressure transfer and a continuous 3D
manufactured solution with independently differentiated Navier forcing.
A restrained thick-cylinder solution checks radial displacement convergence.
That analytical control constrains all axial displacement (plane strain); it
is separate from the example's free-ended 3D housing.
`journal_checks.jl` adds polygonal-frustum volume, complete land/opening
pressure coverage and an independent solid-cylinder external-compression
reference. That reference constrains axial strain and only anchors remaining
rigid motions; it is not the example's two-end-clamped journal.
