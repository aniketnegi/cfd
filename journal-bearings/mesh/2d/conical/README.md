# Conical surface meshes

## Files

- `grid.jl`: shared `background_grid` connectivity, `patch_axis` spacing, and
  `housing_grid` for the complete housing surface.
- `surface_feed.jl`: `surface_feed_grid`, quadrilateral cells around an inlet opening.
- `inlet/seed.jl`: initial background, opening rings and surrounding-node smoothing.
- `inlet/placement.jl`: local node optimization on the requested mesh.
- `inlet/metrics.jl`: cell/face connectivity, geometric measurements and quality.
- `inlet/controlled.jl`: runs that inlet construction and checks its quality limits.
- `../../vtk.jl`: `write_vtk`, shared by both mesh constructors.
- `generate.jl`: command-line input, construction and export.

Load `bearing.jl` and `grid.jl` before loading `surface_feed.jl`. The inlet
imports `BearingParameters` from its caller and `background_grid` / `patch_axis`
from that same `grid.jl`. It does not maintain another parameter type, background
constructor or spacing implementation.
Construction loads neither the exporter nor experiments, Gmsh or fluid solvers.

## Generate a mesh

From `journal-bearings/`:

```sh
julia --project=. mesh/2d/conical/generate.jl \
  --geometry surface-feed --journal-radius-mm 50 --length-mm 100 \
  --clearance-mm 0.05 --cone-angle-deg 10 --eccentricity-ratio 0.6 \
  --n-theta 512 --n-axial 256 --feed-diameter-mm 5 --feed-z-mm 50 \
  --output out/my-inlet.vtk
```

This writes the requested inlet mesh and rejects an existing output file.
For a regular mesh use `--geometry unported` and omit the feed flags. Regular
meshes also support `--experiment NAME`, creating a mesh and settings file.
Named fitted cases are deliberately rejected: the case loader outside `mesh/`
still assumes the old fixed layout. Use `--output` for the new direct meshes.
`--help` lists the options.

Runnable fixed-input examples:
- `mesh/examples/mesh_unported.jl`: uniform spacing, unchanged default.
- `mesh/examples/mesh_adjusted_unported.jl`: the inlet's initial spacing without
  an opening, on a 128×64 background. Refuses to overwrite its output.
- `mesh/examples/mesh_inlet.jl`: the direct fitted inlet.
- `mesh/examples/compare_inlet_counts.jl`: generate six background-count choices
  and a shape-comparison table in a fresh output directory.
- `mesh/examples/compare_inlet_geometry.jl`: vary bearing dimensions and inlet
  geometry, recording accepted meshes and rejected requests.

## Current inlet limits

The constructor uses the supplied `BearingParameters`, inlet diameter and
axial position. The drilling direction remains +Y. Diameter must be positive
and finite; the complete opening must lie strictly between z=0 and z=L.
Eccentricity changes the reported gap, not the housing vertices.

The opening must also leave room for the current transition. Its patch
half-width is 1.5 times the diameter D; the surrounding spacing support spans
z ± 3D and must fit strictly between the ends and away from the angular seam.
Enough grid intervals must remain around that support. A physically possible
opening can therefore still be rejected by this construction. It is never
silently moved or resized. The exported axial end planes are anchored to 0 and
the supplied length in mm, avoiding overshoot from the internal unit conversion.

The cached `p.slope` is carried unchanged through construction, optimization
and measurements. Internal geometry records `R` (housing midpoint radius), `L`
and `feed_z_m` in metres, `ca = cosd(p.semicone_angle_deg)`, and the dimensionless
`slope`. Rebuilding geometry after optimization retains these supplied values.

`n_theta` and `n_axial` now construct the actual background. The inlet replaces
one patch; the resulting cell array is not rectangular. The optimizer runs on
this requested mesh, not on a hidden coarse mesh. There is no count whitelist,
but a background must leave enough room for the transition and pass the quality
checks. Both odd and even counts are supported where the construction succeeds.

The patch count follows the square root of patch area divided by background
cell area at the inlet, rounded to an even count and limited by available space.
The transition uses round(9*patch_cells/8) layers, retaining the earlier chosen
layer/patch ratio. This is a construction rule, not a geometric necessity.
Inspect `surface.patch_cells` and `surface.transition_layers` for the result.

Examples with Rj=50 mm, L=100 mm, clearance=0.05 mm, angle=10°, and a 5 mm
midpoint inlet:

| Background counts | Surface cells | Inlet edges | Vertices |
|-------------------|--------------:|------------:|---------:|
| 128 x 64          |         8,416 |          32 |    8,559 |
| 256 x 128         |        33,664 |          64 |   33,951 |
| 512 x 256         |       134,252 |         120 |  134,823 |
| 640 x 208         |       136,300 |         120 |  136,999 |
| 800 x 256         |       209,892 |         152 |  210,767 |

These intentionally differ from the older coarse-map subdivision results.
The new layout tag is `fitted-land-direct-v1`; it is not relabelled as the old
`controlled-fitted-land-v1` layout.

Rim vertices lie on the drilling/cone intersection. Edges between them are
straight segments. These meshes have an opening, not a resolved 3D passage.
Acceptance requires minimum physical corner score >0.9, minimum physical
face alignment >0.9, and maximum neighbouring **chart-area** ratio <2.
Edge ratio is reported, not gated. Optimization uses chart angles and areas;
acceptance checks do not claim that it minimizes a physical-space objective.
The settings and 180-iteration cap are unchanged; passing does not establish
an optimum. These are surface checks, not 3D Hex8 quality guarantees.

## Calling the inlet constructor

```julia
include("bearing.jl")
include("mesh/2d/conical/grid.jl")
include("mesh/2d/conical/surface_feed.jl")
include("mesh/vtk.jl") # only needed for export

p = BearingParameters(journal_midpoint_radius=60., bearing_length=120.,
    radial_clearance=.08, semicone_angle_deg=15., eccentricity_ratio=.6)
mesh = surface_feed_grid(p, 640, 208; feed_diameter=6., feed_z=70.3)
write_vtk("fresh-mesh.vtk", mesh)
```

For `housing_grid`, load `grid.jl` without `surface_feed.jl`.

## Uniform and adjusted backgrounds

`background_grid(theta_nodes, z_nodes)` creates `(theta,z)` vertices and periodic
quadrilateral connections. Angles are radians; z is mm. The angular vector must
contain at least three strictly increasing values in `[0,2pi)`, with no duplicate
`2pi` endpoint. The axial vector needs at least two strictly increasing values.

`housing_grid(p, n_theta, n_axial)` supplies uniform positions by default. To
adjust spacing, supply either or both position vectors:

```julia
mesh = housing_grid(p, n_theta, n_axial; theta_nodes, z_nodes)
```

Their lengths must be `n_theta` and `n_axial+1`; the axial positions must span
exactly `0` to `p.bearing_length`. The constructor recalculates each cell's gap
and physical area using its own midpoint and widths, rather than assuming a
uniform cell size. Input position vectors are not modified.

`patch_axis` constructs the current piecewise-uniform inlet spacing. Both the
inlet and adjusted-unported example use this function. Sharing it does not make
its spacing mandatory for the unported case. The fitted inlet derives that
spacing from the supplied bearing and inlet geometry.

Returned inlet fields:
- `nodes`: Cartesian housing vertices in mm.
- `quads`: four vertex IDs per solid-surface cell; no cells cover the opening.
- `feed_edges`: pairs of vertex IDs around the opening. There is no `feed_cells` mask.
- `radial_gaps`: cell-centre radial gaps in mm.
- `areas`: physical conical surface areas in mm².
- `quality`: surface corner, face-alignment and spacing measurements.
- `surface`: the internal surface coordinates in metres and face connectivity.
  Its `areas` are coordinate-chart areas, not the physical `areas / 1e6`.

Do not reshape inlet cells into `(n_theta,n_axial)`. The fitted fluid calculation
uses `surface` and its face geometry, not the regular-grid stencil. Fluid
operators outside `mesh/` may still assume reference geometry; this meshing
interface does not make those operators geometry-aware.

VTK exports place inlet `LINES` before surface `POLYGONS`, with cell data in
that order. Rim lines have `feed_rim=1`, zero area and NaN gap/corner values.
They are annotations, not additional surface cells. Regular exports may write
to an existing temporary file, as required by the fluid runner's atomic export;
inlet exports reject existing files. The CLI rejects existing outputs for both.

## Tests

Run each command as a separate Julia process from `journal-bearings/`:

```sh
julia --project=. --check-bounds=yes mesh/tests/meshing.jl
julia --project=. --check-bounds=yes mesh/tests/mesh_cli.jl
julia --project=. --check-bounds=yes mesh/tests/inlet_counts.jl
julia --project=. --check-bounds=yes mesh/tests/inlet_geometry.jl
julia --project=. --check-bounds=yes mesh/tests/chart_geometry.jl
```

The first checks shared connectivity, uniform/adjusted geometry, input validation,
units, export and dependency boundaries.
The second checks command-line inputs, saved settings, hashes and overwrite rules.
The third checks non-whitelisted/odd counts, topology and the cached objective's
analytic gradient against finite differences. The fourth checks supplied
geometry through every stage: cone/rim equations, independently integrated
opening area, Cartesian tangent-based corner/alignment measurements, movable
node selection and eccentricity-independent housing geometry. The fifth checks
signed areas and centroids against 256-bit references, translation behaviour,
seam unwrapping, degenerate/reversed polygons, exact assembly/optimizer area
agreement, and the zero initial area-change penalty.
These tests load no fluid solvers or Gmsh. Objective caching reuses fixed-cell
areas and touches only affected cells/faces; it does not change the objective
or its derivatives.

The optional fitted fluid integration check is separate:

```sh
julia --project=. --check-bounds=yes mesh/tests/meshing_integration.jl
```

That test loads the older experimental fluid operator and its Gmsh dependency;
it is not part of mesh-only checks. Pre-direct inlet snapshots are not a
reference for the direct construction.
Parameterization also adds geometry metadata and consistently uses `p.slope`;
older inlet hashes therefore need a separately reviewed comparison, not a
relaxed snapshot check.

For a future structural refactor, record a direct-mesh baseline before editing,
then compare afterward:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. mesh/tests/meshing_snapshot.jl record out/mesh-baseline.toml
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. mesh/tests/meshing_snapshot.jl compare out/mesh-baseline.toml
```

The snapshot hashes every returned field and exact VTK bytes for the three
example inlet counts and regular cylinders/cones. To compare only the unchanged
regular meshes against an older mixed baseline, append `--regular-only`.
Old inlet snapshots are not a reference answer for the new direct construction.
Generated baselines stay outside Git.

## Compare cell shapes

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. mesh/examples/compare_inlet_counts.jl FRESH_DIRECTORY
```

The script exports accepted VTKs and `comparison.tsv`, including actual counts,
quality, runtime and the physical-area-weighted mean longest/shortest edge ratio.
Ratio 1 means equal edge lengths, not necessarily right angles. The separate
corner/alignment checks still matter. The first timing includes compilation.
Quality-rejected requests are recorded rather than accepted silently. Better
cell shapes do not, by themselves, establish more accurate fluid predictions.

## Compare bearing and inlet geometry

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. mesh/examples/compare_inlet_geometry.jl FRESH_DIRECTORY
```

The checked run exported 19 accepted examples and recorded three rejected
requests. The inputs vary radius (40/50/60 mm), length (80/90/100/120/123.4 mm),
clearance (0.05/0.06/0.08/0.1 mm), angle (0/5/10/12/15°), diameter (3/4/5/6/8 mm),
and axial position. These are individual variations and specific combinations,
not a claim that every combination or every valid bearing can be meshed.
Zero angle is checked as the cylindrical limit of this construction.

Two combined examples:

| Rj / L / clearance (mm) | Angle | D / inlet z (mm) | Background | Cells | Corner min | Alignment min | Chart-area ratio max |
|------------------------|------:|-----------------:|------------|------:|-----------:|--------------:|---------------------:|
| 45 / 90 / 0.06          |   12° | 4 / 35.3         | 127 x 65   | 8,387 | 0.909584 | 0.989899 | 1.630901 |
| 60 / 120 / 0.08         |   15° | 6 / 70.3         | 640 x 208  | 136,704 | 0.917702 | 0.993392 | 1.469237 |

An opening touching an end, a fitting opening whose transition crosses an end,
and an oversized opening are rejected with their reasons. Geometry checks also
verify the cone and drilling equations, physical area and gaps, the supplied
optimization centre, and metric measurements against independent 3D tangents.
This is mesh verification, not fluid-model validation.

Generated evidence is under
`out/cases/conical-inlet-stable-chart-geometry/`; checks and snapshots are under
`out/reviews/stable-chart-geometry/`. These outputs are ignored.

## Parameterization checkpoint

At checkpoint `713cedc`, compared with `a72cf56`, the 128×64, 256×128 and 512×256 reference
inlets retain bitwise-identical chart coordinates, connectivity, radial gaps
and physical areas. Initial objectives and gradients are identical at all three
counts; the full 128×64 optimization history is identical too.
Using the cached slope consistently changes Cartesian coordinates by at most
7.11e-15 mm, with roundoff-level changes in physical quality measurements.
Surface metadata adds `slope` and `feed_z_m`.

Those exact inlet comparisons describe parameterization alone, before the
area/centroid correction below. They are not a claim that today's inlet meshes
match that earlier snapshot.

## Shared chart areas and stable centroids

`inlet/metrics.jl` owns `cell_chart_area`, used by assembly and optimization.
It subtracts the first vertex before forming cross products; the centroid uses
first moments about the same vertex and then restores the origin. Callers
unwrap seam-crossing cells first. Areas remain signed, so assembly still rejects
reversed or collapsed cells rather than hiding their orientation with abs(area).

Assembly areas now exactly match the optimizer's areas at its initial mesh.
The initial area-change anchor is exactly zero. At the previously worst sampled
512×256 cell, relative area error against a 256-bit reference falls from
4.93e-11 to 4.54e-17; centroid error falls from 1.34e-12 m to 4.59e-18 m.

The optimized mesh is not unchanged. Correcting the reference areas changes
the optimizer's numerical inputs even though its formulas, weights, tolerances
and 180-iteration cap stay fixed. Across 21 before/after cases, connectivity and
rim vertices remain exact, but maximum vertex displacement reaches 13.16 µm
and maximum individual physical cell-area change reaches 1.11%. All 19 geometry
examples still pass the same quality gates; the three unavailable requests
still reject. The compared runs do not establish optimizer convergence or
improved fluid accuracy.

In two controlled comparisons (128×64 with D=8 mm and 512×256 with D=5 mm),
changing only the old initial area array reproduces the new optimized chart
coordinates and full iteration history exactly. Both old and new runs reach
the iteration cap. This isolates the layout changes from the centroid fix.
Measured changes are saved in `out/reviews/stable-chart-geometry/checkpoint.tsv`.

All regular-mesh fields and VTK bytes match the earlier baseline. The current
nine-mesh snapshot is also exactly repeatable in all fields and VTK bytes.
Earlier snapshots are retained; their inlet equality checks have not been
relaxed to accommodate the numerical changes.
