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
imports `background_grid` and `patch_axis` from that same `grid.jl`; it does
not maintain another background constructor or spacing implementation.
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

## Current inlet limits

The inlet construction supports only journal radius 50 mm, length 100 mm,
reference radial clearance 0.05 mm, semi-cone angle 10 degrees, and a 5 mm
+Y drilling at z=50 mm. Other layouts are rejected. Eccentricity may vary:
it changes the reported gap, not the fixed housing vertices.

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

Examples from the direct construction:

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
The quality thresholds are surface checks, not 3D Hex8 quality guarantees.

## Calling the inlet constructor

```julia
include("bearing.jl")
include("mesh/2d/conical/grid.jl")
include("mesh/2d/conical/surface_feed.jl")
include("mesh/vtk.jl") # only needed for export

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
mesh = surface_feed_grid(p, 512, 256; feed_diameter=5., feed_z=50.)
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
its spacing mandatory for the unported case. The fitted inlet still uses the
fixed physical dimensions listed above, but no longer restricts the count pair.

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
uses `surface` and its face geometry, not the regular-grid stencil.

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
```

The first checks shared connectivity, uniform/adjusted geometry, input validation,
units, export and dependency boundaries.
The second checks command-line inputs, saved settings, hashes and overwrite rules.
The third checks non-whitelisted/odd counts, topology and the cached objective's
analytic gradient against finite differences. These tests load no fluid solvers
or Gmsh. Objective caching reuses fixed-cell areas and touches only affected
cells/faces; it does not change the objective or its derivatives.

The optional fitted fluid integration check is separate:

```sh
julia --project=. --check-bounds=yes mesh/tests/meshing_integration.jl
```

That test loads the older experimental fluid operator and its Gmsh dependency;
it is not part of mesh-only checks. Historical bitwise geometry comparisons no
longer apply because direct construction intentionally changes the mesh.

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
