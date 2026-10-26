# Conical surface meshes

## Files

- `grid.jl`: `housing_grid`, a regular circumferential/axial mesh without an opening.
- `surface_feed.jl`: `surface_feed_grid`, quadrilateral cells around an inlet opening.
- `inlet/seed.jl`: initial background, opening rings and surrounding-node smoothing.
- `inlet/placement.jl`: local node optimization and subdivision of the coarse cells.
- `inlet/metrics.jl`: cell/face connectivity, geometric measurements and quality.
- `inlet/controlled.jl`: runs that inlet construction and checks its quality limits.
- `../../vtk.jl`: `write_vtk`, shared by both mesh constructors.
- `generate.jl`: command-line input, construction and export.

The constructors expect `bearing.jl` to be loaded by their caller. Neither
constructor loads the exporter, experiments, Gmsh or fluid solvers. The inlet
constructor does not load or call `housing_grid`.

## Generate a mesh

From `journal-bearings/`:

```sh
julia --project=. mesh/2d/conical/generate.jl \
  --geometry surface-feed --journal-radius-mm 50 --length-mm 100 \
  --clearance-mm 0.05 --cone-angle-deg 10 --eccentricity-ratio 0.6 \
  --n-theta 512 --n-axial 256 --feed-diameter-mm 5 --feed-z-mm 50 \
  --experiment my-inlet
```

This creates `out/cases/my-inlet/{mesh.vtk,experiment.toml}`. Existing case
names are rejected. Alternatively use `--output fresh-mesh.vtk` without
`--experiment`. For a regular mesh use `--geometry unported` and omit the two
feed flags. `--help` lists the options.

Runnable fixed-input examples are `examples/mesh_unported.jl` and
`examples/mesh_inlet.jl`, outside the implementation directory.

## Current inlet limits

The inlet construction supports only journal radius 50 mm, length 100 mm,
reference radial clearance 0.05 mm, semi-cone angle 10 degrees, and a 5 mm
+Y drilling at z=50 mm. Other layouts are rejected. Eccentricity may vary:
it changes the reported gap, not the fixed housing vertices.

The two count arguments select one of three subdivisions of the same optimized
coarse layout. They are not dimensions of the returned cell array.

| Background counts | Surface cells | Inlet edges | Vertices |
|---|---:|---:|---:|
| 128 x 64 | 8,416 | 32 | 8,559 |
| 256 x 128 | 33,664 | 64 | 33,951 |
| 512 x 256 | 134,656 | 128 | 135,231 |

Rim vertices lie on the drilling/cone intersection. Edges between them are
straight segments. These meshes have an opening, not a resolved 3D passage.
The quality thresholds are surface checks, not 3D Hex8 quality guarantees.

## Calling the inlet constructor

```julia
include("bearing.jl")
include("mesh/2d/conical/surface_feed.jl")
include("mesh/vtk.jl") # only needed for export

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)
mesh = surface_feed_grid(p, 512, 256; feed_diameter=5., feed_z=50.)
write_vtk("fresh-mesh.vtk", mesh)
```

For `housing_grid`, load `grid.jl` instead of `surface_feed.jl`.

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
julia --project=. --check-bounds=yes test/meshing.jl
julia --project=. --check-bounds=yes test/mesh_cli.jl
```

The first checks geometry, connectivity, units, export and dependency boundaries.
The second checks command-line inputs, saved settings, hashes and overwrite rules.
Neither loads Gmsh or fluid solvers.

Historical mesh equivalence and the fitted fluid integration check are separate:

```sh
julia --project=. --check-bounds=yes test/meshing_integration.jl
```

That test intentionally loads the older experiments and their Gmsh dependency.

For structural refactors, record a baseline before editing, then compare after:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. test/meshing_snapshot.jl record out/mesh-baseline.toml
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=. test/meshing_snapshot.jl compare out/mesh-baseline.toml
```

The snapshot hashes every returned mesh field and exact VTK bytes for all three
inlet resolutions and regular cylinders/cones at the same resolutions. The
baseline stays outside Git; source can be recovered from the local checkpoint.
Task progress and decisions belong in `todos/meshing.org`.
