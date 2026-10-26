# Conical housing surfaces

`surface_feed_grid` now uses the **controlled fitted inlet**, not the old
cell-centre staircase. This promotion was approved after the
[three-way transition comparison](../../../docs/conical-controlled-transition.md).
`housing_grid` remains the unchanged general unported construction.

## Generate a new fitted case

```sh
julia --project=. mesh/2d/conical/generate.jl \
  --geometry surface-feed --journal-radius-mm 50 --length-mm 100 \
  --clearance-mm 0.05 --cone-angle-deg 10 --eccentricity-ratio 0.6 \
  --n-theta 512 --n-axial 256 --feed-diameter-mm 5 --feed-z-mm 50 \
  --experiment my-fitted-inlet
```

Choose a fresh name. The generator writes `out/cases/NAME/mesh.vtk` and
`experiment.toml`. The standalone `surface_feed.jl` example also uses a new
fitted-specific name and refuses to overwrite an existing output.

The supported inlet geometry is currently **Rj=50 mm, L=100 mm, c=0.05 mm,
semi-angle=10°, mouth diameter=5 mm, mouth centre z=50 mm on +Y**. Eccentricity
may vary within the bearing parameter limits: the housing mesh is unchanged,
while the gap diagnostic changes. Other shapes/locations are explicitly
rejected, not generated with hard-coded reference coordinates or a hidden
staircase fallback. This is promotion of the checked case, not a claim of a
generalized mesher. Unported cylinders/cones retain their previous support.

| Background counts | Land quads | Rim edges | Nodes |
|---|---:|---:|---:|
| 128×64 | 8,416 | 32 | 8,559 |
| 256×128 | 33,664 | 64 | 33,951 |
| 512×256 | 134,656 | 128 | 135,231 |

These are the three checked refinement levels. The two counts are **background
resolution labels**, not the dimensions of a rectangular array of fitted cells.

## Public API and units

```julia
include("mesh/2d/conical/surface_feed.jl")
grid = surface_feed_grid(p, 512, 256; feed_diameter=5.0, feed_z=50.0)
write_vtk("fresh-mesh.vtk", grid)
```

- `nodes`: unique, periodic, Cartesian housing vertices, **mm**.
- `quads`: four node IDs per **land** cell. Reservoir interior is omitted.
- `feed_edges`: pairs of node IDs defining the mouth boundary. Nodes lie on the
  exact drilling/cone intersection; the connecting edges approximate it by chords.
- `areas`: physical cone areas for the polygonal chart cells, **mm²**.
- `radial_gaps`: cell-centre radial gaps, **mm**, for the supplied reference `p`.
- `quality`: surface corner, alignment and grading diagnostics—not Hex8 OQ.
- `surface`: the **SI reference chart** and conservative face connectivity used
  by the checked fitted-grid operator. Its `areas` are chart areas in m², before
  applying the cone metric; they are not the public physical `areas / 1e6`.

**Breaking change:** there is no `feed_cells` field on the fitted grid. Do not
reshape these cells into `(n_theta,n_axial)` or pass a fictitious all-false mask.
Inlet conditions belong on `feed_edges` / `surface.boundary == :feed`. This
geometry does not itself prescribe pressure, liquid content or a passage model.

The exported VTK POLYDATA contains land `POLYGONS` plus inlet `LINES`, tagged by
`feed_rim`. As required by VTK, cell data list the **rim lines first**, followed
by land polygons; gap/corner values on lines are NaN, line areas are zero.
The line annotations are not extra fluid/surface cells. `surface_corner_score`
is exported for visual inspection. Existing unported export format and the
solver's atomic temporary-file writing contract are preserved.

Named manifests record `layout = "controlled-fitted-land-v1"`, actual counts,
rim representation and quality. Historical schema-1 files without a layout
remain explicitly `cell-centre-mask` for fed cases. They are not silently
relabelled, rewritten or regenerated. The loader rejects inconsistent fitted
counts, geometry and representation metadata.

## Construction and solver boundary

The implementation is self-contained in `inlet/`:

- `seed.jl`: geometry-only extraction of the reviewed coarse topology; no
  Gmsh, reservoir core or 3D extrusion.
- `placement.jl`: the same bounded angle/log-area placement and fixed-map
  subdivision with the rim correction as the reviewed trial.
- `metrics.jl`: cone-chart face construction and surface-quality checks.
- `controlled.jl`: the checked recipe and regression bounds.

No production code imports experimental files. Experiments remain independent
historical controls. The test suite compares the promoted chart **exactly**
against the reviewed implementation at all three levels: coordinates,
connectivity, centroids, areas, face normals, lengths, neighbours and tags.
It also runs the existing corrected conservative operator on the promoted
coarse chart and recovers the held-pose field's force/flow.

**No solver was silently switched to distorted-grid two-point fluxes.**
The old rectangular-grid Reynolds/JFO/equilibrium routines remain explicitly
rectangular and retain their historical controls. The named Reynolds CLI still
rejects fed cases. A fitted load solve must use the corrected conservative
operator/face geometry, not infer a mask from the background counts. Build the
housing grid once from reference geometry; steady bush pose belongs to the
fluid model, not to mesh-node relocation.

This is a **2D surface** promotion, not a 3D passage or Hex8 sweep. The later
Hex8 minimum OQ ≥0.90 requirement remains separate and untested here. The former
stepped generator and relevant callers are archived at
`out/reviews/pre-surface-inlet-promotion/`; old case outputs are untouched.

## Checks and actual export

```sh
julia --project=. --check-bounds=yes mesh/2d/conical/checks.jl
julia --project=. --check-bounds=yes mesh/2d/conical/generate_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/run_checks.jl
```

**150 mesh + 42 generator/manifest + 22 named-run assertions pass.** All prior
708 numerical assertions were also rerun: **922 passing in total**.

These cover exact three-level promotion equivalence, topology/rim conformity,
units, analytic footprint-area refinement, quality, exports, operator integration,
unsupported inputs, manifest compatibility and the existing cylindrical runner.
The true mouth-area error from polygonal approximation falls
0.6413% → 0.1606% → 0.04018%; no claim of an exact curved edge is made.

Generated example:
`out/cases/conical-surface-feed-controlled-e0p6-d5mm/{mesh.vtk,experiment.toml}`.
It has illustrative **e/c=0.6**, offset **30 µm**, reference clearance **50 µm**;
cell-centre radial gaps are 20.00059–79.99997 µm. It is a geometry export, **not
an operating-point or equilibrium prediction**. The physical regression fixture
instead uses the previous held e/c=0.0186497 pose; these must not be conflated.

[Preview of the actual exported file](../../../out/reviews/surface-inlet-promotion/mesh-preview.png).
At the finest level, minimum corner=0.92980, face alignment=0.95365 and maximum
neighbour-area ratio=1.45097, unchanged from the reviewed candidate.
Generation and export took **24.92 s wall / 24.76 CPU-seconds**, with **0.701 GiB
peak RSS**, including Julia startup/compilation. This is one measured geometry
batch, not a fluid solve or a repeat timing study. Tests and plotting are excluded.
Evidence and source snapshots: `out/reviews/surface-inlet-promotion/`.
