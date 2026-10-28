# Cylindrical and conical surface meshing

This shared mesher builds periodic quadrilateral meshes on cylindrical or
conical housing surfaces, either unported or fitted around a circular +Y drilling.
It does not build a 3D passage, validate fluid predictions, or claim a globally
optimal node layout.

## API

Run from `journal-bearings/`, with lengths in **mm**:

```julia
include("bearing.jl")
include("mesh/2d/grid.jl")
include("mesh/2d/surface_feed.jl") # only needed for an inlet
include("mesh/vtk.jl")             # only needed for export

p = BearingParameters(journal_midpoint_radius=50., bearing_length=100.,
    radial_clearance=.05, semicone_angle_deg=10., eccentricity_ratio=.6)

plain = housing_grid(p, 128, 64)
fitted = surface_feed_grid(p, 128, 64; feed_diameter=5., feed_z=50.)
write_vtk("out/inlet.vtk", fitted)
```

Set `semicone_angle_deg=0.` for a cylinder; the same constructors, fitted-inlet
placement and acceptance checks apply. Nonzero angles produce cones. The CLI
uses `--cone-angle-deg 0` for a cylinder.

Both grids expose `nodes`, `quads`, `radial_gaps`, `areas`, `n_theta`, `n_axial`,
`quality` and `acceptance`. Cartesian points and gaps are mm; physical areas are
mm². Fitted grids additionally expose `feed_edges`, `surface`, `placement` and
`layout="fitted-land-direct-v1"`. The rim is a boundary, not a feed-cell mask.
Unported grids also retain `parameters`, `theta_nodes` (radians) and `z_nodes`
(mm). These are owned coordinate copies for mesh-aware tensor solvers.

`n_theta` and `n_axial` are actual requested background counts. The fitted inlet
replaces a background patch, so its final quads do not form a rectangular array.
There is no resolution whitelist, hidden coarse solve or subdivision.

### Acceptance

A returned public grid has passed **all** applicable checks:

| Check | Requirement |
|---|---|
| Geometry | Finite coordinates, positive gaps/areas and valid orientation |
| Physical corner score | Minimum > 0.9 |
| Physical face alignment | Minimum > 0.9 |
| Neighbouring **chart-area** ratio | Maximum < 2 |
| Fitted placement | Complete physical-motion window, as described below |

Edge ratio and maximum nonorthogonality are reported diagnostics. The edge
ratio is not an additional gate. These are surface metrics, not Hex8 quality.

`acceptance` contains `valid_geometry`, `quality_pass`, `placement_stable`,
`termination_reason`, and `accepted` (the AND of the three checks). Public
constructors reject unsuccessful acceptance rather than return a bad/unsettled
mesh, resize the opening, relocate it or substitute another resolution.

### Fitted placement

L-BFGS moves interior/support nodes while keeping the inlet rim, axial ends,
periodic seam and distant background fixed. It retains the original area
anchors. Current defaults have one owner in `inlet/controlled.jl`:

```julia
surface_feed_grid(p, 128, 64; feed_diameter=5., feed_z=50.,
    movement_tolerance=0.001, stability_steps=20, max_iterations=5000)
```

For **each** movable node, sum its physical 3D chord movements over 20 accepted
steps, divided by its **original shortest adjacent physical edge**. The maximum
per-node path ratio must be at most 0.001 (0.1%). Windows do not overlap. Path
length prevents cancellation; one tiny backtracked move cannot fill a window.
The candidate must then pass the independent quality checks.

The gradient and objective remain diagnostic: a tiny gradient or microscopic
score decrease is not compulsory for practical acceptance. Exact stationary
states need no fictitious iterations. An exactly unchanged trial with nonzero
gradient stops as `numerical_stagnation`, not accepted progress. Work exhaustion,
failed search, stagnation and failed quality are not successful acceptance.
`placement` retains accepted count, reason, window length and measured ratio.

The local-size idea follows established
[CUBIT smoothing practice](https://cubit.sandia.gov/files/cubit/17.08/help_manual/WebHelp/mesh_generation/mesh_modification/mesh_smoothing/mesh_smoothing.htm).
The 0.001 scale is motivated by
[CGAL's documented local-motion tolerance](https://doc.cgal.org/6.1/Mesh_2/group__PkgMesh2Functions.html).
Our cumulative 20-step window and original per-node scales are **engineering
choices**, not a literal copy of CGAL's implementation or a fluid-error bound.
Window stability is not a promise of zero movement under infinite continuation.
The 5000-step cap is only a safeguard, not a convergence criterion.

### Unported spacing

`housing_grid` directly constructs an orthogonal tensor chart. No optimization
or artificial stability iterations are performed; it reports
`direct_construction`. Corner and alignment are analytically 1; actual spacings
determine chart-area and physical edge ratios. These formulas are independently
checked against general quad assembly.

Optional `theta_nodes` (radians) and `z_nodes` (mm) adjust spacing. They must have
lengths `n_theta` and `n_axial+1`; angles increase strictly in `[0,2pi)` without
a duplicate endpoint, each angular cell spans less than pi, and axial positions
increase strictly from exactly 0 to `p.bearing_length`. Inputs are copied, not
mutated. Bad geometry or quality rejects instead of being silently graded.
`patch_axis` supplies the same adjusted spacing used by the inlet seed.

## Geometry and units

The geometry owner is `bearing.jl`: `bearing_geometry`, `journal_radius`,
`housing_radius`, `bearing_radial_geometry` and `bearing_point`. Housing points
are independent of journal eccentricity; eccentricity changes the gaps.

Internal `surface.nodes` are `(R_ref*theta, z/cos(alpha))` in **metres**.
`surface.areas` are chart areas, not physical cone areas. Corner/alignment checks
use the cone metric; public areas integrate its radius scale. `geometry_mm`
retains original dimensions and slope, and axial endpoint conversion preserves
exact supplied end planes. Do not recover analytic gaps/radii with Cartesian
`hypot`, or use chart areas as physical areas.

The inlet diameter must be positive; the entire opening must lie strictly
between the ends. The current seed uses a patch half-width of 1.5D and spacing
support roughly z±3D, which must fit inside the ends and angular seam. Unavailable
requests reject, even when a differently constructed physical inlet might fit.
The +Y drilling direction is fixed; no passage pressure condition is prescribed.

## Files and dependency wiring

- `grid.jl`: background connectivity, adjusted axes and unported construction.
- `quality.jl`: shared metric module, acceptance, physical movement and tensor metrics.
- `metrics.jl`: shared chart polygons, face assembly and physical surface quality.
- `inlet/seed.jl`: fitted rings and harmonic background support.
- `inlet/placement.jl`: L-BFGS, practical stopping and optional diagnostics.
- `inlet/controlled.jl`: inlet defaults, construction and acceptance.
- `surface_feed.jl`: public fitted mapping, physical areas and gaps.
- `../vtk.jl`: shared VTK exporter, separate from construction.
- `generate.jl`: CLI input, construction and export.

Load `bearing.jl`, then `grid.jl`, then the optional `surface_feed.jl`. They share
one parameter/geometry definition and one background implementation. Construction
loads neither exporters, experimental operators, Gmsh nor fluid solvers.

## CLI and examples

```sh
julia --project=. mesh/2d/generate.jl \
  --geometry surface-feed --journal-radius-mm 50 --length-mm 100 \
  --clearance-mm .05 --cone-angle-deg 10 --eccentricity-ratio .6 \
  --n-theta 128 --n-axial 64 --feed-diameter-mm 5 --feed-z-mm 50 \
  --output out/inlet.vtk
```

Use `--geometry unported` without feed flags for a complete surface. CLI output
reports quality and acceptance for both paths. Output files must be fresh.
Unported meshes also support `--experiment NAME`. Named fitted cases deliberately
reject: the case loader does not support `fitted-land-direct-v1`. Use `--output`
for fitted surfaces; see [solver support](../../docs/solver-support.md).

`mesh/examples/mesh_{unported,adjusted_unported,inlet}.jl` are standalone examples
with an optional fresh VTK path. `compare_inlet_counts.jl` and
`compare_inlet_geometry.jl` record accepted and rejected requests.

The VTK writer permits unported overwrites for atomic simulation exports.
Its optional `feed_cells` scalar serves regular-grid simulation exporters. It is rejected for fitted rim meshes and is
not an alternative inlet construction.

## Tests

```sh
julia --project=. mesh/tests/runtests.jl
```

Checks cover geometry, mesh quality, placement convergence, boundary invariants
and command-line exports. Individual suites can also be run directly:

```sh
julia --project=. mesh/tests/practical_acceptance.jl
```

To record and compare mesh snapshots, use a **fresh** baseline:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
julia --project=. --threads=1 --gcthreads=1 mesh/tests/meshing_snapshot.jl record out/snapshot.toml
julia --project=. --threads=1 --gcthreads=1 mesh/tests/meshing_snapshot.jl compare out/snapshot.toml
```

Use the same threading controls for both commands: fitted optimization is
sensitive to floating-point differences. Comparison checks every public mesh
field and the exported VTK bytes against the saved baseline.

## Optional numerical diagnostics

`controlled_layout` is private construction machinery. `trace_iterations`
selects owned search snapshots; `observer` receives owned initial/update/final
states including window metrics and final acceptance. Observation does not
change numerical decisions; callback exceptions propagate.

The fixed-budget examples `compare_inlet_iterations.jl`,
`diagnose_inlet_optimizer.jl` and `inlet_long_run.jl` explicitly select
`stability_steps=0` for gradient/work-budget investigations. This
mode cannot bypass public mesh acceptance. Diagnostic snapshots preserve original
seeds/area anchors; they are not an API for restarting a rebuilt terminal mesh.

`queue_inlet_long_runs.py` is an optional Linux resource-aware study queue. It freezes all numerical sources, limits independent single-threaded
workers, and leaves other workloads untouched. DEBUG writes plain `key=value`
stderr; disk telemetry remains available without DEBUG. Resource/time limits
control work, never numerical success. A finished process is not automatically
an accepted mesh. No queue or long numerical study is part of normal generation.
