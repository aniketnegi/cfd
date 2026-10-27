# Shared cylindrical/conical 2D mesher migration

## Outcome and scope

The canonical surface mesher is now [`mesh/2d/`](../mesh/2d/README.md).
Both `housing_grid` and `surface_feed_grid` use the same implementation:
`semicone_angle_deg=0` produces a cylinder; nonzero angles produce cones.
There is no duplicate cylinder implementation or compatibility wrapper at the
old location. External callers must update their include/CLI paths.

The migration moves all nine mesher files, preserving the `inlet/` subtree.
It updates live includes, CLI-relative dependencies, examples, tests, Julia
source-hash lists, Python frozen-worker lists, experiment source-copy roots,
and documentation. The diameter experiment's already-stale `inlet/metrics.jl`
manifest entry is corrected to the shared metrics/quality dependency locations.

Executable AST comparison confirms all seven construction/quality/placement
files are unchanged after removing comments, documentation and line numbers.
Public function names, parameters, units, mesh fields, acceptance thresholds,
placement defaults, manifest schemas and solver mathematics are unchanged.
The VTK exporter, including its legacy descriptive title, is unchanged so
exports can be compared byte-for-byte.

## Explicit cylindrical generation checks

Both independently varied cylindrical fixtures pass geometry, quality and
placement acceptance. The reference is Rj=50 mm, L=100 mm, c=0.05 mm,
eccentricity ratio=0.6; its fitted inlet has D=5 mm at z=50 mm. The varied
fixture uses Rj=40 mm, L=80 mm, c=0.08 mm, D=4 mm at z=35 mm.

| Cylindrical fixture | Background | Unported cells | Fitted land cells |
|---|---:|---:|---:|
| Reference | 128 × 64 | 8,192 | 8,416 |
| Varied geometry/inlet | 160 × 52 | 8,320 | 8,544 |

The new persistent `mesh/tests/cylindrical.jl` checks constant housing radius,
exact axial ends, independent ray/circle gaps, analytic cylinder area, positive
cell areas, periodic seam connectivity, fitted drill-rim coordinates, complete
boundary classification and an independently integrated drill footprint.
It also runs the real CLI in fresh processes from an unrelated working
directory and requires identical VTK bytes to direct API generation.

Four additional fresh CLI invocations produced retained VTK artifacts for
0°/10° × unported/fitted, all with accepted meshes. The reference cylindrical
fitted mesh has minimum corner score 0.937315, minimum face alignment 0.997809,
and maximum neighbouring chart-area ratio 1.517745; all pass existing gates.

## Verification

| Check | Result |
|---|---|
| Complete isolated mesh runner | **12/12 suites pass; 10,343 unit assertions** |
| New cylindrical API/CLI tests | **66/66 pass**, included in the mesh total |
| Matched-thread original/migrated snapshots | **52/52 pass across nine meshes**, every field and exact VTK bytes |
| Additional cylindrical fitted snapshots | **Both fixtures match exactly**, every field and VTK bytes |
| Executable construction/quality/placement ASTs | **7/7 unchanged** |
| Solver regression suites | **7/7 pass; 645 assertions** |
| Actual frozen-source workers | **2/2 complete**, zero running/queued |
| Affected literal include/source references | **147 checked; all exist** |
| Retained cylinder/cone CLI exports | **4/4 accepted** |

Solver assertions comprise cylindrical full-film (22), conical full-film (136),
conical JFO (171), equilibrium geometry/iteration (87), JFO slider (149), named
runner (22), and the existing Swift–Stieber experiment regressions (58).
The seven maintained files and all snapshot fields are compared without
relaxing tolerances, dropping fields, or changing production arithmetic.

Final numerical checks use fresh Julia processes with bounds checking and
single-threaded Julia/BLAS/OMP controls. Frozen workers use their own archived
source tree; their one-step diagnostic budgets test wiring, **not** mesh
acceptance. Actual accepted cylindrical meshes are checked separately above.

### Two verification issues investigated

1. The initial finest-cone snapshot comparison failed four assertions because
   the initial baseline and comparison had different effective BLAS threading.
   Probes showed `OPENBLAS_NUM_THREADS=1` alone left this Julia build at **eight**
   BLAS threads; adding `OMP_NUM_THREADS=1` gave **one**. Fine optimization
   amplified those floating-point differences. The untouched original source
   from Git checkpoint `5eb6287` was therefore frozen and rerun with exactly the
   same controls as the migrated source. The original baseline and failed logs
   remain preserved; no baseline was overwritten and no field was excluded.
   Matched-thread comparisons are recorded separately as `before-pinned.toml`
   and `after-pinned.log`. Snapshot documentation now makes the controls explicit.
2. The first version of the **new** cylinder test used a long Float64 sum with
   an overly tight total-area bound. Its two area assertions failed at relative
   errors around 1–2e-14 even though cell areas were correct. The independent
   oracle now accumulates those areas in BigFloat at the **same** tolerance.
   No production arithmetic, quality gate or pre-existing test was relaxed.
   The initial failed suite log is retained alongside the final full rerun.
   The redundant first runner was stopped after the corrected full rerun was
   underway; the final complete 12-suite run exits successfully.

## Integration boundaries and preserved work

- The only old mesh-path references remaining in maintained experiment code
  are the two deliberate historical-source checks in
  `expts/conical-inlet-transition/run.jl`. Their old archive paths/hashes stay
  intact. Existing output trees, frozen studies and historical reports were
  not globally rewritten.
- Five pre-existing **untracked** local files also received path-only fixes:
  `docs/conical-controlled-transition.md`, `expts/gangrade-figure6/run.jl`,
  `expts/conical-load-equilibrium/{run,resource_probe}.jl`, and
  `expts/conical-2d-local-quads/README.md`. Their contents are not swept into
  the migration commit. Before-copies and exact local patches are retained.
- The pre-existing root README simulation expansion remains unstaged. The
  commit contains independently reviewed shared-mesher documentation, not that
  unrelated expansion. Other pre-existing work is preserved.
- This does **not** repair or certify the historical fitted fluid/load adapter:
  its missing prerequisite include, obsolete `mesh.level` assumptions and
  reference-case force formulas remain separate simulation-audit issues.
  Named fitted cases remain deliberately unsupported by the external loader.
  Mesh acceptance is not fluid validation or proof of arbitrary-geometry
  solver support. No physical load study or historical force/flow rebaseline
  was performed.

## Evidence and reproduction

Detailed evidence is under
`out/reviews/shared-2d-migration.p69RgH/`: all logs and exit codes, both original
baselines, frozen original source, exact snapshots, accepted CLI VTK files,
threading probes, AST comparison, reference audit, local-only patches, worker
archives and `verification.json`.

From `journal-bearings/`:

```sh
julia --project=. mesh/tests/runtests.jl
```

For a standalone accepted cylindrical fitted mesh (use a fresh output path):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. mesh/2d/generate.jl \
  --geometry surface-feed --journal-radius-mm 50 --length-mm 100 \
  --clearance-mm .05 --cone-angle-deg 0 --eccentricity-ratio .6 \
  --n-theta 128 --n-axial 64 --feed-diameter-mm 5 --feed-z-mm 50 \
  --output out/cylinder-inlet.vtk
```

For an unported cylinder, use `--geometry unported` and omit both feed flags.
For a cone, change `--cone-angle-deg`; no different mesher is needed.
