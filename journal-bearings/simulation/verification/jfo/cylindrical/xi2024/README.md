# Xi2024 finite-cylinder pressure comparison

Compare the pressure curve in Xi, Deng and Li (2024),
[doi10.1007/s40544-023-0791-1](https://doi.org/10.1007/s40544-023-0791-1),
§3.3/Fig.7, using the shared production Reynolds/JFO operators.

**Gap-origin limitation:** §3.3 gives the operating inputs but omits the
journal-gap equation and angular origin, referring to Rom2023. This fixture
uses the fixed conventional `h=c*(1+.8*cos(2pi*x/.08))`, maximum gap at the
supply seam and minimum at40mm. This convention is not fitted. The original
Rom geometry definition remains unrecovered, so agreement is conditional on
this declared assumption, not an exact author-implementation reproduction.

## Run

From `journal-bearings/`, with a fresh study directory:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. simulation/tests/runtests.jl xi2024_jfo_checks.jl
julia --project=. simulation/verification/jfo/cylindrical/xi2024/run.jl out/verification/jfo-xi2024
python simulation/verification/jfo/cylindrical/xi2024/analyse.py out/verification/jfo-xi2024
python simulation/verification/jfo/cylindrical/xi2024/plot.py out/verification/jfo-xi2024
```

Use `--help`, `-v` for diagnostics or `-q` for errors only. The numerical study
needs no PDF, downloaded archive or earlier review output.

## Published inputs and boundaries

`reference.toml` separates published inputs from fixture assumptions.
The unwrapped circumference is80mm (Fig.7 and Eq.18), axial length20mm,
clearance17.5μm, eccentricity ratio.8, viscosity.014Pa·s and wall speed.2m/s.
The derived nominal radius is`.08/(2pi)`m and speed150rpm.
The profile is taken at the actual axial midplane,`y=.01m`.

**All four edges** have0.1MPa absolute pressure; cavity pressure is0absolute.
Equivalently, ambient is0gauge and cavity-0.1MPagauge. They must not both be
set to zero. This is different from the1957 ungrooved, periodic-only benchmark.
There is no drilled hole or finite-width reservoir-cell column.

`cases.jl` constructs the periodic production pressure/transport operators,
retains every internal coefficient and replaces the seam connection with two
half-cell pressure faces. Incoming liquid content is1; outgoing transport uses
the neighbouring film content. This is a verification-only boundary adapter,
not a new public production boundary option or a second interior discretization.

Continuity includes circumferential inlet/outlet AND both axial ends. Net axial
outflow need not vanish: the angular reservoirs supply its difference.
Both signs are preserved when computing gross boundary exchange.

The classical cosine-gap/single-radius fixture and the unchanged production
circle-intersection/separate-radius law are reported separately. Each runs
64×16→128×32→256×64→512×128→1024×256, with exact doubling.

## Reference data and reading limits

`reference.csv` contains graph-derived absolute pressure values:

-24 isolated red markers labelled “Work of Rom”, as reproduced in Xi2024;
-187 visible black-curve points labelled “Current work” in Xi2024.

These are not author raw arrays, and the reproduced Rom markers must not be
called independently recovered original FD results. Two merged red components
are excluded; hidden portions of the black curve are not invented.
Near-zero raster estimates are retained, not clipped to the cavity pressure.

The extraction uses the original477×392pixel Fig.7 image and fixed axis-tick
calibration. The one-pixel coordinate-reading halfwidths are about0.2035mm and
4.067kPa. They describe graph reading, NOT PINN discretization error or
experimental uncertainty. Analysis reports sampled pressure differences and
intersection with these reading boxes; published agreement is not a numerical
pass/fail tolerance.

To regenerate the fixture, obtain the cited PDF, extract page7's embedded
images with `pdfimages -f 7 -l 7 -png`, place its477×392 Fig.7 as `figure7.png`
in a fresh directory, and run:

```sh
python simulation/verification/jfo/cylindrical/xi2024/digitize.py DIRECTORY
```

Review `digitization-overlay.png` before using its CSV. `reference.toml` records
the source PDF/image hashes, tick coordinates and extraction exclusions.

## Numerical checks and outputs

- Independent continuous full-liquid manufactured pressure with four constant
  boundary traces; uniform/graded coordinates and both rotations.
- Arbitrary pressure/content matrix-versus-face balance, half-cell factors,
  boundary transport reversal, stationary limit and actual-plane sampling.
- Independent nodal expanded FD, with a fixed supply NODE representing the two
  identical angular Dirichlet traces. Its separate projected GS solve runs
  through512×128; damping/stopping and same-FD-matrix policy checks test algebra.
- The earlier coupled cavitating manufactured fields remain in the
  [periodic-cylinder study](../README.md). The new full-liquid MMS specifically
  checks the boundary adapter; it is not a new discontinuous-content oracle.

`inputs.toml` fingerprints inputs/source; `results.tsv` records physical
refinement and four-edge continuity. `profiles.tsv` includes prescribed end
traces and actual-midplane pressure/content. `manufactured.tsv`,
`independent.tsv` and `algebraic.tsv` keep different checks separate.
`comparison.tsv`/`analysis.json` report graph differences without altering
physical inputs, targets or reading widths. Finest fields are written as TSV
and nominal-cylinder VTK; plots are `comparison.png`, `verification.png` and
`fields.png`.

No published liquid-content or flow array is available for this case. Those
quantities are numerically checked, not claimed matched to Fig.7. Thermal,
compressible, load-equilibrium and FSI examples are outside this study.
