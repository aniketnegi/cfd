# Gangrade conical-bearing reference

## Published inputs

Gangrade, Phalle and Mantha, *Performance Analysis of a Conical Hydrodynamic
Journal Bearing*, 43:559–573 (2019):
[DOI 10.1007/s40997-018-0217-2](https://doi.org/10.1007/s40997-018-0217-2).

| Quantity | Value | Source |
|---|---:|---|
| Mean journal radius | 50 mm | Table 1 |
| Axial length | 100 mm | Table 1 |
| Radial clearance | 50 µm | Table 1 |
| Semi-cone angle | 10° | Fig. 6 comparison; Fig. 8 |
| Dynamic viscosity | 0.0277 Pa·s | Table 1 |
| Density | 860 kg/m³ | Table 1 |
| Supply/reference pressure | 0.5 MPa | Table 1 |
| Speed for Figs. 6 and 8 | 2000 rpm | Figure-specific text/labels |

Table 1 also lists 2.6 m/s. At 50 mm radius, 2000 rpm gives 10.472 m/s;
these are not interchangeable operating inputs.

Fig. 6(b) is a prescribed-eccentricity numerical comparison. Its ordinate is
the **mid-plane maximum** pressure divided by 0.5 MPa, not the whole-film
maximum. Fig. 8 is a separate load-controlled experiment: 1000 N radial,
200 N axial, V-32 oil and five measured pressure positions.

The rig has a mid-plane oil feed IP1, five pressure ports and two temperature
ports. A 5 mm feed diameter, an experimental eccentricity of 0.6, and a freely
translating bush are **project assumptions**, not recovered rig specifications.
ISO VG 32 alone does not fix dynamic viscosity, density or operating temperature.

## Boundary conditions are not fully specified

The oil paper states circumferential periodicity and atmospheric axial sides,
but prints `p_bar=1` at both axial edges with `p_bar=p/p_s`.
That equation gives 0.5 MPa, not atmosphere. Its rupture conditions specify
zero pressure and zero circumferential pressure gradient at rupture, plus a
zero-pressure start line whose alignment for Fig. 6 is unclear. The numerical
cavitation pressure/reference and actual Fortran boundary implementation are
not recovered.

The companion [water-CFD paper](https://doi.org/10.5958/2321-581X.2018.00016.8)
explicitly imposes **0.5 MPa gauge at one side of the clearance** and atmosphere
at the opposite face. It does not identify which cone end is supplied.
[Publisher article](https://rjetonline.com/HTML_Papers/Research%20Journal%20of%20Engineering%20and%20Technology__PID__2018-9-1-16.html).
This establishes an end-fed water model, not the oil-CFD conditions or pressure
at the experimental IP1 mouth.

The oil paper's equations (14)–(15) balance transverse forces and update
transverse position. They do not specify a third axial-position equilibrium.
Axial seating matters: at 10°, 0.1 mm translation changes radial clearance
by about 17.6 µm.

Missing physical inputs include mouth/passage dimensions, pressure at IP1,
flow, operating oil temperature, axial restraint/seating, bush/shaft structural
dimensions, materials and support compliance. The paper is not a complete EHD
specification.

## Comparison implementation

`expts/gangrade-figure6/reference.toml` records source hashes, figure inputs,
pixel calibration and approximate readings from the native Fig. 6 raster.
The adopted ±1.5 normalized-pressure reading band (±0.75 MPa) is not
experimental uncertainty or an author-supplied confidence interval.

The experiment contains three distinct comparisons:

- `run.jl`: JFO with large-end, small-end or both-end pressure supply;
  ambient-end signed full film is a separate diagnostic.
- `compare_cavitation.jl`: pressure-only Reynolds/Swift–Stieber versus JFO
  with matched boundaries; pinned start-line alignments are separate sensitivities.
- `compare_q4.jl`: independently assembled Q4 Galerkin discretization on
  Fig. 2's 12 sectors × 4 bands (48 elements, 60 unique periodic nodes).
  Nodal gap interpolation, 2×2 quadrature, alignment and boundary interpretation
  are explicit reconstruction assumptions, not author code.

All use 2000 rpm and 0.0277 Pa·s. The JFO comparison assumes zero-gauge
cavitation and pressurized axial replenishment; it does not use the 5 mm port.

The complete published curve is **not matched**. At eccentricity 0.8, graph-read
FEA pressure is about 44.84 MPa. The reconstructed 60-node Q4 result is 50.19 MPa;
refined Q4/FV results are 67.828/67.843 MPa with both ends pressurized and a
pressure-only obstacle. Coarse-grid representation materially affects the
peak, but does not establish the paper's implementation. Refined FE/FV agreement
is numerical evidence, not physical validation.

## Checks and usage

From `journal-bearings/`:

```sh
julia --project=. --check-bounds=yes expts/gangrade-figure6/checks.jl
julia --project=. --check-bounds=yes expts/gangrade-figure6/reynolds_checks.jl
julia --project=. --check-bounds=yes expts/gangrade-figure6/q4_checks.jl
```

`run.jl` requires the local PDF and native raster at the paths recorded in
`reference.toml`, with matching hashes. Comparison scripts also check saved
input/source hashes; they are not standalone runs from geometry parameters:

- `compare_cavitation.jl` reads
  `out/cases/gangrade-figure6/runs/pressure-bc-hypotheses/`.
- `compare_q4.jl` reads the
  `cavitation-closure-comparison/comparison.toml` manifest under that case.

Each runner accepts a fresh output directory and refuses overwrite.
The corresponding `plot.py`, `plot_cavitation.py` and `plot_q4.py` read the
result tables; they require NumPy and Matplotlib. Source snapshots and numeric
tables under `out/` identify the exact inputs behind a saved result.
