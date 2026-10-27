# Mass-conserving cavitation: 1D slider benchmark

## Scope

`jfo_slider.jl` solves steady rigid/isoviscous/incompressible finite-volume
pressure and liquid content. Its tests compare against an independently
evaluated continuum solution for a published sinusoidal slider, including
rupture **and reformation**.

**This is numerical verification, not experimental validation or conical EHD.**
The benchmark's lubricant, pressure, geometry and supply are not choices for
the Gangrade/VG32 bearing. Its open streamwise inlet replenishes the film;
it does not resolve the periodic bearing's supply/inventory ambiguity.

## Audited primary reference

J. Lengiewicz, M. Wichrowski and S. Stupkiewicz, *Mixed formulation and finite
element treatment of the mass-conserving cavitation model*, Tribology
International **72** (2014), 143–155,
[DOI 10.1016/j.triboint.2013.12.012](https://doi.org/10.1016/j.triboint.2013.12.012).

- [Author-hosted accepted manuscript](https://bluebox.ippt.pan.pl/~sstupkie/papers/Lengiewicz_etal_TribInt_2014_accepted.pdf):
  §2 defines entrainment, mass continuity and JFO interfaces; §4.1 gives this
  benchmark; Figs. 2–3 show the analytical solution and FE refinement.
- [Author supplement landing page](https://bluebox.ippt.pan.pl/~sstupkie/files/bearing.html),
  linked from the authors' publication list and the manuscript footnote.
- [Analytical notebook PDF](https://bluebox.ippt.pan.pl/~sstupkie/files/bearing.pdf),
  version 1.3, 2014-02-01: equations, inputs, **printed numerical front values**
  and symbolic solution. This is the quantitative reference, not graph reading.

**Length correction:** the supplement, p. 2, explicitly corrects the paper's
12.5 mm to **125 mm**. The retrieved accepted manuscript already prints
125 mm. Both sources were checked, including rendered parameter/front pages.
Do not silently use 12.5 mm or infer a length from a plot.

Published specification, using the paper's centred x coordinate:

| Input | Value |
| --- | --- |
| Domain | −62.5 ≤ x ≤ 62.5 mm |
| Gap | h(x) = 20 µm − 5 µm cos(2πx/125 mm) |
| Dynamic viscosity | 0.015 Pa·s |
| Entrainment speed u | +2 m/s, **average** of the two wall speeds |
| End pressures | 1 MPa at both ends, relative to the zero cavitation level |
| Cavitation level | p = 0 in the source's pressure coordinate |
| Incoming content | fully liquid: φ=1, implied by p_in > p_cav |

The notebook uses `(4000/2)` mm/s for entrainment. For a stationary wall that
means a moving-wall speed of 4 m/s, not 2 m/s. No density value is necessary:
constant liquid density multiplies every conserved liquid-equivalent flux.
The zero pressure coordinate is not a claim that a physical oil has zero
absolute vapour pressure, nor a selected gauge cavitation pressure for VG32.

Local source snapshots, extracts and checked images:
`out/literature/jfo-slider/`. The evidence manifest records PDF SHA-256 hashes.
These cached sources and run outputs are Git-ignored; source URLs remain here.

## Chosen discretization, distinct from published evidence

Use φ for liquid content (relative density, or fill fraction when gas mass is
neglected), v=1−φ for void, and P=p−p_cav for pressure above cavitation:

```text
Q = u h φ − h³/(12 μ) dP/dx,
dQ/dx = 0,
P ≥ 0,  0 ≤ φ ≤ 1,  P(1−φ) = 0.
```

Q is liquid-equivalent volume flow per unit width (m²/s); mass flow per unit
width is ρ_liquid Q. At rupture the full-film pressure gradient vanishes. At
reformation the pressure gradient and content can jump, but Q is continuous.
This is not a pressure clamp or a pressure-only obstacle model.

Implementation choices (not the paper's finite-element algorithm):

- Uniform cell-centred 1D grid; h evaluated directly at centres and faces.
- One shared face flux, inserted with opposite signs into adjacent cells:
  `Q_face = u h_face φ_upwind − h_face³/(12 μ) (P_right−P_left)/distance`.
- Face-pressure differences are centred; transported content is first-order
  upwind. Pressure boundaries use half-cell distances. The mixture/content
  factor is not applied again to the pressure mobility: in cavitated interiors
  P is identically zero, whereas in full film φ=1.
- Advective inflow content must be explicit. A pressurized inlet requires
  φ=1; underfilled inflow is allowed only at p_cav. Outgoing content comes
  from the last cell, not a forced φ=1 at the outflow.
- Solve `A P + B v = b` with active-set column selection: P unknown and v=0
  in full cells; P=0 and v unknown in cavitated cells. Scale pressure for the
  classification `P/pressure_scale − v > 0`. Julia's sparse linear algebra
  solves each system. All-full initialization; bounded iteration count.
- Do not clip pressure/content after the solve. Failure to converge or failure
  of bounds, complementarity, cell balance, boundary balance or common flux
  checks throws. No claim of global convergence for arbitrary gap functions.
- Nonzero entrainment is required here. A zero-speed underfilled state needs
  a separate inventory contract and is not silently selected by this solver.

The accepted manuscript also studies nonlinear and linear FE weak forms.
We are comparing the shared continuum benchmark, **not** reproducing their FE
algorithm, iteration counts, or formal convergence proof.

## Independent reference and numerical evidence

`jfo_slider_reference.jl` does not use the discrete operator. Starting with
constant flux Q=u h*, it integrates

```text
dp/dx = 12 μ u (h−h*)/h³.
```

At rupture x_c, h*=h(x_c) and p(x_c)=0. Bracketed bisection and 64-point
Gauss–Legendre integration solve the left-reservoir condition. The same h*
then determines reformation x_r by integrating from the right reservoir.
Inside the cavity, φ=h*/h. Outside it, φ=1. This avoids the branch-cut problem
of a trigonometric symbolic antiderivative noted by the authors.

The reference is checked against the supplement's rounded outputs:
`x_c=20.8491 mm`, `x_r=56.7045 mm`, `h*=0.0175034 mm`,
`p(0)=3.31208 MPa`, with tolerances set by printed precision. Increasing
quadrature from 64 to 96 points changes pressure at 201 sample locations by
less than 1e-5 Pa. Thus this is an independent evaluation of the published
solution—not an independently authored second implementation or human review.

Computed continuum values:

- Rupture: **20.849116671 mm**.
- Reformation: **56.704507926 mm**.
- Peak: **6.624157850 MPa**, at x=−20.849116671 mm.
- Common Q: **3.500687225e-5 m²/s**.

Six-grid result, Julia 1.13.2:

| Cells | Pressure relative L2 error | Fill mean absolute error |
| ---: | ---: | ---: |
| 40 | 9.98061e-4 | 6.13292e-3 |
| 80 | 2.49629e-4 | 2.56974e-3 |
| 160 | 6.14542e-5 | 2.60477e-3 |
| 320 | 1.59469e-5 | 7.66791e-4 |
| 640 | 3.76916e-6 | 3.08353e-4 |
| 1280 | 1.04467e-6 | 3.09096e-4 |

Errors compare values at cell centres. Pressure is approximately second order.
**Content error is not monotonic:** the reformation jump moves relative to
cell centres, as the paper also discusses for its scheme. The 1280-cell
content error is slightly larger than at 640; do not hide it or claim smooth
second-order content convergence. The discontinuity dominates these norms.

At 1280 cells (dx=0.09765625 mm):

- Discrete peak: 6.624165301 MPa; approximately 1.1e-6 relative peak error.
- Classified rupture/reformation faces: 20.80078125 / 56.73828125 mm;
  errors 0.04834 / 0.03377 mm, both less than one cell.
- Relative flux error against the continuum: 6.25e-8.
- Maximum cell imbalance: 3.89e-14; boundary imbalance: 6.95e-13;
  whole-domain flux spread: 6.99e-13, each normalized by max|u h|=5e-5 m²/s.
- Seven active-set iterations; no pressure/content clipping.

In addition to the source outputs and six grids, tests cover reverse motion, pressure-reference shifts, speed/viscosity scaling,
constant-gap pressure-driven flow, explicit starved/dry inflow, rejection of
inconsistent inputs and nonconvergence, and a pressure-clamping counterexample.
The §4.2 fluctuating-gap example exercises two cavities and conservation only;
it is not another quantitative full-field benchmark.

Regression guards include pressure refinement ratio <0.4, finest pressure L2
<1e-4, finest fill mean absolute error <1e-3, fronts within two cells, and
continuity residuals <1e-10. These are numerical test guards, not physical
bearing accuracy criteria or error estimates for EHD.

## Reproduce and inspect

From the repository root:

```sh
julia --project=. --check-bounds=yes simulation/reynolds/jfo_slider_checks.jl
julia --project=. --check-bounds=yes simulation/reynolds/jfo_slider_verification.jl
python expts/jfo-slider-verification/plot.py
```

The last two commands refuse to overwrite existing evidence. To repeat, pass
a **fresh output directory** as the sole argument to the verification script,
then that same directory to the plot script. Tests do not require downloaded
PDFs, Python or a plotting package; the optional plot uses NumPy/Matplotlib.

Saved evidence:
`out/cases/lengiewicz-sinusoidal-slider/runs/jfo-verification/`:
`verification.toml`, `convergence.tsv`, `profiles.tsv`, `face-flux.tsv`,
`verification.png`, `checks.log`, `conical-regression.log`,
`cylindrical-regression.log`. The manifest records solver, reference, test,
evidence-script, plot-script and available source-PDF hashes.

## Bearing solver

The [2D periodic conical solver](conical-jfo.md) shares `jfo_active_set.jl`
but verifies its own flux assembly with a curved manufactured cavity,
full-film limits and reservoir accounting. The 1D benchmark alone does not
verify conical metric factors or feed accounting.
