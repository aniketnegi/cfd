# Cylindrical JFO verification

Compare finite, submerged cylinders with the published Jakobsson–Floberg tables,
and check the pressure/content calculation using continuous manufactured fields,
spatial refinement and an independent nodal finite-difference solution.
`reference.csv` contains **24 literal published comparison rows**, not solver
results; `reference.toml` defines their provenance and normalization.

## Run

From `journal-bearings/`, choose a fresh output directory:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. simulation/tests/runtests.jl cylindrical_jfo_checks.jl
julia --project=. simulation/verification/jfo/cylindrical/checks.jl
julia --project=. simulation/verification/jfo/cylindrical/run.jl out/verification/jfo-cylinder
# All 24 rows, both geometry models:
julia --project=. simulation/verification/jfo/cylindrical/run.jl --full out/verification/jfo-cylinder-full
python simulation/verification/jfo/cylindrical/analyse.py out/verification/jfo-cylinder-full
python simulation/verification/jfo/cylindrical/plot.py out/verification/jfo-cylinder-full
```

The default study uses `L/D=1, epsilon=.6` at both cavity pressures. Each case
runs at `64,128,256,512` circumferential cells, doubling the axial count exactly.
Both studies include manufactured fields and independent checks. Use `--help`
for options, `-v` for diagnostics or `-q` for errors only. No old output archive
or downloaded paper is required.

### Numerical models and checks

`cases.jl` supplies a classical cosine-gap/single-radius fixture to the shared
production `assemble_full_film` and `solve_conical_jfo` tensor operators. The
separately labelled `production` rows use the unchanged circle-intersection gap
and distinct journal/housing radii. There is no alternate production JFO solver.
Both models use the dimensional scale `R=50mm, c/R=.00291, mu=.0127Pa*s,
omega=48.1rad/s`; these scale the dimensionless tables, not 24 measured rigs.
The computational seam is at minimum clearance; it is periodic, never pinned.
Plot angles wrap to the source's `[-pi,pi)` convention.

`manufactured.jl` differentiates a continuous pressure/content pair with curved
interfaces. Its artificial distributed source is integrated by 3-by-3 Gauss
quadrature and passed only to the low-level complementarity kernel. Uniform
and graded tensor coordinates and both rotation signs are checked. Pressure
error uses relative area-weighted L2; content error uses area-mean absolute
error. No discrete `A*exact` forcing or prescribed active mask is used.

`independent.jl` includes no production source. It uses nodal expanded pressure
finite differences and a one-sided derivative of `H*(1-fill)`, with an
independent scalar projected Gauss–Seidel iteration in `X=P-(1-fill)`.
It checks the representative cases through `256` angular nodes; `--full` also
checks `epsilon=.4,.8, pv0=-.1` at `L/D=1`. Its end
flow uses second-order one-sided pressure derivatives. This expanded nodal
scheme is not exactly conservative at finite resolution: its signed net end
flow is reported separately, not held to the conservative FV roundoff bound.
Stopping/relaxation controls and a same-FD-matrix production-policy comparison
check algebra independently of the spatial comparison.

All comparisons interpolate at actual coordinates. Self-differences are not
continuum error bounds. Analysis checks physical-state bounds, complementarity,
local/global conservation, manufactured refinement and independent differences;
it assigns **no acceptance threshold to the historical printed values**.

### Outputs

- `inputs.toml`: dimensional scale, resolution and source hashes.
- `results.tsv`: load, attitude, gross inflow/outflow, signed net flow, bounds,
  conservation, reference values and successive-grid field changes.
- `profiles.tsv`: finest-grid pressure/content interpolated at the actual
  midplane, not a substituted axial row.
- `manufactured.tsv`, `independent.tsv`, `algebraic.tsv`: separate continuous,
  spatial and algebraic checks.
- Representative finest fields: JLS, VTK and classical-field TSV files.
  Classical VTK geometry is a nominal reference cylinder, not an eccentric solid.
- `analysis.json`, `verification.png`, `reference.png`, `profiles.png`,
  `fields.png`: derived checks and plots.

Historical comparison, numerical verification and physical validation are
separate conclusions. Pressure/front curves from the sources below have not
been digitized; no experimental field validation is claimed.

## Modern pressure comparison

[Xi2024](xi2024/README.md) provides a separate finite-cylinder pressure-curve
comparison with four prescribed-pressure edges, boundary-adapter checks and
an independent nodal calculation. Its gap-origin convention is an explicit
assumption; its graph-derived values are not author raw arrays. Do not mix
that supply boundary problem with the ungrooved1957 case below.

## Primary reference: Jakobsson–Floberg, 1957

*The Finite Journal Bearing, Considering Vaporization*, Transactions of
Chalmers University of Technology190, is reproduced in Floberg's1961 thesis:
[official record](https://research.chalmers.se/en/publication/177651),
[full text](https://research.chalmers.se/publication/177651/file/177651_Fulltext.pdf).

Use **§5.2**, the finite, ungrooved, submerged cylindrical bearing, with
mass-conserving rupture and reformation. This is the steady Jakobsson–Floberg
part of JFO theory. It is NOT §7's atmospheric-cavity, maximum-gap-grooved
Cameron–Wood comparison used for Swift–Stieber.

The appendix gives load, attitude and oil-flow values for

- `L/D = 1, 0.5, 0.25`;
- eccentricity ratios `0.2, 0.4, 0.6, 0.8`;
- dimensionless cavity pressure `pv0 = -1, -0.1`.

These are the finite positive-eccentricity columns of printed pp106–109
(PDF pp227–230). Some rows are noncavitating limits; the whole set must not be
labelled as24cavitating fields. Zero-eccentricity and singular epsilon1columns
are excluded. The four table pages were visually checked when transcribing.
Printed digits are not a bound on the historical discretization/iteration error.
No numerical acceptance threshold or claim of continuum truth accompanies them.

### Matched model and boundaries

The source assumes a rigid, aligned, isoviscous, incompressible thin film:

```text
H(phi) = h/c = 1 - epsilon*cos(phi)
U = omega*R
```

It uses one nominal radius. Minimum clearance is at `phi=0`; the maximum is
at `phi=-pi,+pi`. Circumference is periodic, with **no imposed angular pressure
line**. Both axial ends are immersed in oil at ambient pressure. The cavity
pressure is below ambient, permitting inward oil flow at the sides as well as
outward flow elsewhere. There is no hole or groove.

Thus the project's provisional drilled-inlet diameter and0.5MPa mouth pressure
are **not inputs to this reference case**. Ambient ends are strictly above the
cavity pressure here; this is different from an unspecified-inventory periodic
film with every reservoir at the cavitation pressure.

The production circle-intersection gap and separate journal/housing radii are
not exactly the source's classical cosine-gap/single-radius model. A matched
classical fixture and a separately labelled production-geometry comparison
must not be silently interchanged.

### Observables and pressure reference

With full axial width `L`, radial clearance `c`, journal radius `R`, dynamic
viscosity `mu`, and angular speed `omega`:

```text
pressure_scale = mu*omega*R^2/c^2
p0  = (p-p_ambient)/pressure_scale
pv0 = (p_cav-p_ambient)/pressure_scale
P0  = W*c^2/(mu*omega*R^3*L)
Q0  = Q_out/(omega*R^2*c)
```

`W` is the total resultant load; the source also uses load per unit width.
`beta` is the angle between the load line and line of centres.

**Q_out is gross outward exchange through both ends**, as defined in §5.23,
not signed net end flow. In a steady submerged bearing without an internal
source, gross outflow equals gross inflow; signed net flow is zero. Report all
three quantities and preserve inward face flows. Comparing the source's
positiveQ0with signed net flow would be comparing different observables.

Dimensional scales must be stated explicitly. A nondimensional negative cavity
pressure is relative to surrounding oil; it is not a negative absolute pressure.
The source experiment's `R=50mm,c/R=.00291,mu=.0127Pa*s,omega=48.1rad/s`
is a possible scale. The nondimensional tables do not imply that every row was
an experiment at these settings.

### Pressure profiles and experimental limits

§5.21 provides calculated isobar/cavity contours in Figs55.2 and56.1–56.3.
§5.27 provides pressure measurements for two finite bearings and their
corresponding theoretical curves. These offer later field/front comparisons,
but the curves have not been digitized into this fixture.

For the first experiment, the target was epsilon.6, but the source says the
2250Nload corresponds to approximately.61. Fig69.1's traces are at axial
positions `(z/R)/(L/D)=3/5` and`1/5`, **not exactly the midplane**. The figure
also labels a vacuum-pump connection. These details must accompany any physical
comparison rather than treating the plotted points as exactepsilon.6midplane data.

## Secondary and cross-check references

### Xu, Zhao, Xu and Li, 2021

[*Reynolds Model versus JFO Theory in Steadily Loaded Journal Bearings*](https://doi.org/10.3390/lubricants9110111)
is open access and provides finite-cylinder pressure profiles, cavity contours,
load/attitude and flow trends. Fig2 gives a fully dimensioned historical case;
Figs4–9 vary cavity pressure, eccentricity andL/D.

Its §3parametric study supplies oil on the maximum-gap line and imposes
ambient pressure on all four unwrapped edges. That is not the same boundary
problem as the primary reference's ungrooved, circumferentially periodic bearing.
The boundary treatment of the separate Fig2validation must be checked rather
than inferred from the later parametric study.

Their Elrod formulation uses a **compressible** pressure–density relation,
`p=pc+beta*ln(theta)`; Fig2 specifies`beta=1e10Pa`. Our incompressible liquid
fraction is not the same variable as their density ratio above one. Pressure
references and the definition of bearing number must also be mapped carefully.
The article's availability statement says the data are in the article; no
separate numerical field dataset was recovered. Fig2 reuses older measurements,
so it is not an independent new experiment. Its “central” profile label needs
mapping to the original source's axial measurement positions.

### Dowson, Miranda and Taylor, 1985

[*Prediction of liquid film journal bearing performance with a consideration of
lubricant film reformation, Part2: experimental results*](https://doi.org/10.1243/PIME_PROC_1985_199_100_02)
has an accessible [university-repository PDF](https://repositorium.uminho.pt/bitstreams/3b08e7db-15ea-4ccf-bf63-81b501c1108b/download).
It reports measured side flow and rupture/reformation boundaries in finite
cylinders with a single hole or finite axial grooves. These are valuable
mass-transport observables beyond peak pressure. Each bush's supply geometry,
pressure and effective-viscosity convention must be matched. Raw numerical
arrays were not recovered; plot digitization is required. The theoretical
Part1 (DOI10.1243/PIME_PROC_1985_199_099_02) remains a lead: its full text was
not accessible in this review.

### Alakhramsing, van Ostayen and Eling, 2015

[DOI10.3390/lubricants3020256](https://www.mdpi.com/2075-4442/3/2/256), §2.2/Fig8,
compares **finite-cylinder** pressure and content profiles with Yu–Keith1994.
This is distinct from its Fig5infinite-width case. Its artificial diffusion,
axial domain/normalization and fully submerged supply must be audited before
using the plot as a matched target. Original Yu–Keith1994
(DOI10.1080/10402009408983286) full text was unavailable here. The later
Ferron-based thermal experiment in the2015paper is not an isoviscous benchmark.

Dynamic NASA Woods–Brewe cases, the Ausas2009oscillating/dynamic cases,
piezoviscous examples and rotor-vibration datasets are not interchangeable
with the steady, rigid, finite-cylinder reference above.
