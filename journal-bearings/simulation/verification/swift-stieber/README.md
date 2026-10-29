# Swift–Stieber verification

Check the [pressure-only solver](../../reynolds/docs/swift-stieber.md) against
independent finite differences, continuum limits and manufactured solutions.
Run from `journal-bearings/`:

```sh
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
julia --project=. --check-bounds=yes simulation/tests/swift_stieber_checks.jl
julia --project=. --check-bounds=yes simulation/verification/swift-stieber/sampling_checks.jl
julia --project=. simulation/verification/swift-stieber/run_cylindrical.jl out/verification/ss-cylinder
julia --project=. simulation/verification/swift-stieber/run_conical.jl out/verification/ss-cone
julia --project=. simulation/verification/swift-stieber/compare_fd.jl out/verification/ss-independent
python simulation/verification/swift-stieber/analyse.py out/verification/ss-cone
python simulation/verification/swift-stieber/plot.py out/verification/ss-cone
```

Use fresh output directories. Runs save inputs, source hashes, comparison tables
and fields there. Python analysis uses the standard library; plotting requires
Matplotlib. Use `--help` for options, `-v` for case diagnostics or `-q`
for errors only.

Analyse and plot cylinder, cone and independent FD comparisons with the same
commands, passing the corresponding run directory. Analysis also accepts the
specialist check tables below; plotting supports cylinder, cone, MMS,
fitted-refinement and FD-comparison tables. Use `plot.py DIRECTORY --table results.tsv` if a directory
contains several supported tables. Postprocessing replaces its derived JSON/PNG.

## Cylinder and independent finite-width calculation

`run_cylindrical.jl` compares the classical cosine gap/single radius with the
unchanged production circle-intersection gap/separate radii. It prescribes
an aligned zero-pressure maximum-gap line and zero axial ends, not an
unspecified submerged or locally fed bearing. Numerical scales are Rj50mm,
c50micrometres,2000rpm,mu0.02752Pa*s. `reference.csv` and `reference.toml`
identify Jakobsson-Floberg1957, reproduced in Floberg1961, its table and
normalization. Printed targets remain literal reference values, not exact
solver acceptance criteria.

The representative case is L/D1,epsilon.8 on64/128/256 angular cells; both
models run. Add `--full` for all11printed finite cases at64/128/256/512.
Each family's axial base doubles exactly with its angular count. The two
geometry models are labelled and are never switched to improve agreement.

`independent_fd.jl` independently constructs the expanded nodal pressure
operator and uses projected SOR, physical boundary nodes and one-sided end
derivatives. It includes NO production geometry/assembly/solver/flux routine.
`compare_fd.jl` loads both implementations to compare the SAME classical
finite-width model. Its default grids are64/128/256; `--full` extends through
1024x512 and adds finer stopping/relaxation controls. It saves full pressure
norms, separate end flows, end-gradient/flux-density norms, rupture markers
and physical-coordinate profiles. No flow-density clipping is applied.

The source defines Q0 as TOTAL two-end leakage/(omega*R^2*c). Its classical
L/D1,epsilon.8 printed Q0 is1.30; independently converged same-model results
are near1.240477. The historical comparison remains unreconciled, not an
identified FV bug or proof the publication is incorrect. Original
Cameron-Wood numerical details are not recovered. Neither an unrelated coarse
mesh description nor infinite-width pressure agreement establishes its
finite-width flow accuracy.

## Unported element-count comparison

```sh
julia --project=. simulation/verification/swift-stieber/q4_checks.jl
julia --project=. simulation/verification/swift-stieber/unported_grid.jl out/verification/ss-unported-grid
python simulation/verification/swift-stieber/analyse.py out/verification/ss-unported-grid
python simulation/verification/swift-stieber/plot.py out/verification/ss-unported-grid
```

`unported_grid.jl` compares the shared FV solver with the verification-only
nodal Galerkin assembly in `q4.jl`. Both consume the same public `housing_grid`
nodes and quads. The12x4 topology has48elements,60periodic Q4 nodes and48FV
cell-centred pressures; the pressure locations are not interchangeable.
Uniform theta/z axes refine through24x8,48x16,96x32,192x64 and384x128.
Eccentricities.6/.8 also get a512x256FV comparison. Add `--full` for all eight
Fig.6(b) eccentricities. The default runs.6/.8.

The geometry, production gap/radii/wall-speed laws,2000rpm andmu0.0277Pa*s
stay fixed. Both ends and the cavity have zero gauge pressure; circumference
is periodic, with no opening or additional angular pin. Journal offset is−Y
(`eccentricity_angle_rad=3pi/2`); the first mesh ray is+X. The0.5MPa reference
normalizes pressure only. These are explicit unported comparison conditions,
not a recovered implementation of the paper's ambiguous axial boundaries.

Q4 bilinearly interpolates the production nodal gaps and uses2x2Gauss
quadrature on the analytical cone metric. It reuses the algebraic pressure
obstacle, not the FV assembly. Actual author coordinates, quadrature and
iteration details are not available. `q4_checks.jl` checks topology, energy,
Couette forcing/sign, scaling, rotation, exact stationary solutions and an
independently differentiated continuous manufactured field.

Tables report actual midplane pressures, nodal/cell counts, residuals and
natural pressure along alpha=0. Fine midplane profiles are compared at common
angles. Separate fine-field sampling at the12coarse nodal angles distinguishes
peak-sampling loss from the coarse solve. Nodal pressures and Fig.2 numbering
maps are exported for.6/.8; Q4 VTK uses POINT_DATA, not fictitious cell pressures.
The small `expts/gangrade-figure6/reference.toml` supplies approximate graph
readings without PDFs or old run directories. Its±0.75MPa visual-reading band
is not experimental uncertainty. A coarse match alone does not establish model
identity or quantify an inlet diameter.

## Fitted pressure-only cavitation

```sh
julia --project=. simulation/verification/swift-stieber/fitted_verification.jl out/verification/ss-fitted
python simulation/verification/swift-stieber/analyse.py out/verification/ss-fitted
python simulation/verification/swift-stieber/plot.py out/verification/ss-fitted --table fitted-mms.tsv
```

`fitted_verification.jl` checks cylinders and10degree cones at128x64 and256x128;
`--full` adds512x256. Known half-film pressure fields use independently derived
forcing, both rotation signs, and fronts away from and crossing the inlet rim.
Analytical end flows, independently integrated rim flows and an analytic
whole-cone force minus the opening integral check boundary flux and force.
The active mask is not supplied. `fitted_reference.jl` also provides projected
SOR for a separate nonlinear-solve comparison on the SAME assembled matrix;
this is an algebraic check, not a second spatial discretization.

Physical-equation refinement uses fixed position, eccentricities.2/.6/.8,
+/-2000rpm,mu0.0277Pa*s,zero-pressure ends and an explicitly prescribed
0.3MPa rim on a5mm midlength opening. These are numerical test conditions,
not adopted Gangrade boundaries. `fitted-refinement.tsv` reports global peaks,
forces and separate boundary exchanges; it is not a Fig.6 midplane comparison.
Use `--table fitted-refinement.tsv` to plot its global peaks. Observable changes
are self-convergence measurements, not absolute error bounds.

## Specialist studies

### Independent end-flow and manufactured checks

```sh
julia --project=. simulation/verification/swift-stieber/fd_checks.jl out/verification/ss-fd-checks
julia --project=. simulation/verification/swift-stieber/fv_flux_checks.jl out/verification/ss-flux-checks
julia --project=. simulation/verification/swift-stieber/long_bearing.jl out/verification/ss-long
julia --project=. simulation/verification/swift-stieber/conical_mms.jl out/verification/ss-mms
julia --project=. simulation/verification/swift-stieber/conical_checks.jl out/verification/ss-controls
```

- `fd_checks.jl`: analytical nodal end derivatives, finite-width manufactured
  PDE/flow and stopping/relaxation controls, without loading the FV code.
- `fv_flux_checks.jl`: independently forced linear/quadratic pressure solutions
  check actual axial conductances, half-cell boundary distances, signs, H^3
  weighting and both-end accounting. Sampling a continuum quadratic and solving
  its FV equations have different cell-centre offsets; both are checked.
- `long_bearing.jl`: independently integrated rupture ODE and infinite-width
  pressure refinement. This does not verify axial leakage. The source's last
  beta entry is retained without claiming it matches the integrated attitude.
- `conical_mms.jl`: known half-film sin^4 pressure, independently differentiated
  eccentric rotating cone coefficients/source, positive cavity reaction,
  exact force and separate-end-flow references. Angles0/10/20, both rotation
  signs and uniform/graded meshes are exercised. The active mask is not supplied
  and the manufactured forcing is NOT assembled as A applied to a desired field.
- `conical_checks.jl`: independent high-precision gaps, cone metric/face fluxes,
  Cartesian surface-normal force quadrature, analytic force convergence,
  continuous zero-angle recovery, rotations/scaling and no-clipping controls.
  Its cylinder checkpoint is a regression fixture, not an independent target.

Manufactured sources are verification fixtures; physical-equation runs invoke
`solve_swift_stieber` without artificial forcing. Mathematical zero-pressure end
conditions are not adopted as Gangrade's inlet conditions.

## Conical refinement and sampling

`run_conical.jl` directly computes every requested level. Its default is the
10degree,epsilon.6 cone with periodic or separately pinned maximum-gap line,
64x32 ->128x64 ->256x128. `--full` runs10degrees withepsilon.2/.6/.8 for both
boundaries, plus20degrees epsilon.6 periodic, through512x256. Conical scales
use mu0.0277Pa*s; input manifests label both cylinder and cone viscosities.

Midplane pressures interpolate actual axial centres at L/2. Whole-field
refinement samples fine pressure at actual coarse pressure-centre coordinates
with periodic bilinear interpolation and physical-area weighting. Consecutive
fine array rows are not assumed to cover a shifted coarse angular cell.
`sampling_checks.jl` exercises this distinction, periodic wrapping and graded
coordinates. Maximum pressure over the entire cone is reported SEPARATELY
from midplane maximum. Self-convergence differences are not absolute error bounds;
no universal second-order Richardson extrapolation is applied to moving peaks.

The force is on the housing sidewall, excluding end faces/shear, not a journal
force or prescribed-load equilibrium. Fully-liquid cavity flux defects are
reaction bookkeeping, NOT JFO liquid-content transport or validated oil
consumption. The fitted study uses the maintained prescribed-rim operator;
it does not model passage flow, thermal physics or elasticity.
