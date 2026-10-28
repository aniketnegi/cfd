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
Matplotlib.

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

## Independent end-flow and manufactured checks

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

Manufactured sources exist only in the test fixture `cases.jl`; physical-equation
cone runs invoke `solve_swift_stieber` itself. Mathematical zero-pressure end
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
consumption. No inlet size, supply patch, fitted-flow support, thermal physics
or elasticity is added by these verification tools.
