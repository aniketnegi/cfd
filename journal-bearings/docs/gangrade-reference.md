# Gangrade conical-bearing reference

## Published inputs

Gangrade, Phalle and Mantha, *Performance Analysis of a Conical Hydrodynamic
Journal Bearing*, 43:559–573 (2019):
[DOI 10.1007/s40997-018-0217-2](https://doi.org/10.1007/s40997-018-0217-2).

| Quantity                  |       Value | Source                      |
|---------------------------|------------:|-----------------------------|
| Mean journal radius       |       50 mm | Table 1                     |
| Axial length              |      100 mm | Table 1                     |
| Radial clearance          |       50 µm | Table 1                     |
| Semi-cone angle           |         10° | Fig. 6 comparison; Fig. 8   |
| Dynamic viscosity         | 0.0277 Pa·s | Table 1                     |
| Density                   |   860 kg/m³ | Table 1                     |
| Supply/reference pressure |     0.5 MPa | Table 1                     |
| Speed for Figs. 6 and 8   |    2000 rpm | Figure-specific text/labels |

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

## Reference data and maintained solvers

`expts/gangrade-figure6/reference.toml` records the figure inputs, primary
source fingerprints, pixel calibration and approximate Fig.6 readings.
The adopted ±1.5 normalized-pressure reading band (±0.75 MPa) is a visual
reading convention, not experimental uncertainty or a confidence interval.
A local PDF/raster is optional evidence for auditing those readings, not a
dependency of the numerical solver or its verification tests.

Use the maintained [Swift–Stieber API](../simulation/reynolds/docs/swift-stieber.md)
for pressure-only calculations and the [verification tools](../simulation/verification/swift-stieber/README.md)
for independently checked cylindrical/conical test problems. The small
`reynolds_cavitation.jl` count adapter delegates to that same shared solver;
`reynolds_checks.jl` exercises it without loading a Gangrade comparison run.

No maintained runner reconstructs the experimental IP1 passage or resolves
the conflicting axial/start-line descriptions. End reservoirs, a pinned
zero-pressure line and a fitted circular feed are distinct boundary models.
The benchmark's maximum-clearance groove is not an implicit Gangrade input.
A diameter supplied by the caller is a study parameter, not a recovered
published dimension.

A future parameter study must state its selected pressure reference, feed
location and boundary interpretation. An exploratory fitted parameter is not
independent validation. The [solver support table](solver-support.md) still
marks fitted fluid simulation unsupported; accepted mesh geometry alone does
not verify a nonorthogonal flow operator.
