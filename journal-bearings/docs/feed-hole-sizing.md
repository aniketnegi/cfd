# Feed-hole size and supply specification

A drilled oil-feed hole, a distribution groove, a pressure-measurement tapping
and an external pipe connection are different dimensions. A pipe thread or a
groove width does not determine the diameter of a standalone film opening.

## Source evidence

| Source | Reported feature | Applicability |
|---|---|---|
| [ISO 3547-3:2017, §6, Table 2](https://cdn.standards.iteh.ai/samples/70440/fc95a90b448d422594acf83989934cf6/ISO-3547-3-2017.pdf) | Minimum hole diameter after forming: 6 mm for inside diameter above 50 through 100 mm; 7 mm above 100 mm | A dimensional standard for wrapped bushes, not a conical or hydraulic-optimization rule |
| [GGB/Timken MBDX catalogue](https://cad.timken.com/viewitems/ggb-mbdx/ggb-mbdx-bushes---metric?pagenum=4&pagesize=25), e.g. MB10050DX and MB10095DX | 100 mm bore; 9.5 mm oil hole; lengths 50 and 95 mm | Product-specific cylindrical bush geometry, not a universal hydrodynamic or conical rule |
| [GGB AuGlide/SY/SP catalogue, §10.1](https://www.ggbearings.com/sites/default/files/2026-06/GGB-Auglide-SY-SP-Bimetal-Bearings-for-Lubricated-Applications-05.2026.pdf) | Tabulated oil-hole sizes for cylindrical SY bushes, including the 9.5 mm group containing 100 mm bores | ISO 3547 dimensions/tolerances and GGB product specifications, not a general conical inlet-sizing law |
| [KC Engineering A-series](https://kceng.com/product/plain-journal-bearings-a-series/) | Two feed holes on the split line, supplied by an external annulus | A different supply layout; the page does not specify a numerical hole diameter |
| [Singh (2015), §4](https://www.journalpressindia.com/download_article.php?aid=1815&jid=49), DOI [10.51976/ijari.321529](https://doi.org/10.51976/ijari.321529) | A 10 mm supply hole branches into two 5 mm holes feeding grooves; reported bearing diameter 70 mm and length 100 mm | A drilled-size precedent, not a standalone-hole or conical benchmark. The paper's accompanying L/D=1 label is inconsistent with those dimensions |
| [ISO 7902-1:2026](https://www.iso.org/standard/93204.html) and the accessible [2020 preview](https://cdn.standards.iteh.ai/samples/78987/a0a866ae9aea45b6a5b40de79d3349d4/ISO-7902-1-2020.pdf) | Steady circular cylindrical hydrodynamic-bearing calculation; distinguishes hole, groove and pocket dimensions | The current standard's scope is cylindrical. The 2020 preview does not expose its §7.3.2 feed-element discussion; no numerical sizing rule is inferred from unavailable clauses |

No cone-specific `d_feed/D` or `d_feed/L` sizing rule is established by these
sources. Catalogue examples are not sufficient to declare a universal standard.
Grease-hole recommendations, gas-bearing restrictors and hydrostatic multi-hole
arrays are not interchangeable with a single oil-fed hydrodynamic opening.

## Working cylinder/cone example

The working feed size is **7 mm**, explicitly supplied by the caller. This
standard-informed example uses ISO 3547-3's wrapped-cylinder dimensions by
analogy, not as conical certification or a recovered Gangrade dimension.
The table specifies minimum formed-hole sizes, not an optimum.

The mean journal diameter is `D_j=100 mm`; the housing bore at midlength is
`2*(50+0.05)=100.1 mm`, derived from the journal radius plus radial clearance.
It is not a separately selected input diameter. Geometric clearance and a
catalogue's nominal bush size must not be conflated: the derived 100.1 mm bore
does not by itself establish an ISO nominal-size classification. For a nominal
100 mm wrapped bush, the table lists a 6 mm minimum. The working 7 mm hole is
an explicit choice above that value, not a requirement inferred from clearance.
No ISO stock-bush nominal/tolerance specification is assigned to this geometry.

For axial length `L=100 mm`,

```text
d_feed / D_j = 0.07
d_feed / L   = 0.07
```

These ratios describe this example; the standard uses bore-size bands, not a
universal percentage. The same physical drill diameter is used for the zero-angle
cylinder and the 10° cone. `surface_feed_grid(...; feed_diameter=7., feed_z=50.)` represents the
rim of the +Y drilling at midlength. The journal offset is −Y, placing the feed
at maximum clearance. There is no distribution groove or resolved oil passage.
Public mesh/solver APIs do not select this size implicitly.

The example retains 2000 rpm, viscosity 0.0277 Pa·s, 0.5 MPa gauge at the feed
rim, zero-gauge axial ends and cavitation, and periodic circumference without
an additional angular pin. Geometry and pressure assumptions are separate:
choosing a diameter does not identify the rig's mouth pressure or flow rate.

## Pressure and liquid supply

For Swift–Stieber, the rim pressure and removed opening change the pressure
boundary problem, but no liquid-content field is transported.

For incompressible JFO, pressure above cavitation implies fully liquid supply.
Incoming Couette transport at the prescribed-pressure rim therefore uses liquid
content one; outgoing transport uses the interior upstream state. The reservoir
is bidirectional, not a check valve. The axial wall velocity has no normal
component at the straight axial ends, so there is no separate axial Couette
inflow-content datum in this reduced model. Pressure-driven end flows still
belong in the global balance. A prescribed rim pressure does not model pipe or
drilling resistance, and is not a prescribed oil-flow rate.

The maintained JFO implementation supports supplied tensor grids with explicit
end/cell reservoirs and fitted land meshes with actual rim transport. The
7 mm study uses the fitted rim, not a pressure-reservoir cell mask. See the
[solver support table](solver-support.md) and [JFO formulation](../simulation/reynolds/docs/conical-jfo.md).
