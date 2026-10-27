"""Static, reproducible mesh/field pictures; no new meshing implementation."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize

out = Path(sys.argv[1])
outputs = ('mesh-comparison.png', 'quality-convergence.png', 'node-placement.png',
           'pressure-fill.png', 'force-flow.png', 'index.html', 'plot.json')
for name in outputs:
    if (out / name).exists():
        raise SystemExit(f'refusing to overwrite plot evidence: {out / name}')
rows = list(csv.DictReader((out / 'comparison.tsv').open(), delimiter='\t'))
kinds = ('original', 'frozen', 'controlled')
labels = ('Harmonic transition', 'Fixed-map subdivision control', 'Controlled placement + fixed blocks')
colours = ('#b33b35', '#807b75', '#087e8b')

def patch(kind, level):
    a = np.genfromtxt(out / f'patch-{kind}-{level}.tsv', names=True, delimiter='\t')
    xy = np.stack([np.stack([a[f'x{i}_mm'], a[f'z{i}_mm']], axis=-1)
                   for i in range(1, 5)], axis=1)
    return a, xy

def draw(ax, kind, level, limits):
    a, xy = patch(kind, level)
    pc = PolyCollection(xy, array=a['corner'], cmap='viridis', norm=Normalize(.65, 1),
                        edgecolors='#2f3639', linewidths=.30, rasterized=True)
    ax.add_collection(pc)
    t = np.linspace(0, 2*np.pi, 200)
    ax.plot(2.5*np.cos(t), 2.5*np.sin(t), color='#dc267f', lw=1.4)
    points = np.genfromtxt(out / f'patch-{kind}-{level}-junctions.tsv', names=True, delimiter='\t')
    ax.scatter(points['x_mm'], points['z_mm'], facecolors='none', edgecolors='#dc267f', s=65, lw=1.4)
    ax.set(xlim=limits[0], ylim=limits[1], xlabel='X across mouth [mm]', ylabel='z − midpoint [mm]', aspect='equal')
    return pc

fig, axes = plt.subplots(2, 3, figsize=(16, 10), layout='constrained')
for c, (kind, title) in enumerate(zip(kinds, labels)):
    draw(axes[0, c], kind, 2, ((-15, 15), (-15, 15)))
    pc = draw(axes[1, c], kind, 4, ((-10, -3), (-10, -3)))
    axes[0, c].set_title(title + '\n33,664 cells over entire bearing', fontsize=12)
    row = next(r for r in rows if r['kind'] == kind and r['level'] == '4')
    axes[1, c].set_title(f"Junction close-up · 134,656 cells total\nmin corner {float(row['corner']):.3f}; max area jump {float(row['area_ratio']):.2f}:1", fontsize=11)
fig.colorbar(pc, ax=axes, shrink=.75, label='Local surface corner score (1 = right angle); NOT Hex8 OQ')
fig.suptitle('Same 5 mm inlet, same cell counts — different transition construction\n'
             'Front (+Y) surface only; magenta circles mark the four five-quad junctions', fontsize=15)
fig.savefig(out / 'mesh-comparison.png', dpi=170); plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
for kind, label, colour in zip(kinds, labels, colours):
    r = [r for r in rows if r['kind'] == kind]
    for ax, key, title in zip(axes, ('corner', 'alignment', 'area_ratio'),
                              ('Minimum corner score ↑', 'Minimum face alignment ↑', 'Maximum neighbour area ratio ↓')):
        ax.plot([int(x['cells']) for x in r], [float(x[key]) for x in r], 'o-', color=colour, label=label)
        ax.set(xscale='log', xlabel='Total surface cells', title=title); ax.grid(alpha=.25)
axes[0].set_ylim(.6, 1.01); axes[1].set_ylim(.8, 1.01)
axes[2].legend(fontsize=8)
fig.suptitle('Worst-cell / worst-face checks, not averages — all are surface metrics')
fig.savefig(out / 'quality-convergence.png', dpi=170); plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 8), layout='constrained')
a, xy = patch('controlled', 1)
ax.add_collection(PolyCollection(xy, facecolors='none', edgecolors='#65828b', linewidths=.6))
a = np.genfromtxt(out / 'node-movement.tsv', names=True, delimiter='\t')
ax.quiver(a['x_old_mm'], a['z_old_mm'], a['x_new_mm']-a['x_old_mm'], a['z_new_mm']-a['z_old_mm'],
          angles='xy', scale_units='xy', scale=1, color='#cc5538', width=.002)
t = np.linspace(0, 2*np.pi, 200)
ax.plot(2.5*np.cos(t), 2.5*np.sin(t), color='#dc267f', lw=2)
ax.set(xlim=(-23, 23), ylim=(-23, 23), xlabel='X across mouth [mm]', ylabel='z − midpoint [mm]', aspect='equal',
       title='Coarse blocks defined once, then subdivided\nArrows: mesh-node relocation, NOT physical wall displacement')
fig.savefig(out / 'node-placement.png', dpi=160); plt.close(fig)

# Read the generated legacy VTK, preserving cell-associated fields and periodic seam.
def vtk(path):
    lines = path.read_text().splitlines()
    p = next(i for i, s in enumerate(lines) if s.startswith('POINTS ')); n = int(lines[p].split()[1])
    xyz = np.array([list(map(float, s.split())) for s in lines[p+1:p+1+n]])
    c = next(i for i, s in enumerate(lines) if s.startswith('CELLS ')); nc = int(lines[c].split()[1])
    cells = np.array([list(map(int, s.split()[1:])) for s in lines[c+1:c+1+nc]])
    fields = {}
    for field in ('pressure_pa', 'fill_fraction'):
        j = next(i for i, s in enumerate(lines) if s.startswith(f'SCALARS {field} '))
        fields[field] = np.array(list(map(float, lines[j+2:j+2+nc])))
    return xyz, cells, fields
xyz, cells, fields = vtk(out / 'controlled-4.vtk')
theta = np.mod(np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])), 360)
angles = theta[cells]
seam = np.ptp(angles, axis=1) > 180
angles[seam] = np.where(angles[seam] < 180, angles[seam]+360, angles[seam])
polygons = np.stack((angles, xyz[cells, 2]), axis=-1)
fig, axes = plt.subplots(2, 1, figsize=(12, 7), layout='constrained')
for ax, key, scale, vmin, vmax, label in zip(axes, ('pressure_pa', 'fill_fraction'), (1e-6, 1),
                                           (0, .96), (.5, 1), ('Pressure [MPa gauge]', 'Liquid content fraction')):
    pc = PolyCollection(polygons, array=fields[key]*scale, cmap='viridis', clim=(vmin, vmax), rasterized=True)
    ax.add_collection(pc)
    ax.set(xlim=(0, 360), ylim=(0, 100), xlabel='Circumferential angle [degrees]', ylabel='Axial coordinate z [mm]')
    ax.annotate('5 mm reservoir mouth', xy=(90, 50), xytext=(115, 68), arrowprops={'arrowstyle': '->'}, color='#dc267f')
    fig.colorbar(pc, ax=ax, label=label)
fig.suptitle('Controlled transition, finest grid — SAME HELD POSE, not new equilibrium\n'
             'e/c = 0.0186497 · offset 0.9325 µm · clearance 50 µm · minimum normal gap 48.322 µm')
fig.savefig(out / 'pressure-fill.png', dpi=160); plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
for kind, label, colour in zip(kinds, labels, colours):
    r = [r for r in rows if r['kind'] == kind]
    for ax, key, title in zip(axes, ('Fz_N', 'radial_residual_N', 'net_feed_L_min'),
                             ('Axial land force [N]', 'Transverse residual [N]', 'Net feed [L/min]')):
        ax.plot([int(x['cells']) for x in r], [float(x[key]) for x in r], 'o-', label=label, color=colour)
        ax.set(xscale='log', xlabel='Total surface cells', title=title); ax.grid(alpha=.25)
axes[0].legend(fontsize=8)
fig.suptitle('Fixed pose: refinement diagnostic, not a 1000 N / 200 N equilibrium\n'
             'Axial target is 200 N, below the displayed axial-force range')
fig.savefig(out / 'force-flow.png', dpi=170); plt.close(fig)

html = '<!doctype html><meta charset="utf-8"><title>Inlet transition comparison</title><style>body{font:17px sans-serif;max-width:1450px;margin:2em auto}img{width:100%}</style>'
html += '<h1>Fitted inlet: controlled transition trial</h1><p>Same inlet and physical inputs. Surface metrics only—not Hex8 OQ, a new equilibrium, or EHD.</p>'
for name in outputs[:5]:
    html += f'<h2>{name[:-4].replace("-", " ").title()}</h2><a href="{name}"><img src="{name}"></a>'
(out / 'index.html').write_text(html)
inputs = [out / 'comparison.tsv', out / 'node-movement.tsv', out / 'controlled-4.vtk'] + list(out.glob('patch-*.tsv'))
(out / 'plot.json').write_text(json.dumps({'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}, indent=2)+'\n')
