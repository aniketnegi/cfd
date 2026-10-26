"""Inspect the actual public POLYDATA export (land polygons plus rim lines)."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection

mesh = Path(sys.argv[1])
output = Path(sys.argv[2])
manifest = output.with_suffix('.json')
if output.exists() or manifest.exists():
    raise SystemExit(f'refusing to overwrite preview: {output}')
lines = mesh.read_text().splitlines()
i = next(i for i, s in enumerate(lines) if s.startswith('POINTS '))
n = int(lines[i].split()[1])
xyz = np.array([list(map(float, s.split())) for s in lines[i+1:i+1+n]])
i = next(i for i, s in enumerate(lines) if s.startswith('LINES '))
nrim = int(lines[i].split()[1])
rim = np.array([list(map(int, s.split()[1:])) for s in lines[i+1:i+1+nrim]])
i = next(i for i, s in enumerate(lines) if s.startswith('POLYGONS '))
nland = int(lines[i].split()[1])
quads = np.array([list(map(int, s.split()[1:])) for s in lines[i+1:i+1+nland]])
i = next(i for i, s in enumerate(lines) if s.startswith('SCALARS surface_corner_score '))
quality = np.array(list(map(float, lines[i+2:i+2+nrim+nland])))[nrim:]
assert np.isfinite(quality).all()
# Never project the opposite (-Y) side onto the inlet side.
front = xyz[quads, 1].mean(axis=1) > 0
polygons = xyz[quads][:, :, (0, 2)].copy()
polygons[:, :, 1] -= 50
boundary = xyz[rim][:, :, (0, 2)].copy()
boundary[:, :, 1] -= 50
fig, axes = plt.subplots(1, 2, figsize=(13, 6), layout='constrained')
for ax, limits, title in zip(axes, ((-15, 15), (-10, -3)),
                             ('Fitted 5 mm inlet', 'Transition-junction close-up')):
    pc = PolyCollection(polygons[front], array=quality[front], cmap='viridis', clim=(.65, 1),
                        edgecolors='#354348', linewidths=.25, antialiased=False, rasterized=True)
    ax.add_collection(pc)
    ax.add_collection(LineCollection(boundary, colors='#dc267f', linewidths=1.4))
    ax.set(xlim=limits, ylim=limits, aspect='equal', title=title,
           xlabel='X across mouth [mm]', ylabel='z − midpoint [mm]')
fig.colorbar(pc, ax=axes, label='Surface corner score — NOT Hex8 OQ', shrink=.8)
fig.suptitle(f'Public surface-feed export: {nland:,} land quads + {nrim} inlet edges\n'
             'Controlled fitted geometry is now the default; no stair-stepped mask')
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=170)
plt.close(fig)
manifest.write_text(json.dumps({'mesh': str(mesh), 'mesh_sha256': hashlib.sha256(mesh.read_bytes()).hexdigest(),
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'surface_cells': nland, 'rim_edges': nrim, 'minimum_surface_corner': float(quality.min())}, indent=2)+'\n')
