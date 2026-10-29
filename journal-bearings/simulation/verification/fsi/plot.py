"""Plot fixed-position housing-FSI coupling history and reference-mesh fields."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require


def table(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(rows, f'empty {path.name}')
    return rows


def plot(directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = table(directory / 'results.tsv')
    with (directory / 'inputs.toml').open('rb') as stream:
        inputs = tomllib.load(stream)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    for r in rows:
        if r['mode'] != 'two-way':
            continue
        h = table(directory / f"coupling-{r['angle_deg']}-{r['epsilon']}-{r['nt']}.tsv")
        x = [int(t['iteration']) for t in h]
        label = f"{r['angle_deg']}°, ε={r['epsilon']}, n={r['nt']}"
        for ax, key in zip(axes, ('gap_defect_m', 'pressure_relative_change', 'solid_relative_residual')):
            y = np.array([float(t[key]) for t in h]); valid = np.isfinite(y) & (y > 0)
            ax.semilogy(np.array(x)[valid], y[valid], label=label)
    for ax, title in zip(axes, ('Normal-gap fixed-point defect (m)', 'Relative pressure change', 'Solid equilibrium residual')):
        ax.set(title=title, xlabel='Coupling iteration'); ax.grid(alpha=.3)
    axes[0].legend(fontsize=8)
    fig.savefig(directory / 'coupling.png', dpi=160); plt.close(fig)
    selected = [r for r in rows if r['mode'] == 'two-way' and float(r['angle_deg']) == 10 and float(r['epsilon']) == .6]
    require(selected, 'missing 10degree/epsilon.6 two-way field')
    r = max(selected, key=lambda row: int(row['nt']))
    stem = f"two-way-{r['angle_deg']}-{r['epsilon']}-{r['nt']}"
    fluid = table(directory / (stem+'-fluid.tsv')); solid = table(directory / (stem+'-solid.tsv'))
    n = int(r['solid_nodes'])//(int(r['n_radial'])+1); solid = solid[:n]
    x, y, z = (np.array([float(t[key]) for t in solid]) for key in ('x_m', 'y_m', 'z_m'))
    theta = np.arctan2(y, x); angle = np.deg2rad(10)
    u = np.array([[float(t[key]) for key in ('ux_m', 'uy_m', 'uz_m')] for t in solid])
    normal = np.column_stack((np.cos(angle)*np.cos(theta), np.cos(angle)*np.sin(theta), np.full(n, np.sin(angle))))
    delta = np.sum(normal*u, axis=1)*1e6
    count = 3 if inputs['elastic_journal'] else 2
    fig, axes = plt.subplots(1, count, figsize=(6*count, 4.5), constrained_layout=True)
    a = axes[0].scatter(np.rad2deg([float(t['theta_rad']) for t in fluid]), [float(t['z_m'])*1000 for t in fluid],
                        c=[float(t['pressure_pa'])/1e6 for t in fluid], s=8, marker='s')
    fig.colorbar(a, ax=axes[0], label='Gauge pressure (MPa)')
    b = axes[1].scatter(np.mod(np.rad2deg(theta), 360), z*1000, c=delta, s=8, marker='s')
    fig.colorbar(b, ax=axes[1], label='Normal displacement (µm)')
    if inputs['elastic_journal']:
        journal = table(directory / (stem+'-journal.tsv'))
        xyz = np.array([[float(t[key]) for key in ('x_m', 'y_m', 'z_m')] for t in journal])
        ju = np.array([[float(t[key]) for key in ('ux_m', 'uy_m', 'uz_m')] for t in journal])
        radius = inputs['journal_midpoint_radius_m']+(inputs['length_m']/2-xyz[:, 2])*np.tan(angle)
        outer = np.abs(np.hypot(xyz[:, 0], xyz[:, 1])-radius) < 1e-10
        jt = np.arctan2(xyz[outer, 1], xyz[outer, 0])
        normal_j = np.column_stack((np.cos(angle)*np.cos(jt), np.cos(angle)*np.sin(jt), np.full(sum(outer), np.sin(angle))))
        jdelta = np.sum(normal_j*ju[outer], axis=1)*1e6
        c = axes[2].scatter(np.mod(np.rad2deg(jt), 360), xyz[outer, 2]*1000, c=jdelta, s=12, marker='s')
        fig.colorbar(c, ax=axes[2], label='Outward journal displacement (µm)')
    for ax, title in zip(axes, ('Two-way JFO cell pressure', 'Inner housing nodal displacement', 'Journal: negative is inward')):
        ax.set(title=title, xlabel='Angle (degrees)', ylabel='Axial coordinate (mm)', xlim=(0, 360), ylim=(0, 100))
    fig.savefig(directory / 'fields.png', dpi=160); plt.close(fig)
    (directory / 'plots.json').write_text(json.dumps(dict(
        plot_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    return directory


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: plot(args.directory), args, label='Fixed-position FSI plots'))
