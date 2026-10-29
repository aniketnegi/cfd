"""Plot conical JFO refinement, manufactured fields and resolved-inlet profiles."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
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
    from matplotlib.ticker import NullLocator
    physical = table(directory / 'physical.tsv'); mms = table(directory / 'manufactured.tsv')
    fd = table(directory / 'independent.tsv')
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for layout, angle, ecc, rpm in sorted({(r['layout'], r['angle_deg'], r['epsilon'], r['rpm']) for r in physical}):
        if layout != 'fitted' or float(rpm) < 0:
            continue
        rows = sorted([r for r in physical if (r['layout'], r['angle_deg'], r['epsilon'], r['rpm']) == (layout, angle, ecc, rpm)], key=lambda r: int(r['nt']))
        x = [int(r['nt']) for r in rows[1:]]
        for ax, key in zip(axes, ('pressure_change', 'fill_change', 'flow_change')):
            ax.loglog(x, [float(r[key]) for r in rows[1:]], 'o-', label=f'{angle}°, ε={ecc}')
    for ax, title in zip(axes, ('Three-plane pressure self-change', 'Three-plane fill self-change', 'Gross boundary-flow vector self-change')):
        ax.set(title=title, xlabel='Circumferential count', ylabel='Fraction')
        ax.set_xticks([128, 256, 512], ['128', '256', '512']); ax.xaxis.set_minor_locator(NullLocator()); ax.grid(alpha=.3)
    axes[0].legend(fontsize=8)
    fig.savefig(directory / 'refinement.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for angle in (0., 10.):
        for rpm in (-2000., 2000.):
            rows = sorted([r for r in mms if r['layout'] == 'fitted' and float(r['angle_deg']) == angle and float(r['rpm']) == rpm], key=lambda r: int(r['nt']))
            for ax, key in zip(axes, ('pressure_relative_L2', 'fill_mean_absolute', 'boundary_flow_relative_L2')):
                ax.loglog([int(r['nt']) for r in rows], [float(r[key]) for r in rows], 'o-', label=f'{angle:g}°, {rpm:g}rpm')
    for ax, title in zip(axes, ('Continuous manufactured pressure L2', 'Manufactured fill absolute mean', 'Manufactured boundary flow error')):
        ax.set(title=title, xlabel='Circumferential count', ylabel='Fraction')
        ax.set_xticks([128, 256, 512], ['128', '256', '512']); ax.xaxis.set_minor_locator(NullLocator()); ax.grid(alpha=.3)
    axes[0].legend(fontsize=8)
    fig.savefig(directory / 'manufactured.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    for angle in (10., 20.):
        rows = sorted([r for r in fd if float(r['angle_deg']) == angle], key=lambda r: int(r['nt']))
        for ax, key in zip(axes, ('pressure_relative_L2', 'fill_mean_absolute')):
            ax.loglog([int(r['nt']) for r in rows], [float(r[key]) for r in rows], 'o-', label=f'{angle:g}°')
    for ax, title in zip(axes, ('Independent nodal FD/FV pressure difference', 'Independent FD/FV fill difference')):
        ax.set(title=title, xlabel='Circumferential count', ylabel='Fraction')
        ax.set_xticks([64, 128, 256], ['64', '128', '256']); ax.xaxis.set_minor_locator(NullLocator()); ax.grid(alpha=.3); ax.legend()
    fig.savefig(directory / 'independent.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    for angle in (0., 10.):
        for ecc in (.6, .8):
            rows = table(directory / f'profiles-{angle}-{ecc}-2000.0.tsv')
            rows = [r for r in rows if float(r['z_m']) == .05]
            theta = np.rad2deg([float(r['theta_rad']) for r in rows])
            axes[0].plot(theta, [float(r['pressure_pa'])/1e6 for r in rows], label=f'{angle:g}°, ε={ecc}')
            axes[1].plot(theta, [float(r['fill_fraction']) for r in rows])
    axes[0].set(ylabel='Gauge pressure (MPa)', title='Actual midplane; 7mm maximum-clearance inlet, 0.5MPa supply')
    axes[1].set(xlabel='Angle (degrees)', ylabel='Reconstructed liquid fraction')
    axes[0].legend(); axes[0].grid(alpha=.3); axes[1].grid(alpha=.3)
    fig.savefig(directory / 'profiles.png', dpi=160); plt.close(fig)
    path = directory / 'fitted-10.0-0.6.tsv'
    field = table(path)
    t = np.rad2deg([float(r['theta_rad']) for r in field]); z = np.array([float(r['z_m']) for r in field])*1000
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    for ax, key, factor, title in zip(axes, ('pressure_pa', 'fill_fraction'), (1e-6, 1.), ('Gauge pressure (MPa)', 'Liquid fraction')):
        # Cell-centre scatter leaves the physical opening empty; no triangulation across it.
        artist = ax.scatter(t, z, c=[float(r[key])*factor for r in field], s=2, marker='s', cmap='viridis', rasterized=True)
        fig.colorbar(artist, ax=ax); ax.set(title=title, ylabel='Axial coordinate (mm)', xlim=(0, 360), ylim=(0, 100))
    axes[1].set_xlabel('Angle (degrees)')
    fig.savefig(directory / 'fields.png', dpi=160); plt.close(fig)
    (directory / 'plots.json').write_text(json.dumps(dict(
        plot_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        command_sha256=hashlib.sha256((Path(__file__).resolve().parents[2] / 'command.py').read_bytes()).hexdigest()), indent=2)+'\n')
    return directory


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: plot(args.directory), args, label='Conical JFO plots'))
