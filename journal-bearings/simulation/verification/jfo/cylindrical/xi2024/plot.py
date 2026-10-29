"""Plot the conditional Xi2024 pressure comparison and numerical controls."""
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from command import execute, parser, require


def rows(path, delimiter='\t'):
    with path.open() as stream:
        result = list(csv.DictReader(stream, delimiter=delimiter))
    require(result, f'empty {path.name}')
    return result


def plot(directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    profiles = rows(directory / 'profiles.tsv')
    reference = rows(Path(__file__).with_name('reference.csv'), ',')
    mms = rows(directory / 'manufactured.tsv'); independent = rows(directory / 'independent.tsv')
    physical = rows(directory / 'results.tsv')
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    for n, color in ((64, '.7'), (256, '.4'), (1024, 'tab:blue')):
        selected = sorted([r for r in profiles if r['model'] == 'classical' and int(r['nt']) == n], key=lambda r: float(r['x_m']))
        for ax, key, scale in ((axes[0], 'pressure_absolute_pa', 1e-6), (axes[1], 'fill_fraction', 1.)):
            ax.plot([1e3*float(r['x_m']) for r in selected], [scale*float(r[key]) for r in selected], label=f'FV classical {n}')
    selected = sorted([r for r in profiles if r['model'] == 'production' and int(r['nt']) == 1024], key=lambda r: float(r['x_m']))
    for ax, key, scale in ((axes[0], 'pressure_absolute_pa', 1e-6), (axes[1], 'fill_fraction', 1.)):
        ax.plot([1e3*float(r['x_m']) for r in selected], [scale*float(r[key]) for r in selected], '--', color='tab:orange', label='FV production1024')
    for trace, color, marker, label in (('xi_visible_curve', 'black', '.', 'Xi2024 visible line'),
                                        ('rom_marker', 'red', 'o', 'Rom markers reproduced in Xi')):
        selected = [r for r in reference if r['trace'] == trace]
        axes[0].errorbar([1e3*float(r['x_m']) for r in selected], [1e-6*float(r['pressure_pa']) for r in selected],
                         xerr=[1e3*float(r['x_reading_halfwidth_m']) for r in selected],
                         yerr=[1e-6*float(r['pressure_reading_halfwidth_pa']) for r in selected],
                         fmt=marker, color=color, ms=3, elinewidth=.4, alpha=.65, label=label)
    axes[0].set(title='Xi2024 Fig.7: fixed cosine-gap/maximum-gap seam convention', ylabel='Absolute pressure (MPa)')
    axes[1].set(title='Computed liquid content — no published content array comparison', ylabel='Liquid fraction')
    for ax in axes:
        ax.set_xlabel('Circumferential coordinate x (mm)'); ax.grid(True, alpha=.3); ax.legend(fontsize=8)
    fig.savefig(directory / 'comparison.png', dpi=170);plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for graded in ('false', 'true'):
        selected = sorted([r for r in mms if r['graded'] == graded and float(r['rpm']) > 0], key=lambda r: int(r['nt']))
        axes[0].loglog([int(r['nt']) for r in selected], [float(r['pressure_relative_L2']) for r in selected], 'o-', label='graded' if graded == 'true' else 'uniform')
    for ax, key in ((axes[1], 'pressure_relative_L2'), (axes[2], 'fill_mean_absolute')):
        ax.loglog([int(r['nt']) for r in independent], [float(r[key]) for r in independent], 'o-', label='Independent nodal FD/FV')
    for ax, title in zip(axes, ('Continuous full-film MMS pressure L2', 'Pressure L2 cross-difference', 'Fill mean absolute cross-difference')):
        ax.set(title=title, xlabel='Circumferential count', ylabel='Fraction');ax.grid(True, which='both', alpha=.3);ax.legend(fontsize=8)
    fig.savefig(directory / 'verification.png', dpi=170);plt.close(fig)

    data = np.genfromtxt(directory / 'fields-classical.tsv', delimiter='\t', names=True)
    x = np.unique(data['x_m']); y = np.unique(data['y_m'])
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    for ax, key, scale, title in ((axes[0], 'pressure_absolute_pa', 1e-6, 'Absolute pressure (MPa)'),
                                 (axes[1], 'fill_fraction', 1., 'Liquid fraction')):
        picture = ax.pcolormesh(x*1e3, y*1e3, scale*data[key].reshape(len(y), len(x)), shading='nearest', cmap='viridis')
        ax.set(title=title, xlabel='x (mm)', ylabel='Axial y (mm)');fig.colorbar(picture, ax=ax)
    fig.savefig(directory / 'fields.png', dpi=170);plt.close(fig)
    return directory


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: plot(args.directory), args, label='Xi2024 JFO plots'))
