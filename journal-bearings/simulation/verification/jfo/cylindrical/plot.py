"""Plot finite-cylinder JFO refinement, historical comparisons and cavity fields."""
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from command import execute, parser, require


def rows(directory, name):
    with (directory / name).open() as stream:
        result = list(csv.DictReader(stream, delimiter='\t'))
    require(result, f'empty {name}')
    return result


def plot(directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    physical = rows(directory, 'results.tsv')
    mms = rows(directory, 'manufactured.tsv')
    independent = rows(directory, 'independent.tsv')
    profiles = rows(directory, 'profiles.tsv')
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for graded in ('false', 'true'):
        selected = sorted([r for r in mms if r['graded'] == graded and float(r['rpm']) > 0], key=lambda r: int(r['nt']))
        for ax, key in zip(axes[0], ('pressure_relative_L2', 'fill_mean_absolute')):
            ax.loglog([int(r['nt']) for r in selected], [float(r[key]) for r in selected], 'o-', label='graded' if graded == 'true' else 'uniform')
    fd_cases = sorted({(float(r['epsilon']), float(r['pv0'])) for r in independent})
    for ecc, pv in fd_cases:
        selected = sorted([r for r in independent if float(r['epsilon']) == ecc and float(r['pv0']) == pv], key=lambda r: int(r['nt']))
        for ax, key in zip(axes[1], ('pressure_relative_L2', 'fill_mean_absolute')):
            ax.loglog([int(r['nt']) for r in selected], [float(r[key]) for r in selected], 'o-', label=f'epsilon={ecc:g}, pv0={pv:g}')
    for ax, title in zip(axes.flat, ('Continuous MMS pressure L2', 'Continuous MMS fill mean absolute',
                                    'Independent FD/FV pressure L2 difference', 'Independent FD/FV fill mean difference')):
        ax.set(title=title, xlabel='Circumferential count', ylabel='Fraction');ax.grid(True, which='both', alpha=.3);ax.legend()
    fig.savefig(directory / 'verification.png', dpi=170);plt.close(fig)

    aspects = sorted({float(r['L_D']) for r in physical}, reverse=True)
    fig, axes = plt.subplots(len(aspects), 3, squeeze=False, figsize=(12, 3.3*len(aspects)), constrained_layout=True)
    for k, aspect in enumerate(aspects):
        for pv, color in ((-1., 'tab:blue'), (-.1, 'tab:orange')):
            for model, style in (('classical', '-'), ('production', '--')):
                selected = sorted([r for r in physical if r['model'] == model and float(r['L_D']) == aspect and float(r['pv0']) == pv and int(r['nt']) == 512], key=lambda r: float(r['epsilon']))
                for ax, key in zip(axes[k], ('P0', 'beta_deg', 'Q0')):
                    x = [float(r['epsilon']) for r in selected]
                    ax.plot(x, [float(r[key]) for r in selected], style, color=color, label=f'{model}, pv0={pv:g}')
                    if model == 'classical':
                        ax.plot(x, [float(r['reference_' + key]) for r in selected], 'x', color=color, ms=8, label=f'published, pv0={pv:g}')
        for ax, key in zip(axes[k], ('Load P0', 'Attitude (degrees)', 'Gross outward Q0')):
            ax.set(title=f'L/D={aspect:g}: {key}', xlabel='Eccentricity ratio');ax.grid(True, alpha=.3)
    axes[0, 0].legend(fontsize=7)
    fig.savefig(directory / 'reference.png', dpi=170);plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for k, pv in enumerate((-1., -.1)):
        for model, style in (('classical', '-'), ('production', '--')):
            selected = sorted([r for r in profiles if r['model'] == model and float(r['L_D']) == 1 and float(r['epsilon']) == .6 and float(r['pv0']) == pv], key=lambda r: float(r['phi_rad']))
            x = np.degrees([float(r['phi_rad']) for r in selected])
            for ax, key in zip(axes[k], ('p0', 'fill_fraction')):
                ax.plot(x, [float(r[key]) for r in selected], style, label=model)
        axes[k, 0].axhline(pv, color='gray', lw=.7)
        for ax, key in zip(axes[k], ('Pressure p0', 'Liquid fill')):
            ax.set(title=f'pv0={pv:g}: {key}, midplane', xlabel='phi (degrees)');ax.grid(True, alpha=.3);ax.legend()
    fig.savefig(directory / 'profiles.png', dpi=170);plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    for k, pv in enumerate((-1., -.1)):
        field = np.genfromtxt(directory / f'field-{pv}.tsv', delimiter='\t', names=True)
        theta = np.unique(field['theta_rad']);z = np.unique(field['z_m'])
        phi = (theta + np.pi) % (2*np.pi) - np.pi;order = np.argsort(phi)
        for ax, key in zip(axes[k], ('p0', 'fill_fraction')):
            values = field[key].reshape(len(z), len(theta))[:, order]
            picture = ax.pcolormesh(np.degrees(phi[order]), (z-.05)/.05, values, shading='nearest', cmap='viridis')
            ax.set(title=f'Classical pv0={pv:g}: {key}', xlabel='phi (degrees)', ylabel='zeta=z/R')
            fig.colorbar(picture, ax=ax)
    fig.savefig(directory / 'fields.png', dpi=170);plt.close(fig)
    return directory


def main():
    p = parser(__doc__)
    args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    return execute(lambda: plot(args.directory), args, label='cylindrical JFO plots')


if __name__ == '__main__':
    sys.exit(main())
