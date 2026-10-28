"""Plot full-film load refinements; output stays in the run directory."""
import sys
import csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(directory):
    with (directory / 'results.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(rows, 'results.tsv contains no solves')
    points = sorted({(float(row['L_D']), float(row['epsilon'])) for row in rows})
    if len(points) > 12:
        points = [point for point in points if point in [(0.7, 0.2), (1.0, 0.8), (0.1, 0.99)]]
    fig, axes = plt.subplots(1, len(points), figsize=(5 * len(points), 4), squeeze=False,
                             layout='constrained')
    for ax, (aspect, ecc) in zip(axes[0], points):
        selected = [row for row in rows if float(row['L_D']) == aspect and float(row['epsilon']) == ecc]
        for model in ['classical', 'production']:
            family = sorted([row for row in selected if row['model'] == model], key=lambda row: float(row['nt']))
            ax.plot([float(row['nt']) for row in family], [float(row['load_ratio']) for row in family], '.-', label=model)
        ax.axhline(float(selected[0]['classical_spectral_load']), color='black', ls='--', label='Classical spectral')
        ax.axhline(float(selected[0]['published']), color='grey', ls=':', label='Printed table')
        ax.set(title=f'L/D={aspect:g}, eccentricity={ecc:g}', xlabel='Circumferential cells', ylabel='Finite / infinite load')
        ax.legend(fontsize=7)
    fig.savefig(directory / 'load-refinement.png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    cli = parser(__doc__)
    args = cli.parse_args()
    if not args.directory.is_dir():
        cli.error(f'not a run directory: {args.directory}')
    raise SystemExit(execute(lambda: plot(args.directory), args,
                             label=f'plot {args.directory}'))
