"""Plot the raw refinements in a Swift-Stieber verification directory."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse import rows


def plot(directory, *, table=None):
    supported = ['fd-comparison.tsv', 'results.tsv', 'mms.tsv', 'cylindrical.tsv']
    available = [name for name in supported if (directory / name).is_file()]
    if table is None:
        require(len(available) == 1,
                'select --table; supported tables found: ' + ', '.join(available))
        table = available[0]
    require(table in available, f'unsupported or missing table: {table}')
    require(rows(directory / table), f'{table}: no result rows')
    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    if table == 'fd-comparison.tsv':
        data = rows(directory / 'fd-comparison.tsv')
        for column in ['FD_Q0', 'FV_Q0', 'FD_Q0_first', 'FD_Q0_third']:
            ax.plot([row['nt'] for row in data], [row[column] for row in data], '.-', label=column)
        ax.set(xlabel='Circumferential cells/intervals', ylabel='Total two-end Q0', title='Independent FD and shared FV')
    elif table == 'results.tsv':
        data = rows(directory / 'results.tsv')
        keys = sorted({(row['angle_deg'], row['epsilon'], row['groove']) for row in data})
        for angle, ecc, groove in keys:
            family = [row for row in data if (row['angle_deg'], row['epsilon'], row['groove']) == (angle, ecc, groove)]
            ax.plot([row['nt'] for row in family], [row['midplane_maximum_pa']/1e6 for row in family], '.-',
                    label=f'{angle:g} deg, eps={ecc:g}, pinned={groove}')
        ax.set(xlabel='Circumferential cells', ylabel='Midplane peak (MPa)', title='Conical numerical boundary cases')
    elif table == 'mms.tsv':
        data = rows(directory / 'mms.tsv')
        keys = sorted({(row['angle_deg'], row['rpm'], row['graded']) for row in data})
        for angle, rpm, graded in keys:
            family = [row for row in data if (row['angle_deg'], row['rpm'], row['graded']) == (angle, rpm, graded)]
            ax.loglog([row['nt'] for row in family], [row['pressure_relative_L2'] for row in family], '.-',
                      label=f'{angle:g} deg, {rpm:g} rpm, graded={graded}')
        ax.set(xlabel='Circumferential cells', ylabel='Known-field relative L2 error', title='Manufactured cavitating cone')
    elif table == 'cylindrical.tsv':
        data = rows(directory / 'cylindrical.tsv')
        keys = sorted({(row['model'], row['L_D'], row['epsilon']) for row in data})
        for model, aspect, ecc in keys:
            family = [row for row in data if (row['model'], row['L_D'], row['epsilon']) == (model, aspect, ecc)]
            ax.plot([row['nt'] for row in family], [row['flow_Q0'] for row in family], '.-', label=f'{model} LD={aspect:g}, eps={ecc:g}')
        ax.set(xlabel='Circumferential cells', ylabel='Total two-end Q0', title='Cylindrical pressure-only flow')
    else:
        raise ValueError('Plot a cylinder, cone, MMS or FD comparison run')
    ax.legend(fontsize=6)
    fig.savefig(directory / 'verification.png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    cli = parser(__doc__)
    cli.add_argument('--table', choices=['fd-comparison.tsv', 'results.tsv', 'mms.tsv', 'cylindrical.tsv'],
                     help='table to plot when the directory contains several supported tables')
    args = cli.parse_args()
    if not args.directory.is_dir():
        cli.error(f'not a run directory: {args.directory}')
    raise SystemExit(execute(lambda: plot(args.directory, table=args.table), args,
                             label=f'plot {args.directory}'))
