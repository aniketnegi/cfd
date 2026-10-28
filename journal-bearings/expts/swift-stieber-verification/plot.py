"""Plot the raw refinements in a Swift-Stieber verification directory."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse import rows


def plot(directory):
    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    if (directory / 'fd-comparison.tsv').exists():
        data = rows(directory / 'fd-comparison.tsv')
        for column in ['FD_Q0', 'FV_Q0', 'FD_Q0_first', 'FD_Q0_third']:
            ax.plot([row['nt'] for row in data], [row[column] for row in data], '.-', label=column)
        ax.set(xlabel='Circumferential cells/intervals', ylabel='Total two-end Q0', title='Independent FD and shared FV')
    elif (directory / 'results.tsv').exists():
        data = rows(directory / 'results.tsv')
        keys = sorted({(row['angle_deg'], row['epsilon'], row['groove']) for row in data})
        for angle, ecc, groove in keys:
            family = [row for row in data if (row['angle_deg'], row['epsilon'], row['groove']) == (angle, ecc, groove)]
            ax.plot([row['nt'] for row in family], [row['midplane_maximum_pa']/1e6 for row in family], '.-',
                    label=f'{angle:g} deg, eps={ecc:g}, pinned={groove}')
        ax.set(xlabel='Circumferential cells', ylabel='Midplane peak (MPa)', title='Conical numerical boundary cases')
    elif (directory / 'mms.tsv').exists():
        data = rows(directory / 'mms.tsv')
        keys = sorted({(row['angle_deg'], row['rpm'], row['graded']) for row in data})
        for angle, rpm, graded in keys:
            family = [row for row in data if (row['angle_deg'], row['rpm'], row['graded']) == (angle, rpm, graded)]
            ax.loglog([row['nt'] for row in family], [row['pressure_relative_L2'] for row in family], '.-',
                      label=f'{angle:g} deg, {rpm:g} rpm, graded={graded}')
        ax.set(xlabel='Circumferential cells', ylabel='Known-field relative L2 error', title='Manufactured cavitating cone')
    elif (directory / 'cylindrical.tsv').exists():
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    plot(parser.parse_args().directory)
