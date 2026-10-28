"""Analyse one freshly generated full-film verification directory."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tomllib
from postprocessing import extrapolate


def analyse(directory):
    with (directory / 'results.tsv').open() as stream:
        rows = [{key: value if key == 'model' else float(value)
                 for key, value in row.items()} for row in csv.DictReader(stream, delimiter='\t')]
    groups = {}
    for row in rows:
        groups.setdefault((row['model'], row['L_D'], row['epsilon']), []).append(row)
        assert row['relative_cell_imbalance'] < 1e-10
        assert row['relative_global_imbalance'] < 1e-10
    summary = []
    for (model, aspect, eccentricity), family in sorted(groups.items()):
        item = dict(model=model, L_D=aspect, epsilon=eccentricity)
        try:
            item['estimated_load_ratio'], item['observed_order'] = extrapolate(family, 'load_ratio')
        except ValueError as error:
            item.update(estimated_load_ratio=None, observed_order=None, limitation=str(error))
        fine = max(family, key=lambda row: row['nt'])
        item.update(finest_ratio=fine['load_ratio'], published=fine['published'],
                    classical_reference=fine['classical_spectral_load'])
        summary.append(item)
    root = Path(__file__).resolve().parents[3]
    metadata = tomllib.loads((directory / 'inputs.toml').read_text())
    for name, expected in metadata['source_sha256'].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Checked {len(rows)} solves in {len(groups)} refinement families')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    analyse(parser.parse_args().directory)
