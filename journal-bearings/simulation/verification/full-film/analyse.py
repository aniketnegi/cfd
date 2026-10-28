"""Analyse one freshly generated full-film verification directory."""
import sys
import csv
import hashlib
import json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require
import tomllib
from postprocessing import extrapolate


def analyse(directory):
    with (directory / 'results.tsv').open() as stream:
        rows = [{key: value if key == 'model' else float(value)
                 for key, value in row.items()} for row in csv.DictReader(stream, delimiter='\t')]
    require(rows, 'results.tsv contains no solves')
    groups = {}
    for row in rows:
        groups.setdefault((row['model'], row['L_D'], row['epsilon']), []).append(row)
        require(row['relative_cell_imbalance'] < 1e-10, 'cell imbalance exceeds tolerance')
        require(row['relative_global_imbalance'] < 1e-10, 'global imbalance exceeds tolerance')
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
        require(hashlib.sha256((root / name).read_bytes()).hexdigest() == expected,
                f'source hash mismatch: {name}')
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Checked {len(rows)} solves in {len(groups)} refinement families')
    return summary


if __name__ == '__main__':
    cli = parser(__doc__)
    args = cli.parse_args()
    if not args.directory.is_dir():
        cli.error(f'not a run directory: {args.directory}')
    raise SystemExit(execute(lambda: analyse(args.directory), args,
                             label=f'analyse {args.directory}'))
