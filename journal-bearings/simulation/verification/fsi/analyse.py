"""Check fixed-position elastic-housing FSI pilot and coupling histories."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tomllib
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require
ROOT = Path(__file__).resolve().parents[3]


def table(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(rows, f'empty {path.name}')
    return rows


def value(row, key):
    result = float(row[key]); require(math.isfinite(result), f'nonfinite {key}')
    return result


def analyse(directory):
    with (directory / 'inputs.toml').open('rb') as stream:
        inputs = tomllib.load(stream)
    for name, expected in inputs['source_sha256'].items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, f'source mismatch: {name}')
    rows = table(directory / 'results.tsv')
    require(len(rows) == (24 if inputs['full'] else 4), 'incomplete FSI study')
    for row in rows:
        require(value(row, 'minimum_gap_m') > 0, 'closed film')
        require(value(row, 'conservative_gap_lower_bound_m') > 0, 'gap positivity not certified by the pilot bound')
        require(0 <= value(row, 'minimum_fill') <= 1, 'invalid liquid content')
        require(value(row, 'liquid_global_balance') < 1e-10, 'liquid imbalance')
        require(value(row, 'maximum_strain') < .01, 'small-strain guard exceeded')
        limit = 1e-9 if row['mode'] == 'one-way' else 1e-5
        require(value(row, 'elastic_residual') <= limit, 'solid force residual')
        require(value(row, 'strain_energy_j') > 0, 'invalid elastic energy')
        if inputs['elastic_journal']:
            require(value(row, 'journal_dofs') > 0 and value(row, 'journal_energy_j') > 0, 'missing elastic journal')
    comparisons = []
    keys = sorted({(r['angle_deg'], r['epsilon'], r['nt']) for r in rows})
    for angle, ecc, count in keys:
        selected = {r['mode']: r for r in rows if (r['angle_deg'], r['epsilon'], r['nt']) == (angle, ecc, count)}
        require(set(selected) == {'one-way', 'two-way'}, 'missing coupling mode')
        one, two = selected['one-way'], selected['two-way']
        history = table(directory / f'coupling-{angle}-{ecc}-{count}.tsv')
        require(len(history) == int(two['coupling_iterations']) and len(history) > 1, 'incomplete coupling history')
        require(value(history[-1], 'gap_defect_m') <= 1e-10, 'gap fixed-point defect')
        require(value(history[-1], 'pressure_relative_change') <= 1e-5, 'pressure not converged')
        require(value(history[-1], 'solid_relative_residual') <= 1e-5, 'solid not converged')
        require(all(value(r, 'minimum_gap_m') > 0 and value(r, 'liquid_global_balance') < 1e-10 for r in history), 'failed iteration')
        comparisons.append(dict(angle_deg=float(angle), epsilon=float(ecc), nt=int(count),
                                one_way_peak_pa=value(one, 'peak_pressure_pa'), two_way_peak_pa=value(two, 'peak_pressure_pa'),
                                peak_change_fraction=value(two, 'peak_pressure_pa')/value(one, 'peak_pressure_pa')-1,
                                two_way_iterations=len(history), two_way_seconds=value(two, 'solve_seconds'),
                                max_normal_displacement_m=value(two, 'max_gap_correction_m'),
                                peak_process_rss_bytes=value(two, 'peak_process_rss_bytes')))
    summary = dict(scope=inputs['scope'], elastic_journal=inputs['elastic_journal'], physical_validation=False,
                   coupled_spatial_convergence_established=False, comparisons=comparisons,
                   explanation='One-way fluid pressure is undeformed; its reported minimum gap is a predicted deformed gap. Two-way fields use the same converged displacement/gap iterate.',
                   analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    path = directory / 'analysis.json'
    path.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return path


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: analyse(args.directory), args, label='Fixed-position FSI analysis'))
