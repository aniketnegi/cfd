"""Analyse Xi2024 pressure comparison; graph agreement is not a pass criterion."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tomllib

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from command import execute, parser, require
ROOT = Path(__file__).resolve().parents[5]


def table(path, delimiter='\t'):
    with path.open() as stream:
        result = list(csv.DictReader(stream, delimiter=delimiter))
    require(result, f'empty {path.name}')
    return result


def finite(row, key):
    value = float(row[key]); require(math.isfinite(value), f'nonfinite {key}')
    return value


def graph_comparison(profile, reference):
    x = np.array([float(r['x_m']) for r in profile])
    p = np.array([float(r['pressure_absolute_pa']) for r in profile])
    require(np.all(np.diff(x) > 0), 'profile coordinates must increase')
    errors = []; targets = []; uncovered = []
    for row in reference:
        point = float(row['x_m']); value = float(row['pressure_pa'])
        dx, dp = float(row['x_reading_halfwidth_m']), float(row['pressure_reading_halfwidth_pa'])
        low, high = max(x[0], point-dx), min(x[-1], point+dx)
        sample_x = np.r_[low, x[(x > low) & (x < high)], high]
        sample_p = np.interp(sample_x, x, p)
        errors.append(float(np.interp(point, x, p))-value); targets.append(value)
        uncovered.append(max(0., sample_p.min()-(value+dp), (value-dp)-sample_p.max()))
    return dict(points=len(errors), sample_relative_L2=float(np.linalg.norm(errors)/np.linalg.norm(targets)),
                sample_rms_pa=float(np.sqrt(np.mean(np.square(errors)))), sample_Linf_pa=float(np.max(np.abs(errors))),
                within_reading_box=int(sum(error <= 1e-8 for error in uncovered)), max_uncovered_pa=float(max(uncovered)))


def analyse(directory):
    with (directory / 'inputs.toml').open('rb') as stream:
        inputs = tomllib.load(stream)
    for name, expected in inputs['source_sha256'].items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, f'source mismatch: {name}')
    physical = table(directory / 'results.tsv'); mms = table(directory / 'manufactured.tsv')
    independent = table(directory / 'independent.tsv'); algebra = table(directory / 'algebraic.tsv')
    profiles = table(directory / 'profiles.tsv')
    reference = table(Path(__file__).with_name('reference.csv'), ',')
    require(len(physical) == 10 and len(mms) == 12 and len(independent) == 4 and len(algebra) == 2, 'incomplete study')
    for row in physical:
        for key in ('local_balance', 'local_face_balance', 'global_balance', 'complementarity'):
            require(abs(finite(row, key)) < 1e-10, f'{key} exceeds tolerance')
        require(0 <= finite(row, 'minimum_fill') <= 1, 'fill outside bounds')
    fine = []
    for model in ('classical', 'production'):
        rows = sorted([r for r in physical if r['model'] == model], key=lambda r: int(r['nt']))
        require([int(r['nt']) for r in rows] == [64, 128, 256, 512, 1024], 'incomplete refinement family')
        require([int(r['nz']) for r in rows] == [16, 32, 64, 128, 256], 'incorrect axial refinement')
        require(finite(rows[-1], 'pressure_mismatch_previous') < .01, 'large final pressure self-change')
        require(finite(rows[-1], 'fill_mismatch_previous') < .01, 'large final fill self-change')
        fine.append(rows[-1])
    groups = {}
    for row in mms:
        require(finite(row, 'local_balance') < 1e-10, 'MMS continuity exceeds tolerance')
        require(finite(row, 'minimum_pressure') > 0 and finite(row, 'fill_error') == 0, 'MMS state invalid')
        groups.setdefault((row['graded'], row['rpm']), []).append(row)
    mms_summary = []
    for key, rows in groups.items():
        rows.sort(key=lambda r: int(r['nt']))
        require([int(r['nt']) for r in rows] == [64, 128, 256], 'incomplete MMS family')
        errors = [finite(r, 'pressure_relative_L2') for r in rows]
        require(errors[-1] < .002 and errors[-1] < errors[0]/2, 'MMS error/refinement exceeds tolerance')
        mms_summary.append(dict(graded=key[0], rpm=float(key[1]), finest_pressure_relative_L2=errors[-1],
                                orders=[math.log2(a/b) for a, b in zip(errors, errors[1:])]))
    for row in independent:
        require(finite(row, 'residual') < 1e-10, 'FD residual exceeds tolerance')
    require([int(r['nt']) for r in independent] == [64, 128, 256, 512], 'incomplete independent refinement')
    require(finite(independent[-1], 'pressure_relative_L2') < .01 and finite(independent[-1], 'fill_mean_absolute') < .01,
            'independent comparison exceeds tolerance')
    for row in algebra:
        require(finite(row, 'pressure_relative_L2') < 1e-7 and finite(row, 'fill_Linf') < 1e-7, 'algebraic sensitivity exceeds tolerance')
    comparisons = []
    for row in physical:
        selected = sorted([r for r in profiles if r['model'] == row['model'] and r['nt'] == row['nt']], key=lambda r: float(r['x_m']))
        require(len(selected) == int(row['nt'])+2, 'incomplete actual-midplane profile')
        for trace in ('rom_marker', 'xi_visible_curve'):
            targets = [r for r in reference if r['trace'] == trace]
            comparisons.append(dict(model=row['model'], nt=int(row['nt']), trace=trace, **graph_comparison(selected, targets)))
    with (directory / 'comparison.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=comparisons[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(comparisons)
    summary = dict(physical_rows=len(physical), geometry_convention_fully_recovered=False,
                   gap_origin_status=inputs['gap_origin_status'], physical_validation=False,
                   published_pressure_arrays_recovered=False, published_content_comparison=False,
                   reference_agreement_pass_criterion=None, manufactured=mms_summary,
                   max_local_balance=max(float(r['local_balance']) for r in physical),
                   max_local_face_balance=max(float(r['local_face_balance']) for r in physical),
                   max_global_balance=max(float(r['global_balance']) for r in physical),
                   finest=fine, independent_finest=independent[-1], algebraic=algebra,
                   graph_comparisons_finest=[r for r in comparisons if r['nt'] == 1024])
    destination = directory / 'analysis.json'
    destination.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return destination


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: analyse(args.directory), args, label='Xi2024 JFO analysis'))
