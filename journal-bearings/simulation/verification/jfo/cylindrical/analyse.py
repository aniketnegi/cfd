"""Check numerical evidence; report historical differences without fitting them."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from command import execute, parser, require

ROOT = Path(__file__).resolve().parents[4]


def table(directory, name):
    with (directory / name).open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(rows, f'empty {name}')
    return rows


def number(row, name):
    value = float(row[name])
    require(math.isfinite(value), f'nonfinite {name}')
    return value


def analyse(directory):
    with (directory / 'inputs.toml').open('rb') as stream:
        inputs = tomllib.load(stream)
    for path, expected in inputs['source_sha256'].items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected,
                f'source mismatch: {path}')
    physical = table(directory, 'results.tsv')
    manufactured = table(directory, 'manufactured.tsv')
    independent = table(directory, 'independent.tsv')
    algebra = table(directory, 'algebraic.tsv')
    expected = 192 if inputs['full'] else 16
    require(len(physical) == expected, f'expected {expected} physical rows, got {len(physical)}')
    independent_cases = 4 if inputs['full'] else 2
    require(len(manufactured) == 12 and len(independent) == 3*independent_cases
            and len(algebra) == 3*independent_cases, 'incomplete numerical controls')
    families = {}
    for row in physical:
        key = (row['model'], number(row, 'L_D'), number(row, 'epsilon'), number(row, 'pv0'))
        families.setdefault(key, []).append(row)
        for name in ('local_balance', 'global_balance', 'complementarity'):
            require(abs(number(row, name)) < 1e-10, f'{name} exceeds tolerance')
        require(0 <= number(row, 'minimum_fill') <= 1, 'fill outside bounds')
        outflow, inflow = number(row, 'gross_out_m3_s'), number(row, 'gross_in_m3_s')
        require(outflow > 0 and inflow > 0, 'missing bidirectional end exchange')
        require(abs(outflow-inflow)/(outflow+inflow) < 1e-10, 'gross-flow balance exceeds tolerance')
        for name in ('P0', 'Q0', 'beta_deg'):
            number(row, name)
    final = []
    for key, rows in families.items():
        rows.sort(key=lambda row: int(row['nt']))
        require([int(r['nt']) for r in rows] == [64, 128, 256, 512], f'incomplete family {key}')
        require(all(int(b['nz']) == 2*int(a['nz']) for a, b in zip(rows, rows[1:])),
                f'axial counts do not double: {key}')
        last = rows[-1]
        require(number(last, 'pressure_mismatch_previous') < .01, f'large pressure self-change: {key}')
        require(number(last, 'fill_mismatch_previous') < .01, f'large fill self-change: {key}')
        final.append(last)
    mms_groups = {}
    for row in manufactured:
        for name in ('local_balance', 'complementarity'):
            require(abs(number(row, name)) < 1e-10, f'MMS {name} exceeds tolerance')
        require(number(row, 'minimum_pressure') >= 0, 'negative MMS pressure')
        require(0 <= number(row, 'minimum_fill') <= number(row, 'maximum_fill') <= 1,
                'MMS fill outside bounds')
        mms_groups.setdefault((row['graded'], row['rpm']), []).append(row)
    mms_summary = []
    for key, rows in mms_groups.items():
        rows.sort(key=lambda row: int(row['nt']))
        require([int(r['nt']) for r in rows] == [64, 128, 256], 'incomplete MMS family')
        require(number(rows[-1], 'pressure_relative_L2') < .01, 'MMS pressure error exceeds tolerance')
        require(number(rows[-1], 'fill_mean_absolute') < .005, 'MMS fill error exceeds tolerance')
        item = {'graded': key[0], 'rpm': float(key[1])}
        for name in ('pressure_relative_L2', 'fill_mean_absolute'):
            errors = [number(row, name) for row in rows]
            require(errors[-1] < errors[0]/2, f'MMS {name} does not refine')
            item[name] = errors[-1]
            item[name + '_orders'] = [math.log2(a/b) for a, b in zip(errors, errors[1:])]
        mms_summary.append(item)
    fd_fine = [row for row in independent if int(row['nt']) == 256]
    require(len(fd_fine) == independent_cases, 'missing fine independent comparisons')
    for row in fd_fine:
        require(number(row, 'pressure_relative_L2') < .01, 'independent pressure difference exceeds tolerance')
        require(number(row, 'fill_mean_absolute') < .01, 'independent fill difference exceeds tolerance')
    for row in independent:
        require(number(row, 'residual') < 1e-10, 'independent residual exceeds tolerance')
    for row in algebra:
        require(number(row, 'pressure_relative_L2') < 1e-7 and number(row, 'fill_Linf') < 1e-7,
                'algebraic sensitivity exceeds tolerance')
    historical = []
    for row in final:
        if row['model'] != 'classical':
            continue
        historical.append({
            'L_D': float(row['L_D']), 'epsilon': float(row['epsilon']), 'pv0': float(row['pv0']),
            'P0': float(row['P0']), 'P0_difference_percent': 100*(float(row['P0'])/float(row['reference_P0'])-1),
            'beta_deg': float(row['beta_deg']), 'beta_difference_deg': float(row['beta_deg'])-float(row['reference_beta_deg']),
            'Q0': float(row['Q0']), 'Q0_difference_percent': 100*(float(row['Q0'])/float(row['reference_Q0'])-1),
        })
    result = {
        'physical_rows': len(physical), 'historical_rows': len(historical),
        'physical_validation': False, 'historical_acceptance_threshold': None,
        'maximum_local_balance': max(float(r['local_balance']) for r in physical),
        'maximum_global_balance': max(float(r['global_balance']) for r in physical),
        'maximum_finest_pressure_self_change': max(float(r['pressure_mismatch_previous']) for r in final),
        'maximum_finest_fill_self_change': max(float(r['fill_mismatch_previous']) for r in final),
        'manufactured': mms_summary, 'independent_finest': fd_fine, 'algebraic': algebra,
        'historical_comparisons': historical,
    }
    destination = directory / 'analysis.json'
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return destination


def main():
    p = parser(__doc__)
    args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    return execute(lambda: analyse(args.directory), args, label='cylindrical JFO analysis')


if __name__ == '__main__':
    sys.exit(main())
