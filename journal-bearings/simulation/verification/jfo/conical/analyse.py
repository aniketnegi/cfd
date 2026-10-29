"""Check conical and fitted-rim JFO transport verification outputs."""
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


def table(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(rows, f'empty table: {path.name}')
    return rows


def number(row, key):
    value = float(row[key])
    require(math.isfinite(value), f'nonfinite {key}')
    return value


def groups(rows, keys):
    result = {}
    for row in rows:
        result.setdefault(tuple(row[k] for k in keys), []).append(row)
    for family in result.values():
        family.sort(key=lambda r: int(r['nt']))
    return result


def analyse(directory):
    with (directory / 'inputs.toml').open('rb') as stream:
        inputs = tomllib.load(stream)
    for name, expected in inputs['source_sha256'].items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, f'source mismatch: {name}')
    full = inputs['full']
    physical = table(directory / 'physical.tsv')
    mms = table(directory / 'manufactured.tsv')
    fd = table(directory / 'independent.tsv')
    algebra = table(directory / 'algebraic.tsv')
    limits = table(directory / 'cylindrical-limit.tsv')
    require(len(physical) == (60 if full else 42), 'incomplete physical study')
    require(len(mms) == (36 if full else 32), 'incomplete manufactured study')
    require(len(fd) == 6 and len(algebra) == 2 and len(limits) == 5, 'incomplete controls')
    for row in physical:
        for key in ('local_balance', 'global_balance', 'complementarity'):
            require(abs(number(row, key)) < 1e-10, f'failed {key}')
        require(0 <= number(row, 'minimum_fill') <= 1 and number(row, 'peak_pa') > 0, 'invalid state')
        out = sum(number(row, f'gross_out_{tag}_m3_s') for tag in ('large', 'small', 'feed'))
        incoming = sum(number(row, f'gross_in_{tag}_m3_s') for tag in ('large', 'small', 'feed'))
        require(abs(out-incoming) / max(out+incoming, 1e-30) < 1e-10, 'gross boundary accounting failed')
    finest = []
    physical_groups = groups(physical, ('layout', 'angle_deg', 'epsilon', 'rpm'))
    require(len(physical_groups) == 18, 'missing physical families')
    for key, rows in physical_groups.items():
        expected = ([64, 128, 256, 512] if full else [64, 128, 256]) if key[0] == 'tensor' else ([128, 256, 512] if full else [128, 256])
        require([int(r['nt']) for r in rows] == expected, f'wrong physical refinement: {key}')
        require(all(int(r['nz']) == int(r['nt'])//2 for r in rows), 'axial refinement is inconsistent')
        row = rows[-1]
        require(number(row, 'pressure_change') < .01 and number(row, 'fill_change') < .01, 'large final field/profile change')
        require(number(row, 'flow_change') < .02 and number(row, 'force_change') < .01, 'large final observable change')
        finest.append(row)
    manufactured = []
    for key, rows in groups(mms, ('layout', 'angle_deg', 'rpm')).items():
        expected = [128, 256, 512] if key[0] == 'fitted' and full else [128, 256] if key[0] == 'fitted' else [64, 128, 256]
        require([int(r['nt']) for r in rows] == expected, f'incomplete MMS family: {key}')
        errors = [number(r, 'pressure_relative_L2') for r in rows]
        fill = [number(r, 'fill_mean_absolute') for r in rows]
        require(errors[-1] < .003 and fill[-1] < .003, 'MMS error exceeds target')
        require(all(b < .8*a for a, b in zip(errors, errors[1:])), 'pressure MMS does not refine')
        require(all(b < .7*a for a, b in zip(fill, fill[1:])), 'fill MMS does not refine')
        for row in rows:
            require(number(row, 'local_balance') < 1e-10 and number(row, 'complementarity') < 1e-10, 'MMS balance/closure failed')
            require(number(row, 'minimum_pressure') >= 0 and 0 <= number(row, 'minimum_fill') <= number(row, 'maximum_fill') <= 1, 'invalid MMS state')
        require(number(rows[-1], 'boundary_flow_relative_L2') < .01 and number(rows[-1], 'force_relative_L2') < .003, 'MMS observable error')
        manufactured.append(dict(layout=key[0], angle_deg=float(key[1]), rpm=float(key[2]),
                                 finest=rows[-1], pressure_orders=[math.log2(a/b) for a, b in zip(errors, errors[1:])],
                                 fill_orders=[math.log2(a/b) for a, b in zip(fill, fill[1:])]))
    require(len(manufactured) == 12, 'missing MMS families')
    independent = []
    for key, rows in groups(fd, ('angle_deg',)).items():
        require([int(r['nt']) for r in rows] == [64, 128, 256], 'incomplete FD refinement')
        require(all(number(r, 'residual') < 1e-11 for r in rows), 'FD convergence failed')
        require(number(rows[-1], 'pressure_relative_L2') < .003 and number(rows[-1], 'fill_mean_absolute') < .006, 'FD/FV disagreement exceeds target')
        independent.append(rows[-1])
    for row in algebra:
        require(number(row, 'pressure_relative_L2') < 1e-6 and number(row, 'fill_Linf') < 1e-5, 'same-matrix solver disagreement')
    require(number(limits[-1], 'pressure_relative_L2') < 1e-6, 'cylinder limit failed')
    require(all(number(b, 'pressure_relative_L2') < number(a, 'pressure_relative_L2') for a, b in zip(limits, limits[1:])), 'cylinder limit does not decrease')
    summary = dict(physical_validation=False, published_conical_JFO_field_matched=False,
                   physical_rows=len(physical), manufactured_rows=len(mms), independent_rows=len(fd),
                   physical_self_change_sampling=dict(tensor='actual coarse cell centres', fitted='common land samples at actual z=25,50,75mm; not a whole-field norm'),
                   max_local_balance=max(number(r, 'local_balance') for r in physical),
                   max_global_balance=max(number(r, 'global_balance') for r in physical),
                   finest=finest, manufactured=manufactured, independent_finest=independent,
                   algebraic=algebra, cylindrical_limit=limits,
                   analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   command_sha256=hashlib.sha256((Path(__file__).resolve().parents[2] / 'command.py').read_bytes()).hexdigest())
    destination = directory / 'analysis.json'
    destination.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return destination


if __name__ == '__main__':
    p = parser(__doc__); args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'run directory does not exist: {args.directory}')
    sys.exit(execute(lambda: analyse(args.directory), args, label='Conical JFO analysis'))
