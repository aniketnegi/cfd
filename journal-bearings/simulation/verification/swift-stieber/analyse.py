"""Check one Swift-Stieber verification run without archived review inputs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import tomllib


def rows(path):
    with path.open() as stream:
        result = []
        for row in csv.DictReader(stream, delimiter='\t'):
            converted = {}
            for key, value in row.items():
                if value in ['true', 'false']:
                    converted[key] = value == 'true'
                else:
                    try:
                        converted[key] = float(value)
                    except ValueError:
                        converted[key] = value
            result.append(converted)
    return result


def analyse(directory):
    metadata = tomllib.loads((directory / 'inputs.toml').read_text())
    root = Path(__file__).resolve().parents[3]
    for name, expected in metadata['source_sha256'].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    sources = [name for name in ['results.tsv', 'cylindrical.tsv', 'mms.tsv',
                                'fd-comparison.tsv', 'fd-mms.tsv', 'fv-flux-checks.tsv', 'long-bearing.tsv',
                                'zero-angle.tsv', 'force-analytic.tsv']
               if (directory / name).exists()]
    assert sources, 'No supported result table in run directory'
    checks = {}
    for name in sources:
        data = rows(directory / name)
        assert data, name
        for row in data:
            for column in ['positive_pressure_residual', 'positive_pressure_relative_residual',
                           'complementarity', 'defect_accounting_relative']:
                if column in row:
                    assert abs(row[column]) < 1e-10, (name, column)
        checks[name] = dict(rows=len(data))
        if name == 'fd-comparison.tsv':
            fine = max(data, key=lambda row: row['nt'])
            checks[name].update(FD_Q0=fine['FD_Q0'], FV_Q0=fine['FV_Q0'],
                pressure_relative_L2=fine['pressure_relative_L2'],
                flow_relative_difference=fine['FD_Q0'] / fine['FV_Q0'] - 1)
        if name == 'results.tsv':
            families = {}
            for row in data:
                families.setdefault((row['angle_deg'], row['epsilon'], row['groove']), []).append(row)
            summary = []
            for key, family in sorted(families.items()):
                family.sort(key=lambda row: row['nt'])
                assert all(b['nt'] == 2*a['nt'] and b['nz'] == 2*a['nz'] for a,b in zip(family, family[1:]))
                last = family[-1]
                assert math.isfinite(last['pressure_common_coordinate_mismatch'])
                summary.append(dict(angle_deg=key[0], epsilon=key[1], groove=key[2],
                    midplane_maximum_pa=last['midplane_maximum_pa'],
                    global_maximum_pa=last['maximum_pressure_pa'],
                    pressure_common_coordinate_mismatch=last['pressure_common_coordinate_mismatch']))
            checks[name]['families'] = summary
    (directory / 'analysis.json').write_text(json.dumps(dict(physical_validation=False, tables=checks), indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    analyse(parser.parse_args().directory)
