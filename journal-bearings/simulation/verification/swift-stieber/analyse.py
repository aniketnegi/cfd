"""Check one Swift-Stieber verification run without archived review inputs."""
import sys
import csv
import hashlib
import json
import math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from command import execute, parser, require
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
        require(hashlib.sha256((root / name).read_bytes()).hexdigest() == expected,
                f'source hash mismatch: {name}')
    sources = [name for name in ['results.tsv', 'cylindrical.tsv', 'mms.tsv',
                                'fd-comparison.tsv', 'fd-mms.tsv', 'fv-flux-checks.tsv', 'long-bearing.tsv',
                                'zero-angle.tsv', 'force-analytic.tsv', 'fitted-mms.tsv',
                                'fitted-refinement.tsv', 'fitted-algebraic.tsv',
                                'unported-grid.tsv', 'unported-agreement.tsv']
               if (directory / name).exists()]
    require(sources, 'No supported result table in run directory')
    checks = {}
    for name in sources:
        data = rows(directory / name)
        require(data, f'{name}: no result rows')
        for row in data:
            for column in ['positive_pressure_residual', 'positive_pressure_relative_residual',
                           'complementarity', 'complementarity_relative',
                           'flux_reconstruction_relative', 'defect_accounting_relative']:
                if column in row:
                    require(abs(row[column]) < 1e-10, f'{name}: {column} exceeds tolerance')
        checks[name] = dict(rows=len(data))
        if name == 'unported-agreement.tsv':
            require(all(row['relative_L2'] < .001 and row['relative_Linf'] < .002
                        and abs(row['relative_peak_difference']) < .001 for row in data),
                    'fine Q4/FV midplane agreement failed')
            checks[name]['comparisons'] = data
        if name == 'unported-grid.tsv':
            reference = root / 'expts/gangrade-figure6/reference.toml'
            require(hashlib.sha256(reference.read_bytes()).hexdigest() == metadata['reference_sha256'],
                    'Gangrade reference hash mismatch')
            summary = []
            for ecc in sorted({row['epsilon'] for row in data}):
                family = [row for row in data if row['epsilon'] == ecc]
                results = {}
                for method in ('Q4', 'FV'):
                    levels = sorted([row for row in family if row['method'] == method], key=lambda row: row['nt'])
                    hierarchy = [row for row in levels if row['nt'] != 512]
                    require(len(hierarchy) == 6 and hierarchy[0]['nt'] == 12 and hierarchy[0]['nz'] == 4,
                            'missing coarse-to-fine hierarchy')
                    require(all(b['nt'] == 2*a['nt'] and b['nz'] == 2*a['nz']
                                for a,b in zip(hierarchy, hierarchy[1:])), 'grid counts must double')
                    for row in levels:
                        expected = row['nt'] * (row['nz'] + 1) if method == 'Q4' else row['nt']*row['nz']
                        require(row['pressure_dofs'] == expected and row['elements'] == row['nt']*row['nz'],
                                'incorrect nodal/cell pressure count')
                        require(math.isfinite(row['midplane_maximum_pa']) and row['midplane_maximum_pa'] > 0
                                and math.isfinite(row['alpha0_max_pa']), 'invalid pressure observable')
                    coarse, fine = levels[0], levels[-1]
                    results[method] = dict(coarse_MPa=coarse['midplane_maximum_pa']/1e6,
                        fine_MPa=fine['midplane_maximum_pa']/1e6,
                        coarse_minus_FEA_MPa=(coarse['midplane_pbar']-coarse['fea_pbar'])/2,
                        coarse_within_visual_reading_band=abs(coarse['midplane_pbar']-coarse['fea_pbar']) <= 1.5,
                        fine_minus_FEA_MPa=(fine['midplane_pbar']-fine['fea_pbar'])/2,
                        alpha0_max_pa=max(row['alpha0_max_pa'] for row in levels))
                summary.append(dict(epsilon=ecc, methods=results))
            checks[name]['families'] = summary
        if name in ('fitted-mms.tsv', 'fitted-refinement.tsv'):
            keys = ('angle_deg', 'front_rad', 'rpm') if name == 'fitted-mms.tsv' else ('angle_deg', 'epsilon', 'rpm')
            families = {}
            for row in data:
                families.setdefault(tuple(row[key] for key in keys), []).append(row)
            summary = []
            for key, family in sorted(families.items()):
                family.sort(key=lambda row: row['nt'])
                require(len(family) >= 2 and all(b['nt'] == 2*a['nt'] and b['nz'] == 2*a['nz']
                        for a,b in zip(family, family[1:])), f'{key}: refinement counts must double')
                fine = family[-1]
                metrics = ('pressure_relative_L2', 'boundary_flow_relative_L2', 'force_relative_L2') if name == 'fitted-mms.tsv' else (
                    'pressure_peak_relative_change', 'force_relative_change', 'boundary_flow_relative_change')
                limit = .003 if name == 'fitted-mms.tsv' else .02
                require(all(math.isfinite(fine[column]) and 0 <= fine[column] < limit for column in metrics),
                        f'{key}: fine comparison exceeds tolerance')
                summary.append(dict(zip(keys, key)) | {column: fine[column] for column in metrics})
            checks[name]['families'] = summary
        if name == 'fitted-algebraic.tsv':
            require(all(row['pressure_relative_L2'] < 1e-7 and row['natural_residual'] < 1e-12
                        and row['control_natural_residual'] < 1e-11 for row in data),
                    'projected-SOR reference comparison failed')
            checks[name].update(data[0])
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
                require(all(b['nt'] == 2*a['nt'] and b['nz'] == 2*a['nz'] for a,b in zip(family, family[1:])),
                        f'{key}: refinement counts must double')
                last = family[-1]
                require(math.isfinite(last['pressure_common_coordinate_mismatch']),
                        f'{key}: missing finite refinement comparison')
                summary.append(dict(angle_deg=key[0], epsilon=key[1], groove=key[2],
                    midplane_maximum_pa=last['midplane_maximum_pa'],
                    global_maximum_pa=last['maximum_pressure_pa'],
                    pressure_common_coordinate_mismatch=last['pressure_common_coordinate_mismatch']))
            checks[name]['families'] = summary
    (directory / 'analysis.json').write_text(json.dumps(dict(physical_validation=False, tables=checks), indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    cli = parser(__doc__)
    args = cli.parse_args()
    if not args.directory.is_dir():
        cli.error(f'not a run directory: {args.directory}')
    raise SystemExit(execute(lambda: analyse(args.directory), args,
                             label=f'analyse {args.directory}'))
