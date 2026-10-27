"""Render saved fields without raster anti-aliasing cracks between cells."""
from pathlib import Path
import hashlib
import json
import runpy
import sys
import matplotlib

source = Path(sys.argv[1]).resolve()
out = source / 'visuals'
if out.exists():
    raise SystemExit(f'refusing to overwrite reviewed visuals: {out}')
out.mkdir()
inputs = [source / name for name in ('comparison.tsv', 'node-movement.tsv', 'controlled-4.vtk')]
inputs += list(source.glob('patch-*.tsv'))
for path in inputs:
    (out / path.name).symlink_to(Path('..') / path.name)
matplotlib.rcParams['patch.antialiased'] = False
script = Path(__file__).with_name('plot.py')
sys.argv = [str(script), str(out)]
runpy.run_path(str(script), run_name='__main__')
(out / 'rendering.json').write_text(json.dumps({
    'purpose': 'Remove white anti-aliasing seams, not change numerical fields',
    'patch.antialiased': False,
    'wrapper_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'plot_script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
}, indent=2)+'\n')
