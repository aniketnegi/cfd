"""Extract visible Fig.7 traces from its original 477x392 pixel image."""
import csv
from collections import deque
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from command import execute, parser, require


def components(mask):
    seen = np.zeros(mask.shape, dtype=bool)
    result = []
    for y, x in zip(*np.where(mask)):
        if seen[y, x]:
            continue
        queue = deque([(y, x)]); seen[y, x] = True; points = []
        while queue:
            yy, xx = queue.popleft(); points.append((yy, xx))
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    a, b = yy+dy, xx+dx
                    if 0 <= a < mask.shape[0] and 0 <= b < mask.shape[1] and mask[a, b] and not seen[a, b]:
                        seen[a, b] = True; queue.append((a, b))
        result.append(np.array(points))
    return result


def extract(directory):
    image_path = directory / 'figure7.png'
    image = Image.open(image_path).convert('RGB')
    pixels = np.array(image)
    require(image.size == (477, 392), 'use the original unscaled PDF Fig.7 image')
    targets = ['digitized.csv', 'digitization.json', 'digitization-overlay.png']
    require(not any((directory / name).exists() for name in targets), 'digitization output already exists')
    # Major tick centres read from the original raster, independently of solutions.
    x_pixels = np.array([64., 162., 261., 359., 457.])
    x_values = np.array([0., .02, .04, .06, .08])
    y_pixels = np.array([325., 276., 227., 178., 129., 79., 30.])
    y_values = np.arange(7)*.2e6
    x_fit = np.polyfit(x_pixels, x_values, 1)
    y_fit = np.polyfit(y_pixels, y_values, 1)
    dx, dp = abs(x_fit[0]), abs(y_fit[0])
    red = (pixels[:, :, 0] > 170) & (pixels[:, :, 1] < 100) & (pixels[:, :, 2] < 100)
    roi = np.zeros(red.shape, dtype=bool); roi[9:328, 61:460] = True
    roi[10:62, 310:460] = False  # Legend, not measured data.
    red &= roi
    selected = []; rejected = []
    for points in components(red):
        height, width = np.ptp(points[:, 0])+1, np.ptp(points[:, 1])+1
        if 12 <= len(points) <= 70 and width <= 9 and height <= 9:
            y, x = points.mean(axis=0)
            selected.append(('rom_marker', x, y))
        else:
            rejected.append({'pixels': len(points), 'bounds_yx': [points.min(axis=0).tolist(), points.max(axis=0).tolist()]})
    black = (pixels.max(axis=2) < 90) & roi
    # Do not infer Xi's hidden line beneath red markers.
    red_columns = np.where(red.any(axis=0))[0]
    blocked = set(i+d for i in red_columns for d in (-1, 0, 1))
    for x in range(65, 457):
        if x in blocked:
            continue
        y = np.flatnonzero(black[:, x])
        if len(y):
            selected.append(('xi_visible_curve', float(x), float(np.median(y))))
    require(sum(kind == 'rom_marker' for kind, _, _ in selected) >= 15, 'too few isolated markers')
    rows = []
    for kind, x, y in sorted(selected):
        rows.append(dict(trace=kind, x_m=float(np.polyval(x_fit, x)), pressure_pa=float(np.polyval(y_fit, y)),
                         x_reading_halfwidth_m=dx, pressure_reading_halfwidth_pa=dp, pixel_x=x, pixel_y=y))
    with (directory / 'digitized.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
    drawing = ImageDraw.Draw(image)
    for kind, x, y in selected:
        color = 'blue' if kind == 'rom_marker' else 'lime'
        drawing.ellipse((x-1, y-1, x+1, y+1), fill=color)
    image.save(directory / 'digitization-overlay.png')
    metadata = dict(source='Xi et al.2024 Fig.7, printed p1171/PDF p7', original_image_size=[477, 392],
                    image_rgb_sha256=hashlib.sha256(pixels.tobytes()).hexdigest(),
                    x_tick_pixels=x_pixels.tolist(), x_tick_m=x_values.tolist(),
                    y_tick_pixels=y_pixels.tolist(), y_tick_pa=y_values.tolist(),
                    x_pixel_to_m=x_fit.tolist(), y_pixel_to_pa=y_fit.tolist(),
                    reading_halfwidth_pixels=1, x_reading_halfwidth_m=dx, pressure_reading_halfwidth_pa=dp,
                    uncertainty_scope='Graph-reading envelope only, not PINN error or experimental uncertainty',
                    rom_markers=sum(row['trace'] == 'rom_marker' for row in rows),
                    xi_visible_curve_points=sum(row['trace'] == 'xi_visible_curve' for row in rows),
                    rejected_red_components=rejected)
    (directory / 'digitization.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return directory


if __name__ == '__main__':
    p = parser(__doc__)
    args = p.parse_args()
    if not args.directory.is_dir():
        p.error(f'directory does not exist: {args.directory}')
    sys.exit(execute(lambda: extract(args.directory), args, label='Xi Fig.7 digitization'))
