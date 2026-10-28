import math


def extrapolate(rows, column):
    rows = sorted(rows, key=lambda row: row['nt'])[-3:]
    if len(rows) != 3:
        raise ValueError('Richardson extrapolation requires three mesh levels')
    for coarse, fine in zip(rows, rows[1:]):
        if any(coarse[axis] <= 0 or fine[axis] != 2*coarse[axis] for axis in ['nt', 'nz']):
            raise ValueError('Richardson requires exact factor-two refinement in both directions')
    a, b, c = [row[column] for row in rows]
    if not all(math.isfinite(value) for value in [a, b, c]) or (a-b)*(b-c) <= 0:
        raise ValueError('A single leading-error term requires nonzero same-sign differences')
    order = math.log((a-b)/(b-c), 2)
    if order <= 0:
        raise ValueError('Richardson estimates require positive observed order')
    return c+(c-b)/(2**order-1), order
