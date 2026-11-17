"""Helpers: load the Fig 1 (img-001) sky mapping; convert Fig 1 page coordinates -> image pixels -> RA/Dec."""
import json
import numpy as np
from geometry_imreg import PolyMap

ROOT = '/Users/mickey/astronomy/silk-nebula'
M929 = (306.3302917480469, -0.036332566291093826, 0.1354925036430359, 1074.4342041015625, 249.78643798828125, -399.9158935546875)
W1, H1 = 1508, 5274


def load_fig1_map():
    js = json.load(open(f'{ROOT}/data/derived/geometry_fig1_wcs.json'))
    pm = PolyMap(js['ra0'], js['de0'], js['order'])
    pm.x0, pm.y0, pm.s = js['x0'], js['y0'], js['s']
    pm.cx, pm.cy = np.array(js['cx']), np.array(js['cy'])
    return pm


def page_to_img1(x, y):
    a, b, c, d, e, f = M929
    A = np.array([[a, c], [b, d]])
    uv = np.linalg.solve(A, np.vstack([np.asarray(x, float) - e, np.asarray(y, float) - f]))
    return uv[0] * W1 - 0.5, uv[1] * H1 - 0.5
