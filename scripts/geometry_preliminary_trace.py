"""Write a PRELIMINARY silk_trace.json (rough spine through P3-P1-P2 extended N-S to ~10 deg).

Will be overwritten by scripts/geometry_build_trace.py with "preliminary": false.
"""
import json, sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import REGIONS, CENTER, ROOT

# Spine: P3 (south) -> P1 -> P2 (north); quadratic in Dec through the three bright regions, extended to
# a ~10 deg total length (paper: 1x10 deg) centred on the paper centre Dec +16.42.
decs_k = np.array([REGIONS['P3'][1], REGIONS['P1'][1], REGIONS['P2'][1]])
ras_k = np.array([REGIONS['P3'][0], REGIONS['P1'][0], REGIONS['P2'][0]])
c = np.polyfit(decs_k, ras_k, 2)
dec = np.arange(11.5, 21.51, 0.5)
ra = np.polyval(c, dec)
spine = []
for r, d in zip(ra, dec):
    inside = 14.7 <= d <= 17.8
    spine.append({'ra': round(float(r), 3), 'dec': round(float(d), 3),
                  'sigma_deg': 0.3 if inside else 0.8,
                  'source': 'prelim: quadratic through P3-P1-P2' + ('' if inside else ' (extrapolated)')})
out = {
    'version': '0.1-prelim',
    'preliminary': True,
    'spine': spine,
    'half_width_deg': 0.5,
    'halpha_patch': [[116.6, 20.8], [116.0, 20.8], [116.0, 21.8], [116.6, 21.8]],
    'lamost_field': {'ra': 116.3, 'dec': 16.3, 'radius_deg': 2.5, 'center_star': 'HIP 37797 (approx; SIMBAD pending)'},
    'regions': {k: {'ra': v[0], 'dec': v[1]} for k, v in REGIONS.items()},
    'curvature': {'center_ra': None, 'center_dec': None, 'radius_deg': None, 'method': 'not yet fitted'},
    'notes': 'PRELIMINARY rough geometry from published numbers only (paper Table 1, 1x10 deg extent). '
             'Halpha patch position is a placeholder near the assumed north end; do not use quantitatively. '
             'Refined version (preliminary=false) will follow from registered paper figures.',
}
(ROOT / 'data/derived/silk_trace.json').write_text(json.dumps(out, indent=1))
print('wrote prelim', len(spine))
