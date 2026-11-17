"""verify_xrayfuv: robustness of the local X-ray step significance to (a) restricting null placements
to high latitude (min b > 8 deg, where N(HI) and absorption are Silk-like) and (b) a realistic
absorption correction exp(sigma f N_HI) with sigma f = 1.5e-21 and 3e-21 cm^2 (f~0.15-0.3).
Output: data/derived/verify_xrayfuv_xray_highb.json
"""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import *
import verify_xrayfuv_xray_edge as X

ra, dec, _ = load_spine()
E = X.load_erass(); pl = X.Placer(E)
cc = CurveCoords(ra, dec); s_kink = float(np.interp(17.6, cc.dec, cc.s))
segs = {'all': (None, None), 'south': (0.0, s_kink), 'north': (s_kink, float(cc.s[-1]))}
out = {}
for sf in (0.0, 1.5e-21, 3e-21):
    img = E['soft'] * np.exp(sf * np.nan_to_num(E['nhi'].astype(float)))
    for seg, (a, b) in segs.items():
        r = X.run_placements(pl, img, E['mask_full'], 0.1, ra, dec, 'soft', 'PSR', a, b)
        rec = {'real': r['real']}
        for nm in ('rot', 'tr'):
            keep = [s for s in r[nm] if s['bmin'] > 8.0]
            for st in ('c', 'cfrac', 'WE'):
                v = np.array([s[st] for s in keep]); o = r['real'][st]
                rec[f'{nm}_{st}'] = {'n': int(v.size), 'mean': float(v.mean()), 'std': float(v.std(ddof=1)),
                                     'z': float((o - v.mean()) / v.std(ddof=1)), 'p_ge': float((np.sum(v >= o) + 1) / (v.size + 1))}
        out[f'sf{sf:.1e}_{seg}'] = rec
        print(sf, seg, 'cfrac %.3f' % r['real']['cfrac'], 'WE %.3f' % r['real']['WE'],
              {k: (v['n'], round(v['z'], 2), round(v['p_ge'], 3)) for k, v in rec.items() if k != 'real'}, flush=True)
(DER / 'verify_xrayfuv_xray_highb.json').write_text(json.dumps(out, indent=1, default=float))
