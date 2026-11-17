"""verify_shockmodel_north: where along the truncated-shock tracks do the NORTHERN Silk ratios fall?
North (wave 2): [O III]/Hb 6.4 (3.2-12, SDSS), [O II]/[O III] 4.4 (DESI, 3.5-5.2) - 7.3 (SDSS).
Scores only these two ratios (+0.1 dex model, +0.3 dex Hb systematics). Output: data/derived/verify_shockmodel_north.json"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_fit as V
ROOT = '/Users/mickey/astronomy/silk-nebula'
fit = json.load(open(ROOT + '/data/derived/verify_shockmodel_fit.json'))
T = np.load(ROOT + '/data/derived/verify_shockmodel_tracks.npz')
OBSN = {'O3_Hb': (6.4, 0.25, 0.32), 'O2_O3': (5.5, 0.12, 0.10)}
res = []
for r in fit['runs']:
    if r['complete']:
        continue
    tag = r['tag']
    tr = {k.split('__', 1)[1]: T[k] for k in T.files if k.startswith(tag + '__')}
    i, c2, rr, terms, _ = V.best_on_track(tr, obs=OBSN)
    res.append(dict(tag=tag, v=r['v'], xn=r['xn'], b=r['b'], chi2=c2, logN=float(np.log10(tr['N'][i])), T=float(tr['T'][i]),
                    logN_south=float(np.log10(r['best']['N'])), O3_Hb=rr['O3_Hb'], O2_O3=rr['O2_O3'], O4363=rr['O4363'],
                    reached_end=bool(i >= len(tr['N']) - 2)))
res.sort(key=lambda d: d['chi2'])
for d in res[:15]:
    print(d)
same = [d for d in res if d['chi2'] < 2]
print('n with chi2<2:', len(same), ' fraction reaching track end (T=1e4):', np.mean([d['reached_end'] for d in same]) if same else None)
print('dlogN (north - south best) for chi2<2:', np.percentile([d['logN'] - d['logN_south'] for d in same], [5, 50, 95]) if same else None)
json.dump(res, open(ROOT + '/data/derived/verify_shockmodel_north.json', 'w'), indent=1)
