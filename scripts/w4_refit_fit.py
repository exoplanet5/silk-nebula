"""w4_refit_fit: refit the truncated / complete radiative-shock grid (w4_refit_grid, verify_shockmodel code) with the
CORRECTED wave-3 observational inputs (w4_refit_common).  Every swept column along every run is a candidate truncation.

Variants: likelihood A (DESI one-sided Hb only) or B (A + SDSS [O III]/Hb ~ 23 weak central); model Hb/Ha native or
x2.3 ('R88-calibrated', the wave-2 correction for the missing radiation field); 4363 nominal (0.075 +- 0.03) or
'floor' (0.050 +- 0.035).  The [O II] doublet is a separate common term (all pressure-balanced models predict ~1.46).
Acceptable = within Delta chi2 = 4 of the global minimum of either Hb variant (the wave-2 definition), per likelihood.
Also: north fit (w3_northsouth ratios), complete-shock scores, UV (thin) per model.
Outputs: data/derived/w4_refit_fit.json
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_fit.py
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_fit as FIT
import verify_altmodel_lines as AL
import w4_refit_common as C

ROOT = '/Users/mickey/astronomy/silk-nebula'
HB_R88 = FIT.HB_R88_FACTOR
UVK = ['OIII5007', 'CIV1549', 'OIII1663', 'NIII1750', 'NIV1486', 'OIV1400', 'SiIV1400', 'NV1240', 'HeII1640',
       'CIII1909', 'CII2326', 'NeIV2423', 'MgII2800', 'SiIII1892', 'OII2470', '2g']


def track(T, tag):
    p = tag + '__'
    return {k[len(p):]: T[k] for k in T.files if k.startswith(p)}


def track_ratios(tr, hbfac):
    G = {g: tr['F_' + g] for g in FIT.GROUPS}
    r = FIT.ratios(G, hbfac)
    with np.errstate(divide='ignore', invalid='ignore'):
        r['Hb_O3'] = 1.0 / r['O3_Hb']
        r['Ha_O3'] = 1.0 / r['O3_Ha']
        r['O2r_3727'] = r['O2r_O3'] / r['O2_O3']
    return r


def uv_at(tr, i):
    L = {k: float(tr['F_uv_' + k][i]) for k in UVK}
    u = AL.uv_CU_per_R5007(L, ebv=0.0)
    return {'FUV_thin': u['FUV']['total'], 'FUV_nonres': u['FUV']['non_resonant'], 'NUV': u['NUV']['total'],
            'FUV_parts': u['FUV']['parts'], 'NUV_parts': u['NUV']['parts']}


def score_track(tr, obs, hbfac, o4363=None):
    r = track_ratios(tr, hbfac)
    tot, terms, n = C.chi2(r, obs, o4363)
    tot = np.where(np.isfinite(tot) & (tr['N'] > 1e15), tot, 1e9)
    return tot, terms, r, n


VARIANTS = {  # name: (obs, hbfac, o4363)
    'A_nat': (C.SOUTH_A, 1.0, None), 'A_R88': (C.SOUTH_A, HB_R88, None),
    'B_nat': (C.SOUTH_B, 1.0, None), 'B_R88': (C.SOUTH_B, HB_R88, None),
    'Afl_nat': (C.SOUTH_A, 1.0, C.O4363_FLOOR), 'Afl_R88': (C.SOUTH_A, HB_R88, C.O4363_FLOOR),
    'N_nat': (C.NORTH, 1.0, None), 'N_R88': (C.NORTH, HB_R88, None),
}
KEEP = ['O3_Hb', 'O2_O3', 'O4363', 'O2r_O3', 'O2r_3727', 'HeII_O3', 'Ne3_O3', 'N2_O3', 'S2_O3', 'O3_Ha', 'N2_Ha',
        'S2_Ha', 'Ha_Hb', 'CIV_O3', 'R23', 'Hbcoll_frac']

if __name__ == '__main__':
    meta = json.load(open(ROOT + '/data/derived/w4_refit_grid.json'))
    T = np.load(ROOT + '/data/derived/w4_refit_tracks.npz')
    runs = []
    for m in meta['runs']:
        tr = track(T, m['tag'])
        rec = dict(m)
        for vn, (obs, hbf, o43) in VARIANTS.items():
            c2, terms, r, n = score_track(tr, obs, hbf, o43)
            if m['complete']:
                i = len(c2) - 1                       # complete: the full column only
            else:
                i = int(np.argmin(c2))
            near = np.nonzero(c2 <= c2[i] + 1.0)[0] if not m['complete'] else np.array([i])
            rec[vn] = dict(i=i, chi2=float(c2[i]), nterms=n, logN=float(np.log10(tr['N'][i])), T=float(tr['T'][i]),
                           logN_1=[float(np.log10(tr['N'][near.min()])), float(np.log10(tr['N'][near.max()]))],
                           ratios={k: float(r[k][i]) for k in KEEP},
                           terms={k: float(v[i]) for k, v in terms.items()})
            if vn == 'A_nat':
                rec['uv'] = uv_at(tr, i)
                rec['frac_O3_of_full'] = float(tr['F_O3'][i] / tr['F_O3'][-1])
        runs.append(rec)
    trn = [r for r in runs if not r['complete']]
    cmp_ = [r for r in runs if r['complete']]
    out = {'inputs': {'SOUTH_A': C.SOUTH_A, 'SOUTH_B': C.SOUTH_B, 'NORTH': C.NORTH, 'DOUBLET': C.DOUBLET,
                      'O4363_FLOOR': C.O4363_FLOOR, 'UV_OBS': C.UV_OBS, 'HB_R88': HB_R88},
           'n_runs': len(runs), 'variants': {}}
    for base in ('A', 'B', 'Afl', 'N'):
        vn, vr = base + '_nat', base + '_R88'
        mn = min(r[vn]['chi2'] for r in trn)
        mr = min(r[vr]['chi2'] for r in trn)
        acc = [r['tag'] for r in trn if r[vn]['chi2'] <= mn + 4 or r[vr]['chi2'] <= mr + 4]
        acc2 = [r['tag'] for r in trn if r[vn]['chi2'] <= mn + 9.7 or r[vr]['chi2'] <= mr + 9.7]
        best = min(trn, key=lambda r: min(r[vn]['chi2'] - mn, r[vr]['chi2'] - mr))
        byv = {}
        for key in ('v', 'xn', 'b'):
            byv[key] = {str(val): {'nat': min(r[vn]['chi2'] for r in trn if r[key] == val),
                                   'R88': min(r[vr]['chi2'] for r in trn if r[key] == val)}
                        for val in sorted(set(r[key] for r in trn))}
        A_ = [r for r in trn if r['tag'] in acc]

        def rng(f):
            v = np.array([f(r) for r in A_])
            return [float(v.min()), float(np.median(v)), float(v.max())]
        pick = lambda r: r[vn] if r[vn]['chi2'] - mn <= r[vr]['chi2'] - mr else r[vr]
        out['variants'][base] = dict(
            chi2_min_nat=mn, chi2_min_R88=mr, nterms=best[vn]['nterms'], n_acc=len(acc), n_acc_2sig=len(acc2),
            acceptable=acc, best=dict(tag=best['tag'], nat=best[vn], R88=best[vr]), chi2_by=byv,
            ranges=dict(v=rng(lambda r: r['v']), xn=rng(lambda r: r['xn']), b=rng(lambda r: r['b']),
                        logN=rng(lambda r: pick(r)['logN']), T=rng(lambda r: pick(r)['T']),
                        **{k: rng(lambda r, k=k: pick(r)['ratios'][k]) for k in KEEP}),
            complete={r['tag']: {'nat': r[vn]['chi2'], 'R88': r[vr]['chi2'], 'O3_Hb': r[vn]['ratios']['O3_Hb'],
                                 'O2_O3': r[vn]['ratios']['O2_O3'], 'terms': r[vn]['terms']} for r in cmp_})
    # doublet: model value at low density is the LDL for every model; the common penalty
    out['doublet_common_chi2'] = {'model_LDL_1.467': C.chi2_doublet(1.467), 'LDL_1.455_T07': C.chi2_doublet(1.455),
                                  'LDL_1.50_Kal09': C.chi2_doublet(1.50)}
    out['runs'] = runs
    json.dump(out, open(ROOT + '/data/derived/w4_refit_fit.json', 'w'), indent=1)
    for base, v in out['variants'].items():
        print(base, f"chi2min nat {v['chi2_min_nat']:.2f} R88 {v['chi2_min_R88']:.2f} nterms {v['nterms']} "
              f"acc {v['n_acc']} (2sig {v['n_acc_2sig']}) best {v['best']['tag']}")
        print('   ranges', {k: [round(x, 3) for x in v['ranges'][k]] for k in ('v', 'xn', 'b', 'logN', 'T', 'O3_Hb', 'O2_O3', 'O4363', 'O2r_O3', 'HeII_O3', 'O3_Ha')})
        print('   chi2 by v', {k: round(min(x['nat'], x['R88']), 1) for k, x in v['chi2_by']['v'].items()})
        print('   chi2 by b', {k: round(min(x['nat'], x['R88']), 1) for k, x in v['chi2_by']['b'].items()})
        print('   chi2 by x', {k: round(min(x['nat'], x['R88']), 1) for k, x in v['chi2_by']['xn'].items()})
        cm = v['complete']
        print('   complete min chi2', round(min(min(c['nat'], c['R88']) for c in cm.values()), 1),
              'O3/Hb range', [round(min(c['O3_Hb'] for c in cm.values()), 1), round(max(c['O3_Hb'] for c in cm.values()), 1)],
              'O2/O3', [round(min(c['O2_O3'] for c in cm.values()), 2), round(max(c['O2_O3'] for c in cm.values()), 2)])
    print('doublet common', out['doublet_common_chi2'])
