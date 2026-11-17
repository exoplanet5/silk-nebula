"""w3_northsouth_shock.py -- wave-3 (tag northsouth): which change of the verify_shockmodel shock family turns the
southern Silk spectrum into the northern one?

Model tracks: data/derived/verify_shockmodel_tracks.npz (verify_shockmodel's own 1-D NEI shock code; cumulative
face-on line fluxes F_X(N) at n0 = 0.1 cm^-3, flux proportional to n0 at fixed b), run metadata from
verify_shockmodel_fit.json, the 25 'consistent' southern models from verify_shockmodel_summary.json.  All ratios are
energy ratios.  Read only; no model is re-run.

Observed (wave-2/3, energy ratios):
  south: [O II]/[O III] 0.9-1.5, [O III]/Hb >= 15
  north: [O II]/[O III] 3.5-9 (DESI 4.0-5.3, SDSS 7.1-7.7), [O III]/Hb 4-30 (SDSS 7-8, > 2.6-2.9 at 2 sigma;
         DESI Hb/[O III] = -0.06 +- 0.09)
  [O III] SB north/south 0.36-0.45 (SDSS / DESI), GALEX FUV per R[O III] north/south 56/92 = 0.61 +- 0.15,
  NUV per R[O III] 81/49 = 1.66 +- 0.30 (verify_xrayfuv_galex_measure.json, NSNS [O III] scale, Dec 17.6-22).
Scenarios, each anchored on one of the 25 consistent southern models (v_s 100-150 km/s, its log N*):
  L  same shock (same v_s, b, X), later truncation N_N > N*: find N_N where the northern ratios hold;
  V  slower shock (v_N < v_S, same b and X) at the same driving pressure (n0 ~ P / v^2): any N where the
     northern ratios hold;
FUV proxy per [O III] = S_CIV * C IV 1549 + O III] 1663 + He II 1640 (S_CIV = 0.21, the verify_shockmodel
resonance-scattering median; the same S for north and south).
Halpha patch (SDSS, verify_sdss): [O II]/Hb 13, [O III]/Hb 1.8, [S II]/Ha 1.06, [N II]/Ha 0.78, Ha/Hb 2.84
compared with the complete-shock runs (*_cmp) and with the ends (T = 1e4 K) of the truncated runs.
Output: data/derived/w3_northsouth_shock.json
"""
import json

import numpy as np

DER = '/Users/mickey/astronomy/silk-nebula/data/derived'
S_CIV = 0.21
OBS = dict(south_o2o3=(0.9, 1.5), south_o3hb_min=15.0, north_o2o3=(3.5, 9.0), north_o3hb=(4.0, 30.0),
           sb_ratio=(0.36, 0.45), fuv_ratio=(0.61, 0.15), nuv_ratio=(1.66, 0.30))
PATCH = dict(O2_Hb=13.0, O3_Hb=1.8, S2_Ha=1.06, N2_Ha=0.78, Ha_Hb=2.84)


def track(t, tag):
    g = lambda k: t[f'{tag}__{k}']  # noqa: E731
    d = {k: g(k) for k in ('N', 'T', 'F_O3', 'F_Hb', 'F_Ha', 'F_O2', 'F_N2', 'F_S2', 'F_CIV', 'F_O3s',
                           'F_HeII1640', 'F_Hb_coll')}
    d['o2o3'] = d['F_O2'] / d['F_O3']
    d['o3hb'] = d['F_O3'] / d['F_Hb']
    d['fuv'] = (S_CIV * d['F_CIV'] + d['F_O3s'] + d['F_HeII1640']) / d['F_O3']
    return d


def north_ok(d, i):
    return (OBS['north_o2o3'][0] <= d['o2o3'][i] <= OBS['north_o2o3'][1] and
            OBS['north_o3hb'][0] <= d['o3hb'][i] <= OBS['north_o3hb'][1])


def best_north(d, idx):
    """index among idx closest (log) to [O II]/[O III] = 5.5, [O III]/Hb = 7."""
    cost = [np.log(d['o2o3'][i] / 5.5) ** 2 + np.log(d['o3hb'][i] / 7.0) ** 2 for i in idx]
    return idx[int(np.argmin(cost))]


def main():
    t = np.load(f'{DER}/verify_shockmodel_tracks.npz')
    fit = json.load(open(f'{DER}/verify_shockmodel_fit.json'))
    summ = json.load(open(f'{DER}/verify_shockmodel_summary.json'))
    runs = {(r['v'], r['xn'], r['b'], r['ab'], r['complete']): r['tag'] for r in fit['runs']}
    south_models = summ['consistent']
    out = dict(obs=OBS, S_CIV=S_CIV, n_south_models=len(south_models), L=[], V=[])
    for sm in south_models:
        tag = runs.get((sm['v'], sm['xn'], sm['b'], 'asplund09', False))
        if tag is None:
            continue
        d = track(t, tag)
        iS = int(np.argmin(np.abs(np.log10(d['N']) - sm['logN'])))
        base = dict(v=sm['v'], xn=sm['xn'], b=sm['b'], logN_S=float(np.log10(d['N'][iS])), T_S=float(d['T'][iS]),
                    o2o3_S=float(d['o2o3'][iS]), o3hb_S=float(d['o3hb'][iS]), fuv_S=float(d['fuv'][iS]))
        # ---- L: later truncation on the same track
        later = [i for i in range(iS + 1, d['N'].size) if north_ok(d, i)]
        rec = dict(base)
        rec['T_end'] = float(d['T'][-1])
        rec['o2o3_end'] = float(d['o2o3'][-1])
        rec['o3hb_end'] = float(d['o3hb'][-1])
        if later:
            iN = best_north(d, later)
            rec.update(found=True, logN_N=float(np.log10(d['N'][iN])), T_N=float(d['T'][iN]),
                       o2o3_N=float(d['o2o3'][iN]), o3hb_N=float(d['o3hb'][iN]),
                       F_O3_ratio=float(d['F_O3'][iN] / d['F_O3'][iS]), F_O2_ratio=float(d['F_O2'][iN] / d['F_O2'][iS]),
                       F_Hb_ratio=float(d['F_Hb'][iN] / d['F_Hb'][iS]), fuv_ratio=float(d['fuv'][iN] / d['fuv'][iS]),
                       n0G_ratio_needed=[float(OBS['sb_ratio'][0] / (d['F_O3'][iN] / d['F_O3'][iS])),
                                         float(OBS['sb_ratio'][1] / (d['F_O3'][iN] / d['F_O3'][iS]))],
                       logN_range_ok=[float(np.log10(d['N'][min(later)])), float(np.log10(d['N'][max(later)]))])
        else:
            rec['found'] = False
        out['L'].append(rec)
        # ---- V: slower shocks, same b and X, equal pressure (n0 ~ 1/v^2), equal G
        for vN in (70, 80, 90, 100, 110):
            if vN >= sm['v']:
                continue
            tagN = runs.get((vN, sm['xn'], sm['b'], 'asplund09', False))
            if tagN is None:
                continue
            dN = track(t, tagN)
            okN = [i for i in range(dN['N'].size) if north_ok(dN, i)]
            r2 = dict(base, v_N=vN)
            if okN:
                iN = best_north(dN, okN)
                sb = (sm['v'] / vN) ** 2 * dN['F_O3'][iN] / d['F_O3'][iS]
                r2.update(found=True, logN_N=float(np.log10(dN['N'][iN])), T_N=float(dN['T'][iN]),
                          o2o3_N=float(dN['o2o3'][iN]), o3hb_N=float(dN['o3hb'][iN]),
                          sb_ratio_equalP=float(sb),
                          P_ratio_needed=[float(OBS['sb_ratio'][0] / sb), float(OBS['sb_ratio'][1] / sb)],
                          fuv_ratio=float(dN['fuv'][iN] / d['fuv'][iS]),
                          F_O2_ratio_equalP=float((sm['v'] / vN) ** 2 * dN['F_O2'][iN] / d['F_O2'][iS]),
                          F_Hb_ratio_equalP=float((sm['v'] / vN) ** 2 * dN['F_Hb'][iN] / d['F_Hb'][iS]))
            else:
                r2['found'] = False
            out['V'].append(r2)
    # ---- summaries
    L = [r for r in out['L'] if r['found']]
    out['L_summary'] = dict(n_found=len(L), n_tried=len(out['L']),
                            F_O3_ratio=[min(r['F_O3_ratio'] for r in L), max(r['F_O3_ratio'] for r in L)] if L else None,
                            fuv_ratio=[min(r['fuv_ratio'] for r in L), max(r['fuv_ratio'] for r in L)] if L else None,
                            F_Hb_ratio=[min(r['F_Hb_ratio'] for r in L), max(r['F_Hb_ratio'] for r in L)] if L else None,
                            dlogN=[min(r['logN_N'] - r['logN_S'] for r in L), max(r['logN_N'] - r['logN_S'] for r in L)]
                            if L else None,
                            T_N=[min(r['T_N'] for r in L), max(r['T_N'] for r in L)] if L else None,
                            n0G_needed=[min(r['n0G_ratio_needed'][0] for r in L), max(r['n0G_ratio_needed'][1] for r in L)]
                            if L else None)
    V = [r for r in out['V'] if r['found']]
    byv = {}
    for r in V:
        byv.setdefault(r['v_N'], []).append(r)
    out['V_summary'] = {str(k): dict(n=len(v), sb_ratio_equalP=[min(x['sb_ratio_equalP'] for x in v),
                                                                  max(x['sb_ratio_equalP'] for x in v)],
                                     P_needed=[min(x['P_ratio_needed'][0] for x in v), max(x['P_ratio_needed'][1] for x in v)],
                                     fuv_ratio=[min(x['fuv_ratio'] for x in v), max(x['fuv_ratio'] for x in v)],
                                     o3hb_N=[min(x['o3hb_N'] for x in v), max(x['o3hb_N'] for x in v)],
                                     o2o3_N=[min(x['o2o3_N'] for x in v), max(x['o2o3_N'] for x in v)])
                        for k, v in sorted(byv.items())}
    out['V_tried'] = len(out['V'])
    print('L later truncation:', json.dumps(out['L_summary']))
    for k, v in out['V_summary'].items():
        print('V slower', k, json.dumps(v))
    # ---- Halpha patch vs complete shocks and truncated-run ends
    comp = []
    for r in fit['runs']:
        tag = r['tag']
        d = track(t, tag)
        e = -1
        rat = dict(O2_Hb=float(d['F_O2'][e] / d['F_Hb'][e]), O3_Hb=float(d['F_O3'][e] / d['F_Hb'][e]),
                   S2_Ha=float(d['F_S2'][e] / d['F_Ha'][e]), N2_Ha=float(d['F_N2'][e] / d['F_Ha'][e]),
                   Ha_Hb=float(d['F_Ha'][e] / d['F_Hb'][e]))
        cost = float(sum(np.log10(rat[k] / PATCH[k]) ** 2 for k in PATCH))
        comp.append(dict(tag=tag, v=r['v'], xn=r['xn'], b=r['b'], complete=r['complete'], T_end=float(d['T'][e]),
                         ratios=rat, log_cost=cost))
    comp.sort(key=lambda x: x['log_cost'])
    out['patch_obs'] = PATCH
    out['patch_best10'] = comp[:10]
    out['patch_complete_runs'] = [c for c in comp if c['complete']]
    print('patch best:', [(c['tag'], round(c['log_cost'], 3), {k: round(v, 2) for k, v in c['ratios'].items()})
                          for c in comp[:5]])
    for c in sorted(out['patch_complete_runs'], key=lambda x: (x['xn'], x['v'])):
        print('  complete', c['tag'], {k: round(v, 2) for k, v in c['ratios'].items()}, round(c['log_cost'], 2))
    json.dump(out, open(f'{DER}/w3_northsouth_shock.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
