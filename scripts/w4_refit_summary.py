"""w4_refit_summary: physical filters on the refitted shocks (per geometry), the north/south pairing (later truncation
vs slower shock), and the key-number collation for notes/w4_refit.md.
Inputs: data/derived/w4_refit_fit.json, w4_refit_geom.json, w4_refit_alt.json, w4_refit_tracks.npz
Filters (south, d = 282 pc, W = 0.28 deg, I_obs 1.62 R):
  curved limb : 1.2 <= I_pk <= 2.5 R; [O III] FWHM <= 38 km/s; FUV within 2 sigma of 80 +- 15 (+20 % model);
                t* <= 1.2e5 yr; cross-profile sigma_W/sigma_E within 0.98 +- 3 x 0.05 (never passed)
  flat sheet  : LOS depth L <= 83 pc (the sky length at 282 pc; 'plausible') or <= 166 pc ('allowed');
                FUV (edge-on or 1-Doppler tilt) within 2 sigma; t* <= 1.2e5 yr; sigma_W/sigma_E 0.83-1.13
  rippled     : FUV within 2 sigma; tangency width <= 0.28 deg
North: later truncation of each south-acceptable run (same v, X, b; N >= N*_S) vs north likelihood; slower shocks from
the north-acceptable set: brightness ratio at equal driving pressure (I_face ~ P eta / (mu m_H v)) and FUV/[O III] N/S.
Output: data/derived/w4_refit_summary.json
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_summary.py
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import w4_refit_common as C
import w4_refit_fit as WF
import verify_shockmodel_fit as FIT

ROOT = '/Users/mickey/astronomy/silk-nebula'
HC = 6.62607015e-27 * 2.99792458e10


def ok_fuv(x):
    s = np.hypot(15.0, 0.2 * x)
    return abs(x - 80.0) <= 2 * s


def rng(v):
    v = np.array([x for x in v if x is not None and np.isfinite(x)])
    if not v.size:
        return None
    return [float(v.min()), float(np.median(v)), float(v.max())]


if __name__ == '__main__':
    fit = json.load(open(ROOT + '/data/derived/w4_refit_fit.json'))
    geo = json.load(open(ROOT + '/data/derived/w4_refit_geom.json'))
    alt = json.load(open(ROOT + '/data/derived/w4_refit_alt.json'))
    T = np.load(ROOT + '/data/derived/w4_refit_tracks.npz')
    runs = {r['tag']: r for r in fit['runs']}
    M = geo['models']
    res = {}
    # ---------------- filters per geometry
    cur, fla, flp, rip = [], [], [], []
    for m in M:
        c, f = m['curved'], m['flat']
        pc = (1.2 <= c['I_pk'] <= 2.5) and c['line_fwhm'] <= 38 and ok_fuv(c['CIV']['FUV']) and c['t_star_yr'] <= 1.2e5
        m['pass_curved_noasym'] = bool(pc)
        m['pass_curved'] = bool(pc and abs(c['asym_WE'] - 0.98) <= 0.15)
        fuv_ok = ok_fuv(f['CIV_span0']['FUV']) or ok_fuv(f['CIV_span1']['FUV'])
        base = fuv_ok and f['t_star_yr'] <= 1.2e5 and 0.83 <= f['asym_WE'] <= 1.13
        m['pass_flat_allowed'] = bool(base and f['L_pc'] <= 166)
        m['pass_flat'] = bool(base and f['L_pc'] <= 83)
        r10 = m['ripple']
        m['pass_ripple'] = bool(any(ok_fuv(r10[f'R{R}_CIV']['FUV']) for R in (10, 20)))
        if m['pass_curved_noasym']:
            cur.append(m)
        if m['pass_flat_allowed']:
            fla.append(m)
        if m['pass_flat']:
            flp.append(m)
        if m['pass_ripple']:
            rip.append(m)
    def describe(S, geom):
        if not S:
            return {'n': 0}
        d = {'n': len(S), 'v': rng([m['v'] for m in S]), 'xn': rng([m['xn'] for m in S]), 'b': rng([m['b'] for m in S]),
             'logN': rng([m['logN'] for m in S]), 'T_star': rng([m['T_star'] for m in S]),
             'O3_Hb': rng([m['ratios']['O3_Hb'] for m in S]), 'O2_O3': rng([m['ratios']['O2_O3'] for m in S]),
             'O4363': rng([m['ratios']['O4363'] for m in S]), 'O2r_O3': rng([m['ratios']['O2r_O3'] for m in S]),
             'O2r_3727': rng([m['ratios']['O2r_3727'] for m in S]), 'HeII_O3': rng([m['ratios']['HeII_O3'] for m in S]),
             'O3_Ha': rng([m['ratios']['O3_Ha'] for m in S]), 'N2_Ha': rng([m['ratios']['N2_Ha'] for m in S]),
             'NUV': rng([m['uv']['NUV'] for m in S]), 'FUV_thin': rng([m['uv']['FUV_thin'] for m in S]),
             'FUV_nonres': rng([m['uv']['FUV_nonres'] for m in S]), 'SiIV': rng([m['uv_parts'].get('SiIV1400') for m in S]),
             'pm_mas_yr_282': rng([m['v'] / (4.74 * 282) * 1e3 for m in S])}
        g = m_ = None
        if geom == 'curved':
            for k in ('n0', 'P_k', 'I_pk', 'G', 'asym_WE', 'line_fwhm', 't_star_yr', 'offset_CIV_peak_deg'):
                d[k] = rng([m['curved'][k] for m in S])
            d['FUV'] = rng([m['curved']['CIV']['FUV'] for m in S]); d['S_eff'] = rng([m['curved']['CIV']['S_eff'] for m in S])
            d['G_CIV'] = rng([m['curved']['CIV']['G_CIV'] for m in S])
        if geom == 'flat':
            for k in ('n0', 'P_k', 'L_pc', 'L_over_len', 'G', 'sstar_pc', 't_star_yr', 't_cool_yr', 'asym_WE',
                      'offset_FUV_deg', 'offset_CIV_deg', 'offset_OII_deg'):
                d[k] = rng([m['flat'][k] for m in S])
            for sp in ('span0', 'span1'):
                d['FUV_' + sp] = rng([m['flat']['CIV_' + sp]['FUV'] for m in S])
                d['S_' + sp] = rng([m['flat']['CIV_' + sp]['S_eff'] for m in S])
                d['G_CIV_' + sp] = rng([m['flat']['CIV_' + sp]['G_CIV'] for m in S])
            d['tau1548'] = rng([m['flat']['CIV_span0']['tau1548'] for m in S])
            for inc in (0, 5, 10):
                d[f'line_fwhm_i{inc}'] = rng([m['flat'][f'line_i{inc}']['fwhm'] for m in S])
                d[f'vshift_i{inc}'] = rng([m['flat'][f'line_i{inc}']['v_shift'] for m in S])
            for dd in (252, 282, 318):
                d[f'P_k_d{dd}'] = rng([m['variants'][f'd{dd}_W0.28']['flat']['P_k'] for m in S])
                d[f'L_pc_d{dd}'] = rng([m['variants'][f'd{dd}_W0.28']['flat']['L_pc'] for m in S])
            for W in (0.15, 0.38):
                d[f'P_k_W{W}'] = rng([m['variants'][f'd282_W{W}']['flat']['P_k'] for m in S])
                d[f'L_pc_W{W}'] = rng([m['variants'][f'd282_W{W}']['flat']['L_pc'] for m in S])
        if geom == 'ripple':
            for k in ('R10_N1', 'R10_N2', 'R10_N4', 'R20_N1', 'R20_N2', 'R20_N4'):
                d['P_k_' + k] = rng([m['ripple'][k]['P_k'] for m in S])
                d['n0_' + k] = rng([m['ripple'][k]['n0'] for m in S])
                d['fwhm1_' + k] = rng([m['ripple'][k]['fwhm1_deg'] for m in S])
            for R in (10, 20):
                d[f'line_fwhm_R{R}'] = rng([m['ripple'][f'R{R}_line_fwhm'] for m in S])
                d[f'FUV_R{R}'] = rng([m['ripple'][f'R{R}_CIV']['FUV'] for m in S])
        return d
    res['n_spectral_acceptable'] = len(M)
    res['curved_without_asym'] = describe(cur, 'curved')
    res['curved_with_asym'] = {'n': sum(m['pass_curved'] for m in M),
                               'asym_all': rng([m['curved']['asym_WE'] for m in M])}
    res['flat_allowed_L166'] = describe(fla, 'flat')
    res['flat_plausible_L83'] = describe(flp, 'flat')
    res['ripple'] = describe(rip, 'ripple')
    res['all_spectral_flat'] = describe(M, 'flat')
    # why flat fails
    fails = {'L>166': 0, 'L>83': 0, 'FUV': 0, 't*': 0, 'asym': 0}
    for m in M:
        f = m['flat']
        fails['L>166'] += f['L_pc'] > 166; fails['L>83'] += f['L_pc'] > 83
        fails['FUV'] += not (ok_fuv(f['CIV_span0']['FUV']) or ok_fuv(f['CIV_span1']['FUV']))
        fails['t*'] += f['t_star_yr'] > 1.2e5; fails['asym'] += not (0.83 <= f['asym_WE'] <= 1.13)
    res['flat_fail_counts'] = {k: int(v) for k, v in fails.items()}
    fc = {'I': 0, 'width': 0, 'FUV': 0, 't*': 0, 'asym': 0}
    for m in M:
        c = m['curved']
        fc['I'] += not (1.2 <= c['I_pk'] <= 2.5); fc['width'] += c['line_fwhm'] > 38
        fc['FUV'] += not ok_fuv(c['CIV']['FUV']); fc['t*'] += c['t_star_yr'] > 1.2e5
        fc['asym'] += abs(c['asym_WE'] - 0.98) > 0.15
    res['curved_fail_counts'] = {k: int(v) for k, v in fc.items()}
    # ---------------- north
    VN = fit['variants']['N']
    mnN = min(VN['chi2_min_nat'], VN['chi2_min_R88'])
    south_tags = [m['tag'] for m in (flp or fla or M)]
    later = []
    for tag in south_tags:
        r = runs[tag]
        tr = WF.track(T, tag)
        iS = r['A_nat']['i']
        best = (1e9, None, None)
        for hb, vn in ((1.0, 'nat'), (FIT.HB_R88_FACTOR, 'R88')):
            c2, terms, rr, n = WF.score_track(tr, C.NORTH, hb)
            c2 = np.where(np.arange(c2.size) >= iS, c2, 1e9)
            j = int(np.argmin(c2))
            if c2[j] < best[0]:
                best = (float(c2[j]), j, vn)
        c2b, j, vn = best
        later.append(dict(tag=tag, chi2N=c2b, ok=bool(c2b <= mnN + 4), logN_N=float(np.log10(tr['N'][j])),
                          T_N=float(tr['T'][j]), O3_ratio_same_n0=float(tr['F_O3'][j] / tr['F_O3'][iS]),
                          O2_O3_N=float(tr['F_O2'][j] / tr['F_O3'][j])))
    res['north_later_truncation'] = dict(n_tested=len(later), n_ok=sum(x['ok'] for x in later),
                                         chi2N=rng([x['chi2N'] for x in later]),
                                         O3_ratio_same_n0_ok=rng([x['O3_ratio_same_n0'] for x in later if x['ok']]),
                                         O2_O3_N_ok=rng([x['O2_O3_N'] for x in later if x['ok']]), chi2N_min=mnN)
    # slower shocks: north-acceptable runs vs south flat survivors at equal pressure
    Sset = flp or fla
    def yield_per_P(tag, i):
        tr = WF.track(T, tag)
        v = runs[tag]['v']
        eta = (tr['F_O3'][i] / (HC / 5006.84e-8)) / (0.1 * v * 1e5)
        return eta / v, tr
    pairs = []
    Ntags = VN['acceptable']
    for tn in Ntags:
        rn = runs[tn]
        iN = rn['N_nat']['i'] if rn['N_nat']['chi2'] - VN['chi2_min_nat'] <= rn['N_R88']['chi2'] - VN['chi2_min_R88'] else rn['N_R88']['i']
        yN, trN = yield_per_P(tn, iN)
        uvN = WF.uv_at(trN, iN)
        for ms in Sset:
            ts = ms['tag']; iS = runs[ts]['A_nat']['i']
            yS, _ = yield_per_P(ts, iS)
            S = ms['flat']['CIV_span0']['S_eff']
            civN = uvN['FUV_parts'].get('CIV1549', 0.0)
            fuvN = uvN['FUV_thin'] - civN * (1 - S)
            fuvS = ms['flat']['CIV_span0']['FUV']
            pairs.append(dict(N=tn, S=ts, vN=rn['v'], vS=ms['v'], I_ratio_equalP=float(yN / yS),
                              PG_needed=[float(0.36 / (yN / yS)), float(0.45 / (yN / yS))],
                              FUV_NS=float(fuvN / fuvS), NUV_NS=float(uvN['NUV'] / ms['uv']['NUV'])))
    if pairs:
        res['north_slower'] = dict(
            n_north_acc=len(Ntags), north_v=rng([runs[t]['v'] for t in Ntags]), north_xn=rng([runs[t]['xn'] for t in Ntags]),
            north_T=rng([runs[t]['N_nat']['T'] for t in Ntags]), north_O3_Hb=rng([runs[t]['N_nat']['ratios']['O3_Hb'] for t in Ntags]),
            north_O2_O3=rng([runs[t]['N_nat']['ratios']['O2_O3'] for t in Ntags]),
            north_O2r_O3=rng([runs[t]['N_nat']['ratios']['O2r_O3'] for t in Ntags]),
            I_ratio_equalP=rng([p['I_ratio_equalP'] for p in pairs]),
            PG_needed=rng([np.mean(p['PG_needed']) for p in pairs]),
            FUV_NS=rng([p['FUV_NS'] for p in pairs]), NUV_NS=rng([p['NUV_NS'] for p in pairs]),
            frac_FUV_NS_within_1sig=float(np.mean([abs(p['FUV_NS'] - 0.61) <= 0.15 for p in pairs])),
            by_vN={str(v): dict(n=sum(1 for p in pairs if p['vN'] == v),
                                FUV_NS=rng([p['FUV_NS'] for p in pairs if p['vN'] == v]),
                                PG=rng([np.mean(p['PG_needed']) for p in pairs if p['vN'] == v]))
                   for v in sorted(set(p['vN'] for p in pairs))})
    # ---------------- alternatives compact table
    from scipy.stats import chi2 as CH
    alts = {}
    for k, f in alt['families'].items():
        b = f['best_A']; s = b['score']
        alts[k] = dict(label=b['label'], chi2A=s['A']['chi2'], n=s['A']['n'], pA=float(CH.sf(s['A']['chi2'], s['A']['n'])),
                       chi2B=s['B']['chi2'], chi2floor=f['best_floor']['score']['Afl']['chi2'],
                       pfloor=float(CH.sf(f['best_floor']['score']['Afl']['chi2'], s['A']['n'])), UV=s.get('UV'),
                       big_terms={t: round(x, 1) for t, x in s['A']['terms'].items() if x > 1.5},
                       ratios={t: b['ratios'].get(t) for t in ('O3_Hb', 'O2_O3', 'O4363', 'O2r_O3', 'HeII_O3')},
                       uv=b['uv'], extra=b['extra'])
    res['alternatives'] = alts
    VA = fit['variants']['A']
    res['shock_spectral'] = dict(chi2_min=min(VA['chi2_min_nat'], VA['chi2_min_R88']), nterms=VA['nterms'],
                                 n_acc=VA['n_acc'], n_trn=sum(1 for r in fit['runs'] if not r['complete']),
                                 ranges=VA['ranges'], chi2_by_v={k: min(x['nat'], x['R88']) for k, x in VA['chi2_by']['v'].items()},
                                 chi2_by_x={k: min(x['nat'], x['R88']) for k, x in VA['chi2_by']['xn'].items()},
                                 complete_min=min(min(c['nat'], c['R88']) for c in VA['complete'].values()),
                                 complete_min_B=min(min(c['nat'], c['R88']) for c in fit['variants']['B']['complete'].values()),
                                 complete_min_floor=min(min(c['nat'], c['R88']) for c in fit['variants']['Afl']['complete'].values()),
                                 doublet_common=fit['doublet_common_chi2'])
    json.dump(res, open(ROOT + '/data/derived/w4_refit_summary.json', 'w'), indent=1,
              default=lambda x: float(x) if np.isscalar(x) else None)
    for k, v in res.items():
        print('==', k, json.dumps(v, default=float)[:1500])
