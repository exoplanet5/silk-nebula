"""verify_shockmodel_summary: condense verify_shockmodel_fit/analysis into verdict-level numbers.
Output: data/derived/verify_shockmodel_summary.json"""
import json, numpy as np
ROOT = '/Users/mickey/astronomy/silk-nebula'
fit = json.load(open(ROOT + '/data/derived/verify_shockmodel_fit.json'))
an = json.load(open(ROOT + '/data/derived/verify_shockmodel_analysis.json'))
runs = [r for r in fit['runs'] if not r['complete']]
out = {}
# chi2 by parameter (min over the others), native and R88-calibrated Hb
for par in ('v', 'xn', 'b'):
    vals = sorted(set(r[par] for r in runs))
    out[f'chi2_by_{par}'] = {str(v): {'native': min(r['best']['chi2'] for r in runs if r[par] == v),
                                      'R88Hb': min(r['best_hbR88']['chi2'] for r in runs if r[par] == v)} for v in vals}
# complete-shock endpoints
cmp_ = [r for r in fit['runs'] if r['complete']]
out['complete_end'] = {r['tag']: r['eta_end'] for r in cmp_}
# fully consistent set
cons = []
for o in an['acceptable']:
    W = o['W']; fuv = o['FUV_RT']['total_CU_per_R']
    ok = dict(P=3e4 <= W['P_ram_k'] <= 4e5, I=1.0 <= W['I_pk'] <= 4.0, width=W['line_fwhm'] <= 38.0,
              age=W['t_star_yr'] <= 1.1e5, fuv=abs(fuv - 92) <= 2 * np.hypot(10, 0.3 * fuv))
    o['consistency'] = ok
    if all(ok.values()):
        cons.append(o)
out['n_acceptable'] = len(an['acceptable'])
out['n_consistent'] = len(cons)
keys = ['v', 'xn', 'b', 'chi2', 'logN', 'T_star', 'eta5007']
out['consistent'] = [{**{k: o[k] for k in keys}, 'n0': o['W']['n0'], 'B0_uG': o['b'] * np.sqrt(o['W']['n0']),
                      'I_pk': o['W']['I_pk'], 'G': o['W']['G'], 'sstar_pc': o['W']['sstar_pc'], 'P_k': o['W']['P_ram_k'],
                      't_star': o['W']['t_star_yr'], 't_cool': o['W']['t_cool_yr'], 'line_fwhm': o['W']['line_fwhm'],
                      'vlos_max': o['W']['vlos_max'], 'ne_OII': o['W']['ne_OII'], 'OII_doublet': o['W']['OII_doublet'],
                      'I_Ha': o['W']['I_Ha_R'], 'FUV_thin': o['FUV_RT']['total_thin_CU_per_R'],
                      'FUV_RT': o['FUV_RT']['total_CU_per_R'], 'S1548': o['FUV_RT']['S1548'],
                      'tau1548': o['FUV_RT']['tau1548_face'], 'NCIV': o['FUV_RT']['NCIV'], 'pm': o['pm_mas_yr'],
                      'I_pk_d250': o['d250']['I_pk'], 'I_pk_d320': o['d320']['I_pk'], 'P_d250': o['d250']['P_ram_k'],
                      'P_d320': o['d320']['P_ram_k'], 'I_Rc15': o['Rc15']['W']['I_pk'], 'P_Rc15': o['Rc15']['W']['P_ram_k'],
                      'ratios': {k: o['ratios'][k] for k in ('O3_Hb', 'O2_O3', 'O3_Ha', 'O4363', 'O2r_O3', 'N2T', 'HeII_O3',
                                                             'CIV_O3', 'N2_Ha', 'S2_Ha', 'Ha_Hb', 'R23')}} for o in cons]
def rng(key, sub=None):
    a = np.array([c[key] if sub is None else c[sub][key] for c in out['consistent']])
    return [float(a.min()), float(np.median(a)), float(a.max())] if a.size else None
out['consistent_ranges'] = {k: rng(k) for k in ('v', 'xn', 'b', 'logN', 'T_star', 'eta5007', 'n0', 'B0_uG', 'I_pk', 'G',
                                                'sstar_pc', 'P_k', 't_star', 't_cool', 'line_fwhm', 'ne_OII', 'OII_doublet',
                                                'I_Ha', 'FUV_thin', 'FUV_RT', 'S1548', 'tau1548', 'NCIV', 'pm', 'chi2',
                                                'I_pk_d250', 'I_pk_d320', 'I_Rc15', 'P_Rc15')}
out['consistent_ratio_ranges'] = {k: rng(k, 'ratios') for k in out['consistent'][0]['ratios']} if cons else None
# which consistency tests fail most often among acceptable
from collections import Counter
fails = Counter(k for o in an['acceptable'] for k, v in o['consistency'].items() if not v)
out['fail_counts'] = dict(fails)
json.dump(out, open(ROOT + '/data/derived/verify_shockmodel_summary.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('chi2_by_v', 'chi2_by_xn', 'chi2_by_b', 'n_acceptable', 'n_consistent', 'fail_counts')}, indent=0))
for c in out['consistent']:
    print(f"v{c['v']} x{c['xn']} b{c['b']} chi2 {c['chi2']:.1f} logN {c['logN']:.2f} T* {c['T_star']:.3g} n0 {c['n0']:.3f} B {c['B0_uG']:.2f} I {c['I_pk']:.2f} P {c['P_k']:.2e} t* {c['t_star']:.2e} tc {c['t_cool']:.2e} FWHMv {c['line_fwhm']:.0f} FUV {c['FUV_thin']:.0f}->{c['FUV_RT']:.0f} S {c['S1548']:.2f} tau {c['tau1548']:.2f} pm {c['pm']:.0f} Ha {c['I_Ha']:.2f} d250/320 I {c['I_pk_d250']:.2f}/{c['I_pk_d320']:.2f} Rc15 I {c['I_Rc15']:.2f} P {c['P_Rc15']:.1e}")
print(json.dumps(out['consistent_ranges'], indent=0))
print(json.dumps(out['consistent_ratio_ranges'], indent=0))
