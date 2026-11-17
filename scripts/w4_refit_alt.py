"""w4_refit_alt: re-score the verify_altmodel non-shock scenarios with the CORRECTED wave-3 inputs (w4_refit_common),
on exactly the same likelihood as the refitted shocks.  Model predictions are read from the verify_altmodel JSON
products (read-only; no altmodel code is re-run).  Each family is scanned over its own free parameter and its best
member is reported:
  warm CIE layer: T (74 values, 1e4 - 3e5 K)          conduction fronts: static (P/k 1.5e4, 2e6), evaporating (R_cl)
  NEI cooling / TML: GS07 isochoric/isobaric Z = 1, 2, start 1.6e6 or 3e5 K, stop 1-5e4 K (+ recombination tail)
  photoionisation: 6 spectra x log U (-5 ... -0.5)    photoionisation + heating: 3 spectra x T_fix
  fossil: 4 (T0, n_H) recombination tracks, every epoch
UV: FUV thin, FUV non-resonant floor (no C IV / Si IV / N V / Mg II), NUV; compared with FUV 80 +- 15, NUV 47 +- 15.
Also derived geometry-dependent budget numbers (flat sheet: G = I_obs / I_face is free).
Outputs: data/derived/w4_refit_alt.json
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_alt.py
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import w4_refit_common as C

ROOT = '/Users/mickey/astronomy/silk-nebula'
LHA, L5007 = 6564.6, 5008.2
HAHB_CASEB = 2.86


def from_alt(d):
    r = {}
    g = lambda k: d.get(k)
    r['O3_Hb'] = g('OIII_Hb'); r['O2_O3'] = g('OII_OIII'); r['O4363'] = g('O4363_5007')
    r['HeII_O3'] = g('HeII_OIII'); r['Ne3_O3'] = g('NeIII_OIII')
    if g('OII7325_OIII') is not None:
        r['O2r_O3'] = g('OII7325_OIII')
    elif g('OII_7325_3727') is not None and r['O2_O3'] is not None:
        r['O2r_O3'] = g('OII_7325_3727') * r['O2_O3']
    if g('OIII_Ha_phot') is not None:
        r['Ha_O3'] = 1.0 / (g('OIII_Ha_phot') * LHA / L5007)
    elif r['O3_Hb'] is not None and r['O3_Hb'] > 0:
        r['Ha_O3'] = HAHB_CASEB / r['O3_Hb']
    if g('NII_Ha') is not None and r.get('Ha_O3') is not None:
        r['N2_O3'] = g('NII_Ha') * r['Ha_O3']
    if g('SII_Ha') is not None and r.get('Ha_O3') is not None:
        r['S2_O3'] = g('SII_Ha') * r['Ha_O3']
    r = {k: v for k, v in r.items() if v is not None and np.isfinite(v) and (k != 'O3_Hb' or v > 0)}
    uv = {'FUV': g('FUV_CU_per_R') or g('FUV'), 'FUV_nonres': g('FUV_nonres_CU_per_R'), 'NUV': g('NUV_CU_per_R') or g('NUV')}
    parts = g('FUV_parts') or {}
    if uv['FUV_nonres'] is None and parts:
        uv['FUV_nonres'] = float(sum(v for k, v in parts.items() if k not in ('CIV1549', 'SiIV1400', 'NV1240')))
    uv['CIV'] = parts.get('CIV1549'); uv['SiIV'] = parts.get('SiIV1400')
    return r, uv


def score(r, uv=None):
    out = {}
    for lab, obs, o43 in (('A', C.SOUTH_A, None), ('B', C.SOUTH_B, None), ('Afl', C.SOUTH_A, C.O4363_FLOOR)):
        tot, terms, n = C.chi2(r, obs, o43)
        out[lab] = dict(chi2=float(tot), n=n, terms={k: float(v) for k, v in terms.items()})
    if uv:
        u = {}
        for k, (mu, sg) in (('FUV', C.UV_OBS['FUV']), ('NUV', C.UV_OBS['NUV'])):
            if uv.get(k):
                u[k] = float((uv[k] - mu)**2 / (sg**2 + (0.3 * uv[k])**2))
        if uv.get('FUV_nonres'):
            mu, sg = C.UV_OBS['FUV']
            x = uv['FUV_nonres']
            u['FUV_floor'] = float(max(x - mu, 0)**2 / (sg**2 + (0.3 * x)**2))   # one-sided: floor above the data
        out['UV'] = u
    return out


def best_of(cands):
    """cands: list of (label, ratios, uv, extra); returns best by chi2 A and the scan."""
    sc = []
    for lab, r, uv, ex in cands:
        s = score(r, uv)
        sc.append(dict(label=lab, ratios=r, uv=uv, extra=ex, score=s))
    ok = [c for c in sc if np.isfinite(c['score']['A']['chi2'])]
    bA = min(ok, key=lambda c: c['score']['A']['chi2'])
    bB = min(ok, key=lambda c: c['score']['B']['chi2'])
    bF = min(ok, key=lambda c: c['score']['Afl']['chi2'])
    return dict(best_A=bA, best_B=bB, best_floor=bF, n=len(ok),
                scan=[dict(label=c['label'], chiA=c['score']['A']['chi2'], chiB=c['score']['B']['chi2'],
                           chiF=c['score']['Afl']['chi2'], **{k: c['ratios'].get(k) for k in ('O3_Hb', 'O2_O3', 'O4363', 'O2r_O3')},
                           FUV=c['uv'].get('FUV') if c['uv'] else None, NUV=c['uv'].get('NUV') if c['uv'] else None)
                      for c in sc])


if __name__ == '__main__':
    cool = json.load(open(ROOT + '/data/derived/verify_altmodel_cooling.json'))
    phot = json.load(open(ROOT + '/data/derived/verify_altmodel_photoion.json'))
    bud = json.load(open(ROOT + '/data/derived/verify_altmodel_budget.json'))
    fam = {}
    # warm CIE layer (own CIE balance)
    c = []
    for d in cool['warm_CIE']['own']:
        r, uv = from_alt(d)
        c.append((f"T={d['T']:.0f}", r, uv, dict(T=d['T'], f5007=d.get('f5007'), NII_Ha=d.get('NII_Ha'), O2r_3727=d.get('OII_7325_3727'))))
    fam['warm_CIE'] = best_of(c)
    # conduction fronts
    c = []
    for k, d in cool['static_fronts'].items():
        r, uv = from_alt(d)
        c.append(('static ' + k, r, uv, dict(I_face=d.get('I5007_faceon_R'), thickness_pc=d.get('thickness_pc'), P_k=d.get('P_k'))))
    fam['conduction_static'] = best_of(c)
    c = []
    runs = cool['evaporation']['runs']
    for k, d in runs.items():
        r, uv = from_alt(d)
        c.append(('evap ' + k, r, uv, dict(I_face=d.get('I5007_faceon_R'), R_cloud_pc=d.get('R_cloud_pc'))))
    fam['conduction_evap'] = best_of(c)
    # NEI cooling flows / TML
    c = []
    for tk, tv in cool['cooling_flows'].items():
        for k, d in tv.items():
            if k.startswith('_'):
                continue
            r, uv = from_alt(d)
            c.append((f'{tk} {k}', r, uv, dict(energy_flux=d.get('energy_flux'), f5007=d.get('f5007'))))
    fam['NEI_cooling_TML'] = best_of(c)
    fam['NEI_from_1.6MK'] = best_of([x for x in c if '1.6e+06' in x[0]])
    fam['NEI_from_3e5K'] = best_of([x for x in c if '3.0e+05' in x[0]])
    # photoionisation equilibrium (free U)
    for sp, rows in phot['models'].items():
        c = []
        for d in rows:
            r, uv = from_alt(d)
            c.append((f"{sp} logU={d['logU']}", r, None, dict(T=d['T'], logU=d['logU'])))
        fam['PIE_' + sp] = best_of(c)
    # PIE + extra heating
    for sp, rows in bud['PIE_plus_heating'].items():
        c = []
        for d in rows:
            r, uv = from_alt(d)
            c.append((f"{sp} Tfix={d['Tfix']}", r, uv, dict(T=d['T'], heat=d.get('heat_extra_over_photo'))))
        fam['PIEheat_' + sp] = best_of(c)
    # fossil
    for k, fv in bud['fossil'].items():
        c = []
        for d in fv['rows']:
            r, uv = from_alt(d)
            c.append((f"{k} t={d['t_yr']:.3g}yr", r, None, dict(t_yr=d['t_yr'], T=d['T'])))
        fam['fossil_' + k] = best_of(c)
        fam['fossil_' + k]['t_fade_yr'] = fv.get('t_fade_yr')
    # geometry-free budget numbers (derived): flat sheet, G free
    E5007 = 6.62607015e-27 * 2.99792458e10 / (5006.84e-8)
    budget = {}
    for G in (10.0, 20.0, 50.0):
        Iface = 1.62 / G
        Fph = Iface * 1e6
        budget[f'G{int(G)}'] = dict(I_face_R=Iface, face_on_5007_erg=Fph * E5007,
                                    warm_heating_flux=Fph * E5007 / 0.0252,
                                    static_front_Pk_for_width=2.0e6 * 0.646 / 1.38,
                                    static_front_Iface_at_that_P=0.3705 * (0.646 / 1.38),
                                    static_front_G_needed=1.62 / (0.3705 * 0.646 / 1.38),
                                    TML_inflow_kms_at_Monogem_P=377.3 * 18.73 / G)
    out = {'families': fam, 'budget_flat_sheet': budget, 'inputs': {'SOUTH_A': C.SOUTH_A, 'UV_OBS': C.UV_OBS}}
    json.dump(out, open(ROOT + '/data/derived/w4_refit_alt.json', 'w'), indent=1,
              default=lambda x: float(x) if np.isscalar(x) else None)
    for k, v in fam.items():
        for lab in ('best_A', 'best_floor'):
            b = v[lab]
            s = b['score']
            print(f"{k:28s} {lab:10s} {b['label']:40s} chiA={s['A']['chi2']:.1f}/{s['A']['n']} chiB={s['B']['chi2']:.1f} "
                  f"chiFl={s['Afl']['chi2']:.1f} UV={s.get('UV')} O3Hb={b['ratios'].get('O3_Hb', np.nan):.1f} "
                  f"O2O3={b['ratios'].get('O2_O3', np.nan):.2f} 4363={b['ratios'].get('O4363', np.nan):.3f} "
                  f"7325={b['ratios'].get('O2r_O3', np.nan):.3f} terms={ {t: round(x, 1) for t, x in s['A']['terms'].items() if x > 1} }")
    print(json.dumps(budget, indent=1))
