"""verify_shockmodel_r88: validate the shockmodel 1-D code against published INCOMPLETE (column-truncated) shocks.

Reference: Raymond, Hester, Cox, Blair, Fesen & Gull 1988, ApJ 324, 869 (R88), Table 5 'Model' columns and Table 7
(read by me from the ADS scan data/raw/verify_xrayfuv_lit/1988ApJ324869R.pdf, pages 17 and 19; the same numbers were
independently hand-transcribed by the xrayfuv agent in scripts/verify_xrayfuv_uv_predict.py -- they agree).
  SE: v_s = 140 km/s, n0 = 4.4 cm^-3, B0 = 11 uG, log N_H(trunc) = 18.20, preshock H I 0, He I 0, He II 0.88
  M : v_s = 120 km/s, n0 = 6.0 cm^-3, B0 = 15 uG, log N_H(trunc) = 18.25, same preionization
  Abundances (R88 sec. III a): He 10.93, C 8.52, N 7.96, O 8.82, Ne 7.92, Mg 7.20, Si 7.50, S 7.20, Fe 7.50.
Also complete shocks: Raymond 1979 (R79) models E-I (n0 = 10, complete; values as transcribed by the xrayfuv agent
from Tables 1/3/4) for [O III]5007/Hb, [O II]/Hb, C IV/Hb.
Output: data/derived/verify_shockmodel_r88.json
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_atomic as A
import verify_shockmodel_flow as F

ROOT = '/Users/mickey/astronomy/silk-nebula'

A.ABUND['r88'] = {'He': 10**(10.93-12), 'C': 10**(8.52-12), 'N': 10**(7.96-12), 'O': 10**(8.82-12),
                  'Ne': 10**(7.92-12), 'Mg': 10**(7.20-12), 'Si': 10**(7.50-12), 'S': 10**(7.20-12),
                  'Fe': 10**(7.50-12)}

# R88 Table 5 model columns (energy, Hb = 100); F(Hb) face-on in 1e-5 erg cm^-2 s^-1 sr^-1
R88 = {
    'SE': dict(v=140, n0=4.4, B=11.0, logN=18.20, lines={
        'CIV1550': 6800, 'HeII1640': 275, 'OIII1664': 537, 'OII3727': 1730, 'NeIII3869': 100, 'OIII4363': 91,
        'Hb': 100, 'OIII4959': 433, 'OIII5007': 1290, 'NII5755': 10, 'HeI5876': 9, 'Ha': 303, 'NII6584': 410,
        'SII6717': 104, 'SII6731': 77, 'OII7319': 32, 'OII7330': 26}, FHb=0.092),
    'M': dict(v=120, n0=6.0, B=15.0, logN=18.25, lines={
        'CIV1550': 4330, 'HeII1640': 166, 'OIII1664': 427, 'OII3727': 1720, 'NeIII3869': 96, 'OIII4363': 74,
        'Hb': 100, 'OIII4959': 353, 'OIII5007': 1050, 'NII5755': 10, 'HeI5876': 10, 'Ha': 304, 'NII6584': 393,
        'SII6717': 116, 'SII6731': 86, 'OII7319': 32, 'OII7330': 26}, FHb=0.107),
}
PRE_R88 = ({'H': {1: 1.0}, 'He': {1: 0.88, 2: 0.12}, 'C': {2: 0.6, 3: 0.4}, 'N': {1: 0.3, 2: 0.6, 3: 0.1},
            'O': {1: 0.5, 2: 0.5}, 'Ne': {1: 0.5, 2: 0.5}, 'Mg': {2: 1.0}, 'Si': {2: 0.5, 3: 0.5},
            'S': {2: 0.6, 3: 0.4}}, 1.0e4)
PRE_R88_LOWO = ({**PRE_R88[0], 'O': {1: 0.9, 2: 0.1}, 'N': {1: 0.8, 2: 0.2}, 'C': {1: 0.3, 2: 0.7}}, 1.0e4)

# Raymond 1979 complete models (n0 = 10), energy relative to Hb = 100 (5007 = 0.749 x doublet)
R79 = {100: dict(O3=547 * 0.749, O2=582, C4=756), 120: dict(O3=572 * 0.749, O2=593, C4=1610),
       141: dict(O3=596 * 0.749, O2=607, C4=2580), 160: dict(O3=495 * 0.749, O2=578, C4=1080),
       200: dict(O3=373 * 0.749, O2=519, C4=612)}


def model_lines(Fc, i):
    hb = Fc['Hb'][i]
    g = lambda *ks: sum(Fc[k][i] for k in ks)
    out = {'CIV1550': g('CIV1548', 'CIV1551'), 'HeII1640': g('HeII1640'), 'OIII1664': g('OIII1666', 'OIII1661'),
           'OII3727': g('OII3726', 'OII3729'), 'NeIII3869': g('NeIII3869'), 'OIII4363': g('OIII4363'), 'Hb': hb,
           'OIII4959': g('OIII4959'), 'OIII5007': g('OIII5007'), 'NII5755': g('NII5755'), 'Ha': g('Ha'),
           'NII6584': g('NII6583'), 'SII6717': g('SII6716'), 'SII6731': g('SII6731'),
           'OII7319': g('OII7320', 'OII7319b'), 'OII7330': g('OII7330', 'OII7330b')}
    return {k: 100 * v / hb for k, v in out.items()}, hb


def run_r88(key, pre=PRE_R88, logN=None):
    p = R88[key]
    net = F.Network('r88')
    sh = F.Shock(net, p['v'], p['n0'], B0_uG=p['B'], pre=pre, complete=True, Tfloor=1e4, Nmax=10**19.5).run()
    j = F.emissivities(sh, thin=1)
    Fc = F.cumulative(sh, j)
    res = {}
    for lab, lN in (('at_R88_N', p['logN'] if logN is None else logN), ('N-0.15dex', p['logN'] - 0.15),
                    ('N+0.15dex', p['logN'] + 0.15)):
        i = int(np.argmin(np.abs(np.log10(np.maximum(sh.N, 1)) - lN)))
        r, hb = model_lines(Fc, i)
        # face-on surface brightness: flux / (4 pi) per sr, in 1e-5 erg cm^-2 s^-1 sr^-1
        res[lab] = dict(logN=float(np.log10(sh.N[i])), T=float(sh.T[i]), ratios={k: round(float(v), 1) for k, v in r.items()},
                        FHb_faceon=float(hb / (4 * np.pi) / 1e-5))
    # [O III]/Hb running maximum and where it occurs
    with np.errstate(divide='ignore', invalid='ignore'):
        rr = Fc['OIII5007'] / Fc['Hb']
    k = int(np.nanargmax(np.where(sh.N > 1e15, rr, np.nan)))
    res['O3Hb_max'] = dict(value=float(rr[k]), logN=float(np.log10(sh.N[k])), T=float(sh.T[k]))
    res['Ts'] = float(sh.Ts)
    return res, sh, Fc


if __name__ == '__main__':
    out = {'R88': {}, 'R79': {}}
    for key in ('SE', 'M'):
        for prelab, pre in (('preO+O++', PRE_R88), ('preO+', PRE_R88_LOWO)):
            res, sh, Fc = run_r88(key, pre)
            out['R88'][f'{key}_{prelab}'] = {'model': res, 'R88_table5': R88[key]['lines'], 'R88_FHb': R88[key]['FHb']}
            m = res['at_R88_N']
            print(f"R88 {key} {prelab}: Ts={res['Ts']:.3g} at logN={m['logN']:.2f} T={m['T']:.3g}  FHb(face)={m['FHb_faceon']:.3f} (R88 {R88[key]['FHb']})")
            for ln in ('OIII5007', 'OII3727', 'OIII4363', 'Ha', 'NII6584', 'SII6717', 'NeIII3869', 'HeII1640',
                       'CIV1550', 'OIII1664', 'NII5755', 'OII7319'):
                print(f"    {ln:10s} mine {m['ratios'][ln]:8.1f}  R88 {R88[key]['lines'][ln]:6d}   "
                      f"(-0.15dex {res['N-0.15dex']['ratios'][ln]:8.1f}, +0.15dex {res['N+0.15dex']['ratios'][ln]:8.1f})")
            print('    O3/Hb max', res['O3Hb_max'])
    # complete shocks vs R79 (n0 = 10, preionized, weak B)
    net = F.Network('asplund09')
    for v in R79:
        sh = F.Shock(net, v, 10.0, B0_uG=1.0, pre='precursor', complete=True, Tfloor=7e3).run()
        Fc = F.cumulative(sh, F.emissivities(sh, thin=2))
        r, hb = model_lines(Fc, -1)
        out['R79'][v] = {'mine': {'O3': r['OIII5007'], 'O2': r['OII3727'], 'C4': r['CIV1550']}, 'R79': R79[v]}
        print(f"R79 v={v}: O3 {r['OIII5007']:.0f} vs {R79[v]['O3']:.0f}; O2 {r['OII3727']:.0f} vs {R79[v]['O2']}; "
              f"C4 {r['CIV1550']:.0f} vs {R79[v]['C4']}")
    # R88 Fig. 8a: 80 km/s, n0 = 8.6, B = 18 uG, preshock H I 0.10, He I 0.2 (Table 7 NW column).  Values read by
    # eye from the scan (log ratio to Hb: [O III]4959+5007 peak ~1.07 at log N ~17.45, [O II] peak ~1.50 at ~17.65,
    # O III]1663 peak ~0.60 at ~17.35; cumulative [O III] 50 % at ~17.05, O III] 50 % at ~16.95; end (log N 18.6):
    # [O III] ~0.48, [O II] ~1.03, O III] ~-0.05).  NB the R88 [O III] curve is the 4959+5007 doublet.
    pre80 = ({'H': {0: 0.10, 1: 0.90}, 'He': {0: 0.2, 1: 0.8}, 'C': {1: 0.5, 2: 0.5}, 'N': {0: 0.1, 1: 0.8, 2: 0.1},
              'O': {0: 0.1, 1: 0.8, 2: 0.1}, 'Ne': {0: 0.3, 1: 0.7}, 'Mg': {1: 0.3, 2: 0.7}, 'Si': {1: 0.3, 2: 0.7},
              'S': {1: 0.6, 2: 0.4}}, 1e4)
    net = F.Network('r88')
    sh = F.Shock(net, 80, 8.6, B0_uG=18.0, pre=pre80, complete=True, Tfloor=8e3, Nmax=10**19).run()
    Fc = F.cumulative(sh, F.emissivities(sh, thin=1))
    lN = np.log10(np.maximum(sh.N, 1))
    o3 = Fc['OIII5007'] + Fc['OIII4959']; o2 = Fc['OII3726'] + Fc['OII3729']; s = Fc['OIII1666'] + Fc['OIII1661']
    with np.errstate(all='ignore'):
        cur = {'O3': o3 / Fc['Hb'], 'O2': o2 / Fc['Hb'], 'O3s': s / Fc['Hb']}
    fig8 = {'O3': (1.07, 17.45, 0.48), 'O2': (1.50, 17.65, 1.03), 'O3s': (0.60, 17.35, -0.05)}
    out['R88_fig8a_80kms'] = {}
    for k, (pk, pN, end) in fig8.items():
        m = np.where(lN > 15.5, cur[k], np.nan)
        i = int(np.nanargmax(m))
        out['R88_fig8a_80kms'][k] = dict(mine_peak_log=float(np.log10(cur[k][i])), mine_peak_logN=float(lN[i]),
                                         mine_end_log=float(np.log10(cur[k][-1])), r88_peak_log=pk, r88_peak_logN=pN,
                                         r88_end_log=end)
        print('R88 Fig8a', k, out['R88_fig8a_80kms'][k])
    c = o3 / o3[-1]
    out['R88_fig8a_80kms']['O3_cum50_logN'] = dict(mine=float(lN[np.argmax(c >= 0.5)]), r88=17.05)
    c = s / s[-1]
    out['R88_fig8a_80kms']['O3s_cum50_logN'] = dict(mine=float(lN[np.argmax(c >= 0.5)]), r88=16.95)
    print(out['R88_fig8a_80kms'])
    json.dump(out, open(ROOT + '/data/derived/verify_shockmodel_r88.json', 'w'), indent=1)
