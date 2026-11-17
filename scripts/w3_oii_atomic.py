"""w3_oii_atomic.py -- (wave 3, tag oii) density / temperature diagnostics and the OH-aware [O II]7325 re-fit.

1. 5-level O+ (PyNeb data files on disk: Kisielius+09 or Tayal 07 collision strengths x Froese Fischer+04 or
   Zeippen 82/Wiese+96 A-values) and S+ (Tayal & Zatsarinny 10 + Irimia & Froese Fischer 05):
   [O II] 3729/3726 and [S II] 6716/6731 versus n_e at T = 0.8-2.0e4 K; inversion of the measured ratios.
2. [O II] (7319+7320+7330+7331)/(3726+3729) versus T at n_e = 1, 30, 100, 300 cm^-3 (density-temperature degeneracy),
   inversion of the measured ratio.
3. DESI [O II]7325 re-fit in the regression spectra (w3_oii_desi_regspec.npz, geometry instrument, F weighting) with
   and without free OH 'residual' Gaussians at the 7303.7, 7316.3, 7329.1, 7340.9, 7358.7 A (air) OH lines, which
   the filament-correlated spectrum shows at comparable amplitude. [O II] internal ratios fixed from the atom.
4. Implications of n_e ~ 100: thermal pressure, emitting path length for the observed [O II] surface brightness,
   cooling and recombination times, filling factor, Halpha expected from the same layer.
Output: data/derived/w3_oii_atomic.json
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_oii_common import (L, C, air2vac, oii_atoms, sii_atom, oii_ratio, oii_7325_3727, sii_ratio, invert,  # noqa
                           design, grid_fit)

D = '/Users/mickey/astronomy/silk-nebula/data/derived/'
TS = [8000.0, 10000.0, 12500.0, 15000.0, 20000.0]
NE = np.logspace(-1, 4, 251)
KB = 1.380649e-16


def curves(ats, s2):
    out = {'ne_grid': NE.tolist()}
    for k, a in ats.items():
        out[f'OII_{k}'] = {f'{T:.0f}': [float(oii_ratio(a, T, n)) for n in NE] for T in TS}
    out['SII_TZ10+IFF05'] = {f'{T:.0f}': [float(sii_ratio(s2, T, n)) for n in NE] for T in TS}
    return out


def ne_of(ats, R, T):
    return {k: invert(lambda n, a=a: oii_ratio(a, T, n), R) for k, a in ats.items()}


def main():
    ats = oii_atoms()
    s2 = sii_atom()
    res = {'curves': curves(ats, s2)}
    res['low_density_limits'] = {k: {f'{T:.0f}': float(oii_ratio(a, T, 0.01)) for T in TS} for k, a in ats.items()}
    res['low_density_limits']['SII'] = {f'{T:.0f}': float(sii_ratio(s2, T, 0.01)) for T in TS}
    print('LDL', json.dumps(res['low_density_limits']))
    # --------------------------------------------------- inversion of the measured doublet ratios
    fits = json.load(open(D + 'w3_oii_desi_fits.json'))
    sd = json.load(open(D + 'w3_oii_sdss.json'))
    wim = json.load(open(D + 'w3_oii_sdss_wim.json'))
    meas = {}
    for sn in ('all', 'south', 'north', 'Zpos', 'Zneg', 'qso_hiz', 'lowz_gal', 'dark', 'bright'):
        k = f'geo|{sn}|F'
        if k in fits:
            o = fits[k]['o2']['kin']
            half = 0.5 * (o['p84'] - o['p16'])
            meas[f'DESI_{sn}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'], sig_boot=half,
                                      sig_inflated_x2=2 * half)
    for rn in ('all', 'south', 'north', 'patch30', 'patch_halfmax'):
        for t in ('fixed', 'vsfree'):
            o = sd['onoff'][rn]['o2'][t]
            meas[f'SDSS_{rn}_{t}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'], sig_boot=0.5 * (o['p84'] - o['p16']))
    for mode in ('between', 'within'):
        for t in ('fixed', 'vsfree'):
            o = wim[mode]['o2'][t]
            meas[f'SDSS_WIM_{mode}_{t}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'])
    inv = {}
    for nm, m in meas.items():
        r = {}
        for T in (10000.0, 12500.0, 20000.0):
            r[f'T{T:.0f}'] = {q: ne_of(ats, m[q], T) for q in ('R', 'p16', 'p84') if np.isfinite(m[q])}
        inv[nm] = r
    res['measured'] = meas
    res['ne_inversion'] = inv
    for nm in ('DESI_all', 'DESI_south', 'DESI_north'):
        q = inv[nm]['T12500']['R']
        print(nm, 'n_e over atomic datasets at 1.25e4 K:', {k: round(v, 0) for k, v in q.items()})
    for nm in ('DESI_all', 'DESI_south', 'DESI_north', 'SDSS_south_vsfree'):
        if nm in inv:
            q = inv[nm]['T12500']
            print(nm, 'R=%.3f [%.3f,%.3f] -> n_e(T=1.25e4, Kal09+FFT04) = %.0f [%s, %s]' % (
                meas[nm]['R'], meas[nm]['p16'], meas[nm]['p84'], q['R']['Kal09+FFT04'],
                f"{q['p84']['Kal09+FFT04']:.0f}", f"{q['p16']['Kal09+FFT04']:.0f}"))
    # DESI all with x2-inflated 68 % interval and 2-sigma upper limit on n_e
    m = meas['DESI_all']
    lo2 = m['R'] - 2 * m['sig_inflated_x2']
    res['DESI_all_ne_2sigma_upper'] = {k: invert(lambda n, a=a: oii_ratio(a, 12500.0, n), lo2) for k, a in ats.items()}
    res['DESI_all_R_2sigma_low_inflated'] = lo2
    # --------------------------------------------------- [S II] inversions
    sii_meas = {}
    for rn in ('all', 'south', 'north', 'patch30', 'patch_halfmax'):
        o = sd['onoff'][rn]['s2']['fixed']
        sii_meas[f'SDSS_{rn}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'])
    for mode in ('between', 'within'):
        o = wim[mode]['s2']['fixed']
        sii_meas[f'SDSS_WIM_{mode}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'])
    for sn in ('all', 'south', 'north'):
        k = f'geo|{sn}|F'
        if k in fits and 's2' in fits[k]:
            o = fits[k]['s2']['kin']
            sii_meas[f'DESI_geo_{sn}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'])
        k = f'ha|{sn}|F'
        if k in fits and 's2' in fits[k]:
            o = fits[k]['s2']['kin']
            sii_meas[f'DESI_Ha_{sn}'] = dict(R=o['ratio'], p16=o['p16'], p84=o['p84'])
    for nm, m in sii_meas.items():
        m['ne_T1e4'] = {q: invert(lambda n: sii_ratio(s2, 1e4, n), m[q]) for q in ('R', 'p16', 'p84')}
        print('[S II]', nm, 'R=%.2f [%.2f,%.2f] n_e(1e4K) %s' % (m['R'], m['p16'], m['p84'],
                                                                {q: round(v, 0) for q, v in m['ne_T1e4'].items()}))
    res['sii'] = sii_meas
    # --------------------------------------------------- 7325/3727 versus T and n_e
    T_grid = np.arange(7000.0, 30001.0, 250.0)
    tab = {}
    for k, a in ats.items():
        tab[k] = {f'ne{n:g}': [float(oii_7325_3727(a, T, n)) for T in T_grid] for n in (1.0, 30.0, 100.0, 300.0)}
    res['r7325_3727_curves'] = dict(T=T_grid.tolist(), **tab)

    def T_of(r, k, n):
        y = np.array(tab[k][f'ne{n:g}'])
        return float(np.interp(r, y, T_grid)) if y.min() < r < y.max() else np.nan
    res['T_OII_for_ratio'] = {f'{r:.3f}': {f'{k}|ne{n:g}': T_of(r, k, n) for k in ats for n in (1.0, 100.0, 300.0)}
                              for r in (0.020, 0.026, 0.035, 0.045, 0.055)}
    print('T(O II) [Kal09 ne=1, ne=100; T07 ne=1]:', {r: (v['Kal09+FFT04|ne1'], v['Kal09+FFT04|ne100'], v['T07+FFT04|ne1'])
                                                       for r, v in res['T_OII_for_ratio'].items()})
    # --------------------------------------------------- DESI 7325 re-fit with OH residual terms
    Z = np.load(D + 'w3_oii_desi_regspec.npz')
    LS = json.load(open(D + 'w3_oii_desi_lsf.json'))
    sR = 0.5 * (LS['7321.0']['median'] + LS['7331.4']['median'])
    oh = list(air2vac(np.array([7303.716, 7316.282, 7329.148, 7340.885, 7358.659])))
    a0 = ats['Kal09+FFT04']
    em = {(j, i): a0.emis(12500.0, 30.0, j, i) for j in (4, 5) for i in (2, 3)}
    # 2P3/2 (4), 2P1/2 (5) -> 2D5/2 (2), 2D3/2 (3): wavelengths from the level energies
    comp = [(1e8 / (a0.E[j - 1] - a0.E[i - 1]), em[(j, i)]) for (j, i) in em]
    tot = sum(e for _, e in comp)
    out7 = {}
    for sn in ('all', 'south', 'north'):
        k = f'S_geo_{sn}_7325_F'
        if k not in Z.files:
            continue
        lam = Z['lam_7325']
        S, bs = Z[k], Z[f'bs_geo_{sn}_7325_F'].astype(float)
        sel = (lam > 7295) & (lam < 7365)
        x = lam[sel]
        Ys = np.vstack([S[sel][None], bs[:, sel]])
        W = np.eye(x.size)
        v = fits[f'geo|{sn}|F']['o3']['v'] if f'geo|{sn}|F' in fits else 11.0
        sint = fits[f'geo|{sn}|F']['o3']['sigma_int_A'] * 7325 / 5008.2 if f'geo|{sn}|F' in fits else 0.2
        s = float(np.sqrt(sR ** 2 + sint ** 2))
        g_o2 = sum(e / tot * np.exp(-0.5 * ((x - lw * (1 + v / C)) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
                   for lw, e in comp)
        s_oh = float(np.sqrt(sR ** 2 + (7330 * 20 / C) ** 2))
        g_oh = [np.exp(-0.5 * ((x - l0) / s_oh) ** 2) / (np.sqrt(2 * np.pi) * s_oh) for l0 in oh]
        base = [((x - 7330) / 35.0) ** p for p in range(2)]
        r = {}
        for tag, cols in (('no_OH', [g_o2] + base), ('OH_free', [g_o2] + g_oh + base),
                          ('OH_free_excl7329', [g_o2] + [g for g, l0 in zip(g_oh, oh) if abs(l0 - 7331.2) > 1] + base)):
            X = np.vstack(cols)
            c2, par, p = grid_fit(Ys, W, [((v, s), X)])
            amp = np.array([q[0] for q in p])
            r[tag] = dict(OII7325_over_OIII=float(amp[0]), p16=float(np.percentile(amp[1:], 16)),
                          p84=float(np.percentile(amp[1:], 84)),
                          OH_amps=[float(q) for q in p[0][1:1 + len(cols) - 3]] if 'OH' in tag else None)
        o2tot = fits[f'geo|{sn}|F']['o2']['kin']['sum_amp']
        for tag in r:
            r[tag]['ratio_7325_3727'] = r[tag]['OII7325_over_OIII'] / o2tot
            r[tag]['ratio_p16'] = r[tag]['p16'] / o2tot
            r[tag]['ratio_p84'] = r[tag]['p84'] / o2tot
            r[tag]['T_OII_ne1_Kal09'] = T_of(r[tag]['ratio_7325_3727'], 'Kal09+FFT04', 1.0)
            r[tag]['T_OII_ne100_Kal09'] = T_of(r[tag]['ratio_7325_3727'], 'Kal09+FFT04', 100.0)
        r['OII3727_over_OIII_kin'] = o2tot
        out7[sn] = r
        print('7325', sn, {t: (round(q['OII7325_over_OIII'], 4), round(q['ratio_7325_3727'], 4),
                               round(q['ratio_p16'], 4), round(q['ratio_p84'], 4), q['T_OII_ne1_Kal09'])
                           for t, q in r.items() if isinstance(q, dict)})
    res['desi_7325_refit'] = out7
    # --------------------------------------------------- implications of n_e = 100
    T = 12500.0
    a = ats['Kal09+FFT04']
    imp = {}
    for ne in (1.0, 10.0, 100.0):
        eps = a.emis(T, ne, 2, 1) + a.emis(T, ne, 3, 1)        # erg cm^3 s^-1 per O+ per n_e
        hnu = 6.62607015e-27 * C * 1e5 / 3728.5e-8
        nO = 4.9e-4 * ne / 1.15                                 # all O as O+, n_H = n_e/1.15
        for I_R in (1.5, 2.5):
            L_cm = I_R * 1e6 / (ne * nO * eps / hnu)
            imp[f'ne{ne:g}_I{I_R}R'] = dict(L_pc=float(L_cm / 3.086e18), L_AU=float(L_cm / 1.496e13),
                                            fill_1p37pc=float(L_cm / 3.086e18 / 1.37),
                                            fill_25pc=float(L_cm / 3.086e18 / 25.3))
        imp[f'ne{ne:g}_P_over_k'] = float(2.2 * ne * T)               # n_tot ~ 2.2 n_e (H+, He, e) for ionised gas
        imp[f'ne{ne:g}_t_rec_yr'] = float(1 / (2.6e-13 * (T / 1e4) ** -0.8 * ne) / 3.156e7)
        # cooling time at ~1.25e4 K: Lambda ~ 2e-24 erg cm^3/s (CIE-ish, solar; order of magnitude)
        imp[f'ne{ne:g}_t_cool_yr'] = float(1.5 * 2.2 * ne * KB * T / (ne * ne / 1.15 * 2e-24) / 3.156e7)
        # Halpha from the same layer if H is ionised (case B, 4 pi j(Ha) = 2.86 x 1.24e-25 T4^-0.94 n_e n_p)
        j_ha = 2.86 * 1.24e-25 * (T / 1e4) ** -0.94 * 0.98
        L_cm = 2.0 * 1e6 / (ne * nO * eps / hnu)
        imp[f'ne{ne:g}_Ha_R_for_2R_OII'] = float(ne * ne / 1.15 * j_ha / (6.62607015e-27 * C * 1e5 / 6564.6e-8)
                                                 * L_cm / 1e6)
    res['implications'] = imp
    print(json.dumps(imp, indent=0)[:1500])
    json.dump(res, open(D + 'w3_oii_atomic.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
