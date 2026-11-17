"""verify_shockmodel_fit: grid of truncated (and complete) radiative shocks vs the Silk's wave-2 spectrum.

Model: verify_shockmodel_flow (1-D steady flow, NEI, no radiative transfer; validated in verify_shockmodel_validate.py
and verify_shockmodel_r88.py).  The flow is scale-free in n0 at fixed b = B0/sqrt(n0) (uG cm^1.5) apart from the
density-sensitive [O II] doublet, so tracks are computed at n0 = 0.1 and rescaled later (verify_shockmodel_analysis).

Grid: v_s = 70..200 km/s; preshock H0 fraction X = 0 ('precursor': H+, He+/He++, O+/O++), 0.05 (WIM-like),
0.3, 1.0 (neutral); b = 1, 3, 5, 7, 10, 14, 20, 30; abundances Asplund+09 (+ refractory-depleted variant for a subset).
Each run integrates until T = 1e4 K; every swept column N along the run is a candidate truncation.

Observed constraints (energy ratios, south-weighted; wave-2 values from notes/verify_desi.md and verify_sdss.md):
  [O III]/Hb 40 (DESI S 46 >28; SDSS S 23), [O II]3727/[O III] 1.35 (DESI S 1.11, all 1.37; SDSS S 1.84),
  [O III]/Ha 7 (DESI S 13.7; SDSS S 2.7, plane fit >6), [N II]6583/[O III] 0.12, [S II]/[O III] 0.05,
  [Ne III]/[O III] 0.10, 4363/5007 0.087, [O II]7325/[O III] 0.034, [N II]5755/6583 0.13,
  He II 4686/[O III] < 0.026, He I 5876/[O III] < 0.005.
Errors: observational sigma (dex) (+) model systematic 0.1 dex (+) 0.3 dex on Hb/Ha terms (R88 comparison:
my no-RT Hb at intermediate columns is 0.3-0.4 dex below Raymond's) (+) 0.1 dex abundance term for N, S, Ne.
Outputs: data/derived/verify_shockmodel_fit.json, verify_shockmodel_tracks.npz
"""
import sys, json, itertools, time, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_atomic as A
import verify_shockmodel_flow as F

ROOT = '/Users/mickey/astronomy/silk-nebula'
N0_REF = 0.1
HB_R88_FACTOR = 2.3   # multiply model Hb/Ha at intermediate columns to emulate R88-type (radiation-field) Hb

OBS = {  # key: (value, sigma_obs_dex, sigma_model_dex) or ('ul', value)
    'O3_Hb': (40.0, 0.20, 0.32), 'O2_O3': (1.35, 0.12, 0.10), 'O3_Ha': (7.0, 0.30, 0.32),
    'N2_O3': (0.12, 0.25, 0.14), 'S2_O3': (0.05, 0.30, 0.20), 'Ne3_O3': (0.10, 0.12, 0.14),
    'O4363': (0.087, 0.10, 0.10), 'O2r_O3': (0.034, 0.10, 0.10), 'N2T': (0.13, 0.25, 0.10),
    'HeII_O3': ('ul', 0.026), 'HeI_O3': ('ul', 0.005),
}
GROUPS = {'O3': ['OIII5007'], 'O4959': ['OIII4959'], 'Hb': ['Hb'], 'Ha': ['Ha'], 'O2': ['OII3726', 'OII3729'],
          'N2': ['NII6583'], 'N2a': ['NII5755'], 'S2': ['SII6716', 'SII6731'], 'Ne3': ['NeIII3869'],
          'O4363': ['OIII4363'], 'HeII': ['HeII4686'], 'HeI': ['HeI5876'], 'CIV': ['CIV1548', 'CIV1551'],
          'O3s': ['OIII1666', 'OIII1661'], 'HeII1640': ['HeII1640'], 'O2r': ['OII7320', 'OII7330', 'OII7319b', 'OII7330b'],
          'O2_3729': ['OII3729'], 'O2_3726': ['OII3726'], 'Hb_coll': ['Hb_coll'], 'Ha_coll': ['Ha_coll']}
LAM5007 = 5006.84


def add_hei(sh, j):
    """He I 5876 recombination (approximate: eps5876 = 1.40 (T/1e4)^-0.15 x eps_Hb(T) per n_e n_He+)."""
    nHep = sh.nH * sh.net.ab['He'] * sh.frac('He', 1)
    T = np.maximum(sh.T, 1e3)
    eps = 1.40 * (T / 1e4)**-0.15 * A.rec_emiss(4861.0, T) * 4861.0 / 5876.0
    j['HeI5876'] = sh.ne * nHep * eps
    return j


def ratios(G, hbfac=1.0):
    with np.errstate(divide='ignore', invalid='ignore'):
        Hb, Ha = G['Hb'] * hbfac, G['Ha'] * hbfac
        r = {'O3_Hb': G['O3'] / Hb, 'O2_O3': G['O2'] / G['O3'], 'O3_Ha': G['O3'] / Ha, 'N2_O3': G['N2'] / G['O3'],
             'S2_O3': G['S2'] / G['O3'], 'Ne3_O3': G['Ne3'] / G['O3'], 'O4363': G['O4363'] / G['O3'],
             'O2r_O3': G['O2r'] / G['O3'], 'N2T': G['N2a'] / G['N2'], 'HeII_O3': G['HeII'] / G['O3'],
             'HeI_O3': G['HeI'] / G['O3'], 'CIV_O3': G['CIV'] / G['O3'], 'O3s_O3': G['O3s'] / G['O3'],
             'HeII1640_O3': G['HeII1640'] / G['O3'], 'R23': (G['O2'] + G['O3'] + G['O4959']) / Hb,
             'Ha_Hb': G['Ha'] / G['Hb'], 'N2_Ha': G['N2'] / Ha, 'S2_Ha': G['S2'] / Ha,
             'O2_doublet': G['O2_3729'] / G['O2_3726'], 'Hbcoll_frac': G['Hb_coll'] / G['Hb']}
    return r


def chi2_terms(r, obs=OBS):
    out = {}
    for k, v in obs.items():
        x = np.maximum(np.nan_to_num(r[k], nan=1e-30, posinf=1e30), 1e-30)
        if v[0] == 'ul':
            out[k] = np.where(x > v[1], (np.log10(x / v[1]) / 0.1)**2, 0.0)
        else:
            s = np.hypot(v[1], v[2])
            out[k] = (np.log10(x / v[0]) / s)**2
    return out


def preshock(xn):
    if xn <= 0:
        return F.PRESHOCK['precursor']
    sw, Tw = F.PRESHOCK['wim']
    sn, Tn = F.PRESHOCK['neutral']
    w = np.clip((xn - 0.05) / 0.95, 0.0, 1.0)
    spec = {}
    for e in set(sw) | set(sn):
        keys = set(sw.get(e, {})) | set(sn.get(e, {}))
        spec[e] = {q: (1 - w) * sw.get(e, {}).get(q, 0.0) + w * sn.get(e, {}).get(q, 0.0) for q in keys}
    return (spec, (1 - w) * Tw + w * Tn)


def run(net, v, xn, b, complete=False):
    B0 = b * np.sqrt(N0_REF)
    sh = F.Shock(net, v, N0_REF, B0_uG=B0, pre=preshock(xn), complete=complete,
                 Tfloor=8e3 if complete else 1e4, Nmax=1e21).run()
    j = add_hei(sh, F.emissivities(sh, thin=1))
    C = F.cumulative(sh, j)
    # per-step quantities for later analysis
    nC4 = sh.nH * net.ab['C'] * sh.frac('C', 3)
    colC4 = np.concatenate([[0.0], np.cumsum(0.5 * (nC4[:-1] + nC4[1:]) * np.diff(sh.s))])
    nH0 = sh.nH * sh.frac('H', 0)
    colH0 = np.concatenate([[0.0], np.cumsum(0.5 * (nH0[:-1] + nH0[1:]) * np.diff(sh.s))])
    j5007 = j['OIII5007']
    return sh, j, C, colC4, colH0, j5007


def resample(sh, j, C, colC4, colH0, nN=220):
    N = sh.N
    Ng = np.logspace(np.log10(max(N[1], 1e14)), np.log10(N[-1]), nN)
    tr = {'N': Ng}
    for k in ('T', 's', 'ne', 't', 'x', 'nH'):
        tr[k] = np.interp(Ng, N, getattr(sh, k))
    for g, ks in GROUPS.items():
        tr['F_' + g] = np.interp(Ng, N, sum(C[k] for k in ks))
    tr['NCIV'] = np.interp(Ng, N, colC4)
    tr['NH0'] = np.interp(Ng, N, colH0)
    # emission-weighted means up to N (5007 and C IV): temperature, flow speed factor (1 - 1/x), density
    for lab, key in (('o3', 'OIII5007'), ('c4', 'CIV1548'), ('o2', 'OII3729')):
        w = j[key] * np.gradient(sh.s)
        cw = np.cumsum(w) + 1e-300
        for q, arr in (('T', sh.T), ('u', 1 - 1 / sh.x), ('ne', sh.ne), ('s', sh.s)):
            tr[f'w{lab}_{q}'] = np.interp(Ng, N, np.cumsum(w * arr) / cw)
    return tr


def best_on_track(tr, hbfac=1.0, obs=OBS):
    G = {g: tr['F_' + g] for g in GROUPS}
    r = ratios(G, hbfac)
    terms = chi2_terms(r, obs)
    c2 = sum(terms.values())
    c2 = np.where(np.isfinite(c2) & (tr['N'] > 1e15), c2, 1e9)
    i = int(np.argmin(c2))
    return i, float(c2[i]), {k: float(v[i]) for k, v in r.items()}, {k: float(v[i]) for k, v in terms.items()}, c2


if __name__ == '__main__':
    t0 = time.time()
    nets = {'asplund09': F.Network('asplund09'), 'ism_dep': F.Network('ism_dep')}
    VS = [70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200]
    XN = [0.0, 0.05, 0.3, 1.0]
    BB = [1.0, 3.0, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0]
    runs = [dict(v=v, xn=x, b=b, ab='asplund09', complete=False) for v, x, b in itertools.product(VS, XN, BB)]
    runs += [dict(v=v, xn=x, b=3.0, ab='ism_dep', complete=False) for v, x in itertools.product((90, 110, 130, 160), (0.0, 0.05))]
    runs += [dict(v=v, xn=x, b=3.0, ab='asplund09', complete=True) for v, x in itertools.product((80, 100, 120, 140, 160, 200), (0.0, 1.0))]
    tracks, summary = {}, []
    for k, rp in enumerate(runs):
        vA = rp['b'] * 1e-6 / np.sqrt(4 * np.pi * nets[rp['ab']].mass_per_H) / 1e5
        if rp['v'] < 1.5 * vA:
            continue
        try:
            sh, j, C, colC4, colH0, _ = run(nets[rp['ab']], rp['v'], rp['xn'], rp['b'], rp['complete'])
        except Exception as e:  # noqa
            print('FAILED', rp, e)
            continue
        tag = f"v{rp['v']}_x{rp['xn']}_b{rp['b']}_{rp['ab']}_{'cmp' if rp['complete'] else 'trn'}"
        tr = resample(sh, j, C, colC4, colH0)
        for key, arr in tr.items():
            tracks[f'{tag}__{key}'] = arr
        i, c2, r, terms, _ = best_on_track(tr)
        i2, c2b, r2, terms2, _ = best_on_track(tr, HB_R88_FACTOR)
        eta = (tr['F_O3'] / (A.H_PL * A.C_L / (LAM5007 * 1e-8))) / (N0_REF * rp['v'] * 1e5)  # 5007 photons per H swept
        rec = dict(tag=tag, **rp, vA=float(vA), Ts=float(sh.Ts), x1=float(sh.x1), N_end=float(sh.N[-1]),
                   t_end_n0=float(sh.t[-1] * N0_REF), s_end_n0=float(sh.s[-1] * N0_REF),
                   eta_end=float(eta[-1]),
                   best={'i': i, 'chi2': c2, 'N': float(tr['N'][i]), 'T': float(tr['T'][i]), 'eta5007': float(eta[i]),
                         's_n0': float(tr['s'][i] * N0_REF), 't_n0': float(tr['t'][i] * N0_REF),
                         'NCIV': float(tr['NCIV'][i]), 'NH0': float(tr['NH0'][i]), 'x': float(tr['x'][i]),
                         'wo3_T': float(tr['wo3_T'][i]), 'wo3_u': float(tr['wo3_u'][i]), 'wc4_T': float(tr['wc4_T'][i]),
                         'wo2_ne_over_n0': float(tr['wo2_ne'][i] / N0_REF), 'ratios': r, 'terms': terms},
                   best_hbR88={'i': i2, 'chi2': c2b, 'N': float(tr['N'][i2]), 'T': float(tr['T'][i2]),
                               'eta5007': float(eta[i2]), 'ratios': r2, 'terms': terms2})
        summary.append(rec)
        print(f"{k+1}/{len(runs)} {tag}: Ts={sh.Ts:.3g} x_end={sh.x[-1]:.1f} chi2={c2:.1f} logN*={np.log10(tr['N'][i]):.2f} "
              f"T*={tr['T'][i]:.3g} O3/Hb={r['O3_Hb']:.0f} O2/O3={r['O2_O3']:.2f} O3/Ha={r['O3_Ha']:.1f} "
              f"4363={r['O4363']:.3f} N2T={r['N2T']:.2f} O2r={r['O2r_O3']:.3f} HeII={r['HeII_O3']:.3f} "
              f"CIV/O3={r['CIV_O3']:.1f} eta={eta[i]:.3f} | R88-Hb chi2={c2b:.1f} logN={np.log10(tr['N'][i2]):.2f} "
              f"({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(ROOT + '/data/derived/verify_shockmodel_tracks.npz', **tracks)
    json.dump({'obs': {k: list(v) for k, v in OBS.items()}, 'n0_ref': N0_REF, 'hb_r88_factor': HB_R88_FACTOR,
               'runs': summary}, open(ROOT + '/data/derived/verify_shockmodel_fit.json', 'w'), indent=1)
