"""verify_shockmodel_grid: grid of truncated / complete radiative-shock models vs the Silk's measured spectrum.

For each (v_s, preshock H0 fraction X, B, abundances) the 1-D model is integrated until T = 1e4 K (truncated
family: every intermediate swept column N is a 'truncated' model), plus complete models (isothermal 7000 K tail,
validated against Allen+08 at 100 km/s) for v <= 140 km/s.  Cumulative line ratios vs N are compared with the
observed ratios (wave-1 + task statement) by a chi^2 built from log ratios; upper limits as one-sided terms.
Outputs: data/derived/verify_shockmodel_grid.npz (tracks), verify_shockmodel_grid.json (best matches).
"""
import sys, json, itertools, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_flow as F

ROOT = '/Users/mickey/astronomy/silk-nebula'

# ---------------- observed constraints (energy ratios; filament-differential where available)
OBS = {
    # name: (value, sigma_dex) for log-normal terms; ('ul', value) for upper limits
    'O3_Hb': (20.0, 0.12),          # DESI regression 20 (15-26); SDSS south > 13.8 (2 sigma)
    'O2_O3': (1.30, 0.10),          # [O II]3726+29/[O III]5007: 1.49 all, 1.12 Dec<19 (DESI wave 1)
    'O3_Ha': (5.3, 0.10),           # SDSS south 4.7+-1.0; DESI positive-residual 6.1 (4.8-7.7)
    'N2_O3': (0.14, 0.15),          # DESI 0.10 (0.09-0.14); SDSS south 0.20+-0.04
    'S2_O3': (0.10, 0.25),          # DESI 0.068; SDSS all-strand 0.19 (+-0.14)
    'Ne3_O3': (0.10, 0.12),         # DESI 0.10 (0.07-0.13)
    'O4363_5007': (0.056, 0.20),    # DESI 0.056 (0.022-0.077), ~2 sigma, Hg-contaminable
    'HeII_O3': ('ul', 0.034),       # 97.5 %
}
LINE_GROUPS = {'O3': ['OIII5007'], 'Hb': ['Hb'], 'Ha': ['Ha'], 'O2': ['OII3726', 'OII3729'], 'N2': ['NII6583'],
               'S2': ['SII6716', 'SII6731'], 'Ne3': ['NeIII3869'], 'O4363': ['OIII4363'], 'HeII': ['HeII4686'],
               'CIV': ['CIV1548', 'CIV1551'], 'O3s': ['OIII1666', 'OIII1661'], 'HeII1640': ['HeII1640'],
               'O2r': ['OII7320', 'OII7330', 'OII7319b', 'OII7330b'], 'O4959': ['OIII4959'],
               'Hb_coll': ['Hb_coll'], 'O2_3729': ['OII3729'], 'O2_3726': ['OII3726']}


def group_cum(C, g):
    return sum(C[k] for k in LINE_GROUPS[g])


def ratios(C):
    g = {k: group_cum(C, k) for k in LINE_GROUPS}
    with np.errstate(divide='ignore', invalid='ignore'):
        r = {'O3_Hb': g['O3'] / g['Hb'], 'O2_O3': g['O2'] / g['O3'], 'O3_Ha': g['O3'] / g['Ha'],
             'N2_O3': g['N2'] / g['O3'], 'S2_O3': g['S2'] / g['O3'], 'Ne3_O3': g['Ne3'] / g['O3'],
             'O4363_5007': g['O4363'] / g['O3'], 'HeII_O3': g['HeII'] / g['O3'],
             'N2_Ha': g['N2'] / g['Ha'], 'S2_Ha': g['S2'] / g['Ha'], 'CIV_O3': g['CIV'] / g['O3'],
             'O3s_O3': g['O3s'] / g['O3'], 'HeII1640_O3': g['HeII1640'] / g['O3'], 'O2r_O3': g['O2r'] / g['O3'],
             'R23': (g['O2'] + g['O3'] + g['O4959']) / g['Hb'], 'Ha_Hb': g['Ha'] / g['Hb'],
             'Hbcoll_frac': g['Hb_coll'] / g['Hb'], 'O2_doublet': g['O2_3729'] / g['O2_3726']}
    return r, g


def chi2(r, obs=OBS):
    tot = 0.0
    for k, v in obs.items():
        x = r[k]
        if v[0] == 'ul':
            tot = tot + np.where(x > v[1], (np.log10(np.maximum(x, 1e-30) / v[1]) / 0.1)**2, 0.0)
        else:
            tot = tot + (np.log10(np.maximum(x, 1e-30) / v[0]) / v[1])**2
    return tot


def preshock_mix(net, xn):
    """Preshock state with H0 fraction xn, linear mix of 'wim' (xH0=0.05) and 'neutral' (xH0=1) states."""
    sw, Tw = F.PRESHOCK['wim']
    sn, Tn = F.PRESHOCK['neutral']
    if xn <= 0.0:
        return F.PRESHOCK['precursor']
    w = np.clip((xn - 0.05) / 0.95, 0.0, 1.0)
    spec = {}
    for e in set(sw) | set(sn):
        keys = set(sw.get(e, {})) | set(sn.get(e, {}))
        spec[e] = {q: (1 - w) * sw.get(e, {}).get(q, 0.0) + w * sn.get(e, {}).get(q, 0.0) for q in keys}
    return (spec, (1 - w) * Tw + w * Tn)


def run_model(net, v, xn, n0=0.1, B0=1.0, complete=False, Tfloor=1e4):
    pre = preshock_mix(net, xn)
    sh = F.Shock(net, v, n0, B0_uG=B0, pre=pre, complete=complete, Tfloor=Tfloor).run()
    j = F.emissivities(sh, thin=2)
    C = F.cumulative(sh, j)
    return sh, C


def resample(sh, C, nN=160):
    """Tracks on a log grid of swept column (from 1e15 to the end)."""
    N = sh.N
    Ng = np.logspace(np.log10(max(N[1], 1e14)), np.log10(N[-1]), nN)
    out = {'N': Ng, 'T': np.interp(Ng, N, sh.T), 's': np.interp(Ng, N, sh.s), 'ne': np.interp(Ng, N, sh.ne),
           't': np.interp(Ng, N, sh.t)}
    for k, a in C.items():
        out['F_' + k] = np.interp(Ng, N, a)
    # C IV column (face-on) behind the shock
    nC4 = sh.nH * sh.net.ab['C'] * sh.frac('C', 3)
    colC4 = np.concatenate([[0.0], np.cumsum(0.5 * (nC4[:-1] + nC4[1:]) * np.diff(sh.s))])
    out['NCIV'] = np.interp(Ng, N, colC4)
    return out


def best_match(tr, obs=OBS):
    C = {k[2:]: v for k, v in tr.items() if k.startswith('F_')}
    r, g = ratios(C)
    c2 = chi2(r, obs)
    c2 = np.where(np.isfinite(c2), c2, 1e9)
    i = int(np.argmin(c2))
    return i, float(c2[i]), {k: float(v[i]) for k, v in r.items()}


if __name__ == '__main__':
    import time
    t0 = time.time()
    VS = [80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200]
    XN = [0.0, 0.05, 0.2, 0.5, 1.0]
    runs = []
    nets = {'asplund09': F.Network('asplund09'), 'ism_dep': F.Network('ism_dep')}
    for v, xn in itertools.product(VS, XN):
        runs.append(dict(v=v, xn=xn, n0=0.1, B0=1.0, ab='asplund09', complete=False))
    for v in (100, 130, 160, 200):
        runs.append(dict(v=v, xn=0.05, n0=0.1, B0=3.0, ab='asplund09', complete=False))
        runs.append(dict(v=v, xn=0.05, n0=0.1, B0=1.0, ab='ism_dep', complete=False))
        runs.append(dict(v=v, xn=1.0, n0=0.1, B0=1.0, ab='ism_dep', complete=False))
    for v in (120, 160):
        runs.append(dict(v=v, xn=0.05, n0=1.0, B0=3.16, ab='asplund09', complete=False))
        runs.append(dict(v=v, xn=0.05, n0=0.03, B0=0.55, ab='asplund09', complete=False))
    for v, xn in itertools.product((80, 90, 100, 110, 120, 130, 140), (0.0, 0.05, 1.0)):
        runs.append(dict(v=v, xn=xn, n0=0.1, B0=1.0, ab='asplund09', complete=True))
    tracks, summary = {}, []
    for k, rp in enumerate(runs):
        net = nets[rp['ab']]
        sh, C = run_model(net, rp['v'], rp['xn'], rp['n0'], rp['B0'], complete=rp['complete'],
                          Tfloor=7e3 if rp['complete'] else 1e4)
        tag = f"v{rp['v']}_x{rp['xn']}_n{rp['n0']}_B{rp['B0']}_{rp['ab']}_{'cmp' if rp['complete'] else 'trn'}"
        tr = resample(sh, C)
        for key, arr in tr.items():
            tracks[f'{tag}__{key}'] = arr
        i, c2, r = best_match(tr)
        Cend = {kk[2:]: vv[-1:] for kk, vv in tr.items() if kk.startswith('F_')}
        rend, gend = ratios(Cend)
        phot5007 = F.photons(tr['F_OIII5007'], 'OIII5007') / (rp['n0'] * rp['v'] * 1e5)
        rec = dict(tag=tag, **rp, Ts=float(sh.Ts), x1=float(sh.x1), N_end=float(sh.N[-1]), t_end=float(sh.t[-1]),
                   s_end=float(sh.s[-1]), best_i=i, best_chi2=c2, best_N=float(tr['N'][i]), best_T=float(tr['T'][i]),
                   best_s_pc=float(tr['s'][i] / 3.0857e18), best_t_yr=float(tr['t'][i] / 3.156e7),
                   best_eta5007=float(phot5007[i]), best_NCIV=float(tr['NCIV'][i]),
                   best_ne=float(tr['ne'][i]), best_ratios=r,
                   end_ratios={kk: float(vv[0]) for kk, vv in rend.items()}, end_eta5007=float(phot5007[-1]))
        summary.append(rec)
        print(f"{k+1}/{len(runs)} {tag}: Ts={sh.Ts:.3g} chi2={c2:.1f} N*={tr['N'][i]:.2e} T*={tr['T'][i]:.3g} "
              f"O3/Hb={r['O3_Hb']:.1f} O2/O3={r['O2_O3']:.2f} O3/Ha={r['O3_Ha']:.1f} N2/O3={r['N2_O3']:.2f} "
              f"S2/O3={r['S2_O3']:.2f} 4363={r['O4363_5007']:.3f} HeII={r['HeII_O3']:.3f} eta={phot5007[i]:.2f} "
              f"| end O3/Hb={rend['O3_Hb'][0]:.1f} O2/O3={rend['O2_O3'][0]:.2f}  ({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(ROOT + '/data/derived/verify_shockmodel_grid.npz', **tracks)
    json.dump({'obs': {k: list(v) if isinstance(v, tuple) else v for k, v in OBS.items()}, 'runs': summary},
              open(ROOT + '/data/derived/verify_shockmodel_grid.json', 'w'), indent=1)
