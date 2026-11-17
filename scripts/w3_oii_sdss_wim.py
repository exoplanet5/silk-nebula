"""w3_oii_sdss_wim.py -- (wave 3, tag oii) low-density control: the Galactic WIM [O II] and [S II] doublets in the
ABSOLUTE total sky of BOSS/SDSS-V 2" sky fibres, isolated by regression on WHAM Halpha.

The WIM (n_e ~ 0.1-0.3 cm^-3) must give [O II] 3729/3726 and [S II] 6716/6731 at their low-density limits
(1.45-1.50 and 1.43-1.45). Zodiacal light (solar Fraunhofer spectrum), O2 Herzberg airglow and OH do not correlate
with Galactic Halpha, so the Halpha-regression spectrum beta(lambda) = d(total sky)/d(I_Halpha,WHAM) isolates the
WIM emission lines with the same instrument, LSF and fitting code as the filament ON-OFF stacks.
Instrument: WHAM-SS DR1 Halpha (Haffner+2003/2016; integrated -80..+80 km/s), mean of beams within 0.55 deg of
each fibre. Sample: BOSS dark groups (wave-2 definition); fibres > 0.45 deg from every strand and outside the
Halpha-patch 30 % contour (filament and patch excluded).
  (a) between groups: group mean spectra regressed on <Halpha>_g with nuisance regressors continuum level (zodi)
      and [O I]5577 (airglow); bootstrap over groups.
  (b) within groups: fibre minus group clipped mean regressed on Halpha_i - <Halpha>_g with the fibre's continuum
      offset as nuisance; bootstrap over groups.
Doublets in beta are fitted as in w3_oii_sdss.py (point-sampled Gaussians, sigma = mean wdisp, linear baseline;
v fixed +12 km/s helio, or free; sigma free). Halpha in beta (as-calibrated R per R) gives the BOSS diffuse-light
calibration as a by-product.
Output: data/derived/w3_oii_sdss_wim.json, data/derived/w3_oii_sdss_wim.npz
"""
import csv
import json
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import load_windows, flam_to_R_per_A  # noqa: E402
from verify_sdss_groups import clipped_mean  # noqa: E402
from w3_oii_common import L, design, grid_fit  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula/'
RNG = np.random.default_rng(20261119)
NB = 500
WINS = {'o2': (3700.0, 3760.0), 's2': (6700.0, 6760.0), 'ha': (6540.0, 6600.0), 'cont': (5040.0, 5100.0)}


def fcol(rows, k):
    out = []
    for r in rows:
        try:
            out.append(float(r[k]))
        except (ValueError, TypeError):
            out.append(np.nan)
    return np.array(out)


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def wham_at(ra, dec, rad=0.55):
    d = fits.open(ROOT + 'data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912-int-grid.fits')[1].data
    wc = SkyCoord(d['GAL-LON'].astype(float) * u.deg, d['GAL-LAT'].astype(float) * u.deg, frame='galactic').icrs
    T = cKDTree(xyz(wc.ra.deg, wc.dec.deg))
    I = np.asarray(d['INTEN'], float)
    idx = T.query_ball_point(xyz(ra, dec), 2 * np.sin(np.radians(rad) / 2))
    return np.array([np.mean(I[j]) if len(j) else np.nan for j in idx])


def ols_slopes(Y, X, w=None):
    """Per-pixel OLS of Y [n, P] on design X [n, k] (with intercept column included by caller)."""
    if w is None:
        w = np.ones(X.shape[0])
    A = (X * w[:, None]).T @ X
    B = (X * w[:, None]).T @ Y
    return np.linalg.solve(A, B)          # [k, P]


def fit_pair(x, S, bs, lines, s0, v0):
    """Doublet fits of a regression spectrum (uniform weights), bootstrap intervals."""
    Ys = np.vstack([S[None], bs])
    W = np.eye(x.size)
    out = {}
    FR = np.arange(0.8, 1.5001, 0.01)
    VG = np.arange(-30.0, 50.01, 1.0)
    for tag, vg, sg in (('fixed', [v0], [s0]), ('s0.90', [v0], [0.9 * s0]), ('s1.10', [v0], [1.1 * s0]),
                        ('s1.20', [v0], [1.2 * s0]), ('vfree', VG, [s0]), ('sfree', [v0], list(FR * s0)),
                        ('vsfree', VG[::2], list(FR[::2] * s0))):
        c2, par, p = grid_fit(Ys, W, (((v, s), design(x, lines, v, [s, s], npoly=2)) for v in vg for s in sg))
        r = np.array([q[1] / q[0] for q in p])
        tot = np.array([q[0] + q[1] for q in p])
        qq = np.nanpercentile(r[1:], [2.5, 16, 50, 84, 97.5])
        out[tag] = dict(ratio=float(r[0]), p2p5=float(qq[0]), p16=float(qq[1]), p84=float(qq[3]), p97p5=float(qq[4]),
                        v=float(par[0][0]), sigma=float(par[0][1]), s_over_wdisp=float(par[0][1] / s0),
                        total=float(tot[0]), total_p16=float(np.percentile(tot[1:], 16)),
                        total_p84=float(np.percentile(tot[1:], 84)),
                        v_p16=float(np.percentile([q[0] for q in par[1:]], 16)),
                        v_p84=float(np.percentile([q[0] for q in par[1:]], 84)))
    return out


def main():
    D = load_windows()
    lam = D['lam']
    F = list(csv.DictReader(open(ROOT + 'data/derived/verify_sdss_fibres.csv')))
    G = list(csv.DictReader(open(ROOT + 'data/derived/verify_sdss_groups.csv')))
    gid = np.array([int(r['gid']) for r in F])
    ra, dec, dmin, diam = fcol(F, 'ra'), fcol(F, 'dec'), fcol(F, 'dmin'), fcol(F, 'diam')
    in30 = np.array([int(r['in_patch30']) for r in F]).astype(bool)
    ins = np.array([r['instrument'] for r in G])
    sb, nf = fcol(G, 'sb5070_AB_cal'), fcol(G, 'n_fib')
    ok = (nf >= 6) & np.isfinite(sb)
    med = np.nanmedian(sb[ok & (ins == 'BOSS')])
    boss = ok & (ins == 'BOSS') & (sb > med - 0.75)
    I5577 = fcol(G, 'I_OI5577_cal')
    Ha = wham_at(ra, dec)
    use = (dmin > 0.45) & ~in30 & np.isfinite(Ha)
    TOT, IV, WD = D['tot'].astype(float), D['ivar'].astype(float), D['wdisp'].astype(float)
    res = {'instrument': 'WHAM-SS DR1 Halpha mean within 0.55 deg', 'n_boss_dark_groups': int(boss.sum())}
    npz = {}
    # ------------------------------------------------------------ per-group means and per-fibre residuals
    gm = {k: [] for k in WINS}
    gsig = {k: [] for k in WINS}
    gH, gC, g5, glist = [], [], [], []
    fib_res = {k: [] for k in WINS}
    fib_z, fib_c, fib_g = [], [], []
    for g in np.where(boss)[0]:
        s = np.where((gid == g) & use)[0]
        if len(s) < 6:
            continue
        dm = diam[s[0]]
        Xw = {}
        for k, (lo, hi) in WINS.items():
            m = (lam > lo) & (lam < hi)
            X = flam_to_R_per_A(TOT[s][:, m], lam[m], dm)
            V = np.where(IV[s][:, m] > 0, flam_to_R_per_A(1.0, lam[m], dm) ** 2 / np.where(IV[s][:, m] > 0,
                                                                                         IV[s][:, m], 1), np.inf)
            X = np.where(np.isfinite(V), X, np.nan)
            M = clipped_mean(X, V)[0]
            Xw[k] = (X, M)
            gm[k].append(M)
            gsig[k].append(np.nanmedian(WD[s][:, m], 0) * lam[m] * np.log(10) * 1e-4)
        cont = np.nanmean(Xw['cont'][0], 1)
        cm = np.nanmean(Xw['cont'][1])
        gH.append(np.mean(Ha[s]))
        gC.append(cm)
        g5.append(I5577[g])
        glist.append(g)
        for k in WINS:
            X, M = Xw[k]
            fib_res[k].append(X - M[None, :])
        fib_z.append(Ha[s] - np.mean(Ha[s]))
        fib_c.append(cont - cm)
        fib_g.append(np.full(len(s), len(glist) - 1))
    gH, gC, g5 = np.array(gH), np.array(gC), np.array(g5)
    ng = len(glist)
    res['n_groups_used'] = ng
    res['group_Ha_percentiles'] = np.percentile(gH, [0, 16, 50, 84, 100]).tolist()
    fz, fc, fg = np.concatenate(fib_z), np.concatenate(fib_c), np.concatenate(fib_g)
    res['n_fibres_within'] = int(fz.size)
    res['within_Ha_std'] = float(np.std(fz))
    boots_g = [RNG.integers(0, ng, ng) for _ in range(NB)]
    for k in WINS:
        lamk = lam[(lam > WINS[k][0]) & (lam < WINS[k][1])]
        npz[f'lam_{k}'] = lamk
        Y = np.array(gm[k])
        sig = np.nanmean(np.array(gsig[k]), 0)
        npz[f'sig_{k}'] = sig
        good = np.all(np.isfinite(Y), 1) & np.isfinite(g5) & np.isfinite(gC)
        # (a) between groups
        Xd = np.c_[np.ones(ng), gH, gC, np.nan_to_num(g5)]
        idx = np.where(good)[0]
        beta = ols_slopes(Y[idx], Xd[idx])[1]
        bb = []
        for bi in boots_g:
            j = idx[np.isin(idx, bi)]
            jj = bi[np.isin(bi, idx)]
            bb.append(ols_slopes(Y[jj], Xd[jj])[1])
        npz[f'between_{k}_S'] = beta
        npz[f'between_{k}_bs'] = np.array(bb, np.float32)
        # (b) within groups (pooled fibres, group-demeaned), nuisance: continuum offset
        R = np.concatenate(fib_res[k])
        okf = np.all(np.isfinite(R), 1) & np.isfinite(fc)
        Xf = np.c_[fz, fc]
        bw = ols_slopes(R[okf], Xf[okf])[0]
        bwb = []
        for bi in boots_g[:300]:
            cnt = np.bincount(bi, minlength=ng)[fg]
            bwb.append(ols_slopes(R[okf], Xf[okf], w=cnt[okf].astype(float))[0])
        npz[f'within_{k}_S'] = bw
        npz[f'within_{k}_bs'] = np.array(bwb, np.float32)
        print(k, 'done', flush=True)
    np.savez_compressed(ROOT + 'data/derived/w3_oii_sdss_wim.npz', **npz)
    # ------------------------------------------------------------ fits
    for mode in ('between', 'within'):
        r = {}
        for k, lines, lo, hi, l0, v0 in (('o2', [L['OII3726'], L['OII3729']], 3708.0, 3750.0, 3728.0, 12.0),
                                          ('s2', [L['SII6731'], L['SII6716']], 6706.0, 6750.0, 6724.0, 12.0)):
            x = npz[f'lam_{k}']
            m = (x > lo) & (x < hi)
            s0 = float(np.interp(l0, x, npz[f'sig_{k}']))
            r[k] = fit_pair(x[m], npz[f'{mode}_{k}_S'][m], npz[f'{mode}_{k}_bs'][:, m].astype(float),
                            lines, s0, v0)
            r[k]['sigma_wdisp'] = s0
        # Halpha slope (as calibrated R per WHAM R): single Gaussian + [N II]
        x = npz['lam_ha']
        m = (x > 6545) & (x < 6595)
        s0 = float(np.interp(6563.0, x, npz['sig_ha']))
        lines = [L['Ha'], L['NII6583']]
        Ys = np.vstack([npz[f'{mode}_ha_S'][m][None], npz[f'{mode}_ha_bs'][:, m].astype(float)])
        c2, par, p = grid_fit(Ys, np.eye(m.sum()), [((12.0, s0), design(x[m], lines, 12.0, [s0, s0], npoly=2))])
        ha = np.array([q[0] for q in p])
        r['ha_slope_cal'] = dict(value=float(ha[0]), p16=float(np.percentile(ha[1:], 16)),
                                 p84=float(np.percentile(ha[1:], 84)),
                                 note='as-calibrated R of Halpha per WHAM R; expected 1/f_eff = 1/0.60 = 1.67 '
                                      'if WHAM and SDSS agree (geocoronal Halpha uncorrelated with WHAM)')
        r['o2_total_per_Ha_photon_cal'] = r['o2']['fixed']['total']
        res[mode] = r
        o, s = r['o2'], r['s2']
        print(f'{mode:8s} [O II] fixed {o["fixed"]["ratio"]:.3f} [{o["fixed"]["p16"]:.3f},{o["fixed"]["p84"]:.3f}] '
              f'vfree {o["vfree"]["ratio"]:.3f} (v={o["vfree"]["v"]:+.0f} [{o["vfree"]["v_p16"]:+.0f},{o["vfree"]["v_p84"]:+.0f}]) '
              f'sfree {o["sfree"]["ratio"]:.3f} (s/wd={o["sfree"]["s_over_wdisp"]:.2f}) vsfree {o["vsfree"]["ratio"]:.3f} '
              f'[{o["vsfree"]["p16"]:.3f},{o["vsfree"]["p84"]:.3f}] | total {o["fixed"]["total"]:.3f} R/R(Ha) cal | '
              f'[S II] fixed {s["fixed"]["ratio"]:.3f} [{s["fixed"]["p16"]:.3f},{s["fixed"]["p84"]:.3f}] vsfree '
              f'{s["vsfree"]["ratio"]:.3f} | Ha slope {r["ha_slope_cal"]["value"]:.2f} '
              f'[{r["ha_slope_cal"]["p16"]:.2f},{r["ha_slope_cal"]["p84"]:.2f}]', flush=True)
    json.dump(res, open(ROOT + 'data/derived/w3_oii_sdss_wim.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
