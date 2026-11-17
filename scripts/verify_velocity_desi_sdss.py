"""verify_velocity: independent re-measurement of the [O III] 5007 heliocentric velocity in SDSS sky fibres and
DESI DR1 residual spectra (own code; wave-1 products used only as input arrays).

SDSS (data/derived/spectra_sdss_lamost_sdss_skyfibres.npz: spec-lite windows of 9,092 SKY fibres; vacuum,
heliocentric log-lambda grid, HELIO_RV = -v_corr):
  * total sky = flux + sky; per plate-MJD-spectrograph group, ON (<0.15 deg from E/W strand of silk_trace v1.0)
    minus OFF (0.45-2.0 deg) mean spectrum;
  * own zero point per group from [O I] 5577 in the group-mean total sky, expected at
    lam0*(1 - HELIO_RV/c) (topocentric rest);
  * stack in velocity space after the zero-point shift; Gaussian (+ tied 4959) fit with the median LSF;
    bootstrap over groups.
DESI (data/derived/spectra_desi_windows/*.npz: SPARCL coadds, vacuum, barycentric; model = Redrock):
  * residual = flux - model, local linear continuum; stack of fibres with predicted positive [O III] residual
    (wave-1 forward-model column model_resid_R > 0.5 R) and of fibres within 0.3 deg of P1/P2/P3;
  * Gaussian fit with LSF from wave_sigma; bootstrap over fibres.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_desi_sdss.py
Output: data/derived/verify_velocity_desi_sdss.json
"""
import glob
import json
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/verify_velocity_desi_sdss.json'
C = 299792.458


def n_air(lv):
    s2 = (1e4 / np.asarray(lv, float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


def air2vac(la):
    v = np.asarray(la, float)
    for _ in range(5):
        v = la * n_air(v)
    return float(v)


O3 = air2vac(5006.843)
O3B = air2vac(4958.911)
OI5577 = air2vac(5577.339)
P = {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72)}


def strand_dist(ra, dec):
    t = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
    d = np.full(np.shape(ra), 99.0)
    which = np.full(np.shape(ra), '', dtype='<U2')
    for name in ('E', 'W'):
        pts = t['strands'][name]['points']
        pr = np.array([p['ra'] for p in pts])
        pd = np.array([p['dec'] for p in pts])
        # densify polyline
        s = np.linspace(0, 1, 20)
        R = np.concatenate([pr[i] + (pr[i + 1] - pr[i]) * s for i in range(len(pr) - 1)])
        D = np.concatenate([pd[i] + (pd[i + 1] - pd[i]) * s for i in range(len(pd) - 1)])
        cosd = np.cos(np.radians(np.asarray(dec)))[:, None]
        dd = np.sqrt(((np.asarray(ra)[:, None] - R[None]) * cosd) ** 2 + (np.asarray(dec)[:, None] - D[None]) ** 2).min(1)
        # only within the strand's Dec range (+0.1)
        inr = (np.asarray(dec) > D.min() - 0.1) & (np.asarray(dec) < D.max() + 0.1)
        dd = np.where(inr, dd, 99.0)
        better = dd < d
        d = np.where(better, dd, d)
        which = np.where(better, name, which)
    return d, which


def gfit_v(u, y, sig, w=None, two=True, mu0=10.0):
    """fit 5007 (+ tied 4959 placed by its own velocity grid handled outside) Gaussian+const in velocity space"""
    def f(x, A, mu, c0, c1):
        return A * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c0 + c1 * x / 1000
    m = np.abs(u) < 400
    p, cov = curve_fit(f, u[m], y[m], p0=[np.nanmax(y[m]), mu0, 0, 0], sigma=None if w is None else 1 / np.sqrt(w[m]),
                       maxfev=20000)
    return p, np.sqrt(np.diag(cov))


# ----------------------------------------------------------------------------- SDSS
def sdss():
    z = np.load(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz')
    lam = 10 ** z['loglam']
    tot = z['flux'] + z['sky']
    ivar = z['ivar']
    ra, dec = z['ra'], z['dec']
    hrv = z['helio_rv']
    wd = z['wdisp']
    d, which = strand_dist(ra, dec)
    key = np.array([f'{p}-{m}-{s}' for p, m, s in zip(z['plate'], z['mjd'], z['spectrograph'])])
    on = d < 0.15
    off = (d > 0.45) & (d < 2.0)
    m5 = (lam > 4985) & (lam < 5030)
    mo = (lam > 5566) & (lam < 5592)
    ug = np.arange(-600, 601, 5.0)
    rows = []
    for k in np.unique(key):
        g = key == k
        if (g & on).sum() < 2 or (g & off).sum() < 4:
            continue
        # zero point from [O I] 5577 in the group-mean total sky
        sky_mean = np.nanmean(tot[g][:, mo], 0)
        x = lam[mo]
        try:
            p, cov = curve_fit(lambda xx, A, mu, s, c0, c1: A * np.exp(-0.5 * ((xx - mu) / s) ** 2) + c0 + c1 * (xx - 5578),
                               x, sky_mean, p0=[sky_mean.max(), 5578.9, 1.3, np.median(sky_mean), 0], maxfev=20000)
        except Exception:
            continue
        hr = np.median(hrv[g])
        lexp = OI5577 * (1 - hr / C)
        dz = C * (p[1] / lexp - 1)
        if not np.isfinite(dz) or abs(dz) > 30:
            continue
        diff = np.nanmean(tot[g & on][:, m5], 0) - np.nanmean(tot[g & off][:, m5], 0)
        xl = lam[m5]
        base = (np.abs(xl - O3) > 4) & (np.abs(xl - O3B) > 4)
        cc = np.polyfit(xl[base], diff[base], 1)
        diff = diff - np.polyval(cc, xl)
        u = C * (xl / O3 - 1) - dz
        ub = C * (xl / O3B - 1) - dz
        sig_v = np.median(wd[g][:, m5]) * 69.03
        nonn, noff = int((g & on).sum()), int((g & off).sum())
        rows.append(dict(key=k, dz=dz, helio_rv=float(hr), n_on=nonn, n_off=noff, dec=float(np.median(dec[g & on])),
                         sig_v=float(sig_v), y=np.interp(ug, u, diff), yb=np.interp(ug, ub, diff),
                         w=nonn * noff / (nonn + noff), inst='sdss' if z['instrument'][g][0] in (b'sdss', 'sdss') else str(z['instrument'][g][0])))
    sig_all = np.hypot(np.median([r["sig_v"] for r in rows]), 20 / 2.3548)

    def stack(sel):
        W = np.array([r['w'] for r in sel])
        Y = np.array([r['y'] for r in sel])
        YB = np.array([r['yb'] for r in sel])
        y = (W[:, None] * Y).sum(0) / W.sum()
        yb = (W[:, None] * YB).sum(0) / W.sum()
        comb = y  # 5007 only for the centroid (4959 reported separately)
        p, e = gfit_v(ug, comb, sig_all)
        pb, eb = gfit_v(ug, yb, sig_all)
        return p, e, pb, eb

    res = {}
    rng = np.random.default_rng(11)
    for name, f in (('all', lambda r: True), ('south_dec<20', lambda r: r['dec'] < 20), ('north_dec>=20', lambda r: r['dec'] >= 20)):
        sel = [r for r in rows if f(r)]
        p, e, pb, eb = stack(sel)
        boots = []
        for _ in range(300):
            bs = [sel[i] for i in rng.integers(0, len(sel), len(sel))]
            try:
                boots.append(stack(bs)[0][1])
            except Exception:
                pass
        res[name] = dict(ngroups=len(sel), V_helio_5007=float(p[1]), err_fit=float(e[1]), err_boot=float(np.std(boots)),
                         amp=float(p[0]), V_helio_4959=float(pb[1]), err_4959=float(eb[1]),
                         median_dz=float(np.median([r['dz'] for r in sel])))
    res['zero_point_dz_all'] = dict(median=float(np.median([r['dz'] for r in rows])),
                                    p16=float(np.percentile([r['dz'] for r in rows], 16)),
                                    p84=float(np.percentile([r['dz'] for r in rows], 84)))
    res['sigma_LSF_used_kms'] = float(sig_all)
    return res, rows, ug


# ----------------------------------------------------------------------------- DESI
def desi_wave():
    lines_air = [3726.03, 3728.82, 3868.76, 4101.74, 4340.47, 4363.21, 4685.71, 4861.33, 4958.91, 5006.84, 5197.90,
                 5200.26, 5754.64, 5875.62, 6300.30, 6363.78, 6548.05, 6562.80, 6583.45, 6678.15, 6716.44, 6730.82,
                 7065.19, 7135.80, 7319.99, 7330.73, 9068.6, 9530.6, 5577.34, 6300.30, 5199.0]
    WAVE = 3600.0 + 0.8 * np.arange(7781)
    c = sorted(set(round(air2vac(l), 2) for l in lines_air))
    segs = []
    for w in c:
        lo, hi = w - 40, w + 40
        if segs and lo <= segs[-1][1]:
            segs[-1][1] = max(segs[-1][1], hi)
        else:
            segs.append([lo, hi])
    idx = np.concatenate([np.where((WAVE >= max(lo, WAVE[0])) & (WAVE <= min(hi, WAVE[-1])))[0] for lo, hi in segs])
    return WAVE[idx]


def desi():
    wave = desi_wave()
    files = sorted(glob.glob(str(ROOT / 'data/derived/spectra_desi_windows/*.npz')))
    ids, F, IV, M, S = [], [], [], [], []
    m5 = (wave > 4975) & (wave < 5040)
    for f in files:
        z = np.load(f)
        if z['flux'].shape[1] != wave.size:
            raise RuntimeError(f'window size mismatch {z["flux"].shape} vs {wave.size}')
        ids.append(z['sparcl_id'])
        F.append(z['flux'][:, m5]); IV.append(z['ivar'][:, m5]); M.append(z['model'][:, m5])
        S.append(z['wave_sigma'][:, m5] if 'wave_sigma' in z.files else np.full((len(z['sparcl_id']), m5.sum()), np.nan))
    ids = np.concatenate(ids); F = np.vstack(F); IV = np.vstack(IV); M = np.vstack(M); S = np.vstack(S)
    x = wave[m5]
    t = Table.read(ROOT / 'data/derived/spectra_desi_fiber_table.csv', format='ascii.csv')
    pos = {s: i for i, s in enumerate(t['sparcl_id'])}
    j = np.array([pos.get(s, -1) for s in ids])
    keep = (j >= 0)
    ids, F, IV, M, S, j = ids[keep], F[keep], IV[keep], M[keep], S[keep], j[keep]
    # unique
    _, ui = np.unique(ids, return_index=True)
    ids, F, IV, M, S, j = ids[ui], F[ui], IV[ui], M[ui], S[ui], j[ui]
    T = t[j]
    R = F - M
    base = (np.abs(x - O3) > 5) & (np.abs(x - O3B) > 5) & (np.abs(x - air2vac(5015.68)) > 3)
    for i in range(R.shape[0]):
        ww = IV[i] > 0
        b = base & ww
        if b.sum() > 10:
            R[i] -= np.polyval(np.polyfit(x[b], R[i][b], 1), x)
    u = C * (x / O3 - 1)
    sig_lsf = np.nanmedian(S[:, np.argmin(np.abs(x - O3))]) / O3 * C
    sig = np.hypot(sig_lsf, 20.0 / 2.3548)

    def stack(sel):
        w = np.where(IV[sel] > 0, IV[sel], 0) * np.clip(T['ffrac'][sel], 0.3, 1)[:, None] ** 0
        y = (w * R[sel]).sum(0) / np.maximum(w.sum(0), 1e-30)
        wy = w.sum(0)
        p, e = gfit_v(u, y, sig, w=wy)
        return p, e, y

    res = {'sigma_LSF_kms': float(sig_lsf), 'n_spectra': int(len(ids))}
    okm = np.asarray(T['ok']) == 1
    sels = {
        'pred_pos_gt0.5R': okm & (np.asarray(T['model_resid_R']) > 0.5),
        'pred_pos_gt1.0R': okm & (np.asarray(T['model_resid_R']) > 1.0),
        'pred_pos_gt0.5R_south_dec<17.6': okm & (np.asarray(T['model_resid_R']) > 0.5) & (np.asarray(T['dec']) < 17.6),
        'pred_neg_lt-0.3R': okm & (np.asarray(T['model_resid_R']) < -0.3),
    }
    for pn, (pra, pdec) in P.items():
        dd = np.hypot((np.asarray(T['ra']) - pra) * np.cos(np.radians(pdec)), np.asarray(T['dec']) - pdec)
        sels[f'within0.3_{pn}'] = okm & (dd < 0.3)
        sels[f'within0.3_{pn}_pos'] = okm & (dd < 0.3) & (np.asarray(T['OIII_R']) > 0)
    rng = np.random.default_rng(7)
    stacks = {}
    for name, sel in sels.items():
        idx = np.where(sel)[0]
        if idx.size < 5:
            res[name] = dict(n=int(idx.size))
            continue
        try:
            p, e, y = stack(idx)
        except Exception as ex:
            res[name] = dict(n=int(idx.size), error=str(ex))
            continue
        boots = []
        for _ in range(300):
            bi = rng.choice(idx, idx.size)
            try:
                boots.append(stack(bi)[0][1])
            except Exception:
                pass
        res[name] = dict(n=int(idx.size), V_helio=float(p[1]), err_fit=float(e[1]), err_boot=float(np.std(boots)),
                         amp=float(p[0]))
        stacks[name] = y
    return res, u, stacks


def sdss_frame_check():
    """Observed [O I] 5577 velocity (vs vacuum rest) in each group-mean total sky against HELIO_RV.
    Topocentric airglow on a heliocentric grid with HELIO_RV = -v_corr should give slope -1."""
    z = np.load(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz')
    lam = 10 ** z['loglam']
    tot = z['flux'] + z['sky']
    key = np.array([f'{p}-{m}-{s}' for p, m, s in zip(z['plate'], z['mjd'], z['spectrograph'])])
    mo = (lam > 5566) & (lam < 5592)
    H, D = [], []
    for k in np.unique(key):
        g = key == k
        sm = np.nanmean(tot[g][:, mo], 0)
        try:
            p, _ = curve_fit(lambda xx, A, mu, s, c0, c1: A * np.exp(-0.5 * ((xx - mu) / s) ** 2) + c0 + c1 * (xx - 5578),
                             lam[mo], sm, p0=[sm.max(), 5578.9, 1.3, np.median(sm), 0], maxfev=20000)
        except Exception:
            continue
        H.append(np.median(z['helio_rv'][g]))
        D.append(C * (p[1] / OI5577 - 1))
    H, D = np.array(H), np.array(D)
    ok = np.abs(D) < 80
    a, b = np.polyfit(H[ok], D[ok], 1)
    return dict(n=int(ok.sum()), slope=float(a), intercept=float(b), rms_if_minus=float(np.std(D[ok] + H[ok])),
                rms_if_plus=float(np.std(D[ok] - H[ok])))


def main():
    s_res, rows, ug = sdss()
    d_res, u, stacks = desi()
    out = dict(sdss=s_res, desi=d_res, sdss_frame_check=sdss_frame_check())
    OUT.write_text(json.dumps(out, indent=1))
    np.savez(ROOT / 'data/derived/verify_velocity_desi_sdss_stacks.npz', desi_u=u,
             **{f'desi_{k}': v for k, v in stacks.items()}, sdss_u=ug,
             sdss_all=np.average([r['y'] for r in rows], axis=0, weights=[r['w'] for r in rows]),
             sdss_south=np.average([r['y'] for r in rows if r['dec'] < 20], axis=0,
                                   weights=[r['w'] for r in rows if r['dec'] < 20]))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
