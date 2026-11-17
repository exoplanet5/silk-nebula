"""judge_obs follow-up to w3_judge_obs_edge.py: robustness checks of the two signals it found.

(a) North [O III] asymmetry (sharp east edge, west tail): split the north into sub-segments that avoid the W2 strand
    (Dec 21.5-24.25, 0.6-1.3 deg west) and mask the Halpha patch 30 % contour; refit with symmetric-truth injections.
(b) GALEX FUV/NUV ridge east of the [O III] ridge: calibrate the centre offset against GALEX sky structure by
    injecting a CO-SPATIAL copy of the observed NSNS ridge shape (scaled to the GALEX amplitude) into off-Silk GALEX
    placements, and refitting exactly as for the real data; also vary the fit window and baseline.

Shares models/estimators with w3_judge_obs_edge.py (copied here so this script runs alone).
Outputs: data/derived/w3_judge_obs_edge2.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_obs_edge2.py
"""
import sys
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter1d
from matplotlib.path import Path as MPath
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import CurveCoords, radec_to_vec, load_spine  # noqa: E402
from w3_figs_common import DER, load_trace  # noqa: E402

RNG = np.random.default_rng(7)
ra_sp, dec_sp, _ = load_spine()
CC = CurveCoords(ra_sp, dec_sp)
TR = load_trace()


def s_at(dec):
    return float(np.interp(dec, CC.dec, CC.s))


BW = 0.02
EDGES = np.arange(-1.6, 1.6001, BW)
XC = 0.5 * (EDGES[1:] + EDGES[:-1])
FITW = 0.9
BLOCK = 0.25
DXF = 0.004
XF = np.arange(-2.5, 2.5001, DXF)


def kern(u):
    out = np.zeros_like(u)
    a = (u > 0) & (u <= 1)
    out[a] = np.sqrt(u[a])
    b = u > 1
    out[b] = np.sqrt(u[b]) - np.sqrt(u[b] - 1)
    return out


def sheet(x, A, xf, w, sp, b0, b1, side):
    h = side * (xf - XF)
    ks = gaussian_filter1d(kern(h / max(w, 1e-4)), max(sp, 1e-3) / DXF, mode='constant')
    return A * np.interp(x, XF, ks) + b0 + b1 * x


def gauss(x, A, x0, s, b0, b1):
    return A * np.exp(-0.5 * ((x - x0) / s) ** 2) + b0 + b1 * x


def sgauss(x, A, x0, sw, se, b0, b1):
    s = np.where(x < x0, sw, se)
    return A * np.exp(-0.5 * ((x - x0) / s) ** 2) + b0 + b1 * x


def fit_all(x, y, e):
    m = np.isfinite(y) & np.isfinite(e) & (np.abs(x) <= FITW) & (e > 0)
    x, y, e = x[m], y[m], e[m]
    A0 = np.nanmax(y) - np.nanmedian(y)
    b0 = np.nanmedian(y[np.abs(x) > 0.6])

    def run(f, p0, lo, hi):
        r = least_squares(lambda p: (f(x, *p) - y) / e, p0, bounds=(lo, hi))
        return r.x, float(np.sum(r.fun ** 2))
    res = {}
    pG, cG = run(gauss, [A0, 0.0, 0.12, b0, 0.0], [0, -0.4, 0.02, -50, -50], [50, 0.4, 0.5, 50, 50])
    res['G'] = dict(p=pG.tolist(), chi2=cG)
    pS, cS = run(sgauss, [pG[0], pG[1], pG[2], pG[2], pG[3], pG[4]], [0, -0.4, 0.02, 0.02, -50, -50],
                 [50, 0.4, 0.6, 0.6, 50, 50])
    res['SG'] = dict(p=pS.tolist(), chi2=cS, ratio_WE=float(pS[2] / pS[3]))
    for side, key in ((+1, 'SHE'), (-1, 'SHW')):
        best = None
        for w0 in (0.06, 0.15):
            for sp0 in (0.03, 0.08):
                xf0 = pG[1] + side * (0.5 * w0 + 0.05)
                p, c = run(lambda xx, A, xf, w, sp, b0_, b1_: sheet(xx, A, xf, w, sp, b0_, b1_, side),
                           [pG[0], xf0, w0, sp0, pG[3], pG[4]], [0, -0.8, 0.01, 0.005, -50, -50],
                           [100, 0.8, 0.8, 0.3, 50, 50])
                if best is None or c < best[1]:
                    best = (p, c)
        res[key] = dict(p=best[0].tolist(), chi2=best[1])
    res['dchi2_SHW_minus_SHE'] = res['SHW']['chi2'] - res['SHE']['chi2']
    res['dchi2_G_minus_SHE'] = res['G']['chi2'] - res['SHE']['chi2']
    return res


# ----------------------------------------------------------------------------- NSNS
M = np.load(DER / 'w3_figs_f1_maps.npz')
RAg, DEg = np.meshgrid(M['ra'], M['dec'])
RAv, DEv, Vv = RAg.ravel(), DEg.ravel(), M['oiii'].astype(float).ravel()
pre = np.isfinite(Vv) & (RAv > 109) & (RAv < 123.5) & (DEv > 9) & (DEv < 30)
RAv, DEv, Vv = RAv[pre], DEv[pre], np.clip(Vv[pre], -5, 15)
patch30 = MPath(TR['patch30'])
PATCH = patch30.contains_points(np.c_[RAv, DEv])


def block_profiles(cc, s0, s1, excl=None):
    x, s = cc.coords(radec_to_vec(RAv, DEv), maxdist_deg=1.8)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    if excl is not None:
        sel &= ~excl
    x, s, v = x[sel], s[sel], Vv[sel]
    nb = max(1, int(round((s1 - s0) / BLOCK)))
    bedge = np.linspace(s0, s1, nb + 1)
    bi = np.clip(np.digitize(s, bedge) - 1, 0, nb - 1)
    xi = np.digitize(x, EDGES) - 1
    P = np.full((nb, XC.size), np.nan)
    for b in range(nb):
        mb = bi == b
        for k in range(XC.size):
            vv = v[mb & (xi == k)]
            if vv.size >= 4:
                P[b, k] = np.median(vv)
    return P


def stack(P, idx=None):
    Q = P if idx is None else P[idx]
    with np.errstate(all='ignore'):
        return np.nanmean(Q, axis=0)


def shifted(dra, ddec=0.0):
    return CurveCoords(ra_sp + dra / np.cos(np.radians(dec_sp)), dec_sp + ddec)


xr, _ = CC.coords(radec_to_vec(RAv, DEv), maxdist_deg=2.0)
EXCL = np.isfinite(xr) & (np.abs(xr) < 1.0)
NULL_SHIFTS = [(a, b) for a in np.r_[-5:-1.9:0.5, 2.0:5.01:0.5] for b in (-2.0, 0.0, 2.0)]

out = {'nsns_north_subsegments': {}, 'galex_offset_null': {}}
for key, d0, d1, lab in (('N1', 20.0, 21.5, 'Dec 20-21.5, Halpha patch 30% contour masked'),
                         ('N2', 24.25, 28.0, 'Dec 24.25-28 (beyond W2)'),
                         ('N3', 21.5, 24.25, 'Dec 21.5-24.25 (W2 0.6-1.3 deg west)')):
    s0, s1 = s_at(d0), s_at(d1)
    P = block_profiles(CC, s0, s1, excl=PATCH)
    nb = P.shape[0]
    y = stack(P)
    boots = np.array([stack(P, RNG.integers(0, nb, nb)) for _ in range(150)])
    e = np.nanstd(boots, axis=0)
    e = np.where(np.isfinite(e) & (e > 0), e, np.nanmedian(e))
    F = fit_all(XC, y, e)
    bs = []
    for q in boots[:60]:
        Fb = fit_all(XC, q, e)
        bs.append((Fb['SG']['ratio_WE'], Fb['dchi2_SHW_minus_SHE']))
    bs = np.array(bs)
    pG = F['G']['p']
    inj = []
    for dra, ddec in NULL_SHIFTS:
        Pn = block_profiles(shifted(dra, ddec), s0, s1, excl=EXCL | PATCH)
        yn = stack(Pn)
        if np.isfinite(yn[np.abs(XC) <= FITW]).mean() < 0.85:
            continue
        Fn = fit_all(XC, yn + pG[0] * np.exp(-0.5 * ((XC - pG[1]) / pG[2]) ** 2), e)
        inj.append((Fn['SG']['ratio_WE'], Fn['dchi2_SHW_minus_SHE']))
    inj = np.array(inj)
    r_obs, d_obs = F['SG']['ratio_WE'], F['dchi2_SHW_minus_SHE']
    res = dict(label=lab, n_blocks=int(nb), ratio_WE=r_obs, sigma_W=F['SG']['p'][2], sigma_E=F['SG']['p'][3],
               dchi2_SHW_minus_SHE=d_obs, dchi2_G_minus_SHE=F['dchi2_G_minus_SHE'], G=F['G']['p'],
               boot_ratio_16_50_84=np.percentile(bs[:, 0], [16, 50, 84]).tolist(),
               boot_frac_SHE_pref=float(np.mean(bs[:, 1] > 0)),
               n_null=int(len(inj)),
               null_ratio_16_50_84=np.percentile(inj[:, 0], [16, 50, 84]).tolist(),
               p_ratio=float((np.sum(inj[:, 0] >= r_obs) + 1) / (len(inj) + 1)),
               p_dchi2=float((np.sum(inj[:, 1] >= d_obs) + 1) / (len(inj) + 1)))
    out['nsns_north_subsegments'][key] = res
    print(key, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res.items()}, flush=True)

# ----------------------------------------------------------------------------- GALEX centre-offset null
import verify_xrayfuv_galex_measure as GM  # noqa: E402
EG = np.arange(-1.2, 1.2001, 0.04)
XG = 0.5 * (EG[1:] + EG[:-1])
GSHIFTS = np.r_[-3.0:-0.79:0.2, 0.8:3.01:0.2]


def gal_rows(D, cc, s0, s1, excl=None, inject=None):
    x, s = cc.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.5)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    if excl is not None:
        sel &= ~excl
    rows = []
    for t in np.unique(D['tile'][sel]):
        m = sel & (D['tile'] == t)
        v = D['v'][m].copy()
        if inject is not None:
            v = v + inject(x[m])
        ref = np.abs(x[m]) > 0.3
        if ref.sum() < 20:
            continue
        v -= np.mean(v[ref])
        bi = np.digitize(x[m], EG) - 1
        ok = (bi >= 0) & (bi < XG.size)
        rows.append((np.bincount(bi[ok], minlength=XG.size), np.bincount(bi[ok], v[ok], XG.size)))
    return rows


def gstack(rows, idx=None):
    rr = rows if idx is None else [rows[i] for i in idx]
    cnt = np.sum([r[0] for r in rr], axis=0); sm = np.sum([r[1] for r in rr], axis=0)
    p = sm / np.maximum(cnt, 1)
    p[cnt < 20] = np.nan
    return p


def cfit(y, e, xmax=0.7, lin=True):
    m = np.isfinite(y) & np.isfinite(e) & (np.abs(XG) <= xmax) & (e > 0)
    f = (lambda p: (gauss(XG[m], *p) - y[m]) / e[m]) if lin else \
        (lambda p: (gauss(XG[m], p[0], p[1], p[2], p[3], 0.0) - y[m]) / e[m])
    p0 = [np.nanmax(y[m]), 0.0, 0.12, 0.0, 0.0] if lin else [np.nanmax(y[m]), 0.0, 0.12, 0.0]
    lo = [0, -0.4, 0.03, -500, -500] if lin else [0, -0.4, 0.03, -500]
    hi = [5000, 0.4, 0.4, 500, 500] if lin else [5000, 0.4, 0.4, 500]
    return least_squares(f, p0, bounds=(lo, hi)).x


def nsns_grid(s0, s1, excl=None):
    x, s = CC.coords(radec_to_vec(RAv, DEv), maxdist_deg=1.5)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    if excl is not None:
        sel &= ~excl
    bi = np.digitize(x[sel], EG) - 1
    ok = (bi >= 0) & (bi < XG.size)
    v = Vv[sel][ok]; bi = bi[ok]
    return np.array([np.median(v[bi == k]) if np.sum(bi == k) > 10 else np.nan for k in range(XG.size)])


for key, d0, d1 in (('S1', 12.0, 15.0), ('S12', 12.0, 17.6)):
    s0, s1 = s_at(d0), s_at(d1)
    yo = nsns_grid(s0, s1)
    eo = np.full_like(yo, np.nanstd(yo[np.abs(XG) > 0.6]) + 1e-3)
    var = {}
    for xmax in (0.5, 0.7, 0.9):
        for lin in (True, False):
            var[f'nsns_w{xmax}_lin{int(lin)}'] = float(cfit(yo, eo, xmax, lin)[1])
    po = cfit(yo, eo)
    # observed NSNS ridge shape (baseline removed), used as the co-spatial template
    tmpl_y = yo - (po[3] + po[4] * XG)
    tmpl_y = np.where(np.abs(XG) <= 0.7, tmpl_y, 0.0)
    tmpl_y = np.nan_to_num(tmpl_y) / max(po[0], 1e-6)
    res = dict(dec=[d0, d1], nsns_x0=float(po[1]), nsns_variants=var)
    for band in ('FUV', 'NUV'):
        D = GM.prep_band(band)
        rows = gal_rows(D, CC, s0, s1)
        yg = gstack(rows)
        bt = np.array([gstack(rows, RNG.integers(0, len(rows), len(rows))) for _ in range(200)])
        eg = np.nanstd(bt, axis=0)
        eg = np.where(np.isfinite(eg) & (eg > 0), eg, np.nanmedian(eg))
        pg = cfit(yg, eg)
        gvar = {}
        for xmax in (0.5, 0.7, 0.9):
            for lin in (True, False):
                gvar[f'w{xmax}_lin{int(lin)}'] = float(cfit(yg, eg, xmax, lin)[1] -
                                                        var[f'nsns_w{xmax}_lin{int(lin)}'])
        # null: co-spatial template (NSNS shape, GALEX amplitude) injected into off-Silk placements
        xg_all, _ = CC.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.5)
        exg = np.isfinite(xg_all) & (np.abs(xg_all) < 0.6)
        amp = float(pg[0])

        def inj(xx, amp=amp, po=po):
            return amp * np.interp(xx, XG - 0.0 + 0.0, tmpl_y, left=0.0, right=0.0)
        dxn = []
        for dra in GSHIFTS:
            rn = gal_rows(D, shifted(dra), s0, s1, excl=exg, inject=inj)
            if len(rn) < 5:
                continue
            yn = gstack(rn)
            try:
                dxn.append(float(cfit(yn, eg)[1] - po[1]))
            except Exception:
                pass
        dxn = np.array(dxn)
        dx = float(pg[1] - po[1])
        res[band] = dict(dx_obs=dx, x0=float(pg[1]), A=amp, window_baseline_variants_dx=gvar,
                         n_null=int(len(dxn)), null_dx_median=float(np.median(dxn)),
                         null_dx_16_84=np.percentile(dxn, [16, 84]).tolist(),
                         null_dx_std=float(np.std(dxn)),
                         p_null_ge_obs=float((np.sum(dxn >= dx) + 1) / (len(dxn) + 1)),
                         z_vs_null=float((dx - np.median(dxn)) / (np.std(dxn) if np.std(dxn) > 0 else np.nan)))
        print(key, band, {k: v for k, v in res[band].items()}, flush=True)
    out['galex_offset_null'][key] = res

(DER / 'w3_judge_obs_edge2.json').write_text(json.dumps(out, indent=1))
print('saved')
