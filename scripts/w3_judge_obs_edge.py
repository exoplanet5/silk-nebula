"""judge_obs: two untested, model-light geometric predictions of an edge-on shock sheet.

(1) Cross-ridge asymmetry of the NSNS DR0.2 [O III] profile. A shock sheet seen edge-on at its tangent has an
    abrupt outer edge (the front) and a 1/sqrt(h) tail on the post-shock side. If the driver is on the concave
    (west, Monogem) side, the front faces east and the tail points west. We fit:
      G   symmetric Gaussian + linear baseline;
      SG  split Gaussian (sigma_W, sigma_E) + linear baseline;
      SHE projected uniform-emissivity sheet, front on the EAST (tail west), Gaussian-smoothed;
      SHW the mirror image (front on the WEST, tail east).
    The kernel for depth w behind the front: K(u) = sqrt(u) (0<u<1), sqrt(u)-sqrt(u-1) (u>1), 0 (u<0), u = h/w.
(2) GALEX FUV/NUV ridge centre relative to the [O III] ridge (verify_shockmodel predicts FUV ~0.1 deg EAST).

Uncertainties: block bootstrap over along-spine blocks (0.25 deg) for NSNS; tile bootstrap for GALEX. A symmetric-
truth calibration injects the best symmetric Gaussian into off-Silk placements of the same estimator (real sky noise
and structure) and refits, giving the null distribution of sigma_W/sigma_E and of chi2(SHW)-chi2(SHE).
Offsets x > 0 = EAST (outside Monogem), from verify_xrayfuv_common.CurveCoords (checked there: P4 -> negative).

Inputs (project files only): w3_figs_f1_maps.npz (NSNS DR0.2 [O III], 1' CAR, Rayleigh),
  verify_xrayfuv_galex_blocks.npz via verify_xrayfuv_galex_measure.prep_band, silk_trace.json v1.0.
Outputs: data/derived/w3_judge_obs_edge.json, figures/w3_judge_obs_edge.png
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_obs_edge.py
"""
import sys
import json
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import CurveCoords, radec_to_vec, load_spine  # noqa: E402
from w3_figs_common import DER, FIG, load_trace  # noqa: E402
import silk_style as ss  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

RNG = np.random.default_rng(20261117)
ra_sp, dec_sp, _ = load_spine()
CC = CurveCoords(ra_sp, dec_sp)
TR = load_trace()


def s_at(dec):
    return float(np.interp(dec, CC.dec, CC.s))


SEGS = [('S1', 12.0, 15.0, 'Dec 12.0-15.0 (single E strand)'),
        ('S2', 15.0, 17.6, 'Dec 15.0-17.6 (W strand 0.1-0.7 deg west)'),
        ('N', 20.0, 28.0, 'Dec 20-28 (north; W2 at Dec 21.5-24)')]
BW = 0.02
EDGES = np.arange(-1.6, 1.6001, BW)
XC = 0.5 * (EDGES[1:] + EDGES[:-1])
FITW = 0.9          # fit window |x| <= FITW (deg)
BLOCK = 0.25        # along-spine block length (deg)


# ----------------------------------------------------------------------------- models
def kern(u):
    out = np.zeros_like(u)
    a = (u > 0) & (u <= 1)
    out[a] = np.sqrt(u[a])
    b = u > 1
    out[b] = np.sqrt(u[b]) - np.sqrt(u[b] - 1)
    return out


DXF = 0.004
XF = np.arange(-2.5, 2.5001, DXF)
from scipy.ndimage import gaussian_filter1d  # noqa: E402


def sheet(x, A, xf, w, sp, b0, b1, side):
    """side=+1: front at xf on the EAST, emission at x < xf (tail west). side=-1: mirror.
    The kernel is normalised to unit peak before smoothing, so A is the unsmoothed peak."""
    h = side * (xf - XF)
    k = kern(h / max(w, 1e-4))
    ks = gaussian_filter1d(k, max(sp, 1e-3) / DXF, mode='constant')
    return A * np.interp(x, XF, ks) + b0 + b1 * x


def gauss(x, A, x0, s, b0, b1):
    return A * np.exp(-0.5 * ((x - x0) / s) ** 2) + b0 + b1 * x


def sgauss(x, A, x0, sw, se, b0, b1):
    s = np.where(x < x0, sw, se)
    return A * np.exp(-0.5 * ((x - x0) / s) ** 2) + b0 + b1 * x


def fit_all(x, y, e):
    m = np.isfinite(y) & np.isfinite(e) & (np.abs(x) <= FITW) & (e > 0)
    x, y, e = x[m], y[m], e[m]
    res = {}
    A0 = np.nanmax(y) - np.nanmedian(y)
    b0 = np.nanmedian(y[np.abs(x) > 0.6]) if np.any(np.abs(x) > 0.6) else 0.0

    def run(f, p0, lo, hi):
        r = least_squares(lambda p: (f(x, *p) - y) / e, p0, bounds=(lo, hi))
        return r.x, float(np.sum(r.fun ** 2))
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
                try:
                    p, c = run(lambda xx, A, xf, w, sp, b0_, b1_: sheet(xx, A, xf, w, sp, b0_, b1_, side),
                               [pG[0], xf0, w0, sp0, pG[3], pG[4]],
                               [0, -0.8, 0.01, 0.005, -50, -50], [100, 0.8, 0.8, 0.3, 50, 50])
                except Exception:
                    continue
                if best is None or c < best[1]:
                    best = (p, c)
        res[key] = dict(p=best[0].tolist(), chi2=best[1])
    res['dchi2_SHW_minus_SHE'] = res['SHW']['chi2'] - res['SHE']['chi2']
    res['dchi2_G_minus_SHE'] = res['G']['chi2'] - res['SHE']['chi2']
    res['n_bins'] = int(m.sum())
    return res


# ----------------------------------------------------------------------------- NSNS profiles
M = np.load(DER / 'w3_figs_f1_maps.npz')
O3 = M['oiii'].astype(float)
RAg, DEg = np.meshgrid(M['ra'], M['dec'])
RAv, DEv, Vv = RAg.ravel(), DEg.ravel(), O3.ravel()
pre = np.isfinite(Vv) & (RAv > 109) & (RAv < 123.5) & (DEv > 9) & (DEv < 30)
RAv, DEv, Vv = RAv[pre], DEv[pre], Vv[pre]
Vv = np.clip(Vv, -5, 15)   # bright-star residuals; medians below make this cosmetic


def block_profiles(cc, ra, dec, val, s0, s1, excl=None):
    x, s = cc.coords(radec_to_vec(ra, dec), maxdist_deg=1.8)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    if excl is not None:
        sel &= ~excl
    x, s, v = x[sel], s[sel], val[sel]
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

out = {'segments': {}, 'notes': {}}
prof_store = {}
for key, d0, d1, lab in SEGS:
    s0, s1 = s_at(d0), s_at(d1)
    P = block_profiles(CC, RAv, DEv, Vv, s0, s1)
    nb = P.shape[0]
    y = stack(P)
    boots = []
    for _ in range(200):
        idx = RNG.integers(0, nb, nb)
        boots.append(stack(P, idx))
    boots = np.array(boots)
    e = np.nanstd(boots, axis=0)
    e = np.where(np.isfinite(e) & (e > 0), e, np.nanmedian(e))
    F = fit_all(XC, y, e)
    # bootstrap distribution of the asymmetry statistics
    bstat = []
    for q in boots[:100]:
        try:
            Fb = fit_all(XC, q, e)
            bstat.append((Fb['SG']['ratio_WE'], Fb['dchi2_SHW_minus_SHE'], Fb['G']['p'][1], Fb['dchi2_G_minus_SHE']))
        except Exception:
            pass
    bstat = np.array(bstat)
    # symmetric-truth calibration: inject the best symmetric Gaussian into off-Silk placements
    pG = F['G']['p']
    inj = []
    injS = []
    pE_best = F['SHE']['p']
    for dra, ddec in NULL_SHIFTS:
        ccn = shifted(dra, ddec)
        Pn = block_profiles(ccn, RAv, DEv, Vv, s0, s1, excl=EXCL)
        yn = stack(Pn)
        if np.isfinite(yn[np.abs(XC) <= FITW]).mean() < 0.85:
            continue
        ysym = yn + pG[0] * np.exp(-0.5 * ((XC - pG[1]) / pG[2]) ** 2)
        ysh = yn + sheet(XC, *pE_best[:4], 0.0, 0.0, +1)
        try:
            Fn = fit_all(XC, ysym, e)
            inj.append((Fn['SG']['ratio_WE'], Fn['dchi2_SHW_minus_SHE'], Fn['dchi2_G_minus_SHE']))
            Fs = fit_all(XC, ysh, e)
            injS.append((Fs['SG']['ratio_WE'], Fs['dchi2_SHW_minus_SHE'], Fs['dchi2_G_minus_SHE']))
        except Exception:
            pass
    inj = np.array(inj)
    injS = np.array(injS)
    r_obs, d_obs = F['SG']['ratio_WE'], F['dchi2_SHW_minus_SHE']
    seg = dict(label=lab, dec=[d0, d1], n_blocks=int(nb), fit=F,
               boot=dict(n=int(len(bstat)),
                         ratio_WE_16_50_84=np.percentile(bstat[:, 0], [16, 50, 84]).tolist(),
                         dchi2_16_50_84=np.percentile(bstat[:, 1], [16, 50, 84]).tolist(),
                         frac_dchi2_pos=float(np.mean(bstat[:, 1] > 0)),
                         x0_16_50_84=np.percentile(bstat[:, 2], [16, 50, 84]).tolist(),
                         dchi2_G_minus_SHE_16_50_84=np.percentile(bstat[:, 3], [16, 50, 84]).tolist()),
               sym_null=dict(n=int(len(inj)),
                             ratio_WE_16_50_84=(np.percentile(inj[:, 0], [16, 50, 84]).tolist() if len(inj) else None),
                             dchi2_16_50_84=(np.percentile(inj[:, 1], [16, 50, 84]).tolist() if len(inj) else None),
                             p_ratio_ge_obs=(float((np.sum(inj[:, 0] >= r_obs) + 1) / (len(inj) + 1)) if len(inj) else None),
                             p_dchi2_ge_obs=(float((np.sum(inj[:, 1] >= d_obs) + 1) / (len(inj) + 1)) if len(inj) else None),
                             dchi2_G_minus_SHE_16_50_84=np.percentile(inj[:, 2], [16, 50, 84]).tolist(),
                             p_GmSHE_le_obs=float((np.sum(inj[:, 2] <= F['dchi2_G_minus_SHE']) + 1) / (len(inj) + 1))),
               sheet_truth=dict(n=int(len(injS)),
                                ratio_WE_16_50_84=np.percentile(injS[:, 0], [16, 50, 84]).tolist(),
                                dchi2_SHW_minus_SHE_16_50_84=np.percentile(injS[:, 1], [16, 50, 84]).tolist(),
                                frac_SHE_preferred=float(np.mean(injS[:, 1] > 0)),
                                dchi2_G_minus_SHE_16_50_84=np.percentile(injS[:, 2], [16, 50, 84]).tolist(),
                                p_GmSHE_le_obs=float((np.sum(injS[:, 2] <= F['dchi2_G_minus_SHE']) + 1) / (len(injS) + 1))))
    out['segments'][key] = seg
    prof_store[key] = (y, e, F)
    print(key, lab, 'nb', nb, 'ratio W/E %.2f' % r_obs, 'dchi2(SHW-SHE) %.1f' % d_obs,
          'dchi2(G-SHE) %.1f' % F['dchi2_G_minus_SHE'], 'boot', seg['boot'], 'null', seg['sym_null'], flush=True)

# ----------------------------------------------------------------------------- GALEX ridge centres
import verify_xrayfuv_galex_measure as GM  # noqa: E402
EG = np.arange(-1.2, 1.2001, 0.04)
XG = 0.5 * (EG[1:] + EG[:-1])


def gal_tile_profiles(D, s0, s1):
    x, s = CC.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.5)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    tiles = np.unique(D['tile'][sel])
    rows = []
    for t in tiles:
        m = sel & (D['tile'] == t)
        ref = m & (np.abs(x) > 0.3)
        if ref.sum() < 20:
            continue
        v = D['v'][m] - np.mean(D['v'][ref])
        bi = np.digitize(x[m], EG) - 1
        ok = (bi >= 0) & (bi < XG.size)
        cnt = np.bincount(bi[ok], minlength=XG.size)
        sm = np.bincount(bi[ok], v[ok], XG.size)
        rows.append((cnt, sm))
    return rows


def gal_stack(rows, idx=None):
    rr = rows if idx is None else [rows[i] for i in idx]
    cnt = np.sum([r[0] for r in rr], axis=0)
    sm = np.sum([r[1] for r in rr], axis=0)
    p = sm / np.maximum(cnt, 1)
    p[cnt < 20] = np.nan
    return p


def centre_fit(x, y, e, xmax=0.7):
    m = np.isfinite(y) & np.isfinite(e) & (np.abs(x) <= xmax) & (e > 0)
    r = least_squares(lambda p: (gauss(x[m], *p) - y[m]) / e[m], [np.nanmax(y[m]), 0.0, 0.12, 0.0, 0.0],
                      bounds=([0, -0.4, 0.03, -500, -500], [5000, 0.4, 0.4, 500, 500]))
    return r.x


# NSNS centre on the same binning, for a like-for-like comparison
def nsns_on_grid(s0, s1):
    x, s = CC.coords(radec_to_vec(RAv, DEv), maxdist_deg=1.5)
    sel = np.isfinite(x) & (s >= s0) & (s < s1)
    bi = np.digitize(x[sel], EG) - 1
    ok = (bi >= 0) & (bi < XG.size)
    v = Vv[sel][ok]; bi = bi[ok]
    return np.array([np.median(v[bi == k]) if np.sum(bi == k) > 10 else np.nan for k in range(XG.size)])


out['galex'] = {}
gal_store = {}
for key, d0, d1 in (('S1', 12.0, 15.0), ('S12', 12.0, 17.6)):
    s0, s1 = s_at(d0), s_at(d1)
    yo = nsns_on_grid(s0, s1)
    eo = np.full_like(yo, np.nanstd(yo[np.abs(XG) > 0.6]) + 1e-3)
    po = centre_fit(XG, yo, eo)
    res = dict(dec=[d0, d1], nsns_x0=float(po[1]), nsns_sigma=float(po[2]))
    for band in ('FUV', 'NUV'):
        D = GM.prep_band(band)
        rows = gal_tile_profiles(D, s0, s1)
        yg = gal_stack(rows)
        bt = np.array([gal_stack(rows, RNG.integers(0, len(rows), len(rows))) for _ in range(300)])
        eg = np.nanstd(bt, axis=0)
        eg = np.where(np.isfinite(eg) & (eg > 0), eg, np.nanmedian(eg))
        pg = centre_fit(XG, yg, eg)
        bx = []
        for q in bt[:200]:
            try:
                bx.append(centre_fit(XG, q, eg)[1])
            except Exception:
                pass
        bx = np.array(bx)
        res[band] = dict(n_tiles=len(rows), A=float(pg[0]), x0=float(pg[1]), sigma=float(pg[2]),
                         x0_boot_16_50_84=np.percentile(bx, [16, 50, 84]).tolist(),
                         dx_vs_nsns=float(pg[1] - po[1]),
                         dx_sigma=float(0.5 * (np.percentile(bx, 84) - np.percentile(bx, 16))))
        gal_store[(key, band)] = (yg, eg, pg)
        print(key, band, res[band], flush=True)
    gal_store[(key, 'NSNS')] = (yo, eo, po)
    out['galex'][key] = res

out['notes'] = {
    'prediction_H1_H2_west_driver': 'front on the east (SHE preferred: dchi2(SHW-SHE) > 0; split-Gaussian W/E > 1 at '
                                    'low levels); FUV ridge ~ +0.1 deg (east) of [O III] (verify_shockmodel)',
    'caveats': 'NSNS DR0.2 PSF and processing unknown (1-arcmin HiPS order 3); the median stack along a wiggly spine '
               'broadens symmetrically; the spine itself was traced from the same images; W strand contaminates the '
               'west wing in S2 and W2 in N; linear baselines can absorb part of a long tail.',
}
(DER / 'w3_judge_obs_edge.json').write_text(json.dumps(out, indent=1))

# ----------------------------------------------------------------------------- figure
fig, axes = ss.new_sky_figure(13, 7.6, 'Silk Nebula · Edge-on shock tests',
                              'Cross-ridge profiles; + offset = east (away from Monogem)', nrows=2, ncols=3)
for j, (key, d0, d1, lab) in enumerate(SEGS):
    ax = axes[0, j]
    y, e, F = prof_store[key]
    ax.errorbar(XC, y, e, fmt='o', ms=2.2, color=ss.OIII_C, ecolor='#2f5f78', lw=0.6, label='NSNS [O III]')
    xx = np.linspace(-FITW, FITW, 600)
    ax.plot(xx, gauss(xx, *F['G']['p']), color=ss.FGDIM, lw=1.0, ls='-', label='Gaussian')
    pE = F['SHE']['p']; pW = F['SHW']['p']
    ax.plot(xx, sheet(xx, *pE, +1), color=ss.OK, lw=1.3, label='sheet, front east')
    ax.plot(xx, sheet(xx, *pW, -1), color=ss.WARN, lw=1.0, ls='--', label='sheet, front west')
    seg = out['segments'][key]
    txt = ('W/E width %.2f\nΔχ² (W−E front) %+.1f\nboot P(>0) %.2f' %
           (F['SG']['ratio_WE'], F['dchi2_SHW_minus_SHE'], seg['boot']['frac_dchi2_pos']))
    ss.label(ax, 0.03, 0.95, txt, size=7.5, transform=ax.transAxes, va='top')
    ax.set_title(lab, fontsize=8.5, loc='left')
    ax.set_xlim(-1.2, 1.2)
    ax.set_xlabel('offset from E spine (deg)')
    ax.set_ylabel('R')
    ax.axvline(0, color='#2a323c', lw=0.8)
    if j == 0:
        ax.legend(loc='upper right', fontsize=6.5)
for j, key in enumerate(('S1', 'S12')):
    ax = axes[1, j]
    for band, col in (('NSNS', ss.OIII_C), ('FUV', '#b388ff'), ('NUV', ss.WARN)):
        yb, eb, pb = gal_store[(key, band)]
        base = pb[3] + pb[4] * XG
        nrm = pb[0]
        ax.errorbar(XG, (yb - base) / nrm, eb / nrm, fmt='o-', ms=2.5, lw=0.8, color=col,
                    label='%s  x0 %+.3f' % (band, pb[1]))
    ax.axvline(0, color='#2a323c', lw=0.8)
    ax.axvline(0.1, color=ss.OK, lw=0.8, ls=':')
    ss.label(ax, 0.105, 1.05, 'predicted FUV', size=7, color=ss.OK)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.6, 1.4)
    ax.set_title('Dec %.1f-%.1f: ridge-normalised' % tuple(out['galex'][key]['dec']), fontsize=8.5, loc='left')
    ax.set_xlabel('offset from E spine (deg)')
    ax.legend(loc='upper left', fontsize=6.5)
ax = axes[1, 2]
ax.axis('off')
lines = ['GALEX ridge offset vs [O III] (deg, + = east)']
for key in ('S1', 'S12'):
    g = out['galex'][key]
    for band in ('FUV', 'NUV'):
        lines.append('%s %s: %+.3f ± %.3f' % (key, band, g[band]['dx_vs_nsns'], g[band]['dx_sigma']))
lines.append('')
lines.append('Shock prediction: FUV +0.1, sharp east edge')
ss.label(ax, 0.02, 0.9, '\n'.join(lines), size=8, transform=ax.transAxes, va='top')
ss.finish(fig, str(FIG / 'w3_judge_obs_edge.png'))
print('saved')
