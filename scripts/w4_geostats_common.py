"""w4_geostats shared helpers: one set of loaders, cross-profile estimators, placements and injections used by
w4_geostats_registration.py, w4_geostats_order.py and w4_geostats_shape.py.

Conventions (identical to verify_xrayfuv_common / w3_figs F6 / w3_judge_obs / w3_judge_skeptic):
  x = signed perpendicular offset from the E spine (silk_trace.json v1.0), deg, + = EAST (outside, away from
      Monogem); s = along-spine arc length from the south end (Dec 12.0).
  Segments: S1 Dec 12.0-15.0 (single E strand, the only clean section), S12 Dec 12.0-17.6 (includes the W
      strand 0.1-0.7 deg west), N Dec 20-28 (W2 strand 0.6-1.3 deg west at Dec 21.5-24; Halpha patch masked),
      N2 Dec 24.25-28 (beyond W2).
Profile estimator ("stack"): per along-spine block of 0.25 deg and per x bin, the median (NSNS) or mean (GALEX,
  DESI, eRASS) of the pixel/fibre values; the stack is the mean over blocks. GALEX values are first de-meaned per
  tile against |x| > 0.3 deg (as F6 / judge_skeptic).
Inputs are project files only (read-only). Downloaded data are never executed; run scripts with python -I.
"""
import sys
import json
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import (CurveCoords, radec_to_vec, vec_to_radec, load_spine, CENTRES, sep_pa,  # noqa
                                   offset_point, rotate_about, DER, ROOT)

FIG = ROOT / 'figures'
ra_sp, dec_sp, TRACE = load_spine()
CC0 = CurveCoords(ra_sp, dec_sp)
HA_PATCH = (115.474, 20.536)


def s_at(dec, cc=CC0):
    return float(np.interp(dec, cc.dec, cc.s))


SEGMENTS = {'S1': (12.0, 15.0), 'S12': (12.0, 17.6), 'N': (20.0, 28.0), 'N2': (24.25, 28.0)}
BLOCK = 0.25


def json_dump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    path.write_text(json.dumps(obj, indent=1, default=conv))


# ------------------------------------------------------------------------------------------------ placements
def shifted_cc(dra, ddec=0.0):
    """Rigid-ish shift of the spine (RA shift at constant angular offset, plus a Dec shift)."""
    return CurveCoords(ra_sp + dra / np.cos(np.radians(dec_sp)), dec_sp + ddec)


def rotated_cc(theta, centre='PSR'):
    r, d = rotate_about(ra_sp, dec_sp, *CENTRES[centre], theta)
    return CurveCoords(r, d)


def near_real(ra, dec, dist):
    x, s = CC0.coords(radec_to_vec(ra, dec), maxdist_deg=dist + 0.5)
    return np.isfinite(x) & (np.abs(x) < dist)


# ------------------------------------------------------------------------------------------------ data
def load_nsns(box=(109.0, 123.5, 9.0, 30.0)):
    M = np.load(DER / 'w3_figs_f1_maps.npz')
    o3 = M['oiii'].astype(float)
    RA, DE = np.meshgrid(M['ra'], M['dec'])
    RA, DE, V = RA.ravel(), DE.ravel(), o3.ravel()
    m = np.isfinite(V) & (RA > box[0]) & (RA < box[1]) & (DE > box[2]) & (DE < box[3])
    return dict(ra=RA[m], dec=DE[m], v=np.clip(V[m], -5, 15), tile=None)


def load_galex(band):
    import verify_xrayfuv_galex_measure as GM   # module-level code reads only project files
    D = GM.prep_band(band)
    return dict(ra=D['ra'], dec=D['dec'], v=D['v'].astype(float), tile=D['tile'])


def load_desi():
    from verify_desi_common import good_fit
    lf = np.load(DER / 'verify_desi_linefluxes.npz')
    cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    det = json.loads((DER / 'verify_desi_detect.json').read_text())
    sc = det['err_scale_from_nulls']
    ok3 = good_fit(lf, 'o3', scale=sc) & cm['ok_contrib'] & np.isfinite(cm['kR_5008'])
    ok2 = good_fit(lf, 'o2sum', scale=sc) & cm['ok_contrib'] & np.isfinite(cm['kR_3727'])
    ok = ok3 & ok2
    return dict(ra=cm['ra'][ok], dec=cm['dec'][ok], o3=(lf['o3__amp'][:, 0] * cm['kR_5008'])[ok],
                o2=(lf['o2sum__amp'][:, 0] * cm['kR_3727'])[ok],
                e3=(lf['o3__err'][:, 0] * sc * cm['kR_5008'])[ok], e2=(lf['o2sum__err'][:, 0] * sc * cm['kR_3727'])[ok],
                petal=cm['group_petal'][ok], n_all=int(len(ok)), n_ok=int(ok.sum()), err_scale=float(sc))


def load_erass():
    z = np.load(DER / 'verify_xrayfuv_erass_pixels.npz')
    soft = (z['b0204'] + z['b0406']).astype(float)
    m = np.isfinite(soft) & ~z['catmask']
    return dict(ra=z['ra'][m].astype(float), dec=z['dec'][m].astype(float), v=soft[m],
                gl=z['gl'][m].astype(float), gb=z['gb'][m].astype(float), nhi=z['nhi'][m].astype(float))


# ------------------------------------------------------------------------------------------------ profiles
class Placed:
    """Pixel/fibre coordinates (x, s, block) relative to one placed curve, cached for repeated injections."""
    def __init__(self, cc, D, seg, maxd=1.8, excl_real=None, mask=None, block=BLOCK):
        d0, d1 = SEGMENTS[seg] if isinstance(seg, str) else seg
        s0, s1 = s_at(d0), s_at(d1)      # arc-length range of the REAL segment (placed curves are rigid copies)
        pre = np.ones(D['ra'].size, bool)
        # cheap pre-cut on a box around the placed segment
        r_seg = cc.ra[(cc.s >= s0) & (cc.s <= s1)]; d_seg = cc.dec[(cc.s >= s0) & (cc.s <= s1)]
        cd = np.cos(np.radians(np.mean(d_seg)))
        pre &= (D['ra'] > r_seg.min() - (maxd + 0.2) / cd) & (D['ra'] < r_seg.max() + (maxd + 0.2) / cd)
        pre &= (D['dec'] > d_seg.min() - maxd - 0.2) & (D['dec'] < d_seg.max() + maxd + 0.2)
        idx = np.nonzero(pre)[0]
        x, s = cc.coords(radec_to_vec(D['ra'][idx], D['dec'][idx]), maxdist_deg=maxd)
        sel = np.isfinite(x) & (s >= s0) & (s < s1)
        if excl_real is not None:
            sel &= ~near_real(D['ra'][idx], D['dec'][idx], excl_real)
        if mask is not None:
            sel &= ~mask[idx]
        self.idx = idx[sel]; self.x = x[sel]; self.s = s[sel]
        nb = max(1, int(round((s1 - s0) / block)))
        self.nblock = nb
        self.blk = np.clip(((self.s - s0) / (s1 - s0) * nb).astype(int), 0, nb - 1)
        self.tile = None if D.get('tile') is None else D['tile'][self.idx]


def stack_profile(P, vals, edges, stat='median', tile_demean=False, min_n=4, demean_out=0.3):
    """vals: values at P.idx (already injected if desired). Returns the block-mean stacked profile."""
    v = vals.astype(float).copy()
    x = P.x
    if tile_demean and P.tile is not None:
        keep = np.ones(v.size, bool)
        ut, inv = np.unique(P.tile, return_inverse=True)
        ref = np.abs(x) > demean_out
        sw = np.bincount(inv[ref], weights=v[ref], minlength=ut.size)
        nn = np.bincount(inv[ref], minlength=ut.size)
        good = nn > 20
        v = v - np.where(good, sw / np.maximum(nn, 1), 0)[inv]
        keep &= good[inv]
        v = np.where(keep, v, np.nan)
    nbin = len(edges) - 1
    xi = np.digitize(x, edges) - 1
    ok = (xi >= 0) & (xi < nbin) & np.isfinite(v)
    key = P.blk[ok] * nbin + xi[ok]
    vv = v[ok]
    nk = P.nblock * nbin
    if stat == 'mean':
        sw = np.bincount(key, weights=vv, minlength=nk)
        n = np.bincount(key, minlength=nk)
        prof = np.where(n >= min_n, sw / np.maximum(n, 1), np.nan)
    else:
        order = np.lexsort((vv, key))
        ks = key[order]; vs = vv[order]
        n = np.bincount(ks, minlength=nk)
        start = np.concatenate([[0], np.cumsum(n)[:-1]])
        lo = start + (n - 1) // 2; hi = start + n // 2
        med = np.full(nk, np.nan)
        g = n >= min_n
        med[g] = 0.5 * (vs[lo[g]] + vs[hi[g]])
        prof = med
    B = prof.reshape(P.nblock, nbin)
    with np.errstate(all='ignore'):
        st = np.nanmean(B, axis=0)
        frac = np.mean(np.isfinite(B), axis=0)
    st[frac < 0.5] = np.nan
    return st, B


# ------------------------------------------------------------------------------------------------ estimators
def splitnorm(p, x):
    A, x0, sw, se, b0, b1 = p
    sig = np.where(x < x0, sw, se)
    return A * np.exp(-0.5 * ((x - x0) / sig) ** 2) + b0 + b1 * x


def fit_split(x, y, xmax=0.9, w=None):
    m = np.isfinite(y) & (np.abs(x) <= xmax)
    if m.sum() < 12:
        return None
    ww = np.ones(m.sum()) if w is None else w[m]
    yy = y[m]; xx = x[m]
    base = np.nanmedian(yy[np.abs(xx) > 0.5]) if np.any(np.abs(xx) > 0.5) else np.nanmedian(yy)
    p0 = [max(np.nanmax(yy) - base, 1e-3), 0.0, 0.15, 0.15, base, 0.0]
    lb = [0, -0.3, 0.02, 0.02, -np.inf, -np.inf]; ub = [np.inf, 0.3, 0.8, 0.8, np.inf, np.inf]
    try:
        r = least_squares(lambda p: (splitnorm(p, xx) - yy) * ww, p0, bounds=(lb, ub))
    except Exception:
        return None
    A, x0, sw, se, b0, b1 = r.x
    return dict(A=float(A), x0=float(x0), sig_w=float(sw), sig_e=float(se), ratio_we=float(sw / se),
                fwhm=float(1.1774 * (sw + se)), b0=float(b0), b1=float(b1))


def fit_gauss_centre(x, y, xmax=0.9):
    """Free-centre, free-width Gaussian + linear baseline; returns the centre."""
    m = np.isfinite(y) & (np.abs(x) <= xmax)
    if m.sum() < 10:
        return np.nan
    xx, yy = x[m], y[m]
    base = np.nanmedian(yy[np.abs(xx) > 0.5]) if np.any(np.abs(xx) > 0.5) else np.nanmedian(yy)
    f = lambda p: p[0] * np.exp(-0.5 * ((xx - p[1]) / p[2]) ** 2) + p[3] + p[4] * xx - yy
    try:
        r = least_squares(f, [max(np.nanmax(yy) - base, 1e-3), 0.0, 0.15, base, 0.0],
                          bounds=([0, -0.35, 0.03, -np.inf, -np.inf], [np.inf, 0.35, 0.6, np.inf, np.inf]))
        return float(r.x[1])
    except Exception:
        return np.nan


DGRID = np.arange(-0.30, 0.3001, 0.005)


def shift_fit(x, y, tmpl_fun, xmax=0.9):
    """y = a T(x - d) + b0 + b1 x, grid over d (5 mdeg) + parabolic refinement. Returns (d, a)."""
    m = np.isfinite(y) & (np.abs(x) <= xmax)
    if m.sum() < 10:
        return np.nan, np.nan
    xx, yy = x[m], y[m]
    chis = np.empty(DGRID.size); amps = np.empty(DGRID.size)
    for i, dl in enumerate(DGRID):
        X = np.vstack([tmpl_fun(xx - dl), np.ones(xx.size), xx]).T
        coef, *_ = np.linalg.lstsq(X, yy, rcond=None)
        chis[i] = np.sum((X @ coef - yy) ** 2); amps[i] = coef[0]
    k = int(np.argmin(chis))
    d = DGRID[k]
    if 0 < k < DGRID.size - 1:
        c0, c1, c2 = chis[k - 1], chis[k], chis[k + 1]
        den = c0 - 2 * c1 + c2
        if den > 0:
            d = d + 0.5 * (c0 - c2) / den * (DGRID[1] - DGRID[0])
    return float(d), float(amps[k])


def moment_centre(x, y, base_fit, half=0.4):
    if base_fit is None:
        return np.nan
    m = np.isfinite(y) & (np.abs(x - base_fit['x0']) <= half)
    yy = y[m] - (base_fit['b0'] + base_fit['b1'] * x[m])
    return float(np.sum(x[m] * yy) / np.sum(yy)) if np.sum(yy) > 0 else np.nan


def template_from_fit(f):
    """Baseline-free split-normal template (unit peak) centred at 0 -> T(x)."""
    sw, se = f['sig_w'], f['sig_e']
    return lambda xx: np.exp(-0.5 * (xx / np.where(xx < 0, sw, se)) ** 2)


def empirical_p(null, obs, side='greater'):
    null = np.asarray(null)[np.isfinite(null)]
    n = null.size
    if n == 0:
        return np.nan
    if side == 'greater':
        k = np.sum(null >= obs)
    elif side == 'less':
        k = np.sum(null <= obs)
    else:
        k = np.sum(np.abs(null) >= abs(obs))
    return float((k + 1) / (n + 1))


def robust_sd(a):
    a = np.asarray(a)[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a)))) if a.size else np.nan
