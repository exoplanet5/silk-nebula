"""verify_xrayfuv: independent test of the claim that the Silk spine marks the outer edge of the
Monogem soft X-ray emission (eRASS1 0.2-0.6 keV; ROSAT R1+R2, R4+R5).

Statistics on the cross-spine profile I(x) (x = signed perpendicular offset, x<0 = west = toward
the Monogem centres), x in [-2.5, +2.5] deg, 0.1-deg bins (eRASS) / 0.2-deg (ROSAT):
  WE    = mean I(-2.5..-0.5) / mean I(+0.5..+2.5)       (wave-1 statistic)
  c     = step at x=0 from I = a + b x + c [x<0]  (step beyond a linear gradient)
  cfrac = c / a
  xbest = location of the best-fitting step (scan -2..+2 deg)
Nulls: (R) rigid rotations of the spine about the Monogem centre (same radius, other PA);
       (T) rigid translations over the surrounding sky; (S) radial scaling scan r -> r + dr.
Placements overlapping the real spine (>10% of band pixels within 1.5 deg of it), with <50%
valid pixels, or reaching l < 180.3 (outside eROSITA-DE) are rejected.
Variants: masks (none / catalogue / catalogue+clip); absorption correction exp(sigma N_HI)
with sigma = 0, 5e-21, 1e-20 cm^2 (f=1, an upper bound); smooth-shell residual maps.
Outputs: data/derived/verify_xrayfuv_xray_edge.json, verify_xrayfuv_xray_nulls.npz
"""
import sys, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import *

rng = np.random.default_rng(20261116)
XB = 2.5

# ------------------------------------------------------------------ data
def load_erass():
    z = np.load(DER / 'verify_xrayfuv_erass_pixels.npz')
    d = {k: z[k] for k in z.files}
    d['soft'] = (d['b0204'] + d['b0406']).astype(float)
    d['mid'] = d['b0610'].astype(float); d['hard'] = d['b1023'].astype(float)
    d['mask_none'] = np.isfinite(d['soft'])
    d['mask_cat'] = d['mask_none'] & ~d['catmask']
    d['mask_full'] = d['good'].astype(bool)
    d['vec'] = radec_to_vec(d['ra'].astype(float), d['dec'].astype(float))
    d['dx'] = 0.1
    return d


def load_rosat():
    S = ROOT / 'data/raw/verify_xrayfuv_skyview'
    hd = fits.open(S / 'wide_rass_r1_l203b13_car.fits')[0]
    w = WCS(hd.header)
    ny, nx = hd.data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    gl, gb = w.pixel_to_world_values(xx, yy)
    gl = gl % 360
    r = {k: fits.getdata(S / f'wide_rass_{k}_l203b13_car.fits').astype(float) for k in ['r1', 'r2', 'r4', 'r5']}
    nh = fits.getdata(S / 'wide_hi4pi_nhi_l203b13_car.fits').astype(float)
    c = SkyCoord(gl.ravel() * u.deg, gb.ravel() * u.deg, frame='galactic').icrs
    d = {'ra': c.ra.deg, 'dec': c.dec.deg, 'gl': gl.ravel(), 'gb': gb.ravel(), 'nhi': nh.ravel(),
         'r12': (r['r1'] + r['r2']).ravel(), 'r45': (r['r4'] + r['r5']).ravel()}
    fin = np.isfinite(d['ra']) & np.isfinite(d['dec'])
    d = {k: v[fin] for k, v in d.items()}
    d['mask_none'] = np.isfinite(d['r12']) & (d['r12'] > 0)
    d['vec'] = radec_to_vec(d['ra'], d['dec'])
    d['dx'] = 0.2
    return d

# ------------------------------------------------------------------ profile statistics
class Placer:
    def __init__(self, data):
        self.d = data
        self.tree = cKDTree(data['vec'])

    def profile(self, ra_c, dec_c, img, mask, dx, smin=None, smax=None):
        cc = CurveCoords(ra_c, dec_c, step=0.02)
        # candidate pixels: within 3.2 deg of any curve point (query the pixel tree with curve pts)
        idx = self.tree.query_ball_point(cc.v[::10], r=2 * np.sin(np.radians(3.2) / 2))
        cand = np.unique(np.concatenate([np.asarray(i, int) for i in idx])) if len(idx) else np.array([], int)
        if cand.size == 0:
            return None
        x, s = cc.coords(self.d['vec'][cand], maxdist_deg=3.2)
        ok = np.isfinite(x) & np.isfinite(img[cand])
        if smin is not None:
            ok &= (s >= smin) & (s <= smax)
        nb = int(round(2 * 3.0 / dx))
        edges = np.linspace(-3.0, 3.0, nb + 1)
        xc = 0.5 * (edges[1:] + edges[:-1])
        bi = np.digitize(x, edges) - 1
        ok_all = ok & (bi >= 0) & (bi < nb)
        okm = ok_all & mask[cand]
        sw = np.bincount(bi[okm], weights=img[cand][okm], minlength=nb)
        n = np.bincount(bi[okm], minlength=nb).astype(float)
        nall = np.bincount(bi[ok_all], minlength=nb).astype(float)
        prof = np.where(n > 0, sw / np.maximum(n, 1), np.nan)
        return xc, prof, n, nall, cand[okm], x[okm]


def stats(xc, prof, xb=XB):
    sel = (np.abs(xc) <= xb) & np.isfinite(prof)
    if sel.sum() < 0.7 * (np.abs(xc) <= xb).sum():
        return None
    x = xc[sel]; y = prof[sel]
    w_in = (x >= -xb) & (x <= -0.5); w_out = (x >= 0.5) & (x <= xb)
    WE = np.nanmean(y[w_in]) / np.nanmean(y[w_out])
    A = np.vstack([np.ones_like(x), x, (x < 0).astype(float)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    # step-location scan
    best = (-np.inf, np.nan)
    for xs in np.arange(-2.0, 2.01, 0.1):
        A2 = np.vstack([np.ones_like(x), x, (x < xs).astype(float)]).T
        c2, *_ = np.linalg.lstsq(A2, y, rcond=None)
        if c2[2] > best[0]:
            best = (c2[2], xs)
    return {'WE': float(WE), 'c': float(coef[2]), 'cfrac': float(coef[2] / coef[0]), 'grad': float(coef[1]),
            'xbest': float(best[1]), 'cbest': float(best[0]), 'level': float(np.nanmean(y))}


def real_overlap(pl, ra_c, dec_c, real_tree_v):
    cc_v = radec_to_vec(*densify(ra_c, dec_c, 0.05)[:2])
    t = cKDTree(real_tree_v)
    dd, _ = t.query(cc_v)
    return float(np.mean(2 * np.degrees(np.arcsin(dd / 2)) < 1.5))


def run_placements(pl, img, mask, dx, ra, dec, label, centre, smin=None, smax=None):
    """Real + rotation + translation + radial-scan placements; returns dict of stats arrays."""
    rv = radec_to_vec(*densify(ra, dec, 0.05)[:2])
    out = {}
    pr = pl.profile(ra, dec, img, mask, dx, smin, smax)
    out['real'] = stats(pr[0], pr[1]); out['real_profile'] = (pr[0].tolist(), pr[1].tolist())
    rc, dc = CENTRES[centre]
    rot = []
    for th in np.arange(-180, 180, 4.0):
        if abs(th) < 8:
            continue
        r2, d2 = rotate_about(ra, dec, rc, dc, th)
        g2 = SkyCoord(r2 * u.deg, d2 * u.deg).galactic
        if np.any(g2.l.deg < 180.3) or np.any(g2.l.deg > 226) or np.any(g2.b.deg < -9) or np.any(g2.b.deg > 35):
            continue
        if real_overlap(pl, r2, d2, rv) > 0.1:
            continue
        p = pl.profile(r2, d2, img, mask, dx, smin, smax)
        if p is None:
            continue
        valid = np.nansum(p[2][np.abs(p[0]) <= XB]) / max(np.nansum(p[3][np.abs(p[0]) <= XB]), 1)
        st = stats(p[0], p[1])
        if st is None or valid < 0.5:
            continue
        st['theta'] = float(th); st['bmin'] = float(g2.b.deg.min())
        rot.append(st)
    out['rot'] = rot
    tr = []
    ra0, dec0 = np.mean(ra), np.mean(dec)
    for k in range(400):
        dra, ddec = rng.uniform(-22, 22), rng.uniform(-18, 18)
        if np.hypot(dra, ddec) < 3:
            continue
        r2, d2 = translate(ra, dec, dra, ddec, ra0, dec0)
        g2 = SkyCoord(r2 * u.deg, d2 * u.deg).galactic
        if np.any(g2.l.deg < 180.3) or np.any(g2.l.deg > 226) or np.any(g2.b.deg < -9) or np.any(g2.b.deg > 35):
            continue
        if real_overlap(pl, r2, d2, rv) > 0.1:
            continue
        p = pl.profile(r2, d2, img, mask, dx, smin, smax)
        if p is None:
            continue
        valid = np.nansum(p[2][np.abs(p[0]) <= XB]) / max(np.nansum(p[3][np.abs(p[0]) <= XB]), 1)
        st = stats(p[0], p[1])
        if st is None or valid < 0.5:
            continue
        st['dra'] = float(dra); st['ddec'] = float(ddec); st['bmin'] = float(g2.b.deg.min())
        tr.append(st)
        if len(tr) >= 150:
            break
    out['tr'] = tr
    scan = []
    sep, pa = sep_pa(rc, dc, ra, dec)
    for dr in np.arange(-6, 6.01, 0.5):
        r2, d2 = offset_point(rc, dc, sep + dr, pa)
        p = pl.profile(r2, d2, img, mask, dx, smin, smax)
        st = stats(p[0], p[1]) if p is not None else None
        if st is not None:
            st['dr'] = float(dr); scan.append(st)
    out['scan'] = scan
    return out


def summarize(out, key='c'):
    obs = out['real'][key] if out['real'] else np.nan
    res = {'obs': obs}
    for nm in ['rot', 'tr']:
        v = np.array([s[key] for s in out[nm]])
        if v.size < 5:
            res[nm] = None; continue
        res[nm] = {'n': int(v.size), 'mean': float(v.mean()), 'std': float(v.std(ddof=1)),
                   'z': float((obs - v.mean()) / v.std(ddof=1)),
                   'p_ge': float((np.sum(v >= obs) + 1) / (v.size + 1))}
    return res

# ------------------------------------------------------------------ smooth Monogem shell model
def shell_model(d, img, mask, centre, rmax=20.0, dr=0.25, pa_bin=15.0, exclude_silk=True, nmin=20):
    """Smooth model I = B(PA) + A(PA) g(r): g(r) = mean normalised radial profile over all
    sectors (excluding the Silk PA range if requested); B, A per 15-deg sector by least squares,
    then smoothed in PA with a 3-sector running mean. Returns model array (nan beyond rmax)."""
    rc, dc = CENTRES[centre]
    r, pa = sep_pa(rc, dc, d['ra'].astype(float), d['dec'].astype(float))
    redges = np.arange(0, rmax + dr, dr); nr = len(redges) - 1
    pedges = np.arange(0, 360 + pa_bin, pa_bin); npa = len(pedges) - 1
    ri = np.digitize(r, redges) - 1; pi = np.digitize(pa, pedges) - 1
    ok = mask & np.isfinite(img) & (ri >= 0) & (ri < nr) & (d['gl'] >= 180.3)
    S = np.zeros((npa, nr)); N = np.zeros((npa, nr))
    np.add.at(S, (pi[ok], ri[ok]), img[ok]); np.add.at(N, (pi[ok], ri[ok]), 1)
    P = np.where(N > nmin, S / np.maximum(N, 1), np.nan)
    sp, spa = sep_pa(rc, dc, *load_spine()[:2])
    silk_pa = (pedges[:-1] + pa_bin / 2 > spa.min() - pa_bin) & (pedges[:-1] + pa_bin / 2 < spa.max() + pa_bin)
    # normalise each sector profile between its 2-6 deg mean and 16-20 deg mean
    rc_ = 0.5 * (redges[1:] + redges[:-1])
    g_list = []
    for k in range(npa):
        if exclude_silk and silk_pa[k]:
            continue
        p = P[k]
        inn = np.nanmean(p[(rc_ > 2) & (rc_ < 6)]); outv = np.nanmean(p[(rc_ > 16) & (rc_ < 20)])
        if not (np.isfinite(inn) and np.isfinite(outv)) or inn <= outv or np.isnan(p).mean() > 0.3:
            continue
        g_list.append((p - outv) / (inn - outv))
    g = np.nanmedian(np.array(g_list), axis=0)
    Bk = np.full(npa, np.nan); Ak = np.full(npa, np.nan)
    for k in range(npa):
        y = P[k]; s = np.isfinite(y) & np.isfinite(g)
        if s.sum() < 0.6 * nr:
            continue
        A = np.vstack([np.ones(s.sum()), g[s]]).T
        c, *_ = np.linalg.lstsq(A, y[s], rcond=None)
        Bk[k], Ak[k] = c
    def circ_smooth(v):
        vv = np.concatenate([v[-1:], v, v[:1]])
        out = np.array([np.nanmean(vv[i:i + 3]) for i in range(len(v))])
        return out
    Bs, As = circ_smooth(Bk), circ_smooth(Ak)
    # interpolate B, A in PA (linear, periodic) and g in r
    pc = pedges[:-1] + pa_bin / 2
    def interp_per(v):
        good = np.isfinite(v)
        return np.interp(pa, np.concatenate([pc[good] - 360, pc[good], pc[good] + 360]), np.tile(v[good], 3))
    model = interp_per(Bs) + interp_per(As) * np.interp(r, rc_, np.nan_to_num(g, nan=0.0))
    model[r > rmax] = np.nan
    return model, {'g': g.tolist(), 'r': rc_.tolist(), 'n_sectors_used': len(g_list),
                   'B': Bk.tolist(), 'A': Ak.tolist(), 'pa_c': pc.tolist(), 'P': np.where(np.isfinite(P), P, -1).tolist()}


def main():
    ra, dec, tr = load_spine()
    dsp = dec
    results = {'date': '2026-11-16', 'notes': 'x<0 = west/inside (toward Monogem centres); c>0 = brighter inside'}
    E = load_erass(); plE = Placer(E)
    cc = CurveCoords(ra, dec); s_all = cc.s
    s_kink = float(np.interp(17.6, cc.dec, cc.s))  # along-spine coordinate of Dec 17.6
    segs = {'all': (None, None), 'south': (0.0, s_kink), 'north': (s_kink, float(s_all[-1]))}
    results['s_kink_deg'] = s_kink
    nulls = {}
    sigmas = {'raw': 0.0, 'nh5e-21': 5e-21, 'nh1e-20': 1e-20}
    # main runs: eRASS soft, mask variants, centres, segments
    runs = []
    for band in ['soft', 'mid', 'hard']:
        for mvar in (['none', 'cat', 'full'] if band == 'soft' else ['full']):
            for sig_k in (sigmas if (band == 'soft' and mvar == 'full') else ['raw']):
                img = E[band] * np.exp(sigmas[sig_k] * np.nan_to_num(E['nhi'].astype(float)))
                for centre in (['PSR', 'K2'] if (band == 'soft' and mvar == 'full' and sig_k == 'raw') else ['PSR']):
                    for seg, (a, b) in segs.items():
                        if band != 'soft' and seg != 'all':
                            continue
                        if mvar != 'full' and seg != 'all':
                            continue
                        runs.append(('eRASS', band, mvar, sig_k, centre, seg, img, E[f'mask_{mvar}'], a, b))
    # smooth-shell residual (eRASS soft full mask, PSR and K2 centres)
    for centre in ['PSR', 'K2']:
        model, info = shell_model(E, E['soft'], E['mask_full'], centre)
        results[f'shell_model_{centre}'] = {k: v for k, v in info.items() if k in ('g', 'r', 'n_sectors_used', 'pa_c', 'A', 'B')}
        resid = E['soft'] - model
        np.save(DER / f'verify_xrayfuv_shellresid_{centre}.npy', resid.astype(np.float32))
        for seg, (a, b) in segs.items():
            runs.append(('eRASS', 'soft_shellresid', 'full', 'raw', centre, seg, resid, E['mask_full'] & np.isfinite(resid), a, b))
    # exposure and N(HI) as 'images' (systematics)
    runs.append(('eRASS', 'exposure', 'none', 'raw', 'PSR', 'all', E['expo'].astype(float), np.isfinite(E['expo']), None, None))
    runs.append(('eRASS', 'nhi', 'none', 'raw', 'PSR', 'all', E['nhi'].astype(float), np.isfinite(E['nhi']), None, None))
    # ROSAT
    R = load_rosat(); plR = Placer(R)
    for band in ['r12', 'r45']:
        for seg, (a, b) in segs.items():
            runs.append(('ROSAT', band, 'none', 'raw', 'PSR', seg, R[band], R['mask_none'], a, b))
        img = R['r12'] * np.exp(1e-20 * R['nhi']) if band == 'r12' else R['r45'] * np.exp(5e-21 * R['nhi'])
        runs.append(('ROSAT', band, 'none', 'nh_corr', 'PSR', 'all', img, R['mask_none'], None, None))
    model, info = shell_model(R, R['r12'], R['mask_none'], 'PSR', dr=0.4, nmin=2)
    results['shell_model_rosat_r12_PSR'] = {k: v for k, v in info.items() if k in ('g', 'r', 'n_sectors_used')}
    resid = R['r12'] - model
    for seg, (a, b) in segs.items():
        runs.append(('ROSAT', 'r12_shellresid', 'none', 'raw', 'PSR', seg, resid, R['mask_none'] & np.isfinite(resid), a, b))

    table = []
    for (survey, band, mvar, sig_k, centre, seg, img, mask, a, b) in runs:
        pl, dd = (plE, E) if survey == 'eRASS' else (plR, R)
        out = run_placements(pl, img, mask, dd['dx'], ra, dec, band, centre, a, b)
        key = f'{survey}|{band}|{mvar}|{sig_k}|{centre}|{seg}'
        row = {'key': key, 'real': out['real'],
               'c': summarize(out, 'c'), 'cfrac': summarize(out, 'cfrac'), 'WE': summarize(out, 'WE'),
               'scan': [{k: s[k] for k in ('dr', 'c', 'cfrac', 'WE', 'xbest')} for s in out['scan']],
               'real_profile': out['real_profile']}
        table.append(row)
        nulls[key] = out
        r = row
        def fz(x):
            return 'NA' if x is None else f"{x['z']:+.2f}/p={x['p_ge']:.3f}(n={x['n']})"
        print(f"{key:55s} WE={r['real']['WE']:.3f} c={r['real']['c']:+.3g} cfrac={r['real']['cfrac']:+.3f} xbest={r['real']['xbest']:+.1f} | "
              f"c: rot {fz(r['c']['rot'])} tr {fz(r['c']['tr'])} | WE: rot {fz(r['WE']['rot'])} tr {fz(r['WE']['tr'])}", flush=True)
    results['runs'] = table
    (DER / 'verify_xrayfuv_xray_edge.json').write_text(json.dumps(results, indent=1, default=float))
    # compact null arrays
    arr = {}
    for key, out in nulls.items():
        for nm in ['rot', 'tr']:
            for st in ['c', 'cfrac', 'WE', 'xbest']:
                arr[f'{key}|{nm}|{st}'] = np.array([s[st] for s in out[nm]])
        arr[f'{key}|rot|theta'] = np.array([s['theta'] for s in out['rot']])
    np.savez_compressed(DER / 'verify_xrayfuv_xray_nulls.npz', **arr)


if __name__ == '__main__':
    main()
