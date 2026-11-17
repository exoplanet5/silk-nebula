"""verify_xrayfuv: re-measure the GALEX FUV/NUV excess on the Silk (independent of wave-1 code).

Input: data/derived/verify_xrayfuv_galex_blocks.npz (1' blocks of the 267 skybg tiles, CU).
Steps per band: (1) per-tile median offset; (2) mean fractional radial (vignetting/flat) pattern
across all tiles of the band, removed; (3) 4-MAD clipping within tiles; (4) two estimators:
  (a) template regression on the NSNS DR0.2 [O III] map (Rayleigh; high-passed with a 1-deg
      median filter, smoothed to 3' FWHM): v = a_tile + beta * T  -> beta in CU per NSNS-Rayleigh;
  (b) trace contrast: mean(|x| <= h) - mean(0.30 <= |x| <= 0.55) within each tile (E spine and
      the W-strand), weighted over tiles, h = 0.15 deg.
Nulls: template / trace shifted in RA by +-0.7..+-3.0 deg (0.1 step) and, for (a), also in
Dec by +-1..3 deg combined with RA shifts. Output: data/derived/verify_xrayfuv_galex_measure.json
"""
import sys, json
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import *

z = np.load(DER / 'verify_xrayfuv_galex_blocks.npz')
G = {k: z[k] for k in z.files}
ns = np.load(DER / 'imaging_optical_nsns_maps.npz')
hdr = fits.Header.fromstring(str(ns['wcs_header']))
wn = WCS(hdr)
T = np.nan_to_num(ns['oiii'].astype(float), nan=0.0)
Tc = np.clip(T, -1, 20)
Tm = ndimage.median_filter(Tc, size=61, mode='nearest')            # ~1 deg
Thp = ndimage.gaussian_filter(Tc - Tm, sigma=3.0 / 2.3548)          # 3' FWHM
# star masks: NSNS variance (bright-star residuals) and Tycho-2 VT<11 (r = max(2', 12' x 10^(-0.2 (VT-5))), cap 15')
SNO = ns['sn_oiii'].astype(float); SNH = ns['sn_halpha'].astype(float)
thrO = np.nanpercentile(SNO, 98); thrH = np.nanpercentile(SNH, 98)
badN = ndimage.binary_dilation((SNO > thrO) | (SNH > thrH), iterations=2)
_ty = read_vizier_tsv(ROOT / 'data/raw/vizier_tycho2/tyc2_VT11_box.tsv')
_f = np.isfinite(_ty['RAmdeg']) & np.isfinite(_ty['DEmdeg']) & np.isfinite(_ty['VTmag'])
TY_RA, TY_DE, TY_VT = _ty['RAmdeg'][_f], _ty['DEmdeg'][_f], _ty['VTmag'][_f]
VC = ndimage.gaussian_filter(np.nan_to_num(ns['vc'].astype(float)) - ndimage.median_filter(np.nan_to_num(ns['vc'].astype(float)), 61, mode='nearest'), 3.0 / 2.3548)
Ha = np.nan_to_num(ns['halpha'].astype(float)); Ha_hp = ndimage.gaussian_filter(np.clip(Ha, -2, 30) - ndimage.median_filter(np.clip(Ha, -2, 30), 61, mode='nearest'), 3.0 / 2.3548)

ra_sp, dec_sp, trace = load_spine()


def strand_coords(name):
    s = trace['strands'][name]
    pts = s['points'] if isinstance(s, dict) and 'points' in s else s
    if isinstance(pts, dict):
        return np.array(pts['ra']), np.array(pts['dec'])
    ra = np.array([p['ra'] if isinstance(p, dict) else p[0] for p in pts])
    de = np.array([p['dec'] if isinstance(p, dict) else p[1] for p in pts])
    return ra, de


def sample_template(img, ra, dec, dra=0.0, ddec=0.0):
    x, y = wn.world_to_pixel_values(ra - dra / np.cos(np.radians(dec)), dec - ddec)
    return ndimage.map_coordinates(img, [y, x], order=1, mode='constant', cval=np.nan)


def prep_band(band):
    s = G['band'] == band
    tile = G['tile'][s]; cu = G['cu'][s].astype(float); rt = G['rtile'][s]
    ra = G['ra'][s]; dec = G['dec'][s]
    med = np.zeros_like(cu)
    for t in np.unique(tile):
        m = tile == t
        med[m] = np.median(cu[m])
    frac = cu / med - 1
    edges = np.arange(0, 0.5501, 0.025)
    ri = np.clip(np.digitize(rt, edges) - 1, 0, len(edges) - 2)
    pat = np.array([np.median(frac[ri == k]) for k in range(len(edges) - 1)])
    resid = cu - med * (1 + pat[ri])
    keep = np.ones_like(resid, bool)
    # star masks
    from scipy.spatial import cKDTree
    rad = np.clip(12.0 * 10 ** (-0.2 * (TY_VT - 5.0)), 2.0, 15.0) / 60.0
    tt = cKDTree(radec_to_vec(TY_RA, TY_DE))
    bv = radec_to_vec(ra, dec)
    near = tt.query_ball_point(bv, r=2 * np.sin(np.radians(15.0 / 60) / 2))
    for i, lst in enumerate(near):
        if lst:
            sep = 2 * np.degrees(np.arcsin(np.linalg.norm(tt.data[lst] - bv[i], axis=1) / 2))
            if np.any(sep < rad[lst]):
                keep[i] = False
    xn, yn = wn.world_to_pixel_values(ra, dec)
    bn_ = ndimage.map_coordinates(badN.astype(float), [yn, xn], order=0, mode='constant', cval=0.0) > 0.5
    star_frac = float((~keep).mean()); keep &= ~bn_
    for t in np.unique(tile):
        m = tile == t
        mad = 1.4826 * np.median(np.abs(resid[m] - np.median(resid[m])))
        keep[m] &= np.abs(resid[m] - np.median(resid[m])) < 4 * mad
    return dict(tile=tile[keep], v=resid[keep], ra=ra[keep], dec=dec[keep], med=med[keep],
                pattern=pat.tolist(), clip_frac=float(1 - keep.mean()), star_frac=star_frac, median_level=float(np.median(cu)))


def regress(D, tmpl):
    ok = np.isfinite(tmpl)
    tile = D['tile'][ok]; v = D['v'][ok]; t = tmpl[ok]
    # demean per tile
    ut, inv = np.unique(tile, return_inverse=True)
    vm = np.bincount(inv, v) / np.bincount(inv); tm = np.bincount(inv, t) / np.bincount(inv)
    dv = v - vm[inv]; dt = t - tm[inv]
    beta = np.sum(dv * dt) / np.sum(dt * dt)
    return float(beta), float(np.sum(dt * dt)), int(ok.sum())


def contrast(D, ra_c, dec_c, h=0.15, off=(0.30, 0.55), dra=0.0):
    cc = CurveCoords(ra_c + dra / np.cos(np.radians(dec_c)), dec_c)
    x, s = cc.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.0)
    on = np.abs(x) <= h; of = (np.abs(x) >= off[0]) & (np.abs(x) <= off[1])
    vals = []; wts = []
    for t in np.unique(D['tile'][on]):
        m = D['tile'] == t
        a = D['v'][m & on]; b = D['v'][m & of]
        if a.size >= 10 and b.size >= 10:
            vals.append(a.mean() - b.mean()); wts.append(a.size)
    if not vals:
        return np.nan, 0
    vals = np.array(vals); wts = np.array(wts, float)
    return float(np.sum(vals * wts) / wts.sum()), len(vals)



def main():
    res = {'date': '2026-11-16'}
    shifts = np.round(np.concatenate([np.arange(-3.0, -0.69, 0.1), np.arange(0.7, 3.01, 0.1)]), 2)
    for band in ['FUV', 'NUV']:
        D = prep_band(band)
        r = {'n_blocks': int(D['v'].size), 'n_tiles': int(np.unique(D['tile']).size), 'clip_frac': D['clip_frac'], 'star_mask_frac': D['star_frac'],
             'median_level_CU': D['median_level'], 'radial_pattern_frac': D['pattern']}
        # (a) template regression
        for tname, img in [('OIII', Thp), ('Halpha', Ha_hp)]:
            for region, sel in [('all', np.ones_like(D['v'], bool)), ('south_dec<17.6', D['dec'] < 17.6),
                                ('north_dec>=17.6', D['dec'] >= 17.6)]:
                Dd = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == D['v'].shape else v) for k, v in D.items()}
                beta, sxx, n = regress(Dd, sample_template(img, Dd['ra'], Dd['dec']))
                nb = []
                for dra in shifts:
                    nb.append(regress(Dd, sample_template(img, Dd['ra'], Dd['dec'], dra=dra))[0])
                for dra in [-2.0, -1.2, 1.2, 2.0]:
                    for ddec in [-2.5, -1.5, 1.5, 2.5]:
                        nb.append(regress(Dd, sample_template(img, Dd['ra'], Dd['dec'], dra=dra, ddec=ddec))[0])
                nb = np.array(nb); nb = nb[np.isfinite(nb)]
                r[f'beta_{tname}_{region}'] = {'beta_CU_per_R': beta, 'null_mean': float(nb.mean()), 'null_std': float(nb.std(ddof=1)),
                                               'n_null': int(nb.size), 'z': float((beta - nb.mean()) / nb.std(ddof=1)), 'n_blocks': n,
                                               'p_ge': float((np.sum(nb >= beta) + 1) / (nb.size + 1))}
                print(band, tname, region, r[f'beta_{tname}_{region}'], flush=True)
        # (a2) joint regression on [O III] + Halpha + vc templates
        for region, sel in [('all', np.ones_like(D['v'], bool)), ('south_dec<17.6', D['dec'] < 17.6), ('north_dec>=17.6', D['dec'] >= 17.6)]:
            Dd = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == D['v'].shape else v) for k, v in D.items()}
            def joint(dra=0.0, ddec=0.0):
                X = np.vstack([sample_template(img, Dd['ra'], Dd['dec'], dra, ddec) for img in (Thp, Ha_hp, VC)]).T
                ok = np.all(np.isfinite(X), axis=1)
                ut, inv = np.unique(Dd['tile'][ok], return_inverse=True)
                v = Dd['v'][ok]; X = X[ok]
                cnt = np.bincount(inv)
                dv = v - (np.bincount(inv, v) / cnt)[inv]
                dX = X - np.vstack([np.bincount(inv, X[:, j]) / cnt for j in range(3)]).T[inv]
                coef, *_ = np.linalg.lstsq(dX, dv, rcond=None)
                return coef
            c0 = joint(); nb = np.array([joint(dra=d) for d in shifts] + [joint(dra=a, ddec=b) for a in (-2.0, -1.2, 1.2, 2.0) for b in (-2.5, -1.5, 1.5, 2.5)])
            r[f'joint_{region}'] = {'beta_OIII': float(c0[0]), 'beta_Ha': float(c0[1]), 'beta_vc': float(c0[2]),
                                    'null_std_OIII': float(nb[:, 0].std(ddof=1)), 'null_mean_OIII': float(nb[:, 0].mean()),
                                    'null_std_Ha': float(nb[:, 1].std(ddof=1)), 'null_std_vc': float(nb[:, 2].std(ddof=1)),
                                    'z_OIII': float((c0[0] - nb[:, 0].mean()) / nb[:, 0].std(ddof=1)),
                                    'z_Ha': float((c0[1] - nb[:, 1].mean()) / nb[:, 1].std(ddof=1))}
            print(band, 'joint', region, r[f'joint_{region}'], flush=True)
        # (b) trace contrast (E spine, W strand), with RA-shift nulls
        for nm, (rc, dc) in [('E', (ra_sp, dec_sp)), ('W', strand_coords('W'))]:
            for region, sel in [('all', np.ones_like(D['v'], bool)), ('south_dec<17.6', D['dec'] < 17.6)]:
                Dd = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == D['v'].shape else v) for k, v in D.items()}
                c0, nt = contrast(Dd, rc, dc)
                nb = np.array([contrast(Dd, rc, dc, dra=d)[0] for d in shifts]); nb = nb[np.isfinite(nb)]
                r[f'contrast_{nm}_{region}'] = {'C_CU': c0, 'n_tiles': nt, 'null_mean': float(nb.mean()),
                                                'null_std': float(nb.std(ddof=1)), 'n_null': int(nb.size),
                                                'z': float((c0 - nb.mean()) / nb.std(ddof=1))}
                print(band, 'contrast', nm, region, r[f'contrast_{nm}_{region}'], flush=True)
        res[band] = r
    (DER / 'verify_xrayfuv_galex_measure.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
