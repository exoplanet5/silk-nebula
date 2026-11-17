"""F6 data: cross-spine profiles of the southern Silk (spine Dec < 17.6) in several bands, with null placements.

Offset x = signed perpendicular distance from the E spine (silk_trace.json v1.0), + = east (outside, away from
Monogem), computed with verify_xrayfuv_common.CurveCoords. "South" = pixels whose nearest spine point has
Dec < 17.6 (along-spine s <= s_split). The same estimator is applied to the real trace and to null placements.

Bands / inputs / nulls:
  NSNS [O III]   w3_figs_f1_maps.npz (NSNS DR0.2, 2'x2' block means)       RA shifts +-1.5..+-5 deg x Dec shifts
                 -2/0/+2 deg (48)
  DESI [O III]   verify_desi_linefluxes.npz x kR_5008 (GFA acceptance), verify_desi good_fit + ok_contrib  (no nulls;
                 per-bin standard errors; fibre residuals are relative to the petal sky mean)
  GALEX FUV/NUV  verify_xrayfuv_galex_measure.prep_band (tile offsets, flat pattern, star masks, clipping); per-tile
                 demeaning with |x| > 0.3 blocks as in verify_xrayfuv_fig_uv          RA shifts +-0.8..+-3.0 (24)
  SHASSA Ha      verify_hiwham_hamaps_part_Sm.npz (Dec 12-15.6; profile + 400 placement nulls, verify_hiwham)
  eRASS1 0.2-0.6 verify_xrayfuv_erass_pixels.npz (catmask applied)          rotations about PSR B0656+14 (same
                 Monogem radius), |theta| = 20..180 deg in 10 deg steps
  GALFA HI       imaging_mw_eq_galfa_nhi_m02p06.fits (N_HI, V_LSR -2..+6)   as NSNS (needs >= 80 % valid bins; the
                 GALFA tiles east of RA 119.33 are missing, so only western placements survive)
  IRIS 100 um    imaging_mw_eq_iris100.fits                                 as NSNS
Null placements ignore pixels within 1.0 deg (GALEX: 0.6 deg) of the real E spine, so the Silk itself does not
enter the null profiles. Ridge amplitude: least squares of a Gaussian (FWHM 0.28 deg, centre 0) + linear baseline over |x| <= 1.2 (GALEX 1.0);
z = (A - median A_null) / (1.4826 MAD of A_null).
Output: data/derived/w3_figs_f6_profiles.npz, w3_figs_f6_profiles.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_f6_profiles_data.py
"""
import sys
import json
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import CurveCoords, radec_to_vec, load_spine, rotate_about, CENTRES
from w3_figs_common import DER, json_dump

FWHM = 0.28
ra_sp, dec_sp, _ = load_spine()
CC0 = CurveCoords(ra_sp, dec_sp)
S_SPLIT = float(np.interp(17.6, CC0.dec, CC0.s))


def shifted(dra, ddec=0.0):
    return CurveCoords(ra_sp + dra / np.cos(np.radians(dec_sp)), dec_sp + ddec)


def rotated(theta):
    r, d = rotate_about(ra_sp, dec_sp, *CENTRES['PSR'], theta)
    return CurveCoords(r, d)


def near_real(ra, dec, dist):
    """True for points within `dist` deg of the real E spine (any s), used to keep null placements off the Silk."""
    x, s = CC0.coords(radec_to_vec(ra, dec), maxdist_deg=dist + 0.5)
    return np.isfinite(x) & (np.abs(x) < dist)


def profile(cc, ra, dec, val, edges, stat='median', smax=S_SPLIT, maxd=3.0, excl=None):
    x, s = cc.coords(radec_to_vec(ra, dec), maxdist_deg=maxd)
    sel = np.isfinite(x) & (s <= smax) & np.isfinite(val)
    if excl is not None:
        sel &= ~excl
    bi = np.digitize(x[sel], edges) - 1
    ok = (bi >= 0) & (bi < len(edges) - 1)
    v = val[sel][ok]; bi = bi[ok]
    out = np.full(len(edges) - 1, np.nan); n = np.zeros(len(edges) - 1, int); e = np.full(len(edges) - 1, np.nan)
    for k in range(len(edges) - 1):
        vv = v[bi == k]
        n[k] = vv.size
        if vv.size >= 5:
            out[k] = np.median(vv) if stat == 'median' else np.mean(vv)
            e[k] = np.std(vv) / np.sqrt(vv.size) * (1.25 if stat == 'median' else 1.0)
    return out, n, e


def ridge_amp(xc, p, xmax=1.2):
    m = np.isfinite(p) & (np.abs(xc) <= xmax)
    if m.sum() < 10:
        return np.nan
    sg = FWHM / 2.3548
    X = np.vstack([np.exp(-0.5 * (xc[m] / sg) ** 2), np.ones(m.sum()), xc[m]]).T
    return float(np.linalg.lstsq(X, p[m], rcond=None)[0][0])


def summarize(name, xc, real, nulls, unit, xmax=1.2, extra=None):
    nulls = np.array([q for q in nulls if np.isfinite(q).mean() >= 0.8]) if len(nulls) else np.zeros((0, xc.size))
    A = ridge_amp(xc, real, xmax)
    An = np.array([ridge_amp(xc, q, xmax) for q in nulls]) if len(nulls) else np.array([])
    An = An[np.isfinite(An)]
    res = dict(name=name, unit=unit, A=A, n_null=int(len(An)))
    if len(An) >= 5:
        med = float(np.median(An)); sig = float(1.4826 * np.median(np.abs(An - med)))
        res.update(A_null_median=med, A_null_sig=sig, z=(A - med) / sig if sig > 0 else np.nan,
                   A_null_std=float(np.std(An)))
    if extra:
        res.update(extra)
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()}, flush=True)
    return res, nulls


def main():
    out, summ = {}, {}
    shifts16 = np.r_[-5:-1.4:0.5, 1.5:5.01:0.5]
    shifts48 = [(a, b) for a in shifts16 for b in (-2.0, 0.0, 2.0)]
    # ---------------- NSNS [O III]
    M = np.load(DER / 'w3_figs_f1_maps.npz')
    o3 = M['oiii'].astype(float)
    ny, nx = o3.shape
    o3b = o3[:ny // 2 * 2, :nx // 2 * 2].reshape(ny // 2, 2, nx // 2, 2).mean(axis=(1, 3))
    rab = M['ra'][:nx // 2 * 2].reshape(-1, 2).mean(1); deb = M['dec'][:ny // 2 * 2].reshape(-1, 2).mean(1)
    RA, DE = np.meshgrid(rab, deb)
    RA, DE, V = RA.ravel(), DE.ravel(), o3b.ravel()
    e05 = np.arange(-2.5, 2.5001, 0.05); xc05 = 0.5 * (e05[1:] + e05[:-1])
    p, n, e = profile(CC0, RA, DE, V, e05)
    ex = near_real(RA, DE, 1.0)
    nl = [profile(shifted(*d), RA, DE, V, e05, excl=ex)[0] for d in shifts48]
    summ['nsns'], nl = summarize('NSNS [O III]', xc05, p, nl, 'R')
    out.update(nsns_x=xc05, nsns=p, nsns_null=nl)
    # ---------------- DESI [O III] fibre residuals
    from verify_desi_common import good_fit
    lf = np.load(DER / 'verify_desi_linefluxes.npz'); cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    det = json.loads((DER / 'verify_desi_detect.json').read_text())
    okd = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & cm['ok_contrib'] & np.isfinite(cm['kR_5008'])
    vR = lf['o3__amp'][:, 0] * cm['kR_5008']
    e10 = np.arange(-2.5, 2.5001, 0.1); xc10 = 0.5 * (e10[1:] + e10[:-1])
    pd_, nd_, ed_ = profile(CC0, cm['ra'][okd], cm['dec'][okd], vR[okd], e10, stat='mean')
    base = np.nanmedian(pd_[(np.abs(xc10) > 0.6) & (np.abs(xc10) < 2.0)])
    out.update(desi_x=xc10, desi=pd_ - base, desi_err=ed_, desi_n=nd_)
    summ['desi'] = dict(name='DESI [O III]', unit='R', baseline_subtracted=float(base), n_fibres=int(nd_.sum()),
                        A=ridge_amp(xc10, pd_ - base))
    print('DESI', summ['desi'], flush=True)
    # ---------------- GALEX FUV / NUV
    import verify_xrayfuv_galex_measure as GM
    eg = np.arange(-1.2, 1.2001, 0.05); xcg = 0.5 * (eg[1:] + eg[:-1])
    gshifts = np.r_[-3.0:-0.79:0.2, 0.8:3.01:0.2]

    def gal_prof(D, cc, excl=None):
        x, s = cc.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.5)
        sel = np.isfinite(x) & (s <= S_SPLIT)
        if excl is not None:
            sel &= ~excl
        v = D['v'].copy()
        for t in np.unique(D['tile'][sel]):
            m = (D['tile'] == t) & sel
            ref = m & (np.abs(x) > 0.3)
            if ref.sum() > 20:
                v[m] -= np.mean(v[ref])
            else:
                sel &= ~m
        bi = np.digitize(x[sel], eg) - 1
        ok = (bi >= 0) & (bi < len(xcg))
        cnt = np.bincount(bi[ok], minlength=len(xcg))
        pr = np.bincount(bi[ok], v[sel][ok], len(xcg)) / np.maximum(cnt, 1)
        pr[cnt < 20] = np.nan
        return pr
    for band in ('FUV', 'NUV'):
        D = GM.prep_band(band)
        pg = gal_prof(D, CC0)
        exg = near_real(D['ra'], D['dec'], 0.6)
        nl = [gal_prof(D, shifted(d), excl=exg) for d in gshifts]
        summ[band.lower()], nl = summarize('GALEX ' + band, xcg, pg, nl, 'CU', xmax=1.0)
        out.update(**{f'{band.lower()}_x': xcg, band.lower(): pg, f'{band.lower()}_null': nl})
    # ---------------- SHASSA (verify_hiwham product)
    H = np.load(DER / 'verify_hiwham_hamaps_part_Sm.npz')
    hs, hn = H['SHASSA_Sm'], H['SHASSA_Sm_nulls']
    summ['shassa'], hn = summarize('SHASSA Ha', H['cen'], hs, list(hn), 'R',
                                   extra=dict(source='verify_hiwham_hamaps_part_Sm.npz', dec_range=[12.0, 15.6]))
    out.update(shassa_x=H['cen'], shassa=hs, shassa_null=hn)
    # ---------------- eRASS1 0.2-0.6 keV
    X = np.load(DER / 'verify_xrayfuv_erass_pixels.npz')
    okx = X['good'] & ~X['catmask']
    xv = (X['b0204'] + X['b0406']).astype(float)[okx]
    xra, xde = X['ra'][okx].astype(float), X['dec'][okx].astype(float)
    pre = (xra > 95) & (xra < 140) & (xde > -5) & (xde < 35)
    xra, xde, xv = xra[pre], xde[pre], xv[pre]
    px, _, _ = profile(CC0, xra, xde, xv, e10, stat='mean')
    thetas = np.r_[-180:-19:10, 20:181:10]
    nl = [profile(rotated(t), xra, xde, xv, e10, stat='mean')[0] for t in thetas]
    summ['erass'], nl = summarize('eRASS1 0.2-0.6 keV', xc10, px, nl, 'cts s-1 deg-2',
                                  extra=dict(mean_level=float(np.nanmean(px)), null='rotations about PSR B0656+14'))
    out.update(erass_x=xc10, erass=px, erass_null=nl)

    # ---------------- EQ grids: GALFA HI, IRIS 100
    def eq_grid(name):
        with fits.open(DER / f'imaging_mw_grids/imaging_mw_eq_{name}.fits') as h:
            d = h[0].data.astype(float); hd = h[0].header
        ii = np.arange(d.shape[1]); jj = np.arange(d.shape[0])
        ra = hd['CRVAL1'] + (ii + 1 - hd['CRPIX1']) * hd['CDELT1']
        de = hd['CRVAL2'] + (jj + 1 - hd['CRPIX2']) * hd['CDELT2']
        R_, D_ = np.meshgrid(ra, de)
        return R_.ravel(), D_.ravel(), d.ravel()
    for key, name, unit, scale in (('hi', 'galfa_nhi_m02p06', '1e19 cm-2', 1e-19), ('iris', 'iris100', 'MJy sr-1', 1.0)):
        R_, D_, V_ = eq_grid(name)
        V_ = V_ * scale
        pq, _, _ = profile(CC0, R_, D_, V_, e10)
        ex = near_real(R_, D_, 1.0)
        nl = [profile(shifted(*d), R_, D_, V_, e10, excl=ex)[0] for d in shifts48]
        summ[key], nl = summarize(name, xc10, pq, nl, unit)
        out.update(**{f'{key}_x': xc10, key: pq, f'{key}_null': nl})
    out['s_split'] = S_SPLIT
    np.savez_compressed(DER / 'w3_figs_f6_profiles.npz', **{k: np.asarray(v) for k, v in out.items()})
    json_dump(dict(s_split_deg=S_SPLIT, fwhm_template=FWHM, bands=summ), DER / 'w3_figs_f6_profiles.json')
    print('saved')


if __name__ == '__main__':
    main()
