"""GALEX FUV/NUV diffuse background along the Silk filament from pipeline coadd skybg maps.

Method
  * each -fd/-nd-skybg map (source-free sigma-clipped background estimate, counts/s per 1.5" pixel) is
    block-averaged to 1' (40x40 px) blocks; only blocks within 0.50 deg of the tile centre with >=80%
    valid pixels are kept (tile edges have known flat-field / vignetting artefacts);
  * counts/s/pix -> continuum units CU (photons cm-2 s-1 sr-1 A-1) with the GALEX flux calibration
    (FUV 1.40e-15, NUV 2.06e-16 erg s-1 cm-2 A-1 per cps; lambda_eff 1528 / 2271 A; 1.5" pixels);
  * every block gets (offset from the E-strand spine, along-spine s) from silk_trace.json;
  * a joint robust fit (alternating medians) of value = tile_offset[t] + profile(offset bin) removes
    tile-to-tile airglow / zodiacal offsets; the profile is normalised to zero at 0.8 < |offset| < 2.5 deg.
Outputs: data/derived/imaging_mw_galex_blocks.npz, imaging_mw_galex_profiles.json,
         data/derived/imaging_mw_grids/imaging_mw_galex_{fuv,nuv}_mosaic_1arcmin.fits
"""
import glob
import json
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import RAW, DER, get_spine, trace_coords, spine_s_of, REGIONS

CPS2CU = {'FUV': 1.40e-15 / (6.626e-27 * 2.998e18 / 1528.0) / (1.5 / 206264.806) ** 2,
          'NUV': 2.06e-16 / (6.626e-27 * 2.998e18 / 2271.0) / (1.5 / 206264.806) ** 2}
BLK = 40
RMAX = 0.50

meta = Table.read(RAW / 'galex_mast' / 'galex_skybg.ecsv')
meta = {u.split('/')[-1]: r for u, r in zip(meta['access_url'], meta)}
spine = get_spine()

rows = {k: [] for k in ['band', 'tile', 'ra', 'dec', 'val', 'exptime', 'survey']}
files = sorted(glob.glob(str(RAW / 'galex_mast' / 'skybg' / '*-skybg.fits.gz')))
for it, f in enumerate(files):
    name = f.split('/')[-1]
    band = 'FUV' if '-fd-' in name else 'NUV'
    try:
        h = fits.open(f)
        d = h[0].data.astype(np.float64)
        w = WCS(h[0].header).celestial
        exptime = float(h[0].header.get('EXPTIME', np.nan))
    except Exception as e:
        print('skip', name, e); continue
    n = d.shape[0] // BLK
    d = d[:n * BLK, :n * BLK]
    good = (d > 0) & np.isfinite(d)
    s = np.where(good, d, 0).reshape(n, BLK, n, BLK).sum(axis=(1, 3))
    c = good.reshape(n, BLK, n, BLK).sum(axis=(1, 3))
    frac = c / BLK ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        m = s / c
    yy, xx = np.mgrid[0:n, 0:n]
    ra, dec = w.wcs_pix2world(xx * BLK + (BLK - 1) / 2, yy * BLK + (BLK - 1) / 2, 0)
    r = np.hypot((ra - h[0].header['CRVAL1']) * np.cos(np.radians(dec)), dec - h[0].header['CRVAL2'])
    keep = (frac >= 0.8) & (r <= RMAX) & np.isfinite(m)
    k = keep.sum()
    rows['band'] += [band] * k; rows['tile'] += [it] * k
    rows['ra'] += list(ra[keep]); rows['dec'] += list(dec[keep]); rows['val'] += list(m[keep] * CPS2CU[band])
    rows['exptime'] += [exptime] * k; rows['survey'] += [name.split('_')[0]] * k
    if it % 40 == 0:
        print(it, name, band, k, np.nanmedian(m[keep]) * CPS2CU[band])

B = {k: np.array(v) for k, v in rows.items()}
off, s_al, beyond, L = trace_coords(B['ra'][None, :], B['dec'][None, :], 'eq', spine)
B['off'] = off[0]; B['s'] = s_al[0]; B['beyond'] = beyond[0]
np.savez_compressed(DER / 'imaging_mw_galex_blocks.npz', **B)
print('blocks', len(B['val']), 'spine length', L)

EDGES = np.arange(-3.0, 3.0001, 0.1)
CEN = 0.5 * (EDGES[1:] + EDGES[:-1])
REF = (np.abs(CEN) > 0.8) & (np.abs(CEN) < 2.5)


def joint_fit(sel, n_iter=15):
    t = B['tile'][sel]; v = B['val'][sel]; o = B['off'][sel]
    ib = np.digitize(o, EDGES) - 1
    ok = (ib >= 0) & (ib < len(CEN))
    t, v, ib = t[ok], v[ok], ib[ok]
    tiles = np.unique(t)
    toff = {tt: np.median(v[t == tt]) for tt in tiles}
    prof = np.zeros(len(CEN))
    for _ in range(n_iter):
        resid = v - np.array([toff[tt] for tt in t])
        for i in range(len(CEN)):
            x = resid[ib == i]
            prof[i] = np.median(x) if x.size >= 5 else np.nan
        prof -= np.nanmedian(prof[REF])
        r2 = v - np.where(np.isfinite(prof[ib]), prof[ib], 0)
        for tt in tiles:
            toff[tt] = np.median(r2[t == tt])
    resid = v - np.array([toff[tt] for tt in t])
    err = np.full(len(CEN), np.nan); n = np.zeros(len(CEN), int); ntile = np.zeros(len(CEN), int)
    for i in range(len(CEN)):
        x = resid[ib == i] - (prof[i] if np.isfinite(prof[i]) else 0)
        n[i] = x.size; ntile[i] = len(np.unique(t[ib == i]))
        if x.size >= 5:
            # 1' blocks are nearly independent (skybg smoothing is a few arcmin) -> inflate by sqrt(4) for safety
            err[i] = 1.2533 * 1.4826 * np.median(np.abs(x)) / np.sqrt(x.size) * 2.0
    return prof, err, n, ntile, np.median(v)


def on_off(prof, err):
    on = np.abs(CEN) <= 0.15
    ref = REF & np.isfinite(prof)
    d = np.nanmean(prof[on]) - np.nanmean(prof[ref])
    e = np.sqrt(np.nansum(err[on] ** 2)) / on.sum()
    return float(d), float(e)


segs = {'all': (None, None)}
for k in ['P1', 'P2', 'P3']:
    s0, _ = spine_s_of(*REGIONS[k], spine=spine)
    segs[k] = (s0 - 0.6, s0 + 0.6)
segs['south_dec12-18'] = (0.0, spine_s_of(116.8, 18.0, spine)[0])
segs['north_dec18-28'] = (spine_s_of(116.8, 18.0, spine)[0], 99)

res = {'method': __doc__, 'edges': EDGES.tolist(), 'centres': CEN.tolist(), 'cps_to_CU': CPS2CU, 'profiles': {}}
for band in ['FUV', 'NUV']:
    for seg, (lo, hi) in segs.items():
        sel = (B['band'] == band) & ~B['beyond']
        if lo is not None:
            sel &= (B['s'] >= lo) & (B['s'] <= hi)
        if sel.sum() < 100:
            continue
        prof, err, n, ntile, lev = joint_fit(sel)
        d, e = on_off(prof, err)
        # east-west step: interior (W, concave) minus exterior (E)
        west = (CEN < -0.5) & (CEN > -2.5); east = (CEN > 0.5) & (CEN < 2.5)
        step = float(np.nanmean(prof[west]) - np.nanmean(prof[east]))
        res['profiles'][f'{band}:{seg}'] = {'profile_CU': np.round(prof, 2).tolist(), 'err_CU': np.round(err, 2).tolist(),
                                           'n_blocks': n.tolist(), 'n_tiles': ntile.tolist(), 'median_level_CU': float(lev),
                                           'on_minus_ref_CU': d, 'on_minus_ref_err_CU': e, 'west_minus_east_CU': step,
                                           's_range': [lo, hi]}
        print(f'{band} {seg:16s} level {lev:7.1f} CU  on-ref {d:+7.1f} +- {e:5.1f} CU  W-E {step:+7.1f}  blocks {sel.sum()}')
(DER / 'imaging_mw_galex_profiles.json').write_text(json.dumps(res, indent=1))

# 1' mosaic with tile offsets removed (for display): EQ 1' lattice RA 108-125, Dec 4-30
for band in ['FUV', 'NUV']:
    sel = (B['band'] == band)
    t = B['tile'][sel]; v = B['val'][sel]
    # tile offsets relative to the global median, from blocks far from the filament (or all if none)
    far = np.abs(B['off'][sel]) > 0.8
    glob_med = np.median(v[far])
    corr = np.zeros_like(v)
    for tt in np.unique(t):
        m = t == tt
        mm = m & far
        corr[m] = (np.median(v[mm]) if mm.sum() > 50 else np.median(v[m])) - glob_med
    vv = v - corr
    NX, NY = 1020, 1560
    i = np.floor((125 - B['ra'][sel]) * 60).astype(int); j = np.floor((B['dec'][sel] - 4) * 60).astype(int)
    ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
    acc = np.zeros((NY, NX)); cnt = np.zeros((NY, NX))
    np.add.at(acc, (j[ok], i[ok]), vv[ok]); np.add.at(cnt, (j[ok], i[ok]), 1)
    with np.errstate(invalid='ignore'):
        img = acc / cnt
    h1 = fits.Header()
    h1['CTYPE1'] = 'RA---CAR'; h1['CTYPE2'] = 'DEC--CAR'; h1['CRVAL1'] = 116.5; h1['CRVAL2'] = 0.0
    h1['CDELT1'] = -1 / 60; h1['CDELT2'] = 1 / 60; h1['CRPIX1'] = NX / 2 + 0.5; h1['CRPIX2'] = 1 - (4 + 1 / 120) * 60
    h1['BUNIT'] = 'CU'; h1['COMMENT'] = f'GALEX {band} skybg coadds, 1 arcmin blocks, tile offsets matched'
    fits.PrimaryHDU(img.astype(np.float32), header=h1).writeto(
        DER / 'imaging_mw_grids' / f'imaging_mw_galex_{band.lower()}_mosaic_1arcmin.fits', overwrite=True)
print('done')
