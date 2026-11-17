"""Mosaic GALFA-HI DR2 wide-product cutouts onto the global 1' RA/Dec lattice and derive HI products.

Outputs (data/derived/imaging_mw_grids/):
  imaging_mw_galfa_cube3.fits     Tb cube block-averaged to 3' on the EQ grid (V_LSR -29.8..+40.1, 0.736 km/s)
  imaging_mw_eq_galfa_*.fits      NHI maps (optically thin, 1.823e18 * int Tb dv) on the EQ grid for velocity windows
Run: ~/.venvs/astro313/bin/python -I scripts/imaging_mw_hi_build.py
"""
import glob
import json
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import RAW, DER, EQ_WCS, EQ_SHAPE, EQ, save_grid

files = sorted(glob.glob(str(RAW / 'galfa_hi_dr2' / 'galfa_dr2_wide_*.fits')))
print(len(files), 'tiles')
# Global 1' lattice: EQ grid at 1' resolution, RA 108-125 (1020 px), Dec 4-30 (1560 px).
# pixel centre RA_i = 125 - (i + 0.5)/60 ; Dec_j = 4 + (j + 0.5)/60 (same lattice as GALFA's: centres at k/60 + 1/120)
NX, NY = 1020, 1560
cube = None
cover = np.zeros((NY, NX), bool)
for f in files:
    h = fits.open(f)
    hd = h[0].header
    d = h[0].data.astype(np.float32)  # astropy applies BSCALE; BLANK -> nan for int data
    nv = d.shape[0]
    if cube is None:
        cube = np.full((nv, NY, NX), np.nan, np.float32)
        v = (np.arange(nv) + 1 - hd['CRPIX3']) * hd['CDELT3'] / 1000.0 + hd['CRVAL3'] / 1000.0
    ra0 = hd['CRVAL1'] + (1 - hd['CRPIX1']) * hd['CDELT1']     # RA of tile pixel 0 centre
    de0 = hd['CRVAL2'] + (1 - hd['CRPIX2']) * hd['CDELT2']
    i0 = int(round((125 - ra0) * 60 - 0.5)); j0 = int(round((de0 - 4) * 60 - 0.5))
    ny, nx = d.shape[1:]
    # clip to mosaic
    xs, ys = slice(max(i0, 0), min(i0 + nx, NX)), slice(max(j0, 0), min(j0 + ny, NY))
    txs = slice(xs.start - i0, xs.stop - i0); tys = slice(ys.start - j0, ys.stop - j0)
    cube[:, ys, xs] = np.where(np.isfinite(d[:, tys, txs]), d[:, tys, txs], cube[:, ys, xs])
    cover[ys, xs] = True
    print(f.split('/')[-1], i0, j0, 'ra0', round(ra0, 4), 'de0', round(de0, 4))
dv = abs(v[1] - v[0])
print('velocity', v[0], v[-1], dv, 'coverage', cover.mean())

# 3' block average onto the EQ grid (340 x 520)
c3 = np.nanmean(cube.reshape(len(v), NY // 3, 3, NX // 3, 3), axis=(2, 4))
# EQ grid x runs with decreasing RA from x=0 at RA 124.975: our 1' lattice i=0 is RA 124.99 -> same orientation
hdr = EQ_WCS.to_header()
hdr['CTYPE3'] = 'VRAD'; hdr['CRVAL3'] = float(v[0]); hdr['CDELT3'] = float(v[1] - v[0]); hdr['CRPIX3'] = 1.0
hdr['CUNIT3'] = 'km/s'; hdr['BUNIT'] = 'K'; hdr['COMMENT'] = 'GALFA-HI DR2 wide Tb, 3 arcmin block average, V_LSR'
out = DER / 'imaging_mw_grids'
fits.PrimaryHDU(c3.astype(np.float32), header=hdr).writeto(out / 'imaging_mw_galfa_cube3.fits', overwrite=True)

WINDOWS = {  # name: (vlo, vhi) km/s LSR
    'galfa_nhi_all': (-29.8, 40.1),
    'galfa_nhi_oiii': (11.0, 16.0),     # [OIII] V_LSR range (Wu+26: +11.2..+15.5)
    'galfa_nhi_red': (-1.0, 6.0),       # Halpha/[NII]/[SII] V_LSR range (Wu+26: -1.0..+6.1)
    'galfa_nhi_m30m10': (-30, -10), 'galfa_nhi_m10m02': (-10, -2), 'galfa_nhi_m02p06': (-2, 6),
    'galfa_nhi_p06p10': (6, 10), 'galfa_nhi_p10p18': (10, 18), 'galfa_nhi_p18p26': (18, 26), 'galfa_nhi_p26p40': (26, 40),
}
meta = {}
for name, (lo, hi) in WINDOWS.items():
    k = (v >= lo) & (v <= hi)
    nhi = 1.823e18 * np.nansum(c3[k], axis=0) * dv
    nhi[~np.isfinite(c3[k]).any(axis=0)] = np.nan
    save_grid(nhi, name, 'eq', 'cm-2', f'GALFA NHI {lo}..{hi} km/s')
    meta[name] = [lo, hi, int(k.sum())]
# moment 1 (intensity weighted velocity) and velocity of peak, for the -10..+30 range
k = (v >= -10) & (v <= 30)
w = np.clip(c3[k], 0, None)
m1 = np.nansum(w * v[k, None, None], axis=0) / np.nansum(w, axis=0)
vpk = v[k][np.nanargmax(np.where(np.isfinite(c3[k]), c3[k], -1), axis=0)]
save_grid(m1, 'galfa_mom1', 'eq', 'km/s', 'GALFA moment-1 V_LSR -10..+30')
save_grid(vpk.astype(float), 'galfa_vpeak', 'eq', 'km/s', 'GALFA velocity of peak Tb, -10..+30')
# also keep a 1' NHI map of the full window for high-resolution profiles
nhi1 = 1.823e18 * np.nansum(cube, axis=0) * dv
nhi1[~cover] = np.nan
h1 = fits.Header()
h1['CTYPE1'] = 'RA---CAR'; h1['CTYPE2'] = 'DEC--CAR'; h1['CRVAL1'] = 116.5; h1['CRVAL2'] = 0.0
h1['CDELT1'] = -1 / 60; h1['CDELT2'] = 1 / 60; h1['CRPIX1'] = NX / 2 + 0.5; h1['CRPIX2'] = 1 - (4 + 1 / 120) * 60
h1['BUNIT'] = 'cm-2'; h1['COMMENT'] = 'GALFA-HI DR2 NHI -29.8..+40.1 km/s at 1 arcmin'
fits.PrimaryHDU(nhi1.astype(np.float32), header=h1).writeto(out / 'imaging_mw_galfa_nhi_1arcmin.fits', overwrite=True)
(out / 'imaging_mw_galfa_meta.json').write_text(json.dumps({'tiles': [f.split('/')[-1] for f in files], 'windows': meta,
                                                             'v0': float(v[0]), 'dv': float(dv), 'nv': int(len(v)),
                                                             'coverage_frac_1arcmin': float(cover.mean())}, indent=1))
print('done')
