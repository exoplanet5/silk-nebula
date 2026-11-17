"""verify_xrayfuv: independent block-averaging of the 267 GALEX GR6/7 pipeline skybg maps
(data/raw/galex_mast/skybg) to 1-arcmin blocks (40x40 px of 1.5"), keeping blocks fully valid
(all pixels > 0) within 0.55 deg of the tile centre. Values in CU via the Morrissey+2007
conversion (FUV 1.40e-15, NUV 2.06e-16 erg s-1 cm-2 A-1 per cps at 1528/2271 A; 1.5" pixels).
Output: data/derived/verify_xrayfuv_galex_blocks.npz
"""
import sys, glob, os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import ROOT, DER

HC = 1.98645e-8  # erg A
OMEGA = (1.5 / 206264.806) ** 2
CU_PER_CPS = {'FUV': 1.40e-15 / (HC / 1528.0) / OMEGA, 'NUV': 2.06e-16 / (HC / 2271.0) / OMEGA}
print('CU per cps/pixel', CU_PER_CPS)
meta = Table.read(ROOT / 'data/raw/galex_mast/galex_skybg.ecsv')
mb = {os.path.basename(str(u)): row for u, row in zip(meta['access_url'], meta)}
files = sorted(glob.glob(str(ROOT / 'data/raw/galex_mast/skybg/*-skybg.fits.gz')))
B = 40
rec = {k: [] for k in ['tile', 'band', 'ra', 'dec', 'cu', 'rtile', 'exptime', 'survey']}
tiles = []
for it, f in enumerate(files):
    bn = os.path.basename(f)
    band = 'FUV' if '-fd-' in bn else 'NUV'
    with fits.open(f) as h:
        d = h[0].data.astype(np.float64); hdr = h[0].header
    w = WCS(hdr)
    ny, nx = d.shape
    nby, nbx = ny // B, nx // B
    blk = d[:nby * B, :nbx * B].reshape(nby, B, nbx, B)
    valid = (blk > 0).all(axis=(1, 3))
    mean = blk.mean(axis=(1, 3))
    yc, xc = np.mgrid[0:nby, 0:nbx]
    px = xc * B + (B - 1) / 2.0; py = yc * B + (B - 1) / 2.0
    cx, cy = (hdr['CRPIX1'] - 1), (hdr['CRPIX2'] - 1)
    rt = np.hypot(px - cx, py - cy) * abs(hdr['CDELT1'])
    keep = valid & (rt <= 0.55)
    ra, dec = w.pixel_to_world_values(px[keep], py[keep])
    row = mb.get(bn)
    expt = float(row['t_exptime']) if row is not None else np.nan
    surv = str(row['survey']) if row is not None else '?'
    n = keep.sum()
    rec['tile'].append(np.full(n, it)); rec['band'].append(np.full(n, band))
    rec['ra'].append(ra); rec['dec'].append(dec); rec['cu'].append(mean[keep] * CU_PER_CPS[band])
    rec['rtile'].append(rt[keep]); rec['exptime'].append(np.full(n, expt)); rec['survey'].append(np.full(n, surv))
    tiles.append(bn)
    if it % 25 == 0:
        print(it, bn, band, n, np.median(mean[keep]) * CU_PER_CPS[band], expt, surv, flush=True)
out = {k: np.concatenate(v) for k, v in rec.items()}
np.savez_compressed(DER / 'verify_xrayfuv_galex_blocks.npz', tiles=np.array(tiles), **out)
print('blocks', len(out['ra']), 'FUV', (out['band'] == 'FUV').sum(), 'NUV', (out['band'] == 'NUV').sum())
