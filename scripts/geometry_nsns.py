"""NSNS (Ziegenbalg, simg.de) image access with the WCS published in the viewer data files
(gem-01-ohs.js: STG, centre 116.0/+18.5, 10"/px full res; gem-01c-ohs.js: STG, centre 83.0/+13.5, 10"/px).
Each WCS is checked/refined against Tycho-2 (small pixel offset fit). Channel mapping: R=[OIII], G=Halpha, B=[SII].
"""
import sys
import numpy as np
from PIL import Image
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from geometry_imreg import load_tycho, blobs

ROOT = '/Users/mickey/astronomy/silk-nebula'
RAW = f'{ROOT}/data/raw/nsns_simg'


def make_wcs(crval, crpix_full, level, H):
    """WCS for a level-`level` downsampled image whose array rows run top->bottom (we flip to FITS y-up)."""
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---STG', 'DEC--STG']
    w.wcs.crval = crval
    w.wcs.crpix = [(crpix_full[0] - 0.5) / level + 0.5, (crpix_full[1] - 0.5) / level + 0.5]
    w.wcs.cd = [[-10 * level / 3600, 0], [0, 10 * level / 3600]]
    return w


def load_narrow():
    img = np.asarray(Image.open(f'{RAW}/gem-01-ohs-2.jpg').convert('RGB')).astype(np.float32)
    w = make_wcs([np.degrees(2.024581932), np.degrees(.3228859116)], (1575.5, 2575.5), 2, img.shape[0])
    return img, w


def load_wide():
    rows, cols = [(0, 1105), (1104, 2210), (2209, 3316)], [(0, 1072), (1071, 2145)]
    big = np.zeros((3317, 2146, 3), np.float32)
    for d, (u, h) in enumerate(rows):
        for p, (g, m) in enumerate(cols):
            t = np.asarray(Image.open(f'{RAW}/wide_l4/4-{d}-{p}.jpg').convert('RGB')).astype(np.float32)
            big[u:u + t.shape[0], g:g + t.shape[1]] = t[:min(t.shape[0], 3317 - u), :min(t.shape[1], 2146 - g)]
    w = make_wcs([np.degrees(1.448623279), np.degrees(.235619449)], (15015.5, 11054.5), 4, 5527)
    return big, w, 5527   # full level-4 image height (rows from top) for FITS y conversion


def pix_to_world(w, col, row, Hfull):
    """col,row: 0-based array indices (row from top). FITS y = Hfull - row."""
    return w.wcs_pix2world(col + 1, Hfull - row, 1)


def world_to_pix(w, ra, de, Hfull):
    x, y = w.wcs_world2pix(ra, de, 1)
    return x - 1, Hfull - y


def check_and_refine(img, w, Hfull, vtmax=9.0, nmax=500, k=8, label=''):
    """Detect stars, match to Tycho-2, fit a constant pixel offset; return (dx, dy, stats)."""
    L = img.mean(2)
    r = L - ndimage.median_filter(L, size=15)
    s = 1.4826 * np.median(np.abs(r - np.median(r)))
    b = blobs(r > k * s, min_area=2, weight=np.clip(r, 0, None))
    b = b[np.argsort(-b[:, 3])][:nmax]
    cat = load_tycho(vtmax)
    cx, cy = world_to_pix(w, cat[:, 0], cat[:, 1], Hfull)
    t = cKDTree(np.c_[cx, cy])
    dx = dy = 0.0
    for tol in (4, 2.5, 1.5):
        d, i = t.query(np.c_[b[:, 0] - dx, b[:, 1] - dy])
        m = d < tol
        dx = float(np.median(b[m, 0] - cx[i[m]])); dy = float(np.median(b[m, 1] - cy[i[m]]))
    d, i = t.query(np.c_[b[:, 0] - dx, b[:, 1] - dy]); m = d < 1.5
    pixsc = abs(w.wcs.cd[0][0]) * 3600
    st = dict(n_det=len(b), n_match=int(m.sum()), dx_px=dx, dy_px=dy, median_resid_arcsec=float(np.median(d[m]) * pixsc),
              rms_resid_arcsec=float(np.sqrt(np.mean(d[m] ** 2)) * pixsc), pixscale_arcsec=pixsc)
    return dx, dy, st


def resample(img, w, Hfull, ra_g, de_g, dx=0.0, dy=0.0, order=1):
    RA, DE = np.meshgrid(ra_g, de_g)
    col, row = world_to_pix(w, RA.ravel(), DE.ravel(), Hfull)
    col = col + dx; row = row + dy
    out = np.stack([ndimage.map_coordinates(img[..., k], [row, col], order=order, cval=np.nan).reshape(RA.shape)
                    for k in range(img.shape[2])], axis=-1)
    bad = (col < 0) | (col > img.shape[1] - 1) | (row < 0) | (row > img.shape[0] - 1)
    out[bad.reshape(RA.shape)] = np.nan
    return out
