"""Register the SWONS deep colour composite (Wu+2026 Fig 1; embedded image img-001, 1508x5274) to the sky by
matching detected stars to Tycho-2 (VizieR I/259, CfA mirror). Grid search over scale/offset, then
iterative polynomial (gnomonic) fit. Output: data/derived/geometry_fig1_wcs.json (+ matched star list).
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_register_fig1.py
"""
import json, sys
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from geometry_imreg import load_tycho, blobs, match_fit, PolyMap, gnomonic, inv_gnomonic

ROOT = '/Users/mickey/astronomy/silk-nebula'
Image.MAX_IMAGE_PIXELS = None


def detect_stars(rgb, hp=9, k=6.0, min_area=2):
    L = rgb.astype(float).mean(axis=2)
    bg = ndimage.median_filter(L, size=hp * 2 + 1)
    r = L - bg
    s = 1.4826 * np.median(np.abs(r - np.median(r)))
    b = blobs(r > k * s, min_area=min_area, weight=np.clip(r, 0, None))
    return b[np.argsort(-b[:, 3])], r, s


def grid_register(det, cat, ra0, de0, scales, parities=(1, -1), ndet=120, ncat=250, tol=0.06):
    """Coarse search: pixel -> tangent plane with xi = par*(x-x0)*sc (east=-x when par=1), eta = -(y-y0)*sc."""
    d = det[:ndet]
    c = cat[np.argsort(cat[:, 2])][:ncat]
    xi_c, eta_c = gnomonic(c[:, 0], c[:, 1], ra0, de0)
    tree = cKDTree(np.c_[xi_c, eta_c])
    best = (0, None)
    x0, y0 = det[:, 0].mean(), det[:, 1].mean()
    for par in parities:
        for sc in scales:
            xi = -par * (d[:, 0] - x0) * sc
            eta = -(d[:, 1] - y0) * sc
            # offsets search via histogram of all pair differences (fast voting)
            dd, ii = tree.query(np.c_[xi, eta], k=60, distance_upper_bound=6)
            ok = np.isfinite(dd)
            rows = np.repeat(np.arange(len(d))[:, None], 60, axis=1)[ok]
            dx = xi_c[ii[ok]] - xi[rows]; dy = eta_c[ii[ok]] - eta[rows]
            H, xe, ye = np.histogram2d(dx, dy, bins=[np.arange(-6, 6.01, tol), np.arange(-6, 6.01, tol)])
            j = np.unravel_index(np.argmax(H), H.shape)
            if H[j] > best[0]:
                best = (H[j], dict(par=par, sc=sc, dx=0.5 * (xe[j[0]] + xe[j[0] + 1]), dy=0.5 * (ye[j[1]] + ye[j[1] + 1]),
                                   x0=x0, y0=y0))
    return best


if __name__ == '__main__':
    rgb = np.asarray(Image.open(f'{ROOT}/data/raw/paper_figures/img-001.jpg').convert('RGB'))
    H, W, _ = rgb.shape
    det, resid, sig = detect_stars(rgb)
    print('detections', len(det), 'noise', sig)
    cat = load_tycho(10.0)
    ra0, de0 = 116.4, 16.4
    scales = 1 / np.arange(330, 520, 4)     # deg per pixel
    nvote, g = grid_register(det, cat, ra0, de0, scales)
    print('grid best votes', nvote, g)
    mapping = lambda x, y: inv_gnomonic(-g['par'] * (x - g['x0']) * g['sc'] + g['dx'],
                                        -(y - g['y0']) * g['sc'] + g['dy'], ra0, de0)
    out = {}
    for order in (1, 2, 3):
        pm, st, (dm, cm) = match_fit(det, cat, mapping, tols=(0.08, 0.05, 0.03, 0.02, 0.015), order=order,
                                     ra0=ra0, de0=de0, ndet=600, ncat=1500)
        print('order', order, st)
        out[order] = (pm, st, dm, cm)
    pm, st, dm, cm = out[3] if out[3][1]['n_match'] > 60 else out[2]
    # leave-one-out style check: residuals for each matched star
    ra, de = pm(dm[:, 0], dm[:, 1])
    res = np.hypot((ra - cm[:, 0]) * np.cos(np.radians(de)), de - cm[:, 1]) * 3600
    corners = {f'{x},{y}': [float(v[0]) for v in pm(np.array([x]), np.array([y]))] for x, y in
               [(0, 0), (W - 1, 0), (0, H - 1), (W - 1, H - 1), (W / 2, H / 2)]}
    # local pixel scale
    (r1, d1), (r2, d2) = pm(np.array([W / 2]), np.array([H / 2])), pm(np.array([W / 2]), np.array([H / 2 + 100]))
    scale = float(np.hypot((r2 - r1) * np.cos(np.radians(d1)), d2 - d1)[0] / 100 * 3600)
    js = dict(image='data/raw/paper_figures/img-001.jpg', shape=[H, W], ra0=ra0, de0=de0, order=pm.order,
              x0=pm.x0, y0=pm.y0, s=pm.s, cx=pm.cx.tolist(), cy=pm.cy.tolist(), stats=st, grid=dict(votes=float(nvote), **{k: float(v) for k, v in g.items()}),
              pixel_scale_arcsec=scale, corners_radec=corners, parity='east-left (normal sky)' if g['par'] == 1 else 'mirrored',
              resid_p50_p90_arcsec=[float(np.percentile(res, 50)), float(np.percentile(res, 90))])
    json.dump(js, open(f'{ROOT}/data/derived/geometry_fig1_wcs.json', 'w'), indent=1)
    np.savetxt(f'{ROOT}/data/derived/geometry_fig1_matched_stars.csv', np.c_[dm[:, :2], cm, res], delimiter=',', fmt='%.5f',
               header='x_px,y_px,tyc_ra,tyc_dec,tyc_vt,resid_arcsec', comments='')
    print(json.dumps({k: js[k] for k in ('pixel_scale_arcsec', 'corners_radec', 'parity', 'resid_p50_p90_arcsec')}, indent=1))
