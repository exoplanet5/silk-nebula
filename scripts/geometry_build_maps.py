"""Resample all registered images onto one plate-carree grid (RA 121->111, Dec 9->30, 0.02 deg) and save
data/derived/geometry_maps.npz. Sources: SWONS composite (Wu+2026 Fig 1, Tycho-registered), Wu+2026 Fig 4
[OIII] background (Tycho-registered), NSNS narrow (20"/px) and wide (40"/px) RGB = [OIII], Halpha, [SII].
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_build_maps.py
"""
import json, sys
import numpy as np
from PIL import Image
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from geometry_imreg import load_tycho, blobs, match_fit
from geometry_fig1_tools import load_fig1_map, W1, H1
import geometry_nsns as ns

ROOT = '/Users/mickey/astronomy/silk-nebula'
STEP = 0.02
ra_g = np.arange(121.0, 110.999, -STEP)
de_g = np.arange(9.0, 30.001, STEP)
RA, DE = np.meshgrid(ra_g, de_g)
out = {'ra_grid': ra_g, 'dec_grid': de_g}
meta = {}


def starless(ch, med=5, op=9):
    return ndimage.grey_opening(ndimage.median_filter(ch, size=med), size=(op, op))


# --- SWONS composite (Fig 1)
pm = load_fig1_map()
rgb = np.asarray(Image.open(f'{ROOT}/data/raw/paper_figures/img-001.jpg').convert('RGB')).astype(np.float32)
sl = np.stack([starless(rgb[..., k]) for k in range(3)], -1)
px, py = pm.inverse(RA.ravel(), DE.ravel(), guess_xy=((118.9 - RA.ravel()) / 4.3 * W1, (26.3 - DE.ravel()) / 13.5 * H1))
ok = (px >= 0) & (px <= W1 - 1) & (py >= 0) & (py <= H1 - 1)
for k, nm in enumerate('RGB'):
    a = ndimage.map_coordinates(sl[..., k], [py, px], order=1, cval=np.nan); a[~ok] = np.nan
    out[f'swons_{nm}'] = a.reshape(RA.shape).astype(np.float32)
meta['swons'] = 'Wu+2026 Fig 1 (img-001), SWONS H-alpha(red)+[OIII](blue) composite, starless (median5+opening9)'

# --- Fig 4 [OIII] filter background (gray, stars dark), flip rows to display orientation, register to Tycho
g4 = np.asarray(Image.open(f'{ROOT}/data/raw/paper_figures/img-004.png')).astype(np.float32)[::-1]
Hh, Ww = g4.shape
inv = 255 - g4
b = blobs(ndimage.uniform_filter(inv, 3) > 200, min_area=6, weight=inv); b = b[np.argsort(-b[:, 2])]
cat = load_tycho(10)
c4 = cat[(cat[:, 0] < 118.8) & (cat[:, 0] > 114.6) & (cat[:, 1] > 12.8) & (cat[:, 1] < 20)]
mp = lambda x, y: (118.80 - (x + 0.5) * 4.2 / Ww, 20.0 - (y + 0.5) * 7.2 / Hh)
pm4, st4, (dm4, cm4) = match_fit(b, c4, mp, ndet=250, ncat=300, order=2)
meta['fig4_registration'] = st4
# nebula map: light pixels = emission (image shows emission as darker? -> check sign later); keep raw grey
g4s = ndimage.grey_closing(ndimage.median_filter(g4, size=5), size=(7, 7))   # remove dark star specks
qx, qy = pm4.inverse(RA.ravel(), DE.ravel(), guess_xy=((118.8 - RA.ravel()) / 4.2 * Ww, (20 - DE.ravel()) / 7.2 * Hh))
ok = (qx >= 0) & (qx <= Ww - 1) & (qy >= 0) & (qy <= Hh - 1)
a = ndimage.map_coordinates(g4s, [qy, qx], order=1, cval=np.nan); a[~ok] = np.nan
out['fig4_grey'] = a.reshape(RA.shape).astype(np.float32)
# axis-implied vs star-implied positions: how misregistered is the paper's background image?
xs, ys = np.meshgrid(np.linspace(0, Ww - 1, 15), np.linspace(0, Hh - 1, 25))
ra_s, de_s = pm4(xs.ravel(), ys.ravel()); ra_a, de_a = mp(xs.ravel(), ys.ravel())
off = np.hypot((ra_s - ra_a) * np.cos(np.radians(de_s)), de_s - de_a) * 60
meta['fig4_axis_vs_star_offset_arcmin'] = {'median': float(np.median(off)), 'max': float(off.max()),
                                           'mean_dRA_deg': float(np.mean(ra_s - ra_a)), 'mean_dDec_deg': float(np.mean(de_s - de_a))}

# --- NSNS narrow and wide
img, w = ns.load_narrow(); dx, dy, st = ns.check_and_refine(img, w, img.shape[0]); meta['nsns_narrow_wcs_check'] = st
a = ns.resample(img, w, img.shape[0], ra_g, de_g, dx, dy)
for k, nm in enumerate(['OIII', 'Ha', 'SII']):
    out[f'nsnsN_{nm}'] = a[..., k].astype(np.float32)
big, w2, H4 = ns.load_wide(); dx, dy, st = ns.check_and_refine(big, w2, H4, vtmax=8); meta['nsns_wide_wcs_check'] = st
# star removal at native resolution before resampling
bigs = np.stack([ndimage.grey_opening(ndimage.median_filter(big[..., k], 3), size=(3, 3)) for k in range(3)], -1)
a = ns.resample(bigs, w2, H4, ra_g, de_g, dx, dy)
for k, nm in enumerate(['OIII', 'Ha', 'SII']):
    out[f'nsnsW_{nm}'] = a[..., k].astype(np.float32)
imgs = np.stack([ndimage.grey_opening(ndimage.median_filter(img[..., k], 3), size=(5, 5)) for k in range(3)], -1)
a = ns.resample(imgs, w, img.shape[0], ra_g, de_g, dx=0, dy=0)
for k, nm in enumerate(['OIII', 'Ha', 'SII']):
    out[f'nsnsNs_{nm}'] = a[..., k].astype(np.float32)
np.savez_compressed(f'{ROOT}/data/derived/geometry_maps.npz', **out)
json.dump(meta, open(f'{ROOT}/data/derived/geometry_maps_meta.json', 'w'), indent=1, default=float)
print(json.dumps(meta, indent=1, default=float))
