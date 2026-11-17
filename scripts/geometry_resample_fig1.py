"""Resample the registered Fig 1 SWONS composite (img-001) onto a plate-carree RA/Dec grid (0.01 deg) and build
star-suppressed colour-index maps: OIII_index = B - R (cyan [OIII] vs red Halpha), HA_index = R - B.
Output: data/derived/geometry_fig1_resampled.npz  (ra_grid, dec_grid, rgb, rgb_starless, valid)
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_resample_fig1.py
"""
import sys
import numpy as np
from PIL import Image
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from geometry_fig1_tools import load_fig1_map, W1, H1

ROOT = '/Users/mickey/astronomy/silk-nebula'
Image.MAX_IMAGE_PIXELS = None
pm = load_fig1_map()
rgb = np.asarray(Image.open(f'{ROOT}/data/raw/paper_figures/img-001.jpg').convert('RGB')).astype(np.float32)
# star suppression: grey-opening removes compact bright features (stars) while keeping diffuse nebulosity
starless = np.stack([ndimage.grey_opening(ndimage.median_filter(rgb[..., k], size=5), size=(9, 9)) for k in range(3)], axis=-1)
step = 0.01
ra_g = np.arange(119.0, 114.2, -step)
de_g = np.arange(12.7, 26.45, step)
RA, DE = np.meshgrid(ra_g, de_g)
# initial guess for inverse: linear interpolation from coarse forward grid
px, py = pm.inverse(RA.ravel(), DE.ravel(), guess_xy=((119.0 - RA.ravel()) / 4.3 * W1, (26.3 - DE.ravel()) / 13.5 * H1))
px = px.reshape(RA.shape); py = py.reshape(RA.shape)
valid = (px >= 0) & (px <= W1 - 1) & (py >= 0) & (py <= H1 - 1)
out = np.zeros(RA.shape + (3,), np.float32); outs = np.zeros_like(out)
for k in range(3):
    out[..., k] = ndimage.map_coordinates(rgb[..., k], [py, px], order=1, cval=np.nan)
    outs[..., k] = ndimage.map_coordinates(starless[..., k], [py, px], order=1, cval=np.nan)
out[~valid] = np.nan; outs[~valid] = np.nan
np.savez_compressed(f'{ROOT}/data/derived/geometry_fig1_resampled.npz', ra_grid=ra_g, dec_grid=de_g, rgb=out,
                    rgb_starless=outs, valid=valid)
# round-trip check
r2, d2 = pm(px[valid][::997], py[valid][::997])
err = np.hypot((r2 - RA[valid][::997]) * np.cos(np.radians(DE[valid][::997])), d2 - DE[valid][::997]) * 3600
print('grid', RA.shape, 'valid frac %.2f' % valid.mean(), 'inverse round-trip max err %.3f arcsec' % err.max())
