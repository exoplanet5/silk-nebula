"""verify_xrayfuv: empirical GALEX FUV/NUV response to dust (dust-scattered light) in the Silk field,
to bound a dust-scattering origin of the on-filament UV excess.
Regress tile-demeaned GALEX 1' blocks (verify_xrayfuv_galex_blocks.npz, same flat/star treatment
as verify_xrayfuv_galex_measure.py) on the high-passed IRIS 100 um map (wave-1 EQ grid, read only)
-> slope (CU per MJy/sr); with RA-shift nulls. Also the median diffuse FUV/NUV colour.
Output: data/derived/verify_xrayfuv_galex_dust.json
"""
import sys, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_xrayfuv_galex_measure as M   # re-uses prep_band (runs its main block on import is avoided below)
from verify_xrayfuv_common import *

ir = fits.open(DER / 'imaging_mw_grids/imaging_mw_eq_iris100.fits')[0]
wi = WCS(ir.header)
I = np.nan_to_num(ir.data.astype(float), nan=np.nanmedian(ir.data))
pix_deg = abs(ir.header.get('CDELT1', 0.05))
Ihp = ndimage.gaussian_filter(I - ndimage.median_filter(I, size=int(round(1.0 / pix_deg)) | 1, mode='nearest'), 0.05 / pix_deg)


def samp(ra, dec, dra=0.0):
    x, y = wi.world_to_pixel_values(ra - dra / np.cos(np.radians(dec)), dec)
    return ndimage.map_coordinates(Ihp, [y, x], order=1, mode='constant', cval=np.nan)


out = {}
for band in ['FUV', 'NUV']:
    D = M.prep_band(band)
    b0, _, n = M.regress(D, samp(D['ra'], D['dec']))
    nb = np.array([M.regress(D, samp(D['ra'], D['dec'], d))[0] for d in M.np.round(np.concatenate([np.arange(-3, -0.69, 0.1), np.arange(0.7, 3.01, 0.1)]), 2)])
    out[band] = {'slope_CU_per_MJysr': b0, 'null_mean': float(nb.mean()), 'null_std': float(nb.std(ddof=1)),
                 'z': float((b0 - nb.mean()) / nb.std(ddof=1)), 'n_blocks': n, 'median_level_CU': D['median_level']}
    print(band, out[band], flush=True)
out['diffuse_FUV_over_NUV_median'] = out['FUV']['median_level_CU'] / out['NUV']['median_level_CU']
# bound on dust-scattered excess for the wave-1 IRIS 100um 2-sigma limit on the filament (0.19 MJy/sr)
for band in ['FUV', 'NUV']:
    s = abs(out[band]['slope_CU_per_MJysr']) + 2 * out[band]['null_std']
    out[band]['max_dust_scattered_CU_for_dI100_0.19'] = float(s * 0.19)
(DER / 'verify_xrayfuv_galex_dust.json').write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
