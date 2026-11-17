"""Point-source masking and smoothing of X-ray grids (EQ and GAL), plus RASS band sums.

eROSITA half-sky rate maps (3' pixels, ~1 count/pixel) and RASS3 intensity maps include point sources:
pixels exceeding the local 21'x21' median by > 4 robust sigma (and their 8 neighbours) are masked, then the
map is smoothed by normalised Gaussian convolution (sigma = 3 px on EQ = 9', FWHM 21'; sigma = 1.5 px on GAL).
RASS diffuse maps (Snowden+97) are already point-source cleaned: R12 = R1+R2, R45 = R4+R5, R67 = R6+R7.
Outputs: imaging_mw_{eq,gal}_<band>_clean.fits (masked, smoothed) and *_r12/_r45/_r67.
"""
import sys

import numpy as np
from scipy.ndimage import median_filter, gaussian_filter, binary_dilation

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import load_grid, save_grid, DER


def clean(a, box=7, nsig=4.0, sigma_px=3.0):
    med = median_filter(np.nan_to_num(a, nan=np.nanmedian(a)), size=box)
    resid = a - med
    mad = median_filter(np.abs(np.nan_to_num(resid)), size=box * 2 + 1) * 1.4826
    src = (resid > nsig * np.maximum(mad, 1e-12)) & np.isfinite(a)
    src = binary_dilation(src, iterations=1)
    good = np.isfinite(a) & ~src
    w = gaussian_filter(good.astype(float), sigma_px)
    s = gaussian_filter(np.where(good, a, 0.0), sigma_px) / np.where(w > 0, w, np.nan)
    s[w < 0.25] = np.nan
    return s, src.mean()


for which, sig in [('eq', 3.0), ('gal', 1.5)]:
    for band in ['ero_0204', 'ero_0406', 'ero_0610', 'ero_1023', 'ero_2350', 'ero_ovii', 'ero_oviii',
                 'rass_int_soft', 'rass_int_hard', 'rass_int_broad']:
        try:
            a = load_grid(band, which)
        except FileNotFoundError:
            continue
        if band in ('ero_ovii', 'ero_oviii'):   # sparse counts maps: no source masking, heavier smoothing
            s, f = clean(a, nsig=1e9, sigma_px=2 * sig)
        else:
            s, f = clean(a, sigma_px=sig)
        save_grid(s, band + '_clean', which, 'as input', f'{band} point sources masked, Gaussian sigma {sig} px')
        print(which, band, f'masked frac {f:.3f}', f'median {np.nanmedian(s):.4g}')
    try:
        r = {i: load_grid(f'rass_r{i}', which) for i in [1, 2, 4, 5, 6, 7]}
    except FileNotFoundError:
        continue
    for nm, (i, j) in {'r12': (1, 2), 'r45': (4, 5), 'r67': (6, 7)}.items():
        save_grid(r[i] + r[j], 'rass_' + nm, which, '1e-6 cts/s/arcmin2', f'RASS diffuse {nm.upper()} (Snowden+97)')
        print(which, nm, np.nanmedian(r[i] + r[j]))
    # eROSITA soft sum 0.2-0.6 keV and 0.6-1.0 already; make 0.2-0.6 sum
    try:
        e = load_grid('ero_0204_clean', which) + load_grid('ero_0406_clean', which)
        save_grid(e, 'ero_0206_clean', which, 'cts/s/deg2', 'eRASS1 0.2-0.6 keV cleaned+smoothed')
    except FileNotFoundError:
        pass
