"""Find [OIII] ridge peaks along RA in Dec bins for each registered image; save peak catalogue.
Output: data/derived/geometry_ridge_peaks.csv (source, dec, ra, amp(norm), prominence, width_deg(FWHM on sky))
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_ridges.py
"""
import numpy as np
from scipy import ndimage, signal

ROOT = '/Users/mickey/astronomy/silk-nebula'
z = np.load(f'{ROOT}/data/derived/geometry_maps.npz')
ra, de = z['ra_grid'], z['dec_grid']
STEP = abs(ra[1] - ra[0])
SRC = {
    'swons': z['swons_B'] - z['swons_R'],                 # [OIII]-blue minus Halpha-red colour index
    'fig4': -z['fig4_grey'],                              # paper Fig 4 [OIII] filter background (emission dark)
    'nsnsN': z['nsnsNs_OIII'],                            # NSNS narrow, [OIII] channel (stars opened)
    'nsnsW': z['nsnsW_OIII'],                             # NSNS wide, [OIII] channel (no continuum)
}
rows = []
DBIN = 0.2
for name, m in SRC.items():
    for d0 in np.arange(9.1, 30.0, DBIN):
        j = (de >= d0 - DBIN / 2) & (de < d0 + DBIN / 2)
        with np.errstate(all='ignore'):
            p = np.nanmean(m[j], axis=0)
        good = np.isfinite(p)
        if good.sum() < 100:
            continue
        # restrict to search box RA 113-120 and require coverage
        sel = good & (ra < 120.0) & (ra > 113.0)
        if sel.sum() < 150:
            continue
        pp = np.interp(np.arange(len(p)), np.where(good)[0], p[good])
        bg = ndimage.median_filter(pp, size=int(1.6 / STEP))            # running median over 1.6 deg (RA)
        r = ndimage.gaussian_filter1d(pp - bg, 1.0)
        noise = 1.4826 * np.median(np.abs(r[sel] - np.median(r[sel])))
        pk, pr = signal.find_peaks(r, prominence=4 * noise, width=2)
        cosd = np.cos(np.radians(d0))
        for k, i in enumerate(pk):
            if not sel[i]:
                continue
            # sub-pixel centre by parabola
            if 0 < i < len(r) - 1:
                a, b, c = r[i - 1], r[i], r[i + 1]
                off = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0
            else:
                off = 0
            rac = ra[i] - off * STEP                      # ra decreases with index
            rows.append((name, round(d0, 2), rac, r[i] / noise, pr['prominences'][k] / noise,
                         pr['widths'][k] * STEP * cosd))
import csv
with open(f'{ROOT}/data/derived/geometry_ridge_peaks.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['source', 'dec_deg', 'ra_deg', 'peak_snr', 'prominence_snr', 'fwhm_sky_deg'])
    for r in rows:
        w.writerow([r[0], f'{r[1]:.2f}', f'{r[2]:.4f}', f'{r[3]:.1f}', f'{r[4]:.1f}', f'{r[5]:.3f}'])
print(len(rows), 'peaks')
