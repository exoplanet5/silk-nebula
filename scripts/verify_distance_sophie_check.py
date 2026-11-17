"""verify_distance: narrow-line (high-pass matched-filter) search for the HD 61717 -16.9 km/s Ca II component toward
the SOPHIE stars (barycentric s1d, wave-1 downloads).  SOPHIE starts at 3872.7 A, so Ca II K and H are covered.
Output: data/derived/verify_distance_sophie_check.json
Run: ~/.venvs/astro313/bin/python -I scripts/verify_distance_sophie_check.py
"""
import glob
import json
import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

ROOT = '/Users/mickey/astronomy/silk-nebula'
C = 299792.458
LINES = {'CaK': (3933.6614, 0.6267), 'CaH': (3968.4673, 0.3116)}
V0, B, FWHM = -16.88, 1.49, C / 75000.0
out = []
for f in sorted(glob.glob(f'{ROOT}/data/raw/spectra_ism_spectra/sophie_*_s1d.fits')):
    with fits.open(f, memmap=False) as h:
        H = h[0].header
        fl = np.array(h[0].data, float)
    w = H['CRVAL1'] + H['CDELT1'] * (np.arange(len(fl)) + 1 - H['CRPIX1'])
    row = dict(file=f.split('/')[-1], star=f.split('/')[-1].split('_')[1])
    for ln, (lam0, fv) in LINES.items():
        v = (w / lam0 - 1) * C               # barycentric (SOPHIE s1d; checked by wave-1 frame test)
        m = (v > -90) & (v < 90) & np.isfinite(fl)
        vv, ff = v[m], fl[m]
        npx = int(round(15.0 / np.median(np.diff(vv)))) | 1
        fn = ff / median_filter(ff, size=npx, mode='nearest')
        sig = np.sqrt(B ** 2 / 2 + (FWHM / 2.3548) ** 2)
        shape = np.exp(-0.5 * ((vv - V0) / sig) ** 2)
        cm = (np.abs(vv - V0) > 3 * FWHM) & ~((vv > -5) & (vv < 30))
        noise = np.std(fn[cm] - 1)
        sel = np.abs(vv - V0) < 4 * sig
        A = np.sum(((1 - fn) * shape)[sel]) / np.sum((shape ** 2)[sel])
        eA = noise / np.sqrt(np.sum((shape ** 2)[sel])) * 1.3
        area = np.sum(shape * np.gradient(vv)) / C * lam0 * 1000
        W, eW = A * area, eA * area
        Nup = 1.13e17 * (max(W, 0) + 3 * eW) / (fv * lam0 ** 2)
        row[ln] = dict(W_mA=round(float(W), 2), eW_mA=round(float(eW), 2), logN_3sig=round(float(np.log10(Nup)), 2))
    out.append(row)
    print(row)
json.dump(dict(note='high-pass (15 km/s running median) matched filter at V_helio -16.88, b 1.49, R 75000 assumed; '
                    'HD 61717 itself gives W(K)=19.7+-0.25 mA with the same filter', rows=out),
          open(f'{ROOT}/data/derived/verify_distance_sophie_check.json', 'w'), indent=1)
