#!/usr/bin/env python
"""spectra_ism: empirical wavelength-frame test for the archival spectra with telluric O2 gamma-band lines.

Reference: UVES HD 61717 red-arm product (SPECSYS = TOPOCENT), so telluric lines sit at v_topo = 0 there.
For every other spectrum the high-pass-filtered 6276-6312 A region (O2 gamma band) is cross-correlated with
the reference; the shift (km/s) at which telluric lines appear is compared with the header BERV:
shift ~ 0 -> topocentric wavelengths; shift ~ +BERV -> barycentric-corrected wavelengths.
Output: data/derived/spectra_ism_frame_test.csv
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_frame_test.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER  # noqa: E402
from spectra_ism_speclib import load, C  # noqa: E402

SP = RAW / 'spectra_ism_spectra'
W1, W2 = 6276.0, 6312.0
LOGW = np.arange(np.log(W1), np.log(W2), 1.0 / C)        # 1 km/s pixels


def prep(sp):
    m = (sp['wave'] > W1 - 5) & (sp['wave'] < W2 + 5)
    if m.sum() < 200:
        return None
    w, f = sp['wave'][m], sp['flux'][m]
    f = np.interp(LOGW, np.log(w), f)
    cont = median_filter(f, size=121)                     # ~120 km/s running median
    x = f / np.where(cont > 0, cont, np.nan) - 1
    x = np.nan_to_num(x - np.nanmean(x))
    return x


ref = prep(load(SP / 'eso_HD61717_ADP.2024-03-09T06:44:09.481.fits'))
lags = np.arange(-60, 61)
rows = []
for p in sorted(SP.glob('*.fits')):
    sp = load(p)
    x = prep(sp)
    if x is None:
        continue
    cc = np.array([np.sum(ref[60:-60] * np.roll(x, -k)[60:-60]) for k in lags])
    cc /= np.sqrt(np.sum(ref[60:-60] ** 2) * np.sum(x[60:-60] ** 2))
    k = np.argmax(cc)
    # parabolic refinement
    if 0 < k < len(lags) - 1:
        a, b, c = cc[k - 1], cc[k], cc[k + 1]
        off = 0.5 * (a - c) / (a - 2 * b + c)
    else:
        off = 0.0
    shift = lags[k] + off
    rows.append(dict(file=p.name, inst=sp['inst'], date=sp['date'], header_frame=sp['frame'], berv=round(sp['berv'], 3),
                     telluric_shift_kms=round(float(shift), 2), ccmax=round(float(cc[k]), 3),
                     inferred=('topo' if abs(shift) < 3 else ('bary' if abs(shift - sp['berv']) < 3 else 'unclear'))))
    print(f"{p.name[:48]:48s} {sp['inst']:7s} hdr={sp['frame']} BERV={sp['berv']:+7.2f} telluric shift={shift:+7.2f} cc={cc[k]:.2f} -> {rows[-1]['inferred']}")
with open(DER / 'spectra_ism_frame_test.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
