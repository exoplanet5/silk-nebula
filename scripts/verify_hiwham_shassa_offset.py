"""verify_hiwham: SHASSA limit on an H-alpha filament that is DISPLACED from the [O III] ridge, or broader than it.

Motivation: every aperture-matched H-alpha limit (DESI template regression, SDSS ON-OFF with r_on = 0.15 deg,
the d = 0 template fit in verify_hiwham_hamaps.py) assumes the H-alpha shares the [O III] cross-profile. In a
radiative shock the H recombination zone lies downstream of the [O III] zone. So here the Gaussian template is
centred at d0 = -0.8..+0.8 deg (0.05 steps), with FWHM 0.28 deg (= [O III] filament (+) 4' beam) and 0.5 deg.

Data: SkyView 'SHASSA Sm' (continuum-subtracted, smoothed ~4', dR -> R), spine piece Dec 12-15.6 (SHASSA
coverage). Profile: mean in 0.05-deg bins of signed offset d, |d| < 2.0. Fit: linear baseline + Gaussian
amplitude over |d| < 1.6. Nulls: the same spine piece at random centres / orientations (0-360 deg) / E-W mirror
within the SHASSA footprint, avoiding the real |d| < 1.5 deg band (< 5 % overlap), >= 85 % finite pixels.
Upper limit on the Gaussian PEAK: 2-sigma UL = max(A - med_null, 0) + 2 * sigma_robust(null), sigma_robust =
(p84 - p16)/2, per (d0, FWHM). The look-elsewhere-free statement is the per-d0 UL; the global significance of
the largest real |z| is calibrated against the same scan on the nulls.

Also merges all verify_hiwham_hamaps_part_*.json into verify_hiwham_hamaps_results.json (the earlier merge was
written before the SHASSA Sm part finished).
Output: data/derived/verify_hiwham_shassa_offset.json, verify_hiwham_shassa_offset.npz
"""
import glob
import json
import sys
import time
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import RAW, DER, strand, Trace

warnings.filterwarnings('ignore')
rng = np.random.default_rng(2026)
t0 = time.time()

# ------------------------------------------------------------------ merge earlier part files
merged = {}
for fpart in sorted(glob.glob(str(DER / 'verify_hiwham_hamaps_part_*.json'))):
    merged.update(json.loads(open(fpart).read()))
(DER / 'verify_hiwham_hamaps_results.json').write_text(json.dumps(merged, indent=1))
print('merged hamaps keys:', list(merged.keys()), flush=True)

# ------------------------------------------------------------------ data
f = RAW / 'skyview/SHASSA_Sm/SHASSA_Sm_wide.fits'
hdr = fits.getheader(f)
assert hdr.get('BUNIT', '').strip() == 'dR', hdr.get('BUNIT')
img = fits.getdata(f).astype(float) / 10.0          # dR -> R
w = WCS(hdr)
ny, nx = img.shape
yy, xx = np.mgrid[0:ny, 0:nx]
rr, dd = w.all_pix2world(xx, yy, 0)
val = img.ravel().copy(); pra = rr.ravel(); pdec = dd.ravel()
val[val == 0] = np.nan                               # SkyView blank outside the survey

ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
m = (REAL.dec >= 12.0) & (REAL.dec <= 15.6)
r_p, d_p = REAL.ra[m], REAL.dec[m]
PIECE = Trace(r_p, d_p)

EDG = np.arange(-2.0, 2.0 + 1e-9, 0.05); CEN = 0.5 * (EDG[1:] + EDG[:-1]); NB = len(CEN)


def profile(tr, pix):
    d, s, ins = tr.coords(pra[pix], pdec[pix], dmax=2.0)
    x = val[pix]
    ok = ins & np.isfinite(d) & np.isfinite(x) & (np.abs(d) < 2.0)
    b = ((d[ok] + 2.0) / 0.05).astype(int)
    n = np.bincount(b, minlength=NB)
    p = np.bincount(b, weights=x[ok], minlength=NB) / np.maximum(n, 1)
    p[n < 5] = np.nan
    return p


D0 = np.round(np.arange(-0.8, 0.8 + 1e-9, 0.05), 3)
FW = (0.284, 0.5)
FITW = np.abs(CEN) < 1.6


def scan(p):
    """Gaussian amplitude at each (FWHM, d0) over a linear baseline; returns array [len(FW), len(D0)]."""
    out = np.full((len(FW), len(D0)), np.nan)
    ok = FITW & np.isfinite(p)
    if ok.sum() < 30:
        return out
    x, y = CEN[ok], p[ok]
    for i, fw in enumerate(FW):
        sg = fw / 2.3548
        for j, d0 in enumerate(D0):
            X = np.c_[np.ones(ok.sum()), x, np.exp(-0.5 * ((x - d0) / sg) ** 2)]
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            out[i, j] = c[2]
    return out


# real
near_real = np.where((np.abs(pra - np.mean(r_p)) * np.cos(np.radians(np.mean(d_p))) < 4.0) & (np.abs(pdec - np.mean(d_p)) < 4.5))[0]
p_real = profile(PIECE, near_real)
A_real = scan(p_real)

# exclusion band of the whole real spine
dall, _, iall = REAL.coords(pra, pdec, dmax=2.5)
realband = iall & (np.abs(dall) < 1.5)
foot = (pdec < 15.6) & (pdec > 5.6) & (pra > 110.6) & (pra < 122.4) & np.isfinite(val)
cen_ok = np.where(foot)[0]

NNULL = 500
A_null, P_null, META = [], [], []
tries = 0
while len(A_null) < NNULL and tries < 40 * NNULL:
    tries += 1
    j = cen_ok[rng.integers(len(cen_ok))]
    rot = rng.uniform(0, 360); mir = bool(rng.integers(2))
    t = Trace(r_p, d_p, rot=rot, mirror=mir, centre=(pra[j], pdec[j]))
    pix = np.where((np.abs(pra - pra[j]) * np.cos(np.radians(pdec[j])) < 4.0) & (np.abs(pdec - pdec[j]) < 4.0))[0]
    dd_, ss_, ii_ = t.coords(pra[pix], pdec[pix], dmax=2.0)
    band = ii_ & (np.abs(dd_) < 1.5)
    if band.sum() < 100 or realband[pix][band].mean() > 0.05:
        continue
    fin = np.isfinite(val[pix][band]).mean()
    if fin < 0.85 or not foot[pix][band & np.isfinite(val[pix])].all():
        continue
    pp = profile(t, pix)
    a = scan(pp)
    if np.isfinite(a).all():
        A_null.append(a); P_null.append(pp); META.append((pra[j], pdec[j], rot, mir))
A_null = np.array(A_null); P_null = np.array(P_null)
print(f'nulls: {len(A_null)} accepted of {tries} tries, {time.time() - t0:.0f} s', flush=True)

med = np.median(A_null, axis=0)
sr = 0.5 * (np.percentile(A_null, 84, axis=0) - np.percentile(A_null, 16, axis=0))
z = (A_real - med) / sr
UL = np.maximum(A_real - med, 0) + 2 * sr
# global: largest z over the d0 scan (per width), calibrated on the nulls' own scans
zn = (A_null - med) / sr
gl = {}
for i, fw in enumerate(FW):
    zmax_real = float(np.nanmax(z[i])); d0_at = float(D0[np.nanargmax(z[i])])
    zmax_null = np.nanmax(zn[:, i, :], axis=1)
    gl[f'fwhm_{fw}'] = {'zmax_real': zmax_real, 'd0_at_zmax': d0_at,
                        'p_global': float((zmax_null >= zmax_real).mean()),
                        'UL2_peak_R_max_over_d0': float(np.nanmax(UL[i])),
                        'UL2_peak_R_at_d0_0': float(UL[i][np.argmin(np.abs(D0))]),
                        'UL2_peak_R_max_abs_d0_le_0.3': float(np.nanmax(UL[i][np.abs(D0) <= 0.3 + 1e-9])),
                        'UL2_peak_R_max_abs_d0_le_0.6': float(np.nanmax(UL[i][np.abs(D0) <= 0.6 + 1e-9]))}
res = {'data': 'SkyView SHASSA Sm (dR/10 -> R); spine piece Dec 12-15.6', 'n_null': int(len(A_null)), 'tries': tries,
       'd0': D0.tolist(), 'fwhm': list(FW), 'A_real': A_real.tolist(), 'null_median': med.tolist(),
       'null_sig_robust': sr.tolist(), 'z': z.tolist(), 'UL2_peak_R': UL.tolist(), 'global': gl,
       'note': ('Peak limits apply to H-alpha + (partially transmitted) [N II] 6548/6583 in the 32 A SHASSA filter, '
                'hence are conservative for H-alpha alone. FWHM 0.284 = 0.276 filament (+) 4 arcmin beam.')}
(DER / 'verify_hiwham_shassa_offset.json').write_text(json.dumps(res, indent=1))
np.savez_compressed(DER / 'verify_hiwham_shassa_offset.npz', cen=CEN, d0=D0, fwhm=np.array(FW), p_real=p_real,
                    p_null=P_null, A_real=A_real, A_null=A_null, meta=np.array(META, dtype=float))
for i, fw in enumerate(FW):
    print(f'FWHM {fw}: ' + ' '.join(f'{d0:+.2f}:{A_real[i, j]:+.2f}/{UL[i, j]:.2f}' for j, d0 in enumerate(D0) if j % 2 == 0))
print(json.dumps(gl, indent=1))
print(f'done {time.time() - t0:.0f} s')
