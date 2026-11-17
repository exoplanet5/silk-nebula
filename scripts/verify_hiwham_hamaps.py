"""verify_hiwham: SHASSA and Finkbeiner (2003) H-alpha on/off limits for the Silk spine, with placement nulls.

* SHASSA (SkyView 'SHASSA Sm', smoothed continuum-corrected H-alpha in dR, ~4' resolution; data only to
  Dec +15.98 here) and 'SHASSA CC' (unsmoothed, star residuals -> binned MEDIANS).
* Finkbeiner 2003 composite (HEALPix Nside 1024, 6' FWHM). Its survey bitmask shows that along the Silk spine
  only Dec <= ~15.75 uses SHASSA; north of that it is WHAM-only (1-deg resolution). Both parts are analysed
  separately, and the wave-1 assumption (6' resolution everywhere) is checked.
Statistic: amplitude A of a Gaussian filament (FWHM_obs = sqrt(0.276^2 + beam^2)) at d=0 over a linear baseline
fitted in |d| < 1.2 deg, on 0.05-deg binned profiles; intrinsic peak = A * FWHM_obs/0.276.
Nulls: the same spine piece placed at random centres/orientations (0-360 deg, random E-W mirror) within the map
footprint, avoiding the real |d| < 1.5 deg band; sigma = std of null amplitudes.
Output: data/derived/verify_hiwham_hamaps_results.json, verify_hiwham_hamaps_profiles.npz
"""
import json
import sys
import warnings

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import RAW, DER, strand, Trace

warnings.filterwarnings('ignore')
rng = np.random.default_rng(42)
ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
EDG = np.arange(-2.0, 2.0 + 1e-9, 0.05); CEN = 0.5 * (EDG[1:] + EDG[:-1])


def piece(dlo, dhi):
    m = (REAL.dec >= dlo) & (REAL.dec <= dhi)
    return REAL.ra[m], REAL.dec[m]


def prof(tr, pra, pdec, val, stat='mean', pix=None):
    if pix is None:
        pix = np.arange(len(pra))
    d, s, ins = tr.coords(pra[pix], pdec[pix], dmax=2.0)
    m = ins & np.isfinite(d) & np.isfinite(val[pix]) & (np.abs(d) < 2.0)
    b = ((d[m] + 2.0) / 0.05).astype(int)
    x = val[pix][m]
    nb = len(CEN)
    n = np.bincount(b, minlength=nb)
    if stat == 'mean':
        p = np.bincount(b, weights=x, minlength=nb) / np.maximum(n, 1)
    else:
        o = np.argsort(b, kind='stable'); bs, xs = b[o], x[o]
        bd = np.searchsorted(bs, np.arange(nb + 1))
        p = np.array([np.median(xs[bd[i]:bd[i + 1]]) if n[i] > 0 else np.nan for i in range(nb)])
    p[n < 5] = np.nan
    return p, n


def amp(p, fwhm_obs):
    w = (np.abs(CEN) < 1.2) & np.isfinite(p)
    sg = fwhm_obs / 2.3548
    X = np.c_[np.ones(w.sum()), CEN[w], np.exp(-0.5 * (CEN[w] / sg) ** 2)]
    if w.sum() < 20:
        return np.nan, np.nan
    c, *_ = np.linalg.lstsq(X, p[w], rcond=None)
    core = np.abs(CEN) < 0.15; ring = (np.abs(CEN) > 0.35) & (np.abs(CEN) < 0.8)
    return c[2], np.nanmean(p[core]) - np.nanmean(p[ring])


def run(name, pra, pdec, val, beam_deg, dec_rng, foot, stat='mean', nnull=400, minfin=0.85):
    """foot(ra, dec) -> bool array of valid footprint for null centres."""
    r, d = piece(*dec_rng)
    tr = Trace(r, d)
    fw = np.hypot(0.276, beam_deg)
    p, n = prof(tr, pra, pdec, val, stat)
    A, C = amp(p, fw)
    dall, sall, iall = REAL.coords(pra, pdec, dmax=2.5)
    realband = iall & (np.abs(dall) < 1.5)
    cen_ok = np.where(foot(pra, pdec) & np.isfinite(val))[0]
    NA, NC, PN = [], [], []
    tries = 0
    while len(NA) < nnull and tries < 30 * nnull:
        tries += 1
        j = cen_ok[rng.integers(len(cen_ok))]
        t = Trace(r, d, rot=rng.uniform(0, 360), mirror=bool(rng.integers(2)), centre=(pra[j], pdec[j]))
        near = (np.abs(pra - pra[j]) * np.cos(np.radians(pdec[j])) < 4.5) & (np.abs(pdec - pdec[j]) < 4.5)
        pix = np.where(near)[0]
        dd, ss, ii = t.coords(pra[pix], pdec[pix], dmax=2.0)
        band = ii & (np.abs(dd) < 1.5)
        if band.sum() < 100 or realband[pix][band].mean() > 0.05:
            continue
        if np.isfinite(val[pix][band]).mean() < minfin or not foot(pra[pix][band], pdec[pix][band]).all():
            continue
        pp, nn = prof(t, pra, pdec, val, stat, pix=pix)
        a, cc = amp(pp, fw)
        if np.isfinite(a):
            NA.append(a); NC.append(cc); PN.append(pp)
    NA, NC = np.array(NA), np.array(NC)
    corr_peak = fw / 0.276
    out = {'name': name, 'dec_range': dec_rng, 'beam_deg': beam_deg, 'stat': stat, 'fwhm_obs_deg': fw,
           'A': float(A), 'C_core': float(C), 'null_n': int(len(NA)), 'null_mean_A': float(NA.mean()),
           'null_std_A': float(NA.std()), 'null_std_C': float(NC.std()), 'null_mean_C': float(NC.mean()),
           'z_A': float((A - NA.mean()) / NA.std()),
           'peak_intrinsic': float((A - NA.mean()) * corr_peak),
           'peak_intrinsic_2sigma_UL': float((max(A - NA.mean(), 0) + 2 * NA.std()) * corr_peak),
           'peak_conversion_factor': corr_peak}
    med = float(np.median(NA)); sr = float(0.5 * (np.percentile(NA, 84) - np.percentile(NA, 16)))
    out.update({'null_median_A': med, 'null_sig_robust_A': sr, 'z_A_robust': float((A - med) / sr),
                'peak_intrinsic_robust': float((A - med) * corr_peak),
                'peak_intrinsic_2sigma_UL_robust': float((max(A - med, 0) + 2 * sr) * corr_peak),
                'null_A_values': [round(float(x), 4) for x in NA]})
    print(f"{name:28s} Dec {dec_rng}: A={A:+.3f} null {NA.mean():+.3f}+-{NA.std():.3f} (n={len(NA)}) "
          f"z={out['z_A']:+.2f}  intrinsic peak 2sig UL {out['peak_intrinsic_2sigma_UL']:.3f} R  (tries {tries}) | robust: "
          f"med {med:+.3f} sig_r {sr:.3f} z {out['z_A_robust']:+.2f} UL {out['peak_intrinsic_2sigma_UL_robust']:.3f} R", flush=True)
    return out, p, np.array(PN)


PARTS = sys.argv[1].split(',') if len(sys.argv) > 1 else ['Sm', 'CC', 'F03']
results = {}
profs = {}
# ------------------------------------------------------------------ SHASSA
for prod, stat in (('Sm', 'mean'), ('CC', 'median')):
    if prod not in PARTS:
        continue
    f = RAW / f'skyview/SHASSA_{prod}/SHASSA_{prod}_wide.fits'
    img = fits.getdata(f).astype(float) / 10.0      # dR -> R
    w = WCS(fits.getheader(f))
    ny, nx = img.shape
    if prod == 'CC':      # 2x2 block to speed up (0.8' -> 1.6' pixels), median later
        pass
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr, dd = w.all_pix2world(xx, yy, 0)
    if prod == 'CC':      # decimate 2x2 (0.8' -> 1.6' sampling) to keep the binned medians affordable
        img, rr, dd = img[::2, ::2], rr[::2, ::2], dd[::2, ::2]
    v = img.ravel().copy(); pr = rr.ravel(); pd_ = dd.ravel()
    v[(v == 0)] = np.nan
    foot = lambda a, b: (b < 15.6) & (b > 5.6) & (a > 110.6) & (a < 122.4)
    o, p, pn = run(f'SHASSA_{prod}', pr, pd_, v, 4.0 / 60, (12.0, 15.6), foot, stat=stat,
                   nnull=400 if prod == 'Sm' else 200)
    results[f'SHASSA_{prod}'] = o; profs[f'SHASSA_{prod}'] = p; profs[f'SHASSA_{prod}_nulls'] = pn

# ------------------------------------------------------------------ Finkbeiner 2003
if 'F03' not in PARTS:
    tag = '_'.join(PARTS)
    (DER / f'verify_hiwham_hamaps_part_{tag}.json').write_text(json.dumps(results, indent=1))
    np.savez_compressed(DER / f'verify_hiwham_hamaps_part_{tag}.npz', cen=CEN, **{k: np.asarray(v) for k, v in profs.items()})
    sys.exit(0)
fk = fits.getdata(RAW / 'lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits', 1)['TEMPERATURE'].astype(float).ravel()
mk = fits.getdata(RAW / 'lambda_finkbeiner_halpha/Halpha_mask_fwhm06_1024.fits', 1)['MASK'].astype(int).ravel()
cen = SkyCoord(116.0 * u.deg, 18.0 * u.deg).galactic
vec = hp.ang2vec(cen.l.deg, cen.b.deg, lonlat=True)
ipix = hp.query_disc(1024, vec, np.radians(16.0), nest=True)
gl, gbb = hp.pix2ang(1024, ipix, nest=True, lonlat=True)
cc = SkyCoord(gl * u.deg, gbb * u.deg, frame='galactic').icrs
fra, fdec = cc.ra.deg, cc.dec.deg
fval = fk[ipix].copy()
fbits = mk[ipix]
shassa_bit = (fbits >> 2) & 1
fval[(fbits >> 3) & 1 == 1] = np.nan     # star-removed pixels
# (i) the SHASSA-based southern part
foot_s = lambda a, b: np.isin(np.ones_like(a, int), [1]) & (b < 15.5) & (b > 4) & (a > 106) & (a < 126)
fv_s = np.where(shassa_bit == 1, fval, np.nan)
print('F03 star-flag fraction in disc:', float(((fbits >> 3) & 1).mean()), flush=True)
o, p, pn = run('Finkbeiner_SHASSApart', fra, fdec, fv_s, 0.1, (12.0, 15.6), lambda a, b: (b < 15.5) & (b > 3) & (a > 106) & (a < 127), minfin=0.6)
results['Finkbeiner_SHASSApart'] = o; profs['Finkbeiner_SHASSApart'] = p
# (ii) the WHAM-only northern part, analysed (a) as if 6' (wave-1 assumption) and (b) with the true ~1 deg beam
fv_w = np.where(shassa_bit == 0, fval, np.nan)
for lab, bm in (('as6arcmin', 0.1), ('as1deg', 1.0)):
    o, p, pn = run(f'Finkbeiner_WHAMpart_{lab}', fra, fdec, fv_w, bm, (15.9, 28.25),
                   lambda a, b: (b > 16.0) & (b < 34) & (a > 100) & (a < 132), nnull=300, minfin=0.6)
    results[f'Finkbeiner_WHAMpart_{lab}'] = o; profs[f'Finkbeiner_WHAMpart_{lab}'] = p
# fraction of the spine (by length) with SHASSA data in the composite
sp_c = SkyCoord(REAL.ra * u.deg, REAL.dec * u.deg).galactic
sp_bits = mk[hp.ang2pix(1024, sp_c.l.deg, sp_c.b.deg, lonlat=True, nest=True)]
results['spine_fraction_with_SHASSA_in_F03'] = float(((sp_bits >> 2) & 1).mean())
results['spine_dec_max_SHASSA_in_F03'] = float(REAL.dec[((sp_bits >> 2) & 1) == 1].max())
results['notes'] = ('Peak values are intrinsic filament H-alpha peak (R) for a Gaussian of FWHM 0.276 deg; '
                    'SHASSA includes [N II] 6548/6583 in its 32 A filter, so its limit applies to Ha+[N II].')
tag = '_'.join(PARTS)
(DER / f'verify_hiwham_hamaps_part_{tag}.json').write_text(json.dumps(results, indent=1))
np.savez_compressed(DER / f'verify_hiwham_hamaps_part_{tag}.npz', cen=CEN, **{k: np.asarray(v) for k, v in profs.items()})
# merge all available part files into the combined results file
import glob
merged = {}
for fpart in sorted(glob.glob(str(DER / 'verify_hiwham_hamaps_part_*.json'))):
    merged.update(json.loads(open(fpart).read()))
(DER / 'verify_hiwham_hamaps_results.json').write_text(json.dumps(merged, indent=1))
print(json.dumps({k: results[k] for k in ['spine_fraction_with_SHASSA_in_F03', 'spine_dec_max_SHASSA_in_F03']}))
