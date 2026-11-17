"""verify_hiwham: significance checks for verify_hiwham_hi_extra.py.

1. PGCC G199.09+19.95 narrow HI line: the same (r < 0.1 deg) minus (0.3-0.8 deg annulus) difference spectrum at
   2000 random GALFA positions; distribution of the peak dT (V_LSR -25..+35) and of its velocity; fraction of
   random positions with a narrow (FWHM < 5 km/s) excess >= the PGCC value; fraction whose peak velocity lies
   within +-1.5 km/s of the W-strand [O III] V_LSR (+6.6).
2. GS202+17+019: elliptical radius rho of the Silk spine vs Dec for two orientation conventions (phi from the
   l axis, and from the b axis).
Output: data/derived/verify_hiwham_hi_extra2.json
"""
import json
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import DER, strand, Trace, galfa_cube, vlsr_minus_vhelio

rng = np.random.default_rng(11)
cube, v, gra, gdec = galfa_cube()
RA, DEC = np.meshgrid(gra, gdec)
pra, pdec = RA.ravel(), DEC.ravel()
C2 = cube.reshape(len(v), -1)
cover = np.isfinite(C2[40])
# unit vectors for fast neighbour search
def uvec(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]


X = uvec(pra[cover], pdec[cover]); idx_cov = np.where(cover)[0]
tree = cKDTree(X)
chord = lambda deg: 2 * np.sin(np.radians(deg) / 2)
k = (v > -25) & (v < 35)


def diffspec(ra0, dec0):
    x0 = uvec(np.array([ra0]), np.array([dec0]))[0]
    on = idx_cov[tree.query_ball_point(x0, chord(0.1))]
    ann_all = idx_cov[tree.query_ball_point(x0, chord(0.8))]
    inner = set(idx_cov[tree.query_ball_point(x0, chord(0.3))])
    ann = np.array([i for i in ann_all if i not in inner])
    if len(on) < 10 or len(ann) < 300:
        return None
    return np.nanmean(C2[:, on], axis=1) - np.nanmean(C2[:, ann], axis=1)


def narrow_peak(dif):
    """peak of a 3-channel (2.2 km/s) smoothed difference and a crude width (channels above half max)."""
    sm = np.convolve(dif, np.ones(3) / 3, mode='same')
    j = np.argmax(np.where(k, sm, -99))
    half = sm[j] / 2
    lo = j
    while lo > 0 and sm[lo] > half:
        lo -= 1
    hi = j
    while hi < len(sm) - 1 and sm[hi] > half:
        hi += 1
    return float(sm[j]), float(v[j]), float((hi - lo) * abs(v[1] - v[0]))


pg = SkyCoord(l=199.09 * u.deg, b=19.95 * u.deg, frame='galactic').icrs
d0 = diffspec(pg.ra.deg, pg.dec.deg)
pk0 = narrow_peak(d0)
ok = np.where(cover & (np.abs(pdec - 18) < 11))[0]
vals = []
while len(vals) < 2000:
    j = ok[rng.integers(len(ok))]
    if SkyCoord(pra[j] * u.deg, pdec[j] * u.deg).separation(pg).deg < 1.0:
        continue
    dd = diffspec(pra[j], pdec[j])
    if dd is None:
        continue
    vals.append(narrow_peak(dd))
vals = np.array(vals)
vW = 18.5 + vlsr_minus_vhelio(115.9, 19.6)[0]
res = {'PGCC': {'peak_dT_K': pk0[0], 'v_peak': pk0[1], 'width_kms_halfmax': pk0[2],
                'random_frac_peak_ge': float(np.mean(vals[:, 0] >= pk0[0])),
                'random_frac_narrow_peak_ge': float(np.mean((vals[:, 0] >= pk0[0]) & (vals[:, 2] <= 5.0))),
                'random_peak_p50_p84_p97.5': [float(np.percentile(vals[:, 0], q)) for q in (50, 84, 97.5)],
                'random_frac_vpeak_within_1.5_of_W': float(np.mean(np.abs(vals[:, 1] - vW) <= 1.5)),
                'random_frac_vpeak_within_1.5_of_W_given_peak_gt_2K': float(np.mean(np.abs(vals[vals[:, 0] > 2, 1] - vW) <= 1.5)),
                'W_strand_vlsr': float(vW), 'n_random': int(len(vals))}}
print(json.dumps(res, indent=1))

# ---------------------------------------------------------------- GS202+17+019 rho along spine
ra_E, dec_E = strand('E')
sp = SkyCoord(ra_E * u.deg, dec_E * u.deg).galactic
sh = dict(l=202.4, b=16.9, a=5.7, bb=4.5, phi=29.6)
out = {}
for conv in ('phi_from_l', 'phi_from_b'):
    x = ((sp.l.deg - sh['l'] + 180) % 360 - 180) * np.cos(np.radians(sh['b'])); y = sp.b.deg - sh['b']
    ph = np.radians(sh['phi'] if conv == 'phi_from_l' else sh['phi'] + 90)
    xr = x * np.cos(ph) + y * np.sin(ph); yr = -x * np.sin(ph) + y * np.cos(ph)
    rho = np.sqrt((xr / sh['a']) ** 2 + (yr / sh['bb']) ** 2)
    out[conv] = {'dec': dec_E[::5].round(2).tolist(), 'rho': rho[::5].round(2).tolist(),
                 'rho_median': float(np.median(rho)), 'frac_0.8_1.2': float(np.mean((rho > 0.8) & (rho < 1.2)))}
    # distance from the Suad HI minimum (203.3, 19.1)
mn = SkyCoord(l=203.3 * u.deg, b=19.1 * u.deg, frame='galactic')
out['sep_spine_to_HI_minimum_deg'] = {'min': float(sp.separation(mn).deg.min()), 'at_dec': float(dec_E[np.argmin(sp.separation(mn).deg)])}
res['GS202+17+019_spine_rho'] = out
print(json.dumps(out, indent=0)[:1500])
(DER / 'verify_hiwham_hi_extra2.json').write_text(json.dumps(res, indent=1))
