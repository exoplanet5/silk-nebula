"""w4_distance: placebo test for the FIMS C IV / south-Silk coincidence on the digitised Kim+07 Fig. 1.

Statistic: contrast = mean C IV along the south spine track (Dec < 17.6, 23 points) minus the mean of the two
parallel tracks offset by +-3.0 deg (perpendicular to the local track). Placebos: the same track rigidly rotated
(random angle) and translated (random centre) anywhere the whole track and its +-3 deg side tracks are on valid
map pixels; (a) all of the map, (b) only placements with all points at b >= 12 deg (the high-latitude part, where
R1 is the only C IV region Kim+07 marked).
Output: data/derived/w4_distance_kim07_placebo.json
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_kim07_placebo.py
"""
import json
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
M = np.load(ROOT / 'data/derived/w4_distance_kim07_map.npz')
val, valid, cx, cy = M['val'], M['valid'], M['cx'], M['cy']


def design(u_, v_, deg=3):
    return np.vstack([u_ ** i * v_ ** j for i in range(deg + 1) for j in range(deg + 1 - i)]).T


def lb2xy(l, b):
    D = design((np.atleast_1d(l) - 200) / 20, (np.atleast_1d(b) - 12) / 20)
    return D @ cx, D @ cy


# local median map (radius 4 px) for fast sampling
from scipy.ndimage import generic_filter
vv = np.where(valid, val, np.nan)
yy, xx = np.mgrid[-4:5, -4:5]
fp = (xx ** 2 + yy ** 2) < 16
med = generic_filter(vv, np.nanmedian, footprint=fp, mode='constant', cval=np.nan)


def samp(l, b):
    x, y = lb2xy(l, b)
    xi, yi = np.round(x).astype(int), np.round(y).astype(int)
    ok = (xi >= 0) & (xi < med.shape[1]) & (yi >= 0) & (yi < med.shape[0])
    out = np.full(len(l), np.nan)
    out[ok] = med[yi[ok], xi[ok]]
    return out


tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
sp = [p for p in tr['spine'] if p['dec'] < 17.6]
c = SkyCoord([p['ra'] for p in sp] * u.deg, [p['dec'] for p in sp] * u.deg).galactic
l0, b0 = c.l.deg, c.b.deg
# local tangent-plane coordinates about the track centre
lc, bc = l0.mean(), b0.mean()
X = (l0 - lc) * np.cos(np.radians(bc)); Y = b0 - bc
t = np.gradient(np.c_[X, Y], axis=0); t /= np.linalg.norm(t, axis=1)[:, None]
N = np.c_[t[:, 1], -t[:, 0]]


def contrast(Xc, Yc, Nc, lc_, bc_):
    def tolb(XX, YY):
        return lc_ + XX / np.cos(np.radians(bc_)), bc_ + YY
    on = samp(*tolb(Xc, Yc))
    s1 = samp(*tolb(Xc + 3 * Nc[:, 0], Yc + 3 * Nc[:, 1]))
    s2 = samp(*tolb(Xc - 3 * Nc[:, 0], Yc - 3 * Nc[:, 1]))
    if np.isnan(on).mean() > 0.2 or np.isnan(s1).mean() > 0.3 or np.isnan(s2).mean() > 0.3:
        return np.nan, None
    bmin = min(tolb(Xc, Yc)[1].min(), tolb(Xc + 3 * Nc[:, 0], Yc + 3 * Nc[:, 1])[1].min(), tolb(Xc - 3 * Nc[:, 0], Yc - 3 * Nc[:, 1])[1].min())
    return np.nanmean(on) - 0.5 * (np.nanmean(s1) + np.nanmean(s2)), bmin


obs, _ = contrast(X, Y, N, lc, bc)
rng = np.random.default_rng(42)
res = []
for k in range(6000):
    th = rng.uniform(0, 2 * np.pi)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    XY = np.c_[X, Y] @ R.T; NN = N @ R.T
    lc_ = rng.uniform(185, 220); bc_ = rng.uniform(-2, 24)
    cval, bmin = contrast(XY[:, 0], XY[:, 1], NN, lc_, bc_)
    if np.isfinite(cval):
        res.append((cval, bmin, lc_, bc_))
res = np.array(res)
allc = res[:, 0]; hi = res[res[:, 1] >= 12, 0]
out = {'observed_contrast_LU': float(obs), 'n_placebo_valid': int(len(allc)), 'n_placebo_b_ge_12': int(len(hi)),
       'p_all': float((allc >= obs).mean()), 'p_b_ge_12': float((hi >= obs).mean()) if len(hi) else None,
       'placebo_all_pct': {q: float(np.percentile(allc, q)) for q in (50, 84, 97.5, 99.7)},
       'placebo_hi_pct': {q: float(np.percentile(hi, q)) for q in (50, 84, 97.5)} if len(hi) else None,
       'note': 'R1 was defined by Kim+07 from the same map, so this tests how unusual a 5.5-deg track this bright is, '
               'not an a-priori detection; a-posteriori choice of the south section (Dec < 17.6) also applies.'}
json.dump(out, open(ROOT / 'data/derived/w4_distance_kim07_placebo.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
