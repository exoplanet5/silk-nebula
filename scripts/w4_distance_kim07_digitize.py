"""w4_distance: digitise Kim et al. 2007 (ApJ 665, L139) Fig. 1, the FIMS/SPEAR C IV 1548/1550 map of the
Monogem ring, and register it to Galactic coordinates so the Silk spine can be overlaid and the C IV map sampled
across the Silk.

Input (downloaded 2026-11-17, logged in data/raw/w4_distance_kim07/fetch_log.tsv):
  data/raw/w4_distance_kim07/kim07_fg1.jpg  (600x611 px journal figure, IOP CDN)
Method
  1. Colour bar (y = 576-590 px) -> linear LUT 0 ... 21091 LU (the bar's printed end labels; ticks 4218 apart).
  2. Dotted white graticule: meridians l = 220..185 (step 5) and parallels b = 25..0 (step 5) are traced by
     iterative robust quadratic fits to white pixels; their 42 intersections fix a 2-D cubic map (x,y)<->(l,b).
  3. Each map pixel -> nearest LUT colour (RGB distance); pixels far from the LUT (graticule, labels, ellipses) or
     with value < 1500 LU (black star-mask diamonds, off-map) are invalid.
  4. The Silk spine (data/derived/silk_trace.json) is projected; C IV is sampled on lines perpendicular to the
     spine (offsets -5..+5 deg) and averaged in the south (Dec < 17.6) and north (Dec >= 17.6) sections.
  5. Cross profiles are fitted with background (linear) + Gaussian of fixed FWHM = effective map beam
     (120' smoothing (+) 55' HEALPix pixel -> ~2.1 deg; 1.9-2.4 tried) -> amplitude A and centre offset.
     Line-integrated C IV across the filament: S_CIV = A * sqrt(2 pi) sigma_b  [LU deg].
Outputs: data/derived/w4_distance_kim07_digitize.json, data/derived/w4_distance_kim07_profiles.csv
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_kim07_digitize.py
"""
import json, sys
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.optimize import curve_fit
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
IMG = ROOT / 'data/raw/w4_distance_kim07/kim07_fg1.jpg'
OUTJ = ROOT / 'data/derived/w4_distance_kim07_digitize.json'
OUTC = ROOT / 'data/derived/w4_distance_kim07_profiles.csv'

a = np.asarray(Image.open(IMG).convert('RGB')).astype(float)
H, W, _ = a.shape

# ---------------- 1. colour bar ----------------
row = a[576:591].mean(0)                       # (W,3)
nonwhite = np.where(row.min(1) < 235)[0]
x0, x1 = nonwhite.min(), nonwhite.max()
lut_rgb = a[576:591, x0:x1 + 1].mean(0)        # (n,3)
lut_val = np.linspace(0, 21091, len(lut_rgb))

# ---------------- 2. graticule ----------------
white = a[:555].min(2) > 175
ys, xs = np.nonzero(white)


def trace(seed, axis, half=3.5, n_iter=8, deg=2):
    """axis='v' -> x(y) meridian; 'h' -> y(x) parallel. seed = constant initial position."""
    coef = np.array([seed]) if True else None
    for it in range(n_iter):
        if axis == 'v':
            pred = np.polyval(coef, ys)
            m = np.abs(xs - pred) < half
            if m.sum() < 20:
                return None
            c = np.polyfit(ys[m], xs[m], min(deg, it))
            res = xs[m] - np.polyval(c, ys[m])
        else:
            pred = np.polyval(coef, xs)
            m = np.abs(ys - pred) < half
            if m.sum() < 20:
                return None
            c = np.polyfit(xs[m], ys[m], min(deg, it))
            res = ys[m] - np.polyval(c, xs[m])
        coef = c
        half = max(2.0, min(half, 3 * np.std(res) + 1))
    return coef


MER_SEED = {220: 40, 215: 117, 210: 193, 205: 268, 200: 342, 195: 417, 190: 492}
PAR_SEED = {25: 70, 20: 155, 15: 236, 10: 318, 5: 398, 0: 479}
mer = {l: trace(x, 'v') for l, x in MER_SEED.items()}
par = {b: trace(y, 'h') for b, y in PAR_SEED.items()}

pts = []   # (x, y, l, b)
for l, cm in mer.items():
    for b, cp in par.items():
        if cm is None or cp is None:
            continue
        # solve x = cm(y), y = cp(x) by fixed-point iteration
        y = PAR_SEED[b]
        for _ in range(50):
            x = np.polyval(cm, y)
            y = np.polyval(cp, x)
        if 0 <= x < W and 0 <= y < 555:
            pts.append((x, y, l, b))
pts = np.array(pts)


def design(u_, v_, deg=3):
    cols = [u_ ** i * v_ ** j for i in range(deg + 1) for j in range(deg + 1 - i)]
    return np.vstack(cols).T


xs_n, ys_n = (pts[:, 0] - 300) / 300, (pts[:, 1] - 280) / 300
ls_n, bs_n = (pts[:, 2] - 200) / 20, (pts[:, 3] - 12) / 20
A_xy = design(xs_n, ys_n)
cl = np.linalg.lstsq(A_xy, pts[:, 2], rcond=None)[0]
cb = np.linalg.lstsq(A_xy, pts[:, 3], rcond=None)[0]
A_lb = design(ls_n, bs_n)
cx = np.linalg.lstsq(A_lb, pts[:, 0], rcond=None)[0]
cy = np.linalg.lstsq(A_lb, pts[:, 1], rcond=None)[0]
res_l = pts[:, 2] - A_xy @ cl
res_b = pts[:, 3] - A_xy @ cb
res_x = pts[:, 0] - A_lb @ cx
res_y = pts[:, 1] - A_lb @ cy


def xy2lb(x, y):
    D = design((np.atleast_1d(x) - 300) / 300, (np.atleast_1d(y) - 280) / 300)
    return D @ cl, D @ cb


def lb2xy(l, b):
    D = design((np.atleast_1d(l) - 200) / 20, (np.atleast_1d(b) - 12) / 20)
    return D @ cx, D @ cy


# ---------------- 3. colour -> value ----------------
img = a[:555]
flat = img.reshape(-1, 3)
d2 = ((flat[:, None, :] - lut_rgb[None, ::2, :]) ** 2).sum(2)     # subsample LUT for speed
k = d2.argmin(1)
val = lut_val[::2][k].reshape(img.shape[:2])
dist = np.sqrt(d2.min(1)).reshape(img.shape[:2])
valid = (dist < 45) & (val > 1500) & ~(img.min(2) > 175)

# ---------------- 4. Silk spine and cross-sections ----------------
tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
sp = tr['spine']
sra = np.array([q['ra'] for q in sp]); sde = np.array([q['dec'] for q in sp])
sc = SkyCoord(sra * u.deg, sde * u.deg).galactic
sl, sb = sc.l.deg, sc.b.deg
wst = tr['strands']['W']
wra = np.array([q['ra'] for q in wst['points']]) if isinstance(wst, dict) and 'points' in wst else None


def sample(l, b, r=4.0):
    x, y = lb2xy(l, b)
    x, y = float(x[0]), float(y[0])
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    m = ((xx - x) ** 2 + (yy - y) ** 2 < r * r) & valid
    if m.sum() < 5:
        return np.nan, int(m.sum())
    return float(np.median(val[m])), int(m.sum())


offs = np.arange(-5.0, 5.01, 0.25)
rows = []
for i in range(len(sl)):
    j0, j1 = max(i - 1, 0), min(i + 1, len(sl) - 1)
    dl = (sl[j1] - sl[j0]) * np.cos(np.radians(sb[i])); db = sb[j1] - sb[j0]
    t = np.array([dl, db]) / np.hypot(dl, db)
    nvec = np.array([t[1], -t[0]])        # perpendicular; sign fixed below so + = toward higher l*cos b (east-ish)
    # East (increasing RA) at this position in Galactic frame:
    e = SkyCoord((sra[i] + 0.1) * u.deg, sde[i] * u.deg).galactic
    ev = np.array([(e.l.deg - sl[i]) * np.cos(np.radians(sb[i])), e.b.deg - sb[i]])
    if np.dot(nvec, ev) < 0:
        nvec = -nvec                       # + offset = east of the spine (convex / away-from-Monogem side)
    for o in offs:
        l = sl[i] + o * nvec[0] / np.cos(np.radians(sb[i])); b = sb[i] + o * nvec[1]
        v, n = sample(l, b)
        rows.append((i, sra[i], sde[i], sl[i], sb[i], o, l, b, v, n))
rows = np.array(rows)

import csv
with open(OUTC, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['ispine', 'ra', 'dec', 'l_spine', 'b_spine', 'offset_deg_east', 'l', 'b', 'civ_LU', 'npix'])
    for r in rows:
        w.writerow([int(r[0])] + [f'{x:.4f}' for x in r[1:8]] + [f'{r[8]:.0f}' if np.isfinite(r[8]) else 'nan', int(r[9])])


def gfun(x, A, mu, c0, c1, fwhm):
    s = fwhm / 2.3548
    return c0 + c1 * x + A * np.exp(-0.5 * ((x - mu) / s) ** 2)


out = {'input': str(IMG.relative_to(ROOT)), 'colourbar_px': [int(x0), int(x1)],
       'graticule_points': int(len(pts)),
       'registration_rms_deg': {'l': float(np.std(res_l)), 'b': float(np.std(res_b))},
       'registration_rms_px': {'x': float(np.std(res_x)), 'y': float(np.std(res_y))},
       'deg_per_px_approx': float(5 / np.median(np.diff(sorted(MER_SEED.values())))),
       'sections': {}}
for name, sel in (('south_dec_lt_17.6', sde < 17.6), ('north_dec_ge_17.6', (sde >= 17.6) & (sde <= 27.0)),
                  ('south_core_dec_12.5_16.5', (sde >= 12.5) & (sde <= 16.5))):
    idx = np.where(sel)[0]
    m = np.isin(rows[:, 0], idx)
    prof, err, nn = [], [], []
    for o in offs:
        v = rows[m & (np.abs(rows[:, 5] - o) < 1e-6), 8]
        v = v[np.isfinite(v)]
        prof.append(np.mean(v) if len(v) else np.nan)
        err.append(np.std(v) / np.sqrt(max(len(v) / 4.0, 1)) if len(v) > 1 else np.nan)   # spine pts correlated over ~1 deg (4 pts)
        nn.append(len(v))
    prof, err = np.array(prof), np.array(err)
    ok = np.isfinite(prof) & np.isfinite(err) & (err > 0)
    fits = {}
    for fw in (1.9, 2.1, 2.4):
        try:
            p, cov = curve_fit(lambda x, A, mu, c0, c1: gfun(x, A, mu, c0, c1, fw), offs[ok], prof[ok],
                               p0=[3000, 0, 6000, 0], sigma=err[ok], absolute_sigma=False, maxfev=20000)
            pe_ = np.sqrt(np.diag(cov))
            sb_ = fw / 2.3548
            fits[f'fwhm_{fw}'] = {'A_LU': float(p[0]), 'e_A': float(pe_[0]), 'mu_deg': float(p[1]), 'e_mu': float(pe_[1]),
                                 'bkg_at0_LU': float(p[2]), 'slope_LU_per_deg': float(p[3]),
                                 'S_line_integral_LU_deg': float(p[0] * np.sqrt(2 * np.pi) * sb_),
                                 'e_S': float(pe_[0] * np.sqrt(2 * np.pi) * sb_)}
        except Exception as e:
            fits[f'fwhm_{fw}'] = {'error': repr(e)}
    # also free-FWHM fit
    try:
        p, cov = curve_fit(gfun, offs[ok], prof[ok], p0=[3000, 0, 6000, 0, 2.1], sigma=err[ok], maxfev=20000)
        fits['free_fwhm'] = {'A_LU': float(p[0]), 'mu_deg': float(p[1]), 'fwhm_deg': float(p[4]),
                             'e': [float(x) for x in np.sqrt(np.diag(cov))]}
    except Exception as e:
        fits['free_fwhm'] = {'error': repr(e)}
    out['sections'][name] = {'n_spine_points': int(len(idx)), 'offsets_deg': offs.tolist(),
                             'mean_civ_LU': [None if not np.isfinite(x) else round(float(x), 1) for x in prof],
                             'err_LU': [None if not np.isfinite(x) else round(float(x), 1) for x in err],
                             'n': nn, 'fits': fits,
                             'spine_value_LU': float(prof[np.argmin(np.abs(offs))]),
                             'spine_l_range': [float(sl[idx].min()), float(sl[idx].max())],
                             'spine_b_range': [float(sb[idx].min()), float(sb[idx].max())]}

# R1 ellipse: white pixels of the ellipse near R1 label region -> l,b of its outline (for the record)
yy, xx = np.nonzero(a[:300, 180:360].min(2) > 200)
xx = xx + 180
sel = (yy > 120) & (yy < 215)
le, be = xy2lb(xx[sel], yy[sel])
out['R1_outline_white_px_lb'] = {'l_min': float(np.percentile(le, 2)), 'l_max': float(np.percentile(le, 98)),
                                 'b_min': float(np.percentile(be, 2)), 'b_max': float(np.percentile(be, 98)),
                                 'note': 'white pixels in the R1 box include graticule dots; ranges are indicative only'}
pg = SkyCoord(115.2847 * u.deg, 20.8334 * u.deg).galactic
out['pgcc_value_LU'] = sample(pg.l.deg, pg.b.deg)[0]
np.savez(ROOT / 'data/derived/w4_distance_kim07_map.npz', val=val, valid=valid, cl=cl, cb=cb, cx=cx, cy=cy, pts=pts)
json.dump(out, open(OUTJ, 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'sections'}, indent=1))
for k, v in out['sections'].items():
    print(k, v['n_spine_points'], 'spine', round(v['spine_value_LU']), {kk: {q: round(w, 3) for q, w in vv.items() if isinstance(w, float)} for kk, vv in v['fits'].items()})
