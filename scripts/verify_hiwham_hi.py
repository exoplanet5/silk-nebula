"""verify_hiwham: GALFA-HI DR2 cross-spine profiles at the revised [O III] velocity, with trace-placement nulls.

Independent of wave-1 code (only the trace JSON and the raw GALFA cutouts are used).

Steps
 1. 2' GALFA mosaic (verify_hiwham_common.galfa_cube) -> N(HI) maps in
      * 20 LSR slices of 2 km/s from -15 to +25 km/s,
      * windows: heliocentric [O III] (V_helio +10.5..+16.5 = V_LSR ~ -2..+6 along the spine, per pixel),
        W strand (V_helio +16..+21), old [O III] (V_LSR +11..+16), red lines (-1..+6), all (-29.8..+40.1).
 2. Cross-spine profiles (0.05 deg bins, |d| < 3 deg) for the real E spine: all, south (spine Dec <= 17.6),
    north, and eleven 1.5-deg segments.
 3. Statistics per profile:
      A_fil   amplitude of a Gaussian HI filament at d = 0 (FWHM = sqrt(0.276^2 + 0.067^2) deg) over a linear
              baseline fitted in |d| < 1.2 deg  [cm^-2]
      C_core  mean(|d|<0.15) - mean(0.35<|d|<0.8)
      S_loc   mean(+0.2..+1.2) - mean(-1.2..-0.2)   (outside minus inside)
      S_wide  mean(+0.5..+2.5) - mean(-2.5..-0.5)
      step    best smoothed-step fit (sigma 0.1 deg) over a linear baseline, d0 scanned -1.5..+2.0 in a +-1.2 deg
              window; returns amplitude (>0 = HI rises outward) and d0
      w1wall  re-implementation of the wave-1 'HI wall' half-rise finder (0.1-deg bins, 3-bin running mean,
              floor = min in -1..+1, peak = max in 0..+3, first mid-level crossing) -> offset and wall/floor
 4. Nulls
      (a) coherent: rigid copies (translation / rotation / E-W mirror) of the whole E spine at positions whose
          |d|<1.5 band overlaps the real band by < 5 % and whose |d|<2.5 band is >= 75 % covered;
      (b) independent segments: each 1.5-deg spine piece placed at random centres and random orientation
          (0-360 deg, random mirror) inside the GALFA footprint, away from the real band.
Outputs: data/derived/verify_hiwham_hi_profiles.npz, data/derived/verify_hiwham_hi_results.json
"""
import json
import sys
import time

import numpy as np
from scipy import sparse

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import (DER, SCR, strand, Trace, galfa_cube, nhi_slice, vlsr_minus_vhelio,
                                  spine_dec_of_s)

rng = np.random.default_rng(20261116)
t0 = time.time()
cube, v, gra, gdec = galfa_cube()
RA, DEC = np.meshgrid(gra, gdec)
ny, nx = RA.shape
pra, pdec = RA.ravel(), DEC.ravel()

# ------------------------------------------------------------------ maps
slices = [(lo, lo + 2.0) for lo in np.arange(-15.0, 25.0, 2.0)]
names = [f'lsr{lo:+.0f}{hi:+.0f}' for lo, hi in slices]
maps = [nhi_slice(cube, v, lo, hi) for lo, hi in slices]
# heliocentric windows -> per-pixel LSR windows
corr = vlsr_minus_vhelio(pra, pdec).reshape(ny, nx)


def helio_window(lo_h, hi_h):
    dv = abs(v[1] - v[0])
    e0, e1 = v - dv / 2, v + dv / 2
    lo = lo_h + corr; hi = hi_h + corr
    out = np.zeros((ny, nx))
    bad = np.zeros((ny, nx), bool)
    for k in range(len(v)):
        w = np.clip(np.minimum(e1[k], hi) - np.maximum(e0[k], lo), 0, None)
        if not np.any(w > 0):
            continue
        t = cube[k]
        out += w * np.nan_to_num(t)
        bad |= (w > 0) & ~np.isfinite(t)
    out *= 1.823e18
    out[bad] = np.nan
    return out


extra = {'oiii_helio': helio_window(10.5, 16.5), 'W_helio': helio_window(16.0, 21.0),
         'old_oiii': nhi_slice(cube, v, 11.0, 16.0), 'red': nhi_slice(cube, v, -1.0, 6.0),
         'all': nhi_slice(cube, v, -29.9, 40.2)}
for k, m in extra.items():
    names.append(k); maps.append(m)
M = np.array([m.ravel() for m in maps]).T.astype(np.float64)      # (npix, nmap)
F = np.isfinite(M)
M0 = np.where(F, M, 0.0)
F = F.astype(np.float64)
cover = np.isfinite(maps[names.index('all')]).ravel()
print(f'maps ready {time.time() - t0:.1f}s; nmap={M.shape[1]}', flush=True)

# ------------------------------------------------------------------ profile machinery
DB = 0.05
EDGES = np.arange(-3.0, 3.0 + 1e-9, DB)
CEN = 0.5 * (EDGES[1:] + EDGES[:-1])
NB = len(CEN)
ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
L = REAL.L
SEG_EDGES = np.arange(0, 16.5 + 1e-9, 1.5)
NSEG = len(SEG_EDGES) - 1
s_split = float(np.interp(17.6, REAL.dec, REAL.s))
GROUPS = [f'seg{i}' for i in range(NSEG)] + ['south', 'north', 'all']


def group_index(s):
    """list of (group id array) memberships for each pixel: returns list of (gid, mask)."""
    out = []
    seg = np.floor(s / 1.5).astype(int)
    for i in range(NSEG):
        out.append((i, seg == i))
    out.append((NSEG, s <= s_split))
    out.append((NSEG + 1, s > s_split))
    out.append((NSEG + 2, np.ones_like(s, bool)))
    return out


def profiles(tr, ngroups_full=True, dmax=3.0, pix=None):
    """Mean profiles (ngroup, NB, nmap) and pixel counts (ngroup, NB) for a Trace copy."""
    if pix is None:
        pix = np.arange(len(pra))
    d, s, inside = tr.coords(pra[pix], pdec[pix], dmax=dmax)
    m = inside & np.isfinite(d) & (np.abs(d) < dmax)
    idx = pix[m]; d = d[m]; s = s[m]
    b = ((d + 3.0) / DB).astype(int)
    rows, cols = [], []
    if ngroups_full:
        for gid, gm in group_index(s):
            rows.append(gid * NB + b[gm]); cols.append(idx[gm])
        ng = len(GROUPS)
    else:
        rows.append(b); cols.append(idx); ng = 1
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(ng * NB, len(pra)))
    S = A @ M0
    N = A @ F
    with np.errstate(invalid='ignore', divide='ignore'):
        P = np.where(N >= 4, S / N, np.nan)
    return P.reshape(ng, NB, -1), N[:, -1].reshape(ng, NB), (d, idx)


# ------------------------------------------------------------------ statistics on one profile (NB,) array
FWHM_FIL = np.hypot(0.276, 4.0 / 60)
SIG_FIL = FWHM_FIL / 2.3548
core = np.abs(CEN) < 0.15
ring = (np.abs(CEN) > 0.35) & (np.abs(CEN) < 0.8)
loc_o = (CEN > 0.2) & (CEN < 1.2); loc_i = (CEN < -0.2) & (CEN > -1.2)
wid_o = (CEN > 0.5) & (CEN < 2.5); wid_i = (CEN < -0.5) & (CEN > -2.5)
D0_GRID = np.arange(-1.5, 2.0001, 0.05)
from scipy.special import erf


def lsq(X, y):
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    if ok.sum() < X.shape[1] + 5:
        return None, None, None
    c, res, rk, sv = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    r = y[ok] - X[ok] @ c
    return c, float(np.sum(r ** 2)), int(ok.sum())


def stats(p):
    out = {}
    w = np.abs(CEN) < 1.2
    X = np.c_[np.ones(w.sum()), CEN[w], np.exp(-0.5 * (CEN[w] / SIG_FIL) ** 2)]
    c, _, _ = lsq(X, p[w])
    out['A_fil'] = np.nan if c is None else c[2]
    out['C_core'] = np.nanmean(p[core]) - np.nanmean(p[ring])
    out['S_loc'] = np.nanmean(p[loc_o]) - np.nanmean(p[loc_i])
    out['S_wide'] = np.nanmean(p[wid_o]) - np.nanmean(p[wid_i])
    best = (np.nan, np.nan, -np.inf)
    for d0 in D0_GRID:
        w = np.abs(CEN - d0) < 1.2
        Xs = np.c_[np.ones(w.sum()), CEN[w] - d0]
        c0, r0, n0 = lsq(Xs, p[w])
        if c0 is None:
            continue
        H = 0.5 * (1 + erf((CEN[w] - d0) / (0.1 * np.sqrt(2))))
        c1, r1, n1 = lsq(np.c_[Xs, H], p[w])
        if c1 is None:
            continue
        gain = (r0 - r1) / max(r1 / max(n1 - 3, 1), 1e-30)
        # keep the strongest POSITIVE (outward-rising) step
        if c1[2] > 0 and gain > best[2]:
            best = (c1[2], d0, gain)
    out['step_amp'], out['step_d0'], out['step_gain'] = best
    # wave-1 wall re-implementation on 0.1-deg bins (mean of pairs of 0.05 bins), 3-bin running mean
    p2 = np.nanmean(p.reshape(-1, 2), axis=1); c2 = CEN.reshape(-1, 2).mean(axis=1)
    k = np.convolve(np.nan_to_num(p2, nan=np.nanmedian(p2)), np.ones(3) / 3, mode='same')
    k[~np.isfinite(p2)] = np.nan
    fl = (c2 >= -1.0) & (c2 <= 1.0); pk = (c2 >= 0.0) & (c2 <= 3.0)
    wall, ratio = np.nan, np.nan
    if np.isfinite(k[pk]).any() and np.isfinite(k[fl]).any():
        floor = np.nanmin(k[fl]); ifl = np.where(fl)[0][np.nanargmin(k[fl])]
        ipk = np.where(pk)[0][np.nanargmax(k[pk])]
        mid = 0.5 * (floor + k[ipk])
        for i in range(ifl, ipk + 1):
            if np.isfinite(k[i]) and k[i] > mid and i > 0:
                x0, x1, y0, y1 = c2[i - 1], c2[i], k[i - 1], k[i]
                wall = x0 + (mid - y0) * (x1 - x0) / (y1 - y0) if np.isfinite(y0) and y1 != y0 else x1
                ratio = k[ipk] / max(floor, 1e10)
                break
    out['w1_wall'], out['w1_ratio'] = wall, ratio
    return out


STATS = ['A_fil', 'C_core', 'S_loc', 'S_wide', 'step_amp', 'step_d0', 'step_gain', 'w1_wall', 'w1_ratio']


def stats_array(P):
    """P: (ngroup, NB, nmap) -> (ngroup, nmap, nstat)"""
    ng, _, nm = P.shape
    out = np.full((ng, nm, len(STATS)), np.nan)
    for g in range(ng):
        for j in range(nm):
            if np.isfinite(P[g, :, j]).sum() > 20:
                st = stats(P[g, :, j])
                out[g, j] = [st[s] for s in STATS]
    return out


# ------------------------------------------------------------------ real trace
P_real, N_real, (d_real_sel, idx_real) = profiles(REAL)
ST_real = stats_array(P_real)
print(f'real done {time.time() - t0:.1f}s', flush=True)
d_all, s_all, in_all = REAL.coords(pra, pdec, dmax=4.0)
real_band = in_all & (np.abs(d_all) < 1.5)

# ------------------------------------------------------------------ (a) coherent full-trace nulls
cands = []
for rot in [-20, -10, 0, 10, 20, 160, 170, 180, 190, 200]:
    for mirror in [False, True]:
        for rc in np.arange(109.0, 120.01, 0.5):
            for dc in np.arange(13.0, 22.01, 0.5):
                cands.append((rot, mirror, rc, dc))
rng.shuffle(cands)
acc = []
sub = rng.choice(len(pra), 40000, replace=False)
acc_cache = SCR / 'hi_coherent_accepted.json'
if acc_cache.exists():
    acc = [tuple(a) for a in json.loads(acc_cache.read_text())]
    cands = []
for rot, mirror, rc, dc in cands:
    tr = Trace(ra_E, dec_E, rot=rot, mirror=mirror, centre=(rc, dc))
    d, s, ins = tr.coords(pra[sub], pdec[sub], dmax=3.0)
    band = ins & (np.abs(d) < 2.5)
    if band.sum() < 50:
        continue
    covf = cover[sub][band].mean()
    inner = ins & (np.abs(d) < 1.5)
    ovl = real_band[sub][inner].mean() if inner.sum() else 1
    if covf >= 0.75 and ovl < 0.05:
        acc.append((rot, mirror, rc, dc, covf, ovl))
print(f'coherent candidates accepted: {len(acc)} of {len(cands)}  {time.time() - t0:.1f}s', flush=True)
if not acc_cache.exists():
    acc_cache.write_text(json.dumps([[float(a[0]), bool(a[1]), float(a[2]), float(a[3]), float(a[4]), float(a[5])] for a in acc]))
rng.shuffle(acc)
# thin: keep placements whose centres differ by >= 0.75 deg or orientation differs
kept = []
for a in acc:
    if all(not (abs(a[2] - b[2]) < 0.75 and abs(a[3] - b[3]) < 0.75 and a[1] == b[1] and abs(a[0] - b[0]) < 15) for b in kept):
        kept.append(a)
kept = kept[:120]
ST_coh = []
for (rot, mirror, rc, dc, covf, ovl) in kept:
    tr = Trace(ra_E, dec_E, rot=rot, mirror=mirror, centre=(rc, dc))
    P, N, _ = profiles(tr)
    ST_coh.append(stats_array(P))
ST_coh = np.array(ST_coh)
print(f'coherent nulls: {len(kept)}  {time.time() - t0:.1f}s', flush=True)

# ------------------------------------------------------------------ (b) independent random segment nulls
seg_pieces = []
for i in range(NSEG):
    m = (REAL.s >= SEG_EDGES[i] - 0.02) & (REAL.s <= SEG_EDGES[i + 1] + 0.02)
    seg_pieces.append((REAL.ra[m], REAL.dec[m]))
ST_real_seg = []          # each real segment as its own Trace piece (same machinery as the null)
P_real_seg = []
for i, (r, d) in enumerate(seg_pieces):
    tr = Trace(r, d)
    P, N, _ = profiles(tr, ngroups_full=False)
    P_real_seg.append(P[0]); ST_real_seg.append(stats_array(P)[0])
ST_real_seg = np.array(ST_real_seg); P_real_seg = np.array(P_real_seg)

NNULL_SEG = 1500
ST_seg = np.full((NNULL_SEG, M.shape[1], len(STATS)), np.nan)
P_seg = np.full((NNULL_SEG, NB, M.shape[1]), np.nan, np.float32)
seg_of = np.zeros(NNULL_SEG, int)
nfill, tries = 0, 0
radec_ok = np.c_[pra[cover], pdec[cover]]
while nfill < NNULL_SEG and tries < 20 * NNULL_SEG:
    tries += 1
    i = rng.integers(NSEG)
    r, d = seg_pieces[i]
    c = radec_ok[rng.integers(len(radec_ok))]
    tr = Trace(r, d, rot=rng.uniform(0, 360), mirror=bool(rng.integers(2)), centre=(c[0], c[1]))
    near = (np.abs(pra - c[0]) * np.cos(np.radians(c[1])) < 3.6) & (np.abs(pdec - c[1]) < 3.6)
    pix = np.where(near)[0]
    dd, ss, ins = tr.coords(pra[pix], pdec[pix], dmax=3.0)
    band = ins & (np.abs(dd) < 2.5)
    if band.sum() < 200:
        continue
    if cover[pix][band].mean() < 0.8:
        continue
    inner = ins & (np.abs(dd) < 1.5)
    if real_band[pix][inner].mean() > 0.05:
        continue
    P, N, _ = profiles(tr, ngroups_full=False, pix=pix)
    P_seg[nfill] = P[0]; ST_seg[nfill] = stats_array(P)[0]; seg_of[nfill] = i
    nfill += 1
    if nfill % 250 == 0:
        print(f'  seg nulls {nfill}  {time.time() - t0:.1f}s', flush=True)
ST_seg = ST_seg[:nfill]; P_seg = P_seg[:nfill]; seg_of = seg_of[:nfill]
print(f'segment nulls: {nfill} (tries {tries})  {time.time() - t0:.1f}s', flush=True)

np.savez_compressed(DER / 'verify_hiwham_hi_profiles.npz', cen=CEN, names=np.array(names), groups=np.array(GROUPS),
                    stats=np.array(STATS), P_real=P_real, N_real=N_real, ST_real=ST_real, ST_coh=ST_coh,
                    coh_placements=np.array([[a[0], a[1], a[2], a[3], a[4], a[5]] for a in kept], float),
                    ST_real_seg=ST_real_seg, P_real_seg=P_real_seg, ST_seg=ST_seg, P_seg=P_seg, seg_of=seg_of,
                    slices=np.array(slices), s_split=s_split, L=L)
print('saved', time.time() - t0)
