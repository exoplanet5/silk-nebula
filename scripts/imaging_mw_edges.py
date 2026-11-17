"""Locate X-ray and HI edges relative to the Silk spine, and the Monogem X-ray ring edge by azimuth.

A) Along-spine edges (EQ grid, trace coordinates): the spine is cut into 1.5 deg segments; in each, the median
   cross-profile (0.1 deg bins, -4..+3 deg) is built and
     X-ray edge  = offset east of the interior peak (searched in -3.5..+0.5) where the profile first drops below
                   the mid-level between that peak and the exterior level (median at +0.5..+3);
     HI wall     = offset (searched -1..+3) where the [OIII]-velocity NHI (+11..+16 km/s) first rises above the
                   mid-level between the cavity floor (min in -1..+1) and the wall peak (max in 0..+3).
B) Monogem ring (GAL grid): radial profiles around PSR B0656+14 in 20-deg azimuth sectors (0.25-deg bins); the
   outer X-ray edge = radius beyond the ring peak (r 4-16 deg) where the profile drops below the mid-level between
   the peak and the outside level (median over r_peak+3..+6). Compared with the radius of the Silk spine points
   that fall in the same sector.
Output: data/derived/imaging_mw_edges.json
"""
import json
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import DER, load_grid, grid_coords, get_spine, MONOGEM_PULSAR, gal

T = np.load(DER / 'imaging_mw_eq_tracecoords.npz')
off, s, beyond, L = T['off'], T['s'], T['beyond'], float(T['L'])
spine = get_spine()
ra_s, dec_s, _ = spine
EDG = np.arange(-4.0, 3.0001, 0.1); CEN = 0.5 * (EDG[1:] + EDG[:-1])


def prof(img, sel):
    v = np.full(len(CEN), np.nan)
    for i in range(len(CEN)):
        m = sel & (off >= EDG[i]) & (off < EDG[i + 1]) & np.isfinite(img)
        if m.sum() >= 5:
            v[i] = np.median(img[m])
    # light smoothing (3-bin running mean) to stabilise edges
    k = np.convolve(np.nan_to_num(v, nan=np.nanmedian(v)), np.ones(3) / 3, mode='same')
    k[~np.isfinite(v)] = np.nan
    return k


def xray_edge(p):
    w = (CEN >= -3.5) & (CEN <= 0.5)
    if np.all(~np.isfinite(p[w])):
        return np.nan, np.nan, np.nan
    ipk = np.where(w)[0][np.nanargmax(p[w])]
    ext = np.nanmedian(p[(CEN > 0.5) & (CEN < 3.0)])
    mid = 0.5 * (p[ipk] + ext)
    for i in range(ipk, len(CEN)):
        if np.isfinite(p[i]) and p[i] < mid:
            # linear interpolation between i-1 and i
            x0, x1, y0, y1 = CEN[i - 1], CEN[i], p[i - 1], p[i]
            return float(x0 + (mid - y0) * (x1 - x0) / (y1 - y0)), float(p[ipk] - ext), float(CEN[ipk])
    return np.nan, float(p[ipk] - ext), float(CEN[ipk])


def hi_wall(p):
    fl = (CEN >= -1.0) & (CEN <= 1.0); pk = (CEN >= 0.0) & (CEN <= 3.0)
    if np.all(~np.isfinite(p[pk])):
        return np.nan, np.nan
    floor = np.nanmin(p[fl]); ipk = np.where(pk)[0][np.nanargmax(p[pk])]
    mid = 0.5 * (floor + p[ipk])
    ifl = np.where(fl)[0][np.nanargmin(p[fl])]
    for i in range(ifl, ipk + 1):
        if np.isfinite(p[i]) and p[i] > mid:
            x0, x1, y0, y1 = CEN[i - 1], CEN[i], p[i - 1], p[i]
            return float(x0 + (mid - y0) * (x1 - x0) / (y1 - y0)), float(p[ipk] / max(floor, 1e10))
    return np.nan, np.nan


out = {'method': __doc__, 'along_spine': {}, 'monogem_ring': {}}
segs = np.arange(0, L - 0.75, 1.5)
for band in ['ero_0206_clean', 'ero_0610_clean', 'rass_r12', 'rass_r45', 'murthy_nuv', 'murthy_fuv', 'halpha_f03']:
    img = load_grid(band, 'eq')
    recs = []
    for s0 in segs:
        sel = (s >= s0) & (s < s0 + 1.5) & ~beyond
        e, amp, xpk = xray_edge(prof(img, sel))
        dec_mid = float(np.interp(s0 + 0.75, np.linspace(0, L, len(dec_s)), dec_s))
        recs.append({'s0': float(s0), 'dec_mid': round(dec_mid, 2), 'edge_off_deg': e, 'peak_off_deg': xpk,
                     'interior_minus_exterior': amp})
    e = np.array([r['edge_off_deg'] for r in recs], float)
    ok = np.isfinite(e)
    out['along_spine'][band] = {'segments': recs, 'edge_median_deg': float(np.nanmedian(e)),
                                'edge_mean_deg': float(np.nanmean(e)), 'edge_sem_deg': float(np.nanstd(e) / np.sqrt(ok.sum())),
                                'n': int(ok.sum())}
    print(f"{band:16s} edge offset median {np.nanmedian(e):+.2f}  mean {np.nanmean(e):+.2f} +- {np.nanstd(e) / np.sqrt(ok.sum()):.2f} deg (n={ok.sum()})")
for band in ['galfa_nhi_oiii', 'galfa_nhi_red', 'galfa_nhi_p10p18', 'galfa_nhi_all']:
    img = load_grid(band, 'eq')
    recs = []
    for s0 in segs:
        sel = (s >= s0) & (s < s0 + 1.5) & ~beyond
        w, contrast = hi_wall(prof(img, sel))
        dec_mid = float(np.interp(s0 + 0.75, np.linspace(0, L, len(dec_s)), dec_s))
        recs.append({'s0': float(s0), 'dec_mid': round(dec_mid, 2), 'wall_off_deg': w, 'wall_over_floor': contrast})
    e = np.array([r['wall_off_deg'] for r in recs], float); ok = np.isfinite(e)
    out['along_spine'][band] = {'segments': recs, 'wall_median_deg': float(np.nanmedian(e)),
                                'wall_mean_deg': float(np.nanmean(e)), 'wall_sem_deg': float(np.nanstd(e) / np.sqrt(max(ok.sum(), 1))),
                                'n': int(ok.sum())}
    print(f"{band:16s} HI wall offset median {np.nanmedian(e):+.2f}  mean {np.nanmean(e):+.2f} (n={ok.sum()})")

# ---------------- B) Monogem ring around PSR B0656+14 (GAL grid)
lon, lat = grid_coords('gal')
pl, pb = gal(*MONOGEM_PULSAR['radec']); pl, pb = float(pl[0]), float(pb[0])
c0 = SkyCoord(l=pl * u.deg, b=pb * u.deg, frame='galactic')
cc = SkyCoord(l=lon.ravel() * u.deg, b=lat.ravel() * u.deg, frame='galactic')
rad = c0.separation(cc).deg.reshape(lon.shape)
paz = c0.position_angle(cc).deg.reshape(lon.shape)          # from Galactic north through increasing l
sl, sb = gal(ra_s, dec_s)
cs = SkyCoord(l=sl * u.deg, b=sb * u.deg, frame='galactic')
srad = c0.separation(cs).deg; spa = c0.position_angle(cs).deg
out['monogem_ring']['centre'] = {'name': 'PSR B0656+14', 'l': pl, 'b': pb}
out['monogem_ring']['silk_radius_deg'] = {'min': float(srad.min()), 'max': float(srad.max()), 'median': float(np.median(srad))}
out['monogem_ring']['silk_pa_deg'] = {'min': float(spa.min()), 'max': float(spa.max())}
redg = np.arange(0, 20.001, 0.25); rc = 0.5 * (redg[1:] + redg[:-1])
for band in ['ero_0206_clean', 'ero_0610_clean', 'rass_r12', 'rass_r45']:
    img = load_grid(band, 'gal')
    recs = []
    for pa0 in range(0, 360, 20):
        dpa = (paz - pa0) % 360
        sec = (dpa < 20) & np.isfinite(img)
        p = np.array([np.median(img[sec & (rad >= redg[i]) & (rad < redg[i + 1])]) if (sec & (rad >= redg[i]) & (rad < redg[i + 1])).sum() > 5 else np.nan
                      for i in range(len(rc))])
        ps = np.convolve(np.nan_to_num(p, nan=np.nanmedian(p)), np.ones(5) / 5, mode='same'); ps[~np.isfinite(p)] = np.nan
        w = (rc >= 4) & (rc <= 16)
        if np.isfinite(ps[w]).sum() < 10:
            continue
        ipk = np.where(w)[0][np.nanargmax(ps[w])]
        outside = np.nanmedian(ps[(rc > rc[ipk] + 3) & (rc < rc[ipk] + 6)])
        mid = 0.5 * (ps[ipk] + outside)
        edge = np.nan
        for i in range(ipk, len(rc)):
            if np.isfinite(ps[i]) and ps[i] < mid:
                edge = float(rc[i - 1] + (mid - ps[i - 1]) * (rc[i] - rc[i - 1]) / (ps[i] - ps[i - 1])); break
        insec = ((spa - pa0) % 360) < 20
        recs.append({'pa0': pa0, 'r_peak': float(rc[ipk]), 'r_edge': edge, 'peak_over_outside': float(ps[ipk] / outside) if outside > 0 else None,
                     'silk_r_in_sector': float(np.median(srad[insec])) if insec.any() else None,
                     'profile_r': rc.tolist(), 'profile': [None if not np.isfinite(x) else float(x) for x in ps]})
    out['monogem_ring'][band] = recs
    for r in recs:
        if r['silk_r_in_sector']:
            print(f"{band:14s} PA {r['pa0']:3d}-{r['pa0'] + 20:3d}: ring peak r={r['r_peak']:.2f} edge r={r['r_edge']:.2f}  Silk r={r['silk_r_in_sector']:.2f}")
(DER / 'imaging_mw_edges.json').write_text(json.dumps(out, indent=1, default=float))
print('Silk radius from PSR B0656+14:', out['monogem_ring']['silk_radius_deg'], 'PA', out['monogem_ring']['silk_pa_deg'])
