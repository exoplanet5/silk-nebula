"""verify_desi_platescale.py -- empirical DESI focal-plane plate scale from desi_dr1.fiberassign (wave-2, desi).

For every fibre in data/raw/datalab_tap/fiberassign_box.csv, the requested focal-plane position (FIBERASSIGN_X/Y, mm)
and the target's gnomonic offset from its tile centre (arcsec) give the radial mapping r_fp(theta). A fit
r_fp = a1 t + a3 t^3 + a5 t^5 + a7 t^7 over all tiles yields the radial (dr/dtheta) and azimuthal (r/theta) plate
scales; the 107-um fibre core then subtends d = 107 / sqrt(ps_rad * ps_az) arcsec. Compares with the 1.52" mean
(DESI-6043) used for surface brightness, and gives the per-fibre solid-angle spread.
Output: data/derived/verify_desi_platescale.json
"""
import csv
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import ROOT, D  # noqa: E402

tiles = {int(r['tileid']): (float(r['tilera']), float(r['tiledec'])) for r in
         csv.DictReader(open(ROOT / 'data/raw/datalab_tap/tiles_box.csv'))}
ra, de, x, y, t = [], [], [], [], []
for r in csv.DictReader(open(ROOT / 'data/raw/datalab_tap/fiberassign_box.csv')):
    try:
        x.append(float(r['fiberassign_x']))
        y.append(float(r['fiberassign_y']))
    except ValueError:
        continue
    ra.append(float(r['target_ra']))
    de.append(float(r['target_dec']))
    t.append(int(r['tileid']))
ra, de, x, y, t = map(np.array, (ra, de, x, y, t))
r0 = np.array([tiles[k][0] for k in t])
d0 = np.array([tiles[k][1] for k in t])
a, d, a0, dd0 = map(np.radians, (ra, de, r0, d0))
cosc = np.sin(dd0) * np.sin(d) + np.cos(dd0) * np.cos(d) * np.cos(a - a0)
theta = np.degrees(np.arccos(np.clip(cosc, -1, 1))) * 3600          # true angular distance (arcsec)
rfp = np.hypot(x, y) * 1000                                         # um
ok = np.isfinite(rfp) & (theta > 60) & (rfp > 0)
tt = theta[ok] / 3600.0
X = np.c_[tt, tt ** 3, tt ** 5, tt ** 7]
c, *_ = np.linalg.lstsq(X, rfp[ok], rcond=None)
res = rfp[ok] - X @ c
th = np.linspace(0.02, 1.62, 161)
rr = np.c_[th, th ** 3, th ** 5, th ** 7] @ c
drr = np.c_[np.ones_like(th), 3 * th ** 2, 5 * th ** 4, 7 * th ** 6] @ c
ps_rad = drr / 3600.0                                               # um / arcsec
ps_az = rr / th / 3600.0
dfib = 107.0 / np.sqrt(ps_rad * ps_az)
# fibre-weighted mean over the actual fibre radii in the box
tf = theta[ok] / 3600.0
d_f = np.interp(tf, th, dfib)
area_mean = np.mean(np.pi * (d_f / 2) ** 2)
out = dict(n_fibres=int(ok.sum()), coef_um_per_deg_powers=c.tolist(), fit_rms_um=float(np.std(res)),
           ps_rad_centre=float(ps_rad[0]), ps_az_centre=float(ps_az[0]), ps_rad_edge=float(ps_rad[-1]),
           ps_az_edge=float(ps_az[-1]), d_fibre_centre_arcsec=float(dfib[0]), d_fibre_edge_arcsec=float(dfib[-1]),
           d_fibre_min=float(dfib.min()), d_fibre_max=float(dfib.max()),
           d_fibre_mean_over_fibres=float(np.mean(d_f)), d_equiv_from_mean_area=float(2 * np.sqrt(area_mean / np.pi)),
           area_rel_spread_rms=float(np.std(np.pi * (d_f / 2) ** 2) / area_mean),
           area_ratio_mean_vs_1p52=float(area_mean / (np.pi * 0.76 ** 2)),
           curve=dict(theta_deg=th[::10].tolist(), d_fibre=dfib[::10].round(4).tolist()))
(D / 'verify_desi_platescale.json').write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != 'curve'}, indent=1))
