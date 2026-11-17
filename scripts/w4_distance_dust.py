"""w4_distance: 3D dust (Edenhofer+2024) around the Silk and the PGCC G199.09+19.95 clump, plus the SFD/IRIS
2-D column toward the clump.

Inputs
  data/raw/w4_distance_eden/eden_subcube.npz  (w4_distance_eden_fetch.py; nside 256, 90-750 pc)
  data/raw/skyview_imaging_mw/eq_SFD_Dust_Map.fits, eq_IRIS_100.fits, eq_Planck_857_I.fits (wave 1, on disk)
  data/derived/silk_trace.json
Outputs
  data/derived/w4_distance_dust.json        (clump on/off profiles, step distance, slices summary)
  data/derived/w4_distance_dust_pixels.npz  (per-pixel E in distance bins, for the figure)
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_dust.py
"""
import json
from pathlib import Path
import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
Z = np.load(ROOT / 'data/raw/w4_distance_eden/eden_subcube.npz')
mean, std = Z['mean'].astype(float), Z['std'].astype(float)
par, spar, ridx = Z['parents16'], Z['std_parents16'], Z['ridx']
rc, rb = Z['r_centers'][ridx], Z['r_bounds']
dr = np.diff(rb)[ridx]
pix = (par[:, None] * 256 + np.arange(256)[None, :]).ravel()
dens = mean.transpose(0, 2, 1).reshape(-1, len(ridx))           # (npix, nrad)
L, B = hp.pix2ang(256, pix, nest=True, lonlat=True)
eq = SkyCoord(L * u.deg, B * u.deg, frame='galactic').icrs
RA, DE = eq.ra.deg, eq.dec.deg

BINS = [(90, 150), (150, 200), (200, 250), (250, 320), (320, 400), (400, 600), (600, 750)]


def ebin(d, lo, hi):
    m = (rc >= lo) & (rc < hi)
    return (d[..., m] * dr[m]).sum(-1)


Ebins = {f'{lo}-{hi}': ebin(dens, lo, hi) for lo, hi in BINS}
np.savez_compressed(ROOT / 'data/derived/w4_distance_dust_pixels.npz', pix=pix, l=L, b=B, ra=RA, dec=DE,
                    **{f'E_{k}': v for k, v in Ebins.items()}, r=rc, dens_mean=dens.astype('f4'))

# ---------- clump on/off ----------
pg = SkyCoord(115.2847 * u.deg, 20.8334 * u.deg)
sep = SkyCoord(RA * u.deg, DE * u.deg).separation(pg).deg
on = sep < 0.20
off = (sep > 0.6) & (sep < 1.5)
cum = np.cumsum(dens * dr, axis=1)
on_prof, off_prof = dens[on].mean(0), np.median(dens[off], 0)
off_sd = np.std(dens[off], 0)
on_cum, off_cum = cum[on].mean(0), np.median(cum[off], 0)
off_cum_sd = np.std(cum[off], 0)
dcum = on_cum - off_cum
# posterior std toward the on pixels (std parent contains the PGCC?)
spix = (spar[:, None] * 256 + np.arange(256)[None, :]).ravel()
sd = std.transpose(0, 2, 1).reshape(-1, len(ridx))
on_std = []
for p in pix[on]:
    j = np.where(spix == p)[0]
    if len(j):
        on_std.append(sd[j[0]])
on_std = np.array(on_std)
# where does the on-off excess build up?
tot = dcum[-1]
d_half = float(np.interp(0.5 * tot, dcum, rc)) if tot > 0 else None
d_16 = float(np.interp(0.16 * tot, dcum, rc)) if tot > 0 else None
d_84 = float(np.interp(0.84 * tot, dcum, rc)) if tot > 0 else None
excess = on_prof - off_prof
kpk = int(np.argmax(excess))
# excess integrated in windows
win = {f'{lo}-{hi}': float(((on_prof - off_prof) * dr)[(rc >= lo) & (rc < hi)].sum()) for lo, hi in BINS}
win_sig = {f'{lo}-{hi}': float(np.sqrt(((off_sd * dr)[(rc >= lo) & (rc < hi)] ** 2).sum())) for lo, hi in BINS}

# ---------- 2-D column toward the clump (SFD, IRIS 100, Planck 857) ----------
def map_on_off(fn):
    h = fits.open(fn)[0]
    w = WCS(h.header)
    ny, nx = h.data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    sky = w.pixel_to_world(xx, yy)
    s = sky.separation(pg).deg
    d = h.data.astype(float)
    return float(np.nanmean(d[s < 0.1])), float(np.nanmedian(d[(s > 0.6) & (s < 1.5)])), \
        float(np.nanstd(d[(s > 0.6) & (s < 1.5)])), float(np.nanmax(d[s < 0.3])), h.header.get('BUNIT')


col = {}
for nm in ('SFD_Dust_Map', 'IRIS_100', 'Planck_857_I'):
    fn = ROOT / f'data/raw/skyview_imaging_mw/eq_{nm}.fits'
    if fn.exists():
        a, b_, c_, mx, un = map_on_off(fn)
        col[nm] = {'on_r<0.1': a, 'off_median_0.6-1.5': b_, 'off_sd': c_, 'max_r<0.3': mx, 'unit': un, 'excess': a - b_}

# ---------- 300-600 pc dust east vs west of the ridge ----------
tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
sra = np.array([p['ra'] for p in tr['spine']]); sde = np.array([p['dec'] for p in tr['spine']])
ridge_ra = np.interp(DE, sde, sra, left=np.nan, right=np.nan)
off_e = (RA - ridge_ra) * np.cos(np.radians(DE))          # + east
sides = {}
for dec0 in range(12, 29, 2):
    m = np.abs(DE - dec0 - 1) <= 1
    row = {}
    for nm, sel in (('W_3-1deg', (off_e > -3) & (off_e < -1)), ('ridge_+-0.5', np.abs(off_e) < 0.5),
                    ('E_1-3deg', (off_e > 1) & (off_e < 3)), ('E_3-5deg', (off_e > 3) & (off_e < 5))):
        mm = m & sel
        if mm.sum() < 3:
            row[nm] = None
            continue
        row[nm] = {'n': int(mm.sum()), **{f'E_{k}': round(float(np.median(Ebins[k][mm])), 4) for k in ('250-320', '320-400', '400-600')},
                   'peak_dist_300_650': round(float(rc[(rc >= 300) & (rc < 650)][np.argmax(np.median(dens[mm][:, (rc >= 300) & (rc < 650)], 0))]), 0)}
    sides[f'Dec{dec0}-{dec0 + 2}'] = row

out = {'date': '2026-11-17', 'source': 'Edenhofer+2024 mean_and_std_healpix.fits (Zenodo 10658339), range reads',
       'units': 'E (Zhang, Green & Rix 2023) ~ E(B-V)',
       'pgcc': {'ra': 115.2847, 'dec': 20.8334, 'l': 199.094, 'b': 19.947,
                'catalogue_distance': 'PGCC (Planck 2016 XXVIII, J/A+A/594/A28): DIST = 0.262 +- 0.121 kpc, DIST_OPTION 2 = optical extinction on SDSS DR7 (DistOEDR7), q_Dist = 2; no size/NH2/temperature (q_Flux = 3)',
                'n_on_pix_r<0.2': int(on.sum()), 'n_off_pix_0.6-1.5': int(off.sum()),
                'E_on_total_90_750': float(on_cum[-1]), 'E_off_total_90_750': float(off_cum[-1]),
                'E_off_total_sd': float(off_cum_sd[-1]), 'on_minus_off_total': float(tot),
                'on_posterior_std_total': (float(np.sqrt(((on_std * dr) ** 2).sum(1)).mean()) if len(on_std) else None),
                'excess_cumulative_16_50_84_pc': [d_16, d_half, d_84],
                'excess_density_peak_pc': float(rc[kpk]),
                'excess_by_window': win, 'off_scatter_by_window': win_sig,
                'profile_r_pc': rc[::4].round(1).tolist(), 'on_density': on_prof[::4].tolist(),
                'off_density_median': off_prof[::4].tolist(), 'off_density_sd': off_sd[::4].tolist()},
       'column_2d': col,
       'ridge_sides_300_600': sides}
json.dump(out, open(ROOT / 'data/derived/w4_distance_dust.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out['pgcc'].items() if not k.startswith(('profile', 'on_dens', 'off_dens'))}, indent=1))
print(json.dumps(col, indent=1))
for k, v in sides.items():
    print(k, {kk: (vv if vv is None else (vv['E_320-400'], vv['E_400-600'], vv['peak_dist_300_650'])) for kk, vv in v.items()})
