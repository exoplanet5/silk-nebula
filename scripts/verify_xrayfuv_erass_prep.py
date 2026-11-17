"""verify_xrayfuv: build an independent eRASS1 pixel table around the Monogem Ring.

Reads the raw eRASS1 DR1 half-sky AIT rate maps (cts/s per 3'x3' pixel, TM1-4+6 summed), cuts
l 180.05-226, b -12..+38, converts to cts/s/deg^2 (x400), masks catalogued eRASS1 point/extended
sources (Merloni+2024 erass1-m via VizieR CfA) with radii set from a stacked source profile, then
masks residual >5 sigma pixels relative to a 1-deg running median. Adds HI4PI N(HI) and an
empirical exposure map (median catalogue MLExp1 in 1-deg cells). Output:
data/derived/verify_xrayfuv_erass_pixels.npz
"""
import sys, glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import ROOT, DER, read_vizier_tsv

RAW = ROOT / 'data/raw/erosita_dr1'
BANDS = {'b0204': 'e01_TM8_0.200_0.400keV_3arcm_c947_AIT_sum_Ratesmap.fits',
         'b0406': 'e01_TM8_0.400_0.600keV_3arcm_c947_AIT_sum_Ratesmap.fits',
         'b0610': 'e01_TM8_0.600_1.000keV_3arcm_c947_AIT_sum_Ratesmap.fits',
         'b1023': 'e01_TM8_1.000_2.300keV_3arcm_c947_AIT_sum_Ratesmap.fits'}
L0, L1, B0, B1 = 180.05, 226.0, -12.0, 38.0

h = fits.getheader(RAW / BANDS['b0204']); w = WCS(h)
L, B = np.meshgrid(np.linspace(L0, L1, 400), np.linspace(B0, B1, 400))
x, y = w.world_to_pixel_values(L, B)
x0, x1 = int(np.nanmin(x)) - 3, min(int(np.nanmax(x)) + 4, h['NAXIS1'])
y0, y1 = int(np.nanmin(y)) - 3, min(int(np.nanmax(y)) + 4, h['NAXIS2'])
print('cut x', x0, x1, 'y', y0, y1)
maps = {}
for k, f in BANDS.items():
    d = fits.getdata(RAW / f, memmap=True)
    maps[k] = np.array(d[y0:y1, x0:x1], dtype=np.float64) * 400.0  # cts/s/deg^2
ny, nx = maps['b0204'].shape
yy, xx = np.mgrid[y0:y1, x0:x1]
gl, gb = w.pixel_to_world_values(xx, yy)
inreg = np.isfinite(gl) & (gl >= L0) & (gl <= L1) & (gb >= B0) & (gb <= B1)
for k in maps:
    inreg &= np.isfinite(maps[k])
c = SkyCoord(gl[inreg] * u.deg, gb[inreg] * u.deg, frame='galactic').icrs
ra = np.full(gl.shape, np.nan); dec = np.full(gl.shape, np.nan)
ra[inreg] = c.ra.deg; dec[inreg] = c.dec.deg

# ---------- catalogue ----------
cat = {}
for f in sorted(glob.glob(str(ROOT / 'data/raw/verify_xrayfuv_erass1cat/erass1m_ra*.tsv'))):
    t = read_vizier_tsv(f)
    for kk, v in t.items():
        cat.setdefault(kk, []).append(v)
cat = {kk: np.concatenate(v) for kk, v in cat.items()}
print('catalogue rows', len(cat['RA_ICRS']))
sx, sy = w.world_to_pixel_values(cat['GLON'], cat['GLAT'])
sx = np.round(sx).astype(int) - x0; sy = np.round(sy).astype(int) - y0
inside = (sx >= 0) & (sx < nx) & (sy >= 0) & (sy < ny)
tot = maps['b0204'] + maps['b0406'] + maps['b0610'] + maps['b1023']

# Median-stacked source profiles (exploration run 2026-11-16, see notes) show the map's response
# to a catalogued source is small and compact: median excess at r=0/1/2/3 px (3') of
# ~1.0/0.7/0.3/0.1 cts/s/deg^2 for MLRate1 0.05-0.1, ~3.5/2.5/1.3/0.6 for 1-5 cts/s, ~18/10/3.7/1.3
# for 5-50 cts/s (total band level ~13.6). Mask radii (pixels) by MLRate1 accordingly:
med_tot = np.nanmedian(tot[inreg])
mask = np.zeros_like(tot, dtype=bool)
rate = cat['MLRate1']; ext = cat['EXT']
use = inside & ((rate >= 0.1) | (ext > 0))
for i, j, r, e in zip(sx[use], sy[use], rate[use], ext[use]):
    a = 2 if r < 0.3 else 3 if r < 1 else 4 if r < 5 else 6 if r < 50 else 15
    if e > 0:
        a = max(a, int(np.ceil(3 * e / 180.0)))
    j0, j1, i0, i1 = max(j - a, 0), min(j + a + 1, ny), max(i - a, 0), min(i + a + 1, nx)
    sub = np.hypot(*np.mgrid[j0 - j:j1 - j, i0 - i:i1 - i]) <= a
    mask[j0:j1, i0:i1] |= sub
prof = np.zeros(1); n = int(use.sum())
print('catalogue mask fraction', mask[inreg].mean())

# residual sigma clip on 0.2-0.6 keV and total, relative to 1 deg running median
soft = maps['b0204'] + maps['b0406']
def resid_mask(img):
    tmp = np.where(inreg & ~mask, img, np.nan)
    med = ndimage.generic_filter(np.nan_to_num(tmp, nan=np.nanmedian(tmp)), np.median, size=7, mode='nearest') if False else None
    # fast robust local level: median filter on a 2x-binned image
    filled = np.where(np.isfinite(tmp), tmp, np.nanmedian(tmp))
    medf = ndimage.median_filter(filled, size=21, mode='nearest')
    madf = ndimage.median_filter(np.abs(filled - medf), size=21, mode='nearest') * 1.4826
    return (img - medf) > 5 * madf
rm = resid_mask(soft)
print('residual-clip mask fraction (extra)', (rm & ~mask)[inreg].mean())
good = inreg & ~mask & ~rm

# HI4PI N(HI) (wide SkyView CAR map) -> pixels
hh = fits.open(ROOT / 'data/raw/verify_xrayfuv_skyview/wide_hi4pi_nhi_l203b13_car.fits')[0]
wh = WCS(hh.header)
px, py = wh.world_to_pixel_values(np.where(inreg, gl, 0), np.where(inreg, gb, 0))
nhi = ndimage.map_coordinates(hh.data.astype(float), [py, px], order=1, mode='nearest')
nhi[~inreg] = np.nan
# SFD E(B-V)
hs = fits.open(ROOT / 'data/raw/verify_xrayfuv_skyview/wide_sfd_ebv_l203b13_car.fits')[0]
ebv = ndimage.map_coordinates(hs.data.astype(float), [py, px], order=1, mode='nearest'); ebv[~inreg] = np.nan

# empirical exposure: median MLExp1 of catalogue sources in 1-deg (20 pix) cells, smoothed
expo_sum = np.zeros_like(tot); expo_n = np.zeros_like(tot)
ok = inside & np.isfinite(cat['MLExp1'])
np.add.at(expo_sum, (sy[ok], sx[ok]), cat['MLExp1'][ok]); np.add.at(expo_n, (sy[ok], sx[ok]), 1)
k = 20
expo = ndimage.uniform_filter(expo_sum, k) / np.maximum(ndimage.uniform_filter(expo_n, k), 1e-9)
expo[ndimage.uniform_filter(expo_n, k) * k * k < 3] = np.nan
expo[~inreg] = np.nan

sel = inreg
np.savez_compressed(DER / 'verify_xrayfuv_erass_pixels.npz',
                    ra=ra[sel].astype(np.float32), dec=dec[sel].astype(np.float32),
                    gl=gl[sel].astype(np.float32), gb=gb[sel].astype(np.float32),
                    b0204=maps['b0204'][sel].astype(np.float32), b0406=maps['b0406'][sel].astype(np.float32),
                    b0610=maps['b0610'][sel].astype(np.float32), b1023=maps['b1023'][sel].astype(np.float32),
                    good=good[sel], catmask=mask[sel], nhi=nhi[sel].astype(np.float32),
                    ebv=ebv[sel].astype(np.float32), expo=expo[sel].astype(np.float32),
                    iy=(yy[sel] - y0).astype(np.int16), ix=(xx[sel] - x0).astype(np.int16),
                    stack_prof=prof, n_stack=n, med_tot=med_tot)
print('saved', sel.sum(), 'pixels; good fraction', good[sel].mean())
print('median soft 0.2-0.6 keV', np.nanmedian(soft[good]), 'cts/s/deg2; median exposure', np.nanmedian(expo[good]))
