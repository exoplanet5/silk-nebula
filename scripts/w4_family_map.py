"""w4_family_map.py -- wave-4 `family` step 1: NSNS DR0.2 [O III] / Halpha mosaic of the whole Monogem Ring
region (order-3 HiPS tiles already on disk, 0.86'/px) on a 3'/px TAN grid, plus spectroscopic coverage
(DESI DR1 tiles, SDSS DR17 plates) and the Monogem X-ray circles.

Inputs (on disk): data/raw/nsns_hips/{oiii,halpha,sii,sn_oiii}/Norder3/*.fits,
                  data/raw/w4_family_tap/tiles_monogem.csv, data/raw/w4_family_sdss/platex_monogem.csv,
                  data/derived/silk_trace.json
Outputs: data/derived/w4_family_nsns_wide.npz (maps + WCS header), figures/w4_family_overview.png
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_map.py
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import csv
import json
import numpy as np
from scipy.ndimage import median_filter, gaussian_filter
from imaging_optical_common import tan_wcs, hips_mosaic, pixel_radec
from silk_style import *

ROOT = '/Users/mickey/astronomy/silk-nebula'
PIX = 3.0
W, SHP = tan_wcs(106.0, 13.0, 42.0, 40.0, PIX)


def circle(ra0, dec0, r, n=721):
    """Small circle of angular radius r (deg) about (ra0, dec0)."""
    t = np.radians(np.linspace(0, 360, n))
    d0, a0, rr = np.radians(dec0), np.radians(ra0), np.radians(r)
    dec = np.arcsin(np.sin(d0) * np.cos(rr) + np.cos(d0) * np.sin(rr) * np.cos(t))
    ra = a0 + np.arctan2(np.sin(t) * np.sin(rr) * np.cos(d0), np.cos(rr) - np.sin(d0) * np.sin(dec))
    return np.degrees(ra), np.degrees(dec)


def main():
    maps = {n: hips_mosaic(n, 3, W, SHP) for n in ['oiii', 'halpha', 'sii', 'sn_oiii']}
    hdr = W.to_header().tostring()
    np.savez_compressed(f'{ROOT}/data/derived/w4_family_nsns_wide.npz', wcs_header=hdr, pix_arcmin=PIX, **maps)

    tr = json.load(open(f'{ROOT}/data/derived/silk_trace.json'))
    sp = tr['spine']
    sra = np.array([p['ra'] for p in sp]); sdec = np.array([p['dec'] for p in sp])

    o3 = maps['oiii'].copy()
    bg = median_filter(np.nan_to_num(o3, nan=np.nanmedian(o3)), size=41)   # 2 deg median background
    hp_o3 = gaussian_filter(np.nan_to_num(o3 - bg), 1.0)
    fig, axes = new_sky_figure(16, 8.2, 'Monogem Ring periphery · NSNS [O III] and spectroscopic coverage',
                               'NSNS DR0.2 order-3 tiles, 2° median-subtracted; DESI DR1 tiles (r 1.6°) and SDSS plates (r 1.49°)',
                               ncols=2, subplot_kw=dict(projection=W))
    for k, (ax, name) in enumerate(zip(axes, ['oiii', 'halpha'])):
        m = maps[name]
        b2 = median_filter(np.nan_to_num(m, nan=np.nanmedian(m)), size=41)
        img = gaussian_filter(np.nan_to_num(m - b2), 1.0)
        img[~np.isfinite(m)] = np.nan
        vmax = 1.5 if name == 'oiii' else 4.0
        ax.imshow(img, origin='lower', cmap='gray', vmin=-0.3 * vmax, vmax=vmax)
        tx = ax.get_transform('world')
        for (ra0, dec0, r, c, lab) in [(104.15, 15.30, 11.0, WARN, 'K24 two-SNR'),
                                       (107.08, 13.50, 11.97, SUN, 'K24 single')]:
            cr, cd = circle(ra0, dec0, r)
            ax.plot(cr, cd, color=c, lw=0.9, ls='--', transform=tx, label=lab)
            ax.plot([ra0], [dec0], '+', color=c, ms=8, transform=tx)
        ax.plot(sra, sdec, color=OIII_C, lw=1.2, transform=tx, label='Silk spine')
        if name == 'oiii':
            for row in csv.DictReader(open(f'{ROOT}/data/raw/w4_family_tap/tiles_monogem.csv')):
                cr, cd = circle(float(row['tilera']), float(row['tiledec']), 1.6, 121)
                ax.plot(cr, cd, color=OK, lw=0.4, alpha=0.6, transform=tx)
            lines = open(f'{ROOT}/data/raw/w4_family_sdss/platex_monogem.csv').read().split('\n', 1)[1]
            seen = set()
            for row in csv.DictReader(lines.splitlines()):
                key = (round(float(row['ra']), 2), round(float(row['dec']), 2))
                if key in seen:
                    continue
                seen.add(key)
                cr, cd = circle(float(row['ra']), float(row['dec']), 1.49, 121)
                ax.plot(cr, cd, color=DANGER, lw=0.4, alpha=0.6, transform=tx)
        # Fesen SNRs (Galactic -> ICRS done offline in the note; positions from Fesen+26 names)
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        for l, b, r, lab in [(206.6, 6.1, 1.8, 'G206.6+6.1'), (203.1, 6.6, 0.8, 'G203.1+6.6'),
                             (191.4, 11.1, 1.4, 'G191.4+11.1'), (190.4, 12.5, 3.3, 'G190.4+12.5'),
                             (191.5, 5.0, 2.85, 'Gemini Hα ring'), (195.9, 2.2, 0.5, 'G195.9+2.2')]:
            c = SkyCoord(l * u.deg, b * u.deg, frame='galactic').icrs
            cr, cd = circle(c.ra.deg, c.dec.deg, r, 181)
            ax.plot(cr, cd, color='#c792ea', lw=0.8, transform=tx)
            label(ax, c.ra.deg, c.dec.deg + r + 0.3, lab, size=7, color='#e0c6ff', transform=tx, ha='center')
        ax.coords[0].set_axislabel('RA'); ax.coords[1].set_axislabel('Dec')
        ax.coords.grid(color='#3a4450', alpha=0.5, lw=0.4)
        ax.set_title('[O III] 5007' if name == 'oiii' else 'Hα', color=FG)
        if k == 0:
            ax.legend(loc='lower left', fontsize=7, ncol=3)
    finish(fig, f'{ROOT}/figures/w4_family_overview.png')


if __name__ == '__main__':
    main()
