"""Figures: multiwavelength panels on the common EQ grid (RA 108-125, Dec 4-30) and the wide Galactic grid
(l 185-225, b 0-35), with the Silk trace (E spine solid, W/W2 strands dashed, Halpha patch red), P1-P4.

Outputs: figures/imaging_mw_eq_panels.png, figures/imaging_mw_gal_panels.png
"""
import json
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FGDIM, WARN, OK, SUN
from imaging_mw_common import (load_grid, overlay_eq, eq_axes, gal_axes, EQ_EXTENT, GAL_EXTENT, get_spine,
                               spine_gal, MONOGEM_PULSAR, GEMINGA, gal, FIG)


def smooth_nan(a, sig):
    m = np.isfinite(a)
    w = gaussian_filter(m.astype(float), sig)
    s = gaussian_filter(np.where(m, a, 0.0), sig) / np.where(w > 0, w, np.nan)
    s[w < 0.3] = np.nan
    return s


def show(ax, a, extent, lo=2, hi=98, cmap='gray', sig=0):
    if sig:
        a = smooth_nan(a, sig)
    v0, v1 = np.nanpercentile(a, [lo, hi])
    im = ax.imshow(a, origin='lower', extent=extent, vmin=v0, vmax=v1, cmap=cmap, interpolation='nearest')
    return im


EQP = [  # grid, title, cmap, smoothing sigma (px), scale, unit
    ('murthy_fuv', 'GALEX FUV · Murthy 2014', 'magma', 1.0, 1, 'CU'),
    ('murthy_nuv', 'GALEX NUV · Murthy 2014', 'magma', 1.0, 1, 'CU'),
    ('rass_r12', 'ROSAT R12 · 1/4 keV', 'inferno', 1.0, 1, '10⁻⁶ c s⁻¹ arcmin⁻²'),
    ('rass_r45', 'ROSAT R45 · 3/4 keV', 'inferno', 2.0, 1, '10⁻⁶ c s⁻¹ arcmin⁻²'),
    ('ero_0206_clean', 'eRASS1 0.2–0.6 keV', 'inferno', 1.0, 1, 'c s⁻¹ deg⁻²'),
    ('ero_0610_clean', 'eRASS1 0.6–1.0 keV', 'inferno', 1.5, 1, 'c s⁻¹ deg⁻²'),
    ('iris100', 'IRIS 100 µm', 'cividis', 0, 1, 'MJy sr⁻¹'),
    ('halpha_f03', 'Hα · Finkbeiner 2003', 'gray', 0, 1, 'R'),
    ('galfa_nhi_oiii', 'GALFA-HI · V_LSR +11…+16', 'viridis', 1.0, 1e-19, '10¹⁹ cm⁻²'),
    ('galfa_nhi_red', 'GALFA-HI · V_LSR −1…+6', 'viridis', 1.0, 1e-19, '10¹⁹ cm⁻²'),
    ('stockert1420', 'Stockert 1420 MHz', 'cubehelix', 0, 1, 'mK'),
    ('haslam408', 'Haslam 408 MHz', 'cubehelix', 0, 1, 'K'),
]


def eq_fig():
    sp = get_spine()
    fig, axes = new_sky_figure(17.0, 15.2, 'Silk Nebula · multiwavelength counterparts',
                               'RA/Dec plate carrée; E spine solid, W strands dashed, Hα patch red, LAMOST field grey',
                               nrows=3, ncols=4)
    for ax, (nm, title, cmap, sig, sc, unit) in zip(axes.ravel(), EQP):
        a = load_grid(nm, 'eq') * sc
        im = show(ax, a, EQ_EXTENT, cmap=cmap, sig=sig)
        overlay_eq(ax, sp, color='#ffffff', lw=0.9, label_regions=(nm == 'murthy_fuv'))
        eq_axes(ax)
        ax.set_title(title, loc='left', fontsize=9.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01)
        cb.ax.tick_params(labelsize=6.5); cb.set_label(unit, fontsize=7)
        ax.tick_params(labelsize=7)
    # mark Monogem pulsar direction (off-grid, to the west)
    finish(fig, FIG / 'imaging_mw_eq_panels.png')


GALP = [
    ('rass_r12', 'ROSAT R12 · 1/4 keV', 'inferno', 1.0, 1),
    ('rass_r45', 'ROSAT R45 · 3/4 keV', 'inferno', 1.5, 1),
    ('ero_0206_clean', 'eRASS1 0.2–0.6 keV', 'inferno', 0.8, 1),
    ('ero_0610_clean', 'eRASS1 0.6–1.0 keV', 'inferno', 1.5, 1),
    ('murthy_fuv', 'GALEX FUV · Murthy 2014', 'magma', 0.7, 1),
    ('halpha_f03', 'Hα · Finkbeiner 2003', 'gray', 0, 1),
    ('iris100', 'IRIS 100 µm', 'cividis', 0, 1),
    ('lab_v+10', 'LAB HI · V_LSR +10…+20', 'viridis', 0, 1),
    ('hi4pi_nhi', 'HI4PI N(HI)', 'viridis', 0, 1),
    ('haslam408', 'Haslam 408 MHz', 'cubehelix', 0, 1),
    ('planck030_pi', 'Planck 30 GHz polarized', 'cubehelix', 1.0, 1),
    ('fermi5', 'Fermi-LAT band 5', 'gray', 1.0, 1),
]


def gal_fig():
    sp = get_spine()
    ls, bs = spine_gal(sp)
    import json as _j
    tr = sp[2].get('raw') or {}
    cur = (tr.get('curvature') or {}).get('fits', {}).get('E_full', {})
    psr_l, psr_b = gal(*MONOGEM_PULSAR['radec']); gem_l, gem_b = gal(*GEMINGA['radec'])
    fig, axes = new_sky_figure(17.0, 13.6, 'Silk Nebula · Galactic context',
                               'l/b plate carrée; Silk E spine (white), PSR B0656+14 (+), Geminga (×), E-spine small-circle fit (dotted)',
                               nrows=3, ncols=4)
    for ax, (nm, title, cmap, sig, sc) in zip(axes.ravel(), GALP):
        a = load_grid(nm, 'gal') * sc
        lo, hi = (2, 98) if nm not in ('hi4pi_nhi', 'lab_v+10', 'iris100', 'fermi5', 'halpha_f03') else (2, 92)
        show(ax, a, GAL_EXTENT, lo=lo, hi=hi, cmap=cmap, sig=sig)
        ax.plot(ls, bs, color='#ffffff', lw=1.2)
        ax.plot(psr_l, psr_b, '+', color=SUN, ms=9, mew=1.4)
        ax.plot(gem_l, gem_b, 'x', color=OK, ms=7, mew=1.4)
        if cur:
            c = SkyCoord(cur['center_ra'] * u.deg, cur['center_dec'] * u.deg)
            pa = np.linspace(0, 360, 721) * u.deg
            ring = c.directional_offset_by(pa, cur['radius_deg'] * u.deg).galactic
            ax.plot(ring.l.deg, ring.b.deg, color='#ffffff', lw=0.7, ls=':')
        if nm == 'rass_r12':
            label(ax, psr_l[0] - 0.8, psr_b[0] - 1.5, 'PSR B0656+14', color=SUN, size=7)
            label(ax, gem_l[0] - 0.8, gem_b[0] + 0.8, 'Geminga', color=OK, size=7)
            label(ax, ls[0] + 1.0, bs[0] - 2.0, 'Silk', color='#ffffff', size=8)
        gal_axes(ax)
        ax.set_title(title, loc='left', fontsize=9.5)
        ax.tick_params(labelsize=7)
    finish(fig, FIG / 'imaging_mw_gal_panels.png')


if __name__ == '__main__':
    eq_fig(); gal_fig(); print('ok')
