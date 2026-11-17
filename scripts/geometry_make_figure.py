"""figures/geometry_trace.png: Silk Nebula spine, strands, Halpha patch, LAMOST field, P1-P4 and curvature,
over the NSNS wide [OIII] map (left) and the registered SWONS composite with digitized LAMOST [OIII] fibres (right).
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_make_figure.py
"""
import json, sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, label, finish, FG, FGDIM, OIII_C, HA_C, WARN, OK, ACCENT, SUN, HALO
from silk_geometry import load_trace
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = '/Users/mickey/astronomy/silk-nebula'
D = f'{ROOT}/data/derived'
tr = load_trace()
z = np.load(f'{D}/geometry_maps.npz')
ra_g, de_g = z['ra_grid'], z['dec_grid']
ext = [ra_g[0], ra_g[-1], de_g[0], de_g[-1]]


def stretch(m, box, lo=3, hi=99.6, gamma=0.75):
    a, b = np.nanpercentile(m[box], [lo, hi])
    return np.clip((m - a) / (b - a), 0, 1) ** gamma


def galgrid(ax, ra_lim, de_lim, lstep=2, bstep=1, lab=True):
    """Dashed constant-l and dotted constant-b lines; labels placed manually along the bottom (l) and the
    east/left edge (b) so they never collide with panel text."""
    RA, DE = np.meshgrid(np.linspace(min(ra_lim) - 1, max(ra_lim) + 1, 300), np.linspace(min(de_lim) - 1, max(de_lim) + 1, 300))
    g = SkyCoord(RA * u.deg, DE * u.deg).galactic
    L, B = g.l.deg, g.b.deg
    ax.contour(RA, DE, L, levels=np.arange(180, 230, lstep), colors=['#b39ddb'], linewidths=0.5, linestyles='--', alpha=0.7)
    ax.contour(RA, DE, B, levels=np.arange(10, 30, bstep), colors=['#b39ddb'], linewidths=0.5, linestyles=':', alpha=0.7)
    if not lab:
        return
    x0, x1 = max(ra_lim), min(ra_lim); y0, y1 = min(de_lim), max(de_lim)
    # l labels along the bottom edge
    xs = np.linspace(x0, x1, 400); ys = np.full_like(xs, y0 + 0.02 * (y1 - y0))
    lb = SkyCoord(xs * u.deg, ys * u.deg).galactic.l.deg
    for lv in np.arange(180, 230, lstep):
        k = np.where(np.diff(np.sign(lb - lv)) != 0)[0]
        if len(k) and 0.06 < (x0 - xs[k[0]]) / (x0 - x1) < 0.94:
            label(ax, xs[k[0]], ys[0], f'l={lv:d}°', color='#d1c4e9', size=6.5, ha='center', va='bottom')
    # b labels along the east (left) edge
    ys = np.linspace(y0, y1, 400); xs = np.full_like(ys, x0 - 0.015 * (x0 - x1))
    bb = SkyCoord(xs * u.deg, ys * u.deg).galactic.b.deg
    for bv in np.arange(10, 30, bstep):
        k = np.where(np.diff(np.sign(bb - bv)) != 0)[0]
        if len(k) and 0.08 < (ys[k[0]] - y0) / (y1 - y0) < 0.9:
            label(ax, xs[0], ys[k[0]], f'b={bv:d}°', color='#d1c4e9', size=6.5, ha='left', va='bottom')


def strand_xy(s):
    p = tr['strands'][s]['points']
    return np.array([q['ra'] for q in p]), np.array([q['dec'] for q in p]), np.array([q['sigma_deg'] for q in p])


def small_circle(ra0, de0, rad, n=720):
    c = SkyCoord(ra0 * u.deg, de0 * u.deg)
    pa = np.linspace(0, 360, n) * u.deg
    p = c.directional_offset_by(pa, rad * u.deg)
    return p.ra.deg, p.dec.deg


fig, axes = new_sky_figure(13.2, 10.6, 'Silk Nebula · geometry of the [O III] arc',
                           'spine and strands measured on 4 star-registered images · RA/Dec J2000 with galactic l,b (lilac)',
                           ncols=2, gridspec_kw={'width_ratios': [0.43, 0.57]})
axL, axR = axes
# ---------------- left: full arc on NSNS wide [OIII]
RA, DE = np.meshgrid(ra_g, de_g)
boxL = (RA < 119) & (RA > 112.5) & (DE > 10) & (DE < 29.5)
axL.imshow(stretch(z['nsnsW_OIII'], boxL), origin='lower', extent=ext, cmap='gray', interpolation='bilinear')
axL.set_xlim(119.6, 112.4); axL.set_ylim(10.4, 29.6)
axL.set_aspect(1 / np.cos(np.radians(20)))
galgrid(axL, (112.4, 119.6), (10.4, 29.6), lstep=2, bstep=2)
# envelope
env = np.array(tr['envelope_polygon'])
axL.add_patch(Polygon(env, closed=True, fc=OIII_C, ec='none', alpha=0.10))
for s, c, lw, nm in (('E', OIII_C, 1.6, 'E spine'), ('W', OK, 1.3, 'W strand'), ('W2', '#80cbc4', 1.0, 'W2')):
    x, y, sg = strand_xy(s)
    axL.plot(x, y, color=c, lw=lw * 0.8, alpha=0.75)
    axL.errorbar(x, y, xerr=sg / np.cos(np.radians(y)), fmt='none', ecolor=c, elinewidth=0.6, alpha=0.6)
# Halpha patch
hp = np.array(tr['halpha_patch'])
axL.add_patch(Polygon(hp, closed=True, fc='none', ec=HA_C, lw=1.3))
label(axL, 114.95, 20.5, 'Hα patch', color=HA_C, size=8)
# LAMOST field
lf = tr['lamost_field']
cr, cd = small_circle(lf['ra'], lf['dec'], lf['radius_deg'])
axL.plot(cr, cd, color='#ffffff', lw=0.9, alpha=0.9)
label(axL, 118.75, 13.0, 'LAMOST\n2.5°', color='#ffffff', size=7.5)
# curvature circle
cv = tr['curvature']
xr, yr = small_circle(cv['center_ra'], cv['center_dec'], cv['radius_deg'], 4000)
k = (xr < 119.6) & (xr > 112.4) & (yr > 10.4) & (yr < 29.6)
axL.plot(np.where(k, xr, np.nan), np.where(k, yr, np.nan), color=SUN, lw=0.9, ls='--', alpha=0.9)
label(axL, 114.0, 17.15, f"to centre\n({cv['center_ra']:.1f}°, {cv['center_dec']:+.1f}°)\nR = {cv['radius_deg']:.1f}°",
      color=SUN, size=7, ha='center', va='bottom')
axL.annotate('', xy=(112.6, 16.75), xytext=(115.0, 16.75), arrowprops=dict(arrowstyle='->', color=SUN, lw=1.0))
# NSNS extensions
label(axL, 117.72, 12.05, 'O200.0b', color=FGDIM, size=6.5, ha='center')
axL.plot([117.72], [12.52], marker='+', color=FGDIM, ms=6)
for n, (r, d) in {k: (v['ra'], v['dec']) for k, v in tr['regions'].items()}.items():
    axL.plot(r, d, marker='o', ms=3, color=WARN, mec='k', mew=0.4)
label(axL, 116.0, 29.15, 'NSNS wide [O III]', color=FGDIM, size=7.5, ha='center')
axL.set_xlabel('RA (deg)'); axL.set_ylabel('Dec (deg)')
axL.grid(alpha=0.25, lw=0.4)
# ---------------- right: LAMOST field on SWONS composite + digitized [OIII] fibres
rgb = np.stack([z['swons_R'], z['swons_G'], z['swons_B']], -1)
boxR = (RA < 118.8) & (RA > 113.5) & (DE > 13) & (DE < 19.5)
img = np.zeros_like(rgb)
for i in range(3):
    img[..., i] = stretch(rgb[..., i], boxR, 1, 99.7, 0.9)
img = np.nan_to_num(img)
axR.imshow(img * 0.85, origin='lower', extent=ext, interpolation='bilinear')
axR.set_xlim(118.9, 113.4); axR.set_ylim(12.9, 19.5)
axR.set_aspect(1 / np.cos(np.radians(16.2)))
galgrid(axR, (113.4, 118.9), (12.9, 19.5), lstep=1, bstep=1)
axR.plot(cr, cd, color='#ffffff', lw=0.9, alpha=0.9)
lam = np.genfromtxt(f'{D}/geometry_lamost_oiii_points_digitized.csv', delimiter=',', skip_header=3, names=True)
sc = axR.scatter(lam['ra_deg'], lam['dec_deg'], c=lam['oiii_rel_intensity'], cmap='magma', s=11, marker='^', lw=0.3,
                 edgecolors='k', vmin=8, vmax=78, zorder=5)
for s, c, lw in (('E', OIII_C, 1.4), ('W', OK, 1.2)):
    x, y, sg = strand_xy(s)
    m = (y < 19.6)
    axR.plot(x[m], y[m], color=c, lw=lw, alpha=0.9, zorder=4)
axR.add_patch(Polygon(hp, closed=True, fc='none', ec=HA_C, lw=1.2))
for n, v in tr['regions'].items():
    axR.plot(v['ra'], v['dec'], marker='X', ms=7, color=WARN, mec='k', mew=0.6, zorder=6)
    dr = v.get('drawn_in_fig1b')
    if dr:
        o = np.array(dr['outline'])
        axR.add_patch(Polygon(o, closed=True, fc='none', ec=WARN, lw=0.8, ls='--', zorder=6))
    label(axR, v['ra'] - 0.12, v['dec'] + 0.08, n, color=WARN, size=8.5, zorder=7)
axR.plot(lf['ra'], lf['dec'], marker='*', ms=9, color='#ffffff', mec='k', mew=0.5, zorder=6)
label(axR, lf['ra'] - 0.1, lf['dec'] - 0.25, 'HIP 37797', color='#ffffff', size=7)
cbx = fig.add_axes([0.905, 0.12, 0.012, 0.25])
cb = fig.colorbar(sc, cax=cbx); cb.set_label('[O III] rel. intensity (digitized)', fontsize=7); cb.ax.tick_params(labelsize=6.5)
axR.set_xlabel('RA (deg)'); axR.set_ylabel('Dec (deg)')
axR.grid(alpha=0.25, lw=0.4)
label(axR, 114.05, 19.1, 'SWONS composite\n(Wu+2026 Fig 1)', color=FGDIM, size=7.5, ha='center')
# legend
hand = [Line2D([], [], color=OIII_C, lw=1.6, label='E spine (main strand)'),
        Line2D([], [], color=OK, lw=1.3, label='W strand'),
        Line2D([], [], color='#80cbc4', lw=1.0, label='W2 strand'),
        Line2D([], [], color=HA_C, lw=1.3, label='Hα patch (half max)'),
        Line2D([], [], color='#ffffff', lw=0.9, label='LAMOST field'),
        Line2D([], [], color=SUN, lw=0.9, ls='--', label='circle fit'),
        Line2D([], [], color=WARN, marker='X', ls='none', label='P1–P4 (Table 1)'),
        Line2D([], [], color=WARN, ls='--', lw=0.8, label='P1–P4 as drawn'),
        Line2D([], [], color='#e0e0e0', marker='^', ls='none', ms=5, label='LAMOST [O III] fibres')]
fig.legend(handles=hand, loc='lower center', ncol=5, frameon=True, fontsize=7.5, bbox_to_anchor=(0.5, 0.012))
fig.subplots_adjust(left=0.05, right=0.895, top=0.93, bottom=0.12, wspace=0.12)
finish(fig, f'{ROOT}/figures/geometry_trace.png', dpi=200, tight=False)
print('ok')
