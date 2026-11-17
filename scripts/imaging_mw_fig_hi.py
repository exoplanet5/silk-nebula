"""Figure: GALFA-HI DR2 channel maps (4 km/s bins, V_LSR -10..+30) and position-velocity diagrams across the Silk.

PV diagram: median Tb in bins of perpendicular offset from the E spine (0.1 deg) for each 0.74 km/s channel,
over the south (Dec 12-18) and north (Dec 18-28) parts of the spine. Output figures/imaging_mw_hi_channels.png
and data/derived/imaging_mw_hi_pv.npz
"""
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM
from imaging_mw_common import DER, FIG, EQ_EXTENT, overlay_eq, eq_axes, get_spine, spine_s_of

h = fits.open(DER / 'imaging_mw_grids' / 'imaging_mw_galfa_cube3.fits')
cube = h[0].data.astype(float); hd = h[0].header
v = hd['CRVAL3'] + (np.arange(cube.shape[0]) + 1 - hd['CRPIX3']) * hd['CDELT3']
T = np.load(DER / 'imaging_mw_eq_tracecoords.npz')
off, s, beyond, L = T['off'], T['s'], T['beyond'], float(T['L'])
sp = get_spine()
s18 = spine_s_of(116.8, 18.0, sp)[0]

# PV diagrams
edges = np.arange(-4.0, 4.0001, 0.1); cen = 0.5 * (edges[1:] + edges[:-1])
vsel = (v >= -15) & (v <= 35)
pv = {}
for seg, (lo, hi) in {'south': (0, s18), 'north': (s18, L)}.items():
    sel = (s >= lo) & (s <= hi) & ~beyond
    img = np.full((vsel.sum(), len(cen)), np.nan)
    for i in range(len(cen)):
        m = sel & (off >= edges[i]) & (off < edges[i + 1])
        if m.sum() > 5:
            img[:, i] = np.nanmedian(cube[vsel][:, m], axis=1)
    pv[seg] = img
np.savez_compressed(DER / 'imaging_mw_hi_pv.npz', v=v[vsel], offset=cen, south=pv['south'], north=pv['north'])

bins = [(-10, -6), (-6, -2), (-2, 2), (2, 6), (6, 10), (10, 14), (14, 18), (18, 22), (22, 26), (26, 30)]
fig, axes = new_sky_figure(17.0, 15.5, 'Silk Nebula · GALFA-HI channel maps',
                           'Tb averaged over 4 km/s (V_LSR); [O III] at +11…+16, Hα/[N II]/[S II] at −1…+6 km/s; blank: RA > 119.3 tiles not obtained',
                           nrows=3, ncols=4)
ax_all = axes.ravel()
for ax, (a, b) in zip(ax_all, bins):
    k = (v >= a) & (v < b)
    m = np.nanmean(cube[k], axis=0)
    v0, v1 = np.nanpercentile(m, [1, 99.5])
    im = ax.imshow(m, origin='lower', extent=EQ_EXTENT, cmap='viridis', vmin=0, vmax=max(v1, 1.0), interpolation='nearest')
    overlay_eq(ax, sp, color='#ffffff', lw=0.8, label_regions=False, lamost=False)
    eq_axes(ax)
    col = OIII_C if (a < 16 and b > 11) else (HA_C if (a < 6 and b > -1) else '#cfd6dd')
    ax.set_title(f'{a:+d} … {b:+d} km/s', loc='left', fontsize=9.5, color=col)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01); cb.ax.tick_params(labelsize=6.5); cb.set_label('K', fontsize=7)
    ax.tick_params(labelsize=7)
for ax, seg in zip(ax_all[10:], ['south', 'north']):
    img = pv[seg]
    vv = v[vsel]
    im = ax.imshow(img, origin='lower', aspect='auto', extent=[cen[0] - 0.05, cen[-1] + 0.05, vv[0], vv[-1]],
                   cmap='viridis', vmin=0, vmax=np.nanpercentile(img, 99.5), interpolation='nearest')
    ax.axvline(0, color='#ffffff', lw=0.8)
    ax.axhspan(11.2, 15.5, color=OIII_C, alpha=0.18, lw=0)
    ax.axhspan(-1.0, 6.1, color=HA_C, alpha=0.15, lw=0)
    label(ax, -3.9, 13.3, '[O III]', color=OIII_C, size=7.5, va='center')
    label(ax, -3.9, 2.5, 'Hα, [N II], [S II]', color=HA_C, size=7.5, va='center')
    ax.set_xlim(-4, 3.0 if seg == 'south' else 4)
    ax.set_xlabel('offset from spine (deg, + east)'); ax.set_ylabel('V_LSR (km/s)')
    ax.set_title(f'PV across the spine · {seg} ({"Dec 12–18" if seg == "south" else "Dec 18–28"})', loc='left', fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01); cb.ax.tick_params(labelsize=6.5); cb.set_label('K', fontsize=7)
finish(fig, FIG / 'imaging_mw_hi_channels.png')
print('ok')
