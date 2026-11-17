"""Figure: the Silk filament at the edge of the Monogem Ring soft X-ray emission.

(a) eRASS1 0.2-0.6 keV (cleaned, smoothed) wide view with the Silk spine and PSR B0656+14, sector boundaries;
(b) radial profiles around PSR B0656+14 in the PA 0-20 and 20-40 deg sectors (which contain the southern Silk)
    for eRASS1 0.2-0.6, ROSAT R12, R45, with the Silk radius marked;
(c) along-spine X-ray half-drop edge and HI-wall offsets per 1.5-deg segment.
Reads data/derived/imaging_mw_edges.json. Output figures/imaging_mw_monogem_edge.png
"""
import json
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FGDIM, OIII_C, WARN, OK, SUN, HA_C, ACCENT
from imaging_mw_common import DER, FIG, load_grid, GAL_EXTENT, gal_axes, get_spine, spine_gal, MONOGEM_PULSAR, gal
from imaging_mw_fig_panels import smooth_nan

E = json.loads((DER / 'imaging_mw_edges.json').read_text())
sp = get_spine(); ls_, bs_ = spine_gal(sp)
pl, pb = E['monogem_ring']['centre']['l'], E['monogem_ring']['centre']['b']
fig = None
fig, axes = new_sky_figure(17.0, 6.6, 'Silk Nebula · Monogem Ring X-ray edge',
                           'soft X-rays fill the region inside (west of) the [O III] spine and fade at the filament; HI at the [O III] velocity rises just outside',
                           nrows=1, ncols=3, gridspec_kw={'width_ratios': [1.15, 1, 1]})
ax = axes[0]
a = load_grid('ero_0206_clean', 'gal')
v0, v1 = np.nanpercentile(a, [2, 98])
ax.imshow(a, origin='lower', extent=GAL_EXTENT, vmin=v0, vmax=v1, cmap='inferno', interpolation='nearest')
ax.plot(ls_, bs_, color='#ffffff', lw=1.3)
ax.plot(pl, pb, '+', color=SUN, ms=10, mew=1.5)
c0 = SkyCoord(l=pl * u.deg, b=pb * u.deg, frame='galactic')
for pa in [0, 20, 40]:
    rr = np.linspace(0, 16, 50)
    pts = c0.directional_offset_by(pa * u.deg, rr * u.deg)
    ax.plot(pts.l.deg, pts.b.deg, color=ACCENT, lw=0.7, ls='--')
for r in [10.0, 11.4]:
    ring = c0.directional_offset_by(np.linspace(0, 360, 361) * u.deg, r * u.deg)
    ax.plot(ring.l.deg, ring.b.deg, color=ACCENT if r > 10.5 else FGDIM, lw=0.7, ls=':')
label(ax, pl - 0.5, pb - 1.8, 'PSR B0656+14', color=SUN, size=7.5)
label(ax, ls_[0] + 0.5, bs_[0] - 1.8, 'Silk', color='#ffffff', size=8)
gal_axes(ax, 'eRASS1 0.2–0.6 keV; PA 0–40° sectors; r = 10°, 11.4°')

ax = axes[1]
cols = {'ero_0206_clean': OIII_C, 'rass_r12': WARN, 'rass_r45': OK}
labs = {'ero_0206_clean': 'eRASS1 0.2–0.6 keV', 'rass_r12': 'ROSAT R12', 'rass_r45': 'ROSAT R45'}
for band in ['ero_0206_clean', 'rass_r12', 'rass_r45']:
    for rec, lsty in zip([r for r in E['monogem_ring'][band] if r['pa0'] in (0, 20)], ['-', '--']):
        r = np.array(rec['profile_r']); p = np.array([np.nan if x is None else x for x in rec['profile']])
        norm = np.nanmedian(p[(r > 14) & (r < 17)])
        ax.plot(r, p / norm, color=cols[band], ls=lsty, lw=1.3, label=f"{labs[band]} PA {rec['pa0']}–{rec['pa0'] + 20}°")
for rec in [r for r in E['monogem_ring']['ero_0206_clean'] if r['pa0'] in (0, 20)]:
    ax.axvline(rec['silk_r_in_sector'], color='#ffffff', lw=1.0, ls='-' if rec['pa0'] == 0 else '--')
label(ax, 11.7, 0.62, 'Silk', color='#ffffff', size=8)
ax.set_xlim(2, 18); ax.set_ylim(0.5, 2.2)
ax.set_xlabel('angular distance from PSR B0656+14 (deg)'); ax.set_ylabel('brightness / median at r = 14–17°')
ax.set_title('radial profiles in the Silk sectors', loc='left', fontsize=9.5)
ax.grid(alpha=0.3); ax.legend(fontsize=6.5, loc='upper left', ncol=1)

ax = axes[2]
for band, col, key, lab_ in [('ero_0206_clean', OIII_C, 'edge_off_deg', 'eRASS1 0.2–0.6 keV edge'),
                             ('rass_r12', WARN, 'edge_off_deg', 'ROSAT R12 edge'),
                             ('galfa_nhi_oiii', OK, 'wall_off_deg', 'HI wall +11…+16 km/s')]:
    segs = E['along_spine'][band]['segments']
    d = [s_['dec_mid'] for s_ in segs]; o = [np.nan if s_[key] is None else s_[key] for s_ in segs]
    ax.plot(o, d, 'o-', color=col, ms=4, lw=1.0, label=lab_)
ax.axvline(0, color='#ffffff', lw=1.2)
ax.axvspan(-0.138, 0.138, color=OIII_C, alpha=0.15, lw=0)
label(ax, 0.08, 27.6, '[O III] spine', color='#ffffff', size=7.5)
ax.set_xlim(-3.5, 3); ax.set_xlabel('offset from spine (deg, + east = outside)'); ax.set_ylabel('Dec of segment centre (deg)')
ax.set_title('edge positions along the spine (1.5° segments)', loc='left', fontsize=9.5)
ax.grid(alpha=0.3); ax.legend(fontsize=7, loc='lower left')
finish(fig, FIG / 'imaging_mw_monogem_edge.png')
print('ok')
