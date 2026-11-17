"""w4_distance: figure - Kim+2007 FIMS/SPEAR C IV map (their Fig. 1) with the Silk overlaid, and C IV cross-profiles.
Reads data/derived/w4_distance_kim07_digitize.json / _map.npz (from w4_distance_kim07_digitize.py).
Output: figures/w4_distance_kim07_civ.png
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_kim07_fig.py
"""
import json, sys
from pathlib import Path
import numpy as np
from PIL import Image
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *  # noqa

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
J = json.load(open(ROOT / 'data/derived/w4_distance_kim07_digitize.json'))
M = np.load(ROOT / 'data/derived/w4_distance_kim07_map.npz')
cx, cy = M['cx'], M['cy']
img = np.asarray(Image.open(ROOT / 'data/raw/w4_distance_kim07/kim07_fg1.jpg').convert('RGB'))


def design(u_, v_, deg=3):
    return np.vstack([u_ ** i * v_ ** j for i in range(deg + 1) for j in range(deg + 1 - i)]).T


def lb2xy(l, b):
    D = design((np.atleast_1d(l) - 200) / 20, (np.atleast_1d(b) - 12) / 20)
    return D @ cx, D @ cy


tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
fig, axes = new_sky_figure(13, 6.4, 'FIMS/SPEAR C IV and the Silk',
                           'Kim et al. 2007 Fig. 1 (C IV 1548/1550, 120′ smoothing, 0–21091 LU) registered to (l, b) · Silk strands from silk_trace.json',
                           ncols=2, gridspec_kw={'width_ratios': [1.05, 1]})
ax = axes[0]
ax.imshow(img[:560], origin='upper')
for key, col, lw in (('E', OIII_C, 1.6), ('W', '#b39ddb', 1.2), ('W2', '#80cbc4', 1.0)):
    pts = tr['spine'] if key == 'E' else tr['strands'][key]['points']
    c = SkyCoord([p['ra'] for p in pts] * u.deg, [p['dec'] for p in pts] * u.deg).galactic
    x, y = lb2xy(c.l.deg, c.b.deg)
    ax.plot(x, y, color=col, lw=lw, path_effects=[HALO], label={'E': 'Silk E strand (spine)', 'W': 'W strand', 'W2': 'W2 strand'}[key])
hp_ = tr['halpha_patch']
c = SkyCoord([p[0] if isinstance(p, list) else p['ra'] for p in hp_] * u.deg,
             [p[1] if isinstance(p, list) else p['dec'] for p in hp_] * u.deg).galactic
x, y = lb2xy(c.l.deg, c.b.deg)
ax.plot(np.r_[x, x[:1]], np.r_[y, y[:1]], color=HA_C, lw=1.0, label='Hα/[S II] patch')
pg = SkyCoord(115.2847 * u.deg, 20.8334 * u.deg).galactic
x, y = lb2xy(pg.l.deg, pg.b.deg)
ax.plot(x, y, marker='*', ms=11, color=SUN, mec='k', ls='none', label='PGCC G199.09+19.95')
# Dec 17.6 junction
sp = tr['spine']; dd = np.array([p['dec'] for p in sp]); rr = np.array([p['ra'] for p in sp])
cj = SkyCoord(np.interp(17.6, dd, rr) * u.deg, 17.6 * u.deg).galactic
x, y = lb2xy(cj.l.deg, cj.b.deg)
ax.plot(x, y, marker='o', ms=6, mfc='none', mec='w', ls='none')
label(ax, x[0] + 8, y[0] + 14, 'Dec 17.6', size=7.5)
ax.set_xlim(0, 600); ax.set_ylim(560, 0)
ax.set_xticks([]); ax.set_yticks([])
ax.legend(loc='lower left', fontsize=7.5, ncol=2, framealpha=0.85)

ax = axes[1]
for name, col, lab in (('south_dec_lt_17.6', OIII_C, 'south (Dec < 17.6)'), ('north_dec_ge_17.6', SII_C, 'north (Dec 17.6–27)')):
    s = J['sections'][name]
    o = np.array(s['offsets_deg']); v = np.array([np.nan if x is None else x for x in s['mean_civ_LU']])
    e = np.array([np.nan if x is None else x for x in s['err_LU']])
    if 'north' in name:
        bad = o > 3.0          # north tracks run off the FIMS map edge here
        v[bad] = np.nan; e[bad] = np.nan
    ax.errorbar(o, v / 1e3, yerr=e / 1e3, color=col, fmt='o', ms=3, lw=0.8, label=lab)
    f = s['fits']['fwhm_2.1']
    sg = 2.1 / 2.3548
    xx = np.linspace(-5, 5, 400)
    ax.plot(xx, (f['bkg_at0_LU'] + f['slope_LU_per_deg'] * xx + f['A_LU'] * np.exp(-0.5 * ((xx - f['mu_deg']) / sg) ** 2)) / 1e3,
            color=col, lw=1.2, alpha=0.9)
    label(ax, -4.8 if 'south' in name else 0.6, (f['bkg_at0_LU'] + f['A_LU']) / 1e3 + (0.5 if 'south' in name else -3.3),
          f"{'S' if 'south' in name else 'N'}: A = {f['A_LU']/1e3:.1f}±{f['e_A']/1e3:.1f}k LU, "
          f"μ = {f['mu_deg']:+.2f}°", color=col, size=8)
ax.axvline(0, color=FGDIM, lw=0.8, ls=':')
ax.set_xlabel('offset from [O III] spine (deg, + = east)')
ax.set_ylabel('C IV (10³ LU, digitised)')
ax.grid(alpha=0.3)
ax.legend(loc='upper right', fontsize=8)
label(ax, -4.8, 2.6, 'fit: linear + Gaussian, FWHM fixed 2.1° (map beam)', size=7.5, color=FGDIM)
ax.set_ylim(2, 13)
finish(fig, ROOT / 'figures/w4_distance_kim07_civ.png')
print('ok')
