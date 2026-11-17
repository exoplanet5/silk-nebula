"""w4_distance: figures for the dust/clump distance and the target list + proper-motion plan.
Outputs: figures/w4_distance_dust_clump.png, figures/w4_distance_targets_pm.png
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_figs.py
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *  # noqa
from scipy.interpolate import griddata

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
SP = np.array([[p['ra'], p['dec']] for p in tr['spine']])
WS = np.array([[p['ra'], p['dec']] for p in tr['strands']['W']['points']])
W2 = np.array([[p['ra'], p['dec']] for p in tr['strands']['W2']['points']])
HP = np.array(tr['halpha_patch'])
PG = (115.2847, 20.8334)


def overlay(ax, lw=1.3):
    ax.plot(SP[:, 0], SP[:, 1], color=OIII_C, lw=lw, path_effects=[HALO])
    ax.plot(WS[:, 0], WS[:, 1], color='#b39ddb', lw=lw * 0.8)
    ax.plot(W2[:, 0], W2[:, 1], color='#80cbc4', lw=lw * 0.7)
    ax.plot(np.r_[HP[:, 0], HP[:1, 0]], np.r_[HP[:, 1], HP[:1, 1]], color=HA_C, lw=0.9)
    ax.plot(*PG, marker='*', ms=12, color=SUN, mec='k', ls='none')


# ------------------------------------------------------------------ figure 1: dust + clump
D = np.load(ROOT / 'data/derived/w4_distance_dust_pixels.npz')
dust = json.load(open(ROOT / 'data/derived/w4_distance_dust.json'))
cs = json.load(open(ROOT / 'data/derived/w4_distance_clump_stars.json'))
fig, axes = new_sky_figure(14, 10, 'Dust toward the Silk and PGCC G199.09+19.95',
                           'Edenhofer+24 3D dust (range-read, nside 256) · SFD tomography of Gaia DR3 GSP-Phot and StarHorse stars',
                           nrows=2, ncols=2)
gx, gy = np.meshgrid(np.linspace(111.5, 121.5, 300), np.linspace(10.5, 29.5, 500))
for ax, key, vmax, ttl in ((axes[0, 0], 'E_320-400', 0.008, 'E(320–400 pc)'), (axes[0, 1], 'E_400-600', 0.018, 'E(400–600 pc)')):
    z = griddata((D['ra'], D['dec']), D[key], (gx, gy), method='linear')
    im = ax.imshow(z, origin='lower', extent=(111.5, 121.5, 10.5, 29.5), cmap='magma', vmin=0, vmax=vmax, aspect=1 / np.cos(np.radians(20)))
    overlay(ax)
    ax.set_xlim(121.5, 111.5); ax.set_ylim(10.5, 29.5)
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label('E (≈ E(B−V), mag)')
    label(ax, 121.2, 28.7, ttl, size=10)
label(axes[0, 0], 120.9, 11.2, 'east ←', size=8, color=FGDIM)
# panel c: Edenhofer on/off profile toward the clump
ax = axes[1, 0]
p = dust['pgcc']
r = np.array(p['profile_r_pc']); on = np.array(p['on_density']); off = np.array(p['off_density_median']); sd = np.array(p['off_density_sd'])
ax.fill_between(r, (off - sd) * 1e5, (off + sd) * 1e5, color=FGDIM, alpha=0.25, lw=0, label='off 0.6–1.5° (median ± sd)')
ax.plot(r, off * 1e5, color=FGDIM, lw=1.2)
ax.plot(r, on * 1e5, color=SUN, lw=1.6, label='on clump (r < 0.2°)')
for d0, c0, t0 in ((160, OK, '160'), (282, OIII_C, '282'), (1000, DANGER, '')):
    ax.axvline(d0, color=c0, lw=0.9, ls=':')
ax.axvspan(252, 318, color=OIII_C, alpha=0.10, lw=0)
ax.set_xscale('log'); ax.set_xlim(90, 750)
ax.set_xlabel('distance (pc)'); ax.set_ylabel('dust density (10⁻⁵ E pc⁻¹)')
ax.legend(loc='upper left', fontsize=7.5)
ax.set_ylim(0, max(on.max(), (off + sd).max()) * 1e5 * 1.15)
label(ax, 95, ax.get_ylim()[1] * 0.55,
      f"3D map on−off to 750 pc: {p['on_minus_off_total']:.4f} ± {p['E_off_total_sd']:.4f}\nSFD clump excess: {dust['column_2d']['SFD_Dust_Map']['excess']:.3f}  → map does not resolve the clump",
      size=7.5, color=FG)
ax.grid(alpha=0.3)
# panel d: SFD tomography slopes
ax = axes[1, 1]
for key, col, lab, dx in (('starhorse_AV_r<0.4', SUN, 'StarHorse A_V (r < 0.4°)', 1.0), ('gaia_A0_r<0.5', ACCENT, 'Gaia GSP-Phot A₀ (r < 0.5°)', 1.04)):
    rows = [q for q in cs['sfd_tomography'][key] if 'slope' in q]
    dm = np.array([np.sqrt(q['d_lo'] * q['d_hi']) for q in rows]) * dx
    ax.errorbar(dm, [q['slope'] for q in rows], yerr=[q['slope_err'] for q in rows],
                xerr=[[m - q['d_lo'] for m, q in zip(dm, rows)], [q['d_hi'] - m for m, q in zip(dm, rows)]],
                fmt='o', color=col, ms=5, lw=1, capsize=0, label=lab)
    for m, q in zip(dm, rows):
        if q['slope'] + q['slope_err'] + 0.4 < 11:
            label(ax, m, q['slope'] + q['slope_err'] + 0.4, f"n={q['n']}", size=6.5, color=col, ha='center')
ax.axhspan(2.8, 3.2, color=OK, alpha=0.18, lw=0)
label(ax, 95, 3.35, 'all clump dust in front (A_V/E ≈ 2.8–3.2)', size=7.5, color=OK)
ax.axhline(0, color=FGDIM, lw=0.8)
label(ax, 95, 0.25, 'clump behind the stars', size=7.5, color=FGDIM)
ax.axvspan(252, 318, color=OIII_C, alpha=0.10, lw=0)
ax.axvline(160, color=OK, lw=0.9, ls=':')
ax.set_xscale('log'); ax.set_xlim(90, 1600); ax.set_ylim(-10, 12)
ax.set_xlabel('stellar distance (pc)'); ax.set_ylabel('dA / dE(B−V)_SFD')
ax.legend(loc='lower right', fontsize=7.5)
ax.grid(alpha=0.3)
fig.subplots_adjust(left=0.06, right=0.97, bottom=0.07, top=0.91, hspace=0.25, wspace=0.18)
finish(fig, ROOT / 'figures/w4_distance_dust_clump.png', tight=False)

# ------------------------------------------------------------------ figure 2: targets + proper motion
T = pd.read_csv(ROOT / 'data/derived/w4_distance_targets.csv')
pm = json.load(open(ROOT / 'data/derived/w4_distance_pm.json'))
fig, axes = new_sky_figure(14, 8.4, 'Silk distance bracket: UV targets and proper motion',
                           'Gaia DR3 blue stars ≤ 0.5° from the strands (plx/err > 5) with GALEX FUV · expected filament motion vs distance',
                           ncols=2, gridspec_kw={'width_ratios': [1.0, 1.15]})
ax = axes[0]
ax.plot(SP[:, 1], np.zeros(len(SP)), color=OIII_C, lw=1.0)
cmap = {'WD': ACCENT, 'hot subdwarf': '#c792ea', 'early-type (B/A)': WARN, 'other blue': FGDIM}
for c, col in cmap.items():
    s = T[T['class'] == c]
    if not len(s):
        continue
    sz = np.where(np.isfinite(s.FUV), np.clip((20 - s.FUV) * 14, 8, 120), 8)
    face = np.where(s.cos_orbits_SN15 <= 5, col, 'none')
    ax.scatter(s.dist_pc, s.dec, s=sz, facecolors=face, edgecolors=col, lw=1.0, label=c)
for _, r in T[(T.cos_orbits_SN15 <= 5) & (T.sep_min <= 0.3) & (T['class'] == 'WD')].iterrows():
    label(ax, r.dist_pc * 1.04, r.dec + 0.15, f"FUV {r.FUV:.1f}, {r.sep_min:.2f}°", size=6.5)
ax.axvspan(252, 318, color=OIII_C, alpha=0.12, lw=0)
ax.axvline(160, color=OK, lw=0.9, ls=':')
ax.axhline(17.6, color=FGDIM, lw=0.8, ls='--')
label(ax, 104, 17.25, 'Dec 17.6 (south below)', size=7, color=FGDIM)
ax.set_xscale('log'); ax.set_xlim(100, 1000); ax.set_ylim(11.3, 29)
ax.set_xlabel('Gaia distance 1/ϖ (pc)'); ax.set_ylabel('Dec (deg)')
ax.legend(loc='upper left', fontsize=7.5, ncol=2, title='filled: COS S/N 15 at C IV in ≤ 5 orbits; size ∝ FUV', title_fontsize=7)
ax.grid(alpha=0.3)
ax = axes[1]
d = np.logspace(np.log10(100), np.log10(1600), 200)
for v, ls in ((100, '--'), (120, '-'), (150, '--')):
    ax.plot(d, v / (4.74047 * d / 1000), color=OIII_C, ls=ls, lw=1.4 if v == 120 else 0.9, label=f'v_s = {v} km/s')
ax.fill_between(d, 100 / (4.74047 * d / 1000), 150 / (4.74047 * d / 1000), color=OIII_C, alpha=0.12, lw=0)
refl = [g['reflex_LSR_gas_mas_yr']['282']['normal_east'] for g in pm['geometry_every_4th_spine_point']]
ax.plot(d, np.abs(np.median(refl)) * 282 / d, color=WARN, lw=1.0, label='|reflex| of LSR-static gas (normal comp.)')
# detection thresholds (5 sigma, two epochs, 12-yr baseline) for per-epoch precision 0.1 / 0.2 / 0.5 arcsec
for sp_, col in ((0.1, OK), (0.2, SUN), (0.5, DANGER)):
    thr = 5 * np.sqrt(2) * sp_ / 12 * 1000
    ax.axhline(thr, color=col, lw=0.9, ls=':')
    label(ax, 105, thr * 1.06, f'5σ in 12 yr, σ_epoch = {sp_}″', size=7, color=col)
for d0 in (160, 282, 1000):
    ax.axvline(d0, color=FGDIM, lw=0.7, ls=':')
    label(ax, d0 * 1.03, 260, f'{d0} pc\n{120 / (4.74047 * d0 / 1000):.0f} mas/yr', size=7.5)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(100, 1600); ax.set_ylim(3, 400)
ax.set_xlabel('distance (pc)'); ax.set_ylabel('proper motion normal to the filament (mas/yr)')
ax.legend(loc='lower left', fontsize=7.5)
ax.grid(alpha=0.3, which='both')
finish(fig, ROOT / 'figures/w4_distance_targets_pm.png')
print('ok')
