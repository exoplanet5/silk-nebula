#!/usr/bin/env python
"""spectra_ism: figures for the WHAM-SS DR1 analysis (run after spectra_ism_wham.py).

figures/spectra_ism_wham_maps.png     total H-alpha, main-component V_LSR, IV (-60..-20) H-alpha, NSNS [O III] beam mean
figures/spectra_ism_wham_spectra.png  WHAM spectra at P1-P4, on-off difference spectrum, cross-filament profiles
"""
import csv
import json
import sys

import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FG, FGDIM, OIII_C, HA_C, SII_C, NII_C, WARN, OK  # noqa: E402
from spectra_ism_common import DER, FIG, REGIONS, draw_trace  # noqa: E402

rows = list(csv.DictReader(open(DER / 'spectra_ism_wham_beams.csv')))


def col(k):
    return np.array([float(r[k]) if r[k] != '' else np.nan for r in rows])


ra, dec = col('ra'), col('dec')
R = json.load(open(DER / 'spectra_ism_wham_results.json'))
S = np.load(DER / 'spectra_ism_wham_stacks.npz')
V = S['vgrid']


def beam_map(ax, val, cmap, norm, title):
    pats = [Ellipse((r, d), 0.98 / np.cos(np.deg2rad(d)), 0.98) for r, d in zip(ra, dec)]
    pc = PatchCollection(pats, cmap=cmap, norm=norm, edgecolor='none')
    pc.set_array(np.ma.masked_invalid(val))
    ax.add_collection(pc)
    ax.set_xlim(128, 104)
    ax.set_ylim(2, 34)
    ax.set_aspect(1 / np.cos(np.deg2rad(18)))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title(title, loc='left', color=FG)
    draw_trace(ax, color='#ffffff', lw=0.9, alpha=0.8, label=False)
    return pc


# ------------------------------------------------------------------------------------------
fig, axs = new_sky_figure(15, 6.6, 'Silk Nebula · WHAM-SS DR1 Hα', '1° beams; V_LSR; white = [O III] spine, dotted = envelope',
                          nrows=1, ncols=4)
c1v = col('c1_v')
iv = col('I_m60_m20')
ob = col('oiii_beam_R_nsns')
panels = [(col('I_tot_R'), 'magma', mcolors.LogNorm(1.0, 12), 'Total Hα (R)'),
          (c1v, 'coolwarm', mcolors.Normalize(-8, 8), 'Main comp. V_LSR (km/s)'),
          (iv, 'viridis', mcolors.Normalize(0, 0.6), 'Hα −60…−20 km/s (R)'),
          (ob, 'cividis', mcolors.Normalize(0, 1.2), 'NSNS [O III], beam mean (R)')]
for ax, (val, cm, nm, tt) in zip(axs, panels):
    pc = beam_map(ax, val, cm, nm, tt)
    cb = fig.colorbar(pc, ax=ax, orientation='horizontal', fraction=0.05, pad=0.08)
    cb.ax.tick_params(labelsize=7, colors=FGDIM)
finish(fig, FIG / 'spectra_ism_wham_maps.png')

# ------------------------------------------------------------------------------------------
fig, axs = new_sky_figure(13.5, 8.2, 'Silk Nebula · WHAM Hα spectra vs LAMOST line velocities',
                          'V_LSR (standard solar motion); vertical lines = Wu+2026 Table 4 V_helio converted to LSR',
                          nrows=2, ncols=3)
axs = axs.ravel()
vl = R['vlsr_recomputed_from_vhelio']
FIP = 12.0
for ax, k in zip(axs[[0, 1, 2, 3]], ['P1', 'P2', 'P3', 'P4']):
    sp, er = S[f'spec_{k}'], S[f'err_{k}']
    ax.fill_between(V, sp - er, sp + er, color='#ffffff', alpha=0.15, lw=0)
    ax.plot(V, sp, color='#ffffff', lw=1.0, label='WHAM')
    nb = R['nearest_beams'][k][0]
    tot = np.zeros_like(V)
    for c in nb['components']:
        s = np.hypot(c['fwhm_int'], FIP) / 2.3548
        g = c['I_R'] / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((V - c['v']) / s) ** 2)
        tot += g
        ax.plot(V, g, ':', color=FGDIM, lw=0.8)
    ax.plot(V, tot, '-', color=OK, lw=0.8, alpha=0.8, ls='--', label='Gaussian fit')
    for ln, cc in (('OIII', OIII_C), ('Ha', HA_C), ('NII', NII_C), ('SII', SII_C)):
        if ln in vl[k]:
            ax.axvline(vl[k][ln], color=cc, lw=1.1, ls='-', alpha=0.9, label=f'LAMOST {ln}')
    ax.set_xlim(-80, 70)
    ax.set_ylim(-0.005, None)
    ax.set_xlabel('V_LSR (km/s)')
    ax.set_ylabel('R / (km/s)')
    c0 = nb['components'][0]
    ax.set_title(f"{k}: beam {nb['sep_deg']:.2f}° away, I = {nb['I_tot_R']:.2f} R", loc='left', color=FG, fontsize=9)
    label(ax, 0.02, 0.92, f"main V = {c0['v']:+.1f}±{c0['ev']:.1f}\nFWHM_int = {c0['fwhm_int']:.0f} km/s",
          transform=ax.transAxes, size=7.5, va='top', color=FG)
    if k == 'P1':
        ax.legend(loc='upper right', fontsize=6.5, ncol=1, frameon=True)

# difference spectrum
ax = axs[4]
d, e = S['diff_mean'], S['diff_err']
ax.fill_between(V, d - e, d + e, color=OK, alpha=0.25, lw=0)
ax.plot(V, d, color=OK, lw=1.1, label='on − off (mean of pairs)')
st = R['onoff_stack']
dO = st['delta_oiii_beam_R_nsns'][0]
vo = R['v_oiii_lsr_mean_P1P3']
s = st['fixed_v%.1f_fwhm_ha_from_same_gas_T1e4' % vo]['fwhm_obs'] / 2.3548
for ratio, ls in ((1.0, '--'), (0.3, ':')):
    g = ratio * dO / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((V - vo) / s) ** 2)
    ax.plot(V, g, ls, color=OIII_C, lw=1.0, label=f'Hα/[O III] = {ratio:g} at [O III] V')
ax.axhline(0, color=FGDIM, lw=0.5)
ax.axvline(vo, color=OIII_C, lw=0.6, alpha=0.6)
ax.set_xlim(-80, 70)
ax.set_xlabel('V_LSR (km/s)')
ax.set_ylabel('R / (km/s)')
ax.set_title(f"On (|d|<{st['on_max_deg']}°) − off (1–3° E+W), {st['n_pairs']} pairs", loc='left', color=FG, fontsize=9)
fx = st['fixed_v%.1f_fwhm_ha_from_same_gas_T1e4' % vo]
label(ax, 0.02, 0.97, f"Hα at {vo:+.1f} km/s:\n{fx['I_R']:.3f}±{fx['eI_R_boot']:.3f} R\nΔ[O III] beam\n= {dO:.2f} R (NSNS)",
      transform=ax.transAxes, size=7.5, va='top', color=FG)
ax.legend(loc='upper right', fontsize=6.5)

# cross-filament profiles
ax = axs[5]
P = R['cross_profiles']
edges = np.array(P['bin_edges_deg'])
xc = 0.5 * (edges[1:] + edges[:-1])
for nm, cc, lab in (('m20_m5', NII_C, '−20…−5'), ('m5_p9', HA_C, '−5…+9'), ('p9_p20', OIII_C, '+9…+20 ([O III] V)'),
                    ('m60_m20', FGDIM, '−60…−20')):
    ax.errorbar(xc + {'m20_m5': -0.06, 'm5_p9': -0.02, 'p9_p20': 0.02, 'm60_m20': 0.06}[nm], P[nm]['mean_R'],
                yerr=P[nm]['sem_R'], fmt='o-', ms=3, lw=0.9, color=cc, label=f'Hα {lab} km/s', capsize=0)
ax2 = ax.twinx()
ax2.plot(xc, P['oiii_beam_nsns']['mean_R'], '-', color='#ffffff', lw=1.6, alpha=0.5)
ax2.set_ylabel('NSNS [O III] beam mean (R)', color=FGDIM)
ax2.tick_params(colors=FGDIM)
ax.axhline(0, color=FGDIM, lw=0.5)
ax.set_xlabel('offset from spine (deg, + = east)')
ax.set_ylabel('Hα beam − local background (R)')
ax.set_title('Cross-filament profiles (Dec 12–28)', loc='left', color=FG, fontsize=9)
ax.legend(loc='upper left', fontsize=6.3)
finish(fig, FIG / 'spectra_ism_wham_spectra.png')
print('done')
