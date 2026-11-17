"""verify_hiwham: figure for HI along the spine, PGCC / H-alpha patch spectra and the catalogued HI shells.
Reads data/derived/verify_hiwham_hi_extra.{json,npz}; writes figures/verify_hiwham_hi_extra.png
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, WARN, OK, ACCENT, FG
from verify_hiwham_common import DER, FIG, vlsr_minus_vhelio, strand

J = json.loads((DER / 'verify_hiwham_hi_extra.json').read_text())
Z = np.load(DER / 'verify_hiwham_hi_extra.npz')
v = Z['A_v']; dec = Z['A_dec']; on = Z['A_on']; off = Z['A_off']
ra_E, dec_E = strand('E')
fig, axes = new_sky_figure(15, 10.5, 'Silk Nebula · HI velocity structure along the spine, at the Hα patch and in the catalogued shells',
                           'GALFA-HI DR2 (2′ mosaic); [O III] V_helio +13.5 (E) and +18.5 (W) converted to LSR per position',
                           nrows=2, ncols=3)
ax = axes.ravel()
vs = (v > -25) & (v < 35)
ext = [v[vs][0], v[vs][-1], dec[0] - 0.25, dec[-1] + 0.25]
im = ax[0].imshow(on[:, vs], origin='lower', aspect='auto', extent=ext, cmap='viridis', vmin=0, vmax=np.nanpercentile(on[:, vs], 99.5))
vo = 13.5 + vlsr_minus_vhelio(np.interp(dec, dec_E, ra_E), dec)
vw = 18.5 + vlsr_minus_vhelio(np.interp(dec, dec_E, ra_E) - 0.4, dec)
ax[0].plot(vo, dec, color=OIII_C, lw=1.6, label='[O III] E (V_helio +13.5)')
ax[0].plot(vw, dec, color=OK, lw=1.2, ls='--', label='[O III] W (V_helio +18.5)')
ax[0].plot(J['A_along_spine']['m1_on_kms'], dec, color='#ffffff', lw=1.0, marker='.', ms=3, label='HI moment-1 (on spine)')
ax[0].axvspan(11.2, 15.5, color=FGDIM, alpha=0.15, lw=0)
label(ax[0], 11.4, dec[-1] - 0.3, 'paper', color=FGDIM, size=7)
ax[0].set_xlabel('V_LSR (km/s)'); ax[0].set_ylabel('spine Dec (deg)')
ax[0].set_title('(a) on-spine HI T_B (|d|<0.15°)', loc='left')
ax[0].legend(fontsize=6.5, loc='lower right')
cb = fig.colorbar(im, ax=ax[0], fraction=0.04, pad=0.01); cb.set_label('K', fontsize=8)
dd = (on - off)[:, vs]
im = ax[1].imshow(dd, origin='lower', aspect='auto', extent=ext, cmap='RdBu_r', vmin=-3, vmax=3)
ax[1].plot(vo, dec, color=OIII_C, lw=1.6); ax[1].plot(vw, dec, color=OK, lw=1.2, ls='--')
ax[1].set_xlabel('V_LSR (km/s)'); ax[1].set_ylabel('spine Dec (deg)')
ax[1].set_title('(b) on − off (0.4–1.2° both sides) HI T_B', loc='left')
cb = fig.colorbar(im, ax=ax[1], fraction=0.04, pad=0.01); cb.set_label('K', fontsize=8)
# (c) PGCC / patch spectra
for key, col, lab in (('C_PGCC_r<0.1', ACCENT, 'PGCC G199.09+19.95 r<0.1°'), ('C_patch_polygon', HA_C, 'Hα patch polygon')):
    ax[2].plot(v, Z[f'{key}_on'], color=col, lw=1.5, label=lab + ' (on)')
    ax[2].plot(v, Z[f'{key}_off'], color=col, lw=0.9, ls=':', label='annulus (off)')
vW_patch = 18.5 + vlsr_minus_vhelio(115.9, 19.6)[0]
vE_patch = 13.5 + vlsr_minus_vhelio(115.9, 19.6)[0]
ax[2].axvspan(vW_patch - 1.5, vW_patch + 1.5, color=OK, alpha=0.2, lw=0)
ax[2].axvline(vE_patch, color=OIII_C, lw=1.0)
label(ax[2], vW_patch + 1.6, 1.0, 'W strand', color=OK, size=7)
label(ax[2], vE_patch - 9, 1.0, 'E', color=OIII_C, size=7)
ax[2].set_xlim(-25, 35); ax[2].set_xlabel('V_LSR (km/s)'); ax[2].set_ylabel('T_B (K)')
fit = J['C_pgcc_patch']['PGCC_r<0.1']['fit_narrow']
ax[2].set_title('(c) PGCC: narrow HI at −0.7 and +6.5 km/s', loc='left')
ax[2].legend(fontsize=6.5, loc='upper left')
# (d, e) shells velocity-radius diagrams
red = Z['D_redges']; rc = 0.5 * (red[1:] + red[:-1])
for k, (nm, ttl) in enumerate((('GS202+17+019', 'GS202+17+019 (Suad+14), v0 +18.6'),
                               ('GS202+19-019_sizes_as_radii', 'GS202+19−019 (E&P05, Δl,Δb as radii), v −18.8'))):
    vr = Z[f'D_{nm}_vr'][:, vs]
    rel = vr - np.nanmean(vr, axis=0, keepdims=True)
    im = ax[3 + k].imshow(rel.T, origin='lower', aspect='auto', extent=[red[0], red[-1], v[vs][0], v[vs][-1]],
                          cmap='RdBu_r', vmin=-2, vmax=2)
    sh = J['D_shells'][nm]['params']
    ax[3 + k].axhline(sh['v0'], color=WARN, lw=1.2)
    ax[3 + k].axhspan(sh['v0'] - sh['dv'] / 2, sh['v0'] + sh['dv'] / 2, color=WARN, alpha=0.08, lw=0)
    ax[3 + k].axhspan(0.3, 4.2, color=OIII_C, alpha=0.25, lw=0)
    r = J['D_shells'][nm]['spine_rho']
    ax[3 + k].axvline(r['min'], color=FGDIM, lw=1.2, ls='--'); ax[3 + k].axvline(r['median'], color=FGDIM, lw=1.2, ls=':')
    label(ax[3 + k], 0.03, 5.2, 'revised [O III]', color=OIII_C, size=7)
    label(ax[3 + k], 0.03, sh['v0'] + 1.0, 'catalogue v', color=WARN, size=7)
    ax[3 + k].set_ylim(-25, 35)
    ax[3 + k].set_xlabel('elliptical radius ρ; grey: spine min (--) and median (:)'); ax[3 + k].set_ylabel('V_LSR (km/s)')
    ax[3 + k].set_title(f'({"de"[k]}) {ttl}', loc='left', fontsize=9)
    cb = fig.colorbar(im, ax=ax[3 + k], fraction=0.04, pad=0.01); cb.set_label('ΔT_B vs ρ-mean (K)', fontsize=7)
# (f) LAB rim-minus-centre z per slice
for nm, col in (('GS202+17+019', WARN), ('GS202+19-019_sizes_as_diameters', ACCENT), ('GS202+19-019_sizes_as_radii', OK)):
    rr = J['D_shells'][nm]['LAB_rim_minus_centre']
    xs = [int(k.split('..')[0]) + 5 for k in rr]
    ax[5].plot(xs, [rr[k]['z_robust'] for k in rr], marker='o', color=col, lw=1.4, label=nm.replace('_sizes_as_', ' Δ as '))
ax[5].axhspan(-2, 2, color='#2a323c', alpha=0.5, lw=0)
ax[5].axvspan(0.3, 4.2, color=OIII_C, alpha=0.25, lw=0)
ax[5].set_xlabel('LAB slice centre V_LSR (km/s, 10 km/s slices)'); ax[5].set_ylabel('rim − centre, z vs l-shifted ellipses')
ax[5].set_title('(f) LAB shell cavity significance', loc='left')
ax[5].legend(fontsize=6.5, loc='upper left')
finish(fig, FIG / 'verify_hiwham_hi_extra.png')
print('saved')
