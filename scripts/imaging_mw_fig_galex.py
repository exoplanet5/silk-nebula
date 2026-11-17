"""GALEX FUV test of radiative-shock excitation: predicted FUV line excess vs measured on-filament contrast.

Prediction (all assumptions stored in data/derived/imaging_mw_fuv_prediction.json):
  I(5007) = 10-13 R (Wu+26 Table 3 converted; notes/paper_review / data/derived/paper_checks.json)
  Raymond et al. (1988) incomplete-shock models for 120-140 km/s, as quoted by Kim et al. (2010, ApJ 722, 388):
     O III] 1661,1666 / [O III] 4959,5007 ~ 0.3 (energy), C IV 1549 / O III] ~ 10-13, He II 1640 / O III] ~ 0.4-0.5
  [O III] 4959+5007 = 1.338 x 5007 (energy); photons = energy x lambda / lambda
  GALEX FUV: CU-equivalent of a line = I_line[ph cm-2 s-1 sr-1] x r(lambda) / W_eff, W_eff = 300 A (+-30%),
     relative response r(1549)=1.0, r(1640)=0.65, r(1663)=0.55 (approximate, read from the published GALEX FUV
     effective-area curve, Morrissey et al. 2007) ; extinction A_FUV = 8.0 E(B-V), E(B-V)=0.04 (SFD along spine; upper bound)
  C IV resonance-scattering suppression factors 1, 3.5 (Sankrit & Wood 2001 typical), 20 (Kim+10 region R2)
Output: figures/imaging_mw_galex.png
"""
import json
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FGDIM, OIII_C, WARN, OK, DANGER, ACCENT
from imaging_mw_common import DER, FIG, overlay_eq, get_spine

R2PH = 1e6 / (4 * np.pi)            # 1 Rayleigh in photons cm-2 s-1 sr-1
A = dict(I5007_R=(10.0, 13.0), r_OIII1663_to_OIII4959_5007_energy=0.3, r_CIV_to_OIII1663=(10.0, 13.0),
         r_HeII_to_OIII1663=(0.4, 0.5), f4959=1.338, W_eff_A=300.0, W_eff_frac_err=0.3,
         resp={'1549': 1.0, '1640': 0.65, '1663': 0.55}, EBV=0.04, AFUV_per_EBV=8.0, civ_suppression=[1.0, 3.5, 20.0])


def predict(i5007_R, civ_supp, ext=True):
    e5007 = i5007_R * R2PH * (6.626e-27 * 2.998e18 / 5007.0)          # erg-ish scale (constant cancels)
    eO3 = A['r_OIII1663_to_OIII4959_5007_energy'] * A['f4959'] * e5007
    nO3 = eO3 / (6.626e-27 * 2.998e18 / 1663.0)
    out = {}
    for tag, lo_hi in [('lo', 0), ('hi', 1)]:
        nC4 = A['r_CIV_to_OIII1663'][lo_hi] * eO3 / (6.626e-27 * 2.998e18 / 1549.0) / civ_supp
        nHe = A['r_HeII_to_OIII1663'][lo_hi] * eO3 / (6.626e-27 * 2.998e18 / 1640.0)
        cu = (nO3 * A['resp']['1663'] + nC4 * A['resp']['1549'] + nHe * A['resp']['1640']) / A['W_eff_A']
        cu_noC = (nO3 * A['resp']['1663'] + nHe * A['resp']['1640']) / A['W_eff_A']
        f = 10 ** (-0.4 * A['AFUV_per_EBV'] * A['EBV']) if ext else 1.0
        out[tag] = (cu * f, cu_noC * f)
    return out


pred = {}
for i5007 in A['I5007_R']:
    for s in A['civ_suppression']:
        p = predict(i5007, s)
        pred[f'I5007={i5007:.0f}R_CIVsupp={s:g}'] = {'total_CU_lo': p['lo'][0], 'total_CU_hi': p['hi'][0],
                                                     'OIII1663+HeII_only_CU_lo': p['lo'][1], 'OIII1663+HeII_only_CU_hi': p['hi'][1]}
T = json.loads((DER / 'imaging_mw_galex_tests.json').read_text())
P = json.loads((DER / 'imaging_mw_profiles.json').read_text())
G = json.loads((DER / 'imaging_mw_galex_profiles.json').read_text())
meas = {
    'tile_local_all': (T['segments']['FUV:all']['median_contrast_CU'], T['segments']['FUV:all']['null_rms_CU']),
    'tile_local_south': (T['segments']['FUV:south_dec12-18']['median_contrast_CU'], T['segments']['FUV:south_dec12-18']['null_rms_CU']),
    'murthy_all': (P['bands']['murthy_fuv']['segments']['all']['contrast'], P['bands']['murthy_fuv']['segments']['all']['contrast_null_rms']),
    'murthy_south': (P['bands']['murthy_fuv']['segments']['south']['contrast'], P['bands']['murthy_fuv']['segments']['south']['contrast_null_rms']),
}
(DER / 'imaging_mw_fuv_prediction.json').write_text(json.dumps({'assumptions': A, 'predictions_CU': pred,
                                                                'measured_CU_value_nullrms': meas, 'doc': __doc__}, indent=1))
for k, v in pred.items():
    print(k, {kk: round(vv) for kk, vv in v.items()})
print(meas)

# ------------------------------------------------------------------ figure
sp = get_spine()
fig, axes = new_sky_figure(16.5, 11.0, 'Silk Nebula · GALEX FUV/NUV test',
                           'pipeline skybg coadds (1′ blocks, tile offsets matched); measured on-filament excess vs radiative-shock expectation',
                           nrows=2, ncols=3)
for ax, b in zip(axes[0, :2], ['fuv', 'nuv']):
    a = fits.getdata(DER / 'imaging_mw_grids' / f'imaging_mw_galex_{b}_mosaic_1arcmin.fits').astype(float)
    m = np.isfinite(a); w = gaussian_filter(m.astype(float), 2.0)
    s = gaussian_filter(np.where(m, a, 0), 2.0) / np.where(w > 0, w, np.nan); s[w < 0.3] = np.nan
    v0, v1 = np.nanpercentile(s, [3, 97])
    im = ax.imshow(s, origin='lower', extent=[125, 108, 4, 30], vmin=v0, vmax=v1, cmap='magma', interpolation='nearest')
    overlay_eq(ax, sp, color='#7fdcff', lw=0.8, label_regions=(b == 'fuv'), lamost=False, halfwidth=False)
    ax.set_xlim(120.5, 112.5); ax.set_ylim(10.5, 23.5); ax.set_aspect(1 / np.cos(np.radians(17)))
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    ax.set_title(f'GALEX {b.upper()} skybg mosaic', loc='left', fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01); cb.set_label('CU', fontsize=7); cb.ax.tick_params(labelsize=6.5)

ax = axes[0, 2]
for i, seg in enumerate(['all', 'south_dec12-18', 'north_dec18-28']):
    for j, band in enumerate(['FUV', 'NUV']):
        r = T['segments'][f'{band}:{seg}']
        x0 = i * 3 + j
        ax.scatter(np.full(len(r['tile_contrasts_CU']), x0) + np.random.default_rng(1).uniform(-0.15, 0.15, len(r['tile_contrasts_CU'])),
                   r['tile_contrasts_CU'], s=9, color=OIII_C if band == 'FUV' else WARN, alpha=0.7, lw=0)
        ax.errorbar(x0 + 0.3, r['median_contrast_CU'], yerr=r['null_rms_CU'], fmt='o', color='#ffffff', ms=4, capsize=2)
        ax.scatter(np.full(len(r['null_medians_CU']), x0 + 0.55), r['null_medians_CU'], s=7, marker='x', color=FGDIM)
ax.axhline(0, color=FGDIM, lw=0.7)
ax.set_xticks([0.5, 3.5, 6.5]); ax.set_xticklabels(['whole spine', 'Dec 12–18', 'Dec 18–28'])
ax.set_ylabel('on − local ref (CU)')
ax.set_title('per-tile contrast; × = null shifts', loc='left', fontsize=9.5)
label(ax, 0.02, 0.97, 'FUV', color=OIII_C, size=8, transform=ax.transAxes, va='top')
label(ax, 0.10, 0.97, 'NUV', color=WARN, size=8, transform=ax.transAxes, va='top')
ax.grid(alpha=0.3)

ax = axes[1, 0]
for band, col in [('FUV', OIII_C), ('NUV', WARN)]:
    for seg, ls in [('south_dec12-18', '-'), ('north_dec18-28', '--')]:
        g = G['profiles'][f'{band}:{seg}']
        x = np.array(G['centres']); y = np.array(g['profile_CU'], float)
        ax.plot(x, y, color=col, ls=ls, lw=1.3, label=f'{band} {seg.split("_")[0]}')
ax.axvspan(-0.138, 0.138, color=OIII_C, alpha=0.12, lw=0); ax.axhline(0, color=FGDIM, lw=0.6)
ax.set_xlim(-3, 3); ax.set_xlabel('offset from spine (deg, + east)'); ax.set_ylabel('CU')
ax.set_title('joint-fit profile (tile offsets removed)', loc='left', fontsize=9.5); ax.grid(alpha=0.3)
ax.legend(fontsize=7, ncol=2, loc='upper left')

ax = axes[1, 1]
labels, lo, hi, cols = [], [], [], []
for s in A['civ_suppression']:
    p13 = pred[f'I5007=13R_CIVsupp={s:g}']; p10 = pred[f'I5007=10R_CIVsupp={s:g}']
    labels.append(f'shock, C IV/{s:g}'); lo.append(p10['total_CU_lo']); hi.append(p13['total_CU_hi']); cols.append(DANGER)
labels.append('shock, O III]+He II only'); lo.append(pred['I5007=10R_CIVsupp=1']['OIII1663+HeII_only_CU_lo'])
hi.append(pred['I5007=13R_CIVsupp=1']['OIII1663+HeII_only_CU_hi']); cols.append(WARN)
for k, nm in [('tile_local_all', 'measured: tiles, whole'), ('tile_local_south', 'measured: tiles, south'),
              ('murthy_all', 'measured: Murthy, whole'), ('murthy_south', 'measured: Murthy, south')]:
    v, e = meas[k]
    labels.append(nm); lo.append(max(v - 2 * e, 1.0)); hi.append(v + 2 * e); cols.append(OK)
y = np.arange(len(labels))[::-1]
for yy, a0, a1, c in zip(y, lo, hi, cols):
    ax.plot([a0, a1], [yy, yy], color=c, lw=6, solid_capstyle='butt', alpha=0.85)
ax.set_xscale('log'); ax.set_xlim(10, 2e4)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
ax.set_xlabel('FUV excess on filament (CU)')
ax.set_title('expected (I₅₀₀₇ = 10–13 R) vs measured (±2σ null)', loc='left', fontsize=9.5)
ax.grid(alpha=0.3, which='both', axis='x')

ax = axes[1, 2]
for band, col in [('murthy_fuv', OIII_C), ('murthy_nuv', WARN)]:
    for seg, ls in [('south', '-'), ('north', '--')]:
        sg = P['bands'][band]['segments'][seg]
        x = np.array(sg['profile_offsets']); yv = np.array(sg['profile'], float)
        ax.plot(x, yv - np.nanmedian(yv[(x > 0.8) & (x < 2.5)]), color=col, ls=ls, lw=1.2, label=f'{band.split("_")[1].upper()} {seg}')
ax.axvspan(-0.138, 0.138, color=OIII_C, alpha=0.12, lw=0); ax.axhline(0, color=FGDIM, lw=0.6)
ax.set_xlim(-4, 4); ax.set_xlabel('offset from spine (deg, + east)'); ax.set_ylabel('CU (relative to east side)')
ax.set_title('Murthy (2014) diffuse maps', loc='left', fontsize=9.5); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2)
finish(fig, FIG / 'imaging_mw_galex.png')
print('ok')
