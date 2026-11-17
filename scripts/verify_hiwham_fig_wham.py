"""verify_hiwham: WHAM summary with robust null widths, bootstrap comparison, and figure.

Reads data/derived/verify_hiwham_wham_spec.npz / _results.json (verify_hiwham_wham.py) and
verify_hiwham_hamaps_results.json (verify_hiwham_hamaps.py, if present).
Writes data/derived/verify_hiwham_wham_summary.json and figures/verify_hiwham_wham.png.

Robust null sigma = (p84 - p16)/2 over the Galactic-longitude-shifted trace copies (heavy tails from copies that
cross bright H II regions are thereby down-weighted); the bootstrap error of the on-off mean (resampling on beams,
the wave-1 approach) is given for comparison.
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, WARN, OK, ACCENT, FG
from verify_hiwham_common import DER, FIG, strand, Trace

Z = np.load(DER / 'verify_hiwham_wham_spec.npz')
R = json.loads((DER / 'verify_hiwham_wham_results.json').read_text())
cols = list(Z['cols'])
VG = Z['vgrid']
meta = Z['null_meta']          # dl, db, rot, mirror
rng = np.random.default_rng(5)


def rob(x):
    x = x[np.isfinite(x)]
    return float(np.median(x)), float(0.5 * (np.percentile(x, 84) - np.percentile(x, 16))), int(len(x))


summ = {'n_null_total': int(len(meta))}
near_mask = (np.abs(meta[:, 0]) <= 36) & (np.abs(meta[:, 1]) <= 6)
for key in ('onoff_south', 'onoff_north', 'reg_south', 'reg_north'):
    arr = Z[f'null_{key}']
    summ[key] = {}
    for q, nm in enumerate(cols):
        real = R[key][nm]['real']
        m, s, n = rob(arr[:, q]); m2, s2, n2 = rob(arr[near_mask, q])
        summ[key][nm] = {'real': real, 'null_median': m, 'null_sig_robust': s, 'z_robust': (real - m) / s, 'n': n,
                         'null_sig_robust_near(|dl|<=36,|db|<=6)': s2, 'n_near': n2,
                         'UL2_robust': max(real - m, 0) + 2 * s}

# bootstrap of the real on-off (beam resampling) for I_win / I_tot
d, s, D = Z['d'], Z['s'], Z['D_real']
meas = Z['meas']; I_tot = Z['I_tot']
s_split = R['s_split']
ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
d2, s2, ins = REAL.coords(Z['ra'], Z['dec'], dmax=4.0)
boot = {}
for nm, y in (('I_win', meas[:, 0]), ('I_tot', I_tot)):
    for part, sel in (('south', s2 <= s_split), ('north', s2 > s_split)):
        on = np.where(ins & sel & (np.abs(d2) < 0.5) & np.isfinite(y))[0]
        off = np.where(ins & (np.abs(d2) > 1.5) & (np.abs(d2) < 3.0) & np.isfinite(y))[0]
        diffs = []
        for j in on:
            mE = off[(d2[off] > 0) & (np.abs(s2[off] - s2[j]) < 0.75)]
            mW = off[(d2[off] < 0) & (np.abs(s2[off] - s2[j]) < 0.75)]
            parts = [np.mean(y[mm]) for mm in (mE, mW) if len(mm)]
            if parts:
                diffs.append(y[j] - np.mean(parts))
        diffs = np.array(diffs)
        bs = [np.mean(rng.choice(diffs, len(diffs))) for _ in range(4000)]
        boot[f'{nm}_{part}'] = {'mean': float(diffs.mean()), 'boot_sig': float(np.std(bs)), 'n_on': int(len(diffs))}
summ['bootstrap_onoff'] = boot

# main-component properties on and off (south / north), real values
on_s = ins & (np.abs(d2) < 0.5) & (s2 <= s_split); on_n = ins & (np.abs(d2) < 0.5) & (s2 > s_split)
off_s = ins & (np.abs(d2) > 1.5) & (np.abs(d2) < 3) & (s2 <= s_split); off_n = ins & (np.abs(d2) > 1.5) & (np.abs(d2) < 3) & (s2 > s_split)
props = {}
for nm, m in (('on_south', on_s), ('off_south', off_s), ('on_north', on_n), ('off_north', off_n)):
    props[nm] = {'n': int(m.sum()), 'I_g_mean': float(np.nanmean(meas[m, 4])), 'I_g_std': float(np.nanstd(meas[m, 4])),
                 'v_g_mean': float(np.nanmean(meas[m, 5])), 'v_g_std': float(np.nanstd(meas[m, 5])),
                 'fwhm_g_mean': float(np.nanmean(meas[m, 6])), 'fwhm_g_std': float(np.nanstd(meas[m, 6])),
                 'fwhm_int_mean': float(np.nanmean(np.sqrt(np.clip(meas[m, 6] ** 2 - 12.0 ** 2, 0, None)))),
                 'I_tot_mean': float(np.nanmean(I_tot[m])), 'D_mean': float(np.mean(D[m]))}
summ['component_props'] = props
summ['dilution'] = {'D_on_median': float(np.median(D[on_s | on_n])), 'D_on_max': float(D[on_s | on_n].max()),
                    'D_on_min': float(D[on_s | on_n].min()),
                    'note': 'beam-averaged unit-peak Gaussian filament (FWHM 0.276 deg) in a 1-deg top-hat beam'}
# limits on the intrinsic filament H-alpha peak (template regression, robust nulls)
lim = {}
cap40 = R['window_capture_fraction']['40.0']
for part in ('south', 'north'):
    for nm in ('I_tot', 'I_win'):
        x = summ[f'reg_{part}'][nm]
        ul = x['UL2_robust'] / (cap40 if nm == 'I_win' else 1.0)
        lim[f'{part}_{nm}'] = {'I_peak': x['real'] - x['null_median'], 'sig': x['null_sig_robust'], 'UL2_peak_R': ul,
                               'OIII_over_Ha_lower_for_2R': 2.0 / ul, 'OIII_over_Ha_lower_for_4R': 4.0 / ul}
summ['wham_peak_limits'] = lim
(DER / 'verify_hiwham_wham_summary.json').write_text(json.dumps(summ, indent=1))
print(json.dumps({k: summ[k] for k in ('bootstrap_onoff', 'component_props', 'dilution', 'wham_peak_limits')}, indent=1))
for key in ('onoff_south', 'onoff_north', 'reg_south', 'reg_north'):
    print(key, {nm: (round(x['real'], 3), round(x['null_sig_robust'], 3), round(x['z_robust'], 2)) for nm, x in summ[key].items()})

# ------------------------------------------------------------------ figure
fig, axes = new_sky_figure(15, 10.5, 'Silk Nebula · WHAM Hα on and off the spine',
                           'WHAM-SS DR1, 1° beams; revised [O III] V_LSR (V_helio +13.5) shaded; nulls = trace shifted in Galactic longitude',
                           nrows=2, ncols=3)
ax = axes.ravel()
vt = (0.3, 4.2)
for k, part in enumerate(('south', 'north')):
    on, off = Z[f'on_{part}'], Z[f'off_{part}']
    ax[k].plot(VG, on, color=HA_C, lw=1.6, label=f'on |d|<0.5° (n={int(Z[f"non_{part}"])})')
    ax[k].plot(VG, off, color=FGDIM, lw=1.4, label=f'off 1.5–3° (n={int(Z[f"noff_{part}"])})')
    ax[k].plot(VG, (on - off) * 5, color=WARN, lw=1.0, label='(on − off) × 5')
    ax[k].axvspan(*vt, color=OIII_C, alpha=0.25, lw=0)
    ax[k].axvspan(11.2, 15.5, color=FGDIM, alpha=0.12, lw=0)
    label(ax[k], 16, ax[k].get_ylim()[1] * 0.92 if ax[k].get_ylim()[1] > 0 else 0.1, 'paper [O III]', color=FGDIM, size=7)
    ax[k].axhline(0, color='#2a323c', lw=0.8)
    ax[k].set_xlim(-80, 80); ax[k].set_xlabel('V_LSR (km/s)'); ax[k].set_ylabel('R / (km/s)')
    ax[k].set_title(f'({"ab"[k]}) mean spectra, {part} (spine Dec {"≤" if part == "south" else ">"} 17.6)', loc='left')
    ax[k].legend(fontsize=7, loc='upper left')
# (c) beam dilution map
sel = ins & (np.abs(d2) < 3.5)
sc = ax[2].scatter(d2[sel], Z['dec'][sel], c=D[sel], s=40, cmap='magma', vmin=0, vmax=0.4, edgecolors='none')
ax[2].axvline(0, color=OIII_C, lw=0.8)
ax[2].axvspan(-0.5, 0.5, color=OIII_C, alpha=0.08, lw=0)
ax[2].set_xlabel('beam offset d (deg)'); ax[2].set_ylabel('beam Dec (deg)')
ax[2].set_title('(c) beam dilution D of a 0.28°-FWHM filament', loc='left')
cb = fig.colorbar(sc, ax=ax[2], fraction=0.04, pad=0.01); cb.set_label('D', fontsize=8)
# (d) template-regression peak vs nulls
for k, (part, col) in enumerate((('south', OIII_C), ('north', WARN))):
    q = cols.index('I_tot')
    arr = Z[f'null_reg_{part}'][:, q]
    ax[3].hist(np.clip(arr, -3, 3), bins=np.arange(-3, 3.01, 0.1), histtype='step', color=col, lw=1.4, density=True,
               label=f'nulls {part} (σ_r={summ[f"reg_{part}"]["I_tot"]["null_sig_robust"]:.2f} R)')
    ax[3].axvline(summ[f'reg_{part}']['I_tot']['real'], color=col, lw=2.0, ls='--')
ax[3].set_xlabel('fitted filament Hα peak I_peak (R), total Hα'); ax[3].set_ylabel('density')
ax[3].set_title('(d) template regression: real (dashed) vs nulls', loc='left')
ax[3].legend(fontsize=7)
# (e) on-off of component centroid / width vs nulls (south+north)
for k, (nm, xl) in enumerate((('v_g', 'Δ centroid (km/s)'), ('fwhm_g', 'Δ FWHM (km/s)'))):
    q = cols.index(nm)
    for part, col in (('south', OIII_C), ('north', WARN)):
        arr = Z[f'null_onoff_{part}'][:, q]
        lo, hi = np.nanpercentile(arr, [2.5, 97.5])
        y = 0.3 + 0.4 * k + (0.1 if part == 'north' else 0)
        ax[4].plot([lo, hi], [y, y], color=col, lw=4, alpha=0.4, solid_capstyle='butt')
        ax[4].plot(np.nanpercentile(arr, [16, 84]), [y, y], color=col, lw=7, alpha=0.6, solid_capstyle='butt')
        ax[4].plot(summ[f'onoff_{part}'][nm]['real'], y, 'o', color='#ffffff', ms=6, mec=col)
    label(ax[4], -7.8, 0.3 + 0.4 * k + 0.17, xl + ' (on − off)', color=FG, size=7.5)
ax[4].set_xlim(-8, 8); ax[4].set_ylim(0.1, 1.0); ax[4].set_yticks([])
ax[4].axvline(0, color='#2a323c')
ax[4].set_title('(e) main component on − off vs null 68/95 % (blue S, orange N)', loc='left', fontsize=9)
# (f) [O III]/Halpha lower limits summary
H = None
try:
    H = json.loads((DER / 'verify_hiwham_hamaps_results.json').read_text())
except Exception:
    pass
print('hamaps available:', None if H is None else list(H.keys()))
rows = [('WHAM 1° template, ±25 km/s', lim['south_I_win']['UL2_peak_R']),
        ('WHAM 1° template, total', lim['south_I_tot']['UL2_peak_R'])]
if H:
    for k2, lab in (('SHASSA_Sm', 'SHASSA Sm 4′'), ('SHASSA_CC', 'SHASSA CC (median)'),
                    ('Finkbeiner_SHASSApart', 'F03, SHASSA part')):
        if k2 in H:
            rows.append((lab, H[k2].get('peak_intrinsic_2sigma_UL_robust', H[k2]['peak_intrinsic_2sigma_UL'])))
try:      # displaced / broadened H-alpha filament limit (verify_hiwham_shassa_offset.py)
    SO = json.loads((DER / 'verify_hiwham_shassa_offset.json').read_text())
    rows.append(('SHASSA Sm, Hα displaced |d0|≤0.8°', max(SO['global'][k]['UL2_peak_R_max_over_d0'] for k in SO['global'])))
except Exception:
    SO = None
LAM = {'OIII': 5008.24, 'Ha': 6564.61, 'Hb': 4862.68}       # vacuum; only ratios of these enter
O3_DESI_REV = 1.6                                           # verify_desi revised south ridge median (R)
# --- aperture-matched spectroscopic limits (other wave-2 verifiers' products; NOT re-derived here)
sdr = json.loads((DER / 'verify_sdss_results_base_r0.15.json').read_text())['regions']['south_dec_lt17p6']
plc = json.loads((DER / 'verify_sdss_placebo_summary.json').read_text())['south'][1]
o3s = sdr['OIII:OIII5007']['I_R']
has, ehas = sdr['HaNII:Ha']['I_R'], np.hypot(sdr['HaNII:Ha']['err_boot'], plc['HaNII:Ha'][1])
hbs, ehbs = sdr['Hb:Hb']['I_R'], np.hypot(sdr['Hb:Hb']['err_boot'], plc['Hb:Hb'][1])
caseB_ph = 2.86 * LAM['Ha'] / LAM['Hb']                     # case-B Ha/Hb in photon units (3.86)
r_sdss_ha = o3s / (max(has, 0) + 2 * ehas)
r_sdss_hb = o3s / (caseB_ph * (max(hbs, 0) + 2 * ehbs))
DR = json.loads((DER / 'verify_desi_ratios.json').read_text())
e2p = LAM['OIII'] / LAM['Ha']                                # energy -> photon for [O III]/Ha
ha_all97 = DR['all']['Ha']['p97p5']
r_desi_ha_cons = e2p / ha_all97                              # all fibres, raw bootstrap 97.5 % (skewed tail)
hs = DR['south_dec<19.5']['Ha']
r_desi_ha_south = e2p / (hs['ratio'] + 1.5 * (hs['p97p5'] - hs['ratio']))   # south, bootstrap x1.5 (fake-line scatter)
o3hb_min = min(DR['derived'][k]['OIII5007/Hb']['lower_97p5'] for k in ('all', 'south_dec<19.5', 'petalFE', 'QSO_or_z>0.5', 'nbr_instrument'))
r_desi_hb = o3hb_min * (LAM['OIII'] / LAM['Hb']) / caseB_ph
r_desi_hb_35 = o3hb_min * (LAM['OIII'] / LAM['Hb']) / (3.5 * LAM['Ha'] / LAM['Hb'])   # if Ha/Hb(energy) up to 3.5
sw1 = json.loads((DER / 'spectra_sdss_lamost_sdss_substacks.json').read_text())['regions']['south_dec_lt_20']['lines']
r_w1 = sw1['OIII']['OIII5007']['I'] / (sw1['HaNII']['Ha']['I'] + 2 * sw1['HaNII']['Ha']['err'])
spec = [(f'SDSS Dec<17.6 Hα, +placebo (verify_sdss)', r_sdss_ha, None, OK, 'D', 'full'),
        (f'SDSS Dec<17.6 Hβ+case B, +placebo', r_sdss_hb, None, OK, 'o', 'full'),
        (f'DESI Hα direct, all 97.5 % (verify_desi)', r_desi_ha_cons, r_desi_ha_south, OIII_C, 'D', 'full'),
        (f'DESI Hβ (min subset >{o3hb_min:.0f}) + case B', r_desi_hb, r_desi_hb_35, OIII_C, 'o', 'full'),
        (f'wave-1 SDSS Dec<20 Hα (superseded)', r_w1, None, FGDIM, 'D', 'none')]
ny_ = len(rows) + len(spec)
for i, (lab, ul) in enumerate(rows):
    y = ny_ - 1 - i
    ax[5].plot([2.0 / ul, 4.0 / ul], [y, y], color=HA_C, lw=6, alpha=0.7, solid_capstyle='butt')
    ax[5].plot([O3_DESI_REV / ul], [y], '|', color='#ffffff', ms=10, mew=1.6)
    label(ax[5], 0.32, y + 0.3, f'{lab}: Hα peak < {ul:.2f} R', color=FG, size=6.5)
for k, (lab, r, r2, col, mk_, fill) in enumerate(spec):
    y = len(spec) - 1 - k
    if r2 is not None:
        ax[5].plot(sorted([r, r2]), [y, y], color=col, lw=1.5, alpha=0.8)
    ax[5].plot([r], [y], mk_, color=col, ms=6.5, mfc=col if fill == 'full' else 'none')
    label(ax[5], 0.32, y + 0.3, f'{lab}: > {r:.1f}' + (f' (→{r2:.1f})' if r2 is not None else ''), color=FG, size=6.5)
ax[5].axhline(len(spec) - 0.5, color='#2a323c', lw=0.8)
ax[5].set_xscale('log'); ax[5].set_xlim(0.3, 30); ax[5].set_ylim(-0.6, ny_ + 0.1); ax[5].set_yticks([])
ax[5].set_xlabel('[O III]5007/Hα 2σ lower limit (photon); bar: [O III] 2–4 R, tick: 1.6 R', fontsize=8)
ax[5].set_title('(f) south strand [O III]/Hα lower limits', loc='left')
summ['oiii_ha_lower_limits'] = {
    'units': 'photon (Rayleigh) ratios; multiply by 1.311 for energy',
    'rows_peak_UL_R': {lab: ul for lab, ul in rows},
    'map_ratio_for_OIII_1p6_2_4R': {lab: [O3_DESI_REV / ul, 2.0 / ul, 4.0 / ul] for lab, ul in rows},
    'sdss_dec_lt17p6_direct_placebo': r_sdss_ha, 'sdss_dec_lt17p6_hb_caseB_placebo': r_sdss_hb,
    'sdss_inputs_dec_lt17p6': {'OIII5007': o3s, 'Ha': [has, float(ehas)], 'Hb': [hbs, float(ehbs)],
                               'placebo_rms_Ha_Hb': [plc['HaNII:Ha'][1], plc['Hb:Hb'][1]]},
    'desi_ha_direct_all_raw97p5': r_desi_ha_cons, 'desi_ha_direct_south_boot1p5': r_desi_ha_south,
    'desi_hb_caseB_min_subset': r_desi_hb, 'desi_hb_HaHb3p5': r_desi_hb_35, 'desi_OIII_Hb_energy_min97p5': o3hb_min,
    'wave1_sdss_dec_lt20_direct': r_w1,
    'note': 'SDSS/DESI inputs are other wave-2 verifiers\' products (verify_sdss, verify_desi), propagated not re-derived.'}
(DER / 'verify_hiwham_wham_summary.json').write_text(json.dumps(summ, indent=1))
finish(fig, FIG / 'verify_hiwham_wham.png')
print('figure saved')
