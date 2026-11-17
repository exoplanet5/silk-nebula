#!/usr/bin/env python
"""spectra_ism: vet the automatic Ca II / Na I absorption fits (spectra_ism_absorption_measure.py) by eye-based
quality flags, flag possible stellar blends, and make the absorption figures.

Quality flags come from inspecting every normalised profile (diagnostic plots, not kept) together with SIMBAD
RV and v sin i (data/raw/spectra_ism_archives/simbad_{basic,rot}_spectra_stars.csv):
  good / ok : interstellar components usable; tentative : residual telluric or low S/N; bad : unusable.
Outputs: data/derived/spectra_ism_abs_curated.csv, spectra_ism_abs_curated.json,
         figures/spectra_ism_abs_profiles.png, figures/spectra_ism_abs_veldist.png
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_absorption_curate.py
"""
import csv
import io
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, FIG, vlsr_minus_vhelio  # noqa: E402
from spectra_ism_speclib import LINES, FVAL, C  # noqa: E402
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, WARN, OK, FG  # noqa: E402

G, O, T, B = 'good', 'ok', 'tentative', 'bad'
VET = {  # (star, instrument): {line: flag}, note
    ('HD 61717', 'UVES'): ({'CaK': G, 'CaH': G, 'NaD2': G, 'NaD1': O}, 'RV* +41.8, vsini 93: photospheric lines broad, removed'),
    ('HD 63432', 'UVES'): ({'CaK': B, 'CaH': B, 'NaD2': O, 'NaD1': T}, 'stellar Ca II residual (no SIMBAD RV/vsini)'),
    ('HD 63774', 'UVES'): ({'CaK': G, 'CaH': O, 'NaD2': G, 'NaD1': G}, ''),
    ('HD 68332', 'UVES'): ({'CaK': O, 'CaH': O, 'NaD2': G, 'NaD1': G}, 'broad stellar Ca II (RV* -9.4) under narrow IS line'),
    ('HD 53588', 'UVES'): ({'CaK': O, 'CaH': O, 'NaD2': O, 'NaD1': O}, 'sharp-lined B9 (RV* +16.4, vsini ~0): +16..+19 possibly stellar'),
    ('HD 52124', 'UVES'): ({'CaK': G, 'CaH': G, 'NaD2': G, 'NaD1': G}, ''),
    ('11 CMi', 'UVES'): ({'CaK': G, 'NaD2': T, 'NaD1': T}, 'vsini 285; Na D telluric residuals'),
    ('BR CMi', 'HARPS'): ({'CaK': B, 'CaH': B, 'NaD2': B, 'NaD1': B}, 'composite B9.5+G8III, low S/N'),
    ('8 Cnc', 'HARPS'): ({'CaK': O, 'CaH': T, 'NaD2': B, 'NaD1': T}, 'Na D telluric residuals'),
    ('HD 67959', 'HARPS'): ({'CaK': B, 'CaH': B, 'NaD2': B, 'NaD1': B}, 'sharp-lined A1V (vsini 18, RV* +23.8): stellar lines at IS velocities'),
    ('HD 52479', 'HARPS'): ({'CaK': T, 'CaH': T, 'NaD2': G, 'NaD1': G}, '963 pc control 13 deg W; stellar Ca II blend'),
    ('BD+17 1832', 'HARPS'): ({'CaK': B, 'CaH': B, 'NaD2': B, 'NaD1': B}, 'low S/N, stellar lines'),
    ('HD 61038', 'SOPHIE'): ({'CaK': O, 'CaH': B, 'NaD2': G, 'NaD1': G}, 'RV* +21.8, vsini 61: broad +17 Ca K residual is stellar'),
    ('HD 64705', 'SOPHIE'): ({'CaK': O, 'CaH': O, 'NaD2': O, 'NaD1': O}, 'sharp-lined B8 (RV* +31.5, vsini ~0): +25..+38 components likely stellar'),
    ('HD 66552', 'SOPHIE'): ({'CaK': O, 'CaH': T, 'NaD2': O, 'NaD1': B}, 'RV* +31.1, vsini 70'),
    ('HD 55869', 'SOPHIE'): ({'CaK': G, 'CaH': G, 'NaD2': G, 'NaD1': G}, 'RV* +28.1, vsini 44'),
    ('HD 57930', 'SOPHIE'): ({'CaK': B, 'CaH': B, 'NaD2': B, 'NaD1': B}, 'blue arm S/N too low; deep Na D at +36 probably stellar (no RV)'),
    ('HD 52710', 'SOPHIE'): ({'CaK': B, 'CaH': B, 'NaD2': T, 'NaD1': T}, 'low S/N'),
    ('64 Gem', 'SOPHIE'): ({'CaK': T, 'CaH': T, 'NaD2': B, 'NaD1': B}, 'low S/N; Na D telluric residuals'),
    ('HD 53588', 'SOPHIE'): ({'CaK': O, 'CaH': O, 'NaD2': O, 'NaD1': O}, 'see UVES; consistent velocities'),
    ('HD 67959', 'SOPHIE'): ({'CaK': B, 'CaH': B, 'NaD2': B, 'NaD1': B}, 'stellar'),
    ('HD 65900', 'SOPHIE'): ({'CaK': B, 'CaH': B, 'NaD2': T, 'NaD1': B}, 'stellar lines at RV* +50.7'),
    ('8 Cnc', 'SOPHIE'): ({'CaK': O, 'CaH': T, 'NaD2': B, 'NaD1': T}, ''),
    ('HD 71557', 'ELODIE'): ({'NaD2': O, 'NaD1': T}, 'ELODIE R=42000, no Ca II coverage'),
    ('HD 59605', 'ELODIE'): ({'NaD2': B, 'NaD1': B}, 'Am star, stellar Na D'),
    ('HD 58144', 'ELODIE'): ({'NaD2': T, 'NaD1': T}, ''),
    ('HD 54217', 'ELODIE'): ({'NaD2': O, 'NaD1': T}, ''),
    ('HD 71262', 'ELODIE'): ({'NaD2': O, 'NaD1': T}, ''),
    ('HD 71310', 'ELODIE'): ({'NaD2': T, 'NaD1': T}, ''),
    ('HD 59473', 'ELODIE'): ({'NaD2': T, 'NaD1': T}, ''),
    ('HD 56359', 'ELODIE'): ({'NaD2': O, 'NaD1': T}, ''),
    ('HD 64934', 'ELODIE'): ({'NaD2': O, 'NaD1': T}, ''),
    ('64 Gem', 'ELODIE'): ({'NaD2': B, 'NaD1': B}, 'telluric residuals'),
}
SHARP = {'HD 53588': 16.4, 'HD 64705': 31.5, 'HD 67959': 23.8}   # sharp-lined stars (vsini < 25 km/s): RV*

comps = list(csv.DictReader(open(DER / 'spectra_ism_abs_components.csv')))
CAND = {' '.join(r['name'].split()): (float(r['ra']), float(r['dec']))
        for r in csv.DictReader(open(DER / 'spectra_ism_archive_candidates.csv'))}
stars = list(csv.DictReader(open(DER / 'spectra_ism_abs_stars.csv')))
P = np.load(DER / 'spectra_ism_abs_profiles.npz')
VG = P['vgrid']
PROF = {k: P[f'p{i}'] for i, k in enumerate(P['keys'])}


def flag(star, inst, line):
    v = VET.get((star, inst))
    return v[0].get(line, B) if v else B


def forced(star, inst, line, v0, fwhm_obs):
    """EW (mA) and 1-sigma error of a Gaussian of given observed FWHM at v0 in a stored normalised profile."""
    y = PROF.get(f'{star}|{inst}|{line}')
    if y is None:
        return np.nan, np.nan
    s = fwhm_obs / 2.3548
    m = np.isfinite(y) & (VG > v0 - 3 * s - 5) & (VG < v0 + 3 * s + 5)
    off = np.isfinite(y) & (np.abs(VG) > 90) & (np.abs(VG) < 250)
    rms = 1.4826 * np.nanmedian(np.abs(y[off] - np.nanmedian(y[off])))
    prof = np.exp(-0.5 * ((VG[m] - v0) / s) ** 2)
    a = np.sum(prof * (1 - y[m])) / np.sum(prof ** 2)
    corr = np.sqrt(max(fwhm_obs / 3.0, 1.0) / 0.5 / 2.0)       # correlated-pixel inflation (~half resolution element)
    ea = rms / np.sqrt(np.sum(prof ** 2)) * corr
    k = s * np.sqrt(2 * np.pi) * LINES[line] / C * 1000
    return float(a * k), float(ea * k)


cur = []
for c in comps:
    q = flag(c['star'], c['inst'], c['line'])
    v, fw = float(c['v_helio']), float(c['fwhm_obs'])
    if q == B or fw > 20 or c['broad_flag'] == 'True' or not (-45 < v < 60):
        continue
    if float(c['EW_mA']) < 2.5 * float(c['eEW_mA']):
        continue
    note = ''
    if c['star'] in SHARP and abs(v - SHARP[c['star']]) < 10:
        note = 'possible stellar blend'
    # doublet confirmation for Na D2 (D1 at same v within 3 km/s) and Ca K (H)
    partner = {'NaD2': 'NaD1', 'NaD1': 'NaD2', 'CaK': 'CaH', 'CaH': 'CaK'}[c['line']]
    pc = [p for p in comps if p['star'] == c['star'] and p['inst'] == c['inst'] and p['line'] == partner
          and abs(float(p['v_helio']) - v) < 3.5]
    conf = bool(pc)
    ra_, de_ = CAND.get(' '.join(c['simbad'].split()), (116.5, 18.0))
    dvl = float(vlsr_minus_vhelio(ra_, de_)[0])
    cur.append(dict(star=c['star'], inst=c['inst'], line=c['line'], quality=q, d_pc=float(c['d_pc']),
                    sep_spine_deg=float(c['sep_spine_deg']), dperp_deg=float(c['dperp_deg']), dec_spine=float(c['dec_spine']),
                    v_helio=round(v, 2), ev=round(float(c['ev']), 2), v_lsr=round(v + dvl, 2), fwhm_obs=round(fw, 1),
                    b_kms=round(float(c['b_kms']), 1), EW_mA=round(float(c['EW_mA']), 2), eEW_mA=round(float(c['eEW_mA']), 2),
                    logN_AOD=c['logN_AOD'], doublet_confirmed=conf, note=note))

# Specific tests
tests = {}
# 1) HD 61717 high-velocity Ca II component: Na D limits at the same velocity / width
hv = [c for c in cur if c['star'] == 'HD 61717' and c['line'] in ('CaK', 'CaH') and c['v_helio'] < -5]
if hv:
    v0 = np.mean([c['v_helio'] for c in hv])
    fw = np.mean([c['fwhm_obs'] for c in hv if c['line'] == 'CaK'])
    tests['HD61717_hv'] = {'v_helio': round(float(v0), 2), 'components': hv}
    for ln in ('NaD2', 'NaD1', 'CaK', 'CaH'):
        ew, e = forced('HD 61717', 'UVES', ln, v0, max(fw, 6.0))
        tests['HD61717_hv'][f'forced_{ln}'] = [round(ew, 2), round(e, 2)]
    ewk, eewk = tests['HD61717_hv']['forced_CaK']
    ewn, eewn = tests['HD61717_hv']['forced_NaD2']
    nca = 1.130e17 * ewk / (FVAL['CaK'] * LINES['CaK'] ** 2)
    nna = 1.130e17 * ewn / (FVAL['NaD2'] * LINES['NaD2'] ** 2)
    rat = nca / nna
    erat = rat * np.hypot(eewk / ewk, eewn / ewn)
    tests['HD61717_hv'].update(logN_CaII_AOD=round(float(np.log10(nca)), 2), logN_NaI_AOD=round(float(np.log10(nna)), 2),
                               CaII_over_NaI=round(float(rat), 1), e_CaII_over_NaI=round(float(erat), 1),
                               K_over_H=round(ewk / tests['HD61717_hv']['forced_CaH'][0], 2),
                               note='Na I weakly detected in D2 (and D1 at ~3 sigma) at the same velocity; Ca II K/H = 1.7 implies mild '
                                    'saturation, so N(Ca II) and the ratio are lower limits')
# 2) forced EW at the [O III] velocity (+25 km/s helio, b = 5 km/s) for good/ok lines
v25 = []
for s in stars:
    q = flag(s['star'], s['inst'], s['line'])
    if q not in (G, O):
        continue
    fw = float(np.hypot(C / float(s['R']), 5.0 * 1.665))
    ew, e = forced(s['star'], s['inst'], s['line'], 25.0, fw)
    v25.append(dict(star=s['star'], inst=s['inst'], line=s['line'], quality=q, d_pc=float(s['d_pc']),
                    sep_spine_deg=float(s['sep_spine_deg']), EW25_mA=round(ew, 2), eEW25_mA=round(e, 2),
                    resolved_comp_22_30=[c['v_helio'] for c in cur if c['star'] == s['star'] and c['inst'] == s['inst']
                                         and c['line'] == s['line'] and 22 <= c['v_helio'] <= 30]))
tests['v25'] = v25

with open(DER / 'spectra_ism_abs_curated.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(cur[0].keys()))
    w.writeheader()
    w.writerows(cur)
json.dump({'vetting': {f'{k[0]}|{k[1]}': {'flags': v[0], 'note': v[1]} for k, v in VET.items()}, 'tests': tests},
          open(DER / 'spectra_ism_abs_curated.json', 'w'), indent=1, default=float)
print('curated components', len(cur))
for c in sorted(cur, key=lambda c: (c['d_pc'], c['star'], c['line'], c['v_helio'])):
    print(f"{c['star']:10s} {c['inst']:6s} {c['line']:5s} {c['quality']:4s} d={c['d_pc']:6.1f} sep={c['sep_spine_deg']:5.2f} "
          f"v={c['v_helio']:+6.1f}+-{c['ev']:.1f} fw={c['fwhm_obs']:5.1f} EW={c['EW_mA']:6.1f}+-{c['eEW_mA']:.1f} conf={c['doublet_confirmed']} {c['note']}")
print(json.dumps(tests.get('HD61717_hv', {}), default=float)[:900])

# ------------------------------------------------------------------------------------------------
# Figure 1: Ca II K and Na I D2 profiles of the vetted sightlines, sorted by distance
# ------------------------------------------------------------------------------------------------
best = [('8 Cnc', 'HARPS'), ('11 CMi', 'UVES'), ('HD 66552', 'SOPHIE'), ('HD 61038', 'SOPHIE'), ('HD 68332', 'UVES'),
        ('HD 63774', 'UVES'), ('HD 53588', 'UVES'), ('HD 63432', 'UVES'), ('HD 52124', 'UVES'), ('HD 61717', 'UVES'),
        ('HD 64705', 'SOPHIE'), ('HD 55869', 'SOPHIE')]
meta = {(s['star'], s['inst']): s for s in stars}
best = [b for b in best if b in meta]
best.sort(key=lambda b: float(meta[b]['d_pc']))
fig, axs = new_sky_figure(12.5, 10.5, 'Silk Nebula · interstellar Ca II K and Na I D2',
                          'archival UVES / HARPS / SOPHIE spectra; sorted by Gaia distance; offsets for clarity',
                          nrows=1, ncols=2, sharey=True)
for ax, ln in zip(axs, ('CaK', 'NaD2')):
    for i, b in enumerate(best):
        y = PROF.get(f'{b[0]}|{b[1]}|{ln}')
        off = 0.55 * i
        q = flag(b[0], b[1], ln)
        if y is not None:
            m = (VG > -80) & (VG < 100)
            yy = np.clip(y[m], 0.45, 1.15)
            ax.plot(VG[m], yy - 1 + off + 1, color=FG if q in (G, O) else FGDIM, lw=0.8)
        s = meta[b]
        txt = f"{b[0]}  {float(s['d_pc']):.0f} pc  {float(s['sep_spine_deg']):.1f}°" + ('' if q in (G, O) else f'  ({q})')
        label(ax, -78, off + 1.07, txt, size=7.2, color='#cfd6dd')
    ax.axvspan(23.8, 26.8, color=OIII_C, alpha=0.25, lw=0)
    ax.axvline(12.5, color=HA_C, lw=0.8)
    ax.axvline(19.2, color=OK, lw=0.8, ls=':')
    ax.axvline(33.9, color=WARN, lw=0.8, ls=':')
    ax.set_xlim(-80, 100)
    ax.set_xlabel('V_helio (km/s)')
    ax.set_title({'CaK': 'Ca II K 3934', 'NaD2': 'Na I D2 5890 (telluric-corrected)'}[ln], loc='left')
axs[0].set_ylabel('normalised flux + offset')
axs[0].set_ylim(0.4, 0.55 * len(best) + 0.75)
axs[0].plot([], [], color=OIII_C, lw=6, alpha=0.4, label='[O III]')
axs[0].plot([], [], color=HA_C, lw=0.8, label='red lines')
axs[0].plot([], [], color=OK, lw=0.8, ls=':', label='LIC')
axs[0].plot([], [], color=WARN, lw=0.8, ls=':', label='Gem cloud')
axs[0].legend(loc='upper right', ncol=2, fontsize=7.5)
finish(fig, FIG / 'spectra_ism_abs_profiles.png')

# ------------------------------------------------------------------------------------------------
# Figure 2: component velocity vs distance, and literature column densities vs distance
# ------------------------------------------------------------------------------------------------
lit = list(csv.DictReader(open(DER / 'spectra_ism_abs_catalog_stars.csv')))
fig, axs = new_sky_figure(13, 5.6, 'Silk Nebula · absorption components vs distance',
                          'left: vetted components (size ∝ EW); right: literature totals (Welsh+2010, Lallement+2003)',
                          nrows=1, ncols=3)
ax = axs[0]
for ln, col, mk in (('CaK', OIII_C, 'o'), ('NaD2', WARN, 's')):
    cc = [c for c in cur if c['line'] == ln and c['quality'] in (G, O)]
    d = np.array([c['d_pc'] for c in cc])
    v = np.array([c['v_helio'] for c in cc])
    ew = np.array([c['EW_mA'] for c in cc])
    near = np.array([c['sep_spine_deg'] < 3 for c in cc])
    st = np.array([c['note'] != '' for c in cc])
    jit = 1.0 if ln == 'CaK' else 1.03
    ax.scatter(d[near & ~st] * jit, v[near & ~st], s=8 + 3 * ew[near & ~st], c=col, marker=mk, edgecolors='white', lw=0.6,
               label=f'{ln[:2]} {ln[2:]} <3° from spine', zorder=3)
    ax.scatter(d[~near & ~st] * jit, v[~near & ~st], s=8 + 3 * ew[~near & ~st], facecolors='none', edgecolors=col, marker=mk, lw=0.8,
               label=f'{ln[:2]} {ln[2:]} 3–11°', zorder=3)
    ax.scatter(d[st] * jit, v[st], s=8 + 3 * ew[st], c='#555', marker='x', lw=0.8, zorder=2)
ax.axhspan(23.8, 26.8, color=OIII_C, alpha=0.2, lw=0)
ax.axhline(12.5, color=HA_C, lw=0.8)
ax.axhline(19.2, color=OK, lw=0.8, ls=':')
ax.axhline(33.9, color=WARN, lw=0.8, ls=':')
ax.set_xscale('log')
ax.set_xlim(40, 700)
ax.set_ylim(-35, 50)
ax.set_xlabel('distance (pc)')
ax.set_ylabel('V_helio (km/s)')
ax.legend(loc='upper left', fontsize=6.5, ncol=2, markerscale=0.6)
if 'HD61717_hv' in tests:
    label(ax, 245, tests['HD61717_hv']['v_helio'], 'HD 61717: Ca II/Na I ≈ 22 ', size=7, ha='right', va='center')
for ax, col_, lab in ((axs[1], 'logN_NaI', 'log N(Na I) (cm⁻²)'), (axs[2], 'logN_CaII', 'log N(Ca II) (cm⁻²)')):
    for r in lit:
        try:
            d = float(r['d_best_pc'])
            y = float(r[col_])
        except ValueError:
            continue
        if d > 1500:
            continue
        lim = r['l_' + col_] == '<'
        sep = float(r['sep_spine_deg'])
        c_ = OIII_C if sep < 3 else (FG if sep < 6 else FGDIM)
        ax.plot(d, y, 'v' if lim else 'o', color=c_, ms=5 if sep < 3 else 4, mfc=c_ if not lim else 'none')
        if sep < 5 and not lim:
            label(ax, d * 1.05, y, r['name'].replace('HD ', ''), size=6.5)
    # our measurements (sum of vetted components, Ca K and Na D2)
    ln = 'NaD2' if 'Na' in col_ else 'CaK'
    for st_ in sorted(set(c['star'] for c in cur)):
        cc = [c for c in cur if c['star'] == st_ and c['line'] == ln and c['quality'] in (G, O) and c['note'] == '']
        if not cc:
            continue
        ew = sum(c['EW_mA'] for c in cc)
        N = 1.130e17 * ew / (FVAL[ln] * LINES[ln] ** 2)
        ax.plot(cc[0]['d_pc'], np.log10(N), 'D', color=WARN, ms=4.5, mfc='none', mew=1.0)
    ax.set_xscale('log')
    ax.set_xlim(8, 1500)
    ax.set_xlabel('distance (pc)')
    ax.set_ylabel(lab)
axs[1].plot([], [], 'o', color=OIII_C, label='literature, <3° from spine')
axs[1].plot([], [], 'o', color=FG, label='3–6°')
axs[1].plot([], [], 'o', color=FGDIM, label='6–12°')
axs[1].plot([], [], 'v', color=FGDIM, mfc='none', label='upper limit')
axs[1].plot([], [], 'D', color=WARN, mfc='none', label='this work (AOD)')
axs[1].legend(loc='upper left', fontsize=6.5)
finish(fig, FIG / 'spectra_ism_abs_veldist.png')
