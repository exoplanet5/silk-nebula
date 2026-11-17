"""w4_refit_fig: summary figure for the w4 refit (shock grid with corrected inputs, alternatives, geometry).
Inputs: data/derived/w4_refit_fit.json, w4_refit_geom.json, w4_refit_alt.json, w4_refit_tracks.npz
Output: figures/w4_refit_summary.png (dpi 200)
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_fig.py
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *
import w4_refit_geom as GEO

ROOT = '/Users/mickey/astronomy/silk-nebula'
fit = json.load(open(ROOT + '/data/derived/w4_refit_fit.json'))
geo = json.load(open(ROOT + '/data/derived/w4_refit_geom.json'))
alt = json.load(open(ROOT + '/data/derived/w4_refit_alt.json'))
T = np.load(ROOT + '/data/derived/w4_refit_tracks.npz')
runs = {r['tag']: r for r in fit['runs']}
VA = fit['variants']

fig, axs = new_sky_figure(15, 9.2, 'Silk Nebula · shock refit with corrected inputs',
                          'w4_refit: truncated-shock grid (verify_shockmodel code) and altmodel scenarios on one likelihood; '
                          'geometry: curved limb vs flat vs rippled sheet', nrows=2, ncols=3)
ax = axs[0, 0]
cols = {'A': ACCENT, 'B': WARN, 'Afl': OK}
labs = {'A': 'A: DESI Hβ limit', 'B': 'B: + SDSS ~23', 'Afl': 'A, 4363 floor'}
for base in ('A', 'B', 'Afl'):
    cb = VA[base]['chi2_by']['v']
    vv = sorted(int(k) for k in cb)
    y = [min(cb[str(v)]['nat'], cb[str(v)]['R88']) for v in vv]
    ax.plot(vv, y, 'o-', color=cols[base], ms=3.5, lw=1.4, label=labs[base])
    cm = VA[base]['complete']
    for tag, c in cm.items():
        v = int(tag.split('_')[0][1:])
        ax.plot(v, min(c['nat'], c['R88']), 'x', color=cols[base], ms=5, alpha=0.8)
mn = min(VA['A']['chi2_min_nat'], VA['A']['chi2_min_R88'])
ax.axhline(mn + 4, color=FGDIM, ls=':', lw=1)
label(ax, 62, mn + 4.6, 'A: χ²min + 4', color=FGDIM, size=7.5)
ax.set_yscale('log'); ax.set_xlabel('v_s (km/s)'); ax.set_ylabel('best χ² (south, 10–11 terms)')
ax.set_title('(a) truncated shocks (lines) vs complete (×)')
ax.legend(loc='upper right', fontsize=7)
ax.grid(alpha=0.3)

# (b) diagnostic plane
ax = axs[0, 1]
acc = set(VA['A']['acceptable'])
for tag, r in runs.items():
    if r['complete']:
        continue
    q = r['A_nat']['ratios']
    c = ACCENT if tag in acc else '#3a4450'
    ax.plot(1 / q['O3_Hb'], q['O2_O3'], '.', color=c, ms=3 if tag in acc else 1.5, alpha=0.9 if tag in acc else 0.5, zorder=3 if tag in acc else 1)
for tag, r in runs.items():
    if r['complete']:
        q = r['A_nat']['ratios']
        ax.plot(1 / q['O3_Hb'], q['O2_O3'], 'x', color=DANGER, ms=5)
accN = set(VA['N']['acceptable'])
for tag in accN:
    q = runs[tag]['N_nat']['ratios']
    ax.plot(1 / q['O3_Hb'], q['O2_O3'], '.', color=WARN, ms=3)
ax.add_patch(plt.Rectangle((0, 1.06), 0.095, 0.16, fc='none', ec=FG, lw=1.2))
label(ax, 0.0055, 0.78, 'south box', size=7.5)
ax.add_patch(plt.Rectangle((0.067, 4.0), 0.27, 4.0, fc='none', ec=WARN, lw=1.2, ls='--'))
label(ax, 0.07, 8.4, 'north', color=WARN, size=7.5)
mk = {'warm_CIE': ('s', SUN, 'warm CIE'), 'conduction_static': ('^', NII_C, 'cond. front'),
      'NEI_cooling_TML': ('D', HA_C, 'NEI / TML'), 'PIE_hotgas_kT140eV': ('v', FGDIM, 'PIE hot gas'),
      'PIE_BB81kK_WDJ0736': ('v', '#888888', 'PIE WD'), 'PIEheat_hotgas140eV': ('P', OK, 'PIE + heat'),
      'fossil_T0=4e+04|nH=0.5': ('*', '#ff9de2', 'fossil 4e4 K')}
for k, (m, c, lab) in mk.items():
    b = alt['families'][k]['best_A']['ratios']
    ax.plot(1 / b['O3_Hb'], b['O2_O3'], m, color=c, ms=7, mec='k', mew=0.5, label=lab, zorder=5)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(5e-3, 2); ax.set_ylim(0.3, 30)
ax.set_xlabel('Hβ / [O III]5007'); ax.set_ylabel('[O II]3727 / [O III]5007')
ax.set_title('(b) shocks (blue acc., grey rest; red × complete; orange north)')
ax.legend(loc='lower right', fontsize=6.5, ncol=2)
ax.grid(alpha=0.3)

# (c) cross-profiles of the best model, three geometries, blurred to the stacked FWHM
ax = axs[0, 2]
best = min(geo['models'], key=lambda m: m['chi2'])
r = dict(runs[best['tag']]); r['pick'] = {'i': runs[best['tag']]['A_nat']['i']}
tr = GEO.track(T, best['tag'])
i = r['pick']['i']; d = GEO.D_FID
n0 = best['flat']['n0']
xg, jg = GEO.flat_profile(tr, i, n0)
xdeg = np.degrees(xg / (d * GEO.PC))
xp, yp = GEO.blur_to(xdeg, jg / jg.max(), GEO.W_STACK)
ax.plot(-(xp - xdeg[np.argmax(jg)]), yp / yp.max(), color=OK, lw=1.6, label=f"flat sheet  σW/σE {best['flat']['asym_WE']:.2f}")
Rc = d * np.radians(11.9) * GEO.PC
pr = GEO.SA.profile(tr, i, best['curved']['n0'], 0.1, Rc)
hu = np.linspace(0, pr['hs'].max() * 0.999, 4000); Iu = np.interp(hu, pr['hs'], pr['I'])
xd = np.degrees(hu / (d * GEO.PC)); m = xd < 1.5
xp, yp = GEO.blur_to(xd[m], Iu[m] / Iu.max(), GEO.W_STACK)
ax.plot(-(xp - xd[m][np.argmax(Iu[m])]), yp / yp.max(), color=DANGER, lw=1.6, label=f"curved 58-pc limb  σW/σE {best['curved']['asym_WE']:.2f}")
ax.axvspan(-0.02, 0.02, color=FGDIM, alpha=0.15)
label(ax, -0.55, 0.92, 'observed σW/σE 0.98 ± 0.05', size=7.5)
label(ax, 0.35, 0.05, 'east (front) →', size=7.5, color=FGDIM)
label(ax, -0.62, 0.05, '← west (Monogem)', size=7.5, color=FGDIM)
ax.set_xlim(-0.7, 0.7); ax.set_ylim(0, 1.08)
ax.set_xlabel('offset across ridge (deg, + = east)'); ax.set_ylabel('normalised [O III] (blurred to 0.41°)')
ax.set_title(f"(c) best model {best['tag'].replace('_trn', '')}")
ax.legend(loc='upper right', fontsize=7)
ax.grid(alpha=0.3)

# (d) FUV per R vs v_s (flat edge-on sheet, C IV scattered; curved limb), NUV
ax = axs[1, 0]
S = geo['summary']
summ = json.load(open(ROOT + '/data/derived/w4_refit_summary.json'))
rng_ = np.random.default_rng(3)
xc = {0.0: ACCENT, 0.05: OK, 0.3: WARN, 1.0: DANGER}
for m in geo['models']:
    jit = rng_.uniform(-2.5, 2.5)
    ax.plot(m['v'] + jit, m['flat']['CIV_span0']['FUV'], 'o', color=xc[m['xn']], ms=2.8, alpha=0.75)
    ax.plot(m['v'] + jit, m['uv']['FUV_nonres'], '_', color=FGDIM, ms=5, alpha=0.6)
    ax.plot(m['v'] + jit, m['uv']['NUV'], '.', color=HB_C, ms=2, alpha=0.5)
for X, c in xc.items():
    ax.plot([], [], 'o', color=c, ms=4, label=f'FUV flat, X(H0)={X}')
ax.plot([], [], '_', color=FGDIM, ms=6, label='FUV non-resonant floor')
ax.plot([], [], '.', color=HB_C, ms=5, label='NUV (all models)')
ax.axhspan(68, 92, color=ACCENT, alpha=0.18); label(ax, 197, 70, 'GALEX FUV 68–92', color=ACCENT, size=7.5, ha='right')
ax.axhspan(36, 57, color=HB_C, alpha=0.12); label(ax, 197, 39, 'NUV 36–57', color=HB_C, size=7.5, ha='right')
wc = alt['families']['warm_CIE']['best_A']['uv']
ax.axhline(wc['NUV'], color=SUN, ls=':', lw=1); label(ax, 72, wc['NUV'] * 1.05, 'warm CIE NUV', color=SUN, size=7)
ax.set_yscale('log'); ax.set_ylim(25, 500); ax.set_xlim(65, 205)
ax.set_xlabel('v_s (km/s)'); ax.set_ylabel('CU per R of [O III]5007')
ax.set_title('(d) UV per R, spectrally acceptable shocks')
ax.legend(loc='upper left', fontsize=6.3)
ax.grid(alpha=0.3)

# (e) alternatives: spectral chi2 (A) and the budget/X-ray verdicts
ax = axs[1, 1]
fams = [('shock (best)', None, ACCENT), ('warm CIE 6e4 K', 'warm_CIE', SUN), ('NEI cool 1.6 MK', 'NEI_from_1.6MK', HA_C),
        ('NEI cool 3e5 K', 'NEI_from_3e5K', HA_C), ('cond. static', 'conduction_static', NII_C),
        ('cond. evap.', 'conduction_evap', NII_C), ('PIE hot gas', 'PIE_hotgas_kT140eV', FGDIM),
        ('PIE WD 81 kK', 'PIE_BB81kK_WDJ0736', FGDIM), ('PIE 150 kK', 'PIE_BB150kK', FGDIM),
        ('PIE + heat', 'PIEheat_hotgas140eV', OK), ('fossil 4e4 K', 'fossil_T0=4e+04|nH=0.5', '#ff9de2'),
        ('fossil 2e4 K', 'fossil_T0=2e+04|nH=0.5', '#ff9de2')]
from scipy.stats import chi2 as CH
for k, (lab, key, c) in enumerate(fams):
    if key is None:
        b = VA['A']['best']; x2 = min(b['nat']['chi2'], b['R88']['chi2']); n = b['nat']['nterms']
        x2f = min(VA['Afl']['chi2_min_nat'], VA['Afl']['chi2_min_R88'])
    else:
        b = alt['families'][key]
        x2, n = b['best_A']['score']['A']['chi2'], b['best_A']['score']['A']['n']
        x2f = b['best_floor']['score']['Afl']['chi2']
    ax.barh(k, x2, color=c, alpha=0.85)
    ax.plot(x2f, k, '|', color='w', ms=9, mew=1.5)
    label(ax, 98, k, f'{n} terms', size=6.5, color=FGDIM, va='center', ha='right')
ax.set_yticks(range(len(fams))); ax.set_yticklabels([f[0] for f in fams], fontsize=7.5)
ax.invert_yaxis(); ax.set_xscale('log'); ax.set_xlim(0.2, 100)
ax.axvline(CH.ppf(0.95, 9), color=FGDIM, ls=':', lw=1)
label(ax, CH.ppf(0.95, 9) * 1.05, len(fams) - 0.6, 'p = 0.05 (9 dof)', size=6.5, color=FGDIM)
ax.set_xlabel('optical χ² (likelihood A; | = 4363 floor)')
ax.set_title('(e) optical spectrum alone no longer excludes most')
ax.grid(alpha=0.3, axis='x')

# (f) predicted [O III] FWHM for the survivors of each geometry
ax = axs[1, 2]
summ = json.load(open(ROOT + '/data/derived/w4_refit_summary.json'))
fp, cs = summ['flat_plausible_L83'], summ['curved_without_asym']
okf = lambda x: abs(x - 80) <= 2 * np.hypot(15, 0.2 * x)
def rip(R, N):
    v = [m['ripple'][f'R{R}_line_fwhm'] for m in geo['models']
         if m['ripple'][f'R{R}_N{N}']['fwhm1_deg'] <= 0.28 and okf(m['ripple'][f'R{R}_CIV']['FUV'])]
    return [min(v), float(np.median(v)), max(v)]
W = [('flat 0°', fp['line_fwhm_i0'], OK), ('flat 5°', fp['line_fwhm_i5'], OK), ('flat 10°', fp['line_fwhm_i10'], OK),
     ('limb 58', cs['line_fwhm'], DANGER), ('rip. 20', rip(20, 2), SII_C), ('rip. 10', rip(10, 2), SII_C)]
for k, (lab, q, c) in enumerate(W):
    ax.plot([k, k], [q[0], q[2]], color=c, lw=5, alpha=0.85, solid_capstyle='butt')
    ax.plot(k, q[1], 'o', color='w', ms=3)
for k, sh in ((1, 'vshift_i5'), (2, 'vshift_i10')):
    label(ax, k, fp['line_fwhm_i0'][2] + 3, f"Δv {fp[sh][0]:.0f}–{fp[sh][2]:.0f}", size=6.3, color=OK, ha='center')
k = len(W)
for x, lab, c in ((13.2, 'warm\n6e4 K', SUN), (6.1, 'PIE\n1.2e4 K', FGDIM)):
    ax.plot(k, x, 's', color=c, ms=7, mec='k'); label(ax, k, x + 3.5, lab, size=6.3, color=c, ha='center')
    k += 1
ax.axhline(33, color=ACCENT, ls='--', lw=1); label(ax, -0.3, 34.5, 'DESI ≤ 33', color=ACCENT, size=7.5)
ax.set_xticks(range(len(W) + 2)); ax.set_xticklabels([w[0] for w in W] + ['warm', 'PIE'], fontsize=6.8, rotation=35)
ax.set_xlim(-0.6, k - 0.4)
ax.set_ylabel('[O III]5007 FWHM, thermal + flow only (km/s)'); ax.set_ylim(0, 60)
ax.set_title('(f) line width tests geometry (survivors)')
ax.grid(alpha=0.3, axis='y')
finish(fig, ROOT + '/figures/w4_refit_summary.png')
print('ok')
