"""verify_shockmodel_figures: figures for the shockmodel verification (dpi 200, silk_style).
  figures/verify_shockmodel_tracks.png   cumulative line ratios vs swept column (truncation) + complete-shock values
  figures/verify_shockmodel_closure.png  brightness / pressure / line-width / age closure over (v_s, b)
  figures/verify_shockmodel_fuv.png      GALEX FUV per R(5007): optically thin vs resonance-scattered; profiles
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *
import verify_shockmodel_fit as V
import verify_shockmodel_analysis as AN

ROOT = '/Users/mickey/astronomy/silk-nebula'
fit = json.load(open(ROOT + '/data/derived/verify_shockmodel_fit.json'))
an = json.load(open(ROOT + '/data/derived/verify_shockmodel_analysis.json'))
summ = json.load(open(ROOT + '/data/derived/verify_shockmodel_summary.json'))
TR = np.load(ROOT + '/data/derived/verify_shockmodel_tracks.npz')
runs = {r['tag']: r for r in fit['runs']}


def tr(tag):
    return {k.split('__', 1)[1]: TR[k] for k in TR.files if k.startswith(tag + '__')}


def ratios_track(t):
    G = {g: t['F_' + g] for g in V.GROUPS}
    return V.ratios(G)


# ------------------------------------------------------------------ figure 1: tracks
fig, axs = new_sky_figure(12, 8.4, 'Silk · truncated radiative shocks: cumulative line ratios vs swept column',
                          'own 1-D NEI steady-flow model (Asplund09, b = B0/√n0 = 10 µG cm^1.5); grey bands = observed (wave 2); '
                          '◆ = best truncation; ▶ = complete-shock end', nrows=2, ncols=2)
axs = axs.ravel()
cols = {100: '#7ad97a', 120: OIII_C, 140: '#c792ea', 160: WARN}
show = [(v, 0.05, 10.0, '-') for v in (100, 120, 140, 160)] + [(120, 1.0, 10.0, '--')]
obsband = {'O3_Hb': (24, 60, 46), 'O2_O3': (1.07, 1.84, 1.35), 'O4363': (0.059, 0.096, 0.087), 'CIV_O3': None}
for (v, x, b, ls) in show:
    tag = f'v{v}_x{x}_b{b}_asplund09_trn'
    if tag not in runs:
        continue
    t = tr(tag)
    r = ratios_track(t)
    lN = np.log10(t['N'])
    i = runs[tag]['best']['i']
    lab = f'{v} km/s, X(H0)={x}'
    for ax, key in zip(axs[:3], ('O3_Hb', 'O2_O3', 'O4363')):
        ax.plot(lN, r[key], ls, color=cols[v], lw=1.4, label=lab)
        ax.plot(lN[i], r[key][i], 'D', color=cols[v], ms=5, mec='k')
    axs[3].plot(lN, r['CIV_O3'], ls, color=cols[v], lw=1.4, label=lab)
    axs[3].plot(lN[i], r['CIV_O3'][i], 'D', color=cols[v], ms=5, mec='k')
    ax2 = axs[3]
# complete endpoints (b = 3) at right edge
for tag, rr in runs.items():
    if not rr['complete']:
        continue
    t = tr(tag)
    G = {g: t['F_' + g][-1:] for g in V.GROUPS}
    r = V.ratios(G)
    vv = rr['v']
    c = cols.get(vv, FGDIM)
    for ax, key in zip(axs, ('O3_Hb', 'O2_O3', 'O4363', 'CIV_O3')):
        ax.plot(20.15, r[key][0], '>', color=c, ms=5, alpha=0.8)
for ax, key, yl in zip(axs, ('O3_Hb', 'O2_O3', 'O4363', 'CIV_O3'),
                       ('[O III]5007 / Hβ', '[O II]3727 / [O III]5007', '[O III]4363 / 5007', 'C IV 1549 / [O III]5007 (energy, thin)')):
    ax.set_yscale('log')
    ax.set_xlim(15.0, 20.3)
    ax.set_xlabel('log swept column N_H behind the shock (cm⁻²)')
    ax.set_ylabel(yl)
    ax.grid(alpha=0.3)
    if obsband.get(key):
        lo, hi, c = obsband[key]
        ax.axhspan(lo, hi, color='#9aa4ae', alpha=0.18, lw=0)
        ax.axhline(c, color='#cfd6dd', lw=0.8, ls=':')
axs[0].set_ylim(1, 3000)
axs[0].axhspan(1.2, 4.5, xmin=0.93, xmax=1.0, color=DANGER, alpha=0.35, lw=0)
label(axs[0], 19.2, 1.3, 'published\ncomplete\n1.2–4.5', color=DANGER, size=7)
for nm, lN_, val in (('R88 SE', 18.20, 12.9), ('R88 M', 18.25, 10.5)):
    axs[0].plot(lN_, val, 's', color=SUN, ms=5)
    label(axs[0], lN_ + 0.05, val * 1.25, nm, color=SUN, size=7)
label(axs[0], 15.1, 1.6, 'obs: DESI-S 46 (>28), SDSS-S 23', color=FG, size=7)
axs[1].set_ylim(0.01, 30)
axs[2].set_ylim(0.01, 0.3)
axs[3].set_ylim(0.1, 300)
label(axs[3], 15.1, 0.15, 'GALEX needs C IV suppressed ≥3–5× (see fuv figure)', color=FG, size=7)
axs[0].legend(loc='upper right', fontsize=7)
finish(fig, ROOT + '/figures/verify_shockmodel_tracks.png')

# ------------------------------------------------------------------ figure 2: closure maps over (v, b)
acc = an['acceptable']
vs = sorted(set(o['v'] for o in acc))
bs = sorted(set(o['b'] for o in acc))


def grid(fn):
    Z = np.full((len(bs), len(vs)), np.nan)
    for o in acc:
        i, j = bs.index(o['b']), vs.index(o['v'])
        val = fn(o)
        # keep the best-chi2 preshock state per cell
        if np.isnan(Z[i, j]) or o['chi2'] < grid.best.get((i, j), 1e9):
            Z[i, j] = val
            grid.best[(i, j)] = o['chi2']
    return Z


def fresh(fn):
    grid.best = {}
    return grid(fn)


fig, axs = new_sky_figure(13, 8.6, 'Silk · truncated-shock closure: n0 fixed by the observed 0.28° width (R_c = 58 pc, d = 280 pc)',
                          'cells = (v_s, b) with acceptable line ratios (best preshock X); hatched = consistent on all five tests',
                          nrows=2, ncols=2)
axs = axs.ravel()
panels = [(lambda o: o['W']['I_pk'], 'predicted [O III] ridge peak (R)  [obs 1.2–2.3, ≤4]', (0.1, 50), 'log'),
          (lambda o: o['W']['P_ram_k'], 'implied ρ0 v² / k (K cm⁻³)  [Monogem 3e4–4e5]', (3e3, 2e6), 'log'),
          (lambda o: o['W']['line_fwhm'], 'predicted [O III] FWHM (km/s)  [obs ≤ 33 (26–38)]', (10, 80), 'lin'),
          (lambda o: o['W']['t_star_yr'], 'time since shock entered the gas t* (yr)  [age 1.1e5]', (1e4, 2e6), 'log')]
import matplotlib.colors as mcolors
cons = {(o['v'], o['b']) for o in acc if all(o.get('consistency', {}).values()) and o.get('consistency')}
# recompute consistency flags here (analysis json has none; summary json has the list)
cons = {(c['v'], c['b']) for c in summ['consistent']}
for ax, (fn, ttl, lim, sc) in zip(axs, panels):
    Z = fresh(fn)
    norm = mcolors.LogNorm(*lim) if sc == 'log' else mcolors.Normalize(*lim)
    X, Y = np.meshgrid(np.arange(len(vs) + 1) - 0.5, np.arange(len(bs) + 1) - 0.5)
    m = ax.pcolormesh(X, Y, Z, norm=norm, cmap='viridis')
    for (vv, bb) in cons:
        j, i = vs.index(vv), bs.index(bb)
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch='///', ec='#ffffff', lw=0.0))
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=SUN, lw=1.2))
    for i in range(len(bs)):
        for j in range(len(vs)):
            if np.isfinite(Z[i, j]):
                txt = f'{Z[i, j]:.0f}' if sc == 'lin' else (f'{Z[i, j]:.1f}' if Z[i, j] < 10 else f'{Z[i, j]:.0e}'.replace('e+0', 'e'))
                ax.text(j, i, txt, ha='center', va='center', fontsize=5.5, color='w', path_effects=[HALO])
    ax.set_xticks(range(len(vs))); ax.set_xticklabels([str(int(v)) for v in vs], fontsize=7)
    ax.set_yticks(range(len(bs))); ax.set_yticklabels([f'{b:g}' for b in bs], fontsize=7)
    ax.set_xlabel('shock velocity v_s (km/s)')
    ax.set_ylabel('b = B0/√n0 (µG cm^1.5)')
    ax.set_title(ttl, fontsize=8.5, color=FG, loc='left')
    plt.colorbar(m, ax=ax, fraction=0.04, pad=0.01)
finish(fig, ROOT + '/figures/verify_shockmodel_closure.png')

# ------------------------------------------------------------------ figure 3: FUV + profiles
fig, axs = new_sky_figure(13, 5.6, 'Silk · GALEX FUV per R(5007) and predicted cross-filament / line profiles',
                          'FUV = C IV 1548/51 + O III] 1661/66 + He II 1640 only (N IV], Si IV, C II, 2γ omitted → lower bound before scattering)',
                          nrows=1, ncols=3)
ax = axs[0]
for o in acc:
    th, rt = o['FUV_RT']['total_thin_CU_per_R'], o['FUV_RT']['total_CU_per_R']
    c = SUN if (o['v'], o['b']) in cons else '#5a6570'
    ax.plot(o['v'] - 1.5, th, 'o', ms=3, color=c, alpha=0.8, mfc='none')
    ax.plot(o['v'] + 1.5, rt, 'o', ms=3, color=c, alpha=0.9)
ax.axhspan(92 - 10, 92 + 10, color=OIII_C, alpha=0.25, lw=0)
label(ax, 72, 98, 'GALEX south 92±10 CU/R (xrayfuv)', color=OIII_C, size=7)
ax.set_yscale('log'); ax.set_ylim(10, 1000)
ax.set_xlabel('v_s (km/s)'); ax.set_ylabel('FUV CU per R of [O III]5007')
ax.plot([], [], 'o', mfc='none', color=FGDIM, label='optically thin'); ax.plot([], [], 'o', color=FGDIM, label='edge-on C IV scattering (MC)')
ax.plot([], [], 'o', color=SUN, label='consistent models'); ax.legend(fontsize=7, loc='upper left')
ax.grid(alpha=0.3)
# best consistent model profiles
best = min(summ['consistent'], key=lambda c: c['chi2']) if summ['consistent'] else None
if best:
    tag = f"v{best['v']}_x{best['xn']}_b{best['b']}_asplund09_trn"
    rec, t = runs[tag], tr(tag)
    n0 = best['n0']
    Rc = 280 * np.radians(11.9) * AN.PC
    pr = AN.profile(t, rec['best']['i'], n0, fit['n0_ref'], Rc)
    hdeg = np.degrees(pr['hs'] / (280 * AN.PC))
    axs[1].plot(-hdeg, pr['I'], color=OIII_C, lw=1.6, label=f"model {best['v']} km/s, n0={n0:.3f}")
    prC = AN.profile(t, rec['best']['i'], n0, fit['n0_ref'], Rc, key='F_CIV', lam=1549.5)
    axs[1].plot(-np.degrees(prC['hs'] / (280 * AN.PC)), prC['I'] * best['S1548'], color=WARN, lw=1.2,
                label=f"C IV (×S={best['S1548']:.2f} scattering)")
    axs[1].axvspan(-0.28, 0, color='#9aa4ae', alpha=0.12)
    label(axs[1], -0.27, max(pr['I']) * 0.9, 'obs FWHM 0.28°', color=FG, size=7)
    axs[1].set_xlim(-1.2, 0.05)
    axs[1].set_xlabel('offset from shock front (deg, − = toward Monogem)')
    axs[1].set_ylabel('surface brightness (R)')
    axs[1].legend(fontsize=7, loc='upper left'); axs[1].grid(alpha=0.3)
    lp = AN.line_profile(t, rec['best']['i'], n0, fit['n0_ref'], Rc, best['v'], pr['h_pk'])
    axs[2].plot(lp['v'], lp['p'], color=OIII_C, lw=1.6, label=f"intrinsic, FWHM {lp['fwhm']:.0f} km/s")
    lsf = 98.0 / 2.3548
    conv = np.convolve(lp['p'], np.exp(-0.5 * (lp['v'] / lsf)**2), mode='same')
    axs[2].plot(lp['v'], conv / conv.max(), color=FGDIM, lw=1.0, ls='--', label='⊗ DESI LSF (98 km/s)')
    axs[2].axvspan(-16.5, 16.5, color='#9aa4ae', alpha=0.15)
    label(axs[2], -60, 0.95, 'obs intrinsic FWHM ≤ 33', color=FG, size=7)
    axs[2].set_xlim(-90, 90); axs[2].set_xlabel('V − V(preshock gas) (km/s)'); axs[2].set_ylabel('normalised')
    axs[2].legend(fontsize=7, loc='lower center'); axs[2].grid(alpha=0.3)
finish(fig, ROOT + '/figures/verify_shockmodel_fuv.png')
print('figures written')
