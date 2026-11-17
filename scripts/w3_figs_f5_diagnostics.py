"""F5 excitation diagnostics: observed Silk constraints vs model loci (energy ratios).

(a) [O III]5007/Hb vs [O II]3727/[O III]5007;  (b) [O III]4363/5007 vs [O III]5007/Hb (right axis: T_e).
Models (own wave-2 codes, read from disk):
  verify_shockmodel_tracks.npz / _fit.json : truncated NEI shock tracks (cumulative ratios vs swept column; v = 100-160,
      X(H0) = 0.05, b = 10), best truncation (diamond), complete-shock end points (b = 3; triangles)
  Raymond+88 SE / M incomplete models (verify_xrayfuv_uv_predictions.json; 4363 for SE from verify_xrayfuv.md)
  verify_altmodel_photoion.json : one-zone PIE curves vs log U (81 kK WD, 150 kK, 0.14 keV hot gas, 35 kK DIG, 28 eV)
  verify_altmodel_cooling.json  : GS07 NEI cooling from 1.6 MK truncated at 5,3,2,1e4 K (+H-recombination tail);
      warm CIE layer vs T; static / evaporating conduction fronts
  verify_altmodel_scenarios.json: PIE + extra heating to 4e4 K; recombining fossil
Observed (energy units; wave-2): DESI south [O II]/[O III] 1.11-1.16 (+-0.08), all 1.37-1.41 (1.21-2.07), north 4.2-4.4
  (3.5-5.2); [O III]/Hb > 11 (box) / > 18 (Gaussian) at 2 sigma; 4363/5007 0.06-0.09 (0.083 +- 0.029).
  SDSS south (Dec<20; verify_sdss_results_base_r0.15.json, photon -> energy): [O III]/Hb 22.8 (12.4-46.2),
  > 7.3 (2 sigma with placebo); [O II]/[O III] 1.84 (1.46-2.37); 4363/5007 0.038 +- 0.035.
  SDSS north: [O III]/Hb 6.2 (3.1-12.0), [O II]/[O III] 7.3 (6.1-8.8). Ha patch (SDSS 30 % contour): 1.73 / 11.8.
Output: figures/w3_figs_f5_diagnostics.png
"""
import sys
import json
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, panel_tag

PH2E_HB = 4861.3 / 5006.8      # [O III]5007/Hb photon -> energy
PH2E_43 = 5006.8 / 4363.2      # 4363/5007 photon -> energy


def r43(T):
    return 1.3356 / (7.90 * np.exp(32900.0 / T))


def T_of_r(r):
    r = np.asarray(r, float)
    with np.errstate(all='ignore'):
        return 32900.0 / np.log(1.3356 / (7.90 * r))


def main():
    TR = np.load(DER / 'verify_shockmodel_tracks.npz')
    fit = json.loads((DER / 'verify_shockmodel_fit.json').read_text())
    runs = {r['tag']: r for r in fit['runs']}
    pio = json.loads((DER / 'verify_altmodel_photoion.json').read_text())
    coo = json.loads((DER / 'verify_altmodel_cooling.json').read_text())
    scn = json.loads((DER / 'verify_altmodel_scenarios.json').read_text())['scenarios']
    sd = json.loads((DER / 'verify_sdss_results_base_r0.15.json').read_text())['regions']

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.2, 7.4), gridspec_kw=dict(wspace=0.30))
    fig.text(0.012, 0.985, 'Silk Nebula · excitation diagnostics vs models', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.952, 'energy ratios · observed south (DESI Dec < 19.5, SDSS Dec < 20) and north · '
             'model loci from the wave-2 shock and alternative-mechanism codes', color=FGDIM, fontsize=9.3, va='top')

    # ---------------- models
    shock_cols = {100: '#7ad97a', 120: OIII_C, 140: '#c792ea', 160: WARN}
    for v, c in shock_cols.items():
        tag = f'v{v}_x0.05_b10.0_asplund09_trn'
        t = {k.split('__', 1)[1]: TR[k] for k in TR.files if k.startswith(tag + '__')}
        lN = np.log10(t['N'])
        m = t['F_O3'] >= 0.2 * t['F_O3'][-1]     # show the track once >= 20 % of the final [O III] has formed
        o3hb = t['F_O3'] / t['F_Hb']; o2o3 = t['F_O2'] / t['F_O3']; r4363 = t['F_O4363'] / t['F_O3']
        axA.plot(o2o3[m], o3hb[m], color=c, lw=1.6, alpha=0.95, zorder=3)
        axB.plot(o3hb[m], r4363[m], color=c, lw=1.6, alpha=0.95, zorder=3)
        i = runs[tag]['best']['i']
        axA.plot(o2o3[i], o3hb[i], 'D', color=c, ms=6, mec='k', zorder=6)
        axB.plot(o3hb[i], r4363[i], 'D', color=c, ms=6, mec='k', zorder=6)
        axA.plot(o2o3[m][-1], o3hb[m][-1], 'o', color=c, ms=3.5, zorder=5)
    # R88-calibrated Hb shifts the shock [O III]/Hb down by 2.3 (verify_shockmodel)
    t = {k.split('__', 1)[1]: TR[k] for k in TR.files if k.startswith('v120_x0.05_b10.0_asplund09_trn__')}
    i = runs['v120_x0.05_b10.0_asplund09_trn']['best']['i']
    xx, yy = (t['F_O2'] / t['F_O3'])[i], (t['F_O3'] / t['F_Hb'])[i]
    axA.annotate('', xy=(xx * 0.93, yy / 2.3), xytext=(xx * 0.93, yy * 0.92),
                 arrowprops=dict(arrowstyle='-|>', color=FGDIM, lw=0.9, mutation_scale=8), zorder=4)
    label(axA, xx * 0.88, yy / 1.5, 'Hβ ×2.3\n(R88 cal.)', color=FGDIM, size=6.5, ha='right', va='center')
    # complete-shock end points
    for tag, rr in runs.items():
        if not rr['complete']:
            continue
        t = {k.split('__', 1)[1]: TR[k][-1] for k in TR.files if k.startswith(tag + '__')}
        c = shock_cols.get(rr['v'], FGDIM)
        axA.plot(t['F_O2'] / t['F_O3'], t['F_O3'] / t['F_Hb'], '^', color=c, ms=6, mec='k', mew=0.5, zorder=5)
        axB.plot(t['F_O3'] / t['F_Hb'], t['F_O4363'] / t['F_O3'], '^', color=c, ms=6, mec='k', mew=0.5, zorder=5)
    axA.axhspan(1.2, 4.5, color=DANGER, alpha=0.08, lw=0, zorder=0)
    label(axA, 0.06, 1.3, 'published complete shocks: [O III]/Hβ 1.2–4.5', color=DANGER, size=7, ha='left',
          va='bottom')
    # Raymond+88 incomplete
    for nm, x_, y_, r_ in (('R88 SE', 1.341, 12.9, 0.071), ('R88 M', 1.638, 10.5, None)):
        axA.plot(x_, y_, 's', color=SUN, ms=6, mec='k', zorder=6)
        if nm == 'R88 M':
            label(axA, x_ * 1.07, y_, 'R88 SE, M', color=SUN, size=7, ha='left', va='center')
        if r_ is not None:
            axB.plot(y_, r_, 's', color=SUN, ms=6, mec='k', zorder=6)
            label(axB, y_ * 0.93, r_ * 1.0, nm, color=SUN, size=7, ha='right', va='center')
    # photoionisation curves
    pcols = {'BB81kK_WDJ0736': ('PIE · WD J0736 (81 kK)', '#9aa4ae', '-'), 'BB150kK': ('PIE · 150 kK', '#9aa4ae', '--'),
             'hotgas_kT140eV': ('PIE · 0.14 keV hot gas', '#9aa4ae', ':'), 'BB35kK_DIG': ('PIE · DIG (35 kK)', '#6b7480', '-.')}
    for k, (lab, c, ls) in pcols.items():
        rows = pio['models'][k]
        o3hb = np.array([r['OIII_Hb'] for r in rows]); o2o3 = np.array([r['OII_OIII'] for r in rows])
        r43_ = np.array([r['O4363_5007'] for r in rows])
        ok = (o3hb > 0.05) & (o2o3 < 200)
        axA.plot(o2o3[ok], o3hb[ok], color=c, lw=1.1, ls=ls, zorder=2)
        axB.plot(o3hb[ok], r43_[ok], color=c, lw=1.1, ls=ls, zorder=2)
    # GS07 NEI cooling flows truncated at decreasing T (connected), plus the complete (tail) point
    for fam, c, mk in (('isobaric_Z1', '#ff9e80', 'o'), ('isochoric_Z1', '#ff9e80', 's')):
        F = coo['cooling_flows'][fam]
        seq = ['T1.6e+06-5e+04', 'T1.6e+06-3e+04', 'T1.6e+06-2e+04', 'T1.6e+06-1e+04', 'T1.6e+06-1e+04+tail']
        xs = [F[k]['OII_OIII'] for k in seq]; ys = [F[k]['OIII_Hb'] for k in seq]; zs = [F[k]['O4363_5007'] for k in seq]
        axA.plot(xs, ys, color=c, lw=0.9, ls='--', marker=mk, ms=3.5, zorder=3)
        axB.plot(ys, zs, color=c, lw=0.9, ls='--', marker=mk, ms=3.5, zorder=3)
    # warm CIE vs T
    W = coo['warm_CIE']['own'][:len(coo['warm_CIE']['T'])]
    wx = np.array([r['OII_OIII'] for r in W]); wy = np.array([r['OIII_Hb'] for r in W]); wz = np.array([r['O4363_5007'] for r in W])
    ok = np.isfinite(wx) & np.isfinite(wy) & (wy > 0.3) & (wx < 100) & (wx > 0.02)
    axA.plot(wx[ok], wy[ok], color='#e8eaed', lw=0.9, ls=(0, (1, 1.5)), zorder=2)
    axB.plot(wy[ok], wz[ok], color='#e8eaed', lw=0.9, ls=(0, (1, 1.5)), zorder=2)
    # conduction fronts, PIE + heating, fossil
    cf = [coo['static_fronts']['P1.5e+04']] + [coo['evaporation']['runs'][k] for k in coo['evaporation']['runs']]
    axA.plot([c['OII_OIII'] for c in cf], [c['OIII_Hb'] for c in cf], 'X', color='#ff5252', ms=6.5, mec='k', zorder=6)
    axB.plot([c['OIII_Hb'] for c in cf], [c['O4363_5007'] for c in cf], 'X', color='#ff5252', ms=6.5, mec='k', zorder=6)
    for k_, mk_, c_ in (('4d PIE + extra heating to 4e4 K', 'P', '#9aa4ae'),
                        ('6 fossil, recombining (hot start 4e4 K, n_H = 100)', 'v', '#9aa4ae')):
        p = scn[k_]['pred']
        axA.plot(p['OII_OIII'], p['OIII_Hb'], mk_, color=c_, ms=7, mec='k', zorder=6)
        axB.plot(p['OIII_Hb'], p['O4363'], mk_, color=c_, ms=7, mec='k', zorder=6)

    # ---------------- observations
    def obsbox(ax, x0, x1, y0, y1, c, alpha=0.28, hatch=None):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=c, alpha=alpha, edgecolor=c, lw=1.2,
                               hatch=hatch, zorder=7))
    # DESI south / all: [O III]/Hb lower limits
    obsbox(axA, 1.03, 1.24, 11, 300, OIII_C, 0.30)
    axA.plot([1.13, 1.13], [18, 18], '_', color=OIII_C, ms=14, mew=2.2, zorder=8)
    axA.errorbar(1.39, 8.2, xerr=[[1.39 - 1.21], [2.07 - 1.39]], fmt='|', color=OIII_C, ms=8, lw=1.0, alpha=0.9,
                 zorder=7)
    label(axA, 2.2, 8.2, 'DESI all [O II]', color=OIII_C, size=7, ha='left', va='center')
    label(axA, 0.98, 125, 'DESI S\n[O III]/Hβ > 11–18', color=OIII_C, size=7.5, ha='right', va='center')
    # DESI north ([O III]/Hb not measured by DESI north; drawn at the SDSS-north range)
    s_ = sd['south_dec_lt20']['ratios_photon']
    yS = s_['OIII5007/Hb']['value'] * PH2E_HB
    yS16, yS84 = s_['OIII5007/Hb']['p16'] * PH2E_HB, s_['OIII5007/Hb']['p84'] * PH2E_HB
    xS = s_['OII/OIII5007_energy_BOSS']
    axA.errorbar(xS['value'], yS, xerr=[[xS['value'] - xS['p2.5']], [xS['p97.5'] - xS['value']]],
                 yerr=[[yS - yS16], [yS84 - yS]], fmt='o', color='#ffffff', ms=7, mec='k', lw=1.2, capsize=2, zorder=9)
    label(axA, xS['value'] * 1.1, yS * 1.18, 'SDSS S', color='#ffffff', size=7.5, ha='left', va='bottom')
    n_ = sd['north_dec_ge20']['ratios_photon']
    yN = n_['OIII5007/Hb']['value'] * PH2E_HB
    yN16, yN84 = n_['OIII5007/Hb']['p16'] * PH2E_HB, n_['OIII5007/Hb']['p84'] * PH2E_HB
    xN = n_['OII/OIII5007_energy_BOSS']
    axA.errorbar(xN['value'], yN, xerr=[[xN['value'] - xN['p2.5']], [xN['p97.5'] - xN['value']]],
                 yerr=[[yN - yN16], [yN84 - yN]], fmt='o', color='#b8ecff', mfc='none', ms=7, lw=1.1, capsize=2, zorder=9)
    label(axA, xN['value'] * 1.08, yN * 0.72, 'SDSS N', color='#b8ecff', size=7.5, ha='left', va='top')
    obsbox(axA, 3.5, 5.2, yN16, yN84, '#b8ecff', 0.18)
    label(axA, 4.25, yN84 * 1.08, 'DESI N [O II]', color='#b8ecff', size=7, ha='center', va='bottom')
    p_ = sd['halpha_patch_30pct']['ratios_photon']
    yP = p_['OIII5007/Hb']['value'] * PH2E_HB
    xP = p_['OII/OIII5007_energy_BOSS']['value']
    axA.plot(xP, yP, '*', color=HA_C, ms=12, mec='k', zorder=9)
    label(axA, xP * 1.0, yP * 0.78, 'Hα patch', color=HA_C, size=7.5, ha='center', va='top')

    # panel B observations
    obsbox(axB, 11, 300, 0.057, 0.091, OIII_C, 0.30)
    axB.plot([11, 300], [0.083, 0.083], color=OIII_C, lw=1.0, ls='-', alpha=0.9, zorder=8)
    axB.plot([18, 18], [0.054, 0.112], color=OIII_C, lw=1.0, ls=':', zorder=8)
    label(axB, 290, 0.053, 'DESI S: 0.06–0.09 (2–3σ)', color=OIII_C, size=7.5, ha='right', va='top')
    y43 = 0.038
    axB.errorbar(yS, y43, xerr=[[yS - yS16], [yS84 - yS]], yerr=0.035 * 0.97, fmt='o', color='#ffffff', ms=7, mec='k',
                 lw=1.1, capsize=2, zorder=9)
    label(axB, yS * 1.1, y43 * 0.82, 'SDSS S', color='#ffffff', size=7.5, ha='left', va='top')

    # ---------------- axes
    axA.set_xscale('log'); axA.set_yscale('log')
    axA.set_xlim(0.05, 30); axA.set_ylim(0.3, 300)
    axA.set_xlabel('[O II] 3726+29 / [O III] 5007'); axA.set_ylabel('[O III] 5007 / Hβ')
    axB.set_xscale('log'); axB.set_yscale('log')
    axB.set_xlim(0.5, 300); axB.set_ylim(0.001, 0.25)
    axB.set_xlabel('[O III] 5007 / Hβ'); axB.set_ylabel('[O III] 4363 / 5007')
    for ax in (axA, axB):
        ax.grid(alpha=0.22, which='major')
        ax.tick_params(labelsize=8.5)
    sec = axB.twinx(); sec.set_yscale('log'); sec.set_ylim(axB.get_ylim())
    Ts = np.array([8e3, 1e4, 1.5e4, 2e4, 3e4, 5e4, 1e5])
    sec.set_yticks(r43(Ts)); sec.set_yticklabels(['8k', '10k', '15k', '20k', '30k', '50k', '100k'])
    sec.minorticks_off()
    sec.set_ylabel('T_e([O III]) (K, low-density)', color=FGDIM); sec.tick_params(labelsize=8, colors=FGDIM)
    panel_tag(axA, 'a  [O III]/Hβ vs [O II]/[O III]')
    panel_tag(axB, 'b  T_e indicator vs [O III]/Hβ')

    hand = [Line2D([], [], color=shock_cols[v], lw=1.6, label=f'truncated shock {v} km/s') for v in shock_cols] + \
           [Line2D([], [], marker='D', ls='', color=FGDIM, mec='k', label='best truncation'),
            Line2D([], [], marker='^', ls='', color=FGDIM, mec='k', label='complete shock (own code)'),
            Line2D([], [], marker='s', ls='', color=SUN, mec='k', label='Raymond+88 incomplete'),
            Line2D([], [], color='#9aa4ae', lw=1.1, ls='-', label='PIE · WD 81 kK (vs U)'),
            Line2D([], [], color='#9aa4ae', lw=1.1, ls='--', label='PIE · 150 kK'),
            Line2D([], [], color='#9aa4ae', lw=1.1, ls=':', label='PIE · 0.14 keV gas'),
            Line2D([], [], color='#6b7480', lw=1.1, ls='-.', label='PIE · DIG 35 kK'),
            Line2D([], [], color='#ff9e80', lw=0.9, ls='--', marker='o', ms=3, label='NEI cooling flow (GS07)'),
            Line2D([], [], color='#e8eaed', lw=0.9, ls=(0, (1, 1.5)), label='warm CIE layer (vs T)'),
            Line2D([], [], marker='X', ls='', color='#ff5252', mec='k', label='conduction fronts'),
            Line2D([], [], marker='P', ls='', color='#9aa4ae', mec='k', label='PIE + heating to 4e4 K'),
            Line2D([], [], marker='v', ls='', color='#9aa4ae', mec='k', label='recombining fossil')]
    fig.legend(handles=hand, loc='lower center', ncol=8, fontsize=7.3, frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, 0.03), columnspacing=1.2, handlelength=1.8)
    ss.credit(fig)
    fig.subplots_adjust(left=0.055, right=0.95, top=0.90, bottom=0.21)
    fig.savefig(FIG / 'w3_figs_f5_diagnostics.png', dpi=200)
    plt.close(fig)
    print('saved; SDSS S [O III]/Hb energy %.1f (%.1f-%.1f); patch %.2f / %.1f; north %.1f / %.1f'
          % (yS, yS16, yS84, yP, xP, yN, xN['value']))


if __name__ == '__main__':
    main()
