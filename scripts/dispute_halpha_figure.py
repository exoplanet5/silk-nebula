"""dispute_halpha_figure.py -- Figure 4 of the Comment on Wu et al. 2026 (AJ 172, 298; arXiv:2610.13329v1).

(a) Forward model for region P1 at the adopted joint (centroid + width, self-consistent sky/nebula) fit Fsc: WIM
    Halpha, true geocoronal Halpha,
    the subtracted MRS-N sky template, the residual spectrum and the fitted Gaussian, in V_helio.
(b) Fitted Halpha centroid shift vs geocoronal over-subtraction (model/true - 1) for P1-P4, with the
    observed Table 4 Halpha - [N II] offsets.
(c) Fitted intrinsic Halpha FWHM vs over-subtraction, with the Table 2 Halpha widths and WHAM widths.
(d) Per-region velocities: LAMOST Halpha, [N II], [S II] (Table 4), WHAM-SS Halpha, model-corrected Halpha.

Input: data/derived/dispute_halpha_results.json (scripts/dispute_halpha_analysis.py),
       data/derived/verify_velocity_summary.json (adopted [O III] velocity, Sect. 3 of the Comment).
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_halpha_figure.py
Output: figures/dispute_halpha_fig4.png (dpi 200)
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
from silk_style import (plt, new_sky_figure, label, finish, FG, FGDIM, BG2, HA_C, NII_C, SII_C, OK, ACCENT,  # noqa: E402
                        SUN, OIII_C)
from matplotlib.gridspec import GridSpec  # noqa: E402

D = json.load(open(ROOT / 'data/derived/dispute_halpha_results.json'))
VS = json.load(open(ROOT / 'data/derived/verify_velocity_summary.json'))
REGS = ['P1', 'P2', 'P3', 'P4']
RCOL = {'P1': '#80deea', 'P2': '#f48fb1', 'P3': '#ffe082', 'P4': '#b0bec5'}
GRID_KW = dict(color='#2a323c', lw=0.6, alpha=0.6)


def panel_title(ax, text):
    ax.set_title(text, loc='left', color=FG, fontsize=10.5, pad=6)


def main():
    fig = plt.figure(figsize=(16, 8.6))
    fig.text(0.012, 0.985, 'LAMOST Hα offset: a plausible geocoronal-subtraction residual', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.952, 'Wu et al. 2026 Tables 2, 4 · WHAM-SS DR1 · MRS-N sky model (Zhang et al. 2021) · '
                           '2024-02-01 plate · forward model, not a fit to spectra', color=FGDIM, fontsize=9.5, va='top')
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.12, 1.0, 1.0], left=0.045, right=0.985, top=0.885, bottom=0.085,
                  wspace=0.27, hspace=0.34)
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, 1], sharex=axb)
    axd = fig.add_subplot(gs[:, 2])

    # ------------------------------------------------------------------ (a) profiles, P1, joint fit
    a = D['fig4a']
    v = np.array(a['v'])
    norm = np.max(a['wim'])
    wim, geo, tm = (np.array(a[k]) / norm for k in ('wim', 'geo', 'tmpl'))
    res, fitl = np.array(a['resid']) / norm, np.array(a['fit']) / norm
    axa.fill_between(v, 0, wim, color=HA_C, alpha=0.18, lw=0)
    axa.plot(v, wim, color=HA_C, lw=1.6, label='WIM Hα (true)')
    axa.plot(v, geo, color=SUN, lw=1.6, label='geocoronal Hα (true)')
    axa.plot(v, tm, color=ACCENT, lw=1.6, ls='--', label='sky template subtracted')
    axa.plot(v, res, color=FG, lw=2.2, label='residual (what is fitted)')
    axa.plot(v, fitl, color=OK, lw=1.4, ls=(0, (4, 2)), label='Gaussian + quadratic fit')
    axa.axhline(0, color=FGDIM, lw=0.6)
    # centroid ticks below the zero line, labels next to each curve's peak
    for x, c in ((a['geo_centre'], SUN), (a['tmpl_centre'], ACCENT), (a['wim_centre'], HA_C), (a['fit_centre'], OK)):
        axa.plot([x, x], [-0.15, -0.03], color=c, lw=2.0, solid_capstyle='butt')
    shift = a['fit_centre'] - a['wim_centre']
    label(axa, -21, 1.04, f"geocoronal\n{a['geo_centre']:+.1f}", color=SUN, size=8.4, ha='right', va='center')
    label(axa, -1, 1.22, f"template\n{a['tmpl_centre']:+.1f}", color=ACCENT, size=8.4, ha='left', va='center')
    label(axa, 22, 1.08, f"true WIM\n{a['wim_centre']:+.1f}", color=HA_C, size=8.4, ha='left', va='center')
    label(axa, 36, 0.80, f"fitted {a['fit_centre']:+.1f}\n(shift {shift:+.1f} km/s)", color=OK, size=8.4,
          ha='left', va='center')
    over = 100 * (1 / a['kt'] - 1)
    label(axa, -72, -0.30, f"P1 · subtracted sky = {a['g_selfcons']:.2f}× true WIM Hα\n"
                           f"template {over:.0f}% too bright, at RV(OH 6554) {a['dm']:+.1f} km/s\n"
                           f"ticks: centroids (V$_{{helio}}$, km/s)",
          color=FGDIM, size=8.2, ha='left', va='center')
    axa.set_xlim(-75, 80)
    axa.set_ylim(-0.48, 1.92)
    axa.set_xlabel('V$_{helio}$ (km/s)')
    axa.set_ylabel('flux (WIM peak = 1, LSF-convolved)')
    axa.grid(**GRID_KW)
    axa.legend(loc='upper right', fontsize=8, framealpha=0.9)
    panel_title(axa, '(a) What the fit sees after sky subtraction')

    # ------------------------------------------------------------------ (b) centroid shift
    b = D['fig4b']
    f = 100 * np.array(b['oversub_frac'])
    off = D['offsets_table4']
    F = D['fits']['Fsc_joint_selfcons']
    B = D['fits']['Bsc_centroid_only_selfcons']
    nan = lambda a_: np.array([np.nan if x is None else x for x in a_], float)   # unconverged points -> gap
    lo, hi = F['oversub_pct_1sigma']
    for ax in (axb, axc):
        ax.axvspan(-30, 0, color='#ffffff', alpha=0.035, lw=0)
        ax.axvspan(lo, hi, color=OK, alpha=0.16, lw=0)
        ax.axvline(F['oversub_pct'], color=OK, lw=1.1)
        ax.axvline(B['oversub_pct'], color=FGDIM, lw=1.0, ls=':')
        ax.axvline(0, color=FGDIM, lw=0.6)
        ax.grid(**GRID_KW)
    mobs, eobs = off['Ha_minus_NII']['wmean'], off['Ha_minus_NII']['err']
    vals = [off['per_region'][r]['Ha_minus_NII'] for r in REGS]
    axb.axhspan(min(vals), max(vals), color=HA_C, alpha=0.10, lw=0)
    axb.axhspan(mobs - eobs, mobs + eobs, color=HA_C, alpha=0.22, lw=0)
    axb.axhline(mobs, color=HA_C, lw=1.0)
    for r in REGS:
        axb.plot(f, nan(b['shift_F'][r]), color=RCOL[r], lw=1.6, label=r)
    axb.axhline(0, color=FGDIM, lw=0.6)
    label(axb, -28, max(vals) + 0.3, f'observed Hα − [N II]\n{mobs:+.1f} ± {eobs:.1f} km/s', color=HA_C, size=8.2,
          ha='left', va='bottom')
    label(axb, -28, 1.6, 'under-\nsubtracted\n→ blueshift', color=FGDIM, size=8, ha='left', va='center')
    label(axb, 78, 11.6, 'over-subtracted\n→ redshift', color=FGDIM, size=8, ha='right', va='top')
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch as _Patch
    hb = [(_Patch(facecolor=OK, alpha=0.3, lw=0), Line2D([], [], color=OK, lw=1.1)),
          Line2D([], [], color=FGDIM, lw=1.0, ls=':')]
    lb = [f"joint fit {F['oversub_pct']:.0f}% ({lo:.0f}–{hi:.0f}%), centroid + width",
          f"centroid-only fit {B['oversub_pct']:.0f}%"]
    axb.legend(hb, lb, loc='lower right', fontsize=7.8, framealpha=0.92, handlelength=1.8)
    axb.set_ylim(-14, 12)
    axb.set_ylabel('Hα centroid shift (km/s)')
    plt.setp(axb.get_xticklabels(), visible=False)
    panel_title(axb, '(b) Centroid shift vs over-subtraction')

    # ------------------------------------------------------------------ (c) fitted width
    t2 = D['fits']['table2_Ha_fwhm']
    wv = np.array([t2[r][0] for r in REGS]); we = np.array([t2[r][1] for r in REGS])
    wm = np.sum(wv / we ** 2) / np.sum(1 / we ** 2); wme = np.sqrt(1 / np.sum(1 / we ** 2))
    whw = [D['wham_per_region'][r]['fwhm_int_nearest'] for r in REGS]
    axc.axhspan(wm - wme, wm + wme, color=HA_C, alpha=0.22, lw=0)
    axc.axhline(wm, color=HA_C, lw=1.0)
    axc.axhspan(min(whw), max(whw), facecolor='none', edgecolor=OK, hatch='////', lw=0, alpha=0.5)
    for r in REGS:
        axc.plot(f, nan(b['fwhm_F'][r]), color=RCOL[r], lw=1.6, label=r)
    label(axc, -24, wm + wme + 0.6, f'Table 2 Hα {wm:.1f} ± {wme:.1f}', color=HA_C, size=8.3, ha='left',
          va='bottom')
    label(axc, -28, min(whw) - 0.8, f'WHAM WIM {min(whw):.0f}–{max(whw):.0f}\n(true width, model input)', color=OK,
          size=8.3, ha='left', va='top')
    axc.set_ylim(0, 68)
    axc.set_xlim(-30, 80)
    axc.legend(loc='upper right', ncol=2, fontsize=7.8, framealpha=0.9, handlelength=1.4, columnspacing=0.9,
               title='region', title_fontsize=7.8)
    axc.set_xlabel('geocoronal over-subtraction, model/true − 1 (%)')
    axc.set_ylabel('fitted intrinsic Hα FWHM (km/s)')
    panel_title(axc, '(c) Fitted Hα width vs over-subtraction')

    # ------------------------------------------------------------------ (d) per-region velocities
    T4 = {'P1': (16.3, 1.2, 11.8, 1.6, 11.8, 2.5), 'P2': (17.0, 1.2, 12.0, 1.6, 12.0, 1.6),
          'P3': (18.4, 1.4, 11.7, 1.9, 11.7, 1.9), 'P4': (19.7, 1.4, 13.3, 2.3, 12.7, 2.1)}   # Table 4 V_helio
    o3 = VS['combined']['V_helio']; o3e = VS['combined']['err']
    axd.axvspan(o3 - o3e, o3 + o3e, color=OIII_C, alpha=0.13, lw=0)
    axd.axvline(o3, color=OIII_C, lw=0.8, alpha=0.8)
    yl = []
    for i, r in enumerate(REGS):
        y0 = -i
        ha, hae, n2, n2e, s2, s2e = T4[r]
        wh = D['wham_per_region'][r]
        corr = ha - F['pred_shift'][r]
        rows = [(ha, hae, 'o', HA_C, True, 'LAMOST Hα (Table 4)'),
                (corr, hae, 'o', HA_C, False, 'Hα − model shift (joint fit)'),
                (n2, n2e, 'D', NII_C, True, '[N II] (Table 4)'),
                (s2, s2e, 'X', SII_C, True, '[S II] (Table 4)'),
                (wh['v_helio_nearest'], max(wh['v_lsr_within1deg_std'], 0.5), 's', OK, True, 'WHAM Hα (nearest beam)')]
        for k, (x, e, mk, c, filled, lab) in enumerate(rows):
            yy = y0 + 0.30 - 0.15 * k
            axd.errorbar(x, yy, xerr=e, fmt=mk, ms=6.5, color=c, mfc=c if filled else 'none', mec=c, ecolor=c,
                         elinewidth=1.1, capsize=2.5, label=lab if i == 0 else None, zorder=3)
        label(axd, 24.7, y0 + 0.43, f"Hα − [N II] {ha - n2:+.1f}", color=HA_C, size=8.2, ha='right', va='center')
        yl.append(y0)
    axd.set_yticks(yl)
    axd.set_yticklabels([r + ('\n(no nebula)' if r == 'P4' else '') for r in REGS])
    for i in range(1, 4):
        axd.axhline(-i + 0.5, color='#2a323c', lw=0.8)
    axd.set_xlim(5, 25)
    axd.set_ylim(-4.3, 0.62)
    axd.set_xlabel('V$_{helio}$ (km/s)')
    axd.grid(axis='x', **GRID_KW)
    sec = axd.secondary_xaxis('top', functions=(lambda x: x - 12.2, lambda x: x + 12.2))
    sec.set_xlabel('V$_{LSR}$ (km/s, ≈ V$_{helio}$ − 12.2)', color=FGDIM, fontsize=8.5)
    sec.tick_params(colors=FGDIM, labelsize=8)
    from matplotlib.patches import Patch
    h, l = axd.get_legend_handles_labels()
    h.append(Patch(facecolor=OIII_C, alpha=0.3, edgecolor=OIII_C)); l.append(f'[O III] {o3:+.1f} ± {o3e:.1f} (Sect. 3)')
    axd.legend(h, l, loc='lower right', fontsize=7.8, framealpha=0.92, bbox_to_anchor=(1.0, 0.0), ncol=1)
    panel_title(axd, '(d) Velocities per region')

    finish(fig, str(ROOT / 'figures/dispute_halpha_fig4.png'), dpi=200, tight=False)


if __name__ == '__main__':
    main()
