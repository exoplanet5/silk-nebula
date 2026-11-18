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

Paper version (comment manuscript): SILK_THEME=light SILK_PAPER=1 ... -> paper/figures/dispute_halpha_fig4.{png,pdf}
    light palette, no figure title/subtitle or credit line, 7.1 in wide layout with (a)+(d) on top and (b)+(c)
    below, print-size fonts (>= 6.5 pt) and line widths, "W26" wording. No annotation lies on a hatch, grid or
    reference line: annotations sit on solid boxes clear of the reference lines, the grid is blanked next to
    them, the Halpha - [N II] band of panel (b) is named in the legend, and the legends are opaque. Without the
    environment variables the original dark figure is reproduced unchanged.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import silk_style as ss  # noqa: E402
from silk_style import (plt, new_sky_figure, label, finish, FG, FGDIM, BG2, HA_C, NII_C, SII_C, OK, ACCENT,  # noqa: E402
                        SUN, OIII_C)
from matplotlib.gridspec import GridSpec  # noqa: E402

PAPER, LIGHT = ss.PAPER, ss.LIGHT

D = json.load(open(ROOT / 'data/derived/dispute_halpha_results.json'))
VS = json.load(open(ROOT / 'data/derived/verify_velocity_summary.json'))
REGS = ['P1', 'P2', 'P3', 'P4']
if LIGHT:
    # saturated region colours that read on white (the dark-theme pastels vanish there)
    RCOL = {'P1': '#0097a7', 'P2': '#c2185b', 'P3': '#c08a00', 'P4': '#546e7a'}
    GRID_KW = dict(color=ss.GRID, lw=0.5, alpha=1.0)
    SHADE_C, SHADE_A = ss.FG, 0.045          # "under-subtracted" side of panels (b), (c)
    SEP_C = '#b9c1c9'                        # region separators in panel (d)
else:
    RCOL = {'P1': '#80deea', 'P2': '#f48fb1', 'P3': '#ffe082', 'P4': '#b0bec5'}
    GRID_KW = dict(color='#2a323c', lw=0.6, alpha=0.6)
    SHADE_C, SHADE_A = '#ffffff', 0.035
    SEP_C = '#2a323c'


def pk(dark, paper):
    """Value for the original (dark, screen) figure or for the print-size paper figure."""
    return paper if PAPER else dark


if PAPER:
    # the figure is drawn at its print size (7.1 in wide), so point sizes below are the printed sizes
    plt.rcParams.update({'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.labelsize': 7.4,
                         'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'xtick.major.size': 2.6, 'ytick.major.size': 2.6, 'xtick.major.pad': 2.0,
                         'ytick.major.pad': 2.0, 'axes.labelpad': 2.5, 'hatch.linewidth': 0.6,
                         'legend.borderpad': 0.35, 'legend.labelspacing': 0.3, 'legend.handletextpad': 0.5,
                         'legend.borderaxespad': 0.4})

FS_T = pk(10.5, 8.0)     # panel titles
FS_A = pk(8.4, 6.9)      # annotations
FS_L = pk(7.8, 6.6)      # legends
W26 = 'W26 ' if PAPER else ''
# Paper version: annotations get a solid background box, so that no grid line, reference line or hatch shows
# through the gaps between their words (at print size such a line reads like an underscore or an extra stroke).
# BOX is white; BOX_SH has the colour of the grey "under-subtracted" strip of panels (b), (c), for labels inside it.
if PAPER:
    from matplotlib.colors import to_rgb
    _SH = tuple((1 - SHADE_A) * w + SHADE_A * c for w, c in zip(to_rgb(ss.AXBG), to_rgb(SHADE_C)))
    BOX = dict(bbox=dict(boxstyle='square,pad=0.15', facecolor=ss.AXBG, edgecolor='none'))
    BOX_SH = dict(bbox=dict(boxstyle='square,pad=0.15', facecolor=_SH, edgecolor='none'))
else:
    BOX = BOX_SH = {}


def sgn(x):
    """Signed one-decimal number; typographic minus sign in the paper version."""
    t = f'{x:+.1f}'
    return t.replace('-', '−') if PAPER else t


def panel_title(ax, text):
    ax.set_title(text, loc='left', color=FG, fontsize=FS_T, pad=pk(6, 4))


def blank_sectors(fig, jobs):
    """Paper version only. jobs: (ax, boxed annotation, x0, x1, colour, top or None).

    Blanks the grid lines between x0 and x1 (the frame or a vertical reference line) over the height of the
    annotation (up to `top` if given), so that no stub of a grid line is left between the text box and that line
    (next to a word such a stub reads like a dash, in front of an arrow like a longer shaft)."""
    from matplotlib.patches import Rectangle
    fig.canvas.draw()
    for ax, t, x0, x1, col, top in jobs:
        bb = t.get_bbox_patch().get_window_extent()
        y0, y1 = ax.transData.inverted().transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])[:, 1]
        y1 = y1 if top is None else top
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=col, edgecolor='none', zorder=1.7))


def main():
    if PAPER:
        fig = plt.figure(figsize=(7.1, 6.3))
        top = GridSpec(1, 2, figure=fig, width_ratios=[1.22, 1.0], left=0.077, right=0.985, top=0.918,
                       bottom=0.455, wspace=0.30)
        bot = GridSpec(1, 2, figure=fig, left=0.077, right=0.985, top=0.375, bottom=0.075, wspace=0.22)
        axa = fig.add_subplot(top[0, 0])
        axd = fig.add_subplot(top[0, 1])
        axb = fig.add_subplot(bot[0, 0])
        axc = fig.add_subplot(bot[0, 1])
    else:
        fig = plt.figure(figsize=(16, 8.6))
        fig.text(0.012, 0.985, 'LAMOST Hα offset: a plausible geocoronal-subtraction residual', color='#cfd6dd',
                 fontsize=13, fontweight='bold', va='top')
        fig.text(0.012, 0.952, 'Wu et al. 2026 Tables 2, 4 · WHAM-SS DR1 · MRS-N sky model (Zhang et al. 2021) · '
                               '2024-02-01 plate · forward model, not a fit to spectra', color=FGDIM, fontsize=9.5,
                 va='top')
        gs = GridSpec(2, 3, figure=fig, width_ratios=[1.12, 1.0, 1.0], left=0.045, right=0.985, top=0.885,
                      bottom=0.085, wspace=0.27, hspace=0.34)
        axa = fig.add_subplot(gs[:, 0])
        axb = fig.add_subplot(gs[0, 1])
        axc = fig.add_subplot(gs[1, 1], sharex=axb)
        axd = fig.add_subplot(gs[:, 2])
    sect = []                                         # paper: jobs for blank_sectors()

    # ------------------------------------------------------------------ (a) profiles, P1, joint fit
    a = D['fig4a']
    v = np.array(a['v'])
    norm = np.max(a['wim'])
    wim, geo, tm = (np.array(a[k]) / norm for k in ('wim', 'geo', 'tmpl'))
    res, fitl = np.array(a['resid']) / norm, np.array(a['fit']) / norm
    axa.fill_between(v, 0, wim, color=HA_C, alpha=pk(0.18, 0.12), lw=0)
    axa.plot(v, wim, color=HA_C, lw=pk(1.6, 1.1), label='WIM Hα (true)')
    axa.plot(v, geo, color=SUN, lw=pk(1.6, 1.1), label='geocoronal Hα (true)')
    axa.plot(v, tm, color=ACCENT, lw=pk(1.6, 1.1), ls='--', label='sky template subtracted')
    axa.plot(v, res, color=FG, lw=pk(2.2, 1.5), label='residual (what is fitted)')
    axa.plot(v, fitl, color=OK, lw=pk(1.4, 1.1), ls=(0, (4, 2)), label='Gaussian + quadratic fit')
    axa.axhline(0, color=FGDIM, lw=0.6)
    # centroid ticks below the zero line, labels next to each curve's peak
    for x, c in ((a['geo_centre'], SUN), (a['tmpl_centre'], ACCENT), (a['wim_centre'], HA_C), (a['fit_centre'], OK)):
        axa.plot([x, x], [-0.15, -0.03], color=c, lw=pk(2.0, 1.3), solid_capstyle='butt')
    shift = a['fit_centre'] - a['wim_centre']
    # paper: the three peak labels sit clear above the y = 1 grid line (no grid line under or through them)
    label(axa, -21, pk(1.04, 1.15), f"geocoronal\n{sgn(a['geo_centre'])}", color=SUN, size=FS_A, ha='right',
          va='center', **BOX)
    label(axa, -1, pk(1.22, 1.24), f"template\n{sgn(a['tmpl_centre'])}", color=ACCENT, size=FS_A, ha='left',
          va='center', **BOX)
    label(axa, pk(22, 26), pk(1.08, 1.15), f"true WIM\n{a['wim_centre']:+.1f}", color=HA_C, size=FS_A, ha='left',
          va='center', **BOX)
    fit_txt = (f"fitted {a['fit_centre']:+.1f}\nshift {shift:+.1f} km/s" if PAPER
               else f"fitted {a['fit_centre']:+.1f}\n(shift {shift:+.1f} km/s)")
    label(axa, pk(36, 34), 0.80, fit_txt, color=OK, size=FS_A, ha='left', va='center', **BOX)
    over = 100 * (1 / a['kt'] - 1)
    label(axa, -72, pk(-0.30, -0.375), f"P1 · subtracted sky = {a['g_selfcons']:.2f}× true WIM Hα\n"
                                       f"template {over:.0f}% too bright, at RV(OH 6554) {a['dm']:+.1f} km/s\n"
                                       + ("ticks: line centroids" if PAPER else
                                          f"ticks: centroids (V$_{{helio}}$, km/s)"),
          color=FGDIM, size=pk(8.2, 6.7), ha='left', va='center')
    axa.set_xlim(-75, 80)
    axa.set_ylim(pk(-0.48, -0.56), pk(1.92, 2.0))
    if PAPER:
        axa.set_yticks(np.arange(0, 2.01, 0.5))       # no -0.5 grid line through the note
        # ... and no vertical grid lines through it either: blank strip below the centroid ticks
        axa.axhspan(-0.56, -0.20, facecolor=ss.AXBG, edgecolor='none', zorder=1.7)
    axa.set_xlabel('V$_{helio}$ (km/s)')
    axa.set_ylabel('flux (WIM peak = 1, LSF-convolved)')
    axa.grid(**GRID_KW)
    axa.legend(loc='upper right', fontsize=pk(8, FS_L), framealpha=pk(0.9, 1.0), handlelength=pk(2.0, 2.2))
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
        ax.axvspan(-30, 0, color=SHADE_C, alpha=SHADE_A, lw=0)
        ax.axvspan(lo, hi, color=OK, alpha=0.16, lw=0)
        ax.axvline(F['oversub_pct'], color=OK, lw=pk(1.1, 0.9))
        ax.axvline(B['oversub_pct'], color=FGDIM, lw=pk(1.0, 0.8), ls=':')
        ax.axvline(0, color=FGDIM, lw=0.6)
        ax.grid(**GRID_KW)
    mobs, eobs = off['Ha_minus_NII']['wmean'], off['Ha_minus_NII']['err']
    YT_B = 15                                         # paper: top of the y axis of panel (b)
    vals = [off['per_region'][r]['Ha_minus_NII'] for r in REGS]
    axb.axhspan(min(vals), max(vals), color=HA_C, alpha=pk(0.10, 0.08), lw=0)
    axb.axhspan(mobs - eobs, mobs + eobs, color=HA_C, alpha=pk(0.22, 0.18), lw=0)
    axb.axhline(mobs, color=HA_C, lw=pk(1.0, 0.8))
    for r in REGS:
        axb.plot(f, nan(b['shift_F'][r]), color=RCOL[r], lw=pk(1.6, 1.1), label=r)
    axb.axhline(0, color=FGDIM, lw=0.6)
    if PAPER:
        # Lines of at most 11 characters fit between the frame and the x = 0 line (left) or right of the dotted
        # line, so the two side labels cross neither. The 15-character band label fits in neither slot (it would
        # run over the x = 0 or the dotted line), so the band is named in the legend instead.
        t = label(axb, -28.5, YT_B - 0.75, 'under-\nsubtracted\n→ blueshift', color=FGDIM, size=FS_A, ha='left',
                  va='top', **BOX_SH)
        t.set_path_effects([])                        # no white halo on the grey strip
        sect.append((axb, t, -30, 0, _SH, YT_B))
        t = label(axb, 78.5, YT_B - 0.75, 'over-\nsubtracted\n→ redshift', color=FGDIM, size=FS_A, ha='right',
                  va='top', **BOX)
        sect.append((axb, t, B['oversub_pct'], 80, ss.AXBG, YT_B))
    else:
        label(axb, -28, max(vals) + 0.3, f'observed Hα − [N II]\n{mobs:+.1f} ± {eobs:.1f} km/s',
              color=HA_C, size=8.2, ha='left', va='bottom')
        label(axb, -28, 1.6, 'under-\nsubtracted\n→ blueshift', color=FGDIM, size=8, ha='left', va='center')
        label(axb, 78, 11.6, 'over-subtracted\n→ redshift', color=FGDIM, size=8, ha='right', va='top')
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch as _Patch
    hb = [(_Patch(facecolor=OK, alpha=0.3, lw=0), Line2D([], [], color=OK, lw=pk(1.1, 0.9))),
          Line2D([], [], color=FGDIM, lw=pk(1.0, 0.8), ls=':')]
    if PAPER:
        hb.insert(0, (_Patch(facecolor=HA_C, alpha=0.3, lw=0), Line2D([], [], color=HA_C, lw=0.8)))
        lb = [f'W26 Hα − [N II]\n{mobs:+.1f} ± {eobs:.1f} km/s',
              f"joint fit {F['oversub_pct']:.0f}% ({lo:.0f}–{hi:.0f}%)\n(centroid + width)",
              f"centroid-only fit {B['oversub_pct']:.0f}%"]
    else:
        lb = [f"joint fit {F['oversub_pct']:.0f}% ({lo:.0f}–{hi:.0f}%), centroid + width",
              f"centroid-only fit {B['oversub_pct']:.0f}%"]
    axb.legend(hb, lb, loc='lower right', fontsize=FS_L, framealpha=pk(0.92, 1.0), handlelength=1.8)
    axb.set_ylim(-14, pk(12, YT_B))
    axb.set_ylabel('Hα centroid shift (km/s)')
    if PAPER:
        axb.set_xlim(-30, 80)
        axb.set_xlabel('geocoronal over-subtraction, model/true − 1 (%)')
    else:
        plt.setp(axb.get_xticklabels(), visible=False)
    panel_title(axb, '(b) Centroid shift vs over-subtraction')

    # ------------------------------------------------------------------ (c) fitted width
    t2 = D['fits']['table2_Ha_fwhm']
    wv = np.array([t2[r][0] for r in REGS]); we = np.array([t2[r][1] for r in REGS])
    wm = np.sum(wv / we ** 2) / np.sum(1 / we ** 2); wme = np.sqrt(1 / np.sum(1 / we ** 2))
    whw = [D['wham_per_region'][r]['fwhm_int_nearest'] for r in REGS]
    axc.axhspan(wm - wme, wm + wme, color=HA_C, alpha=pk(0.22, 0.18), lw=0)
    axc.axhline(wm, color=HA_C, lw=pk(1.0, 0.8))
    axc.axhspan(min(whw), max(whw), facecolor='none', edgecolor=OK, hatch='////', lw=0, alpha=pk(0.5, 0.6))
    for r in REGS:
        axc.plot(f, nan(b['fwhm_F'][r]), color=RCOL[r], lw=pk(1.6, 1.1), label=r)
    if PAPER:
        # right of x = 40 all curves lie below the band, so the label sits there clear of them (and, in three
        # lines, right of the dotted line)
        t = label(axc, 78.5, wm + wme + 1.2, f'W26 Hα\n(Table 2)\n{wm:.1f} ± {wme:.1f}', color=HA_C, size=FS_A,
                  ha='right', va='bottom', **BOX)
        sect.append((axc, t, B['oversub_pct'], 80, ss.AXBG, None))
        # clear below the hatched band and left of the x = 0 line; the box (in the colour of the grey strip it
        # lies in) keeps the y = 20 grid line out of the gaps of the text
        t = label(axc, -28.5, min(whw) - 1.8, f'WHAM {min(whw):.0f}–{max(whw):.0f}\n(true width,\nmodel input)',
                  color=OK, size=FS_A, ha='left', va='top', **BOX_SH)
        t.set_path_effects([])                        # no white halo on the grey strip
        sect.append((axc, t, -30, 0, _SH, None))
    else:
        label(axc, -24, wm + wme + 0.6, f'Table 2 Hα {wm:.1f} ± {wme:.1f}', color=HA_C, size=8.3, ha='left',
              va='bottom')
        label(axc, -28, min(whw) - 0.8, f'WHAM WIM {min(whw):.0f}–{max(whw):.0f}\n(true width, model input)',
              color=OK, size=8.3, ha='left', va='top')
    axc.set_ylim(0, 68)
    axc.set_xlim(-30, 80)
    axc.legend(loc='upper right', ncol=2, fontsize=FS_L, framealpha=pk(0.9, 1.0), handlelength=1.4,
               columnspacing=0.9, title='region', title_fontsize=FS_L)
    axc.set_xlabel('geocoronal over-subtraction, model/true − 1 (%)')
    axc.set_ylabel('fitted intrinsic Hα FWHM (km/s)')
    panel_title(axc, '(c) Fitted Hα width vs over-subtraction')

    # ------------------------------------------------------------------ (d) per-region velocities
    T4 = {'P1': (16.3, 1.2, 11.8, 1.6, 11.8, 2.5), 'P2': (17.0, 1.2, 12.0, 1.6, 12.0, 1.6),
          'P3': (18.4, 1.4, 11.7, 1.9, 11.7, 1.9), 'P4': (19.7, 1.4, 13.3, 2.3, 12.7, 2.1)}   # Table 4 V_helio
    o3 = VS['combined']['V_helio']; o3e = VS['combined']['err']
    axd.axvspan(o3 - o3e, o3 + o3e, color=OIII_C, alpha=pk(0.13, 0.12), lw=0)
    axd.axvline(o3, color=OIII_C, lw=0.8, alpha=0.8)
    yl = []
    for i, r in enumerate(REGS):
        y0 = -i
        ha, hae, n2, n2e, s2, s2e = T4[r]
        wh = D['wham_per_region'][r]
        corr = ha - F['pred_shift'][r]
        rows = [(ha, hae, 'o', HA_C, True, f'{W26 or "LAMOST "}Hα (Table 4)'),
                (corr, hae, 'o', HA_C, False, 'W26 Hα − model shift' if PAPER else 'Hα − model shift (joint fit)'),
                (n2, n2e, 'D', NII_C, True, f'{W26}[N II] (Table 4)'),
                (s2, s2e, 'X', SII_C, True, f'{W26}[S II] (Table 4)'),
                (wh['v_helio_nearest'], max(wh['v_lsr_within1deg_std'], 0.5), 's', OK, True, 'WHAM Hα (nearest beam)')]
        for k, (x, e, mk, c, filled, lab) in enumerate(rows):
            yy = y0 + pk(0.30, 0.22) - pk(0.15, 0.16) * k
            axd.errorbar(x, yy, xerr=e, fmt=mk, ms=pk(6.5, 3.6), color=c, mfc=c if filled else 'none', mec=c,
                         ecolor=c, elinewidth=pk(1.1, 0.8), capsize=pk(2.5, 1.6), label=lab if i == 0 else None,
                         zorder=3, **({'mew': 0.8, 'capthick': 0.8} if PAPER else {}))
        # paper: the box (drawn below the error bars, above the grid) keeps the x = 20 grid line out of the label
        label(axd, 24.7, y0 + pk(0.43, 0.365), f"Hα − [N II] {ha - n2:+.1f}", color=HA_C, size=pk(8.2, FS_A), ha='right',
              va='center', **(dict(zorder=2.5, **BOX) if PAPER else {}))
        yl.append(y0 + pk(0.0, -0.10))
    axd.set_yticks(yl)
    axd.set_yticklabels([r + ('\n(no nebula)' if r == 'P4' else '') for r in REGS])
    for i in range(1, 4):
        axd.axhline(-i + 0.5, color=SEP_C, lw=pk(0.8, 0.6), **({'zorder': 3.2} if PAPER else {}))
    axd.set_xlim(5, 25)
    axd.set_ylim(pk(-4.3, -4.98), pk(0.62, 0.55))
    axd.set_xlabel('V$_{helio}$ (km/s)')
    axd.grid(axis='x', **GRID_KW)
    sec = axd.secondary_xaxis('top', functions=(lambda x: x - 12.2, lambda x: x + 12.2))
    sec.set_xlabel('V$_{LSR}$ (km/s, ≈ V$_{helio}$ − 12.2)', color=FGDIM, fontsize=pk(8.5, 7.0),
                   labelpad=pk(4.0, 3.2))
    sec.tick_params(colors=FGDIM, labelsize=pk(8, 7))
    if PAPER:
        axd.set_xticks([5, 10, 15, 20, 25])
        sec.set_xticks([-5, 0, 5, 10])
    from matplotlib.patches import Patch
    h, l = axd.get_legend_handles_labels()
    h.append(Patch(facecolor=OIII_C, alpha=0.3, edgecolor=OIII_C))
    l.append(f'[O III] {o3:+.1f} ± {o3e:.1f}' + ('' if PAPER else ' (Sect. 3)'))
    axd.legend(h, l, loc='lower right', fontsize=FS_L, framealpha=pk(0.92, 1.0), bbox_to_anchor=(1.0, 0.0), ncol=1)
    panel_title(axd, '(d) Velocities per region')

    if PAPER:
        blank_sectors(fig, sect)
    finish(fig, str(ROOT / 'figures/dispute_halpha_fig4.png'), dpi=200, tight=False)


if __name__ == '__main__':
    main()
