"""dispute_brightness_fig.py -- Figure 5 of the dispute document (Table 3 [O III] surface brightness).

Panels: (a) mag arcsec^-2 -> rayleigh for 30-70 A bands and AB/Vega +-0.1 mag, with Table 3;
        (b) brightness ladder of all measurements with their angular scale;
        (c) additive background error (fraction of the sky) needed to turn I_true into the Table 3 values;
        (d) multiplicative systematics vs the factor needed.
Inputs: data/derived/dispute_brightness_calc.json (run dispute_brightness_calc.py first).
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_brightness_fig.py
Paper version (light theme, 7.1 in wide, PNG + PDF in paper/figures):
     SILK_THEME=light SILK_PAPER=1 ~/.venvs/astro313/bin/python -I scripts/dispute_brightness_fig.py
"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
from silk_style import (new_sky_figure, label, finish, FG, FGDIM, ACCENT, WARN, DANGER, OK, SUN, BG2, HALO,  # noqa: E402
                        LIGHT, PAPER, EDGE)
from verify_paper_energetics import mag_to_R  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib import transforms as mtransforms  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

C = json.load(open(ROOT / 'data/derived/dispute_brightness_calc.json'))
inp, conv, rat, mult, add = C['inputs'], C['conversion'], C['ratios'], C['multiplicative'], C['additive']
P = ('P1', 'P2', 'P3')
MK = {'P1': 'o', 'P2': 's', 'P3': '^'}
GRAY = '#6b7580' if LIGHT else '#7f8a95'
W26 = 'W26' if PAPER else 'Wu+26'
if PAPER:
    # printed at 7.1 in (two-column figure*): sizes below are final print sizes
    plt.rcParams.update({'font.size': 7, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5, 'xtick.major.size': 2.5,
                         'ytick.major.size': 2.5, 'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
                         'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2.5, 'hatch.linewidth': 0.5})
else:
    plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10.5, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5})

_FS = {13: 9, 10.5: 7.5, 10: 7, 9.5: 7, 9: 6.8, 8.8: 6.6, 8.5: 6.6, 8: 6.6}


def fs(s):
    """Font size: the original value, or its print size in paper mode."""
    return _FS[s] if PAPER else s


def LW(x):
    return max(0.55, 0.65 * x) if PAPER else x


def MS(x):
    return 0.6 * x if PAPER else x


def MEW(x):
    return max(0.5, 0.75 * x) if PAPER else x


# ------------------------------------------------------------------------------------------------ layout
if PAPER:
    FW, FH = 7.1, 7.15                         # inches = print size
    A_X, A_X1, B_LX, B_X, X1 = 0.50, 2.72, 2.88, 4.47, 7.02  # (a) axes; (b) label block left, axes left; right edge
    TOP1, BOT1 = FH - 0.25, FH - 0.25 - 2.62   # top row axes
    C_X, C_X1, D_LX, D_X = 0.56, 3.06, 3.22, 4.95
    TOP2, BOT2 = 3.30, 0.40                    # bottom row axes
    fig, ax0 = new_sky_figure(FW, FH, '', '')
    ax0.remove()

    def _ax(x0, x1, y0, y1):
        return fig.add_axes([x0 / FW, y0 / FH, (x1 - x0) / FW, (y1 - y0) / FH])
    axa = _ax(A_X, A_X1, BOT1, TOP1)
    axb = _ax(B_X, X1, BOT1, TOP1)
    axc = _ax(C_X, C_X1, BOT2, TOP2)
    axd = _ax(D_X, X1, BOT2, TOP2)
    TAGX = {'a': 0.04, 'b': B_LX, 'c': 0.04, 'd': D_LX}
    TAGY = {'a': TOP1, 'b': TOP1, 'c': TOP2, 'd': TOP2}
else:
    fig, ax0 = new_sky_figure(16, 12.6, 'Table 3 [O III] brightness check',
                              'Wu+26 Table 3 (SWONS) vs DESI DR1, SDSS/BOSS fibres, SDSS imaging sky, NSNS DR0.2 · '
                              'λ5006.84 Å, AB')
    ax0.remove()
    gs_top = fig.add_gridspec(1, 2, left=0.06, right=0.835, top=0.895, bottom=0.555, width_ratios=[1.0, 1.35], wspace=0.62)
    gs_bot = fig.add_gridspec(1, 2, left=0.06, right=0.975, top=0.445, bottom=0.075, width_ratios=[1.0, 1.05], wspace=0.42)
    axa = fig.add_subplot(gs_top[0])
    axb = fig.add_subplot(gs_top[1])
    axc = fig.add_subplot(gs_bot[0])
    axd = fig.add_subplot(gs_bot[1])


def ptag(ax, s, txt):
    if PAPER:   # letter at the upper left of the panel block (labels included), short title after it
        y = (TAGY[s] + 0.05) / FH
        fig.text(TAGX[s] / FW, y, f'({s})', fontsize=fs(13), fontweight='bold', color=FG, va='bottom', ha='left')
        fig.text((TAGX[s] + 0.30) / FW, y + 0.004 / FH, txt, fontsize=fs(10.5), color=FGDIM, va='bottom', ha='left')
        return
    ax.text(0.0, 1.015, s, transform=ax.transAxes, fontsize=13, fontweight='bold', color=FG, va='bottom', ha='left')
    ax.text(0.045, 1.015, txt, transform=ax.transAxes, fontsize=10.5, color=FGDIM, va='bottom', ha='left')


# ------------------------------------------------------------------------------------------------ (a) conversion
m = np.linspace(23.12, 26.75, 400)
if LIGHT:   # sequential blues, darker with width; the adopted 60 A line in the text colour
    cols = {30: '#6aa6d1', 40: '#3f86bd', 50: '#1f659f', 60: FG, 70: '#0d3d6b'}
    SHADE_C, SHADE_A = FG, 0.13
else:
    cols = {30: '#3b6e8f', 40: '#4f8fb3', 50: '#6fb3d9', 60: '#ffffff', 70: '#9fd6f2'}
    SHADE_C, SHADE_A = '#ffffff', 0.16
axa.axhspan(1.4, 2.5, color=OK, alpha=0.16, lw=0)
axa.fill_between(m, mag_to_R(m, 60, 0.1), mag_to_R(m, 60, -0.1), color=SHADE_C, alpha=SHADE_A, lw=0)
for W in (30, 40, 50, 70, 60):
    axa.plot(m, mag_to_R(m, W), color=cols[W], lw=LW(2.2) if W == 60 else LW(1.2), zorder=3)
for W, dy in ((30, 1.0), (40, 1.0), (50, 1.0), (60, 1.0), (70, 1.0)):
    yv = mag_to_R(23.12, W)
    label(axa, 23.08, yv * dy, f'{W} Å', color=cols[W], size=fs(9) if not PAPER else 6.6, va='center', ha='left')
for k in P:
    axa.plot(inp['table3_mag'][k], conv['AB_W60'][k], MK[k], color=DANGER, ms=MS({'P1': 8, 'P2': 10.5, 'P3': 7}[k]),
             mec='white', mew=MEW(0.7), zorder=5)
label(axa, 24.05, 13.5, 'Table 3 P1–P3', color=DANGER, size=fs(9.5), ha='right', va='center')
q = inp['qhyccd_mag']
axa.plot([q, q], conv['qhyccd_R_W30_W60'], color=WARN, lw=LW(2.2), zorder=4)
axa.plot(q, mag_to_R(q, 50), 'D', mfc='none', mec=WARN, ms=MS(7), mew=MEW(1.6), zorder=5)
label(axa, q + 0.08, 8.2, 'QHYCCD 24.3', color=WARN, size=fs(9), ha='right', va='center')
for I in (1.4, 2.5):
    mm = conv['mag_for_R_W60'][str(I)]
    axa.plot([mm, mm], [0.8, I], color=OK, lw=LW(0.9), ls=':', zorder=2)
YDM = 2.1 if PAPER else 2.0   # height of the delta-m arrow (paper: lifted above the two-line band label)
if PAPER:
    label(axa, 22.68, 1.44, 'DESI · SDSS · NSNS\n1.4–2.5 R', color=OK, size=fs(8.5), va='bottom', ha='right',
          multialignment='right', linespacing=1.15)
else:
    label(axa, 24.75, 1.47, 'DESI · SDSS · NSNS 1.4–2.5 R', color=OK, size=9, va='bottom', ha='left')
axa.annotate('', xy=(23.9, YDM), xytext=(25.6, YDM),
             arrowprops=dict(arrowstyle='<->', color=FG, lw=LW(1.0), shrinkA=0, shrinkB=0) if PAPER else
             dict(arrowstyle='<->', color=FG, lw=1.0))
label(axa, 24.55, YDM + 0.18 if not PAPER else YDM * 1.06,
      f'Δm {rat["delta_mag_region_matched"][0]:.1f}–{rat["delta_mag_region_matched"][1]:.1f} mag',
      color=FG, size=fs(9), ha='center', va='bottom')
if PAPER:
    label(axa, 22.68, 0.86, 'shade: ±0.1 mag\nAB/Vega at 60 Å', color=FGDIM, size=fs(8.5), ha='right', va='bottom',
          multialignment='right', linespacing=1.15)
else:
    label(axa, 25.05, 0.88, 'shade: ±0.1 mag AB/Vega at 60 Å', color=FGDIM, size=8.5, ha='left', va='bottom')
axa.set_xlim(26.75, 22.6 if PAPER else 22.75)
axa.set_yscale('log')
axa.set_ylim(0.8, 30)
axa.set_yticks([1, 2, 3, 5, 10, 20])
axa.set_yticklabels(['1', '2', '3', '5', '10', '20'])
if PAPER:
    axa.set_xticks([26, 25, 24, 23])
    axa.set_xticks(np.arange(23.0, 26.75, 0.5), minor=True)
    axa.set_xticklabels(['26', '25', '24', '23'])
    axa.set_xlabel('[O III] band (mag arcsec$^{-2}$)')
else:
    axa.set_xlabel('[O III]-band surface brightness (mag arcsec$^{-2}$)')
axa.set_ylabel('I(5007) (R)')
axa.grid(alpha=0.35, lw=0.5)
ptag(axa, 'a', 'mag → rayleigh, line-only band')

# ------------------------------------------------------------------------------------------------ (b) ladder
rows = [f'{W26} Table 3', 'Table 3, lowering extremes', 'QHYCCD note (SWONS)', 'MRS-N rel. × sky line',
        'NSNS DR0.2', 'SDSS imaging sky (g)', 'SDSS/BOSS fibres', 'SDSS/BOSS south', 'DESI DR1']
scales = ['subregions, 50″ px', 'panel d stack', 'brightest region', 'sky ratio ±×5',
          '≤0.25° of ridge, 51″ px', '17′ ridge, 19 runs', '2–3″ fibres', 'ON−OFF, 0.15°',
          '1.5″ fibres']
ny = len(rows)
Y = {r_: ny - 1 - i for i, r_ in enumerate(rows)}
if LIGHT:
    axb.axvspan(*inp['fesen_limit_R_as_quoted'], color=FG, alpha=0.03, lw=0)
    axb.axvspan(*inp['fesen_limit_R_as_quoted'], fill=False, hatch='//', ec='#c3cad1', lw=0, zorder=0.5)
else:
    axb.axvspan(*inp['fesen_limit_R_as_quoted'], color='#ffffff', alpha=0.05, lw=0)
    axb.axvspan(*inp['fesen_limit_R_as_quoted'], fill=False, hatch='//', ec='#2e363f', lw=0, zorder=0.5)
axb.axvspan(1.4, 2.5, color=OK, alpha=0.16, lw=0)
label(axb, 0.52, ny - 0.35, f'Fesen+24 limits 0.5–5 R ({W26} quote)', color=FGDIM, size=fs(8.5), va='center', ha='left')
label(axb, 1.43, -0.65, 'adopted 1.4–2.5 R', color=OK, size=fs(8.5), va='center', ha='left')


def pts(row, vals, filled=True, color=ACCENT, errs=None, ms=8):
    ms = MS(ms)
    msz = {'P1': ms, 'P2': ms * 1.3, 'P3': ms * 0.9}
    y = Y[row]
    for k in P:
        if vals.get(k) is None:
            continue
        kw = dict(mfc=color if filled else 'none', mec='white' if filled else color, mew=MEW(0.7) if filled else MEW(1.5))
        e = None if errs is None else errs[k]
        axb.errorbar(vals[k], y, xerr=e, fmt=MK[k], color=color, ms=msz[k], capsize=0, elinewidth=LW(1.2), zorder=5, **kw)


def bar(row, lo, hi, color, ls='-', lw=5, alpha=0.55):
    axb.plot([lo, hi], [Y[row], Y[row]], color=color, lw=LW(lw), alpha=alpha, solid_capstyle='butt', ls=ls, zorder=3)


bar(rows[0], *conv['all_variants_range_R'], DANGER, lw=3, alpha=0.45)
pts(rows[0], conv['AB_W60'], color=DANGER)
pts(rows[1], mult['stacked_R'], filled=False, color=DANGER)
bar(rows[2], *conv['qhyccd_R_W30_W60'], WARN)
axb.plot(mag_to_R(inp['qhyccd_mag'], 50), Y[rows[2]], 'D', mfc='none', mec=WARN, ms=MS(7), mew=MEW(1.6), zorder=5)
lm = inp['lamost_mrsn_scaled_R']
bar(rows[3], *lm['APO_16_84'], GRAY, ls=(0, (2, 2)), lw=1.4, alpha=0.9)
pts(rows[3], dict(zip(P, lm['APO_ratio'])), color=GRAY, ms=7)
pts(rows[3], dict(zip(P, lm['UVES_ratio'])), filled=False, color=GRAY, ms=7)
label(axb, 11.5, Y[rows[3]] + 0.33, 'UVES ratio', color=GRAY, size=fs(8), ha='center', va='bottom')
label(axb, 2.9, Y[rows[3]] + 0.33, 'APO ratio', color=GRAY, size=fs(8), ha='center', va='bottom')
pts(rows[4], inp['nsns_onridge_excess_R'])
pts(rows[4], inp['nsns_peak_excess_R'], filled=False)
im = inp['sdss_imaging_R']
axb.errorbar(im['g_PSP'][0], Y[rows[5]], xerr=im['g_PSP'][1], fmt='D', color=ACCENT, mec='white', mew=MEW(0.7), ms=MS(7),
             elinewidth=LW(1.2), zorder=5)
for kk in ('g_FramesSub', 'g_Frames'):
    axb.plot(im[kk][0], Y[rows[5]], 'D', mfc='none', mec=ACCENT, mew=MEW(1.4), ms=MS(7), zorder=5)
pts(rows[6], inp['sdss_P_mean_R'], errs=inp['sdss_P_err_R'])
axb.plot(inp['sdss_brightest_fibre_R'], Y[rows[6]], '*', mfc='none', mec=ACCENT, mew=MEW(1.3), ms=MS(11), zorder=5)
eon = float(np.hypot(inp['sdss_south_onoff_err_stat_R'], inp['sdss_south_onoff_err_placebo_R']))
axb.errorbar(inp['sdss_south_onoff_R'], Y[rows[7]], xerr=eon, fmt='D', color=ACCENT, mec='white', mew=MEW(0.7), ms=MS(7),
             elinewidth=LW(1.2), zorder=5)
nodes = [v[0] for k, v in inp['desi_nodes_R'].items() if k.startswith('E@1') and 14 <= float(k[2:]) <= 18]
bar(rows[8], min(nodes), max(nodes), ACCENT, lw=4, alpha=0.4)
d0 = inp['desi_south_median_R']
axb.errorbar(d0, Y[rows[8]], xerr=d0 * inp['desi_sys_frac'], fmt='D', color=ACCENT, mec='white', mew=MEW(0.7), ms=MS(7),
             elinewidth=LW(1.2), zorder=5)
axb.plot(inp['desi_nodes_R']['W@18.00'][0] if 'W@18.00' in inp['desi_nodes_R'] else 3.88, Y[rows[8]], 's', mfc='none',
         mec=ACCENT, mew=MEW(1.4), ms=MS(7), zorder=5)
axb.set_xscale('log')
axb.set_xlim(0.5, 20)
axb.set_ylim(-1.05, ny - 0.1)
axb.set_xticks([0.5, 1, 2, 3, 5, 10, 20])
axb.set_xticklabels(['0.5', '1', '2', '3', '5', '10', '20'])
axb.set_yticks([Y[r_] for r_ in rows])


def _rowcol(r_):
    return (DANGER if r_.startswith('Wu') or r_.startswith('W26') or r_.startswith('Table')
            else (WARN if 'QHY' in r_ else (GRAY if 'MRS' in r_ else FG)))


if PAPER:   # name over angular scale, both left of the axis (no right-hand twin axis at print width)
    axb.set_yticklabels([''] * ny)
    axb.tick_params(axis='y', length=0)
    tr = mtransforms.blended_transform_factory(axb.transAxes, axb.transData) + \
        mtransforms.ScaledTranslation(-4 / 72, 0, fig.dpi_scale_trans)
    for r_, sc in zip(rows, scales):
        axb.text(0, Y[r_] + 0.2, r_, transform=tr, fontsize=fs(9.5), color=_rowcol(r_), ha='right', va='center')
        axb.text(0, Y[r_] - 0.24, sc, transform=tr, fontsize=fs(8.5), color=FGDIM, ha='right', va='center')
else:
    axb.set_yticklabels(rows, fontsize=9.5)
    for t, r_ in zip(axb.get_yticklabels(), rows):
        t.set_color(_rowcol(r_))
    axr = axb.twinx()
    axr.set_ylim(axb.get_ylim())
    axr.set_yticks([Y[r_] for r_ in rows])
    axr.set_yticklabels(scales, fontsize=8.5, color=FGDIM)
    axr.tick_params(axis='y', length=0)
axb.set_xlabel('I(5007) (R)')
axb.grid(axis='x', alpha=0.35, lw=0.5)
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
h = [Line2D([], [], ls='', marker=MK[k], color=FG, mec=FG, ms=MS(7), label=k) for k in P]
h += [Line2D([], [], ls='', marker='D', mfc=FG, mec=FG, ms=MS(6), label='ridge-wide')]
filled = tuple(Line2D([], [], ls='', marker=mk, mfc=FG, mec=FG, ms=MS(6)) for mk in ('o', 's', '^', 'D'))
opened = tuple(Line2D([], [], ls='', marker=mk, mfc='none', mec=FG, mew=MEW(1.2), ms=MS(6)) for mk in ('o', 's', '^', 'D'))
star = Line2D([], [], ls='', marker='*', mfc='none', mec=FG, mew=MEW(1.1), ms=MS(9))
if PAPER:
    labs = ['P1', 'P2', 'P3', 'ridge-wide', 'filled: mean/median', 'open: peak/alt.', 'brightest fibre']
    axb.legend(handles=h + [filled, opened, star], labels=labs, ncol=4, loc='upper center',
               bbox_to_anchor=((B_LX + X1) / 2 / FW, (BOT1 - 0.36) / FH), bbox_transform=fig.transFigure,
               frameon=False, fontsize=fs(8.8), handletextpad=0.5, columnspacing=1.0, handlelength=3.0,
               borderpad=0.2, labelspacing=0.35, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.1)})
else:
    labs = ['P1', 'P2', 'P3', 'ridge-wide', 'filled: mean / median', 'open: peak / alt.', 'brightest fibre']
    axb.legend(handles=h + [filled, opened, star], labels=labs, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.115),
               frameon=False, fontsize=8.8, handletextpad=0.7, columnspacing=1.6, handlelength=3.4,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)})
ptag(axb, 'b', 'all measurements on one scale')

# ------------------------------------------------------------------------------------------------ (c) additive map
It = np.linspace(0.5, 4.0, 300)
sky = np.linspace(20.8, 22.6, 300)
X, S = np.meshgrid(It, sky)
I3m = float(np.mean([conv['AB_W60'][k] for k in P]))
Z = 100 * (I3m - X) / mag_to_R(S, 60)
lev = [4, 6, 8, 10, 12.5, 15, 20, 25]
if LIGHT:   # light sequential map (pale = small error), dark contour lines and labels
    cmap_c = ListedColormap(plt.get_cmap('PuBu')(np.linspace(0.03, 0.62, 256)))
    cf = axc.contourf(X, S, Z, levels=[3] + lev + [32], cmap=cmap_c, extend='max')   # fill the >32% corner
    cs = axc.contour(X, S, Z, levels=lev, colors='#3c4650', linewidths=LW(0.8))
    CLAB_C = FG
else:
    cf = axc.contourf(X, S, Z, levels=[3] + lev + [32], cmap='cividis', alpha=0.85)
    cs = axc.contour(X, S, Z, levels=lev, colors='#0b0e12', linewidths=0.8)
    CLAB_C = '#ffffff'
man = []
for L in lev:
    for xx in (3.6, 0.8):
        sk = 21.5 - 2.5 * np.log10((I3m - xx) / (L / 100) / mag_to_R(21.5, 60))
        if 20.85 < sk < 22.55:
            man.append((xx, sk))
            break
lbl = axc.clabel(cs, fmt=lambda v: f'{v:g}%', fontsize=fs(9), inline=True, inline_spacing=6 if not PAPER else 3,
                 manual=man)
for t in lbl:
    t.set_color(CLAB_C)
    t.set_path_effects([HALO])
from matplotlib.patches import Rectangle  # noqa: E402
axc.add_patch(Rectangle((1.4, 21.1), 1.1, 0.9, fill=False, ec=OK, lw=LW(2.0), zorder=5))
fr = add['frac_of_sky_range_21.1_22.0']
bx = [100 * (I3m - 2.5) / mag_to_R(21.1, 60), 100 * (I3m - 1.4) / mag_to_R(22.0, 60)]   # box corners, mean of P1-P3
BOXTXT = '{:.0f}–{:.0f}% of sky\nP1–P3: {:.0f}–{:.0f}%' if PAPER else '{:.0f}–{:.0f}% of sky\n(P1–P3: {:.0f}–{:.0f}%)'
label(axc, 1.95, 21.72, BOXTXT.format(bx[0], bx[1], 100 * fr[0], 100 * fr[1]),
      color=FG if LIGHT else '#ffffff', size=fs(10), ha='center', va='center', zorder=6, fontweight='bold')
axc.axhline(inp['sqm_mag'], color=SUN, lw=LW(1.0) if not PAPER else 0.9, ls='--', zorder=4)
label(axc, 2.62, inp['sqm_mag'] + 0.035, f'SQM 21.5 ({W26})', color=SUN, size=fs(8.5), va='bottom', ha='left')
label(axc, 1.42, 21.05, 'measured I_true', color=OK, size=fs(8.5), va='top', ha='left')
if PAPER:   # the Table 3 mean moves from the panel title into the panel
    label(axc, 3.95, 20.83, f'Table 3 mean {I3m:.1f} R', color=FG, size=fs(8.5), va='bottom', ha='right')
axc.set_xlabel('true [O III] brightness I_true (R)')
axc.set_ylabel('sky in the [O III] band (mag arcsec$^{-2}$)')
axc.set_xlim(0.5, 4.0)
axc.set_ylim(20.8, 22.6)
if PAPER:
    ptag(axc, 'c', 'background error needed, % of sky')
else:
    ptag(axc, 'c', f'background error needed, % of sky (Table 3 mean {I3m:.1f} R)')

# ------------------------------------------------------------------------------------------------ (d) multiplicative budget
F = mult['factors_true_over_T3']
fil = conv['filter_true_over_paper']
NEU = EDGE if LIGHT else FG     # neutral bars: mid grey on white, light grey on black
items = [('AB vs Vega, ±0.1 mag', F['AB_Vega_pm0.1mag'], FG),
         ('effective width 50–70 Å', F['width_50_70A'], FG),
         ('[O III]4959 leakage', F['4959_leak_0_to_full'], FG),
         ('filter at stop, field 0–5°', fil['collimated_0_5deg'], FG),
         ('filter at stop, field 10–15°', fil['collimated_10_15deg'], FG),
         ('filter in f/2.8 beam', fil['f2.8_converging'], FG),
         ('untuned filter, f/0.95 beam', fil['f0.95_converging_untuned'], WARN),
         ('P1 starlight, ~0.3 mag', [F['P1_star_0.3mag'], F['P1_star_0.3mag']], FG),
         ('lowering extremes stacked', [mult['stack_lowering_P1'], mult['stack_lowering_P2P3']], DANGER)]
need = rat['required_factor_true_over_T3']
axd.axvspan(need[0], need[1], color=OK, alpha=0.22, lw=0)
axd.axvline(1.0, color=FGDIM, lw=LW(0.9))
n = len(items)
for i, (nm, (lo, hi), col) in enumerate(items):
    y = n - 1 - i
    cb = NEU if col == FG else col
    if abs(hi - lo) < 1e-3:
        axd.plot(lo, y, 'o', color=cb, ms=MS(8), mec='white', mew=MEW(0.7), zorder=5)
    else:
        axd.plot([lo, hi], [y, y], color=cb, lw=LW(7), solid_capstyle='butt', alpha=0.85, zorder=4)
axd.set_yticks([n - 1 - i for i in range(n)])
axd.set_yticklabels([it[0] for it in items], fontsize=fs(9.5))
for t, it in zip(axd.get_yticklabels(), items):
    t.set_color(it[2] if it[2] != FG else FG)
yl = n - 1 - [it[0] for it in items].index('lowering extremes stacked')
axd.annotate('', xy=(need[1] * 1.03, yl), xytext=(mult['stack_lowering_P1'] * 0.97, yl),
             arrowprops=dict(arrowstyle='->', color=FG, lw=LW(1.1)))
g = mult['stacked_over_independent_range']
label(axd, 0.235, yl - 0.55, f'still {g[0]:.1f}–{g[1]:.1f}× too bright', color=FG, size=fs(9), ha='left', va='center')
label(axd, np.sqrt(need[0] * need[1]), n - 0.15, f'needed\n{need[0]:.2f}–{need[1]:.2f}', color=OK, size=fs(9),
      ha='center', va='top')
label(axd, 2.9, n - 1 - [it[0] for it in items].index('untuned filter, f/0.95 beam') + 0.42, 'raises I', color=WARN,
      size=fs(8.5), ha='center', va='bottom')
axd.set_xscale('log')
axd.set_xlim(0.08, 6)
axd.set_ylim(-1.0, n)
axd.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
axd.set_xticklabels(['0.1', '0.2', '0.5', '1', '2', '5'])
axd.set_xlabel('factor I_true / I_Table3')
axd.grid(axis='x', alpha=0.35, lw=0.5)
ptag(axd, 'd', 'multiplicative systematics vs factor needed')

out = ROOT / 'figures/dispute_brightness_5.png'
out = finish(fig, out, dpi=200, tight=False)
print('wrote', out)
