"""F6 cross-spine profiles of the southern Silk (spine Dec < 17.6) from data/derived/w3_figs_f6_profiles.npz/.json.

Row 1: [O III] in Rayleigh (NSNS DR0.2 median profile with its null band; DESI DR1 fibre residuals, petal-relative,
baseline-subtracted, mean +- s.e.). Other rows: profile minus its mean over the wings 0.5 <= |x| <= x_max, divided by the null
scatter sigma_null (median over x of the robust std of the mean-subtracted null profiles), i.e. in units of the
placement-null sigma; bands = 16-84 % and 2.5-97.5 % of the null profiles on the same scale.
Output: figures/w3_figs_f6_profiles.png

Paper mode (SILK_THEME=light SILK_PAPER=1, see silk_style.py): light palette, journal layout at print size
7.1 x 7.4 in ([O III] full width on top, the six normalised profiles in a 3 x 2 grid below), panel letters (a)-(g),
null-band key in panel (a), no figure title or credit; PNG (300 dpi) + PDF go to paper/figures/ via ss.save_paper.
Same data, normalisation and numbers as the default run; the default (dark) run is unchanged.
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, panel_tag

PAPER = ss.PAPER
LIGHT = ss.LIGHT
NUV_C = NII_C if LIGHT else '#c792ea'      # light theme: saturated purple instead of the pastel
XR_C = '#ad1457' if LIGHT else '#ff9e80'   # light theme: deep pink instead of the pale salmon

ROWS = [('fuv', 'GALEX FUV', OIII_C, 1.2, 'CU', '%.0f'),
        ('nuv', 'GALEX NUV', NUV_C, 1.2, 'CU', '%.0f'),
        ('shassa', 'SHASSA Hα+[N II] · Dec 12–15.6', HA_C, 2.0, 'R', '%.2f'),
        ('erass', 'eRASS1 0.2–0.6 keV', XR_C, 2.5, 'cts s⁻¹ deg⁻²', '%.2f'),
        ('hi', 'GALFA HI · V_LSR −2…+6', HB_C, 2.0, '10¹⁹ cm⁻²', '%.1f'),
        ('iris', 'IRIS 100 µm', SII_C, 2.0, 'MJy sr⁻¹', '%.2f')]
NULLDESC = {'fuv': 'RA shifts', 'nuv': 'RA shifts', 'shassa': 'random placements', 'erass': 'rotations about PSR',
            'hi': 'RA/Dec shifts', 'iris': 'RA/Dec shifts'}


def norm_rows(x, p, nl, xmax):
    m = np.abs(x) <= xmax + 1e-9
    off = m & (np.abs(x) >= 0.5)                 # reference level = mean of the off-ridge wings
    pm = p - np.nanmean(p[off])
    nm = nl - np.nanmean(nl[:, off], axis=1, keepdims=True)
    med = np.nanmedian(nm, 0)
    sig = 1.4826 * np.nanmedian(np.abs(nm - med), 0)
    s = float(np.nanmedian(sig[m]))
    return m, pm / s, nm / s, s


def main():
    P = np.load(DER / 'w3_figs_f6_profiles.npz')
    J = json.loads((DER / 'w3_figs_f6_profiles.json').read_text())['bands']
    if PAPER:
        main_paper(P, J)
        return
    fig, axs = plt.subplots(len(ROWS) + 1, 1, figsize=(9.2, 13.2), sharex=True,
                            gridspec_kw=dict(hspace=0.10, height_ratios=[1.45] + [1] * len(ROWS)))
    fig.text(0.012, 0.992, 'Silk Nebula · cross-spine profiles (south, spine Dec < 17.6)', color=TITLE,
             fontsize=13, fontweight='bold', va='top')
    fig.text(0.012, 0.975, 'offset from the E spine, + = east (outside, away from Monogem); grey = null placements '
             '(16–84 %, 2.5–97.5 %)', color=FGDIM, fontsize=9, va='top')
    # ---- [O III] row in R
    ax = axs[0]
    x = P['nsns_x']; p = P['nsns']; nl = P['nsns_null']
    m = np.abs(x) <= 2.0
    base = np.nanmedian(p[m & (np.abs(x) > 0.6)])
    nb = nl - np.nanmedian(np.where((np.abs(x) > 0.6) & m, nl, np.nan), axis=1, keepdims=True)
    ax.fill_between(x[m], np.nanpercentile(nb, 2.5, 0)[m], np.nanpercentile(nb, 97.5, 0)[m], color=FGDIM, alpha=0.15, lw=0)
    ax.fill_between(x[m], np.nanpercentile(nb, 16, 0)[m], np.nanpercentile(nb, 84, 0)[m], color=FGDIM, alpha=0.25, lw=0)
    ax.plot(x[m], p[m] - base, color=OIII_C, lw=1.8, label='NSNS DR0.2 (2′ blocks, median)')
    xd = P['desi_x']; md = np.abs(xd) <= 2.0
    ax.errorbar(xd[md], P['desi'][md], yerr=P['desi_err'][md], fmt='o', ms=3.5, color=FG if LIGHT else '#ffffff', mec=BG if LIGHT else 'k', mew=0.4,
                lw=0.8, capsize=0, label='DESI DR1 fibre residuals (GFA-corrected)')
    ax.set_ylabel('[O III] (R)')
    ax.legend(loc='upper right', fontsize=7.6, framealpha=0.85)
    panel_tag(ax, '[O III] 5007 · NSNS ridge %.2f R (z = %.0f) · DESI %.2f R' %
              (J['nsns']['A'], J['nsns']['z'], J['desi']['A']), size=8.6)
    ax.set_ylim(-0.6, 3.2)
    # ---- normalised rows
    for ax, (key, name, c, xmax, unit, fmt) in zip(axs[1:], ROWS):
        x = P[f'{key}_x']; p = P[key]; nl = np.asarray(P[f'{key}_null'], float)
        m, pz, nz, s = norm_rows(x, p, nl, xmax)
        ax.fill_between(x[m], np.nanpercentile(nz, 2.5, 0)[m], np.nanpercentile(nz, 97.5, 0)[m], color=FGDIM, alpha=0.15,
                        lw=0)
        ax.fill_between(x[m], np.nanpercentile(nz, 16, 0)[m], np.nanpercentile(nz, 84, 0)[m], color=FGDIM, alpha=0.25,
                        lw=0)
        ax.plot(x[m], pz[m], color=c, lw=1.6)
        b = J[key]
        txt = f'{name} · ridge {fmt % b["A"]} ± {fmt % b["A_null_sig"]} {unit} (z = {b["z"]:+.1f}, {b["n_null"]} {NULLDESC[key]})'
        panel_tag(ax, txt, size=8.0)
        ax.text(0.995, 0.06, f'σ_null = {fmt % s} {unit}', transform=ax.transAxes, ha='right', va='bottom',
                fontsize=7.3, color=FGDIM, path_effects=[HALO])
        ax.set_ylabel('(I − ⟨I⟩_wings)/σ_null', fontsize=7.5)
        lo = np.nanmin([np.nanpercentile(nz, 2.5, 0)[m].min(), np.nanmin(pz[m])])
        hi = np.nanmax([np.nanpercentile(nz, 97.5, 0)[m].max(), np.nanmax(pz[m])])
        ax.set_ylim(lo - 0.6, hi + 0.35 * (hi - lo) + 0.8)
    for ax in axs:
        ax.axvline(0, color=FGDIM, lw=0.7, ls=':')
        ax.axhline(0, color=FGDIM, lw=0.5, alpha=0.6)
        ax.grid(alpha=0.18)
        ax.tick_params(labelsize=8)
    axs[-1].set_xlim(-2.0, 2.0)
    axs[-1].set_xlabel('offset from E spine (deg)')
    axs[-1].text(0.0, -0.33, '← west · inside, toward Monogem', transform=axs[-1].transAxes, fontsize=8, color=FGDIM,
                 ha='left', va='top')
    axs[-1].text(1.0, -0.33, 'east · outside →', transform=axs[-1].transAxes, fontsize=8, color=FGDIM, ha='right',
                 va='top')
    ss.credit(fig)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.962, bottom=0.065)
    ss.save_paper(fig, FIG / 'w3_figs_f6_profiles.png', dpi=200)
    plt.close(fig)
    print('saved')


# ------------------------------------------------------------------------------------------------ paper version
# Print size 7.1 x 7.4 in (two-column figure*; with its caption it fits one AASTeX page unscaled), drawn 1:1:
# every font >= 6.5 pt, every line >= 0.5 pt.
P_W, P_H = 7.1, 7.4
P_FS = dict(tick=7.0, label=7.5, tag=7.5, text=7.0, note=6.8, legend=7.0)
P_BAND = ('#e2e6eb', '#b4bdc7') if LIGHT else ((0.60, 0.64, 0.68, 0.15), (0.60, 0.64, 0.68, 0.25))  # 2.5-97.5, 16-84 %
P_NULLDESC = {'fuv': 'RA shifts', 'nuv': 'RA shifts', 'shassa': 'random placements', 'erass': 'rotations about PSR',
              'hi': 'RA/Dec shifts', 'iris': 'RA/Dec shifts'}
P_NAME = {'fuv': 'GALEX FUV', 'nuv': 'GALEX NUV', 'shassa': 'SHASSA Hα+[N II] · Dec 12–15.6',
          'erass': 'eRASS1 0.2–0.6 keV', 'hi': 'GALFA HI · V_LSR −2…+6 km s⁻¹', 'iris': 'IRIS 100 µm'}
P_DIR_W, P_DIR_E = '← west · inside, toward Monogem', 'east · outside →'


def _p_axes(fig, x0, y0, w, h):
    """Axes at (x0, y0, w, h) in inches from the lower left of the P_W x P_H figure."""
    return fig.add_axes([x0 / P_W, y0 / P_H, w / P_W, h / P_H])


def _p_style(ax, major):
    from matplotlib.ticker import MultipleLocator
    ax.set_xlim(-2.0, 2.0)
    ax.xaxis.set_major_locator(MultipleLocator(major))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.axvline(0, color=FGDIM, lw=0.7, ls=':', zorder=1.5)
    ax.axhline(0, color=FGDIM, lw=0.5, zorder=1.5)
    ax.grid(color=GRID, lw=0.5, alpha=1.0 if LIGHT else 0.18, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=P_FS['tick'], length=3, width=0.6, pad=2)
    ax.tick_params(which='minor', length=1.8, width=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def _p_bands(ax, x, nb):
    ax.fill_between(x, np.nanpercentile(nb, 2.5, 0), np.nanpercentile(nb, 97.5, 0), color=P_BAND[0], lw=0, zorder=1)
    ax.fill_between(x, np.nanpercentile(nb, 16, 0), np.nanpercentile(nb, 84, 0), color=P_BAND[1], lw=0, zorder=1.2)


def _p_xlabel(ax, dy, two_lines):
    """x label and the west/east reading aids below an axes (offsets in points)."""
    kw = dict(xycoords='axes fraction', textcoords='offset points', va='top', fontsize=P_FS['label'])
    ax.annotate('offset from E spine (deg)', xy=(0.5, 0), xytext=(0, dy), ha='center', color=FG, **kw)
    dy2 = dy - 10 if two_lines else dy
    kw['fontsize'] = P_FS['note']
    ax.annotate(P_DIR_W, xy=(0, 0), xytext=(0, dy2), ha='left', color=FGDIM, **kw)
    ax.annotate(P_DIR_E, xy=(1, 0), xytext=(0, dy2), ha='right', color=FGDIM, **kw)


def main_paper(P, J):
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator
    fig = plt.figure(figsize=(P_W, P_H))
    # vertical budget (in): 0.05 top margin + top_h + 0.44 (x label of (a), gap) + 3 x (tit_h + h) + 2 x gap_r + bottom
    left, right, gapc, bottom = 0.50, 0.13, 0.40, 0.48
    colw = (P_W - left - right - gapc) / 2
    top_h = 1.60
    top_y = P_H - 0.05 - top_h
    tit_h, gap_r = 0.40, 0.12                                       # three title lines above each grid panel
    grid_top = top_y - 0.44
    h = (grid_top - bottom - 3 * tit_h - 2 * gap_r) / 3             # 1.13 in

    # ---- (a) [O III] in R, full width
    ax = _p_axes(fig, left, top_y, P_W - left - right, top_h)
    x = P['nsns_x']; p = P['nsns']; nl = P['nsns_null']
    m = np.abs(x) <= 2.0
    base = np.nanmedian(p[m & (np.abs(x) > 0.6)])
    nb = nl - np.nanmedian(np.where((np.abs(x) > 0.6) & m, nl, np.nan), axis=1, keepdims=True)
    _p_bands(ax, x[m], nb[:, m])
    ax.plot(x[m], p[m] - base, color=OIII_C, lw=1.4, zorder=3, label='NSNS DR0.2 (2′ blocks, median)')
    xd = P['desi_x']; md = np.abs(xd) <= 2.0
    ax.errorbar(xd[md], P['desi'][md], yerr=P['desi_err'][md], fmt='o', ms=3.0, color=FG, mec=BG, mew=0.5,
                elinewidth=0.7, capsize=0, zorder=4, label='DESI DR1 fibre residuals (GFA-corrected)')
    h0, l0 = ax.get_legend_handles_labels()
    h0 += [Patch(color=P_BAND[1], lw=0), Patch(color=P_BAND[0], lw=0)]
    l0 += ['null placements, 16–84 %', 'null placements, 2.5–97.5 %']
    lg = ax.legend(h0, l0, loc='upper right', fontsize=P_FS['legend'], framealpha=0.92, borderaxespad=0.4,
                   handlelength=1.8, labelspacing=0.35, borderpad=0.45)
    lg.get_frame().set_linewidth(0.5)
    _p_style(ax, 0.5)
    ax.set_ylim(-0.6, 3.2)
    ax.set_ylabel('[O III] 5007 (R)', fontsize=P_FS['label'], labelpad=3)
    ypt = 1 / (72 * top_h)                                          # 1 pt in axes-fraction units
    panel_tag(ax, '(a) [O III] 5007', x=0.008, y=1 - 2.8 * ypt, size=P_FS['tag'])
    ax.text(0.008, 1 - 11.8 * ypt, 'NSNS ridge %.2f R (z = %.0f) · DESI %.2f R' % (J['nsns']['A'], J['nsns']['z'], J['desi']['A']),
            transform=ax.transAxes, ha='left', va='top', fontsize=P_FS['text'], color=FG, path_effects=[HALO])
    _p_xlabel(ax, -15, two_lines=False)

    # ---- (b)-(g) normalised profiles, 3 x 2 grid
    axs = []
    for i, (key, name, c, xmax, unit, fmt) in enumerate(ROWS):
        r, k = divmod(i, 2)
        y0 = grid_top - (r + 1) * (tit_h + h) - r * gap_r
        ax = _p_axes(fig, left + k * (colw + gapc), y0, colw, h)
        axs.append(ax)
        x = P[f'{key}_x']; p = P[key]; nl = np.asarray(P[f'{key}_null'], float)
        m, pz, nz, s = norm_rows(x, p, nl, xmax)
        _p_bands(ax, x[m], nz[:, m])
        ax.plot(x[m], pz[m], color=c, lw=1.2, zorder=3)
        b = J[key]
        lo = np.nanmin([np.nanpercentile(nz, 2.5, 0)[m].min(), np.nanmin(pz[m])])
        hi = np.nanmax([np.nanpercentile(nz, 97.5, 0)[m].max(), np.nanmax(pz[m])])
        ax.set_ylim(lo - 0.08 * (hi - lo), hi + 0.10 * (hi - lo))
        _p_style(ax, 1.0)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
        if r < 2:
            ax.tick_params(labelbottom=False)
        else:
            _p_xlabel(ax, -15, two_lines=True)
        kw = dict(xy=(0, 1), xycoords='axes fraction', textcoords='offset points', ha='left', va='bottom',
                  annotation_clip=False)
        ax.annotate(f'({"bcdefg"[i]}) {P_NAME[key]}', xytext=(0, 20.5), fontsize=P_FS['tag'], fontweight='bold',
                    color=TITLE, **kw)
        ax.annotate(f'ridge {fmt % b["A"]} ± {fmt % b["A_null_sig"]} {unit} · z = {b["z"]:+.1f}', xytext=(0, 11),
                    fontsize=P_FS['text'], color=FG, **kw)
        ax.annotate(f'σ_null = {fmt % s} {unit} · {b["n_null"]} {P_NULLDESC[key]}', xytext=(0, 2),
                    fontsize=P_FS['note'], color=FGDIM, **kw)
    # one y label for the six normalised panels
    yc = (bottom + grid_top) / 2 / P_H
    fig.text(0.07 / P_W, yc, '(I − ⟨I⟩_wings) / σ_null', rotation=90, ha='left', va='center',
             fontsize=P_FS['label'], color=FG)
    out = ss.save_paper(fig, FIG / 'w3_figs_f6_profiles.png', dpi=300)
    plt.close(fig)
    print('saved', out, out.replace('.png', '.pdf'))


if __name__ == '__main__':
    main()
