"""F1 overview map: NSNS [O III] background, Silk strands, Ha patch, paper P1-P4, LAMOST field,
eRASS1 0.2-0.6 keV contours, Monogem/Geminga directions, SDSS and DESI spectroscopic coverage.

Inputs: data/derived/w3_figs_f1_maps.npz (w3_figs_f1_data.py), silk_trace.json v1.0,
        verify_sdss_fibres.csv (SDSS/BOSS sky fibres used by verify_sdss),
        verify_desi_coadd_meta.npz (DESI DR1 coadds used by verify_desi).
Output: figures/w3_figs_f1_overview.png
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_f1_overview.py
"""
import sys
import csv
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, load_trace, small_circle, ang_sep

PSR = (104.9506, 14.2393)         # PSR B0656+14 (ATNF, via silk_trace.json candidate centres)
GEM = (98.4756, 17.7703)          # Geminga (now)
K2 = (104.15, 15.3)               # Knies+2024 two-SNR Monogem centre (via silk_trace.json)


def main():
    M = np.load(DER / 'w3_figs_f1_maps.npz')
    T = load_trace()
    ra, de, o3 = M['ra'], M['dec'], M['oiii'].astype(float)
    cmap = LinearSegmentedColormap.from_list('o3', ['#05070a', '#0b2233', '#125a7a', '#4fc3f7', '#e9f8ff'])
    sm = ndimage.gaussian_filter(np.nan_to_num(o3, nan=0.0), 0.8)
    stretch = np.arcsinh(np.clip(sm, -0.3, 6) / 0.6) / np.arcsinh(4.0 / 0.6)

    fig = plt.figure(figsize=(9.6, 12.0))
    fig.text(0.012, 0.992, 'Silk Nebula · overview', color='#cfd6dd', fontsize=13, fontweight='bold', va='top')
    fig.text(0.012, 0.972, 'NSNS DR0.2 [O III] 5007 · eRASS1 0.2–0.6 keV contours · spectroscopic coverage',
             color=FGDIM, fontsize=9.5, va='top')
    ax = fig.add_axes([0.085, 0.07, 0.86, 0.885])
    ext = [ra.max() + 0.5 / 60, ra.min() - 0.5 / 60, de.min() - 0.5 / 60, de.max() + 0.5 / 60]
    im = ax.imshow(stretch, origin='lower', extent=ext, cmap=cmap, vmin=0.0, vmax=1.0, interpolation='bilinear')
    ax.set_aspect(1 / np.cos(np.radians(19.0)))
    ax.set_xlim(124.0, 108.0); ax.set_ylim(8.0, 30.0)

    # spectroscopic coverage (light)
    fib = list(csv.DictReader(open(DER / 'verify_sdss_fibres.csv')))
    fra = np.array([float(r['ra']) for r in fib]); fde = np.array([float(r['dec']) for r in fib])
    ax.scatter(fra, fde, s=1.2, color='#ffd54f', alpha=0.35, lw=0, zorder=2)
    cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    dsel = cm['ok_contrib']
    ax.scatter(cm['ra'][dsel], cm['dec'][dsel], s=0.8, color='#c792ea', alpha=0.35, lw=0, zorder=2)

    # eRASS1 0.2-0.6 keV contours
    X = M['xray']
    Xs = ndimage.gaussian_filter(np.nan_to_num(X, nan=np.nanmedian(X)), 1.5)   # extra 0.15 deg smoothing
    Xs[~np.isfinite(X)] = np.nan
    lev = [6.0, 7.0, 8.0]
    cs = ax.contour(M['x_ra'], M['x_dec'], Xs, levels=lev, colors=['#ff9e80'], linewidths=[0.7, 1.0, 1.3],
                    alpha=0.8, zorder=3)
    ax.clabel(cs, levels=lev, fmt='%.0f', fontsize=6.5, colors='#ffb59e', inline_spacing=2)

    # LAMOST MRS-N plate field
    lr, ld, lrad = T['lamost']
    cr, cd = small_circle(lr, ld, lrad)
    ax.plot(cr, cd, color=OK, lw=1.0, ls=(0, (5, 3)), zorder=4)
    label(ax, 119.15, 13.05, 'LAMOST MRS-N\nplate (r = 2.5°)', color=OK, size=7.5, ha='left', va='top')

    # strands
    for k, c, lw in (('E', '#ffffff', 0.9), ('W', '#ffffff', 0.9), ('W2', '#ffffff', 0.8)):
        P = T[k]
        ax.plot(P[:, 0] + 0.45 / np.cos(np.radians(P[:, 1])) * (1 if k == 'E' else -1), P[:, 1],
                color=c, lw=0.6, alpha=0.0)  # placeholder (no offset guides)
    # draw strands as thin dotted lines so the [O III] ridge stays visible
    ax.plot(T['E'][:, 0], T['E'][:, 1], color='#ffffff', lw=0.7, ls=(0, (1, 2.2)), alpha=0.75, zorder=5)
    ax.plot(T['W'][:, 0], T['W'][:, 1], color='#ffffff', lw=0.7, ls=(0, (1, 2.2)), alpha=0.75, zorder=5)
    ax.plot(T['W2'][:, 0], T['W2'][:, 1], color='#ffffff', lw=0.7, ls=(0, (1, 2.2)), alpha=0.75, zorder=5)
    e = T['E']
    label(ax, e[-1, 0] - 0.25, e[-1, 1] + 0.35, 'E strand (spine)', size=8, ha='left')
    label(ax, e[0, 0] + 0.35, e[0, 1] - 0.45, 'S end (116.54, +12.00)', size=7, color=FGDIM, ha='left')
    w = T['W']
    label(ax, w[-1, 0] - 0.55, w[-1, 1] + 0.05, 'W', size=8, ha='left')
    w2 = T['W2']
    label(ax, w2[len(w2) // 2, 0] - 0.35, w2[len(w2) // 2, 1], 'W2', size=8, ha='left')

    # Halpha patch
    pp = T['patch']
    ax.fill(pp[:, 0], pp[:, 1], facecolor=HA_C, alpha=0.18, edgecolor=HA_C, lw=1.0, zorder=5)
    label(ax, 114.55, 20.95, 'Hα patch', color=HA_C, size=8, ha='left')

    # paper regions
    for k, (r, d) in T['regions'].items():
        ax.plot(r, d, marker='s', ms=6.5, mfc='none', mec='#ffd54f', mew=1.3, zorder=7)
        if k == 'P4':      # label to the west (display right)
            label(ax, r - 0.28, d, k, color='#ffd54f', size=8.5, ha='left', va='center')
        else:              # label to the east (display left), away from the bright ridge
            label(ax, r + 0.30, d, k, color='#ffd54f', size=8.5, ha='right', va='center')

    # directions to the Monogem pulsar and Geminga (both west of the map)
    # separations are the mean distances of the southern spine (verify_distance / silk_trace.json)
    for (rr, dd_), txt in ((PSR, 'PSR B0656+14 · 11.4°\nMonogem two-SNR\ncentre · 12.1°'),
                           (GEM, 'Geminga · 17.6°')):
        y0 = dd_
        ax.annotate('', xy=(108.05, y0), xytext=(109.7, y0),
                    arrowprops=dict(arrowstyle='-|>', color=WARN, lw=1.2, mutation_scale=11), zorder=8)
        label(ax, 108.15, y0 + 0.3, txt, color=WARN, size=7.3, ha='right', va='bottom')

    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    ax.grid(color='#9aa4ae', alpha=0.15, lw=0.5)
    ax.tick_params(labelsize=8.5)

    # colour bar (inset, top-left corner of the map)
    cax = ax.inset_axes([0.035, 0.955, 0.30, 0.013])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    tk = np.array([0, 0.5, 1, 2, 4])
    cb.set_ticks(np.arcsinh(tk / 0.6) / np.arcsinh(4.0 / 0.6)); cb.set_ticklabels([str(t) for t in tk])
    cb.ax.tick_params(labelsize=7, colors=FGDIM, length=2)
    cb.set_label('NSNS [O III] (R)', fontsize=7.5, color=FGDIM, labelpad=1)

    hand = [Line2D([], [], color='#ffffff', ls=(0, (1, 2.2)), lw=1.0, label='strands E / W / W2'),
            Patch(facecolor=HA_C, alpha=0.35, edgecolor=HA_C, label='Hα patch'),
            Line2D([], [], marker='s', ls='', mfc='none', mec='#ffd54f', mew=1.3, label='Wu+26 P1–P4'),
            Line2D([], [], color=OK, ls=(0, (5, 3)), label='LAMOST plate'),
            Line2D([], [], color='#ff9e80', lw=1.0, label='eRASS1 0.2–0.6 keV (cts s⁻¹ deg⁻²)'),
            Line2D([], [], marker='o', ls='', ms=3, color='#ffd54f', alpha=0.7, label='SDSS/BOSS sky fibres'),
            Line2D([], [], marker='o', ls='', ms=3, color='#c792ea', alpha=0.7, label='DESI DR1 coadds')]
    lg = ax.legend(handles=hand, loc='lower left', ncol=2, fontsize=7.3, framealpha=0.85, borderpad=0.5,
                   handlelength=1.8, columnspacing=1.0, bbox_to_anchor=(0.0, 0.0))
    lg.get_frame().set_facecolor(BG2)
    ss.credit(fig)
    fig.savefig(FIG / 'w3_figs_f1_overview.png', dpi=200)
    plt.close(fig)
    print('saved', FIG / 'w3_figs_f1_overview.png')


if __name__ == '__main__':
    main()
