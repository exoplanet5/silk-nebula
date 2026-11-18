"""F1 overview map: NSNS [O III] background, Silk strands, Ha patch, W26 P1-P4, LAMOST field,
eRASS1 0.2-0.6 keV contours, Monogem/Geminga directions, SDSS and DESI spectroscopic coverage.

Inputs: data/derived/w3_figs_f1_maps.npz (w3_figs_f1_data.py), silk_trace.json v1.0,
        verify_sdss_fibres.csv (SDSS/BOSS sky fibres used by verify_sdss),
        verify_desi_coadd_meta.npz (DESI DR1 coadds used by verify_desi).
Output: figures/w3_figs_f1_overview.png (default: dark theme, 9.6 x 12 in, 200 dpi)
        SILK_THEME=light: negative [O III] map (white sky, dark-blue emission), overlays re-coloured for white
        SILK_PAPER=1: 7.1 x 8.3 in print layout (legend and colour bar below the map, no title or credit),
                      PNG (300 dpi) + PDF written to paper/figures/ via silk_style.save_paper
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_f1_overview.py
     SILK_THEME=light SILK_PAPER=1 ~/.venvs/astro313/bin/python -I scripts/w3_figs_f1_overview.py
"""
import sys
import csv
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import patheffects as pe
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, load_trace, small_circle, ang_sep

PSR = (104.9506, 14.2393)         # PSR B0656+14 (ATNF, via silk_trace.json candidate centres)
GEM = (98.4756, 17.7703)          # Geminga (now)
K2 = (104.15, 15.3)               # Knies+2024 two-SNR Monogem centre (via silk_trace.json)


def pfx(effects):
    """path_effects keyword only when set (keeps the dark default calls unchanged)."""
    return {} if not effects else dict(path_effects=effects)


def theme():
    """Overlay colours. Dark values are the original figure's; light values are chosen for a white page."""
    if ss.LIGHT:
        wh = pe.withStroke(linewidth=2.2, foreground=(1.0, 1.0, 1.0, 0.85))
        return dict(
            # negative [O III] map: white sky, dark-blue emission; the top end stays blue so stars do not dominate
            cmap=LinearSegmentedColormap.from_list('o3neg', [
                (0.00, '#ffffff'), (0.12, '#f3f7fb'), (0.30, '#c6dcef'), (0.50, '#6fa6d3'),
                (0.70, '#2b6cad'), (1.00, '#0c3a6e')]),
            title=TITLE, sdss='#a8790a', sdss_a=0.55, desi=NII_C, desi_a=0.55,
            xray='#e65100', xlab='#b23c00', xhalo=[ss.HALO], lamost=OK,
            strand=FG, strand_a=0.9, strand_lw=0.8, strand_pe=[wh],
            region=FG, region_pe=[pe.withStroke(linewidth=3.0, foreground=(1.0, 1.0, 1.0, 0.9))],
            arrow='#37474f', sdim=FGDIM, grid=FGDIM, grid_a=0.22, ha_a=0.25, legbg=BG)
    return dict(
        cmap=LinearSegmentedColormap.from_list('o3', ['#05070a', '#0b2233', '#125a7a', '#4fc3f7', '#e9f8ff']),
        title='#cfd6dd', sdss='#ffd54f', sdss_a=0.35, desi='#c792ea', desi_a=0.35,
        xray='#ff9e80', xlab='#ffb59e', xhalo=None, lamost=OK,
        strand='#ffffff', strand_a=0.75, strand_lw=0.7, strand_pe=None,
        region='#ffd54f', region_pe=None, arrow=WARN, sdim=FGDIM, grid='#9aa4ae', grid_a=0.15, ha_a=0.18, legbg=BG2)


def layout():
    """Sizes and limits. The paper layout is set at its 7.1 in print width, so point sizes are final."""
    if ss.PAPER:
        return dict(figsize=(7.1, 8.3), xlim=(124.5, 107.5), ylim=(10.6, 29.0),
                    ax_left_in=0.52, ax_width_in=6.50, ax_bottom_in=0.80,
                    tick=7.5, axlab=8.0, strand=7.5, send=7.0, ha=7.5, lam=7.0, reg=8.0, ms=6.0,
                    arrow=7.0, arrow_x=(109.3, 107.65), arrow_lab_x=107.75, clab=7.0, sdss_s=1.0, desi_s=0.8,
                    leg=7.0, cb_tick=7.0, cb_lab=7.5, regions_label='W26 P1–P4', xray_label='eRASS1 0.2–0.6 keV',
                    strands_label='strands E / W / W2', lam_xy=(120.45, 13.3),
                    cbar_in=(5.35, 0.25, 1.50, 0.075))
    return dict(figsize=(9.6, 12.0), xlim=(124.0, 108.0), ylim=(8.0, 30.0),
                tick=8.5, axlab=None, strand=8, send=7, ha=8, lam=7.5, reg=8.5, ms=6.5,
                arrow=7.3, arrow_x=(109.7, 108.05), arrow_lab_x=108.15, clab=6.5, sdss_s=1.2, desi_s=0.8,
                leg=7.3, cb_tick=7, cb_lab=7.5, regions_label='Wu+26 P1–P4',
                xray_label='eRASS1 0.2–0.6 keV (cts s⁻¹ deg⁻²)', strands_label='strands E / W / W2',
                lam_xy=(119.15, 13.05))


def main():
    th, L = theme(), layout()
    M = np.load(DER / 'w3_figs_f1_maps.npz')
    T = load_trace()
    ra, de, o3 = M['ra'], M['dec'], M['oiii'].astype(float)
    cmap = th['cmap']
    sm = ndimage.gaussian_filter(np.nan_to_num(o3, nan=0.0), 0.8)
    stretch = np.arcsinh(np.clip(sm, -0.3, 6) / 0.6) / np.arcsinh(4.0 / 0.6)
    asp = 1 / np.cos(np.radians(19.0))

    fig = plt.figure(figsize=L['figsize'])
    if not ss.PAPER:
        fig.text(0.012, 0.992, 'Silk Nebula · overview', color=th['title'], fontsize=13, fontweight='bold', va='top')
        fig.text(0.012, 0.972, 'NSNS DR0.2 [O III] 5007 · eRASS1 0.2–0.6 keV contours · spectroscopic coverage',
                 color=FGDIM, fontsize=9.5, va='top')
        ax = fig.add_axes([0.085, 0.07, 0.86, 0.885])
    else:
        # axes box sized to the data aspect so nothing is shrunk; legend + colour bar sit in a strip below
        W, H = L['figsize']
        (x0, x1), (y0, y1) = L['xlim'], L['ylim']
        h_in = L['ax_width_in'] * (y1 - y0) * asp / abs(x1 - x0)
        ax = fig.add_axes([L['ax_left_in'] / W, L['ax_bottom_in'] / H, L['ax_width_in'] / W, h_in / H])
    ext = [ra.max() + 0.5 / 60, ra.min() - 0.5 / 60, de.min() - 0.5 / 60, de.max() + 0.5 / 60]
    im = ax.imshow(stretch, origin='lower', extent=ext, cmap=cmap, vmin=0.0, vmax=1.0, interpolation='bilinear')
    ax.set_aspect(asp)
    ax.set_xlim(*L['xlim']); ax.set_ylim(*L['ylim'])

    # spectroscopic coverage (light)
    fib = list(csv.DictReader(open(DER / 'verify_sdss_fibres.csv')))
    fra = np.array([float(r['ra']) for r in fib]); fde = np.array([float(r['dec']) for r in fib])
    ax.scatter(fra, fde, s=L['sdss_s'], color=th['sdss'], alpha=th['sdss_a'], lw=0, zorder=2,
               rasterized=ss.PAPER)
    cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    dsel = cm['ok_contrib']
    ax.scatter(cm['ra'][dsel], cm['dec'][dsel], s=L['desi_s'], color=th['desi'], alpha=th['desi_a'], lw=0,
               zorder=2, rasterized=ss.PAPER)

    # eRASS1 0.2-0.6 keV contours
    X = M['xray']
    Xs = ndimage.gaussian_filter(np.nan_to_num(X, nan=np.nanmedian(X)), 1.5)   # extra 0.15 deg smoothing
    Xs[~np.isfinite(X)] = np.nan
    lev = [6.0, 7.0, 8.0]
    cs = ax.contour(M['x_ra'], M['x_dec'], Xs, levels=lev, colors=[th['xray']], linewidths=[0.7, 1.0, 1.3],
                    alpha=0.8, zorder=3)
    cl = ax.clabel(cs, levels=lev, fmt='%.0f', fontsize=L['clab'], colors=th['xlab'], inline_spacing=2)
    if th['xhalo']:
        for t in cl:
            t.set_path_effects(th['xhalo'])

    # LAMOST MRS-N plate field
    lr, ld, lrad = T['lamost']
    cr, cd = small_circle(lr, ld, lrad)
    ax.plot(cr, cd, color=th['lamost'], lw=1.0, ls=(0, (5, 3)), zorder=4)
    label(ax, *L['lam_xy'], 'LAMOST MRS-N\nplate (r = 2.5°)', color=th['lamost'], size=L['lam'], ha='left', va='top')

    # strands
    if not ss.PAPER:
        for k, c, lw in (('E', '#ffffff', 0.9), ('W', '#ffffff', 0.9), ('W2', '#ffffff', 0.8)):
            P = T[k]
            ax.plot(P[:, 0] + 0.45 / np.cos(np.radians(P[:, 1])) * (1 if k == 'E' else -1), P[:, 1],
                    color=c, lw=0.6, alpha=0.0)  # placeholder (no offset guides)
    # draw strands as thin dotted lines so the [O III] ridge stays visible
    for k in ('E', 'W', 'W2'):
        ax.plot(T[k][:, 0], T[k][:, 1], color=th['strand'], lw=th['strand_lw'], ls=(0, (1, 2.2)),
                alpha=th['strand_a'], zorder=5, **pfx(th['strand_pe']))
    e = T['E']
    label(ax, e[-1, 0] - 0.25, e[-1, 1] + 0.35, 'E strand (spine)', size=L['strand'], ha='left')
    label(ax, e[0, 0] + 0.35, e[0, 1] - 0.45, 'S end (116.54, +12.00)', size=L['send'], color=th['sdim'], ha='left')
    w = T['W']
    label(ax, w[-1, 0] - 0.55, w[-1, 1] + 0.05, 'W', size=L['strand'], ha='left')
    w2 = T['W2']
    label(ax, w2[len(w2) // 2, 0] - 0.35, w2[len(w2) // 2, 1], 'W2', size=L['strand'], ha='left')

    # Halpha patch
    pp = T['patch']
    ax.fill(pp[:, 0], pp[:, 1], facecolor=HA_C, alpha=th['ha_a'], edgecolor=HA_C, lw=1.0, zorder=5)
    label(ax, 114.55, 20.95, 'Hα patch', color=HA_C, size=L['ha'], ha='left')

    # W26 regions P1-P4
    for k, (r, d) in T['regions'].items():
        ax.plot(r, d, marker='s', ms=L['ms'], mfc='none', mec=th['region'], mew=1.3, zorder=7,
                **pfx(th['region_pe']))
        if k == 'P4':      # label to the west (display right)
            label(ax, r - 0.28, d, k, color=th['region'], size=L['reg'], ha='left', va='center')
        else:              # label to the east (display left), away from the bright ridge
            label(ax, r + 0.30, d, k, color=th['region'], size=L['reg'], ha='right', va='center')

    # directions to the Monogem pulsar and Geminga (both west of the map)
    # separations are the mean distances of the southern spine (verify_distance / silk_trace.json)
    for (rr, dd_), txt in ((PSR, 'PSR B0656+14 · 11.4°\nMonogem two-SNR\ncentre · 12.1°'),
                           (GEM, 'Geminga · 17.6°')):
        y0 = dd_
        ax.annotate('', xy=(L['arrow_x'][1], y0), xytext=(L['arrow_x'][0], y0),
                    arrowprops=dict(arrowstyle='-|>', color=th['arrow'], lw=1.2, mutation_scale=11), zorder=8)
        label(ax, L['arrow_lab_x'], y0 + 0.3, txt, color=th['arrow'], size=L['arrow'], ha='right', va='bottom')

    akw = {} if L['axlab'] is None else dict(fontsize=L['axlab'])
    ax.set_xlabel('RA (deg)', **akw); ax.set_ylabel('Dec (deg)', **akw)
    ax.grid(color=th['grid'], alpha=th['grid_a'], lw=0.5)
    ax.tick_params(labelsize=L['tick'])

    hand = [Line2D([], [], color=th['strand'], ls=(0, (1, 2.2)), lw=1.0, label=L['strands_label'],
                   **pfx(th['strand_pe'])),
            Patch(facecolor=HA_C, alpha=0.35, edgecolor=HA_C, label='Hα patch'),
            Line2D([], [], marker='s', ls='', mfc='none', mec=th['region'], mew=1.3, label=L['regions_label']),
            Line2D([], [], color=th['lamost'], ls=(0, (5, 3)), label='LAMOST plate',
                   **({'lw': 1.0} if ss.PAPER else {})),
            Line2D([], [], color=th['xray'], lw=1.0, label=L['xray_label']),
            Line2D([], [], marker='o', ls='', ms=3, color=th['sdss'], alpha=0.8 if ss.LIGHT else 0.7,
                   label='SDSS/BOSS sky fibres'),
            Line2D([], [], marker='o', ls='', ms=3, color=th['desi'], alpha=0.8 if ss.LIGHT else 0.7,
                   label='DESI DR1 coadds')]

    if not ss.PAPER:
        # colour bar (inset, top-left corner of the map)
        cax = ax.inset_axes([0.035, 0.955, 0.30, 0.013])
    else:
        # colour bar in the strip below the map, right of the legend
        W, H = L['figsize']
        cx, cy, cw, ch = L['cbar_in']
        cax = fig.add_axes([cx / W, cy / H, cw / W, ch / H])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    tk = np.array([0, 0.5, 1, 2, 4])
    cb.set_ticks(np.arcsinh(tk / 0.6) / np.arcsinh(4.0 / 0.6)); cb.set_ticklabels([str(t) for t in tk])
    cb.ax.tick_params(labelsize=L['cb_tick'], colors=FGDIM if not ss.LIGHT else FG, length=2)
    cb.set_label('NSNS [O III] (R)', fontsize=L['cb_lab'], color=FGDIM if not ss.LIGHT else FG, labelpad=1)
    if ss.PAPER:
        cb.ax.xaxis.set_label_position('top')
        cb.outline.set_linewidth(0.5)

    if not ss.PAPER:
        lg = ax.legend(handles=hand, loc='lower left', ncol=2, fontsize=L['leg'], framealpha=0.85, borderpad=0.5,
                       handlelength=1.8, columnspacing=1.0, bbox_to_anchor=(0.0, 0.0))
        lg.get_frame().set_facecolor(th['legbg'])
        ss.credit(fig)
        fig.savefig(FIG / 'w3_figs_f1_overview.png', dpi=200)
        out = FIG / 'w3_figs_f1_overview.png'
    else:
        W, H = L['figsize']
        # compact two-row legend in the strip below the map (outside the data)
        fig.legend(handles=hand, loc='lower left', bbox_to_anchor=(0.06 / W, 0.0), ncol=4, fontsize=L['leg'],
                   frameon=False, borderpad=0.3, handlelength=2.0, handletextpad=0.5, columnspacing=1.0,
                   labelspacing=0.45)
        plt.rcParams['savefig.dpi'] = 300          # resolution of the image / rasterised layers in the PDF
        out = ss.save_paper(fig, FIG / 'w3_figs_f1_overview.png', dpi=300)
    plt.close(fig)
    print('saved', out)


if __name__ == '__main__':
    main()
