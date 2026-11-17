"""spectra_desi_fig_coverage.py -- DESI DR1 spectra density + tiles vs the Silk Nebula trace.
Output: figures/spectra_desi_coverage.png"""
import csv
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, label, finish, FG, FGDIM, OIII_C, WARN, SUN, HA_C, OK  # noqa: E402
from silk_geometry import ROOT, REGIONS, load_trace  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402


def draw_trace(ax, t, lw=1.2, alpha=0.9, envelope=True, color=OIII_C):
    sp = t['spine']
    ax.plot([p['ra'] for p in sp], [p['dec'] for p in sp], '-', color=color, lw=lw, alpha=alpha, label='E spine')
    for k, ls in (('W', '--'), ('W2', ':')):
        if k in t.get('strands', {}):
            pts = t['strands'][k]['points']
            ax.plot([p['ra'] for p in pts], [p['dec'] for p in pts], ls, color=color, lw=lw, alpha=alpha,
                    label=f'{k} strand')
    if envelope and 'envelope_polygon' in t:
        e = np.array(t['envelope_polygon'])
        ax.plot(np.r_[e[:, 0], e[0, 0]], np.r_[e[:, 1], e[0, 1]], '-', color=color, lw=0.5, alpha=0.4)


def draw_regions(ax, size=8):
    for k, (ra, dec) in REGIONS.items():
        ax.plot(ra, dec, 'x' if k != 'P4' else '+', color=SUN if k != 'P4' else FGDIM, ms=size, mew=1.5)
        label(ax, ra - 0.12, dec, k, color=SUN if k != 'P4' else FGDIM, size=8.5, va='center')


if __name__ == '__main__':
    t = load_trace()
    ra, dec, prog = [], [], []
    for r in csv.DictReader(open(ROOT / 'data/derived/spectra_desi_catalog.csv')):
        ra.append(float(r['ra']))
        dec.append(float(r['dec']))
        prog.append(r['program'])
    ra, dec = np.array(ra), np.array(dec)
    tiles = list(csv.DictReader(open(ROOT / 'data/raw/datalab_tap/tiles_box.csv')))
    fig, ax = new_sky_figure(8.2, 11.0, 'DESI DR1 spectra · Silk Nebula field',
                             f'{len(ra):,} SPARCL/zpix coadds · tiles (r = 1.6°) · [O III] trace')
    H, xe, ye = np.histogram2d(ra, dec, bins=[np.arange(111, 122.01, 0.2), np.arange(6, 28.01, 0.2)])
    dens = H.T / (0.2 * 0.2 * np.cos(np.radians(0.5 * (ye[:-1] + ye[1:])))[:, None])
    im = ax.imshow(np.where(dens > 0, dens, np.nan), origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]),
                   cmap='magma', vmin=0, vmax=4000, aspect='auto', interpolation='nearest')
    for tl in tiles:
        if tl['survey'] != 'main' or tl['program'] not in ('dark', 'bright'):
            continue
        c = Circle((float(tl['tilera']), float(tl['tiledec'])), 1.6 / np.cos(np.radians(float(tl['tiledec']))),
                   fill=False, lw=0.35, ec=('#9aa4ae' if tl['program'] == 'dark' else '#ffb84f'), alpha=0.35)
        ax.add_patch(c)
    draw_trace(ax, t, lw=1.6)
    draw_regions(ax)
    pil = list(csv.DictReader(open(ROOT / 'data/derived/spectra_desi_sel_pilot.csv')))
    ax.plot([float(r['ra']) for r in pil], [float(r['dec']) for r in pil], '.', ms=0.6, color=OK, alpha=0.6,
            label='pilot spectra')
    ax.set_xlim(122, 111)
    ax.set_ylim(6, 28)
    ax.set_aspect(1 / np.cos(np.radians(17)))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.grid(alpha=0.25, lw=0.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('spectra per deg²', color=FG)
    ax.plot([], [], '-', color='#9aa4ae', lw=0.8, label='main dark tile')
    ax.plot([], [], '-', color='#ffb84f', lw=0.8, label='main bright tile')
    ax.legend(loc='lower left', ncol=2, fontsize=7, markerscale=8)
    finish(fig, ROOT / 'figures/spectra_desi_coverage.png')
