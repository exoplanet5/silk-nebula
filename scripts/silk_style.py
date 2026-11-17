"""Shared figure style for the Silk Nebula project (dark palette, Menlo, haloed labels).

Palette and halo follow ~/sda/starship-ITF-14/work/ift14_style.py so sky figures match the
user's map figures. Usage:

    import sys; sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
    from silk_style import *
    fig, ax = new_sky_figure(10, 8, 'Silk Nebula · WHAM Hα', 'V_LSR −10…+30 km/s')
    ...
    finish(fig, 'figures/whams_halpha.png')
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe, font_manager

try:
    font_manager.fontManager.addfont('/System/Library/Fonts/Menlo.ttc')
    plt.rcParams['font.family'] = 'Menlo'
except Exception:
    pass

BG, BG2, BG3, FG, FGDIM = '#101418', '#171c22', '#1e242c', '#e8eaed', '#9aa4ae'
ACCENT, WARN, DANGER, OK, SUN = '#4fc3f7', '#ffb84f', '#ff5252', '#7ad97a', '#ffd54f'
OIII_C, HA_C, SII_C, NII_C, HB_C = '#4fc3f7', '#ff6f6f', '#ffb84f', '#c792ea', '#7ad97a'
HALO = pe.withStroke(linewidth=3, foreground=(5/255, 8/255, 12/255, 0.85))

plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': '#05070a', 'savefig.facecolor': BG,
    'axes.edgecolor': '#2a323c', 'axes.labelcolor': FG, 'text.color': FG,
    'xtick.color': FGDIM, 'ytick.color': FGDIM, 'grid.color': '#2a323c',
    'legend.facecolor': BG2, 'legend.edgecolor': '#2a323c', 'legend.fontsize': 8,
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
})


def new_sky_figure(w, h, title, subtitle='', nrows=1, ncols=1, **kw):
    """Figure with the title + subtitle at upper left; returns (fig, axes)."""
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h), **kw)
    fig.text(0.012, 0.985, title, color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    if subtitle:
        fig.text(0.012, 0.955, subtitle, color=FGDIM, fontsize=9.5, va='top')
    return fig, axes


def label(ax, x, y, text, color='#ffffff', size=8.5, **kw):
    return ax.text(x, y, text, color=color, fontsize=size, path_effects=[HALO], **kw)


def credit(fig, text='Illustration: Mickey'):
    fig.text(0.99, 0.006, text, color=FGDIM, fontsize=7, ha='right', va='bottom')


def finish(fig, path, dpi=200, tight=True):
    credit(fig)
    if tight:
        fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path
