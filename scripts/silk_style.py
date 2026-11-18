"""Shared figure style for the Silk Nebula project (Menlo, haloed labels; dark or light theme).

Default (no environment variables): the original dark palette, which follows
~/sda/starship-ITF-14/work/ift14_style.py so sky figures match the user's map figures.

Paper mode, used for the comment manuscript in paper/:
    SILK_THEME=light   light palette (white background, dark text, saturated line colours)
    SILK_PAPER=1       no credit line, no figure-level title/subtitle (the caption carries them),
                       finish() writes PNG (300 dpi) + PDF into $SILK_FIGDIR instead of figures/
    SILK_FIGDIR=...    output directory for paper figures (default: <project>/paper/figures)

Usage:

    import sys; sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
    from silk_style import *
    fig, ax = new_sky_figure(10, 8, 'Silk Nebula · WHAM Hα', 'V_LSR −10…+30 km/s')
    ...
    finish(fig, 'figures/whams_halpha.png')
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe, font_manager

try:
    font_manager.fontManager.addfont('/System/Library/Fonts/Menlo.ttc')
    plt.rcParams['font.family'] = 'Menlo'
except Exception:
    pass

THEME = os.environ.get('SILK_THEME', 'dark').lower()
LIGHT = THEME == 'light'
PAPER = os.environ.get('SILK_PAPER', '') == '1'
PROJECT = Path('/Users/mickey/astronomy/silk-nebula')
PAPER_FIGDIR = Path(os.environ.get('SILK_FIGDIR', str(PROJECT / 'paper/figures')))

if LIGHT:
    # light theme: white page, near-black text, colours darkened for contrast on white
    BG, BG2, BG3, FG, FGDIM = '#ffffff', '#f4f6f8', '#e9edf1', '#1b1f24', '#59636e'
    AXBG, EDGE, GRID = '#ffffff', '#8a949e', '#e1e6eb'
    TITLE = '#1b1f24'
    ACCENT, WARN, DANGER, OK, SUN = '#0277bd', '#d17a00', '#d32f2f', '#2e8b3a', '#b8860b'
    OIII_C, HA_C, SII_C, NII_C, HB_C = '#0277bd', '#d32f2f', '#d17a00', '#7b3fa0', '#2e8b3a'
    HALO = pe.withStroke(linewidth=3, foreground=(1.0, 1.0, 1.0, 0.85))
    LABEL_C = '#1b1f24'
else:
    BG, BG2, BG3, FG, FGDIM = '#101418', '#171c22', '#1e242c', '#e8eaed', '#9aa4ae'
    AXBG, EDGE, GRID = '#05070a', '#2a323c', '#2a323c'
    TITLE = '#cfd6dd'
    ACCENT, WARN, DANGER, OK, SUN = '#4fc3f7', '#ffb84f', '#ff5252', '#7ad97a', '#ffd54f'
    OIII_C, HA_C, SII_C, NII_C, HB_C = '#4fc3f7', '#ff6f6f', '#ffb84f', '#c792ea', '#7ad97a'
    HALO = pe.withStroke(linewidth=3, foreground=(5/255, 8/255, 12/255, 0.85))
    LABEL_C = '#ffffff'

plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': AXBG, 'savefig.facecolor': BG,
    'axes.edgecolor': EDGE, 'axes.labelcolor': FG, 'text.color': FG,
    'xtick.color': FGDIM if not LIGHT else FG, 'ytick.color': FGDIM if not LIGHT else FG,
    'grid.color': GRID,
    'legend.facecolor': BG2, 'legend.edgecolor': EDGE if LIGHT else '#2a323c', 'legend.fontsize': 8,
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
})


if PAPER:
    # Journal figures: embed TrueType (Type 42) fonts so PDF text stays extractable.
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # AJ house style uses American spelling; normalize every figure string in paper mode.
    import re as _re
    from matplotlib import text as _mtext
    _US = {
        'fibre': 'fiber', 'fibres': 'fibers', 'centre': 'center', 'centres': 'centers',
        'centred': 'centered', 'anticentre': 'anticenter', 'colour': 'color', 'colours': 'colors',
        'colourbar': 'colorbar', 'normalised': 'normalized', 'normalisation': 'normalization',
        'labelled': 'labeled', 'modelled': 'modeled', 'artefact': 'artifact', 'artefacts': 'artifacts',
        'grey': 'gray', 'favoured': 'favored', 'favours': 'favors', 'analysed': 'analyzed',
        'catalogue': 'catalog', 'catalogued': 'cataloged',
    }
    _US_RE = _re.compile(r'\b(' + '|'.join(sorted(_US, key=len, reverse=True)) + r')\b', _re.I)

    def _american(s):
        def rep(m):
            w = m.group(0)
            r = _US[w.lower()]
            return r.upper() if w.isupper() else (r.capitalize() if w[0].isupper() else r)
        return _US_RE.sub(rep, s)

    _orig_set_text = _mtext.Text.set_text

    def _set_text_us(self, s):
        if isinstance(s, str):
            s = _american(s)
        return _orig_set_text(self, s)

    _mtext.Text.set_text = _set_text_us


def new_sky_figure(w, h, title, subtitle='', nrows=1, ncols=1, **kw):
    """Figure with the title + subtitle at upper left (omitted in paper mode); returns (fig, axes)."""
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h), **kw)
    if not PAPER:
        fig.text(0.012, 0.985, title, color=TITLE, fontsize=12, fontweight='bold', va='top')
        if subtitle:
            fig.text(0.012, 0.955, subtitle, color=FGDIM, fontsize=9.5, va='top')
    return fig, axes


def label(ax, x, y, text, color=None, size=8.5, **kw):
    return ax.text(x, y, text, color=LABEL_C if color is None else color, fontsize=size,
                   path_effects=[HALO], **kw)


def credit(fig, text='Illustration: Mickey'):
    if PAPER:
        return None
    fig.text(0.99, 0.006, text, color=FGDIM, fontsize=7, ha='right', va='bottom')


def paper_path(path):
    """In paper mode, map an output path to $SILK_FIGDIR/<basename>; otherwise return it unchanged."""
    if not PAPER:
        return str(path)
    PAPER_FIGDIR.mkdir(parents=True, exist_ok=True)
    return str(PAPER_FIGDIR / Path(path).name)


def save_paper(fig, path, dpi=300, **kw):
    """Save PNG (dpi) and a PDF twin in paper mode; plain savefig otherwise. Returns the PNG path."""
    out = paper_path(path)
    fig.savefig(out, dpi=dpi, **kw)
    if PAPER:
        fig.savefig(str(Path(out).with_suffix('.pdf')), **kw)
    return out


def finish(fig, path, dpi=200, tight=True):
    credit(fig)
    if tight:
        fig.tight_layout(rect=(0, 0, 1, 1) if PAPER else (0, 0.02, 1, 0.94))
    out = save_paper(fig, path, dpi=300 if PAPER else dpi)
    plt.close(fig)
    return out
