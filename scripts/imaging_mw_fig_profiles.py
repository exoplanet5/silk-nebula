"""Figure: cross-filament profiles (offset from the E-strand spine) for key bands, P1/P2/P3/all, with P4 level.

Reads data/derived/imaging_mw_profiles.json. Output figures/imaging_mw_profiles.png
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FGDIM, OIII_C, HA_C, OK, WARN, ACCENT
from imaging_mw_common import DER, FIG

P = json.loads((DER / 'imaging_mw_profiles.json').read_text())
G = json.loads((DER / 'imaging_mw_galex_profiles.json').read_text())
KEY = [('murthy_fuv', 1, 'CU'), ('murthy_nuv', 1, 'CU'), ('rass_r12', 1, '10$^{-6}$ c s$^{-1}$ arcmin$^{-2}$'),
       ('rass_r45', 1, '10$^{-6}$ c s$^{-1}$ arcmin$^{-2}$'), ('ero_0206_clean', 1, 'c s$^{-1}$ deg$^{-2}$'),
       ('ero_0610_clean', 1, 'c s$^{-1}$ deg$^{-2}$'), ('iris100', 1, 'MJy sr$^{-1}$'), ('sfd_ebv', 1e3, 'mmag'),
       ('galfa_nhi_oiii', 1e-19, '10$^{19}$ cm$^{-2}$'), ('galfa_nhi_red', 1e-19, '10$^{19}$ cm$^{-2}$'),
       ('halpha_f03', 1, 'R'), ('stockert1420', 1, 'mK'), ('haslam408', 1, 'K'), ('planck030_pi', 1e6, '$\\mu$K')]
SEGC = {'all': '#e8eaed', 'P1': OIII_C, 'P2': OK, 'P3': WARN}

fig, axes = new_sky_figure(17, 13.5, 'Silk Nebula · cross-filament profiles',
                           'median in offset bins from the E-strand spine (+ = east); segments: whole spine and ±0.6° around P1–P3; dashed = P4 disk median',
                           nrows=4, ncols=4)
axes = axes.ravel()
for ax, (name, sc, unit) in zip(axes, KEY):
    b = P['bands'].get(name)
    if b is None:
        ax.axis('off'); continue
    for seg, col in SEGC.items():
        sg = b['segments'].get(seg)
        if not sg:
            continue
        x = np.array(sg['profile_offsets']); y = np.array(sg['profile'], float) * sc; e = np.array(sg['profile_err'], float) * sc
        lw = 1.6 if seg == 'all' else 1.0
        ax.plot(x, y, color=col, lw=lw, label=seg if name == KEY[0][0] else None)
        if seg == 'all':
            ax.fill_between(x, y - e, y + e, color=col, alpha=0.15, lw=0)
    ax.axhline(b['regions']['P4']['median'] * sc, color=FGDIM, ls='--', lw=0.8)
    ax.axvspan(-0.138, 0.138, color=OIII_C, alpha=0.12, lw=0)
    sa = b['segments']['all']
    txt = f"C = {sa['contrast'] * sc:+.3g} ({sa['contrast_sig']:+.1f}σ)\nW−E = {sa['step_W_minus_E'] * sc:+.3g} ({sa['step_sig']:+.1f}σ)"
    label(ax, 0.03, 0.04, txt, size=7, transform=ax.transAxes, va='bottom', color='#cfd6dd')
    ax.set_title(b['label'], loc='left', fontsize=9)
    ax.set_xlim(-4, 4)
    ax.set_xlabel('offset from spine (deg, + east)', fontsize=8)
    ax.set_ylabel(unit, fontsize=8)
    ax.grid(alpha=0.3)
# GALEX tile-level FUV & NUV (joint fit) in the last two panels
for ax, band in zip(axes[len(KEY):len(KEY) + 2], ['FUV', 'NUV']):
    for seg, col in [('all', '#e8eaed'), ('P1', OIII_C), ('P2', OK), ('P3', WARN)]:
        g = G['profiles'].get(f'{band}:{seg}')
        if not g:
            continue
        x = np.array(G['centres']); y = np.array(g['profile_CU'], float); e = np.array(g['err_CU'], float)
        ax.plot(x, y, color=col, lw=1.6 if seg == 'all' else 1.0)
        if seg == 'all':
            ax.fill_between(x, y - e, y + e, color=col, alpha=0.15, lw=0)
    ax.axvspan(-0.138, 0.138, color=OIII_C, alpha=0.12, lw=0)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.set_title(f'GALEX {band} tile skybg (offsets removed)', loc='left', fontsize=9)
    ax.set_xlim(-3, 3); ax.set_xlabel('offset from spine (deg, + east)', fontsize=8); ax.set_ylabel('CU', fontsize=8)
    ax.grid(alpha=0.3)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc='upper right', ncol=4, fontsize=8, bbox_to_anchor=(0.99, 0.985))
finish(fig, FIG / 'imaging_mw_profiles.png')
print('ok')
