"""w3_canon_fig.py -- canonical [O III] brightness and velocity ladders by dataset (wave 3, tag `canon`).

Values are copied from data/derived/w3_canon_values.json entries (sb.*, v.*), which quote the wave-2 notes.
Output: figures/w3_canon_ladder.png (dpi 200).
Run:  ~/.venvs/astro313/bin/python -I scripts/w3_canon_fig.py
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *  # noqa: E402,F401,F403
import numpy as np  # noqa: E402

OUT = '/Users/mickey/astronomy/silk-nebula/figures/w3_canon_ladder.png'

# ---- brightness rows: (label, centre, lo, hi, colour, filled, region)
SB = [
    ('Wu+26 SWONS P1-P3', 11.5, 9.94, 12.75, DANGER, True, 'S'),
    ('QHYCCD note', 5.0, 3.3, 6.7, FGDIM, True, 'S'),
    ('NSNS on-ridge P1-P3', 2.1, 2.04, 2.18, WARN, True, 'S'),
    ('NSNS 1\' peaks P1-P3', 3.3, 2.95, 3.64, WARN, False, 'S'),
    ('SDSS imaging 17\' (g)', 2.31, 1.93, 2.69, OK, True, 'S'),
    ('SDSS fibres plate peaks', 2.05, 1.6, 2.5, ACCENT, False, 'S'),
    ('SDSS fibres ON-OFF', 1.39, 1.29, 1.49, ACCENT, True, 'S'),
    ('DESI ridge (GFA cal.)', 1.62, 1.49, 1.75, OIII_C, True, 'S'),
    ('DESI north', 0.63, 0.56, 0.71, OIII_C, True, 'N'),
    ('SDSS fibres north peaks', 0.6, 0.4, 0.8, ACCENT, False, 'N'),
    ('NSNS northern arc', 0.60, 0.58, 0.62, WARN, True, 'N'),
]
SB_OLD = [('DESI wave 1 (bug)', 2.4, 'DESI ridge (GFA cal.)'), ('SDSS wave 1 (bug)', 0.79, 'SDSS fibres ON-OFF')]

# ---- velocity rows (V_helio): (label, centre, lo, hi, colour, filled)
V = [
    ('LAMOST Table 4 P1-P3', 25.8, 23.8, 26.8, DANGER, False),
    ('LAMOST corrected P1-P3', 15.7, 13.7, 16.7, DANGER, True),
    ('LAMOST Fig 3a vs [N I]', 13.6, 12.4, 14.8, DANGER, True),
    ('DESI on-ridge stack', 14.4, 13.2, 15.6, OIII_C, True),
    ('DESI joint regression', 11.3, 10.2, 13.9, OIII_C, True),
    ('DESI south Dec<19.5', 10.6, 10.2, 12.0, OIII_C, True),
    ('DESI north Dec>=19.5', 20.9, 18.5, 23.1, OIII_C, False),
    ('SDSS south', 10.0, 8.0, 12.0, ACCENT, True),
    ('SDSS north', 16.4, 10.5, 22.3, ACCENT, False),
    ('Paper [N II]/[S II]', 12.5, 11.7, 13.3, NII_C, True),
    ('WHAM WIM Halpha', 12.8, 11.9, 13.6, HA_C, True),
]


def main():
    fig, (a1, a2) = new_sky_figure(12.6, 6.2, 'Silk Nebula · canonical [O III] values', '', ncols=2)
    fig.subplots_adjust(left=0.165, right=0.985, bottom=0.10, top=0.83, wspace=0.62)
    # --- brightness
    y = np.arange(len(SB))[::-1]
    nN = sum(1 for r in SB if r[6] == 'N')
    a1.fill_betweenx([nN - 0.5, len(SB) - 0.5], 1.4, 2.5, color=OIII_C, alpha=0.13, lw=0)
    a1.fill_betweenx([-0.5, nN - 0.5], 0.4, 0.9, color=FGDIM, alpha=0.16, lw=0)
    a1.set_ylim(-0.5, len(SB) - 0.5)
    for yi, (lab, c, lo, hi, col, filled, reg) in zip(y, SB):
        a1.plot([lo, hi], [yi, yi], color=col, lw=2.2, solid_capstyle='butt')
        a1.plot(c, yi, 'o', ms=6, mec=col, mfc=col if filled else 'none', mew=1.5)
    labels = [r[0] for r in SB]
    for lab, val, row in SB_OLD:
        yi = y[labels.index(row)]
        a1.plot(val, yi, 'x', ms=6, color=FGDIM, mew=1.3)
        label(a1, val * 1.07, yi + 0.18, 'wave 1', color=FGDIM, size=6.5)
    a1.set_yticks(y)
    a1.set_yticklabels(labels, fontsize=7.5)
    a1.set_xscale('log')
    a1.set_xlim(0.3, 20)
    a1.set_xticks([0.3, 0.5, 1, 2, 5, 10, 20])
    a1.set_xticklabels(['0.3', '0.5', '1', '2', '5', '10', '20'])
    a1.set_xlabel('I(5007)  [R]')
    a1.axhline(y[labels.index('DESI north')] + 0.5, color='#2a323c', lw=0.8)
    label(a1, 2.55, len(SB) - 0.85, 'adopted south 1.4-2.5', color=OIII_C, size=7, va='center')
    label(a1, 0.92, -0.25, 'adopted north 0.4-0.9', color=FGDIM, size=7, va='center')
    a1.set_title('Surface brightness', loc='left', color=FG)
    a1.grid(axis='x', lw=0.4)

    # --- velocity
    y2 = np.arange(len(V))[::-1]
    a2.axvspan(13.3 - 1.0, 13.3 + 1.0, color=OK, alpha=0.18, lw=0)
    a2.axvline(13.3, color=OK, lw=0.8)
    for yi, (lab, c, lo, hi, col, filled) in zip(y2, V):
        a2.plot([lo, hi], [yi, yi], color=col, lw=2.2, solid_capstyle='butt')
        a2.plot(c, yi, 'o', ms=6, mec=col, mfc=col if filled else 'none', mew=1.5)
    a2.set_yticks(y2)
    a2.set_yticklabels([r[0] for r in V], fontsize=7.5)
    a2.set_xlim(5, 30)
    a2.set_xlabel('V_helio  [km/s]')
    label(a2, 12.4, y2[0] + 0.55, 'adopted +13.3 ± 1.0', color=OK, size=7.5, ha='right')
    a2.set_ylim(y2[-1] - 0.7, y2[0] + 0.9)
    a2.set_title('[O III] velocity', loc='left', color=FG)
    a2.grid(axis='x', lw=0.4)
    sec = a2.secondary_xaxis('top', functions=(lambda v: v - 12.1, lambda v: v + 12.1))
    sec.set_xlabel('V_LSR at P1  [km/s]', fontsize=8)
    finish(fig, OUT, tight=False)
    print(OUT)


if __name__ == '__main__':
    main()
