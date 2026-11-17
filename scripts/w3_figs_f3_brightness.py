"""F3 [O III] 5007 brightness reconciliation (Rayleigh, log axis), with aperture / scale notes.

Values (wave-2-corrected; wave-1 only where never re-tested, and marked so):
  Wu+26 Table 3 -> R: notes/paper_review.md (12.8 / 9.9 / 10.1 R; 8.3-14.9 R for 50-70 A effective width)
  QHYCCD press note 24.3 mag/arcsec2 -> 3.3-6.7 R (3-6 nm): notes/discovery.md
  Wu+26 MRS-N relative x sky-line scale: notes/spectra_sdss_lamost.md (wave 1; 2.8-3.2 R, 16-84 % 1.2-8.6)
  NSNS DR0.2: notes/imaging_optical.md (on-ridge |d|<=0.25 deg medians 2.04-2.18 R; 1' peaks 2.95-3.64 R)
  SDSS imaging sky (g): notes/verify_precovery.md (PSP 2.31+-0.38; Frames 2.55+-0.53; all confounders 1.65+-0.31)
  DESI: notes/verify_desi.md + verify_desi_brightness.json (south median 1.62, variants 1.60-1.71, +-8 % sys;
        E nodes Dec 14-18 1.44-2.30; W@18 3.88+-0.20; north 0.56-0.71)
  SDSS fibres (total sky): notes/verify_sdss.md (south ON-OFF 1.39+-0.07+-0.07; W 1.93+-0.11; per-plate peaks
        1.6-2.5 S, 0.4-0.8 N; P1/P2/P3 means 1.49/1.78/1.87 +-0.2)
  LAMOST LRS stellar residuals: notes/spectra_sdss_lamost.md (wave 1, not re-verified; 0.88+-0.13)
Output: figures/w3_figs_f3_brightness.png
"""
import sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import FIG

PAPER_C = DANGER
ROWS = [
    # group, label, scale note, items: list of (kind, values, colour, marker, filled)
    ('Imaging photometry', 'Wu+26 Table 3 · P1 / P2 / P3', 'star-free sub-regions',
     [('range', (8.3, 14.9), PAPER_C, None, True), ('pts', (12.8, 9.9, 10.1), PAPER_C, 'o', True)]),
    ('Imaging photometry', 'QHYCCD note · brightest region', 'V-calibrated, 3–6 nm',
     [('range', (3.3, 6.7), WARN, None, True), ('pts', (5.6,), WARN, 'o', False)]),
    ('Imaging photometry', 'NSNS DR0.2 · ridge median / 1′ peaks', '51″ pixels',
     [('range', (2.04, 2.18), OIII_C, None, True), ('pts', (2.10,), OIII_C, 'o', True),
      ('pts', (2.95, 3.64), OIII_C, '^', False)]),
    ('Imaging photometry', 'SDSS imaging sky g · 2002–07', '17′ ridge fit',
     [('range', (1.65, 2.55), OIII_C, None, True), ('err', (2.31, 0.38), OIII_C, 'D', True)]),
    ('Fibre spectroscopy', 'DESI DR1 · south ridge median', '1.5″ fibres',
     [('err', (1.62, 0.13), OIII_C, 'o', True)]),
    ('Fibre spectroscopy', 'DESI DR1 · E nodes Dec 14–18 / W @ Dec 18', '1.5″ fibres',
     [('range', (1.44, 2.30), OIII_C, None, True), ('err', (3.88, 0.20), OIII_C, 's', True)]),
    ('Fibre spectroscopy', 'DESI DR1 · north ridge', '1.5″ fibres',
     [('range', (0.56, 0.71), '#b8ecff', None, True), ('pts', (0.63,), '#b8ecff', 'o', True)]),
    ('Fibre spectroscopy', 'SDSS/BOSS · south ON−OFF / W strand', '2–3″ fibres, r_on 0.15°',
     [('err', (1.39, 0.10), OIII_C, 'o', True), ('err', (1.93, 0.11), OIII_C, 's', True)]),
    ('Fibre spectroscopy', 'SDSS/BOSS · P1 / P2 / P3 (0.3° means)', '2–3″ fibres',
     [('pts', (1.49, 1.78, 1.87), OIII_C, 'o', True)]),
    ('Fibre spectroscopy', 'SDSS/BOSS · per-plate peaks S / N', '2–3″ fibres',
     [('range', (1.6, 2.5), OIII_C, None, True), ('range', (0.4, 0.8), '#b8ecff', None, True)]),
    ('Fibre spectroscopy', 'LAMOST LRS star residuals (wave 1)', '3.3″, |d| < 0.2°, Dec 12–28',
     [('err', (0.88, 0.13), FGDIM, 'o', True)]),
    ('Fibre spectroscopy', 'Wu+26 MRS-N × sky-line scale (wave 1)', 'relative flux → R via [N I]',
     [('range', (1.2, 8.6), FGDIM, None, False), ('pts', (2.8, 3.1, 3.2), FGDIM, 'o', True)]),
]


def main():
    fig = plt.figure(figsize=(11.8, 7.6))
    fig.text(0.012, 0.985, 'Silk Nebula · [O III] 5007 brightness reconciliation', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.952, 'peak / ridge intensity in Rayleigh; bars = ranges or ±1σ; scale of each measurement at right',
             color=FGDIM, fontsize=9.5, va='top')
    ax = fig.add_axes([0.33, 0.095, 0.47, 0.79])
    n = len(ROWS)
    y = np.zeros(n); off = 0; prev = None
    for i, r in enumerate(ROWS):
        if prev is not None and r[0] != prev:
            off += 0.9
        y[i] = -(i + off); prev = r[0]
    ax.axvspan(1.4, 2.5, color=OK, alpha=0.13, lw=0, zorder=0)
    ax.axvspan(10, 13, color=DANGER, alpha=0.13, lw=0, zorder=0)
    for i, (g, lab, scale, items) in enumerate(ROWS):
        yy = y[i]
        for kind, vals, c, mk, filled in items:
            if kind == 'range':
                ax.plot(vals, [yy, yy], color=c, lw=5 if filled else 1.2, alpha=0.45 if filled else 0.8,
                        solid_capstyle='butt', zorder=2, ls='-' if filled else (0, (3, 2)))
            elif kind == 'pts':
                for v in vals:
                    ax.plot(v, yy, marker=mk, ms=6.5, color=c if filled else 'none', mec=c, mew=1.3, zorder=4,
                            ls='')
            elif kind == 'err':
                v, e = vals
                ax.errorbar(v, yy, xerr=e, fmt=mk, color=c, ms=6.5, lw=1.4, capsize=2.5, mec='k', mew=0.4, zorder=5)
        ax.text(-0.012, yy, lab, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=8.7, color=FG)
        ax.text(1.012, yy, scale, transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=8.0,
                color=FGDIM)
    prev = None
    for i, r in enumerate(ROWS):
        if r[0] != prev:
            ax.text(-0.012, y[i] + 0.78, r[0], transform=ax.get_yaxis_transform(), ha='right', va='center',
                    fontsize=9.2, color=WARN, fontweight='bold')
            prev = r[0]
    ax.text(1.012, y[0] + 0.78, 'scale', transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=9.2,
            color=WARN, fontweight='bold')
    ax.set_xscale('log')
    ax.set_xlim(0.3, 20)
    ax.set_xticks([0.3, 0.5, 1, 2, 3, 5, 10, 20]); ax.set_xticklabels(['0.3', '0.5', '1', '2', '3', '5', '10', '20'])
    ax.set_ylim(y.min() - 0.8, y.max() + 1.5)
    ax.set_yticks([])
    ax.set_xlabel('[O III] 5007 intensity (R)')
    ax.grid(axis='x', alpha=0.25, which='both')
    label(ax, 1.87, y.max() + 1.05, 'south consensus\n1.4–2.5 R', color=OK, size=8, ha='center', va='center')
    label(ax, 11.4, y.max() + 1.05, 'paper\n10–13 R', color=DANGER, size=8, ha='center', va='center')
    yg = 0.5 * (y[3] + y[4]) + 0.15
    ax.annotate('', xy=(2.6, yg), xytext=(9.8, yg),
                arrowprops=dict(arrowstyle='<->', color=FGDIM, lw=0.9, mutation_scale=9))
    label(ax, 5.05, yg - 0.32, 'paper / consensus ≈ 5–8', color=FGDIM, size=8, ha='center', va='center')
    ss.credit(fig)
    fig.savefig(FIG / 'w3_figs_f3_brightness.png', dpi=200)
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    main()
