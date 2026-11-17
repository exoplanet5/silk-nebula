"""F7 precovery and discovery timeline 2000-2026.

Inputs:
  verify_sdss_precovery.csv / _placebo.csv : per plate-MJD [O III] ON-OFF S/N including a per-plate placebo (verify_sdss)
  verify_precovery_table.json : SDSS imaging g-sky per-run S/N (PSP, gauss028 template), per-plate fibre S/N from
      wave-1 residual fluxes (used only for plates verify_sdss could not placebo-test: 927, 1263), LAMOST LRS per
      plan-date S/N, DESI DR1 per-month S/N (partly circular: regressed on a forward model fitted to all epochs)
  discovery_timeline.json / notes/discovery.md / notes/verify_precovery.md : dated events
Output: figures/w3_figs_f7_timeline.png
"""
import sys
import csv
import json
import datetime as dt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG

CAP = 19.0


def dyear(s):
    p = [int(x) for x in s.split('-')]
    d = dt.date(p[0], p[1] if len(p) > 1 else 7, p[2] if len(p) > 2 else 15)
    y0 = dt.date(d.year, 1, 1); y1 = dt.date(d.year + 1, 1, 1)
    return d.year + (d - y0).days / (y1 - y0).days


def main():
    T = json.loads((DER / 'verify_precovery_table.json').read_text())
    pc = {(r['plate'], r['mjd']): r for r in csv.DictReader(open(DER / 'verify_sdss_precovery.csv'))}
    pl = list(csv.DictReader(open(DER / 'verify_sdss_precovery_placebo.csv')))
    fig = plt.figure(figsize=(14.5, 10.4))
    fig.text(0.012, 0.988, 'Silk Nebula · precovery and discovery timeline', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.963, '[O III] signal on the strand per epoch (S/N) · data epochs · discovery and publication record',
             color=FGDIM, fontsize=9.5, va='top')
    ax1 = fig.add_axes([0.06, 0.50, 0.92, 0.42])
    ax2 = fig.add_axes([0.06, 0.345, 0.92, 0.12], sharex=ax1)
    ax3 = fig.add_axes([0.06, 0.06, 0.92, 0.215])

    # ---------------- top: S/N per epoch
    ax1.axhspan(-3, 3, color=FGDIM, alpha=0.08, lw=0)
    ax1.axhline(3, color=FGDIM, lw=0.8, ls='--'); ax1.axhline(5, color=FGDIM, lw=0.6, ls=':')
    ax1.text(2000.15, 3.15, '3σ', color=FGDIM, fontsize=8, va='bottom'); ax1.text(2000.15, 5.15, '5σ', color=FGDIM,
                                                                                    fontsize=8, va='bottom')
    im = T['sdss_imaging']['g_psp_gauss028_OIIIunits']
    for r in im:
        core = 12 <= r['dec_cross'] <= 20.5
        ax1.plot(dyear(r['date']), r['snr'], 'D', ms=5.5 if core else 4.5, color=SUN if core else 'none',
                 mec=SUN if core else '#8a7a40', mew=1.0, zorder=5)
    # SDSS/BOSS plates (verify_sdss, placebo-included S/N)
    for r in pl:
        s = float(r['snr_placebo'])
        if not np.isfinite(s):
            continue
        d = pc[(r['plate'], r['mjd'])]['date_first_ut']
        ax1.plot(dyear(d), min(s, CAP), 'o', ms=6, color=OIII_C, mec='k', mew=0.4, zorder=6)
        if s > CAP:
            ax1.annotate('', xy=(dyear(d), CAP + 1.2), xytext=(dyear(d), CAP),
                         arrowprops=dict(arrowstyle='-|>', color=OIII_C, lw=1, mutation_scale=8))
    # plates without placebo test (verify_precovery, wave-1 per-fibre residuals)
    for r in T['sdss_fibres']:
        if r['plate_mjd'] in ('927-52577', '1263-52708'):
            ax1.plot(dyear(r['date']), r['snr'], 'o', ms=6.5, color='none', mec=OIII_C, mew=1.3, zorder=6)
    # LAMOST LRS
    for r in T['lamost_lrs']:
        if r['snr'] is not None:
            ax1.plot(dyear(r['date']), r['snr'], 's', ms=4.8, color=OK, alpha=0.85, mec='k', mew=0.3, zorder=4)
    # DESI DR1 months
    for r in T['desi']:
        y = dyear(r['month'] + '-15')
        ax1.plot(y, min(r['snr'], CAP), 'P', ms=8, color='#c792ea', mec='k', mew=0.4, zorder=6)
        if r['snr'] > CAP:
            ax1.annotate('', xy=(y, CAP + 1.2), xytext=(y, CAP),
                         arrowprops=dict(arrowstyle='-|>', color='#c792ea', lw=1, mutation_scale=8))
            label(ax1, y, CAP + 1.35, '%.0f' % r['snr'], color='#c792ea', size=7, ha='center', va='bottom')
    # precovery milestones
    ms = [(dyear('2002-10-30'), 5.0, 'plate 927 · 2002-10-30\n(3 fibres, no placebo)', 'right'),
          (dyear('2002-12-31'), 3.26, 'imaging run 3560\n2002-12-31', 'left'),
          (dyear('2003-01-30'), 3.38, 'plate 1203 · 2003-01-30\n(two analyses)', 'left'),
          (dyear('2005-03-03'), 9.68, 'plate 2074 · 2005-03-03\nfirst ≥ 5σ', 'left')]
    offs = [(-0.25, 4.0), (2.4, -5.5), (1.0, 2.3), (0.45, 2.6)]
    for (x, y, t, ha), (dx, dy) in zip(ms, offs):
        ax1.annotate(t, xy=(x, y), xytext=(x + dx, y + dy), fontsize=7.4, color=FG, ha=ha, va='center',
                     arrowprops=dict(arrowstyle='-', color=FGDIM, lw=0.7), path_effects=[HALO], zorder=8)
    ax1.set_ylim(-3.5, CAP + 3.0)
    ax1.set_xlim(2000.0, 2027.0)
    ax1.set_ylabel('S/N of on-strand [O III]')
    ax1.grid(alpha=0.2)
    plt.setp(ax1.get_xticklabels(), visible=False)
    hand = [Line2D([], [], marker='D', ls='', color=SUN, mec=SUN, label='SDSS imaging g sky, run (Dec 12–20.5)'),
            Line2D([], [], marker='D', ls='', color='none', mec='#8a7a40', label='same, north (broadband confounder)'),
            Line2D([], [], marker='o', ls='', color=OIII_C, mec='k', label='SDSS/BOSS plate (placebo-included)'),
            Line2D([], [], marker='o', ls='', color='none', mec=OIII_C, mew=1.3, label='plate, no placebo test'),
            Line2D([], [], marker='s', ls='', color=OK, mec='k', label='LAMOST LRS plan (star residuals)'),
            Line2D([], [], marker='P', ls='', color='#c792ea', mec='k', ms=8, label='DESI DR1 month (partly circular)')]
    ax1.legend(handles=hand, loc='upper left', ncol=2, fontsize=7.6, framealpha=0.9, bbox_to_anchor=(0.0, 1.0))

    # ---------------- middle: data epochs
    bars = [('SDSS imaging', '2001-11-12', '2007-02-21', SUN),
            ('SDSS-I/II spectra', '2001-11-19', '2008-02-10', OIII_C),
            ('BOSS / eBOSS plates', '2010-12-06', '2019-02-01', OIII_C),
            ('LAMOST LRS', '2011-11-09', '2024-01-13', OK),
            ('DESI DR1', '2021-01-03', '2022-02-28', '#c792ea'),
            ('NSNS DR0.2 data', '2022-12-12', '2025-04-04', '#b8ecff')]
    for i, (nm, a, b, c) in enumerate(bars):
        y = len(bars) - 1 - i
        ax2.plot([dyear(a), dyear(b)], [y, y], color=c, lw=6, solid_capstyle='butt', alpha=0.85)
        if nm.startswith('NSNS'):
            ax2.text(dyear(a) - 0.15, y, nm, color=c, fontsize=7.6, va='center', ha='right', path_effects=[HALO])
        else:
            ax2.text(dyear(b) + 0.15, y, nm, color=c, fontsize=7.6, va='center', ha='left', path_effects=[HALO])
    ax2.set_ylim(-0.8, len(bars) - 0.2)
    ax2.set_yticks([])
    ax2.axvspan(2022.75, 2027.0, color=FGDIM, alpha=0.08, lw=0)
    ax2.grid(axis='x', alpha=0.2)
    ax2.set_xticks(range(2000, 2028, 2))
    ax2.tick_params(labelsize=8.5)

    # ---------------- bottom: discovery record (zoom 2022.75-2027)
    # (date, text, colour, label level, label x)
    ev = [('2022-12', 'NSNS first [O III]\ncamera', OK, 1.2, 2022.98),
          ('2023-11', 'SWONS detects the\nfilament (self-reported)', OIII_C, 2.6, 2023.50),
          ('2024-01', 'NSNS full array:\nindependent detection', OK, 4.4, 2023.62),
          ('2024-01-14', 'AstroBin zh47r9:\nfirst public image', OIII_C, 4.4, 2024.15),
          ('2024-02-01', 'LAMOST MRS-N\nplate (paper data)', OIII_C, 1.2, 2024.20),
          ('2024-06-11', 'Cloudy Nights post\n(Ziegenbalg)', OK, 2.8, 2024.48),
          ('2024-07-12', 'QHYCCD note:\n"Silk Nebula"', OIII_C, 1.2, 2024.76),
          ('2025-05-22', 'NSNS DR0.2', OK, 1.8, 2025.38),
          ('2025-08-28', 'RNAAS 9, 227\n(Ziegenbalg)', OK, 3.6, 2025.66),
          ('2026-08-27', 'Fesen+ arXiv:\nMonogem limb', WARN, 2.4, 2026.42),
          ('2026-10-09', 'Wu+ arXiv\n2610.13329', OIII_C, 4.2, 2026.62),
          ('2026-11-04', 'AJ 172, 298', OIII_C, 1.2, 2026.80)]
    for d, t, c, lev, lx in ev:
        x = dyear(d)
        ax3.plot([x, x, lx], [0, lev - 0.35, lev - 0.05], color=c, lw=0.9, alpha=0.8)
        ax3.plot(x, 0, 'o', color=c, ms=5, mec='k', mew=0.4, zorder=5)
        ax3.text(lx, lev, f'{d}\n{t}', color=c, fontsize=7.3, ha='center', va='bottom', path_effects=[HALO])
    ax3.axhline(0, color=FGDIM, lw=0.8)
    ax3.set_xlim(2022.75, 2027.0); ax3.set_ylim(-0.4, 5.9)
    ax3.set_yticks([])
    ax3.set_xticks([2023, 2024, 2025, 2026, 2027])
    ax3.tick_params(labelsize=8.5)
    ax3.grid(axis='x', alpha=0.2)
    hand3 = [Line2D([], [], color=OIII_C, marker='o', label='SWONS / Wu+ team'),
             Line2D([], [], color=OK, marker='o', label='NSNS / Ziegenbalg'),
             Line2D([], [], color=WARN, marker='o', label='Fesen / Drechsler / Schaeffer')]
    ax3.legend(handles=hand3, loc='upper center', fontsize=7.6, framealpha=0.9)
    for xa, xb in ((2022.75, 2022.75), (2027.0, 2027.0)):
        con = ConnectionPatch(xyA=(xa, -0.8), coordsA=ax2.transData, xyB=(xb, 5.9), coordsB=ax3.transData,
                              color=FGDIM, lw=0.7, ls=':')
        fig.add_artist(con)
    ss.credit(fig)
    fig.savefig(FIG / 'w3_figs_f7_timeline.png', dpi=200)
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    main()
