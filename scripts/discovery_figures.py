"""Figures for the discovery-history note (notes/discovery.md).

1. figures/discovery_timeline.png        dated events by lane (2023-2026), from data/derived/discovery_timeline.json
2. figures/discovery_extents_map.png     RA/Dec map of all published geometries vs Monogem circles and neighbours
3. figures/discovery_quoted_extents.png  quoted lengths over time + [O III] surface-brightness estimates

Run after scripts/discovery_geometry.py:
    ~/.venvs/astro313/bin/python -I scripts/discovery_figures.py
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, FG, FGDIM, ACCENT, WARN, DANGER, OK, SUN,  # noqa: E402
                        OIII_C, HA_C, SII_C, NII_C, BG2)
from silk_geometry import CENTER, REGIONS, MONOGEM_PULSAR, GEMINGA, load_trace  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
TL = json.loads((ROOT / 'data/derived/discovery_timeline.json').read_text())
GEO = json.loads((ROOT / 'data/derived/discovery_extents.json').read_text())

LANES = [('swons', 'SWONS / LAMOST (Wu+)', OIII_C),
         ('ziegenbalg', 'Ziegenbalg / NSNS', SUN),
         ('fesen_nhz', 'Fesen / NHZ / Drechsler', HA_C),
         ('other_amateur', 'Other amateurs', OK),
         ('desi', 'DESI', NII_C),
         ('press', 'Press', FGDIM)]
SHORT = {  # compact labels for the timeline (keyed by date + lane); '' = marker only
    ('2023-11', 'swons'): 'SWONS candidate',
    ('2023-12', 'swons'): '',
    ('2024-01', 'ziegenbalg'): 'arc in full-array previews',
    ('2024-01-10', 'swons'): 'ε-180',
    ('2024-01-14', 'swons'): 'zh47r9 (1st image)',
    ('2024-02-01', 'swons'): 'LAMOST',
    ('2024-04', 'swons'): 'CNA magazine No.4',
    ('2024-05-21', 'swons'): '',
    ('2024-06-11', 'ziegenbalg'): 'CN post (14°)',
    ('2024-06-12', 'fesen_nhz'): '"part of Monogem"',
    ('2024-06-15', 'ziegenbalg'): '',
    ('2024-06-28', 'ziegenbalg'): 'NSNS page/K1a',
    ('2024-07-12', 'swons'): 'QHYCCD PR (15°)',
    ('2024-07-13', 'fesen_nhz'): '',
    ('2024-07-16', 'swons'): '',
    ('2024-10-01', 'ziegenbalg'): 'DR0.1',
    ('2025-02-05', 'other_amateur'): 'y32271 (other object)',
    ('2025-03-04', 'other_amateur'): 'visual (Bartels)',
    ('2025-03-18', 'other_amateur'): '',
    ('2025-05-02', 'other_amateur'): 'zwfplw (15°)',
    ('2025-05-22', 'ziegenbalg'): 'DR0.2',
    ('2025-06-11', 'ziegenbalg'): '17° (O200.0)',
    ('2025-08-28', 'ziegenbalg'): 'RNAAS 9, 227',
    ('2025-09-27', 'ziegenbalg'): 'Zenodo',
    ('2026-02-28', 'other_amateur'): '"AstroBin first?"',
    ('2026-07-30', 'swons'): 'AJ received',
    ('2026-08-27', 'fesen_nhz'): 'Fesen+ arXiv',
    ('2026-08-28', 'desi'): 'DESI sky-fibre paper',
    ('2026-10-05', 'swons'): '',
    ('2026-10-09', 'swons'): 'arXiv',
    ('2026-11-04', 'swons'): 'AJ (+credit)',
    ('2026-11-16', 'press'): 'Chinese press',
    ('2023..2024', 'desi'): 'Koposov? (hearsay)',
}
CHAR_DAYS = 6.6          # approx. width of one 7-pt Menlo character in days on the timeline axis
LEVELS = (0.20, -0.24, 0.40, -0.44)
REL_MARK = {'high': 'o', 'medium': 's', 'low': 'X'}


def parse_date(s):
    if '..' in s:
        a, b = s.split('..')
        return parse_date(a), parse_date(b)
    p = [int(x) for x in s.split('-')]
    if len(p) == 1:
        return date(p[0], 7, 1)
    if len(p) == 2:
        return date(p[0], p[1], 15)
    return date(*p)


def fig_timeline():
    fig, ax = new_sky_figure(13, 7.0, 'Silk Nebula · discovery timeline',
                             'dated public records 2023–2026 · marker: ● primary/timestamped  ■ self-reported/approx.  ✕ hearsay')
    lane_y = {k: i for i, (k, _, _) in enumerate(LANES[::-1])}
    t0, t1 = date(2023, 9, 1), date(2026, 12, 15)
    for k, name, col in LANES:
        y = lane_y[k]
        ax.axhline(y, color=col, lw=0.6, alpha=0.35)
    placed = {}  # lane -> list of (level, x0, x1)
    for ev in TL['events']:
        lane = ev['lane']
        if lane not in lane_y:
            continue
        d = parse_date(ev['date'])
        col = dict((k, c) for k, _, c in LANES)[lane]
        y = lane_y[lane]
        if isinstance(d, tuple):
            a, b = max(d[0], t0), min(d[1], t1)
            ax.plot([a, b], [y, y], color=col, lw=5, alpha=0.35, solid_capstyle='butt')
            d = a + (b - a) / 2
        elif d < t0:
            continue
        ax.scatter([d], [y], marker=REL_MARK[ev['reliability']], s=46, color=col, edgecolor='k', lw=0.5, zorder=5)
        txt = SHORT.get((ev['date'], lane), ev['event'][:24])
        if not txt:
            continue
        x0 = mdates.date2num(d) - 3
        x1 = x0 + CHAR_DAYS * len(txt) + 6
        used = placed.setdefault(lane, [])
        for lev in LEVELS:
            if all(not (lv == lev and x0 < b1 and a1 < x1) for lv, a1, b1 in used):
                break
        used.append((lev, x0, x1))
        ax.plot([d, d], [y, y + lev * 0.75], color=col, lw=0.5, alpha=0.7)
        label(ax, x0, y + lev, txt, color=FG, size=7.0, ha='left', va='center')
    ax.set_xlim(t0, t1)
    ax.set_ylim(-0.7, len(LANES) - 0.4)
    ax.set_yticks([lane_y[k] for k, _, _ in LANES])
    ax.set_yticklabels([n for _, n, _ in LANES], fontsize=8.5, fontweight='bold')
    for tl, (_, _, c) in zip(ax.get_yticklabels(), LANES):
        tl.set_color(c)
    ax.tick_params(axis='y', length=0)
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.grid(axis='x', alpha=0.4)
    ax.axvline(date(2026, 11, 16), color=WARN, lw=0.8, ls='-')
    label(ax, date(2026, 11, 16), len(LANES) - 0.55, 'today', color=WARN, size=7.5, ha='right')
    for s in ('left', 'right', 'top'):
        ax.spines[s].set_visible(False)
    return finish(fig, ROOT / 'figures/discovery_timeline.png')


def ellipse_pts(ra, dec, maj_arcmin, mnr_arcmin, pa_deg, n=240):
    """Sky ellipse sampled with directional offsets (major axis along PA, E of N)."""
    c = SkyCoord(ra * u.deg, dec * u.deg)
    a, b = maj_arcmin / 120.0, mnr_arcmin / 120.0
    th = np.linspace(0, 2 * np.pi, n)
    # radius of ellipse in direction phi relative to major axis
    r = a * b / np.sqrt((b * np.cos(th)) ** 2 + (a * np.sin(th)) ** 2)
    p = c.directional_offset_by((pa_deg + np.degrees(th)) * u.deg, r * u.deg)
    return p.ra.deg, p.dec.deg


def circle_pts(ra, dec, r_deg, n=360):
    c = SkyCoord(ra * u.deg, dec * u.deg)
    p = c.directional_offset_by(np.linspace(0, 360, n) * u.deg, r_deg * u.deg)
    return p.ra.deg, p.dec.deg


def fig_map():
    fig, ax = new_sky_figure(9.2, 11.4, 'Silk Nebula · published geometries',
                             'NSNS ellipses · measured [O III] spine · Monogem (Knies+2024) · Fesen+2026 SNRs · ◆ rim crossing')
    dec0 = 17.0
    # Monogem circles
    mono = GEO['monogem_crossings']
    styles = {'single (Knies+2024)': (FGDIM, '-', 'Monogem 24° (single)'),
              'two-SNR: Monogem (Knies+2024)': (FGDIM, ':', 'Monogem / G205.6+12.4 (two-SNR)'),
              'two-SNR: G205.6+12.4 (Knies+2024)': (FGDIM, ':', None)}
    for k, (col, ls, lab) in styles.items():
        m = mono[k]
        x, y = circle_pts(m['center'][0], m['center'][1], m['radius_deg'])
        ax.plot(x, y, color=col, ls=ls, lw=1.0, alpha=0.8, label=lab)
    label(ax, 110.3, 22.2, 'G205.6+12.4 (two-SNR)', color=FGDIM, size=7, ha='center')
    # neighbours
    for k, v in GEO['neighbours'].items():
        if v['kind'] == 'pulsar':
            ax.scatter(v['ra'], v['dec'], marker='*', s=80, color=SUN, edgecolor='k', lw=0.4, zorder=6)
            label(ax, v['ra'] - 0.4, v['dec'] + 0.3, k, color=SUN, size=7.5)
            continue
        r = max(v['size_deg']) / 2
        x, y = circle_pts(v['ra'], v['dec'], r)
        ax.plot(x, y, color=WARN, lw=0.9)
        label(ax, v['ra'] - r * 0.4, v['dec'] + r + 0.25, k, color=WARN, size=7.5, ha='center')
    # NSNS ellipses
    E = GEO['ellipses']
    for k, col, ls, lw, lab in [('K1a_2024', SUN, '--', 1.0, 'NSNS K1a 2024 (14°)'),
                                ('O200.0_2025', SUN, '-', 1.5, 'NSNS O200.0 2025 (17°)'),
                                ('O200.0a_Ha_SII', HA_C, '-', 1.2, 'Hα+[S II]-only part'),
                                ('O200.0b_faint_south', OIII_C, ':', 1.2, 'faint S extension'),
                                ('O205.0_east', NII_C, '-.', 1.0, 'E filaments (Monogem?)')]:
        e = E[k]
        x, y = ellipse_pts(e['ra'], e['dec'], e['maj'], e['mnr'], e['pa'])
        ax.plot(x, y, color=col, ls=ls, lw=lw, label=lab)
    e = E['O200.0_2025']
    ax.plot([e['north_end_radec'][0], e['south_end_radec'][0]], [e['north_end_radec'][1], e['south_end_radec'][1]],
            color=SUN, lw=0.6, alpha=0.7)
    # Monogem crossing point(s) on the 2025 axis
    src = GEO.get('measured_trace', {}).get('monogem_crossings') or mono  # prefer crossings on the measured spine
    for k in ('single (Knies+2024)', 'two-SNR: G205.6+12.4 (Knies+2024)'):
        for (xr, yd) in src[k]['crossings_radec']:
            ax.scatter(xr, yd, marker='D', s=26, color=FG, edgecolor='k', zorder=8)
            label(ax, xr + 0.35, yd, f'rim {yd:.1f}°', color=FG, size=6.5, ha='right', va='center')
    # zwfplw approximate footprint (long axis ~ N-S)
    z = GEO['astrobin_solutions']['zwfplw']
    L, W = z['field_deg']
    c = SkyCoord(z['ra'] * u.deg, z['dec'] * u.deg)
    pa = (z['orientation'] - 270.0) % 360  # image x-axis ~ along this PA (approximate)
    corners = []
    for sl, sw in [(1, 1), (1, -1), (-1, -1), (-1, 1), (1, 1)]:
        p = c.directional_offset_by(pa * u.deg, sl * L / 2 * u.deg)
        p = p.directional_offset_by((pa + 90) * u.deg, sw * W / 2 * u.deg)
        corners.append((p.ra.deg, p.dec.deg))
    corners = np.array(corners)
    ax.plot(corners[:, 0], corners[:, 1], color=OK, lw=0.9, ls='--', label='AstroBin zwfplw field (approx.)')
    # LAMOST plate and paper regions
    h = GEO['hip37797']
    x, y = circle_pts(h['ra'], h['dec'], 2.5)
    ax.plot(x, y, color=ACCENT, lw=1.0, label='LAMOST plate (HIP 37797, Ø5°)')
    for k, (ra, dec) in REGIONS.items():
        ax.scatter(ra, dec, s=14, color=ACCENT if k != 'P4' else FGDIM, edgecolor='k', lw=0.3, zorder=7)
        if k == 'P4':
            label(ax, ra - 0.25, dec, k, color=FGDIM, size=7.5, ha='left', va='center')
        else:
            label(ax, ra + 0.25, dec, k, color=ACCENT, size=7.5, ha='right', va='center')
    ax.scatter(*CENTER, marker='+', s=60, color=FG, zorder=7)
    label(ax, CENTER[0] - 0.25, CENTER[1] + 0.5, 'Wu+ centre', color=FG, size=7, ha='left')
    # measured trace (geometry teammate), if present
    tr = load_trace()
    if tr:
        if tr.get('preliminary', True):
            pts = np.array([(p['ra'], p['dec']) for p in tr['spine'] if p.get('sigma_deg', 1) <= 0.3])
            lab = f"trace {tr.get('version', '')} (prelim.)"
        else:
            pts = np.array([(p['ra'], p['dec']) for p in tr['spine']])
            lab = f"[O III] spine (silk_trace v{tr.get('version', '')})"
            w = tr.get('strands', {}).get('W', {}).get('points', [])
            if w:
                w = np.array([(p['ra'], p['dec']) for p in w])
                ax.plot(w[:, 0], w[:, 1], color=OIII_C, lw=1.0, alpha=0.8)
            hp = np.array(tr.get('halpha_patch') or [])
            if len(hp):
                ax.fill(hp[:, 0], hp[:, 1], color=HA_C, alpha=0.8, lw=0, label='Hα patch (measured)')
        if len(pts):
            ax.plot(pts[:, 0], pts[:, 1], color=OIII_C, lw=2.2, alpha=0.95, label=lab)
    label(ax, 112.3, 13.3, 'Sh2-274', color=FGDIM, size=7)
    ax.scatter(112.18, 13.35, s=8, color=FGDIM)
    ax.set_xlim(124.5, 95.5)
    ax.set_ylim(2.5, 31.5)
    ax.set_aspect(1 / np.cos(np.radians(dec0)))
    ax.set_xlabel('RA (deg, J2000)')
    ax.set_ylabel('Dec (deg, J2000)')
    ax.grid(alpha=0.35)
    ax.legend(loc='upper center', bbox_to_anchor=(0.42, -0.06), ncol=3, fontsize=7, framealpha=0.85)
    # galactic latitude guides
    for bb in (10, 15, 20, 25):
        g = SkyCoord(l=np.linspace(180, 220, 200) * u.deg, b=np.full(200, bb) * u.deg, frame='galactic').icrs
        ax.plot(g.ra.deg, g.dec.deg, color='#2f3a46', lw=0.6, ls='--')
        sel = (g.ra.deg > 96) & (g.ra.deg < 124) & (g.dec.deg > 3) & (g.dec.deg < 30.5)
        if sel.any():
            i = np.where(sel)[0][np.argmax(g.dec.deg[sel])]
            label(ax, g.ra.deg[i] - 0.2, g.dec.deg[i] - 0.3, f'b=+{bb}°', color='#7c8a98', size=6.5, ha='left')
    return finish(fig, ROOT / 'figures/discovery_extents_map.png')


def fig_quoted():
    fig, axes = new_sky_figure(12, 5.9, 'Silk Nebula · quoted size and brightness',
                               'angular length claimed by each source vs. publication date; [O III] surface brightness → Rayleigh',
                               ncols=2, gridspec_kw={'width_ratios': [1.25, 1]})
    ax = axes[0]
    # (date, length, width, colour, label, dy, ha)
    pts = [('2024-01-14', 10, 1.0, OIII_C, 'SWONS AstroBin (1.0°)', 0.5, 'left'),
           ('2024-06-11', 14, 2.5, SUN, 'Ziegenbalg CN/NSNS (2.5°)', -0.6, 'left'),
           ('2024-07-12', 15, 1.5, OIII_C, 'QHYCCD PR (1.5°)', 0.5, 'left'),
           ('2025-05-02', 15, None, OK, 'Jackson AstroBin', -0.6, 'center'),
           ('2025-06-11', 17, 2.5, SUN, 'NSNS O200.0 (2.5°)', 0.55, 'right'),
           ('2025-08-28', 17, 3.0, SUN, 'RNAAS 2025 (3.0°)', -0.65, 'left'),
           ('2026-08-27', 15, None, HA_C, 'Fesen+ 2026', 0.5, 'center'),
           ('2026-10-09', 10, 1.0, OIII_C, 'Wu+ 2026 (1.0°)', 0.5, 'right')]
    for d, L, W, col, lab, dy, ha in pts:
        dd = parse_date(d)
        ax.scatter(dd, L, s=60 if W is None else 40 + 25 * W, color=col, edgecolor='k', lw=0.5, zorder=5)
        dx = {'left': -10, 'right': 10, 'center': 0}[ha]
        label(ax, mdates.date2num(dd) + dx, L + dy, lab, color=FG, size=7.2, ha=ha, va='center')
    ax.set_ylim(8, 19)
    ax.set_xlim(date(2023, 11, 1), date(2026, 12, 31))
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.set_ylabel('quoted length (deg)   [marker size ∝ width]')
    ax.grid(alpha=0.35)
    ax2 = ax.twinx()
    ax2.set_ylim(*[2 * 160 * np.tan(np.radians(v / 2)) for v in ax.get_ylim()])
    ax2.set_ylabel('length at 160 pc (pc)', color=FGDIM)

    ax = axes[1]
    sb = GEO['surface_brightness']
    rows = [('Wu+ P1 (23.60)', sb['Wu+2026 P1'], OIII_C), ('Wu+ P2 (23.87)', sb['Wu+2026 P2'], OIII_C),
            ('Wu+ P3 (23.85)', sb['Wu+2026 P3'], OIII_C), ('QHYCCD PR (24.3)', sb['QHYCCD PR brightest'], OIII_C)]
    for i, (lab, v, col) in enumerate(rows):
        lo, mid, hi = v['R_dlam3nm'], v['R_dlam5nm'], v['R_dlam6nm']
        ax.plot([lo, hi], [i, i], color=col, lw=3, alpha=0.6, solid_capstyle='butt')
        ax.scatter([hi], [i], color=col, edgecolor='k', zorder=5, s=36)
        label(ax, hi + 0.3, i, f'{hi:.1f} R (6 nm) … {lo:.1f} R (3 nm)', color=FG, size=7, va='center')
    ax.scatter([2.04, 1.53], [4, 4], color=SUN, edgecolor='k', s=36, zorder=5)
    label(ax, 2.6, 4, 'NSNS K1a / K1a′ (est.)', color=FG, size=7, va='center')
    ax.axvspan(0.5, 5, color=HA_C, alpha=0.12)
    label(ax, 0.6, 5.0, 'Fesen+2024 [O III] limit 0.5–5 R', color=HA_C, size=7)
    ax.set_yticks(range(5))
    ax.set_yticklabels([r[0] for r in rows] + ['Ziegenbalg NSNS'], fontsize=7.5)
    ax.set_ylim(-0.6, 5.5)
    ax.set_xlim(0, 22)
    ax.set_xlabel('[O III] 5007 surface brightness (Rayleigh)')
    ax.grid(axis='x', alpha=0.35)
    return finish(fig, ROOT / 'figures/discovery_quoted_extents.png')


if __name__ == '__main__':
    for f in (fig_timeline, fig_map, fig_quoted):
        print(f())
