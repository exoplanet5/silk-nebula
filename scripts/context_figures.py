"""Figures for the 'context' topic.

figures/context_map.png                  Galactic map (l 232-178, b -13..+35): ROSAT 1/4 keV background, Monogem
                                         outlines, SNRs, pulsars + LSR back-tracks, hot stars, PNe, Silk trace
figures/context_curvature_pulsars.png    RA/Dec: Silk spine, small-circle fits, centre-of-curvature region, pulsar
                                         back-tracks and MC birth loci, candidate centres from the literature
figures/context_hotstars.png             RA/Dec zoom: hot ionising candidates (Teff, Q), Stromgren radii, LSR PMs
Inputs: data/derived/context_*.csv/json, data/derived/silk_trace.json, data/raw/context_skyview/*.fits
"""
import sys, json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *
from silk_style import *

D = ROOT / 'data/derived'
TR = load_trace()
SP = TR['spine']
SRA = np.array([q['ra'] for q in SP]); SDE = np.array([q['dec'] for q in SP])
WRA = np.array([q['ra'] for q in TR['strands']['W']['points']]); WDE = np.array([q['dec'] for q in TR['strands']['W']['points']])
HAP = np.array(TR['halpha_patch'])
PSR_C = '#ffd54f'; HOT_C = '#c792ea'; PN_C = OK; SNR_C = '#ff8a65'; MON_C = '#e8eaed'


def circle_radec(ra0, dec0, rho, n=361):
    """Small circle of angular radius rho (deg) around (ra0, dec0) -> (ra, dec) arrays."""
    c = unitvec(ra0, dec0)
    a = np.cross(c, [0, 0, 1.0]); a /= np.linalg.norm(a) if np.linalg.norm(a) > 0 else 1
    b = np.cross(c, a)
    t = np.linspace(0, 2 * np.pi, n)
    r = np.radians(rho)
    v = np.cos(r) * c[None] + np.sin(r) * (np.cos(t)[:, None] * a[None] + np.sin(t)[:, None] * b[None])
    return radec_of(v)


def to_lb(ra, dec):
    l, b = gal(np.asarray(ra), np.asarray(dec))
    return l, b


def plot_lb(ax, ra, dec, *a, **k):
    l, b = to_lb(ra, dec)
    # break lines that wrap
    return ax.plot(l, b, *a, **k)


def fig_map():
    h1 = fits.open(ROOT / 'data/raw/context_skyview/rass_bkg1_l205b5_70x60_car.fits')[0]
    h2 = fits.open(ROOT / 'data/raw/context_skyview/rass_bkg2_l205b5_70x60_car.fits')[0]
    img = gaussian_filter(np.nan_to_num(h1.data) + np.nan_to_num(h2.data), 0.8 / 2.355 / 0.2)
    w = WCS(h1.header)
    ny, nx = img.shape
    H = h1.header   # CAR with CRVAL1=0: compute l linearly to avoid the +-180 deg native-longitude wrap
    lc = (H['CRVAL1'] + H['CDELT1'] * (np.arange(nx) + 1 - H['CRPIX1'])) % 360
    bc = H['CRVAL2'] + H['CDELT2'] * (np.arange(ny) + 1 - H['CRPIX2'])
    fig, ax = new_sky_figure(12.5, 10.5, 'Silk Nebula · astrophysical neighbourhood',
                             'ROSAT 1/4 keV (R1+R2) background · Monogem outlines · pulsars with LSR back-tracks (0–1 Myr) · hot stars · SNRs')
    vmin, vmax = np.percentile(img, [3, 99.5])
    ax.imshow(img, origin='lower', extent=[lc[0], lc[-1], bc[0], bc[-1]], cmap='bone', vmin=vmin, vmax=vmax,
              alpha=0.85, aspect='auto', interpolation='bilinear')
    # Monogem outlines
    for (ra0, dec0, rho), ls, c, lab in [((107.08, 13.50, 11.97), '--', MON_C, 'Monogem (Knies+24, 1 SNR)'),
                                         ((104.15, 15.30, 11.00), ':', MON_C, 'Monogem (Knies+24, 2 SNRs)'),
                                         ((110.33, 12.94, 8.82), ':', OK, 'G205.6+12.4 (Knies+24)')]:
        r, d = circle_radec(ra0, dec0, rho)
        plot_lb(ax, r, d, ls=ls, color=c, lw=1.3, alpha=0.9)
    l0, b0 = to_lb(107.08, 13.50)
    label(ax, 207.5, -3.9, 'Monogem Ring', color=MON_C, size=9, ha='center')
    lg, bg = to_lb(110.33, 12.94)
    label(ax, lg[0] + 4.5, bg[0] - 7.0, 'G205.6+12.4', color=OK, size=8)
    # SNRs and candidates
    snr = pd.read_csv(D / 'context_snrs.csv')
    for _, s in snr.iterrows():
        if 'Monogem' in s['name'] or 'G205.6' in s['name'] or 'Ziegenbalg' in s['name'] or not np.isfinite(s.ra):
            continue
        if s['name'].startswith(('GS2', 'R20', 'W06')):          # HI shells / filaments: not SNR outlines
            continue
        rad = max(s.radius_deg, 0.25) if np.isfinite(s.radius_deg) else 0.5
        r, d = circle_radec(s.ra, s.dec, rad, 120)
        plot_lb(ax, r, d, color=SNR_C, lw=1.0, alpha=0.9)
    for nm, dl, db in [('G190.4+12.5', 0, 3.8), ('IC 443', -1.0, -2.0), ('Monoceros', 2.6, -2.6),
                       ('G206.7+5.9', -2.2, 1.4), ('S147', 1.5, -2.0), ('G203.1+6.6', 0.3, 1.6)]:
        m = snr[snr['name'].str.contains(nm.split()[0] if nm != 'IC 443' else 'IC443')]
        if len(m):
            label(ax, m.l.values[0] + dl, m.b.values[0] + db, nm, color=SNR_C, size=7.5, ha='center')
    # Monogem-related optical filaments R20 (Reich+20) and W06 (Weinberger+06), as small markers
    for nm, l_, b_ in [('R20', 203.7, 11.5), ('W06', 212.54, 10.62)]:
        ax.plot(l_, b_, 'x', color=HA_C, ms=6, mew=1.4, zorder=7)
        label(ax, l_ - 0.6, b_ - 0.4, nm, color=HA_C, size=7.5, va='top')
    # Origem Loop (Gao & Han 2013; SIMBAD centre 94.134,+16.072, diameter 400')
    r, d = circle_radec(94.1342, 16.0724, 400 / 120.0)
    plot_lb(ax, r, d, color='#90a4ae', lw=1.0, ls='-.')
    lo, bo = to_lb(94.1342, 16.0724)
    label(ax, lo[0] - 1.0, bo[0] - 4.2, 'Origem Loop', color='#90a4ae', size=7.5, ha='left')
    # lambda Ori, Barnard's Loop, Gem OB1 MC (SIMBAD positions)
    for nm, ra, dec in [('λ Ori', 83.775, 9.933), ('Gem OB1 MC', 94.75, 23.0)]:   # Barnard's Loop (b~-20) is off the map
        l_, b_ = to_lb(ra, dec)
        ax.plot(l_, b_, 'x', color='#90a4ae', ms=6, mew=1.2)
        label(ax, l_[0] - 0.8, b_[0] + 0.8, nm, color='#90a4ae', size=7.5)
    # fitted curvature circle and centre
    cv = TR['curvature']
    r, d = circle_radec(cv['center_ra'], cv['center_dec'], cv['radius_deg'])
    plot_lb(ax, r, d, color=OIII_C, lw=0.9, ls=(0, (6, 3, 1, 3)), alpha=0.8)
    lc_, bc_ = to_lb(cv['center_ra'], cv['center_dec'])
    ax.plot(lc_, bc_, '+', color=OIII_C, ms=13, mew=2)
    label(ax, lc_[0] - 1.2, bc_[0] - 0.4, 'centre of\ncurvature', color=OIII_C, size=8, ha='left', va='top')
    # hot stars
    hs = pd.read_csv(D / 'context_hotstars.csv')
    hot = hs[(hs.dist_pc < 1000) & (hs.teff >= 40000)]
    ax.scatter(hot.l, hot.b, s=14, marker='D', facecolor='none', edgecolor=HOT_C, lw=0.9, zorder=5)
    for nm, lab_, dl, db in [('UCAC4 524-042847', 'WD J0736+1444\n81 kK, 151 pc', -0.6, -2.2),
                             ('GALEX J074320.3+160752', 'PN G203.8+18.5 (DO, 1.2 kpc)', -0.9, -1.7)]:
        m = hs[hs.name == nm]
        if len(m):
            ax.plot(m.l, m.b, 'D', color=HOT_C, ms=5, zorder=6)
            label(ax, m.l.values[0] + dl, m.b.values[0] + db, lab_, color=HOT_C, size=7.5, ha='right' if dl < 0 else 'left')
    # PNe (SIMBAD)
    s = pd.read_csv(ROOT / 'data/raw/context_catalogs/simbad_hotstars_cone25.csv')
    pn = s[s.otype == 'PN']
    lp, bp = to_lb(pn.ra.values, pn.dec.values)
    ax.scatter(lp, bp, s=16, marker='o', facecolor='none', edgecolor=PN_C, lw=0.9, zorder=5)
    # pulsars < 1 kpc
    ps = pd.read_csv(D / 'context_pulsars_30deg.csv')
    near = ps[ps.dist_best_kpc < 1.0]
    ax.scatter(near.l, near.b, s=22, marker='^', color=PSR_C, edgecolor='k', lw=0.4, zorder=7)
    tr = pd.read_csv(D / 'context_pulsar_tracks.csv')
    mc = pd.read_csv(D / 'context_pulsar_birth_mc.csv')
    for nm, c in [('B0656+14', PSR_C), ('Geminga', '#ff7043')]:
        t = tr[tr['name'] == nm]
        ax.plot(t.l, t.b, '-', color=c, lw=1.6, zorder=6)
        tk = t[np.isclose((t.t_myr * 10) % 1, 0) | np.isclose((t.t_myr * 10) % 1, 1)]
        ax.plot(tk.l, tk.b, 'o', color=c, ms=2.5, zorder=6)
        m = mc[mc['name'] == nm]
        ax.scatter(m.l, m.b, s=3, color=c, alpha=0.25, lw=0, zorder=4)
    for nm, lab_, dl, db in [('J0659+1414', 'PSR B0656+14', 0.6, 1.0), ('J0633+1746', 'Geminga', 0.8, 0.9)]:
        m = ps[ps.psrj == nm]
        ax.plot(m.l, m.b, '*', color=PSR_C, ms=12, mec='k', mew=0.6, zorder=8)
        label(ax, m.l.values[0] + dl, m.b.values[0] + db, lab_, color=PSR_C, size=8.5, ha='right' if dl < 0 else 'left')
    # Silk trace
    ls_, bs_ = to_lb(SRA, SDE)
    ax.plot(ls_, bs_, '-', color=OIII_C, lw=3.0, solid_capstyle='round', zorder=9)
    lw_, bw_ = to_lb(WRA, WDE)
    ax.plot(lw_, bw_, '-', color=OIII_C, lw=1.4, zorder=9)
    lh, bh = to_lb(HAP[:, 0], HAP[:, 1])
    ax.add_patch(Polygon(np.c_[lh, bh], closed=True, color=HA_C, alpha=0.9, zorder=9))
    for k, (ra, dec) in REGIONS.items():
        l_, b_ = to_lb(ra, dec)
        ax.plot(l_, b_, 's' if k == 'P4' else 'o', ms=3.5, color='w', mfc='none' if k == 'P4' else 'w', mew=0.8, zorder=10)
    label(ax, 211.5, 22.6, 'Silk [O III]', color=OIII_C, size=10, fontweight='bold')
    ax.set_xlim(232, 178); ax.set_ylim(-13, 35)
    ax.set_xlabel('Galactic longitude l (deg)'); ax.set_ylabel('Galactic latitude b (deg)')
    ax.set_aspect('equal'); ax.grid(alpha=0.18)
    handles = [Line2D([], [], color=OIII_C, lw=3, label='Silk spine'),
               Line2D([], [], color=HA_C, lw=6, label='Hα patch'),
               Line2D([], [], color=MON_C, ls='--', label='Monogem X-ray'),
               Line2D([], [], color=OIII_C, ls=(0, (6, 3, 1, 3)), label='curvature fit'),
               Line2D([], [], color=SNR_C, label='SNR / cand.'),
               Line2D([], [], color=PSR_C, marker='^', ls='', label='pulsar <1 kpc'),
               Line2D([], [], color=PSR_C, lw=1.6, marker='o', ms=2.5, label='back-track 0–1 Myr'),
               Line2D([], [], color=HOT_C, marker='D', mfc='none', ls='', label='hot star ≥40 kK'),
               Line2D([], [], color=PN_C, marker='o', mfc='none', ls='', label='PN')]
    ax.legend(handles=handles, ncol=5, loc='lower left', fontsize=7.5, framealpha=0.85, handletextpad=0.6, columnspacing=1.2)
    finish(fig, ROOT / 'figures/context_map.png')


def fig_curvature():
    reg = np.load(D / 'context_curvature_regionA.npz')
    cj = json.loads((D / 'context_curvature.json').read_text())
    fig, ax = new_sky_figure(11, 7.3, 'Centre of curvature vs pulsar birth sites',
                             'small-circle fits to the Silk [O III] strands · pulsar LSR back-tracks (MC: v_r, distance, PM)')
    fits_ = TR['curvature']['fits']
    cols = {'E_full': OIII_C, 'E_north_ge17p75': OK, 'E_south_le17p5': WARN, 'W': '#c792ea', 'E_full_plus_W': FGDIM}
    for k, c in cols.items():
        f = fits_[k]
        r, d = circle_radec(f['center_ra'], f['center_dec'], f['radius_deg'])
        m = (r > 40) & (r < 130)
        ax.plot(np.where(m, r, np.nan), d, color=c, lw=0.8, alpha=0.6)
        ax.plot(f['center_ra'], f['center_dec'], '+', color=c, ms=10, mew=1.8)
        b = f['boot_16_84']
        ax.add_patch(Rectangle((b['center_ra'][0], b['center_dec'][0]), b['center_ra'][1] - b['center_ra'][0],
                               b['center_dec'][1] - b['center_dec'][0], fill=False, ec=c, lw=1.0, ls='--'))
        off = {'E_full': (-0.2, 1.1), 'E_north_ge17p75': (-3.6, 0.2), 'E_south_le17p5': (-0.6, 0.6), 'W': (-0.6, -1.4),
               'E_full_plus_W': (-0.6, -1.7)}[k]
        label(ax, f['center_ra'] + off[0], f['center_dec'] + off[1],
              f"{k.replace('_ge17p75','').replace('_le17p5','').replace('_plus_','+')} R={f['radius_deg']:.0f}°",
              color=c, size=7.5, ha='right' if off[0] > 0 else 'left')
    GR, GD = np.meshgrid(reg['ra'], reg['dec'])
    ax.contour(GR, GD, reg['dchi2'], levels=[2.30, 6.18], colors=[OIII_C, OIII_C], linewidths=[1.4, 0.8])
    mc = pd.read_csv(D / 'context_pulsar_birth_mc.csv')
    tr = pd.read_csv(D / 'context_pulsar_tracks.csv')
    for nm, c in [('B0656+14', PSR_C), ('Geminga', '#ff7043'), ('Geminga_ATNF', '#ffab91')]:
        m = mc[mc['name'] == nm]
        ax.scatter(m.ra, m.dec, s=3, color=c, alpha=0.3, lw=0)
        t = tr[tr['name'] == nm]
        ax.plot(t.ra, t.dec, '-', color=c, lw=1.5)
        tk = t[np.isclose((t.t_myr * 10) % 1, 0) | np.isclose((t.t_myr * 10) % 1, 1)]
        ax.plot(tk.ra, tk.dec, 'o', color=c, ms=2.5)
    for nm, ra, dec, c, dx, dy in [('PSR B0656+14', 104.9506, 14.2393, PSR_C, 0.5, -1.5), ('Geminga', 98.4756, 17.7703, '#ff7043', 0.7, 0.3)]:
        ax.plot(ra, dec, '*', color=c, ms=13, mec='k', mew=0.6, zorder=8)
        label(ax, ra + dx, dec + dy, nm, color=c, size=8.5, ha='right' if nm == 'Geminga' else 'left')
    for nm, ra, dec, mk in [('Monogem (K24 1-SNR)', 107.08, 13.50, 'x'), ('Monogem (K24 2-SNR)', 104.15, 15.30, 'x'),
                            ('G205.6+12.4', 110.33, 12.94, 'x'), ('Origem Loop', 94.1342, 16.0724, 'P'),
                            ('NSNS Geminga-SNR', 90.1412, 11.8861, 'v'), ('NSNS C0', 73.1139, 6.5749, 'v'),
                            ('Geminga birth (Pellizza+)', 84.7088, 7.6303, 's')]:
        ax.plot(ra, dec, mk, color=MON_C, ms=6, mew=1.2)
        dy = {'Monogem (K24 2-SNR)': 0.8, 'Origem Loop': -1.6}.get(nm, -1.3)
        label(ax, ra - 0.5, dec + dy, nm, color=FGDIM, size=7)
    ax.plot(SRA, SDE, '-', color=OIII_C, lw=3)
    ax.plot(WRA, WDE, '-', color=OIII_C, lw=1.2)
    ax.add_patch(Polygon(HAP, closed=True, color=HA_C, alpha=0.9))
    label(ax, SRA[40] + 0.6, SDE[40], 'Silk', color=OIII_C, size=10, fontweight='bold')
    ax.set_xlim(125, 65); ax.set_ylim(-2, 34)
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)'); ax.grid(alpha=0.2)
    ax.set_aspect(1 / np.cos(np.radians(16)))
    handles = [Line2D([], [], color=OIII_C, lw=1.4, label='Δχ²=2.3, 6.2 (E strand)'),
               Line2D([], [], color=OIII_C, ls='--', label='bootstrap 16–84% boxes'),
               Line2D([], [], color=PSR_C, lw=1.5, marker='o', ms=2.5, label='B0656+14 track (100 kyr ticks)'),
               Line2D([], [], color='#ff7043', lw=1.5, marker='o', ms=2.5, label='Geminga track'),
               Line2D([], [], color=FGDIM, marker='.', ls='', label='MC loci at τc')]
    ax.legend(handles=handles, ncol=3, loc='lower right', fontsize=7.5)
    finish(fig, ROOT / 'figures/context_curvature_pulsars.png')


def fig_hotstars():
    hs = pd.read_csv(D / 'context_hotstars.csv')
    fig, ax = new_sky_figure(10, 8.6, 'Hot ionising candidates near the Silk',
                             'Teff ≥ 40 kK, d < 1 kpc (filled) · circle = Strömgren radius at n = 0.1 cm⁻³ · arrows = LSR proper motion × 50 kyr')
    m = (hs.teff >= 40000) & (hs.dist_pc < 1000) & (hs.ra > 100) & (hs.ra < 134) & (hs.dec > 2) & (hs.dec < 31)
    h = hs[m]
    sc = ax.scatter(h.ra, h.dec, c=h.teff / 1e3, s=12 + 8 * np.log10(h.Q_H / 1e43).clip(0, 6) ** 2, cmap='plasma',
                    vmin=40, vmax=110, edgecolor='w', lw=0.4, zorder=6)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.01); cb.set_label('Teff (kK)')
    for _, r in h.iterrows():
        rho = np.degrees(r['Rs_pc_n0.1'] / r.dist_pc)
        if rho < 20 and r.sep_trace_deg < 6:
            cr, cd = circle_radec(r.ra, r.dec, rho, 180)
            ax.plot(cr, cd, color=HOT_C, lw=0.6, alpha=0.6)
        if np.isfinite(r.pm_lsr_masyr):
            ang = np.radians(r.pa_pm_lsr); L = r.pm_lsr_masyr * 50e3 / 3.6e6
            ax.annotate('', xy=(r.ra + L * np.sin(ang) / np.cos(np.radians(r.dec)), r.dec + L * np.cos(ang)), xytext=(r.ra, r.dec),
                        arrowprops=dict(arrowstyle='->', color=FGDIM, lw=0.8))
    oth = hs[(hs.dist_pc >= 1000) & (hs.teff >= 30000) & (hs.ra > 100) & (hs.ra < 134) & (hs.dec > 2) & (hs.dec < 31)]
    ax.scatter(oth.ra, oth.dec, s=14, facecolor='none', edgecolor=FGDIM, lw=0.6, zorder=5)
    # (label x, label y, ha) in data coords; RA axis is inverted so ha='left' extends toward lower RA
    for nm, lab_, x, y, ha in [('UCAC4 524-042847', 'WD J0736+1444 (DOZ, 81 kK, 151 pc)', 113.8, 14.15, 'left'),
                               ('GALEX J074507.5+144939', 'GALEX J0745+1449 (DO, 594 pc)', 116.9, 14.2, 'right'),
                               ('GALEX J074320.3+160752', 'PN G203.8+18.5 (DO, 1.2 kpc)', 115.0, 16.85, 'left'),
                               ('HD  60848', 'HD 60848 (O8Vpe, 1.9 kpc)', 113.6, 17.65, 'left'),
                               ('PN A66   21', 'Medusa (Sh2-274)', 112.0, 12.55, 'left'),
                               ('NGC  2392', 'NGC 2392', 112.0, 21.4, 'left')]:
        r = hs[hs.name == nm]
        if len(r):
            ax.plot([r.ra.values[0], x], [r.dec.values[0], y + 0.15], '-', color='#e1bee7', lw=0.5, alpha=0.7)
            label(ax, x, y, lab_, color='#e1bee7', size=7.5, ha=ha)
    ax.plot(SRA, SDE, '-', color=OIII_C, lw=3, zorder=4)
    ax.plot(WRA, WDE, '-', color=OIII_C, lw=1.2, zorder=4)
    ax.add_patch(Polygon(HAP, closed=True, color=HA_C, alpha=0.9, zorder=4))
    for k, (ra, dec) in REGIONS.items():
        ax.plot(ra, dec, 's' if k == 'P4' else 'o', color='w', ms=4, mfc='none' if k == 'P4' else 'w', zorder=7)
        label(ax, ra - 0.25, dec, k, color='w', size=7.5, va='center')
    ax.set_xlim(134, 100); ax.set_ylim(2, 31)
    ax.set_aspect(1 / np.cos(np.radians(16)))
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)'); ax.grid(alpha=0.2)
    finish(fig, ROOT / 'figures/context_hotstars.png')


if __name__ == '__main__':
    which = sys.argv[1:] or ['map', 'curv', 'hot']
    if 'map' in which: fig_map()
    if 'curv' in which: fig_curvature()
    if 'hot' in which: fig_hotstars()
