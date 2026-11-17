"""verify_precovery figures (silk_style, dpi 200).
  figures/verify_precovery_shiftscan.png  combined core amplitude vs RA shift (g, u, r, i)
  figures/verify_precovery_dec.png        per-run amplitude vs crossing Dec, compared with DESI / SDSS fibres / paper
  figures/verify_precovery_frames4458.png run 4458 frame-level residual sky maps (g, i) and cross-ridge profiles
  figures/verify_precovery_timeline.png   per-epoch detection significance vs date, all datasets
"""
import sys, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from datetime import datetime
from astropy.time import Time
from silk_style import *
from verify_precovery_lib import DER, FIG, load_strands, desi_nodes, fibre_nodes

S = json.loads((DER / 'verify_precovery_scan.json').read_text())
ST = json.loads((DER / 'verify_precovery_stats.json').read_text())
TB = json.loads((DER / 'verify_precovery_table.json').read_text())
SED = json.loads((DER / 'verify_precovery_frames_sed.json').read_text())


def fig_scan():
    fig, axes = new_sky_figure(12, 4.8, 'SDSS imaging sky · RA-shift scan of the ridge template',
                               '19 runs crossing Dec 12–20.5 (2002–2007); template and ±2.5° window shifted together',
                               ncols=2)
    ax = axes[0]
    sty = {'g_psp_gauss028_OIIIunits': (OIII_C, '-', 'g PSP  [O III]'), 'g_framessub_gauss028_OIIIunits': (OIII_C, '--', 'g FramesSub  [O III]'),
           'u_psp_gauss028_OIIunits': (SII_C, '-', 'u PSP  [O II]'), 'u_framessub_gauss028_OIIunits': (SII_C, '--', 'u FramesSub  [O II]')}
    for k, (c, ls, lab) in sty.items():
        sc = S[k]['scan']; x = [p['shift'] for p in sc]
        ax.plot(x, [p['A'] for p in sc], color=c, ls=ls, lw=1.6, marker='o', ms=2.5, label=lab)
        if ls == '-':
            ax.plot(x, [p['A_expected_from_ridge'] for p in sc], color=c, ls=':', lw=1.0, alpha=0.8)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.axvspan(-0.6, 0.6, color='#ffffff', alpha=0.05)
    ax.set_xlabel('template + window shift in RA (deg)'); ax.set_ylabel('combined peak amplitude (R)')
    ax.legend(loc='upper right', ncol=2)
    label(ax, 0.02, 0.04, 'dotted: response expected from the s = 0 ridge alone', size=7.5, transform=ax.transAxes, color=FGDIM)
    ax = axes[1]
    for k, c, lab in [('g_psp_gauss028_OIIIunits', OIII_C, 'g PSP'), ('g_framessub_gauss028_OIIIunits', '#9be3ff', 'g FramesSub'),
                      ('u_psp_gauss028_OIIunits', SII_C, 'u PSP'), ('u_framessub_gauss028_OIIunits', '#ffd9a0', 'u FramesSub'),
                      ('r_framessub_gauss028', HA_C, 'r FramesSub'), ('i_framessub_gauss028', NII_C, 'i FramesSub')]:
        sc = S[k]['scan']; sd = S[k]['offpeak_std']
        ax.plot([p['shift'] for p in sc], [p['A'] / sd for p in sc], color=c, lw=1.3, label=lab)
    ax.axhspan(-1, 1, color='#ffffff', alpha=0.06); ax.axhline(0, color=FGDIM, lw=0.6)
    ax.axhline(3, color=FGDIM, lw=0.6, ls=':')
    ax.set_xlabel('shift in RA (deg)'); ax.set_ylabel('amplitude / off-peak scatter (|s| ≥ 2°)')
    ax.legend(loc='upper right', ncol=2)
    finish(fig, FIG / 'verify_precovery_shiftscan.png')


def fig_dec():
    fig, ax = new_sky_figure(11, 5.6, 'SDSS imaging sky · ridge amplitude along the filament',
                             'per-run peak of unit E+W Gaussian ridges (FWHM 0.28°); errors from off-window nulls; DESI and fibre curves: wave-1 calibration')
    pr = ST['per_run']
    for k, c, mk, lab, dx in [('g_framessub_gauss028_OIIIunits', OIII_C, 'o', 'SDSS g FramesSub [O III]5007', 0.0),
                              ('g_psp_gauss028_OIIIunits', '#9be3ff', 's', 'SDSS g PSP [O III]5007', 0.08),
                              ('u_framessub_gauss028_OIIunits', SII_C, 'D', 'SDSS u FramesSub [O II]3727', -0.08)]:
        rows = pr[k]
        ax.errorbar([r['dec_cross'] + dx for r in rows], [r['A'] for r in rows], [r['e'] for r in rows], fmt=mk, ms=3.5,
                    color=c, ecolor=c, elinewidth=0.8, alpha=0.9, label=lab)
    dn = desi_nodes()
    for kk, c, lab in (('E', '#7ad97a', 'DESI E (wave 1)'), ('W', '#c6f0a8', 'DESI W (wave 1)')):
        v = dn[kk]; ax.plot(v[:, 0], v[:, 1], color=c, lw=1.4, marker='^', ms=4, label=lab)
    fn = fibre_nodes()
    for kk, c, lab in (('E', '#c792ea', 'fibres E (wave 1)'), ('W', '#e7c6ff', 'fibres W (wave 1)')):
        v = fn[kk]; ax.errorbar(v[:, 0], v[:, 1], v[:, 2], color=c, lw=1.2, marker='v', ms=4, label=lab)
    ax.axhspan(10, 13, color=DANGER, alpha=0.15)
    label(ax, 22.6, 11.2, 'Wu+2026 imaging: 10–13 R', color=DANGER, size=8.5)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.set_xlim(11.8, 27.2); ax.set_ylim(-6, 17.5)
    ax.axvline(20.5, color=FGDIM, lw=0.6, ls=':')
    label(ax, 20.6, -5.3, 'north: u, r, i excess without g (broadband confounder)', size=7.5, color=FGDIM)
    ax.set_xlabel('Dec of the ridge crossing (deg)'); ax.set_ylabel('peak intensity (R)')
    ax.legend(loc='upper left', ncol=4, fontsize=7.2, bbox_to_anchor=(0.0, 1.0))
    finish(fig, FIG / 'verify_precovery_dec.png')


def fig_frames():
    from scipy.spatial import cKDTree
    strands, _ = load_strands()
    J = json.loads((DER / 'verify_precovery_frames.json').read_text())
    fig, axes = new_sky_figure(13, 6.4, 'SDSS run 4458 (2004-02-16) · frame-level sky (image + global sky model)',
                               '25″ block medians, common mode / column offsets / camcol drifts removed, smoothed 4′; units: R of [O III]5007 (g) or i-band equivalent',
                               ncols=3, gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
    for ax, band, vmax in ((axes[0], 'g', 6.0), (axes[1], 'i', 120.0)):
        z = np.load(DER / f'verify_precovery_frames_resid_{band}.npz')
        ra, dec, r, m = z['ra'], z['dec'], z['resid'], z['mask']
        ok = m & np.isfinite(r)
        # smooth with a 4' Gaussian using a KD-tree in tangent coordinates
        x = (ra - 116.5) * np.cos(np.radians(17)); y = dec
        t = cKDTree(np.c_[x[ok], y[ok]])
        gx, gy = np.meshgrid(np.arange(x[ok].min(), x[ok].max(), 0.03), np.arange(y[ok].min(), y[ok].max(), 0.03))
        sig = 4 / 60
        vals = np.full(gx.shape, np.nan)
        rr = r[ok]
        for i in range(gx.shape[0]):
            q = t.query_ball_point(np.c_[gx[i], gy[i]], 2.5 * sig)
            for j, lst in enumerate(q):
                if len(lst) > 20:
                    d2 = (x[ok][lst] - gx[i, j]) ** 2 + (y[ok][lst] - gy[i, j]) ** 2
                    w = np.exp(-0.5 * d2 / sig ** 2)
                    vals[i, j] = np.sum(w * rr[lst]) / w.sum()
        im = ax.pcolormesh(gx / np.cos(np.radians(17)) + 116.5, gy, vals, cmap='magma' if band == 'g' else 'viridis',
                           vmin=-vmax / 2, vmax=vmax, shading='auto')
        for k, c in (('E', '#ffffff'), ('W', '#9be3ff')):
            s = strands[k]
            ax.plot(s.ra_nodes, s.dec_nodes, color=c, lw=0.8, ls='--', alpha=0.8)
        ax.set_xlim(gx.max() / np.cos(np.radians(17)) + 116.5, gx.min() / np.cos(np.radians(17)) + 116.5)
        ax.set_ylim(gy.min(), gy.max())
        ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
        cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
        cb.set_label('R' if band == 'g' else 'i-band, R5007-equivalent')
        sb = SED['bands'][band]; un = sb['unit_nmgy_per_R'] * 1e3
        label(ax, 0.03, 0.96, f"{band}: A = {sb['base']['A'] / un:.1f} R; re-fitted shift nulls {sb['null_mean'] / un:.1f} ± {sb['null_std'] / un:.1f}",
              transform=ax.transAxes, size=8, va='top')
    ax = axes[2]
    e = np.array(J['analysis']['g']['tot']['profile_edges']); c = 0.5 * (e[1:] + e[:-1])
    for comp, col, lab in (('tot', OIII_C, 'g total sky'), ('img', '#9be3ff', 'g sky-subtracted image'), ('skymodel', SUN, 'g global sky model')):
        p = J['analysis']['g'][comp]['profile']['14.5-19.5']
        ax.plot(c, p['p'], color=col, lw=1.5, label=lab, drawstyle='steps-mid')
    pi = J['analysis']['i']['tot']['profile']['14.5-19.5']
    ax.plot(c, np.array(pi['p']) * J['analysis']['i']['conv'] / J['analysis']['g']['conv'], color=NII_C, lw=1.0,
            label='i total (same nmgy scale)', drawstyle='steps-mid')
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.set_xlabel('distance from E spine (deg, + east)'); ax.set_ylabel('g-band R of [O III]5007')
    ax.legend(loc='upper left', fontsize=7.5)
    ax.set_xlim(-1.8, 1.8)
    finish(fig, FIG / 'verify_precovery_frames4458.png')


def fig_timeline():
    fig, ax = new_sky_figure(12, 5.4, 'Silk Nebula · dated detections before and after discovery',
                             'single-epoch significance (run, plate, plan or month); line-specific spectra vs broadband sky')

    def dt(s):
        return datetime.strptime(s[:10] if len(s) >= 10 else s + '-15', '%Y-%m-%d')
    im = TB['sdss_imaging']
    for k, c, mk, lab in [('g_psp_gauss028_OIIIunits', OIII_C, 'o', 'SDSS imaging sky g ([O III])'),
                          ('u_framessub_gauss028_OIIunits', SII_C, 'D', 'SDSS imaging sky u ([O II])')]:
        rows = im[k]
        so = [r for r in rows if r['dec_cross'] < 20.5]; no = [r for r in rows if r['dec_cross'] >= 20.5]
        ax.scatter([dt(r['date']) for r in so], [r['snr'] for r in so], s=16, marker=mk, color=c, label=lab + ', Dec < 20.5', zorder=3)
        ax.scatter([dt(r['date']) for r in no], [r['snr'] for r in no], s=18, marker=mk, facecolors='none', edgecolors=c,
                   label=lab + ', north (confounded)', zorder=3, alpha=0.6)
    fb = [r for r in TB['sdss_fibres'] if r['n_on'] >= 2]
    ax.scatter([dt(r['date']) for r in fb], [r['snr'] for r in fb], s=16, marker='v', color=NII_C, label='SDSS/BOSS sky fibres [O III]', zorder=3)
    la = [r for r in TB['lamost_lrs'] if r['snr'] is not None]
    ax.scatter([dt(r['date']) for r in la], [r['snr'] for r in la], s=14, marker='s', color=HA_C, label='LAMOST LRS stellar residuals', zorder=3)
    de = TB['desi']
    ax.scatter([dt(r['month']) for r in de], [min(r['snr'], 14.5) for r in de], s=22, marker='^', color=OK, label='DESI DR1 (clipped at 14.5)', zorder=3)
    ax.axhline(3, color=FGDIM, lw=0.7, ls=':')
    for d, t, c, ha in [('2022-12-12', 'NSNS DR0.2 data start', ACCENT, 'right'), ('2023-11-15', 'SWONS detection', DANGER, 'right'),
                        ('2024-01-14', 'AstroBin zh47r9', WARN, 'left')]:
        ax.axvline(dt(d), color=c, lw=0.9, alpha=0.8)
        label(ax, dt(d), 13.8, t, color=c, size=7.5, rotation=90, va='top', ha=ha)
    ax.set_ylim(-4, 15); ax.set_ylabel('detection significance (σ)')
    ax.set_xlim(datetime(2001, 6, 1), datetime(2025, 1, 1))
    ax.legend(loc='upper left', ncol=2, fontsize=7.2, bbox_to_anchor=(0.0, 0.99))
    label(ax, datetime(2001, 8, 1), -3.4, 'DSS plates 1950–1996: not detected (A < 14 × NSNS, 3σ)', size=7.5, color=FGDIM)
    finish(fig, FIG / 'verify_precovery_timeline.png')


if __name__ == '__main__':
    fig_scan(); fig_dec(); fig_frames(); fig_timeline()
    print('ok')
