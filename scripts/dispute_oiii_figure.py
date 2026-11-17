"""dispute_oiii: Figure 3 of the Comment on Wu et al. 2026 (AJ 172, 298; arXiv:2610.13329v1).

(a) [N I] 5198/5200 airglow in the paper's own Fig. 3a spectrum (PDF vector polyline, re-extracted by
    scripts/verify_velocity_fig3_vectors.py), with the positions at which the airglow must appear if the plotted
    wavelength scale were air or vacuum, topocentric or heliocentric.
(b) [O III] 5007 in the same spectrum (and the Fig. 3b P1-P3 spectra, normalised), with air/vacuum rest and the
    plotted-scale positions of a line with true V_helio = +13.3 and +25 km/s.
(c) Velocity ladder: the paper's published Table 4, the paper's Fig. 3 on its plotted scale, corrected values, DESI DR1, SDSS,
    WHAM and the paper's red lines.

Inputs: data/derived/verify_velocity_fig3_spectra.npz, verify_velocity_fig3_measure.json,
data/derived/dispute_oiii_values.json (from scripts/dispute_oiii_values.py).
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_oiii_figure.py
Output: figures/dispute_oiii_fig3.png (dpi 200)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *  # noqa: E402,F403
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
OUT = ROOT / 'figures/dispute_oiii_fig3.png'
C = 299792.458

V = json.loads((D / 'dispute_oiii_values.json').read_text())
FM = json.loads((D / 'verify_velocity_fig3_measure.json').read_text())
SP = np.load(D / 'verify_velocity_fig3_spectra.npz')

# frame-marker styles (shared by panels a and b)
FR = {
    'air_topocentric': dict(color=OK, ls='-', lw=1.3, label='air, topocentric'),
    'air_heliocentric': dict(color=WARN, ls='--', lw=1.3, label='air, heliocentric'),
    'vacuum_topocentric': dict(color=ACCENT, ls=':', lw=1.5, label='vac, topocentric'),
    'vacuum_heliocentric_LAMOST_native': dict(color=NII_C, ls='-.', lw=1.3, label='vac, helio (LAMOST)'),
}
DATA_C = '#dfe4ea'


def gauss(x, a, mu, s):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def panel_a(ax):
    x, y = SP['a_wave'], SP['a_flux']
    m = (x > 5196.6) & (x < 5202.3)
    r = V['NI_doublet_refit']
    c0, c1 = r['c0'], r['c1']
    yn = (y - (c0 + c1 * (x - 5199))) / r['A1']
    ax.step(x[m], yn[m], where='mid', color=DATA_C, lw=1.1, zorder=3)
    xx = np.linspace(5196.6, 5202.3, 1200)
    mod = (gauss(xx, r['A1'], r['mu1'], r['sig']) + gauss(xx, r['A2'], r['mu2'], r['sig'])) / r['A1']
    ax.plot(xx, mod, color=OIII_C, lw=1.0, alpha=0.75, zorder=4)
    fp = V['airglow_frames']['frame_positions_A']
    for k, st in FR.items():
        for lam in fp[k]:
            ax.axvline(lam, ymin=0.0, ymax=0.60, color=st['color'], ls=st['ls'], lw=st['lw'], alpha=0.95, zorder=2)
    ax.set_xlim(5196.6, 5202.3)
    ax.set_ylim(-0.3, 2.0)
    ax.set_xlabel('wavelength on the paper’s plotted scale (Å)')
    ax.set_ylabel('flux − continuum (norm.)')
    ax.set_title('(a) [N I] airglow in paper Fig. 3a', loc='left', color=FG)
    ft = V['frame_test']
    label(ax, 5196.72, 1.86, '[N I] sits at air, topocentric rest', color=OK, size=8.4, va='center')
    label(ax, 5196.72, 1.69, f"offset {ft['air_topocentric']['measured']:+.1f} ± {ft['air_topocentric']['err']:.1f} km/s "
          '(5 sky lines)', color=FG, size=8.4, va='center')
    label(ax, 5196.72, 1.52, f"air-helio excluded at {ft['air_heliocentric']['n_sigma']:.0f}σ, "
          f"vacuum at {ft['vacuum_heliocentric_LAMOST_native']['n_sigma']:.0f}σ", color=FGDIM, size=8.4, va='center')
    h = [Line2D([], [], **{k: v for k, v in st.items() if k != 'label'}) for st in FR.values()]
    ax.legend(h, [st['label'] for st in FR.values()], loc='upper right', ncol=1, fontsize=8.2,
              handlelength=2.6, framealpha=0.92)
    ax.grid(alpha=0.25)


def panel_b(ax):
    o3a = V['oiii_rest']['air']
    o3v = V['oiii_rest']['vac']
    x, y = SP['a_wave'], SP['a_flux']
    m = (x > 5005.6) & (x < 5008.75)
    o = V['oiii_fig3a']
    ax.step(x[m], (y[m] - o['cont']) / o['A'], where='mid', color=DATA_C, lw=1.2, zorder=4, label='Fig. 3a')
    xx = np.linspace(5005.6, 5008.75, 800)
    ax.plot(xx, gauss(xx, 1, o['mu'], o['sig']), color=OIII_C, lw=1.1, alpha=0.8, zorder=5)
    cols = {'P1': '#c792ea', 'P2': '#ff8a80', 'P3': '#ffd54f'}
    for r, c in cols.items():
        b = FM['panel_b'][f'b_{r}']
        xb, yb = SP[f'b_{r}_wave'], SP[f'b_{r}_flux']
        mb = (xb > 5005.6) & (xb < 5008.75)
        ax.step(xb[mb], (yb[mb] - b['cont']) / b['A'], where='mid', color=c, lw=0.8, alpha=0.55, zorder=3,
                label=f'Fig. 3b {r}')
    # rest wavelengths (same styles as panel a)
    ax.axvline(o3a, ymax=0.62, **{k: v for k, v in FR['air_topocentric'].items() if k != 'label'})
    ax.axvline(o3v, ymax=0.62, **{k: v for k, v in FR['vacuum_topocentric'].items() if k != 'label'})
    label(ax, o3a - 0.03, -0.27, 'air rest', color=OK, size=8.4, ha='right')
    label(ax, o3v - 0.03, -0.27, 'vacuum rest', color=ACCENT, size=8.4, ha='right')
    # where a line of given true V_helio appears on the plotted scale
    pp = V['plotted_position_for_true_vhelio']
    l13, l25 = pp['13.3']['lam'], pp['25.0']['lam']
    ax.axvline(l13, ymax=0.76, color=OK, lw=2.2, alpha=0.55)
    ax.axvline(l25, ymax=0.76, color=DANGER, lw=1.6, ls='--')
    label(ax, l13 - 0.03, 1.27, 'true V$_{helio}$ +13.3', color=OK, size=8.4, ha='right')
    label(ax, l25 + 0.03, 1.27, 'true V$_{helio}$ +25', color=DANGER, size=8.4, ha='left')
    ax.set_xlim(5005.6, 5008.75)
    ax.set_ylim(-0.35, 1.75)
    ax.set_xlabel('wavelength on the paper’s plotted scale (Å)')
    ax.set_ylabel('flux − continuum (norm.)')
    ax.set_title('(b) [O III] 5007 on the same scale', loc='left', color=FG)
    # top axis: velocity read off the plotted scale (the paper's Table 4 reading)
    sec = ax.secondary_xaxis('top', functions=(lambda l: C * (l / o3a - 1), lambda v: o3a * (1 + v / C)))
    sec.set_xlabel('V read from plotted scale (km/s, air rest)', color=FGDIM, fontsize=8.8)
    sec.tick_params(labelsize=8.5)
    b = V['budget_fig3a']
    label(ax, 5005.68, 1.61, f"fit {o['mu']:.3f} Å → {o['V_plot']:+.1f} km/s as read (paper Table 4: +23.8…+26.8)",
          color=FG, size=8.4, va='center')
    label(ax, 5005.68, 1.46, f"vs own [N I]: V$_{{helio}}$ {b['V_helio']:+.1f} ± {b['err']:.1f} km/s",
          color=OK, size=8.4, va='center')
    ax.legend(loc="upper left", fontsize=8, ncol=1, framealpha=0.92, bbox_to_anchor=(0.0, 0.67))
    ax.grid(alpha=0.25)


def panel_c(ax):
    R = V['regions']
    I = V['independent']
    RL = V['red_lines']
    S = dict(pub=dict(color=DANGER, marker='o', ms=6, mfc=DANGER),
             plot=dict(color=WARN, marker='s', ms=5.5, mfc='none'),
             corr=dict(color=OK, marker='D', ms=5.5, mfc=OK),
             desi=dict(color=ACCENT, marker='^', ms=6, mfc=ACCENT),
             sdss=dict(color='#80deea', marker='v', ms=6, mfc='none'),
             wim=dict(color=NII_C, marker='o', ms=5, mfc=NII_C))
    rows = ['P1', 'P2', 'P3', 'Fig. 3a', 'DESI DR1', 'SDSS', '[N II]/[S II]', 'WHAM Hα', 'paper Hα']
    ypos = {r: len(rows) - i for i, r in enumerate(rows)}

    def pt(row, dy, v, e, st, xerr_asym=None):
        y = ypos[row] + dy
        xe = e if xerr_asym is None else np.array([[v - xerr_asym[0]], [xerr_asym[1] - v]])
        ax.errorbar(v, y, xerr=xe, fmt=st['marker'], color=st['color'], mfc=st['mfc'], mec=st['color'],
                    ms=st['ms'], elinewidth=1.2, capsize=2.5, zorder=5)

    # adopted band and budget
    A = I['adopted']
    ax.axvspan(A['v'] - A['e'], A['v'] + A['e'], color=OK, alpha=0.13, zorder=0)
    ax.axvline(A['v'], color=OK, lw=0.8, alpha=0.6, zorder=0)
    for r in ('P1', 'P2', 'P3'):
        q = R[r]
        pt(r, +0.27, *q['pub_vhelio'], S['pub'])
        pt(r, +0.09, *q['fig3b_plotted'], S['plot'])
        pt(r, -0.09, *q['corr_vhelio'], S['corr'])
        pt(r, -0.27, *q['desi_0p5'][:2], S['desi'])
    o = V['oiii_fig3a']
    b = V['budget_fig3a']
    pt('Fig. 3a', +0.12, o['V_plot'], o['e_V_plot'], S['plot'])
    pt('Fig. 3a', -0.12, b['V_helio'], b['err'], S['corr'])
    ax.annotate('', xy=(b['V_helio'] + 1.5, ypos['Fig. 3a'] - 0.12), xytext=(o['V_plot'] - 1.0, ypos['Fig. 3a'] + 0.12),
                arrowprops=dict(arrowstyle='->', color=FGDIM, lw=1.0), zorder=4)
    label(ax, 28.0, ypos['Fig. 3a'] + 0.02,
          f"−{-b['sky_zero_point']:.1f} sky ref.\n−{-b['solar_blend']:.1f} solar blend\n−{-b['v_corr']:.1f} missing v$_{{corr}}$",
          color=FGDIM, size=8.2, va="center", ha="left")
    pt('DESI DR1', +0.12, I['desi_stack']['v'], I['desi_stack']['e_tot'], S['desi'])
    pt('DESI DR1', -0.12, I['desi_rev2']['v'], None, S['desi'], xerr_asym=(I['desi_rev2']['p16'], I['desi_rev2']['p84']))
    pt('SDSS', +0.12, I['sdss_all']['v'], I['sdss_all']['e'], S['sdss'])
    pt('SDSS', -0.12, I['sdss_south']['v'], float(np.hypot(I['sdss_south']['e'], I['sdss_south']['e_zp'])), S['sdss'])
    # red lines and WIM: per-region points
    dys = np.linspace(0.27, -0.27, 4)
    for dy, r in zip(dys, ('P1', 'P2', 'P3', 'P4')):
        q = R[r]
        ax.errorbar(q['nii'][0], ypos['[N II]/[S II]'] + dy, xerr=q['nii'][1], fmt='o', color=NII_C, mfc=NII_C, ms=4,
                    elinewidth=1.0, capsize=2, zorder=5)
        ax.errorbar(q['sii'][0], ypos['[N II]/[S II]'] + dy, xerr=0, fmt='x', color=SII_C, ms=5, zorder=6)
        ax.errorbar(q['wham_nearest']['v_helio'], ypos['WHAM Hα'] + dy, xerr=q['wham_nearest']['ev'], fmt='o',
                    color=HA_C, mfc='none', ms=4.5, elinewidth=1.0, capsize=2, zorder=5)
        ax.errorbar(q['ha'][0], ypos['paper Hα'] + dy, xerr=q['ha'][1], fmt='o', color='#8a939c', mfc='#8a939c',
                    ms=4, elinewidth=1.0, capsize=2, zorder=5)
    ax.set_yticks([ypos[r] for r in rows])
    ax.set_yticklabels(['P1', 'P2', 'P3', 'Fig. 3a', 'DESI DR1', 'SDSS', 'paper [N II]/[S II]', 'WHAM Hα (WIM)',
                        'paper Hα'], fontsize=9)
    for t in ax.get_yticklabels():
        t.set_color(FG)
    for yy in (ypos['Fig. 3a'] - 0.5, ypos['SDSS'] - 0.5):
        ax.axhline(yy, color='#2a323c', lw=0.9)
    ax.set_xlim(2, 36)
    ax.set_ylim(0.4, len(rows) + 0.75)
    ax.set_xlabel('V$_{helio}$ (km/s)')
    ax.set_title('(c) [O III] velocity ladder', loc='left', color=FG)
    sec = ax.secondary_xaxis('top', functions=(lambda v: v - 12.19, lambda v: v + 12.19))
    sec.set_xlabel('V$_{LSR}$ (km/s, ≈ V$_{helio}$ − 12.2)', color=FGDIM, fontsize=8.8)
    sec.tick_params(labelsize=8.5)
    label(ax, A['v'] + 0.25, len(rows) + 0.5, f"adopted {A['v']:+.1f} ± {A['e']:.1f}", color=OK, size=8.4, va='center')
    h = [Line2D([], [], ls='none', marker=s['marker'], color=s['color'], mfc=s['mfc'], mec=s['color'], ms=s['ms'])
         for s in (S['pub'], S['plot'], S['corr'], S['desi'], S['sdss'])]
    h += [Line2D([], [], ls='none', marker='x', color=SII_C, ms=5)]
    ax.legend(h, ['paper Table 4', 'paper Fig. 3, as plotted', 'corrected (−10.1)', 'DESI DR1', 'SDSS sky fibres',
                  '[S II]'], loc='center right', bbox_to_anchor=(0.995, 0.445), fontsize=8.2, ncol=2, framealpha=0.92,
              columnspacing=1.0, handletextpad=0.4)
    ax.grid(axis='x', alpha=0.25)


def main():
    fig = plt.figure(figsize=(15, 9.4))
    fig.text(0.012, 0.985, '[O III] velocity: a wavelength-reference artefact', color='#cfd6dd', fontsize=12.5,
             fontweight='bold', va='top')
    fig.text(0.012, 0.955, 'Wu et al. 2026 Fig. 3 re-extracted from the PDF vectors · LAMOST MRS-N, 2024-02-01 · '
             'v$_{corr}$ = −8.32 km/s', color=FGDIM, fontsize=9.5, va='top')
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.0], height_ratios=[1, 1], left=0.055, right=0.985,
                  top=0.875, bottom=0.07, wspace=0.2, hspace=0.42)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[1, 0]))
    panel_c(fig.add_subplot(gs[:, 1]))
    finish(fig, OUT, dpi=200, tight=False)
    print(OUT)


if __name__ == '__main__':
    main()
