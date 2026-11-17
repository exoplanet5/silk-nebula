"""verify_velocity: summary figure.

(a) The paper's own Fig. 3a spectrum (vector-extracted): [N I] 5198/5200 airglow and [O III] 5007 on the
    published wavelength scale, converted to velocity w.r.t. air rest.  Sky lines sit at ~0 (topocentric scale).
(b) [O III] velocity comparison (heliocentric).
(c) LAMOST DR12 stellar RV blue-arm minus red-arm per MRS plate in the Silk box and on nearby nights.
(d) Halpha centroid bias from geocoronal subtraction vs (true / subtracted) geocoronal intensity.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_figure.py
Output: figures/verify_velocity_summary.png
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, finish, label, FG, FGDIM, ACCENT, WARN, DANGER, OK, OIII_C, HA_C,
                        NII_C, SII_C)  # noqa: E402
import verify_velocity_halpha_sim as HS  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
C = 299792.458
VCORR = -8.32


def main():
    z = np.load(ROOT / 'data/derived/verify_velocity_fig3_spectra.npz')
    w, f = z['a_wave'], z['a_flux']
    o = np.argsort(w)
    w, f = w[o], f[o]
    fs = json.loads((ROOT / 'data/derived/verify_velocity_fig3_solar.json').read_text())
    la = json.loads((ROOT / 'data/derived/verify_velocity_lamost_arms.json').read_text())
    ds = json.loads((ROOT / 'data/derived/verify_velocity_desi_sdss.json').read_text())
    hs = json.loads((ROOT / 'data/derived/verify_velocity_halpha_sim.json').read_text())

    fig, axs = new_sky_figure(13, 9.2, 'Silk Nebula [O III] velocity forensics',
                              'Paper Fig. 3 spectra recovered from the PDF vectors; LAMOST arms; DESI/SDSS re-measured; '
                              'geocoronal Hα simulation', nrows=2, ncols=2)
    # ---------------------------------------------------------------- (a)
    ax = axs[0, 0]
    lines = [('[N I] 5197.90', 5197.902, ACCENT, 5195.5, 5199.0), ('[N I] 5200.26', 5200.257, '#80deea', 5199.2, 5202.0),
             ('[O III] 5006.84', 5006.843, WARN, 5004.5, 5010.5)]
    for name, l0, col, lo, hi in lines:
        m = (w > lo) & (w < hi)
        v = C * (w[m] / l0 - 1)
        y = f[m]
        base = np.median(np.concatenate([y[:5], y[-5:]]))
        y = (y - base) / (y.max() - base)
        ax.plot(v, y, color=col, lw=1.4, drawstyle='steps-mid', label=name)
    ax.axvline(0, color=FGDIM, lw=0.8, ls=':')
    vo = fs['V_helio_calibration_free']['5007_plain_plot_scale']
    ax.axvline(vo, color=WARN, lw=0.8, ls='--')
    ax.axvspan(23.8, 26.8, color=DANGER, alpha=0.18, lw=0)
    label(ax, 52, 1.04, 'paper Table 4 [O III]\n+23.8…+26.8 (red band)', color=DANGER, size=8)
    label(ax, -106, 1.04, 'sky lines at ≈0:\nscale is topocentric (air)', color=ACCENT, size=8)
    label(ax, -106, 0.66, '[O III] vs [N I]: V_topo ≈ +22\n+ v_corr (−8.3) → V_helio ≈ +14', color=FG, size=8)
    ax.set_xlim(-110, 110)
    ax.set_ylim(-0.35, 1.25)
    ax.set_xlabel('velocity on the paper\'s plotted wavelength scale (km/s, air rest)')
    ax.set_ylabel('normalised flux')
    ax.legend(loc='lower right', fontsize=7.5)
    ax.set_title('(a) paper Fig. 3a spectrum: airglow vs [O III]', loc='left')

    # ---------------------------------------------------------------- (b)
    ax = axs[0, 1]
    items = []
    items.append(('paper Table 4 [O III] P1', 26.8, 4.3, DANGER))
    items.append(('paper Table 4 [O III] P2', 26.8, 4.1, DANGER))
    items.append(('paper Table 4 [O III] P3', 23.8, 3.3, DANGER))
    corr = 8.32 + 1.8
    items.append(('Table 4 − 10.1 (frame fix) P1', 26.8 - corr, 4.4, WARN))
    items.append(('Table 4 − 10.1 (frame fix) P2', 26.8 - corr, 4.2, WARN))
    items.append(('Table 4 − 10.1 (frame fix) P3', 23.8 - corr, 3.4, WARN))
    v3a = 24.07 - 1.8 + VCORR
    items.append(('paper Fig. 3a vs own [N I] airglow', v3a, 1.4, OK))
    d = ds['desi']['pred_pos_gt0.5R']
    items.append(('DESI on-ridge stack (this work)', d['V_helio'], d['err_boot'], OIII_C))
    items.append(('DESI wave 1 (flux-weighted)', 14.5, 1.5, '#5f9ea0'))
    s = ds['sdss']['all']
    items.append(('SDSS sky fibres ON−OFF (this work)', s['V_helio_5007'], s['err_boot'], OIII_C))
    items.append(('SDSS wave 1 (south)', 12.5, 2.6, '#5f9ea0'))
    items.append(('WHAM Hα (V_LSR≈0)', 12.2, 1.0, HA_C))
    items.append(('paper [N II] P1–P4', 12.2, 0.8, NII_C))
    items.append(('paper Hα P1–P4', 17.9, 1.5, HA_C))
    y = np.arange(len(items))[::-1]
    for yi, (name, v, e, col) in zip(y, items):
        ax.errorbar(v, yi, xerr=e, fmt='o', color=col, ms=4.5, capsize=2, lw=1.2)
        label(ax, -14.5, yi - 0.1, name, color=col, size=7.5)
    ax.axvspan(13.3 - 1.1, 13.3 + 1.1, color=OK, alpha=0.15, lw=0)
    ax.set_xlim(-15, 35)
    ax.set_yticks([])
    ax.set_xlabel('V_helio (km/s)')
    ax.set_title('(b) [O III] heliocentric velocity: best +13.3 ± 1.1', loc='left')

    # ---------------------------------------------------------------- (c)
    ax = axs[1, 0]
    bp = la['box_plates_b0_minus_r0']
    dates = [r['obsdate'] for r in bp]
    yy = [r['median'] for r in bp]
    ee = [r['mad_sigma'] / np.sqrt(r['n']) for r in bp]
    xx = np.array([int(dd[:4]) + (int(dd[5:7]) - 0.5) / 12 for dd in dates])
    ax.errorbar(xx, yy, yerr=ee, fmt='o', color=ACCENT, ms=3.5, lw=0.8, label='MRS plates in the Silk box')
    nn = la['nearby_nights_plates_b0_minus_r0']
    xn = np.array([2024.08 + 0.002 * i for i in range(len(nn))])
    ax.errorbar(xn, [r['median'] for r in nn], yerr=[r['mad_sigma'] / np.sqrt(r['n']) for r in nn], fmt='s',
                color=FGDIM, ms=3, lw=0.8, label='plates 2024-01-26…02-07')
    mn = la['MN_plate_coadd_b0_minus_r0_all']
    ax.errorbar([2024.09], [mn['median']], yerr=[mn['median_boot_sigma']], fmt='*', color=WARN, ms=13, lw=1.2,
                label='MN074209N160100N01 (paper plate)')
    ax.axhspan(10, 15, color=DANGER, alpha=0.15, lw=0)
    label(ax, 2018.4, 12.0, 'offset needed to explain\npaper [O III] (+10…+15)', color=DANGER, size=8)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.set_ylim(-6, 16)
    ax.set_xlabel('observation date')
    ax.set_ylabel('stellar RV blue arm − red arm (km/s)')
    ax.legend(loc='lower right', fontsize=7.5)
    ax.set_title('(c) LAMOST DR12 per-arm stellar RVs: no blue-arm offset', loc='left')

    # ---------------------------------------------------------------- (d)
    ax = axs[1, 1]
    kts = np.arange(0.4, 1.61, 0.05)
    for dm, col, nm in ((0.0, ACCENT, 'model at RV(OH 6554)'), (-8.24, WARN, 'model at λ6554×6563/6554')):
        for g, ls in ((0.9, '-'), (1.3, '--')):
            b = [HS.simulate(g=g * kt, k=1 / kt, dm=dm)[0] for kt in kts]
            ax.plot(kts, b, color=col, ls=ls, lw=1.4, label=f'{nm}, g={g}')
    ax.axhspan(4.5, 6.7, color=HA_C, alpha=0.2, lw=0)
    label(ax, 1.15, 5.2, 'paper Hα − [N II]\n+4.5…+6.7 km/s', color=HA_C, size=8)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.axvline(1, color=FGDIM, lw=0.6, ls=':')
    b = hs['best_fit_dm_kt']
    ax.plot([b['k_true_over_model']], [5.6], '*', color=OK, ms=12)
    label(ax, b['k_true_over_model'] + 0.03, 6.9, f"best fit: true/model = {b['k_true_over_model']:.2f}", color=OK, size=8)
    ax.set_xlabel('true / subtracted geocoronal Hα intensity')
    ax.set_ylabel('Hα centroid bias (km/s)')
    ax.set_ylim(-10, 14)
    ax.legend(loc='upper right', fontsize=7)
    ax.set_title('(d) geocoronal over-subtraction pushes Hα redward', loc='left')
    finish(fig, ROOT / 'figures/verify_velocity_summary.png', dpi=200)


if __name__ == '__main__':
    main()
