"""Figure: SDSS sky_g amplitude vs RA shift of the NSNS [OIII] template (observed vs injected signal).
Input: data/derived/imaging_optical_sdss_placebo_w4.0.json -> figures/imaging_optical_sdss_shift.png"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
import silk_style as ss
from imaging_optical_common import DER, FIG

d = json.loads((DER / 'imaging_optical_sdss_placebo_w4.0.json').read_text())
s = np.array(d['shifts_deg']); x = np.array(d['A_over_expected']); e = np.array(d['err_over_expected'])
xi = np.array(d['injected_A_over_expected']); k = d['curve_scale_obs_over_injected']
fig, ax = ss.new_sky_figure(8.5, 4.8, 'Silk Nebula · SDSS template-shift test',
                            f"{d['n_runs']} runs, Dec 12–20.5, fields ≤{d['window_deg']:.0f}° from ridge")
fig.texts[1].set_y(0.935)
ax.fill_between(s, x - e, x + e, color=ss.OIII_C, alpha=0.25, lw=0)
ax.plot(s, x, 'o-', color=ss.OIII_C, ms=3.5, lw=1.2, label='observed g sky')
ax.plot(s, xi, '--', color=ss.FGDIM, lw=1.0, label='injected NSNS × 1.0')
ax.plot(s, k * xi, '-', color='#ffffff', lw=1.0, label=f'injected × {k:.2f}')
ax.plot(s, 2.7 * xi, ':', color=ss.DANGER, lw=1.2, label='injected × 2.7 (paper, min.)')
ax.axhspan(d['placebo_mean_absshift_ge2'] - d['placebo_std_absshift_ge2'],
           d['placebo_mean_absshift_ge2'] + d['placebo_std_absshift_ge2'], color=ss.FGDIM, alpha=0.12, lw=0)
ss.label(ax, 2.25, -0.95, 'null ±1σ (|shift| ≥ 2°)', size=7, color=ss.FGDIM, ha='center')
ax.axhline(0, color='#2a323c', lw=0.8)
ax.set_xlabel('RA shift of template [deg on sky]'); ax.set_ylabel('A / A(NSNS calibration)')
ax.set_xlim(-3.05, 3.05); ax.set_ylim(-1.2, 3.0); ax.grid(alpha=0.3)
ax.legend(fontsize=7.5, loc='upper left', ncol=2)
ss.finish(fig, FIG / 'imaging_optical_sdss_shift.png')
