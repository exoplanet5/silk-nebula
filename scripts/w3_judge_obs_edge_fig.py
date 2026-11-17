"""Summary figure for the judge_obs edge tests: observed statistics against their null bands.

Inputs: data/derived/w3_judge_obs_edge.json, w3_judge_obs_edge2.json (written by w3_judge_obs_edge*.py)
Output: figures/w3_judge_obs_edge_summary.png
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_obs_edge_fig.py
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import DER, FIG  # noqa: E402
import silk_style as ss  # noqa: E402

E1 = json.loads((DER / 'w3_judge_obs_edge.json').read_text())
E2 = json.loads((DER / 'w3_judge_obs_edge2.json').read_text())

fig, (a1, a2) = ss.new_sky_figure(12, 5.2, 'Silk Nebula · Edge and UV-offset tests',
                                  'Observed vs off-Silk null placements (grey: 16–84 % of nulls)', ncols=2)
# left: W/E width ratio
rows = [('S1 Dec 12–15', E1['segments']['S1']['fit']['SG']['ratio_WE'], E1['segments']['S1']['boot']['ratio_WE_16_50_84'],
         E1['segments']['S1']['sym_null']['ratio_WE_16_50_84'], E1['segments']['S1']['sym_null']['p_ratio_ge_obs'])]
for k, lab in (('N1', 'N1 Dec 20–21.5'), ('N3', 'N3 Dec 21.5–24.25'), ('N2', 'N2 Dec 24.25–28')):
    r = E2['nsns_north_subsegments'][k]
    rows.append((lab, r['ratio_WE'], r['boot_ratio_16_50_84'], r['null_ratio_16_50_84'], r['p_ratio']))
for i, (lab, obs, bt, nl, p) in enumerate(rows):
    y = len(rows) - 1 - i
    a1.fill_betweenx([y - 0.25, y + 0.25], nl[0], nl[2], color='#3a434d', lw=0)
    a1.errorbar(obs, y, xerr=[[obs - bt[0]], [bt[2] - obs]], fmt='o', color=ss.OIII_C, ms=6, capsize=3)
    ss.label(a1, max(bt[2], nl[2]) + 0.08, y, 'p = %.2f' % p, size=8, va='center')
a1.axvline(1.0, color=ss.FGDIM, lw=0.8, ls=':')
a1.set_yticks(range(len(rows)))
a1.set_yticklabels([r[0] for r in rows][::-1])
a1.set_xlim(0.6, 3.6)
a1.set_xlabel('[O III] width ratio σ_W / σ_E (>1: sharp east edge)')
a1.set_title('NSNS [O III] cross-ridge asymmetry', loc='left', fontsize=9)
# right: GALEX centre offsets
rows = []
for seg, lab in (('S1', 'Dec 12–15'), ('S12', 'Dec 12–17.6')):
    for band in ('FUV', 'NUV'):
        g = E2['galex_offset_null'][seg][band]
        v = list(g['window_baseline_variants_dx'].values())
        rows.append(('%s %s' % (band, lab), g['dx_obs'], (min(v), max(v)), g['null_dx_16_84'], g['z_vs_null']))
for i, (lab, obs, vr, nl, z) in enumerate(rows):
    y = len(rows) - 1 - i
    a2.fill_betweenx([y - 0.25, y + 0.25], nl[0], nl[1], color='#3a434d', lw=0)
    col = '#b388ff' if lab.startswith('FUV') else ss.WARN
    a2.errorbar(obs, y, xerr=[[obs - vr[0]], [vr[1] - obs]], fmt='o', color=col, ms=6, capsize=3)
    ss.label(a2, 0.125, y, 'z = %.1f' % z, size=8, va='center')
a2.axvline(0.0, color=ss.FGDIM, lw=0.8, ls=':')
a2.axvline(0.1, color=ss.OK, lw=1.0, ls='--')
ss.label(a2, 0.102, len(rows) - 0.55, 'shock prediction', size=7.5, color=ss.OK)
a2.set_yticks(range(len(rows)))
a2.set_yticklabels([r[0] for r in rows][::-1])
a2.set_xlim(-0.07, 0.16)
a2.set_xlabel('UV ridge − [O III] ridge (deg, + = east); bars: fit variants')
a2.set_title('GALEX ridge offset', loc='left', fontsize=9)
ss.finish(fig, str(FIG / 'w3_judge_obs_edge_summary.png'))
print('saved')
