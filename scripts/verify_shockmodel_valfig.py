"""verify_shockmodel_valfig: validation figure for the shockmodel 1-D code.
  (a) isobaric non-equilibrium cooling from 5e6 K: O2+ fraction and cooling efficiency vs Gnat & Sternberg 2007
  (b) model vs published line ratios (energy, rel. Hb): complete shocks (Allen+08 Table 6; Raymond 1979 E-I) and
      truncated shocks (Raymond+88 Table 5 SE/M), plus R88 Fig. 8a peak values (80 km/s)
Output: figures/verify_shockmodel_validation.png
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *

ROOT = '/Users/mickey/astronomy/silk-nebula'
val = json.load(open(ROOT + '/data/derived/verify_shockmodel_validation.json'))
r88 = json.load(open(ROOT + '/data/derived/verify_shockmodel_r88.json'))

fig, axs = new_sky_figure(12.5, 5.2, 'Silk · shockmodel code validation',
                          'left: isobaric NEI cooling vs Gnat & Sternberg 2007 · right: line ratios vs published shock models',
                          nrows=1, ncols=2)
ib = val['isobaric']
ax = axs[0]
ax.plot(ib['curve_T'], ib['curve_O2'], color=OIII_C, lw=1.6, label='O²⁺/O this code')
ax.plot(ib['gs07_T'], ib['gs07_O2'], color=OIII_C, lw=1.0, ls='--', label='O²⁺/O GS07')
ax.set_xscale('log'); ax.set_xlim(1e4, 2e6); ax.set_ylim(0, 1)
ax.set_xlabel('T (K)'); ax.set_ylabel('ion fraction')
ax2 = ax.twinx()
ax2.plot(ib['curve_T'], ib['curve_L'], color=WARN, lw=1.4, label='Λ this code')
ax2.plot(ib['gs07_cefT'], ib['gs07_cef'], color=WARN, lw=1.0, ls='--', label='Λ GS07')
ax2.set_yscale('log'); ax2.set_ylim(1e-23, 1e-21); ax2.set_ylabel('Λ = L/(n_e n_H) (erg cm³ s⁻¹)', color=WARN)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=7, loc='upper right')
ax.grid(alpha=0.3)

ax = axs[1]
pts = []
mk = {'O': ('o', OIII_C), 'C': ('s', WARN), 'T': ('D', HA_C)}
names = {'OIII5007': '[O III]', 'OII3727': '[O II]', 'CIV1549': 'C IV', 'CIV1550': 'C IV', 'Ha': 'Hα', 'NII6583': '[N II]',
         'NII6584': '[N II]', 'OIII4363': '4363', 'HeII1640': 'He II', 'OIII1664': 'O III]', 'SII6716_31': '[S II]'}
for key, rec in val['allen'].items():
    if not key.endswith('_7000'):
        continue
    for ln in ('OIII5007', 'OII3727', 'CIV1549', 'Ha', 'NII6583', 'OIII4363', 'HeII1640', 'OIII1664'):
        ax.plot(rec['allen08'][ln], rec['mine'][ln], 'o', color=OIII_C, ms=4, alpha=0.8)
for v, rec in r88['R79'].items():
    for ln, k in (('O3', '[O III]'), ('O2', '[O II]'), ('C4', 'C IV')):
        ax.plot(rec['R79'][ln] / 100, rec['mine'][ln] / 100, 's', color=WARN, ms=4, alpha=0.8)
for key in ('SE_preO+O++', 'M_preO+O++'):
    rec = r88['R88'][key]
    m = rec['model']['at_R88_N']['ratios']
    for ln in ('OIII5007', 'OII3727', 'OIII4363', 'Ha', 'NII6584', 'NeIII3869', 'HeII1640', 'CIV1550', 'OIII1664'):
        ax.plot(rec['R88_table5'][ln] / 100, m[ln] / 100, 'D', color=HA_C, ms=4, alpha=0.9)
        if ln in ('OIII5007', 'OII3727', 'CIV1550'):
            label(ax, rec['R88_table5'][ln] / 100 * 1.08, m[ln] / 100, f"{key.split('_')[0]} {names[ln]}", color=HA_C, size=6)
f8 = r88['R88_fig8a_80kms']
for k in ('O3', 'O2', 'O3s'):
    ax.plot(10**f8[k]['r88_peak_log'], 10**f8[k]['mine_peak_log'], '^', color=OK, ms=6)
    ax.plot(10**f8[k]['r88_end_log'], 10**f8[k]['mine_end_log'], 'v', color=OK, ms=6)
xx = np.logspace(-2, 2.2, 10)
ax.plot(xx, xx, color=FGDIM, lw=0.8)
ax.fill_between(xx, xx / 2, xx * 2, color='#9aa4ae', alpha=0.12, lw=0)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.02, 120); ax.set_ylim(0.02, 120)
ax.set_xlabel('published ratio to Hβ'); ax.set_ylabel('this code, ratio to Hβ')
for lab, (m_, c_) in (('Allen+08 complete 100/200 km/s', ('o', OIII_C)), ('Raymond 1979 complete 100-200', ('s', WARN)),
                      ('Raymond+88 truncated SE/M (at their N)', ('D', HA_C)),
                      ('R88 Fig.8a 80 km/s: peak ▲ / complete ▼', ('^', OK))):
    ax.plot([], [], m_, color=c_, label=lab)
ax.legend(fontsize=7, loc='upper left')
label(ax, 6, 0.03, 'grey band = factor 2', color=FGDIM, size=7)
ax.grid(alpha=0.3)
finish(fig, ROOT + '/figures/verify_shockmodel_validation.png')
print('ok')
