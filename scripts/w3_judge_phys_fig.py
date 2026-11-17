"""w3_judge_phys figure: hypothesis posteriors and the distance-independence of the shock closure.

Reads data/derived/w3_judge_phys_checks.json (run w3_judge_phys_checks.py first).
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_phys_fig.py
"""
import json
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *  # noqa: F401,F403
import numpy as np

ROOT = '/Users/mickey/astronomy/silk-nebula/'
ck = json.load(open(ROOT + 'data/derived/w3_judge_phys_checks.json'))

# Posteriors as adopted in notes/w3_judge_phys.md (judgement, not computed)
H = [('H1 Monogem limb, truncated shock (250-320 pc)', 0.64, 'favored'),
     ('H2 other remnant/bubble, truncated shock', 0.26, 'viable'),
     ('H4 superbubble wall at 1-2 kpc', 0.04, 'disfavored'),
     ('H6 conduction/TML/NEI/WD or EUV photoion./fossil', 0.025, 'disfavored'),
     ('H9 unanticipated', 0.012, 'viable'),
     ('H5 DIG / photoionised bubble wall', 0.01, 'ruled_out'),
     ('H3 paper model: 160 pc, 2e51 erg, PDS', 0.005, 'ruled_out'),
     ('H7 stellar/pulsar bow shock or wake', 0.005, 'ruled_out'),
     ('H8 IVC/HVC infall shock', 0.003, 'ruled_out')]
col = {'favored': OK, 'viable': ACCENT, 'disfavored': WARN, 'ruled_out': DANGER}

fig, axes = new_sky_figure(13.2, 5.6, 'Silk Nebula · origin hypotheses, physical-consistency judge',
                           'Posterior (judgement) and the shock-closure distance scan (derived from verify_shockmodel runs)',
                           nrows=1, ncols=2, gridspec_kw={'width_ratios': [1.25, 1]})
ax = axes[0]
y = np.arange(len(H))[::-1]
for yi, (name, p, v) in zip(y, H):
    ax.barh(yi, p, color=col[v], alpha=0.85, height=0.62)
    label(ax, max(p, 0.004) * 1.12, yi, f'{p:.3f}', size=8, va='center')
ax.set_yticks(y)
ax.set_yticklabels([h[0] for h in H], fontsize=8)
ax.set_xscale('log')
ax.set_xlim(2e-3, 1.4)
ax.set_xlabel('posterior probability (log)')
ax.grid(axis='x', alpha=0.4)
for v, c in col.items():
    ax.barh([-5], [0], color=c, label=v.replace('_', ' '))
ax.set_ylim(-0.6, len(H) - 0.4)
ax.legend(loc='lower right', ncol=2, fontsize=7.5)

ax = axes[1]
sc = ck['closureW_distance_scan']['by_distance_pc']
ds = np.array(sorted(float(k) for k in sc))
lo = np.array([sc[str(int(d))]['Pram_k_range'][0] for d in ds])
hi = np.array([sc[str(int(d))]['Pram_k_range'][1] for d in ds])
med = np.array([sc[str(int(d))]['Pram_k_median'] for d in ds])
nm = [sc[str(int(d))]['n_match'] for d in ds]
bv = [sc[str(int(d))]['b_values'] for d in ds]
ax.fill_between(ds, lo, hi, color=ACCENT, alpha=0.25, label='ram P/k of matching models')
ax.plot(ds, med, 'o-', color=ACCENT, ms=4, label='median')
ax.axhspan(1.5e4, 1.6e5, color=SUN, alpha=0.08, label='Monogem P/k (definition spread)')
for d, n, b, m in zip(ds, nm, bv, med):
    label(ax, d, m * 1.9, f'n={n}\nb={int(min(b))}-{int(max(b))}', size=7, ha='center', va='bottom', color=FG)
for d0, d1, txt, c in ((150, 170, 'LB wall', FGDIM), (252, 318, 'Monogem limb', OK), (1000, 2000, 'H4', WARN)):
    ax.axvspan(d0, d1, color=c, alpha=0.12)
    label(ax, np.sqrt(d0 * d1), 1.15e4, txt, size=7.5, ha='center', color=c)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(130, 2100); ax.set_ylim(9e3, 5e5)
ax.set_xlabel('assumed distance (pc)')
ax.set_ylabel('driving pressure P/k (K cm$^{-3}$)')
ax.grid(alpha=0.35)
ax.legend(loc='upper right', fontsize=7.2, ncol=2)
finish(fig, ROOT + 'figures/w3_judge_phys_posteriors.png')
print('wrote figures/w3_judge_phys_posteriors.png')
