"""w4_geostats figure: (a) spine distance from the Monogem centres with the eRASS1 sector edges; (b) eRASS1 radial
profiles in the Silk sector; (c) chance-coincidence probabilities under the unified null; (d) south cross-profiles;
(e) ionisation-order offsets against their injection nulls; (f) profile asymmetry vs sheet curvature parameter q.
Inputs: data/derived/w4_geostats_{registration,order,shape}.json + _registration_erass.npz, _order_profiles.npz.
Output: figures/w4_geostats_summary.png (dpi 200)
Run: ~/.venvs/astro313/bin/python -I scripts/w4_geostats_fig.py
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w4_geostats_common import DER, FIG, ra_sp, dec_sp, CENTRES, sep_pa  # noqa
from silk_style import new_sky_figure, finish, label, OIII_C, WARN, OK, FG, FGDIM, ACCENT, DANGER, SII_C  # noqa

R = json.loads((DER / 'w4_geostats_registration.json').read_text())
O = json.loads((DER / 'w4_geostats_order.json').read_text())
S = json.loads((DER / 'w4_geostats_shape.json').read_text())
Z = np.load(DER / 'w4_geostats_registration_erass.npz')
PZ = np.load(DER / 'w4_geostats_order_profiles.npz')
VIOLET = '#b388ff'
COL = {'K1': WARN, 'K2': OIII_C, 'PSR': OK}

fig, axes = new_sky_figure(15, 9.6, 'Silk Nebula · geometry statistics',
                           'one null per question: Monogem-radius coincidence · ionisation order · cross-profile shape',
                           nrows=2, ncols=3)

# (a) spine distance from each centre with sector edges
ax = axes[0, 0]
for k in ('K1', 'K2', 'PSR'):
    d, _ = sep_pa(*CENTRES[k], ra_sp, dec_sp)
    ax.plot(dec_sp, d, color=COL[k], lw=1.5, label=k)
    e = R['erass_radius'][k]
    rh = e['sector']['r_half']; rp = e['sector']['r_peak']
    ax.hlines(rh, 12, 17.6, color=COL[k], lw=3.5, alpha=0.35)
    ax.hlines(rp, 12, 17.6, color=COL[k], lw=1.0, ls=(0, (2, 2)), alpha=0.8)
ax.axhspan(11.0 - 1.1, 11.0 + 1.1, color=OIII_C, alpha=0.07)
ax.axvline(17.6, color=FGDIM, lw=0.6)
ax.set_xlabel('spine Dec (deg)'); ax.set_ylabel('distance from centre (deg)')
ax.set_ylim(7.5, 17); ax.set_xlim(12, 28.3)
ax.legend(loc='upper left', fontsize=7.5, ncol=3)
label(ax, 18.0, 8.0, 'thick: eRASS1 outer half-max (Silk sector)\ndashed: shell peak · band: Knies K2 R ± σ', size=7, color=FGDIM)
ax.set_title('(a) south spine sits on the local X-ray half-max edge', loc='left', fontsize=9)

# (b) eRASS1 radial profiles in the Silk sector
ax = axes[0, 1]
for k in ('K1', 'K2', 'PSR'):
    rc = Z[f'{k}_rc']; p = Z[f'{k}_prof_s']
    ax.plot(rc, p, color=COL[k], lw=1.3, label=f'{k} sector')
    ax.axvline(R['erass_radius'][k]['silk_south_mean_dist'], color=COL[k], lw=0.8, ls=':')
ax.set_xlabel('radius from centre (deg)'); ax.set_ylabel('eRASS1 0.2–0.6 keV (cts s⁻¹ deg⁻²)')
ax.set_xlim(4, 18)
ax.legend(loc='upper right', fontsize=7.5)
le = R['local_edge_null']
txt = '\n'.join(f"{k}: Silk − r½ = {le[k]['delta_obs']:+.2f}°, rotation p = {le[k]['p_abs_le_obs']:.3f}" for k in ('K1', 'K2', 'PSR'))
label(ax, 4.3, ax.get_ylim()[0] + 0.05 * np.diff(ax.get_ylim())[0], txt, size=7, va='bottom')
ax.set_title('(b) X-ray shell in the Silk sector (dotted: Silk)', loc='left', fontsize=9)

# (c) coincidence probabilities
ax = axes[0, 2]
rows = [('K2 only', R['nulls']['K2_only']), ('K1/K2/PSR', R['nulls']['primary']), ('  + mirror', R['nulls']['mirror']),
        ('  band ±6°', R['nulls']['band6']), ('  PSR σ wide', R['nulls']['psr_wide']), ('  eRASS radii', R['nulls']['erass_radii']),
        ('  + G205.6', R['nulls']['with_G205'])]
y = np.arange(len(rows))[::-1]
p1 = [r[1]['r0=0.15|z0=1.0']['p_given_touch1.5'] for r in rows]
p3 = [r[1]['r0=0.15|z0=1.0']['p_given_touch3.0'] for r in rows]
ax.barh(y + 0.18, p1, height=0.34, color=OIII_C, label='rim zone ±1.5°')
ax.barh(y - 0.18, p3, height=0.34, color=VIOLET, label='rim zone ±3°')
for yy, a in zip(y, p1):
    label(ax, a + 0.003, yy + 0.18, f'{a:.3f}', size=7, va='center')
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
ax.set_xlabel('P(E | H0, trace reaches rim)')
ax.set_xlim(0, 0.17)
rp = R['reproduce_judges']
label(ax, 0.004, -1.0, f"judge_phys def.: P = {rp['phys_theta20']['p']:.4f}–{rp['phys_theta30']['p']:.4f} → BF {rp['phys_theta20']['bf_with_limb_031']:.0f}–{rp['phys_theta30']['bf_with_limb_031']:.0f}\n"
      f"judge_skeptic def. (K2): P = {rp['skeptic_K2']['p_given_touch']:.3f}", size=6.8, color=FGDIM, va='center')
ax.legend(loc='upper right', fontsize=7)
ax.set_ylim(-1.6, len(rows) - 0.4)
ax.set_title('(c) chance coincidence, one null', loc='left', fontsize=9)

# (d) south S1 profiles
ax = axes[1, 0]
x25 = PZ['x25']; x10 = PZ['x10']


def norm(p, xx):
    base = np.nanmedian(p[(np.abs(xx) > 0.6) & (np.abs(xx) < 1.2)])
    q = p - base
    return q / np.nanmax(q[np.abs(xx) < 0.3])


ax.plot(x25, norm(PZ['S1_nsns'], x25), color=OIII_C, lw=1.8, label='NSNS [O III]')
ax.plot(x25, norm(PZ['S1_FUV'], x25), color=VIOLET, lw=1.2, label='FUV')
ax.plot(x25, norm(PZ['S1_NUV'], x25), color=OK, lw=1.0, alpha=0.85, label='NUV')
ax.plot(x10, norm(PZ['S1_desi_o3'], x10), 'o', color=OIII_C, ms=3.5, mfc='none', label='DESI [O III]')
ax.plot(x10, norm(PZ['S1_desi_o2'], x10), 's', color=WARN, ms=3.5, mfc='none', label='DESI [O II]')
ax.axvline(0, color=FGDIM, lw=0.6)
ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.5, 1.4)
ax.set_xlabel('offset from spine (deg; + = east, outside)'); ax.set_ylabel('normalised')
ax.legend(loc='upper left', fontsize=6.8, ncol=2)
f = O['S1']['FUV']['fit']
label(ax, 0.32, 0.80, f"FUV split normal: mode {f['x0'] - O['S1']['nsns_fit']['x0']:+.3f}°\nσW/σE = {f['ratio_we']:.2f}", size=7, va='top')
ax.set_title('(d) Dec 12–15 cross-profiles', loc='left', fontsize=9)

# (e) offsets vs nulls
ax = axes[1, 1]
items = []
for seg, tag in (('S1', 'S1'), ('S12', 'S12')):
    for b, c in (('FUV', VIOLET), ('NUV', OK)):
        r = O[seg][b]
        t = r['test_coshape_delta']
        items.append((f'{b} δ {tag}', t['obs'], t['median'], t['sd'], r['inj_pred']['median'], c, t['p_one']))
        m = r['test_mode_displacement']
        items.append((f'{b} mode {tag}', m['obs'], m['median'], m['sd'], r['inj_pred']['delta_true'], c, m['p_one']))
    d = O[seg]['desi']
    items.append((f'[O II]−[O III] {tag}', d['d_O2_vs_O3']['obs'], d['d_O2_vs_O3']['orient_null_median'],
                  d['d_O2_vs_O3']['orient_null_sd'], d['d_O2_vs_O3']['inj_pred']['median'], WARN, d['d_O2_vs_O3']['orient_p_two']))
d = O['S1']['desi']
items.append(('DESI−NSNS [O III] S1', d['d_O3_vs_NSNS']['obs'], d['d_O3_vs_NSNS']['orient_null_median'],
              d['d_O3_vs_NSNS']['orient_null_sd'], np.nan, OIII_C, d['d_O3_vs_NSNS']['orient_p_two']))
yy = np.arange(len(items))[::-1]
for y0, (nm, o, med, sd, pr, c, pp) in zip(yy, items):
    ax.plot([med - sd, med + sd], [y0, y0], color=FGDIM, lw=5, alpha=0.35, solid_capstyle='butt')
    ax.plot([med - 2 * sd, med + 2 * sd], [y0, y0], color=FGDIM, lw=1, alpha=0.6)
    ax.plot(o, y0, 'o', color=c, ms=6)
    if np.isfinite(pr):
        ax.plot(pr, y0, 'x', color=DANGER, ms=6, mew=1.5)
    label(ax, 0.148, y0, f'p {pp:.2f}', size=6.8, va='center')
ax.set_yticks(yy); ax.set_yticklabels([i[0] for i in items], fontsize=7)
ax.axvline(0, color=FGDIM, lw=0.6)
ax.set_xlim(-0.16, 0.2)
ax.set_xlabel('offset (deg; + = east)   p: empirical, placement/orientation null')
ax.plot([], [], 'x', color=DANGER, label='shock prediction (injected)')
ax.plot([], [], color=FGDIM, lw=5, alpha=0.35, label='co-spatial null ±1σ')
ax.legend(loc='lower left', fontsize=6.5)
ax.set_title('(e) ionisation order vs injection nulls', loc='left', fontsize=9)

# (f) asymmetry vs q
ax = axes[1, 2]
qs = np.array(S['S1']['likelihood_vs_q']['sb0.03']['q'])
qplot = np.where(qs == 0, 0.04, qs)
for seg, c in (('S1', OIII_C), ('N', WARN), ('N2', VIOLET)):
    for sb, ls in (('sb0.03', '-'), ('sb0.06', '--')):
        L = S[seg]['likelihood_vs_q'][sb]
        m = np.array(L['ratio_mean']); sd = np.array(L['ratio_sd'])
        ax.plot(qplot, m, color=c, lw=1.1, ls=ls, alpha=0.9)
        if seg == 'S1' and sb == 'sb0.03':
            ax.fill_between(qplot, m - sd, m + sd, color=c, alpha=0.12)
    o = S[seg]['observed']['primary']
    r = o['fit']['ratio_we']; b = o['boot_ratio']
    ax.axhspan(b['p16'], b['p84'], color=c, alpha=0.18)
    ax.axhline(r, color=c, lw=1.4)
    label(ax, 0.045 if seg != 'N2' else 0.3, r + 0.06, f"{seg} observed {r:.2f}", size=7, color=c)
ax.set_xscale('log')
ax.set_xlabel('q = LOS sagitta / layer depth  (0 = flat ribbon)')
ax.set_ylabel('σW/σE, injected into real sky (mean)')
ax.set_ylim(0.6, 4.5)
G = S['S1']['likelihood_vs_q']['sb0.03']['G']
for qq, gg in zip(qplot[1:], G[1:]):
    ax.text(qq, 4.38, f'{gg:.0f}', color=FGDIM, fontsize=6.5, ha='center', va='top')
ax.text(0.04, 4.38, 'G(S1):', color=FGDIM, fontsize=6.5, ha='left', va='top')
ax.plot([], [], color=FGDIM, ls='-', label='model, blur σ 0.03°'); ax.plot([], [], color=FGDIM, ls='--', label='model, blur σ 0.06°')
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.93), fontsize=7)
ax.set_title('(f) cross-profile asymmetry: ribbon → limb', loc='left', fontsize=9)

out = finish(fig, str(FIG / 'w4_geostats_summary.png'))
print('saved', out)
