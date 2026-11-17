"""dispute_distance_kd.py -- Section 1 of the Comment on Wu et al. 2026 (AJ 172, 298):
"The 160 pc distance: Table 4 is in kpc, not '100 pc'".

Everything here is derived from the paper's OWN numbers (Table 1 positions, Table 4 velocities and
distances, Fig. 5b [O III] V_LSR histogram, Sect. 3.3 / 4.1 / 4.2 text) plus the rotation curve the
paper names (Reid et al. 2014 URC, as used in the kd package of Wenger et al. 2018).
No new observational data are used.

Re-used (imported, not modified):
  scripts/verify_paper_kd.py  -- Table 4 transcription (T4, REG), Reid+14 URC in the kd form, solar motion,
                                 traditional-distance solver
  scripts/paper_review_lib.py -- further rotation-curve variants (flat, Reid+14/Reid+19 A5 linear, kd Reid+19 URC)
  data/derived/verify_paper_kd.json -- kd.pdf_kd emulation (MC modes/HPD), cross-checked here

Outputs: data/derived/dispute_distance_kd.json, figures/dispute_distance_fig1.png
Run:     ~/.venvs/astro313/bin/python -I scripts/dispute_distance_kd.py
"""
import sys
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm, chi2 as chi2dist

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import verify_paper_kd as vk          # noqa: E402  (T4, REG, vlsr_model, helio_to_lsr)
import paper_review_lib as prl        # noqa: E402  (MODELS)
from astropy.coordinates import SkyCoord   # noqa: E402
import astropy.units as u             # noqa: E402

OUT_JSON = ROOT / 'data/derived/dispute_distance_kd.json'
OUT_FIG = ROOT / 'figures/dispute_distance_fig1.png'
FIG5B = ROOT / 'sources/arxiv/figs/fig5b_crop.png'

LINES = ['OIII', 'SII', 'NII', 'Ha']
LNAME = {'OIII': '[O III]', 'SII': '[S II]', 'NII': '[N II]', 'Ha': 'Hα'}
A_OORT, A_OORT_ERR = 15.3, 0.4        # Bovy 2017, MNRAS 468, L63 (km/s/kpc)
D = np.arange(0.0, 6.0001, 0.001)     # kpc


def lb_of(region):
    ra, dec = vk.REG[region]
    c = SkyCoord(ra * u.deg, dec * u.deg).galactic
    return float(c.l.deg), float(c.b.deg)


LB = {k: lb_of(k) for k in vk.REG}

# rotation-curve variants: name -> f(l, b, d[kpc]) = V_LSR (km/s)
MODELS = {
    'Reid+14 URC, kd form (main)': lambda l, b, d: vk.vlsr_model(l, b, d),
    'Reid+14 URC, x cos b': lambda l, b, d: vk.vlsr_model(l, b, d, proj='3d'),
    'flat 240, x cos b': lambda l, b, d: vk.vlsr_model(l, b, d, curve='flat', proj='3d'),
    'flat 236 (R0 8.15)': prl.MODELS['flat_R0_8.15_T236_stdsun'],
    'Reid+14 A5 linear': prl.MODELS['reid14_A5_linear'],
    'Reid+19 A5 linear': prl.MODELS['reid19_A5_linear'],
    'kd Reid+19 URC (kd default)': prl.MODELS['kd_reid19_URC'],
}
MAIN = 'Reid+14 URC, kd form (main)'


def vrot(region, d, model=MAIN):
    l, b = LB[region]
    return np.asarray(MODELS[model](l, b, np.asarray(d, float)), float)


def solve_d(region, v, model=MAIN, dmin=0.001):
    """Outer-Galaxy traditional kinematic distance: V(d) is monotonic increasing here (l = 202-206 deg)."""
    g = D[D >= dmin]
    vg = vrot(region, g, model)
    assert np.all(np.diff(vg) > 0), model
    if v < vg[0] or v > vg[-1]:
        return None
    return float(np.interp(v, vg, g))


def hpd_truncnorm(mu, sig, lo, frac=0.683):
    """68.3 % highest-density interval of N(mu, sig) truncated to V > lo (kd drops samples with no solution)."""
    v = np.linspace(lo, mu + 8 * sig, 40001)
    p = norm.pdf(v, mu, sig)
    p /= p.sum()
    order = np.argsort(p)[::-1]
    k = np.searchsorted(np.cumsum(p[order]), frac) + 1
    sel = v[order[:k]]
    return float(sel.min()), float(sel.max())


def digitize_fig5b():
    """Read the [O III] (brown) step histogram of the paper's Fig. 5b from the 300-dpi page crop.
    Axis calibration from the tick marks; bar levels from the topmost brown pixel per column."""
    import matplotlib.image as mpimg
    im = mpimg.imread(FIG5B)
    rgb = (im[..., :3] * 255).astype(int)
    blk = rgb.sum(-1) < 200
    xt = np.nonzero(blk[978:984].sum(0) >= 4)[0]
    xt = np.array([g.mean() for g in np.split(xt, np.nonzero(np.diff(xt) > 2)[0] + 1)])
    yt = np.nonzero(blk[:, 124:131].sum(1) >= 5)[0]
    xv = np.array([-10, 0, 10, 20, 30.])
    yv = np.array([60, 50, 40, 30, 20, 10, 0.])
    ax_, bx_ = np.polyfit(xt, xv, 1)
    ay_, by_ = np.polyfit(yt, yv, 1)
    brown = (np.abs(rgb[..., 0] - 166) < 30) & (np.abs(rgb[..., 1] - 86) < 30) & (np.abs(rgb[..., 2] - 40) < 35)
    cols = np.arange(int(xt[0]) - 70, 832)
    top = np.full(cols.size, np.nan)
    for i, x in enumerate(cols):
        ys = np.nonzero(brown[200:974, x])[0]
        if ys.size:
            top[i] = ys.min() + 200
    lw_half = 2.5   # half line width (px) of the drawn step outline
    level = ay_ * (top + lw_half) + by_
    vel = ax_ * cols + bx_
    # bins: matplotlib hist with 9 equal bins between the outer vertical edges (line centres = first/last
    # column above the baseline -/+ half the line width)
    ok = np.isfinite(top)
    up = ok & (level > 0.5)
    c0, c1 = cols[up][0] + lw_half, cols[up][-1] - lw_half
    v0, v1 = ax_ * c0 + bx_, ax_ * c1 + bx_
    edges = np.linspace(v0, v1, 10)
    pct = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = ok & (vel > a + 0.6) & (vel < b - 0.6) & ~((vel > 8.4) & (vel < 9.1))   # skip the dashed reference line
        pct.append(float(np.median(level[m])))
    pct = np.array(pct)
    # bin heights are integer multiples of 100/N
    best = min(range(40, 300), key=lambda n: np.sum((pct * n / 100 - np.round(pct * n / 100)) ** 2) + 1e-4 * n)
    counts = np.round(pct * best / 100).astype(int)
    return dict(edges_kms=[round(float(e), 2) for e in edges], percent=[round(float(p), 2) for p in pct],
                N_inferred=int(best), counts=counts.tolist(),
                px_per_kms=round(float(1 / ax_), 3), px_per_percent=round(float(-1 / ay_), 3))


def main():
    res = {'about': 'Section 1 of the Comment on Wu+2026: Table 4 distance units, re-derived from the paper\'s own numbers',
           'lb_deg': {k: [round(v[0], 3), round(v[1], 3)] for k, v in LB.items()},
           'models': list(MODELS)}
    kdref = json.loads((ROOT / 'data/derived/verify_paper_kd.json').read_text())

    # ---- 1) rotation velocity at small distances; distance needed for the [O III] velocities
    near = {}
    for k in ['P1', 'P2', 'P3', 'P4']:
        row = {}
        for m in MODELS:
            v = vrot(k, np.array([0.0, 0.16, 0.30, 1.0, 1.6]), m)
            row[m] = dict(V_0=round(float(v[0]), 2), V_160pc=round(float(v[1]), 2), V_300pc=round(float(v[2]), 2),
                          V_1kpc=round(float(v[3]), 2), V_1p6kpc=round(float(v[4]), 2))
        near[k] = row
    res['V_rot_near'] = near
    v160 = [near[k][MAIN]['V_160pc'] for k in ['P1', 'P2', 'P3']]
    v160_all = [near[k][m]['V_160pc'] - near[k][m]['V_0'] for k in ['P1', 'P2', 'P3'] for m in MODELS]
    v300_max = max(near[k][MAIN]['V_300pc'] for k in ['P1', 'P2', 'P3'])
    res['summary_near'] = dict(V160_main_P1P2P3=v160, V160_minus_V0_all_models=[round(min(v160_all), 2), round(max(v160_all), 2)],
                               V300_main_max_P1P2P3=v300_max)

    # ---- 2) local rotation gradient dV/dd along P1-P4, 0.05-2.5 kpc, all models; Oort first-order slope
    grad = {}
    gmin, gmax = np.inf, -np.inf
    sel = (D >= 0.05) & (D <= 2.5)
    for m in MODELS:
        grad[m] = {}
        for k in ['P1', 'P2', 'P3', 'P4']:
            g = np.gradient(vrot(k, D, m), D)[sel]
            grad[m][k] = [round(float(g.min()), 2), round(float(g.max()), 2)]
            gmin, gmax = min(gmin, g.min()), max(gmax, g.max())
    oort = {k: round(A_OORT * np.sin(np.radians(2 * LB[k][0])) * np.cos(np.radians(LB[k][1])) ** 2, 2) for k in LB}
    res['gradient_kms_per_kpc'] = dict(per_model=grad, band_all=[round(float(gmin), 2), round(float(gmax), 2)],
                                       band_main=[round(min(v[0] for v in grad[MAIN].values()), 2),
                                                  round(max(v[1] for v in grad[MAIN].values()), 2)],
                                       oort_A_sin2l_cos2b=oort, oort_note='A = 15.3 +- 0.4 km/s/kpc (Bovy 2017)')
    G_LO = res['gradient_kms_per_kpc']['band_all'][0]
    G_HI = max(res['gradient_kms_per_kpc']['band_all'][1], max(oort.values()))
    res['gradient_kms_per_kpc']['band_used'] = [G_LO, G_HI]

    # ---- 3) row-by-row test of the two unit hypotheses
    rows = []
    zk = {'OIII': [], 'Ha': [], 'red': []}
    z1 = {'OIII': [], 'Ha': [], 'red': []}
    for k, lines in vk.T4.items():
        for ln in LINES:
            if ln not in lines:
                continue
            vh, ve, vl, mcm, mcn, mcp, tr = lines[ln]
            l, b = LB[k]
            vl_h = float(vk.helio_to_lsr(l, b, vh))
            d_rec = solve_d(k, vl)
            d_rec_h = solve_d(k, vl_h)
            d_models = [solve_d(k, vl, m) for m in MODELS]
            d_models = [x for x in d_models if x is not None]
            v_kpc, v_100 = float(vrot(k, tr)), float(vrot(k, tr / 10))
            zkp, z10 = (vl - v_kpc) / ve, (vl - v_100) / ve
            grp = ln if ln in ('OIII', 'Ha') else 'red'
            zk[grp].append(zkp)
            z1[grp].append(z10)
            # implied gradient from the printed MC interval (velocity-only minimum width)
            vfloor = float(vrot(k, 0.01))                      # V at kd's 0.01-kpc grid floor
            lo, hi = hpd_truncnorm(vl, ve, vfloor)
            wv = hi - lo
            wd = mcp + mcn
            g_kpc, g_100 = wv / wd, 10 * wv / wd
            # cross-check with the kd.pdf_kd emulation in verify_paper_kd.json
            em = kdref['monte_carlo'][k][ln]['mine']
            w_em = em[2] - em[1]                                 # emulated 68.3 % HPD width (kpc), tabulated sigma_v
            w_ratio_kpc, w_ratio_100 = wd / w_em, (wd / 10) / w_em
            rows.append(dict(region=k, line=ln, V_helio=vh, sigma_v=ve, V_LSR=vl, V_LSR_from_Vhelio=round(vl_h, 2),
                             MC_mode=mcm, MC_minus=mcn, MC_plus=mcp, MC_lower=round(mcm - mcn, 2), trad=tr,
                             d_recomputed_kpc=None if d_rec is None else round(d_rec, 3),
                             d_recomputed_from_Vhelio_kpc=None if d_rec_h is None else round(d_rec_h, 3),
                             d_recomputed_model_range=None if not d_models else [round(min(d_models), 3), round(max(d_models), 3)],
                             V_rot_at_trad_as_kpc=round(v_kpc, 2), V_rot_at_trad_as_100pc=round(v_100, 2),
                             z_kpc=round(zkp, 2), z_100pc=round(z10, 2),
                             W_V_kms=round(wv, 2), W_d_printed=round(wd, 2),
                             g_implied_kpc=round(g_kpc, 2), g_implied_100pc=round(g_100, 1),
                             emulated_MC_mode_HPD_kpc=em, emulated_HPD_width_kpc=round(w_em, 3),
                             printed_over_emulated_width_kpc=round(w_ratio_kpc, 2),
                             printed_over_emulated_width_100pc=round(w_ratio_100, 3),
                             trad_as_100pc_pc=round(tr * 100, 0), MC_mode_as_100pc_pc=round(mcm * 100, 0)))
    res['rows'] = rows

    def chi(zz, n):
        c = float(np.sum(np.square(zz)))
        return dict(chi2=round(c, 2), dof=n, p=float(chi2dist.sf(c, n)))
    tests = {}
    for grp, n in [('OIII', 3), ('Ha', 4), ('red', 8)]:
        tests[grp] = dict(kpc=chi(zk[grp], n), x100pc=chi(z1[grp], n))
    allk = sum((zk[g] for g in zk), [])
    all1 = sum((z1[g] for g in z1), [])
    tests['all15'] = dict(kpc=chi(allk, 15), x100pc=chi(all1, 15))
    res['unit_tests_velocity_space'] = tests
    # consistency with verify_paper_kd.json
    for r in rows:
        ref = kdref['traditional'][r['region']][r['line']]
        mine_ref = ref['my_trad_from_table_V']['URC_kd']
        if mine_ref is not None and r['d_recomputed_kpc'] is not None:
            assert abs(mine_ref - r['d_recomputed_kpc']) < 0.003, (r, mine_ref)
    oiii_t = [r for r in rows if r['line'] == 'OIII']
    res['OIII_trad_printed_vs_recomputed'] = [(r['region'], r['trad'], r['d_recomputed_kpc'], r['d_recomputed_from_Vhelio_kpc'],
                                              r['d_recomputed_model_range']) for r in oiii_t]
    res['OIII_trad_per_model_kpc'] = {r['region']: {m: round(solve_d(r['region'], r['V_LSR'], m), 3) for m in MODELS}
                                      for r in oiii_t}

    # ---- 4) implied-gradient summary and floor
    gk = {r['region'] + ' ' + r['line']: r['g_implied_kpc'] for r in rows}
    g1 = {r['region'] + ' ' + r['line']: r['g_implied_100pc'] for r in rows}
    res['implied_gradient_summary'] = dict(
        physical_band=[G_LO, G_HI],
        kpc_range=[min(gk.values()), max(gk.values())], x100pc_range=[min(g1.values()), max(g1.values())],
        OIII_kpc=[gk[f'P{i} OIII'] for i in (1, 2, 3)], OIII_100pc=[g1[f'P{i} OIII'] for i in (1, 2, 3)],
        n_rows_kpc_below_band_max=int(sum(v <= G_HI for v in gk.values())),
        n_rows_100pc_above_band_max=int(sum(v > G_HI for v in g1.values())), n_rows=len(rows),
        OIII_100pc_over_band_max=[round(g1[f'P{i} OIII'] / G_HI, 2) for i in (1, 2, 3)],
        Ha_100pc_over_band_max=[round(g1[f'P{i} Ha'] / G_HI, 2) for i in (1, 2, 3, 4)],
        OIII_printed_over_emulated_width_kpc=[r['printed_over_emulated_width_kpc'] for r in rows if r['line'] == 'OIII'],
        OIII_emulated_over_printed_width_100pc=[round(1 / r['printed_over_emulated_width_100pc'], 1)
                                                for r in rows if r['line'] == 'OIII'],
        width_note='emulated = kd.pdf_kd emulation (verify_paper_kd.json) with the tabulated sigma_v; printed width = '
                   'MC_plus + MC_minus (kpc reading) or /10 (x100 pc reading)',
        logic='Rotation-curve resampling, KDE smoothing and any extra error terms can only widen the MC interval, '
              'so W_V / W_d <= true dV/dd if the MC used the tabulated velocity errors or larger.')
    red = [r for r in rows if r['line'] != 'OIII']
    sn = [r for r in rows if r['line'] in ('SII', 'NII')]
    res['red_line_MC'] = dict(lower_bounds=[r['MC_lower'] for r in red],
                              emulated_lower_edge_SII_NII_kpc=[min(r['emulated_MC_mode_HPD_kpc'][1] for r in sn),
                                                               max(r['emulated_MC_mode_HPD_kpc'][1] for r in sn)],
                              modes_as_100pc_pc=[min(r['MC_mode_as_100pc_pc'] for r in red), max(r['MC_mode_as_100pc_pc'] for r in red)],
                              trad_as_100pc_pc=[min(r['trad_as_100pc_pc'] for r in red), max(r['trad_as_100pc_pc'] for r in red)],
                              kd_dist_min_kpc=0.01, note='lower bound = mode - err; all 12 red-line rows (Hα, [S II], [N II])')

    # ---- 5) fiducial = median of six [O III] values
    six = [r['MC_mode'] for r in oiii_t] + [r['trad'] for r in oiii_t]
    res['fiducial'] = dict(six_values=six, median=float(np.median(six)), as_100pc_pc=round(float(np.median(six)) * 100, 1),
                           as_kpc_pc=round(float(np.median(six)) * 1000, 0),
                           OIII_range_units=[min(six), max(six)], paper_text_range_pc='~100-190 pc (Sect. 3.3)')

    # ---- 6) paper-internal: "V_LSR differences of 10-16 km/s ... distance discrepancy of 100-150 pc" (Sect. 4.2)
    res['paper_internal_gradient'] = dict(dV=[10, 16], dd_pc=[100, 150],
                                          implied_kms_per_kpc=[round(10 / 0.15, 0), round(16 / 0.10, 0)],
                                          matched_ends=[round(10 / 0.10, 0), round(16 / 0.15, 0)],
                                          in_kpc=[round(10 / 1.5, 1), round(16 / 1.0, 1)])
    dif = []
    for k in ['P1', 'P2', 'P3']:
        for col in (6, 3):
            o = vk.T4[k]['OIII'][col]
            for ln in ('SII', 'NII', 'Ha'):
                dif.append(o - vk.T4[k][ln][col])
    res['paper_internal_gradient']['OIII_minus_red_distance_units'] = [round(min(dif), 2), round(max(dif), 2)]

    # ---- 7) downstream scalings
    d_old, d_new = 0.160, float(np.median(six))
    res['downstream'] = dict(
        length_10deg_pc={'160 pc': round(2 * 160 * np.tan(np.radians(5)), 1), f'{d_new:.3f} kpc': round(2 * d_new * 1000 * np.tan(np.radians(5)), 0)},
        width_1deg_pc={'160 pc': round(160 * np.radians(1), 2), f'{d_new:.3f} kpc': round(d_new * 1000 * np.radians(1), 1)},
        Rs_30pc_scaled_pc=round(30 * d_new / d_old, 0),
        z_kpc_at_new=round(d_new * np.sin(np.radians(LB['P1'][1])), 3),
        z_pc_at_160pc=round(160 * np.sin(np.radians(LB['P1'][1])), 1),
        note='sizes scale as d; energetics handled in Section 6')

    # ---- 8) the paper's own Fig. 5b [O III] V_LSR histogram
    f5 = digitize_fig5b()
    e = np.array(f5['edges_kms'])
    c = np.array(f5['counts'])
    f5['n_first_bin'] = int(c[0])
    f5['frac_above_first_bin'] = round(float(c[1:].sum() / c.sum()), 3)
    f5['frac_above_second_bin'] = round(float(c[2:].sum() / c.sum()), 3)
    p = np.array(f5['percent'])
    f5['raw_percent_sum'] = round(float(p.sum()), 2)
    f5['raw_frac_above_first_bin'] = round(float(p[1:].sum() / p.sum()), 3)
    f5['raw_frac_above_second_bin'] = round(float(p[2:].sum() / p.sum()), 3)
    f5['N_note'] = ('N is not printed in the paper; bin heights are integer multiples of 100/81 % '
                    '(N = 81 or a multiple). Fractions do not depend on N.')
    f5['binned_mean_kms'] = round(float(np.sum(0.5 * (e[1:] + e[:-1]) * c) / c.sum()), 2)
    f5['d_at_edges_P1_kpc'] = [None if x is None else round(x, 2) for x in (solve_d('P1', v) for v in e)]
    res['fig5b_OIII'] = f5

    OUT_JSON.write_text(json.dumps(res, indent=1, default=lambda o: o.item() if isinstance(o, np.generic) else str(o)))
    print(json.dumps({k: res[k] for k in ['summary_near', 'unit_tests_velocity_space', 'implied_gradient_summary', 'red_line_MC',
                                          'fiducial', 'paper_internal_gradient', 'downstream', 'fig5b_OIII', 'OIII_trad_printed_vs_recomputed']},
                     indent=1, default=str))
    print('gradient band all', res['gradient_kms_per_kpc']['band_all'], 'main', res['gradient_kms_per_kpc']['band_main'],
          'oort', oort)
    for r in rows:
        print(f"{r['region']} {r['line']:4s} V={r['V_LSR']:5.1f}±{r['sigma_v']:.1f} MC {r['MC_mode']:.2f}+{r['MC_plus']:.2f}-{r['MC_minus']:.2f} "
              f"tr {r['trad']:.2f} rec {r['d_recomputed_kpc']} (Vh {r['d_recomputed_from_Vhelio_kpc']}) z {r['z_kpc']:+.2f}/{r['z_100pc']:+.2f} "
              f"W_V {r['W_V_kms']:.2f} g {r['g_implied_kpc']:.2f}/{r['g_implied_100pc']:.1f} em {r['emulated_MC_mode_HPD_kpc']}")
    make_figure(res)


# ------------------------------------------------------------------------------------------------ figure
def make_figure(res):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from silk_style import new_sky_figure, finish, label, FG, FGDIM, OK, WARN, DANGER, BG2

    REGC = {'P1': '#4fc3f7', 'P2': '#7ad97a', 'P3': '#c792ea', 'P4': '#ffb84f'}
    MK = {'OIII': 'o', 'Ha': 's', 'SII': '^', 'NII': 'v'}
    rows = res['rows']
    tests = res['unit_tests_velocity_space']

    fig, ax0 = new_sky_figure(13.2, 10.6, 'Wu+2026 Table 4: the distances are in kpc',
                              "Table 4 and Fig. 5b of the paper vs the Reid+14 rotation curve it names "
                              "(kd form), P1–P4 sightlines")
    ax0.remove()
    gs = fig.add_gridspec(2, 12, left=0.06, right=0.948, top=0.905, bottom=0.075, hspace=0.30, wspace=1.7)
    axa = fig.add_subplot(gs[0, :8])
    axb = fig.add_subplot(gs[0, 8:])
    axc = fig.add_subplot(gs[1, :5])
    axd = fig.add_subplot(gs[1, 5:])

    def ptitle(ax, letter, text):
        ax.text(0.0, 1.015, letter, transform=ax.transAxes, color=FG, fontsize=12, fontweight='bold', va='bottom')
        ax.text(0.04 if ax is not axb else 0.075, 1.015, text, transform=ax.transAxes, color='#cfd6dd', fontsize=9.5, va='bottom')

    # ---------------- (a) V_LSR(d) along P1-P3 with Table 4 values placed under both readings
    ax = axa
    d = np.linspace(0, 2.6, 800)
    env = np.array([[vrot(k, d, m) for m in MODELS] for k in ['P1', 'P2', 'P3']]).reshape(-1, d.size)
    ax.fill_between(d, env.min(0), env.max(0), color='#9aa4ae', alpha=0.13, lw=0)
    for k, ls in [('P1', '-'), ('P2', '--'), ('P3', '-.')]:
        ax.plot(d, vrot(k, d), color=REGC[k], lw=1.6, ls=ls)
    for r in rows:
        if r['region'] == 'P4' or r['line'] not in ('OIII', 'Ha'):
            continue
        c, mk = REGC[r['region']], MK[r['line']]
        ax.errorbar(r['trad'], r['V_LSR'], yerr=r['sigma_v'], fmt=mk, color=c, ms=7, mec='#05070a', mew=0.6,
                    capsize=0, elinewidth=1.1, zorder=5)
        x1 = r['trad'] / 10
        ax.errorbar(x1, r['V_LSR'], yerr=r['sigma_v'], fmt=mk, mfc='none', mec=c, color=c, ms=7, mew=1.3,
                    capsize=0, elinewidth=0.9, alpha=0.95, zorder=5)
        ax.plot([x1, x1], [float(vrot(r['region'], x1)), r['V_LSR'] - r['sigma_v']], color=c, lw=0.8, ls=':', zorder=3)
    ax.set_xlim(0, 2.6)
    ax.set_ylim(-3, 24.5)
    ax.axhline(0, color='#2a323c', lw=0.8)
    ax.set_xlabel('heliocentric distance (kpc)')
    ax.set_ylabel('V_LSR (km/s)')
    ax.grid(alpha=0.25)
    t = tests['OIII']
    gaps = [r['V_LSR'] - r['V_rot_at_trad_as_100pc'] for r in rows if r['line'] == 'OIII']
    label(ax, 0.245, 21.4, f"open: Table 4 read as ×100 pc\n[O III] {min(gaps):.0f}–{max(gaps):.0f} km/s above rotation\n"
          f"χ² = {t['x100pc']['chi2']:.1f} (3 dof), p = {t['x100pc']['p']:.0e}".replace('e-0', 'e-'), color=WARN, size=8.3, va='top')
    label(ax, 1.93, 12.6, f"filled: read as kpc\non the curves\n[O III] χ² = {t['kpc']['chi2']:.2f} (3 dof)", color=OK, size=8.3, va='top')
    hand = [Line2D([], [], color=REGC['P1'], ls='-', lw=1.6, label='P1 l 203.8 b 19.5'),
            Line2D([], [], color=REGC['P2'], ls='--', lw=1.6, label='P2 l 202.5 b 19.6'),
            Line2D([], [], color=REGC['P3'], ls='-.', lw=1.6, label='P3 l 205.5 b 18.5'),
            Line2D([], [], color='#9aa4ae', lw=6, alpha=0.3, label='7 curve variants'),
            Line2D([], [], color=FG, marker='o', ls='none', ms=6, label='[O III] trad.'),
            Line2D([], [], color=FG, marker='s', ls='none', ms=6, label='Hα trad.')]
    ax.legend(handles=hand, loc='upper left', bbox_to_anchor=(0.375, 0.995), ncol=2, fontsize=7.6, framealpha=0.9)
    # inset 0-300 pc
    ins = ax.inset_axes([0.585, 0.10, 0.39, 0.36])
    di = np.linspace(0, 0.3, 300)
    for k, ls in [('P1', '-'), ('P2', '--'), ('P3', '-.')]:
        ins.plot(di * 1000, vrot(k, di), color=REGC[k], lw=1.4, ls=ls)
    ins.axvline(160, color=WARN, lw=1.0, ls=':')
    v16 = res['summary_near']['V160_main_P1P2P3']
    ins.set_xlim(0, 300)
    ins.set_ylim(0, 4.2)
    ins.set_facecolor('#0b0f14')
    ins.tick_params(labelsize=7)
    ins.set_xlabel('d (pc)', fontsize=7.5, labelpad=1)
    ins.set_ylabel('V_rot (km/s)', fontsize=7.5, labelpad=1)
    ins.grid(alpha=0.25)
    label(ins, 168, 0.35, f"160 pc: {min(v16):.1f}–{max(v16):.1f} km/s", color=WARN, size=7.5)
    label(ins, 8, 3.55, '≤ 300 pc: ≤ %.1f km/s' % res['summary_near']['V300_main_max_P1P2P3'], color=FG, size=7.5)
    ins.annotate('', xy=(160, 4.15), xytext=(160, 3.35), arrowprops=dict(arrowstyle='->', color=WARN, lw=1.0))
    label(ins, 168, 3.0, '[O III] 11–16\noff scale', color=WARN, size=7, va='bottom')
    ptitle(ax, 'a', 'Rotation V_LSR(d) vs Table 4 [O III] and Hα')

    # ---------------- (b) the paper's Fig. 5b [O III] histogram on the same velocity scale
    ax = axb
    f5 = res['fig5b_OIII']
    e = np.array(f5['edges_kms'])
    pct = np.array(f5['percent'])
    ax.barh(0.5 * (e[1:] + e[:-1]), pct, height=np.diff(e) * 0.96, color='#4fc3f7', alpha=0.55, edgecolor='#4fc3f7', lw=0.8)
    vmax300 = res['summary_near']['V300_main_max_P1P2P3']
    ax.axhspan(-5, vmax300, color=WARN, alpha=0.18, lw=0)
    ax.axhline(np.mean(v16), color=WARN, lw=1.0, ls=':')
    ax.set_ylim(-5, 39)
    ax.set_xlim(0, 26)
    ax.set_xlabel('share of [O III] fibres (%)')
    ax.set_ylabel('[O III] V_LSR (km/s)')
    ax.grid(alpha=0.25)
    label(ax, 10.5, -3.8, 'rotation at d ≤ 300 pc', color=WARN, size=7.8)
    label(ax, 5.0, 33.0, f"{100 * f5['frac_above_first_bin']:.1f} % above {e[1]:.1f} km/s\n"
          f"{100 * f5['frac_above_second_bin']:.0f} % above {e[2]:.1f} km/s", color=FG, size=8.3)
    # secondary axis: P1 kinematic distance for each velocity
    vg = vrot('P1', D)

    s0 = (vg[1] - vg[0]) / (D[1] - D[0])

    def v2d(v):
        v = np.asarray(v, float)
        return np.where(v < vg[0], v / s0, np.interp(v, vg, D))

    def d2v(x):
        x = np.asarray(x, float)
        return np.where(x < 0, x * s0, np.interp(x, D, vg))
    sec = ax.secondary_yaxis('right', functions=(v2d, d2v))
    sec.set_yticks([0.16, 0.5, 1, 2, 3, 4])
    sec.set_yticklabels(['0.16', '0.5', '1', '2', '3', '4'])
    sec.tick_params(labelsize=7.5, colors=FGDIM)
    sec.set_ylabel('P1 kinematic distance (kpc)', fontsize=8.5, color=FGDIM)
    ptitle(ax, 'b', "Paper's Fig. 5b [O III]")

    # ---------------- (c) velocity-space parity: V_rot at the printed distance vs printed V_LSR
    ax = axc
    ax.plot([-6, 20], [-6, 20], color=FGDIM, lw=1.0, zorder=1)
    for r in rows:
        c, mk = REGC[r['region']], MK[r['line']]
        dim = r['line'] in ('SII', 'NII')
        al = 0.55 if dim else 1.0
        ax.errorbar(r['V_LSR'], r['V_rot_at_trad_as_kpc'], xerr=r['sigma_v'], fmt=mk, color=c, ms=6.5, mec='#05070a', mew=0.5,
                    capsize=0, elinewidth=0.9, alpha=al, zorder=4)
        ax.errorbar(r['V_LSR'], r['V_rot_at_trad_as_100pc'], xerr=r['sigma_v'], fmt=mk, mfc='none', mec=c, color=c, ms=6.5,
                    mew=1.2, capsize=0, elinewidth=0.6, alpha=al * 0.9, zorder=3)
    ax.set_xlim(-5, 20)
    ax.set_ylim(-2, 18)
    ax.set_xlabel('Table 4 V_LSR (km/s)')
    ax.set_ylabel('rotation V_LSR at the Table 4 distance (km/s)')
    ax.grid(alpha=0.25)
    ta, th, tr = tests['all15'], tests['Ha'], tests['red']
    label(ax, -4.4, 16.9, f"filled = kpc: all 15 rows χ² = {ta['kpc']['chi2']:.1f}, p = {ta['kpc']['p']:.2f}", color=OK, size=7.8, va='top')
    label(ax, -4.4, 15.6, f"open = ×100 pc: χ² = {ta['x100pc']['chi2']:.0f}, p = {ta['x100pc']['p']:.0e}".replace('e-0', 'e-'),
          color=WARN, size=7.8, va='top')
    label(ax, 11.0, -1.55, 'open: all ≤ 1.9 km/s', color=WARN, size=7.6)
    label(ax, 7.3, 8.5, '1:1', color=FGDIM, size=7.6, rotation=44)
    label(ax, -4.4, 8.4, '[S II]/[N II] (faded):\nV ≈ 0, not\ndiscriminating', color=FGDIM, size=7.2, va='top')
    hand = [Line2D([], [], color=REGC[k], marker='o', ls='none', ms=6, label=k) for k in ['P1', 'P2', 'P3', 'P4']]
    hand += [Line2D([], [], color=FG, marker=MK[ln], ls='none', ms=6, mfc='none' if False else FG, label=LNAME[ln]) for ln in LINES]
    ax.legend(handles=hand, loc='lower right', bbox_to_anchor=(1.0, 0.25), ncol=2, fontsize=7.2, framealpha=0.9,
              handletextpad=0.3, columnspacing=0.8)
    ptitle(ax, 'c', 'Parity in velocity, both readings')

    # ---------------- (d) implied gradient from the printed MC widths
    ax = axd
    glo, ghi = res['implied_gradient_summary']['physical_band']
    ax.axhspan(glo, ghi, color=OK, alpha=0.20, lw=0)
    ax.axhspan(ghi, 1000, color=DANGER, alpha=0.09, lw=0)
    x = np.arange(len(rows))
    for i, r in enumerate(rows):
        c, mk = REGC[r['region']], MK[r['line']]
        ax.plot(i, r['g_implied_kpc'], mk, color=c, ms=7, mec='#05070a', mew=0.5, zorder=4)
        ax.plot(i, r['g_implied_100pc'], mk, mfc='none', mec=c, ms=7, mew=1.3, zorder=4)
        ax.plot([i, i], [r['g_implied_kpc'], r['g_implied_100pc']], color=c, lw=0.6, alpha=0.5, zorder=2)
    ax.set_yscale('log')
    ax.set_ylim(0.6, 330)
    ax.set_xlim(-0.7, len(rows) - 0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([LNAME[r['line']] for r in rows], rotation=60, fontsize=7.3, ha='right', rotation_mode='anchor')
    for j, k in enumerate(['P1', 'P2', 'P3', 'P4']):
        idx = [i for i, r in enumerate(rows) if r['region'] == k]
        ax.text(np.mean(idx), 0.64, k, color=REGC[k], fontsize=8.5, ha='center', va='bottom', fontweight='bold')
        if j:
            ax.axvline(idx[0] - 0.5, color='#2a323c', lw=0.8)
    ax.set_ylabel('implied dV/dd = W_V / W_d (km/s/kpc)')
    ax.grid(alpha=0.25, axis='y')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:g}'))
    from matplotlib.patches import Patch
    hand = [Patch(facecolor=OK, alpha=0.35, edgecolor=OK, lw=0.8, label=f'rotation dV/dd, {glo:.1f}–{ghi:.1f}'),
            Patch(facecolor=DANGER, alpha=0.30, edgecolor=DANGER, lw=0.8, label='forbidden: narrower than σ_v allows'),
            Line2D([], [], color=FG, marker='o', ls='none', ms=6, mfc='none', mew=1.3, label='Table 4 read as ×100 pc'),
            Line2D([], [], color=FG, marker='o', ls='none', ms=6, label='Table 4 read as kpc')]
    ax.legend(handles=hand, loc='upper center', ncol=2, fontsize=7.6, framealpha=0.92, handletextpad=0.5,
              columnspacing=1.4, borderpad=0.5)
    ptitle(ax, 'd', 'Gradient implied by the Monte-Carlo widths')

    finish(fig, str(OUT_FIG), dpi=200, tight=False)
    print('wrote', OUT_FIG)


if __name__ == '__main__':
    main()
