"""dispute_kinematic_kd.py -- Section 2 of the dispute comment (tag: kinematic).

Why no kinematic distance is possible toward the Silk Nebula.

Re-uses (imports, does not modify) the rotation-curve code of scripts/verify_paper_kd.py:
  Persic+96 universal rotation curve with the Reid+14 A5 parameters as tabulated in the `kd` package
  (a1=241, a2=0.90, a3=1.46, R0=8.34 kpc), the `kd` projection V = R0 sin l (Theta/R - Theta0/R0),
  and the IAU standard solar motion used by `kd` (10.27, 15.32, 7.74) km/s.

Computes
  1. gradient map dV_LSR/dd over l 180-230, b 0-30 (mean slope 0-2 kpc; also local slope at d->0 and 1 kpc),
     values at P1-P4 and along the Silk spine (data/derived/silk_trace.json), position of the spine relative
     to the paper's 20-deg Galactic-anticentre (GAC) exclusion zone (both readings: |l-180|<20 and a 20-deg radius);
  2. distance shift produced by sigma_pec = 3-10 km/s;
  3. P(d | V_LSR) for the paper's P1 [O III] V_LSR = +15.5 +- 4.3 (Table 4) and the corrected
     [O III] V_LSR = +1.2 +- 1.0 (verify_velocity_summary.json), sigma_pec = 0, 3, 7, 10 km/s:
       (A) Bayesian posterior, prior uniform in d on [0, D_MAX], rotation-curve parameters marginalised;
       (B) kd-style Monte Carlo inversion (Wenger+18 method C) with a peculiar-velocity term added,
           samples below V(d=0) kept at d = 0 (kd itself drops them);
  4. formal kinematic distances of gas whose true distance is known to be < 15 pc (Redfield & Linsky 2008
     LIC and Gem clouds, kinematic-calculator output saved in data/raw/spectra_ism_lism/dynlism_P1.html)
     and of the northern [O III] (w3_northsouth: V_LSR +5.7...+10.4).

Outputs data/derived/dispute_kinematic_kd.json and figures/dispute_kinematic_fig2.png.
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_kinematic_kd.py
"""
import sys
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import verify_paper_kd as vk  # noqa: E402  (rotation curve, projection, solar motion, parameter errors)

rng = np.random.default_rng(20261117)

SIG_PEC = [0.0, 3.0, 7.0, 10.0]
D_MAX = 6.0                     # kpc, upper edge of the flat distance prior (z = 2.0 kpc at b = 19.5)
D_MAX_ALT = [4.0, 10.0]         # sensitivity of (A) to the prior edge

# ---------------- inputs ----------------
trace = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
REGIONS = {k: (v['l'], v['b']) for k, v in trace['regions'].items()}
spine = trace['spine']
SL = np.array([p['l'] for p in spine])
SB = np.array([p['b'] for p in spine])
SDEC = np.array([p['dec'] for p in spine])

vsum = json.loads((ROOT / 'data/derived/verify_velocity_summary.json').read_text())
V_CORR = round(vsum['combined']['V_LSR_standard'], 1)          # +1.2 (P1-local LSR correction)
V_CORR_ERR = round(vsum['combined']['err'], 1)                  # 1.0
V_PAPER, V_PAPER_ERR = vk.T4['P1']['OIII'][2], vk.T4['P1']['OIII'][1]   # 15.5, 4.3 (Table 4, P1 [O III])
PAPER_MC = vk.T4['P1']['OIII'][3:6]                             # 1.63 (+1.01/-0.84), Table 4 Distance^1
NORTH_VLSR = (5.7, 10.4)                                        # w3_northsouth.md Sect. 1 item 2 (DESI)
NORTH_DEC = 24.0                                                # representative northern spine point
PAPER_FIG1_DEC_MAX = 26.29     # Wu+26 Fig. 1 covers Dec 12.82-26.29 (notes/geometry.md)
CNC8 = dict(star='8 Cnc', d_pc=63.0, sep_deg=4.7, V_helio_CaII=[19.1, 19.4],
            source='notes/spectra_ism.md vetted-component table (HARPS, SOPHIE): LIC at 63 pc')


def lism_from_html(path):
    """Parse the Redfield & Linsky (2008) kinematic-calculator output (heliocentric v_r)."""
    txt = re.sub(r'<[^>]*>', ' ', Path(path).read_text(errors='replace'))
    out = {}
    for name in ['LIC', 'Gem', 'Aur', 'G']:
        m = re.search(r'\n\s*' + name + r'\s+(-?\d+\.\d+)\s+(\d+\.\d+)', txt)
        if m:
            out[name] = (float(m.group(1)), float(m.group(2)))
    trav = re.search(r'traverses the following clouds:(.*?)Sight line', txt, re.S)
    out['traverses'] = trav.group(1).split() if trav else []
    return out


# ---------------- helpers ----------------
def vrot(l, b, d, **kw):
    return vk.vlsr_model(l, b, d, **kw)


def grad_mean(l, b, **kw):
    return vrot(l, b, 2.0, **kw) / 2.0


def grad_local(l, b, d0, h=0.01, **kw):
    return (vrot(l, b, d0 + h, **kw) - vrot(l, b, max(d0 - h, 0.0), **kw)) / (d0 + h - max(d0 - h, 0.0))


def angsep(l1, b1, l2, b2):
    l1, b1, l2, b2 = map(np.radians, (l1, b1, l2, b2))
    c = np.sin(b1) * np.sin(b2) + np.cos(b1) * np.cos(b2) * np.cos(l1 - l2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def draw_theta(n):
    return {k: rng.normal(vk.R14[k], vk.R14_ERR[k], n) for k in ['a1', 'a2', 'a3', 'R0']}


def trad(l, b, v, **kw):
    """Traditional (single-valued, outer Galaxy) kinematic distance by bisection; NaN if v <= V(0)."""
    if v <= vrot(l, b, 0.0, **kw):
        return np.nan
    lo, hi = 0.0, 30.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if vrot(l, b, mid, **kw) < v:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mc_invert(l, b, v, verr, spec, n=40000):
    """(B) kd-style MC: resample rotation curve + V (+ peculiar term); invert V(d) by vectorised bisection."""
    th = draw_theta(n)
    vs = v + rng.normal(0, verr, n) + rng.normal(0, spec, n) if spec > 0 else v + rng.normal(0, verr, n)
    lo, hi = np.zeros(n), np.full(n, 30.0)
    v0 = vk.vlsr_model(l, b, 0.0, **th)
    vmax = vk.vlsr_model(l, b, 30.0, **th)
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        below = vk.vlsr_model(l, b, mid, **th) < vs
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    d = 0.5 * (lo + hi)
    at_zero = vs <= v0
    d[at_zero] = 0.0
    beyond = vs >= vmax
    return d, float(at_zero.mean()), float(beyond.mean())


def posterior(l, b, v, verr, spec, dmax=D_MAX, nth=800, step=0.0025):
    """(A) Bayesian P(d | V): flat prior on [0, dmax], Gaussian velocity likelihood with
    sigma^2 = verr^2 + spec^2, marginalised over Reid+14 parameter draws."""
    d = np.arange(0.0, dmax + step / 2, step)
    th = draw_theta(nth)
    sig = np.hypot(verr, spec)
    like = np.zeros_like(d)
    for j in range(nth):
        vm = vk.vlsr_model(l, b, d, a1=th['a1'][j], a2=th['a2'][j], a3=th['a3'][j], R0=th['R0'][j])
        like += np.exp(-0.5 * ((v - vm) / sig) ** 2)
    p = like / np.trapezoid(like, d)
    return d, p


def summarise_pdf(d, p):
    c = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(d))])
    c /= c[-1]
    q = lambda f: float(np.interp(f, c, d))
    return dict(median=round(q(0.5), 3), p16=round(q(0.16), 3), p84=round(q(0.84), 3),
                p2_5=round(q(0.025), 3), p97_5=round(q(0.975), 3), mode=round(float(d[np.argmax(p)]), 3),
                P_lt_0p3kpc=round(float(np.interp(0.3, d, c)), 3), P_gt_1kpc=round(1 - float(np.interp(1.0, d, c)), 3))


def summarise_samples(s):
    return dict(median=round(float(np.median(s)), 3), p16=round(float(np.percentile(s, 16)), 3),
                p84=round(float(np.percentile(s, 84)), 3), P_lt_0p3kpc=round(float(np.mean(s < 0.3)), 3),
                P_gt_1kpc=round(float(np.mean(s > 1.0)), 3))


# ---------------- main ----------------
def main():
    res = {'script': 'scripts/dispute_kinematic_kd.py', 'date': '2026-11-17',
           'rotation_curve': 'Persic+96 URC, Reid+14 A5 parameters as in kd (a1 241+-8, a2 0.90+-0.06, a3 1.46+-0.16, R0 8.34+-0.16)',
           'projection': "kd form V = R0 sin l (Theta/R - Theta0/R0) (no cos b); '3d' = x cos b",
           'solar_motion': 'IAU standard (10.27, 15.32, 7.74) km/s, as in kd and the paper'}

    # 1. gradient map
    lg = np.arange(180.0, 230.0 + 0.001, 0.2)
    bg = np.arange(0.0, 30.0 + 0.001, 0.2)
    LL, BB = np.meshgrid(lg, bg)
    G = grad_mean(LL, BB)
    res['gradient_map'] = dict(l_range=[180, 230], b_range=[0, 30], definition='V_rot(2 kpc)/2 kpc, kd form',
                               max=round(float(G.max()), 2),
                               argmax_lb=[round(float(LL.flat[np.argmax(G)]), 1), round(float(BB.flat[np.argmax(G)]), 1)])
    reg = {}
    for k, (l, b) in REGIONS.items():
        row = {}
        for name, kw in {'URC_kd': {}, 'URC_3d': dict(proj='3d'), 'flat_3d': dict(curve='flat', proj='3d')}.items():
            row[name] = dict(mean_0_2kpc=round(float(grad_mean(l, b, **kw)), 2),
                             local_d0=round(float(grad_local(l, b, 0.0, **kw)), 2),
                             local_1kpc=round(float(grad_local(l, b, 1.0, **kw)), 2),
                             V_0p16=round(float(vrot(l, b, 0.16, **kw)), 2), V_0p3=round(float(vrot(l, b, 0.3, **kw)), 2),
                             V_1p6=round(float(vrot(l, b, 1.6, **kw)), 2))
        reg[k] = dict(l=l, b=b, sep_from_GAC_deg=round(float(angsep(l, b, 180.0, 0.0)), 2),
                      dl_from_180=round(abs(l - 180.0), 2), models=row)
    res['regions'] = reg
    allg = [reg[k]['models'][m]['mean_0_2kpc'] for k in reg for m in reg[k]['models']]
    allg0 = [reg[k]['models'][m]['local_d0'] for k in reg for m in reg[k]['models']]
    allg1 = [reg[k]['models'][m]['local_1kpc'] for k in reg for m in reg[k]['models']]
    allv3 = [reg[k]['models'][m]['V_0p3'] for k in reg for m in reg[k]['models']]
    res['P1_P4_summary'] = dict(mean_0_2kpc=[min(allg), max(allg)], local_d0=[min(allg0), max(allg0)],
                                local_1kpc=[min(allg1), max(allg1)], V_rot_at_0p3kpc=[min(allv3), max(allv3)])
    gs = grad_mean(SL, SB)
    sep = angsep(SL, SB, 180.0, 0.0)
    south = SDEC < 19.5
    inside_band = np.abs(SL - 180.0) < 20.0
    # arc length of the spine inside the longitude band
    seg = angsep(SL[:-1], SB[:-1], SL[1:], SB[1:])
    segin = 0.5 * (inside_band[:-1].astype(float) + inside_band[1:].astype(float))
    i_cross = int(np.argmax(inside_band))
    res['spine'] = dict(n_points=int(SL.size), l_range=[round(float(SL.min()), 2), round(float(SL.max()), 2)],
                        b_range=[round(float(SB.min()), 2), round(float(SB.max()), 2)],
                        grad_mean_0_2kpc_kd=[round(float(gs.min()), 2), round(float(gs.max()), 2)],
                        grad_south_DecLT19p5=[round(float(gs[south].min()), 2), round(float(gs[south].max()), 2)],
                        grad_north_DecGE19p5=[round(float(gs[~south].min()), 2), round(float(gs[~south].max()), 2)],
                        V_rot_0p3kpc_max=round(float(vrot(SL, SB, 0.3).max()), 2),
                        sep_from_GAC_deg=[round(float(sep.min()), 2), round(float(sep.max()), 2)],
                        dl_from_180=[round(float(np.abs(SL - 180).min()), 2), round(float(np.abs(SL - 180).max()), 2)],
                        frac_points_inside_lband=round(float(inside_band.mean()), 3),
                        arc_inside_lband_deg=round(float((seg * segin).sum()), 2), arc_total_deg=round(float(seg.sum()), 2),
                        lband_crossing_dec=round(float(SDEC[i_cross]), 2),
                        arc_inside_lband_deg_within_paper_fig1=round(float(
                            (seg * segin * (0.5 * (SDEC[:-1] + SDEC[1:]) <= PAPER_FIG1_DEC_MAX)).sum()), 2),
                        paper_fig1_dec_max=PAPER_FIG1_DEC_MAX,
                        wenger18_mask='Wenger+18 Figs. 5, 8, 19 captions: "160 < l < 200 ... masked"',
                        note='"within 20 deg of the GAC" read as |l-180|<20 (Wenger+18 define GC/GAC by l) or as a 20-deg radius about (180,0)')

    # 2. sigma_pec -> distance shift
    gmin, gmax = min(allg), max(allg)
    res['dd_from_sigma_pec'] = {f'{s:.0f}': [round(s / gmax, 2), round(s / gmin, 2)] for s in [3, 5, 7, 10]}
    res['dd_note'] = 'Delta d = sigma_pec / (dV/dd) with the P1-P4 mean 0-2 kpc slope range (all three models)'

    # 3. posteriors
    l1, b1 = REGIONS['P1']
    cases = {'paper_P1_OIII': (V_PAPER, V_PAPER_ERR), 'corrected_OIII': (V_CORR, V_CORR_ERR)}
    post, mc, curves = {}, {}, {}
    for cname, (v, ve) in cases.items():
        post[cname], mc[cname] = {}, {}
        for s in SIG_PEC:
            d, p = posterior(l1, b1, v, ve, s)
            post[cname][f'{s:.0f}'] = summarise_pdf(d, p)
            post[cname][f'{s:.0f}']['prior_edge_sensitivity'] = {
                f'dmax_{dm:.0f}': {kk: vv for kk, vv in summarise_pdf(*posterior(l1, b1, v, ve, s, dmax=dm, nth=300)).items()
                                   if kk in ('median', 'p16', 'p84')} for dm in D_MAX_ALT}
            curves[f'{cname}_{s:.0f}'] = (d, p)
            smp, fz, fb = mc_invert(l1, b1, v, ve, s)
            mc[cname][f'{s:.0f}'] = summarise_samples(smp)
            mc[cname][f'{s:.0f}'].update(frac_at_d0=round(fz, 3), frac_beyond_30kpc=round(fb, 4))
        post[cname]['input'] = dict(V_LSR=v, err=ve, l=l1, b=b1)
    res['posterior_flat_prior'] = post
    res['posterior_note'] = (f'(A) flat prior on [0, {D_MAX:.0f}] kpc, Gaussian likelihood sigma^2 = err^2 + sigma_pec^2, '
                             'marginalised over 800 Reid+14 parameter draws; percentiles of the normalised posterior')
    res['mc_inversion'] = mc
    res['mc_note'] = '(B) 40000 samples; V(d) inverted by bisection; samples with V <= V_rot(0) kept at d = 0 (kd drops them)'
    res['paper_table4_P1_OIII'] = dict(V_LSR=V_PAPER, err=V_PAPER_ERR, MC_mode=PAPER_MC[0],
                                       MC_interval=[round(PAPER_MC[0] - PAPER_MC[1], 2), round(PAPER_MC[0] + PAPER_MC[2], 2)],
                                       traditional=vk.T4['P1']['OIII'][6])
    res['traditional_check'] = dict(paper_V=round(trad(l1, b1, V_PAPER), 3), corrected_V=round(trad(l1, b1, V_CORR), 3))

    # 4. gas of known (tiny) distance and the northern [O III]
    lism = lism_from_html(ROOT / 'data/raw/spectra_ism_lism/dynlism_P1.html')
    lres = {}
    for name in ['LIC', 'Gem']:
        vh, e = lism[name]
        vl = float(vk.helio_to_lsr(l1, b1, vh))
        lres[name] = dict(V_helio=vh, err=e, V_LSR=round(vl, 2), true_distance='< 15 pc (Redfield & Linsky 2008)',
                          formal_traditional_kd_kpc=round(trad(l1, b1, vl), 3),
                          V_rot_at_15pc=round(float(vrot(l1, b1, 0.015)), 3))
        lres[name]['formal_over_15pc'] = round(lres[name]['formal_traditional_kd_kpc'] / 0.015, 1)
    lres['traverses'] = lism['traverses']
    lres['LIC_observed'] = dict(CNC8, formal_LIC_kd_over_star_distance=round(
        lres['LIC']['formal_traditional_kd_kpc'] / (CNC8['d_pc'] / 1000.0), 1))
    res['lism_clouds_P1'] = lres
    i_n = int(np.argmin(np.abs(SDEC - NORTH_DEC)))
    ln_, bn_ = float(SL[i_n]), float(SB[i_n])
    res['north_OIII'] = dict(V_LSR_range=NORTH_VLSR, l=round(ln_, 3), b=round(bn_, 3), dec=float(SDEC[i_n]),
                             grad_mean_0_2kpc=round(float(grad_mean(ln_, bn_)), 2),
                             formal_traditional_kd_kpc=[round(trad(ln_, bn_, v), 3) for v in NORTH_VLSR],
                             V_rot_north_minus_P1_at=dict({f'{d}kpc': round(float(vrot(ln_, bn_, d) - vrot(l1, b1, d)), 2)
                                                           for d in [0.3, 1.0, 1.6]}),
                             source='w3_northsouth.md Sect. 1 item 2 (DESI; south V_LSR -2.8...-0.2)')
    res['south_OIII_DESI'] = dict(V_LSR_range=[-2.8, -0.2], V_rot_at_d0_P1=round(float(vrot(l1, b1, 0.0)), 4),
                                  formal_traditional_kd_kpc=[None if np.isnan(x) else round(x, 3)
                                                             for x in (trad(l1, b1, -2.8), trad(l1, b1, -0.2))],
                                  note='V_LSR below V_rot(0): no outer-Galaxy solution (formal d = 0)',
                                  source='w3_northsouth.md Sect. 1 item 2')

    out = ROOT / 'data/derived/dispute_kinematic_kd.json'
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ['P1_P4_summary', 'spine', 'dd_from_sigma_pec', 'traditional_check',
                                          'lism_clouds_P1', 'north_OIII', 'south_OIII_DESI', 'gradient_map']}, indent=1))
    for c in post:
        for s in SIG_PEC:
            a, m = post[c][f'{s:.0f}'], mc[c][f'{s:.0f}']
            print(c, s, 'post', a['median'], a['p16'], a['p84'], 'P<0.3', a['P_lt_0p3kpc'], 'P>1', a['P_gt_1kpc'],
                  'edge', a['prior_edge_sensitivity'], '| mc', m['median'], m['p16'], m['p84'], 'P<0.3', m['P_lt_0p3kpc'],
                  'at0', m['frac_at_d0'])
    make_fig(LL, BB, G, curves, res)


# ---------------- figure ----------------
def make_fig(LL, BB, G, curves, res):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from silk_style import new_sky_figure, finish, label, FG, FGDIM, WARN, NII_C, ACCENT, BG2, HALO

    RAMP = ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab']   # ordinal ramp, validated (dataviz validator, --ordinal, dark)
    LS = ['-', (0, (6, 2)), (0, (3, 1.5)), (0, (1, 1.2))]
    SPINE_C = '#ff8a50'

    fig, axs = new_sky_figure(13.0, 9.8, 'Kinematic distance toward the Silk Nebula',
                              'Reid+14 rotation curve (kd form) · IAU standard solar motion · σ_pec = cloud peculiar velocity',
                              nrows=2, ncols=2)
    # (a) gradient map
    ax = axs[0, 0]
    cmap = LinearSegmentedColormap.from_list('grad', ['#05070a', '#0d366b', '#256abf', '#5598e7', '#9ec5f4'])
    pm = ax.pcolormesh(LL, BB, G, cmap=cmap, vmin=0, vmax=14, shading='auto', rasterized=True)
    cs = ax.contour(LL, BB, G, levels=[2, 4, 6, 8, 10, 12], colors=['#cfd6dd'], linewidths=0.6, alpha=0.7)
    lab_b = 25.0   # place contour labels at b = 25 deg, clear of the spine, P1-P4 and the 20-deg circle
    ib = int(np.argmin(np.abs(BB[:, 0] - lab_b)))
    manual = [(float(np.interp(lev, G[ib, :], LL[ib, :])), lab_b) for lev in cs.levels]   # G rises with l here
    for tx in ax.clabel(cs, fmt='%d', fontsize=7.5, colors='#ffffff', inline_spacing=2, manual=manual):
        tx.set_path_effects([HALO])
    # exclusion zone: |l-180| < 20 band and 20-deg radius circle
    ax.axvspan(180, 200, color='#ffffff', alpha=0.07, lw=0)
    ax.axvline(200, color='#e8eaed', lw=1.0, ls=(0, (4, 2)))
    t = np.linspace(0, 2 * np.pi, 721)
    # small circle of radius 20 deg about (180, 0) in (l, b)
    lc0, bc0, r = np.radians(180.0), 0.0, np.radians(20.0)
    bb = np.arcsin(np.sin(bc0) * np.cos(r) + np.cos(bc0) * np.sin(r) * np.cos(t))
    ll = lc0 + np.arctan2(np.sin(t) * np.sin(r) * np.cos(bc0), np.cos(r) - np.sin(bc0) * np.sin(bb))
    ax.plot(np.degrees(ll), np.degrees(bb), color='#e8eaed', lw=1.0, ls=(0, (1, 1.5)))
    ax.plot(SL, SB, color=SPINE_C, lw=2.4, solid_capstyle='round')
    for k, (l, b) in REGIONS.items():
        ax.plot(l, b, 'o', ms=4.5, mfc='#ffffff', mec='#05070a', mew=0.8, zorder=5)
    label(ax, 211.0, 22.6, 'Silk spine', color=SPINE_C, size=8.5)
    label(ax, 207.6, 15.2, 'P1–P4', color='#ffffff', size=8)
    label(ax, 198.8, 28.4, '|l−180°| < 20°', color=FG, size=8, ha='left')
    label(ax, 195.5, 9.0, '20° from GAC', color=FG, size=8, ha='left')
    sp = res['spine']
    label(ax, 229.0, 1.2, f"spine: {sp['grad_mean_0_2kpc_kd'][0]:.1f}–{sp['grad_mean_0_2kpc_kd'][1]:.1f} km/s/kpc",
          color=FG, size=8)
    ax.set_xlim(230, 180)
    ax.set_ylim(0, 30)
    ax.set_xlabel('Galactic longitude l (deg)')
    ax.set_ylabel('Galactic latitude b (deg)')
    cb = fig.colorbar(pm, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('km/s/kpc', color=FG)
    cb.ax.tick_params(colors=FGDIM)
    ax.set_title('(a) Rotation gradient dV_LSR/dd, mean over 0–2 kpc', loc='left', color=FG)

    # (b) V_LSR vs distance
    ax = axs[0, 1]
    l1, b1 = REGIONS['P1']
    nn = res['north_OIII']
    d = np.linspace(0, 3.0, 400)
    v1 = vrot(l1, b1, d)
    vn = vrot(nn['l'], nn['b'], d)
    ax.fill_between(d, v1 - 7, v1 + 7, color='#9aa4ae', alpha=0.16, lw=0, label='P1 rotation ± 7 km/s')
    ax.plot(d, v1, color=FG, lw=1.6, label='rotation, P1 (l 203.8°, b 19.5°)')
    ax.plot(d, vn, color=FG, lw=1.2, ls=(0, (4, 2)), label='rotation, north spine (Dec 24°)')
    ax.axhspan(V_PAPER - V_PAPER_ERR, V_PAPER + V_PAPER_ERR, color=WARN, alpha=0.20, lw=0)
    ax.axhline(V_PAPER, color=WARN, lw=1.2)
    ax.axhspan(V_CORR - V_CORR_ERR, V_CORR + V_CORR_ERR, color=ACCENT, alpha=0.30, lw=0)
    ax.axhline(V_CORR, color=ACCENT, lw=1.2)
    ax.axhspan(NORTH_VLSR[0], NORTH_VLSR[1], color=NII_C, alpha=0.18, lw=0)
    lis = res['lism_clouds_P1']
    for name in ['LIC', 'Gem']:
        vv, dk = lis[name]['V_LSR'], lis[name]['formal_traditional_kd_kpc']
        ax.plot([0.06, dk], [vv, vv], color='#ffffff', lw=0.9, ls=(0, (1, 1.6)), zorder=5)
        ax.plot(dk, vv, 'o', ms=5, mfc='none', mec='#ffffff', mew=1.1, zorder=6)
        ax.plot(0.03, vv, marker='D', ms=6.5, mfc='#ffffff', mec='#05070a', mew=0.8, zorder=6, clip_on=False)
    label(ax, 0.09, lis['LIC']['V_LSR'] + 1.5,
          f"LIC {lis['LIC']['V_LSR']:+.1f} → {lis['LIC']['formal_traditional_kd_kpc']:.2f} kpc", color='#ffffff', size=8)
    label(ax, 0.09, lis['Gem']['V_LSR'] - 1.9,
          f"Gem {lis['Gem']['V_LSR']:+.1f} → {lis['Gem']['formal_traditional_kd_kpc']:.2f} kpc  (both < 15 pc)",
          color='#ffffff', size=8)
    label(ax, 0.30, V_PAPER + 0.8, f'paper P1 [O III] {V_PAPER:+.1f} ± {V_PAPER_ERR}', color=WARN, size=8)
    label(ax, 0.30, V_CORR - 2.4, f'corrected [O III] {V_CORR:+.1f} ± {V_CORR_ERR}', color=ACCENT, size=8)
    label(ax, 1.95, NORTH_VLSR[1] - 2.3, 'north [O III] +5.7…+10.4', color=NII_C, size=8)
    ax.legend(loc='lower right', fontsize=7.5, frameon=True, facecolor=BG2)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(-8, 32)
    ax.axhline(0, color='#2a323c', lw=0.8)
    ax.set_xlabel('distance d (kpc)')
    ax.set_ylabel('V_LSR (km/s)')
    ax.grid(alpha=0.25)
    ax.set_title('(b) Circular-rotation velocity and measured velocities', loc='left', color=FG)

    # (c), (d) posteriors
    pm_ = res['paper_table4_P1_OIII']
    for ax, cname, ttl, vtxt in [(axs[1, 0], 'paper_P1_OIII', '(c) P(d | V_LSR = +15.5 ± 4.3), paper P1 [O III]', None),
                                 (axs[1, 1], 'corrected_OIII', '(d) P(d | V_LSR = +1.2 ± 1.0), corrected [O III]', None)]:
        post = res['posterior_flat_prior'][cname]
        for i, s in enumerate(SIG_PEC):
            dd, p = curves[f'{cname}_{s:.0f}']
            ax.plot(dd, p / p.max(), color=RAMP[i], ls=LS[i], lw=1.8, label=f'σ_pec = {s:.0f} km/s')
            q = post[f'{s:.0f}']
            y = -0.08 - 0.075 * i
            ax.plot([q['p16'], q['p84']], [y, y], color=RAMP[i], lw=3.2, solid_capstyle='butt')
            ax.plot(q['median'], y, '|', color='#ffffff', ms=7, mew=1.4)
            label(ax, q['p84'] + 0.07, y - 0.022, f'{s:.0f}', color=FGDIM, size=7.5)
        if cname == 'paper_P1_OIII':
            y = -0.08 - 0.075 * 4
            ax.plot(pm_['MC_interval'], [y, y], color=WARN, lw=3.2, solid_capstyle='butt')
            ax.plot(pm_['MC_mode'], y, '|', color='#ffffff', ms=7, mew=1.4)
            label(ax, pm_['MC_interval'][1] + 0.08, y - 0.025, 'paper Table 4 MC (kpc)', color=WARN, size=7.5)
        for xv, txt in [(0.16, '160 pc'), (1.6, '1.6 kpc')]:
            ax.axvline(xv, color=FGDIM, lw=0.8, ls=(0, (2, 2)))
            label(ax, xv + 0.05, 1.085, txt, color=FGDIM, size=7.5)
        if cname == 'corrected_OIII':
            label(ax, 4.95, -0.10, 'bars: 16–84 % (σ_pec), tick = median', color=FGDIM, size=7.5, ha='right')
        ax.set_xlim(0, 5.0)
        ax.set_ylim(-0.42, 1.16)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.axhline(0, color='#2a323c', lw=0.8)
        ax.set_xlabel('distance d (kpc)')
        ax.set_ylabel('posterior (peak = 1)')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper right', fontsize=8, frameon=True, facecolor=BG2)
        ax.set_title(ttl, loc='left', color=FG)
    finish(fig, str(ROOT / 'figures/dispute_kinematic_fig2.png'))


if __name__ == '__main__':
    main()
