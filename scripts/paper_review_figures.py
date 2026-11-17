#!/usr/bin/env python
"""paper_review figures for the critical review of Wu et al. 2026 (arXiv:2610.13329).
Reads data/derived/paper_checks.json and data/derived/paper_review_mc_samples.npz written by
scripts/paper_checks.py; writes figures/paper_review_*.png (dpi 200, silk_style).

Run: ~/.venvs/astro313/bin/python -I scripts/paper_review_figures.py
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, plt, BG2, FG, FGDIM, ACCENT, WARN, DANGER, OK, SUN,  # noqa: E402
                        OIII_C, HA_C, SII_C, NII_C)
from silk_geometry import ROOT, REGIONS, load_trace  # noqa: E402
import paper_review_lib as L  # noqa: E402

R = json.loads((ROOT / 'data/derived/paper_checks.json').read_text())
MC = np.load(ROOT / 'data/derived/paper_review_mc_samples.npz')
FIG = ROOT / 'figures'
T2, T4 = R['inputs']['table2'], R['inputs']['table4']
LB = R['regions_lb']
REGC = {'P1': '#e6d84f', 'P2': '#7ad97a', 'P3': '#e070e0', 'P4': '#ff9f43'}


def style_ax(ax):
    ax.grid(True, lw=0.4, alpha=0.5)
    for s in ax.spines.values():
        s.set_color('#2a323c')


# ---------------------------------------------------------------------------------------------
# A. Kinematic distance
# ---------------------------------------------------------------------------------------------
def fig_kinematics():
    fig, axs = new_sky_figure(16, 5.8, 'Silk Nebula · kinematic distance check',
                              'Wu+2026 Table 4 V_LSR vs Galactic rotation at l≈204°, b≈19.5°', ncols=3)
    ax = axs[0]
    l, b = LB['P1']
    d = np.logspace(np.log10(0.02), np.log10(5), 400)
    curves = [('flat_R0_8.15_T236_stdsun', 'flat 236 km/s', FGDIM, '-'),
              ('reid14_A5_linear', 'Reid+14 linear', ACCENT, '-'),
              ('kd_reid14_URC', 'kd Reid+14 URC', SUN, '--'),
              ('kd_reid19_URC', 'kd Reid+19 URC', OK, ':')]
    for key, lab, col, ls in curves:
        ax.plot(d, L.MODELS[key](l, b, d), color=col, ls=ls, lw=1.6, label=lab)
    v14 = L.MODELS['reid14_A5_linear'](l, b, d)
    ax.fill_between(d, v14 - 7, v14 + 7, color=ACCENT, alpha=0.10, lw=0, label='±7 km/s cloud motions')
    ax.axhspan(11.2, 15.5, color=OIII_C, alpha=0.25, lw=0)
    ax.axhspan(-1.0, 6.1, color=HA_C, alpha=0.15, lw=0)
    label(ax, 0.022, 13.4, '[O III] V_LSR 11–16', color=OIII_C, size=8.5, va='center')
    label(ax, 0.022, 2.6, 'Hα/[N II]/[S II] −1…+6', color=HA_C, size=8.5, va='center')
    for dd, txt, col, yy in ((0.16, '160 pc\n(paper)', DANGER, 28.5), (0.288, '288 pc\nPSR B0656', WARN, 22.5), (1.6, '1.6 kpc\n(kd kpc)', OK, 28.5)):
        ax.axvline(dd, color=col, lw=1.0, ls='--')
        label(ax, dd * 1.06, yy, txt, color=col, size=8, va='top')
    ax.set_xscale('log'); ax.set_xlim(0.02, 5); ax.set_ylim(-12, 30)
    ax.set_xlabel('distance (kpc)'); ax.set_ylabel('V_LSR (km/s)')
    ax.set_title('(a) rotation model toward P1', loc='left')
    ax.legend(loc='lower right', ncol=1, fontsize=7.5)
    style_ax(ax)

    ax = axs[1]
    bins = np.linspace(0, 6, 61)
    for sp, col in (('0', SUN), ('5', ACCENT), ('10', DANGER)):
        s = MC[f'kin_P1_OIII_pec{sp}']
        kin = R['kinematic_distance']['monte_carlo']['reid14_linear']['P1_OIII'][f'sig_pec_{sp}']
        ax.hist(s, bins=bins, density=True, histtype='step', lw=1.8, color=col,
                label=f'σ_pec={sp:>2} km/s: {kin["p16"]:.1f}–{kin["p84"]:.1f} kpc; P(<0.35)={kin["P_d_lt_0.35kpc"]:.2f}')
    ax.axvline(0.16, color=DANGER, ls='--', lw=1); ax.axvline(1.63, color=OK, ls='--', lw=1)
    label(ax, 0.2, 0.62, '160 pc', color=DANGER, size=8)
    label(ax, 1.68, 0.62, 'Table 4: "1.63"', color=OK, size=8)
    label(ax, 0.12, 0.93, '← d=0 pile-up\n   (V below rotation)', color=DANGER, size=7.5, va='top')
    ax.set_xlim(0, 6); ax.set_ylim(0, 1.1)
    ax.set_xlabel('kinematic distance (kpc)'); ax.set_ylabel('probability density')
    ax.set_title('(b) Monte Carlo, P1 [O III] 15.5±4.3 km/s', loc='left')
    ax.legend(loc='upper right', fontsize=7.2, title='16–84 % range; P(d<0.35 kpc)', title_fontsize=7.2)
    style_ax(ax)

    ax = axs[2]
    lines = [('OIII', '[O III]', OIII_C, 'o'), ('Ha', 'Hα', HA_C, 's'), ('NII', '[N II]', NII_C, 'D'), ('SII', '[S II]', SII_C, '^')]
    for i, reg in enumerate(('P1', 'P2', 'P3', 'P4')):
        for j, (ln, lab, col, mk) in enumerate(lines):
            if ln in T4[reg]:
                vh, sv, vl = T4[reg][ln][:3]
                ax.errorbar(i + (j - 1.5) * 0.12, vl, yerr=sv, fmt=mk, color=col, ms=6, capsize=2,
                            label=lab if i == 0 else None)
        w = R['sii_density']['wham_ss_dr1'][reg][0]
        ax.errorbar(i + 0.42, w['V_LSR_gauss'], yerr=max(w['V_LSR_gauss_err'], 1.0), fmt='*', color=OK, ms=10, capsize=2,
                    label='WHAM Hα (1° beam)' if i == 0 else None)
    ev = R['kinematic_distance']['expected_vlsr']['P1']['reid14_A5_linear']
    for key, txt, col in (('0.16kpc', 'rotation @160 pc', DANGER), ('0.3kpc', '@300 pc', WARN), ('1.6kpc', '@1.6 kpc', OK)):
        ax.axhline(ev[key], color=col, ls='--', lw=0.9)
        label(ax, -0.45, ev[key] + (0.5 if key != '0.16kpc' else -1.4), txt, color=col, size=7.5, ha='left')
    ax.set_xticks(range(4)); ax.set_xticklabels(['P1', 'P2', 'P3', 'P4 (bkg)'])
    ax.set_xlim(-0.5, 3.7); ax.set_ylim(-6, 32)
    ax.set_ylabel('V_LSR (km/s)')
    ax.set_title('(c) measured velocities by region', loc='left')
    ax.legend(loc='upper right', fontsize=7.5, ncol=2)
    style_ax(ax)
    return finish(fig, FIG / 'paper_review_kinematics.png')


# ---------------------------------------------------------------------------------------------
# B. [S II] density and emission-measure self-consistency
# ---------------------------------------------------------------------------------------------
def fig_density():
    fig, axs = new_sky_figure(13.5, 5.8, 'Silk Nebula · [S II] density vs emission measure',
                              'Table 2 ratios with standard atomic data; path length L = EM / n_e²', ncols=2)
    ax = axs[0]
    rr = np.linspace(0.46, 1.445, 400)
    ax.plot(rr, L.ne_sanders(rr), color=ACCENT, lw=1.8, label='Sanders+16 (T=10⁴ K)')
    rp = np.linspace(0.45, 1.43, 400)
    ax.plot(rp, L.ne_proxauf(rp), color=SUN, lw=1.4, ls='--', label='Proxauf+14 eq. 5 (paper)')
    ax.axhspan(64, 150, color=DANGER, alpha=0.18, lw=0)
    label(ax, 0.5, 98, 'paper: n_e 64–150', color=DANGER, size=8.5, va='center')
    ax.axvspan(1.253, 1.356, color=DANGER, alpha=0.10, lw=0)
    ax.axvline(1.45, color=FGDIM, ls=':', lw=1)
    label(ax, 1.435, 1.3, 'low-density\nlimit 1.45', color=FGDIM, size=7.5, ha='right', va='bottom')
    for reg in ('P1', 'P2', 'P3', 'P4'):
        dr = R['sii_density']['regions'][reg]
        Rt = dr['R_table']
        lo, md, hi = dr['R_MC_p16_50_84']
        ne = L.ne_sanders(Rt)
        ax.errorbar(Rt, ne, xerr=[[Rt - lo], [hi - Rt]], fmt='o', color=REGC[reg], ms=7, capsize=3, label=f'{reg}: R={Rt:.2f} → {ne:.0f} cm⁻³')
    ax.set_yscale('log'); ax.set_xlim(0.45, 1.5); ax.set_ylim(1, 2e4)
    ax.set_xlabel('[S II] 6716/6731'); ax.set_ylabel('n_e (cm⁻³)')
    ax.set_title('(a) ratio → density', loc='left')
    ax.legend(loc='upper right', fontsize=7.5)
    style_ax(ax)

    ax = axs[1]
    ne = np.logspace(-1, 3.3, 200)
    em_o3 = (L.em_from_oiii(10, 1e4), L.em_from_oiii(13, 1e4))
    em_o3_hot = (L.em_from_oiii(10, 3e4), L.em_from_oiii(13, 3e4))
    ems = [R['self_consistency']['EM_from_WHAM_Ha'][r]['EM_Ha_WHAM_T1e4'] for r in ('P1', 'P2', 'P3', 'P4')]
    ax.fill_between(ne, em_o3[0] / ne ** 2, em_o3[1] / ne ** 2, color=OIII_C, alpha=0.45, lw=0, label='[O III] 10–13 R, T=10⁴ K')
    ax.fill_between(ne, em_o3_hot[0] / ne ** 2, em_o3_hot[1] / ne ** 2, color=OIII_C, alpha=0.18, lw=0, label='[O III] 10–13 R, T=3×10⁴ K')
    ax.fill_between(ne, min(ems) / ne ** 2, max(ems) / ne ** 2, color=HA_C, alpha=0.30, lw=0, label='WHAM Hα 2.9–4.4 R (total)')
    tr = R['sizes']['trace']['per_distance']
    for key, txt, col in (('160pc_paper', '160 pc', DANGER), ('300pc_Monogem', '300 pc', WARN), ('1.6kpc_kd_units', '1.6 kpc', OK)):
        w = tr[key]['main_fwhm_pc']
        ax.axhline(w, color=col, ls='--', lw=1)
        em_c = L.em_from_oiii(11.5, 1e4)
        ne_x = np.sqrt(em_c / w)
        ax.plot(ne_x, w, 'o', color=col, ms=7, mec='k', zorder=5)
        label(ax, 5.5, w * (1.25 if key != '160pc_paper' else 0.55),
              f'strand FWHM 0.28° @ {txt} = {w:.2g} pc → n_e≈{ne_x:.1f}', color=col, size=7.5, va='center')
    ax.axvspan(64, 150, color=DANGER, alpha=0.18, lw=0)
    ax.axvspan(190, 440, color=SUN, alpha=0.12, lw=0)
    label(ax, 64, 3e-6, ' paper', color=DANGER, size=8)
    label(ax, 195, 3e-6, ' Table 2 R', color=SUN, size=8)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.1, 2000); ax.set_ylim(1e-6, 1e3)
    ax.set_xlabel('n_e (cm⁻³)'); ax.set_ylabel('required path length L (pc)')
    ax.set_title('(b) emission measure → path length', loc='left')
    ax.legend(loc='upper right', fontsize=7.5)
    style_ax(ax)
    return finish(fig, FIG / 'paper_review_density_em.png')


# ---------------------------------------------------------------------------------------------
# C. Line ratios vs WIM, and H-alpha budget
# ---------------------------------------------------------------------------------------------
def fig_ratios():
    fig, axs = new_sky_figure(13.5, 5.8, 'Silk Nebula · red-arm line ratios vs Galactic WIM',
                              'WHAM convention: [N II]6583/Hα, [S II]6716/Hα (energy units)', ncols=2)
    ax = axs[0]
    m6 = R['classification']['wham_literature']['Madsen2006_Table5_high_latitude_filaments_l217-236']
    first = True
    for k, v in m6.items():
        if k == 'columns':
            continue
        ax.plot(v[0], v[1], 'o', color=FGDIM, ms=6, mfc='none', label='Madsen+06 high-|b| filaments' if first else None)
        if k.startswith('Arc'):
            label(ax, v[0] - 0.01, v[1] - 0.035, 'Arc (I_Hα 1.3 R)', color=FGDIM, size=7, ha='right')
        first = False
    ax.plot(0.27, 0.11, '*', color=WARN, ms=13, label='Madsen+06 O-star H II mean')
    ax.axvspan(0.4, 0.6, ymin=0, ymax=1, color=FGDIM, alpha=0.08, lw=0)
    label(ax, 0.5, 0.06, 'WIM typical\n[N II]/Hα≈0.5', color=FGDIM, size=7.5, ha='center')
    for reg in ('P1', 'P2', 'P3', 'P4'):
        c = R['classification']['regions'][reg]
        ax.errorbar(c['NII6583_Ha'][0], c['SII6716_Ha'][0], xerr=c['NII6583_Ha'][1], yerr=c['SII6716_Ha'][1],
                    fmt='o', color=REGC[reg], ms=8, capsize=3, label=f'LAMOST {reg}' + (' (nebula-free)' if reg == 'P4' else ''))
    ax.axhline(0.4 / 1.704, color=DANGER, ls='--', lw=1)
    label(ax, 0.98, 0.4 / 1.704 + 0.01, '"shock" cut [S II]/Hα=0.4 (sum)', color=DANGER, size=7.5, ha='right', va='bottom')
    ax.set_xlim(0, 1.0); ax.set_ylim(0, 0.65)
    ax.set_xlabel('[N II] 6583 / Hα'); ax.set_ylabel('[S II] 6716 / Hα')
    ax.set_title('(a) nebula = background = WIM', loc='left')
    ax.legend(loc='upper left', fontsize=7.5)
    style_ax(ax)

    ax = axs[1]
    x = np.arange(4)
    regs = ('P1', 'P2', 'P3', 'P4')
    wham = [R['sii_density']['wham_ss_dr1'][r][0]['I_Ha_R'] for r in regs]
    cal = R['self_consistency']['Ha_excess_on_filament']['R_per_LAMOST_unit_mean']
    lam = [T2[r]['Ha'][0] * cal for r in regs]
    ax.bar(x - 0.27, wham, 0.25, color=HA_C, alpha=0.75, label='WHAM Hα total (R)')
    ax.bar(x, lam, 0.25, color=HA_C, alpha=0.35, label=f'LAMOST Hα × {cal:.1f} R/unit')
    exc = R['self_consistency']['Ha_excess_on_filament']['regions']
    for i, r in enumerate(regs[:3]):
        e = exc[r]
        ax.errorbar(i + 0.0, e['Ha_excess_R'], yerr=e['Ha_excess_err_R'], fmt='o', color=SUN, ms=6, capsize=3, mec='k',
                    label='Hα excess over P4 (±1σ)' if i == 0 else None)
    label(ax, 3.0, 4.0, 'no [O III]', color=FGDIM, size=8, ha='center')
    o3 = [R['surface_brightness']['regions'][r]['R_AB_60A'] for r in regs[:3]]
    ax.bar(x[:3] + 0.27, o3, 0.25, color=OIII_C, alpha=0.8, label='[O III] 5007 from Table 3 (R)')
    for i, r in enumerate(regs[:3]):
        lim = R['self_consistency']['oiii_ha_limits'][r]
        label(ax, i + 0.27, o3[i] + 0.4, f'>{lim["OIII_over_Ha_photon_2sigma_LL"]:.0f}×', color=OIII_C, size=8, ha='center')
    ax.set_xticks(x); ax.set_xticklabels(['P1', 'P2', 'P3', 'P4 (bkg)'])
    ax.set_ylabel('intensity (Rayleigh)'); ax.set_ylim(-0.5, 15.5)
    ax.set_title('(b) Hα budget vs [O III] (2σ [O III]/Hα limits)', loc='left')
    ax.legend(loc='upper right', fontsize=7.5)
    style_ax(ax)
    return finish(fig, FIG / 'paper_review_line_ratios.png')


# ---------------------------------------------------------------------------------------------
# D. Line-width and doublet-ratio simulations
# ---------------------------------------------------------------------------------------------
def fig_fits():
    fig, axs = new_sky_figure(13.5, 5.6, 'Silk Nebula · what S/N≈6 [O III] fits can measure',
                              'Gaussian+linear fits, R=7500 (FWHM_inst 40 km/s), paper\'s FWHM cut applied', ncols=2)
    ax = axs[0]
    bins = np.arange(0, 94, 4.0)
    ymax = 0
    for w, col in (('5.4', FGDIM), ('20', ACCENT), ('30', OK), ('45', SUN)):
        a = MC[f'fit_SN6.1_Wint{w}']
        pas = a[:, 2].astype(bool)
        f = R['line_widths']['fit_mc'][f'SN6.1_Wint{w}']
        h, _, _ = ax.hist(a[pas, 3], bins=bins, histtype='step', lw=1.8, color=col, density=True,
                          label=f'true {w:>4} km/s: pass {f["frac_pass_cut"]:.0%}, median {f["Wrec_pass_p16_50_84"][1]:.0f}')
        ymax = max(ymax, h.max())
    ax.axvspan(43.1, 48.9, color=OIII_C, alpha=0.25, lw=0)
    label(ax, 46, ymax * 0.62, 'paper\n43–49', color=OIII_C, size=8.5, ha='center', va='center')
    label(ax, 2, ymax * 1.05, 'true width 5–45 km/s →\nsame 36–47 km/s', color=FG, size=8, va='top')
    ax.set_xlim(0, 92); ax.set_ylim(0, ymax * 1.15)
    ax.set_xlabel('recovered intrinsic FWHM of fits passing the cut (km/s)'); ax.set_ylabel('density')
    ax.set_title('(a) [O III] width at S/N 6.1', loc='left')
    ax.legend(loc='upper right', fontsize=7.5)
    style_ax(ax)

    ax = axs[1]
    rat = MC['ratio_mc']; rat = rat[np.isfinite(rat)]
    ax.hist(np.clip(rat, 0, 10), bins=np.linspace(0, 10, 51), histtype='stepfilled', color=OIII_C, alpha=0.35, density=True,
            label='simulated, true 2.98')
    ax.axvline(2.98, color=SUN, lw=1.5, label='atomic value 2.98')
    for i, reg in enumerate(('P1', 'P2', 'P3')):
        o = R['oiii_doublet']['regions'][reg]
        ax.errorbar(o['ratio'], 0.52 - 0.05 * i, xerr=o['sigma'], fmt='o', color=REGC[reg], capsize=3, label=f'{reg}: {o["ratio"]:.2f}±{o["sigma"]:.2f}')
    m = R['oiii_doublet']['mc_independent_fits']
    label(ax, 9.8, 0.30, f'sim. 16–84%: {m["ratio_p16_50_84"][0]:.1f}–{m["ratio_p16_50_84"][2]:.1f}\n'
          f'S/N(4959)<3 in {m["frac_SN4959_below_3"]:.0%} of trials', color=FG, size=8, ha='right')
    label(ax, 9.75, 0.16, '≥10', color=FGDIM, size=7.5, ha='center')
    ax.set_xlim(0, 10); ax.set_ylim(0, 0.6)
    ax.set_xlabel('[O III] 5007/4959 (independent fits)'); ax.set_ylabel('density')
    ax.set_title('(b) doublet ratio at S/N 6.1 / ~2', loc='left')
    ax.legend(loc='upper right', fontsize=7.5)
    style_ax(ax)
    return finish(fig, FIG / 'paper_review_fit_sims.png')


# ---------------------------------------------------------------------------------------------
# E. Energetics and physical scales vs distance
# ---------------------------------------------------------------------------------------------
def fig_energetics():
    fig, axs = new_sky_figure(13.5, 5.6, 'Silk Nebula · size and SNR energetics vs distance',
                              'Cioffi+88 PDS energy (Boumis+22 form) with R_s scaled from 30 pc at 160 pc', ncols=2)
    ax = axs[0]
    d = np.logspace(np.log10(0.08), np.log10(2.5), 200)
    rs = 30.0 * d / 0.16
    for n0, col in ((1.0, SUN), (0.3, ACCENT), (0.1, OK), (0.01, FGDIM)):
        vs = L.vs_dopita(64.0, n0)
        ax.plot(d, L.e_boumis(n0, vs, rs), color=col, lw=1.8, label=f'n0={n0:g}, V_s(n_e=64)={vs:.0f} km/s')
    ax.plot(d, L.e_boumis(1.0, 119.0, rs), color=DANGER, lw=1.2, ls='--', label='paper inputs (n0=1, 119 km/s)')
    ax.plot(0.16, 2.0e51, 'o', color=DANGER, ms=8, mec='k', zorder=5)
    label(ax, 0.17, 3.2e51, 'paper: 2×10⁵¹ erg', color=DANGER, size=8, va='bottom')
    ax.plot(1.6, L.e_boumis(1.0, 119.0, 300.0), 'o', color=OK, ms=8, mec='k', zorder=5)
    label(ax, 1.45, 5e54, '2.9×10⁵⁴ erg', color=OK, size=8, ha='right')
    ax.axhspan(5e50, 3e51, color=OK, alpha=0.12, lw=0)
    label(ax, 0.085, 1.3e51, 'single SN', color=OK, size=8.5, va='center')
    for dd, txt, col in ((0.16, '160 pc', DANGER), (0.288, '288 pc', WARN), (1.6, '1.6 kpc', OK)):
        ax.axvline(dd, color=col, ls=':', lw=1)
        label(ax, dd * 1.05, 3e55, txt, color=col, size=8, va='top')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.08, 2.5); ax.set_ylim(1e49, 1e56)
    ax.set_xlabel('distance (kpc)'); ax.set_ylabel('E (erg)')
    ax.set_title('(a) explosion energy', loc='left')
    ax.legend(loc='lower right', fontsize=7.3)
    style_ax(ax)

    ax = axs[1]
    tr = R['sizes']['trace']
    arc, fw = tr['spine_arc_length_deg'], tr['main_strand_fwhm_deg']
    bc = 19.44
    ax.plot(d, np.deg2rad(arc) * d * 1e3, color=OIII_C, lw=1.8, label=f'spine length ({arc:.1f}°)')
    ax.plot(d, 2 * d * 1e3 * np.tan(np.deg2rad(5)), color=OIII_C, lw=1.2, ls='--', label='paper extent (10°)')
    ax.plot(d, np.deg2rad(fw) * d * 1e3, color=ACCENT, lw=1.4, ls=':', label=f'main-strand FWHM ({fw:.2f}°)')
    ax.plot(d, d * 1e3 * np.sin(np.deg2rad(bc)) + 20.8, color=WARN, lw=1.8, label='height above midplane')
    for dd, txt, col in ((0.16, '160 pc', DANGER), (0.288, '288 pc', WARN), (1.6, '1.6 kpc', OK)):
        ax.axvline(dd, color=col, ls=':', lw=1)
        label(ax, dd * 1.05, 1.5e3, txt, color=col, size=8, va='top')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.08, 2.5); ax.set_ylim(0.2, 2e3)
    ax.set_xlabel('distance (kpc)'); ax.set_ylabel('pc')
    ax.set_title('(b) physical scales', loc='left')
    ax.legend(loc='lower right', fontsize=7.5)
    style_ax(ax)
    return finish(fig, FIG / 'paper_review_energetics.png')


# ---------------------------------------------------------------------------------------------
# F. Sky context: H-alpha (Finkbeiner 2003) with trace, regions, LAMOST field, WHAM beams
# ---------------------------------------------------------------------------------------------
def fig_sky():
    import healpy as hp
    from astropy.io import fits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    fk = fits.open(ROOT / 'data/raw/lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits')[1].data.field(0).ravel().astype(float)
    ra = np.linspace(110.5, 122.5, 361); dec = np.linspace(8, 30, 661)
    RA, DE = np.meshgrid(ra, dec)
    g = SkyCoord(RA.ravel() * u.deg, DE.ravel() * u.deg).galactic
    img = hp.get_interp_val(fk, g.l.deg, g.b.deg, nest=True, lonlat=True).reshape(RA.shape)
    fig, ax = new_sky_figure(7.6, 11.0, 'Silk Nebula · Hα context of the LAMOST regions',
                             'Finkbeiner (2003) Hα (WHAM + VTSS patchwork); [O III] trace v1.0')
    im = ax.imshow(img, origin='lower', extent=(ra[0], ra[-1], dec[0], dec[-1]), cmap='magma', vmin=1.5, vmax=6.5,
                   aspect=1 / np.cos(np.deg2rad(19.0)))
    ax.plot([], [], ':', color=FG, lw=0.8, label='WHAM 1° beams used (I_Hα)')
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label('Hα (R)')
    tr = load_trace()
    if tr:
        ax.plot([p['ra'] for p in tr['spine']], [p['dec'] for p in tr['spine']], color=OIII_C, lw=1.6, label='[O III] spine')
        env = np.array(tr['envelope_polygon'])
        if env.ndim == 2 and env.shape[1] >= 2:
            ax.plot(np.r_[env[:, 0], env[0, 0]], np.r_[env[:, 1], env[0, 1]], color=OIII_C, lw=0.7, alpha=0.6, label='[O III] envelope')
        hpch = np.array(tr.get('halpha_patch', []))
        if hpch.size:
            ax.plot(np.r_[hpch[:, 0], hpch[0, 0]], np.r_[hpch[:, 1], hpch[0, 1]], color=HA_C, lw=1.0, label='Hα patch (N end)')
        lf = tr.get('lamost_field', {})
        if lf:
            t = np.linspace(0, 2 * np.pi, 200)
            cen = SkyCoord(lf['ra'] * u.deg, lf['dec'] * u.deg)
            circ = cen.directional_offset_by(t * u.rad, lf['radius_deg'] * u.deg)
            ax.plot(circ.ra.deg, circ.dec.deg, color=FG, lw=1, ls='--', label='LAMOST plate (r=2.5°)')
    for reg, (r0, d0) in REGIONS.items():
        ax.plot(r0, d0, 'o', color=REGC[reg], ms=7, mec='k')
        label(ax, r0 - 0.25, d0, reg, color=REGC[reg], size=10, va='center', fontweight='bold')
        w = R['sii_density']['wham_ss_dr1'][reg][0]
        c = SkyCoord(l=w['l'] * u.deg, b=w['b'] * u.deg, frame='galactic').icrs
        t = np.linspace(0, 2 * np.pi, 100)
        circ = c.directional_offset_by(t * u.rad, 0.5 * u.deg)
        ax.plot(circ.ra.deg, circ.dec.deg, color=REGC[reg], lw=0.8, ls=':')
        label(ax, c.ra.deg + 0.55, c.dec.deg - 0.55, f'{w["I_Ha_R"]:.1f} R', color=REGC[reg], size=7.5)
    ax.set_xlim(122.5, 110.5); ax.set_ylim(8, 30)
    ax.set_xlabel('RA (deg, J2000)'); ax.set_ylabel('Dec (deg, J2000)')
    ax.legend(loc='lower left', fontsize=7.5)
    return finish(fig, FIG / 'paper_review_sky_context.png')


if __name__ == '__main__':
    for f in (fig_kinematics, fig_density, fig_ratios, fig_fits, fig_energetics, fig_sky):
        print(f())
