"""dispute_other_fig6.py -- Comment section 6 ("further problems"): numbers for Table 7 and Figure 6.

Recomputes, from the values printed in Wu et al. 2026 (arXiv:2610.13329v1 = AJ 172, 298; Table 2 as transcribed in
data/derived/paper_checks.json -> inputs.table2), the red-arm line ratios, the [S II] densities, the Eq. 2/3
temperature and turbulence, the N2 / N2S2Ha abundances, the energetics arithmetic, the [O III] doublet S/N and the
emission-measure path lengths.  Reuses (read-only imports / JSON reads):
  scripts/verify_paper_sii.py         5-level S+ atom (TZ10 / RBS96 Omega; TZ10, MZ82, RGJ19 A) and Eq. 5 (Proxauf+14)
  scripts/verify_paper_energetics.py  Dopita 1979 / Cioffi+88 / Boumis+22 formulae
  data/derived/paper_checks.json      Table 2/4 transcription, WHAM-SS DR1 beams, [O III] width Monte Carlo
  data/derived/verify_paper_energetics.json   distance grid, Balmer-bounded n0 and E
  data/derived/w3_oii_summary.json    WIM [S II] doublet on this sightline (SDSS, DESI regressions)
Literature reference values (Madsen, Reynolds & Haffner 2006 Table 5 and Sect. 4-5 text; Pettini & Pagel 2004;
Marino et al. 2013; Haffner et al. 2003 EM conversion) are typed in below with their source.
Outputs: data/derived/dispute_other_values.json, figures/dispute_other_fig6.png
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_other_fig6.py
"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import verify_paper_sii as S                      # noqa: E402  (reads PyNeb text tables only)
import verify_paper_energetics as EN              # noqa: E402  (functions only)

DER = ROOT / 'data/derived'
PC_AU = 206264.806
REG = ['P1', 'P2', 'P3', 'P4']
rng = np.random.default_rng(20261117)

pc = json.load(open(DER / 'paper_checks.json'))
ven = json.load(open(DER / 'verify_paper_energetics.json'))
w3o = json.load(open(DER / 'w3_oii_summary.json'))
T2 = pc['inputs']['table2']
T4_VT = pc['inputs']['table4_Vt']
out = {'source': 'scripts/dispute_other_fig6.py', 'paper': 'Wu et al. 2026, AJ 172, 298 (arXiv:2610.13329v1)'}


def ratio_err(a, b):
    """a, b = (value, err); independent errors."""
    r = a[0] / b[0]
    return r, r * np.hypot(a[1] / a[0], b[1] / b[0])


def sum_err(*xs):
    return sum(x[0] for x in xs), float(np.sqrt(sum(x[1] ** 2 for x in xs)))


# ---------------------------------------------------------------- (1) red-arm ratios
rat = {}
for k in REG:
    t = T2[k]
    rat[k] = {
        'SII_sum_Ha': ratio_err(sum_err(t['S6716'], t['S6731']), t['Ha']),
        'SII6716_Ha': ratio_err(t['S6716'], t['Ha']),
        'NII6583_Ha': ratio_err(t['N6583'], t['Ha']),
        'NII_sum_Ha': ratio_err(sum_err(t['N6548'], t['N6583']), t['Ha']),
        'table2_SII_Ha': t['ratio_SII_Ha'], 'table2_NII_Ha': t['ratio_NII_Ha'],
        'table2_SII_doublet': t['ratio_SII'],
        'SII_doublet_from_intensities': t['S6716'][0] / t['S6731'][0]}
for k in REG[:3]:
    rat[k]['delta_vs_P4_sigma'] = {q: (rat[k][q][0] - rat['P4'][q][0]) / np.hypot(rat[k][q][1], rat['P4'][q][1])
                                   for q in ('SII_sum_Ha', 'SII6716_Ha', 'NII6583_Ha', 'NII_sum_Ha')}
out['red_ratios'] = rat
# literature reference values (energy ratios, WHAM; Madsen, Reynolds & Haffner 2006, ApJ 652, 401)
M06_T5 = {  # Table 5, high-latitude filaments: ([N II]6583/Ha, [S II]6716/Ha)
    'A1': (0.49, 0.32), 'A2': (0.44, 0.29), 'B1': (0.63, 0.33), 'B2': (0.58, 0.26), 'C': (0.27, 0.26),
    'D': (0.46, 0.29), 'E1': (0.50, 0.40), 'E2': (0.40, 0.36), 'F': (0.46, 0.23), 'Arc': (0.72, 0.46)}
n2 = [v[0] for v in M06_T5.values()]; s2 = [v[1] for v in M06_T5.values()]
R_LDL = 1.454   # TZ10 5-level, 1e4 K (verify_paper_sii.json limits)
out['literature_refs'] = {
    'M06_Table5_NII6583_Ha_range': [min(n2), max(n2)], 'M06_Table5_SII6716_Ha_range': [min(s2), max(s2)],
    'M06_Table5_SII_sum_Ha_range_derived_LDL': [min(s2) * (1 + 1 / R_LDL), max(s2) * (1 + 1 / R_LDL)],
    'M06_faint_WIM_SII_Ha': '>~0.4 in regions of faint Ha (Sect. 5)',
    'M06_faintest_WIM_NII_Ha': '>~0.7 for the WIM at l > 150 deg (Sect. 5)',
    'M06_typical_WIM_NII_Ha': '~0.5, in some cases > 1.0 (Sect. 4)',
    'M06_Ostar_HII_mean': {'NII_Ha': 0.27, 'SII_Ha': 0.11},
    'paper_SNR_threshold_SII_Ha': {'abstract_sect5': 0.4, 'sect1': 0.45, 'sect4.1': [0.45, 0.5]}}

# ---------------------------------------------------------------- (2) [S II] density
combos = [('TZ10', 'TZ10'), ('TZ10', 'MZ82'), ('TZ10', 'RGJ19'), ('RBS96', 'MZ82')]
Rtab = {k: T2[k]['ratio_SII'] for k in REG}
dens = {'five_level': {}, 'eq5_1e4K': {}, 'eq5_8000K_sqrtT': {}}
for k, R in Rtab.items():
    dens['five_level'][k] = {f'{c}+{a}@{int(T)}': float(S.ne_from_R(R, T, coll=c, aval=a))
                             for c, a in combos for T in (8000.0, 1e4)}
    dens['eq5_1e4K'][k] = float(S.proxauf(R))
    dens['eq5_8000K_sqrtT'][k] = float(S.proxauf(R) * np.sqrt(0.8))
allv = [v for k in REG for v in dens['five_level'][k].values()]
tzv = [v for k in REG for key, v in dens['five_level'][k].items() if key.startswith('TZ10+')]
dens['five_level_range_all'] = [min(allv), max(allv)]
dens['five_level_range_TZ10_Omega'] = [min(tzv), max(tzv)]
dens['eq5_range'] = [min(list(dens['eq5_8000K_sqrtT'].values())), max(list(dens['eq5_1e4K'].values()))]
dens['R_LDL_5level'] = {f'{c}+{a}@{int(T)}': float(S.ratio(1e-3, T, coll=c, aval=a)) for c, a in combos for T in (8000.0, 1e4)}
rr = np.linspace(1.0, 1.44, 44001); pp = S.proxauf(rr)
dens['R_for_ne'] = {str(n): {'eq5': float(rr[np.argmin(np.abs(pp - n))]),
                             '5level_TZ10_1e4': float(S.ratio(n, 1e4)), '5level_TZ10_8000': float(S.ratio(n, 8000.0))}
                    for n in (64, 150)}
# uncertainty of the Table 2 doublet ratio: Table 2 gives no error on R.  Two readings:
#  (i) per-line S/N from the Table 2 S/N column (44-65 for [S II]): sigma_R = R*sqrt(2)/SN
#  (ii) the quoted +-0.02 intensity errors taken as independent per line (common-mode sky-line scaling would cancel)
sig = {}
for k in REG:
    t = T2[k]; R = Rtab[k]; sn = t['SN_red']['S']
    s_sn = R * np.sqrt(2) / sn
    s_int = (t['S6716'][0] / t['S6731'][0]) * np.hypot(t['S6716'][1] / t['S6716'][0], t['S6731'][1] / t['S6731'][0])
    sig[k] = {'SN_SII': sn, 'sigma_R_from_SN': s_sn, 'sigma_R_from_intensity_errors': s_int,
              'nsigma_to_R_ne64_eq5_SN': (dens['R_for_ne']['64']['eq5'] - R) / s_sn,
              'nsigma_to_R_ne64_eq5_int': (dens['R_for_ne']['64']['eq5'] - R) / s_int}
for k in REG[:3]:
    sig[k]['P4_minus_region_nsigma_SN'] = (Rtab[k] - Rtab['P4']) / np.hypot(sig[k]['sigma_R_from_SN'], sig['P4']['sigma_R_from_SN'])
dens['R_uncertainty'] = sig
dens['wim_sightline_SII_doublet'] = {'SDSS_WIM_between': [w3o['sdss_wim_control']['between'][x] for x in ('sii_fixed', 'sii_p16', 'sii_p84')],
                                     'SDSS_WIM_within': [w3o['sdss_wim_control']['within'][x] for x in ('sii_fixed', 'sii_p16', 'sii_p84')],
                                     'DESI_Ha_component_all': [w3o['atomic']['sii']['DESI_Ha_all'][x] for x in ('R', 'p16', 'p84')]}
# emission-measure argument: WHAM-SS DR1 nearest beam, EM = 2.25 pc cm^-6 per R at 8000 K (Haffner et al. 2003)
em = {}
for k in REG:
    I = pc['sii_density']['wham_ss_dr1'][k][0]['I_Ha_R']
    EM = 2.25 * I
    ne5 = dens['five_level'][k]['TZ10+TZ10@8000']
    em[k] = {'I_Ha_WHAM_R': I, 'EM_8000K': EM,
             'L_AU_at_ne_5level_TZ10_8000K': EM / ne5 ** 2 * PC_AU,
             'L_AU_at_ne_eq5_8000K': EM / dens['eq5_8000K_sqrtT'][k] ** 2 * PC_AU,
             'L_pc_at_ne64': EM / 64 ** 2, 'L_pc_at_ne150': EM / 150 ** 2,
             'ne_rms_over_1kpc': np.sqrt(EM / 1000.0), 'ne_rms_over_3kpc': np.sqrt(EM / 3000.0)}
dens['emission_measure'] = em
# a 1.37-pc filament path (0.28 deg at 280 pc) relative to the EM-implied paths at the Table 2 densities
_paths = [em[k][q] for k in REG for q in ('L_AU_at_ne_5level_TZ10_8000K', 'L_AU_at_ne_eq5_8000K')]
dens['filament_path_over_EM_path_range'] = [1.37 * PC_AU / max(_paths), 1.37 * PC_AU / min(_paths)]
# [O III] gas: brightness that n_e = 64 / 150 would give if it filled the 0.28-deg strand at 280 pc (L = 1.37 pc), x(O++) = 1,
# solar O/H; I(5007) per EM(n_e n_H) from verify_paper_energetics.json (own 5-level O++ atom), n_e/n_H = 1.15
osc = ven['oiii_selfconsistency']
L_str = osc['geometry']['paths_pc']['filament_L=w']
dens['oiii_gas_if_dense'] = {f'ne{n}_T{T}': n ** 2 / osc['ne_over_nH'] * L_str * osc['emissivity'][T]['I5007_per_EM_x1']
                             for n in (64, 150) for T in ('10000', '20000')}
dens['oiii_gas_if_dense']['L_pc'] = L_str
dens['oiii_gas_ne_from_EM_280pc_canonical'] = 'n_e 0.1-0.7 (edge-on sheet/shell), 0.5-3 (filament) for 2-4 R; verify_paper d1 / w3_canon Sect. 7'
out['sii_density'] = dens

# ---------------------------------------------------------------- (3) T_e (Eq. 2) and V_t (Eq. 3), paper coefficients
tev = {}
for k in REG:
    WH, WN = T2[k]['W_H'], T2[k]['W_N']
    T = 23.5 * (WH[0] ** 2 - WN[0] ** 2)
    Vt = 1.04 * WN[0] * np.sqrt(1 - 0.071 * WH[0] ** 2 / WN[0] ** 2)
    h = rng.normal(WH[0], WH[1], 400000); n = rng.normal(WN[0], WN[1], 400000)
    Ts = 23.5 * (h ** 2 - n ** 2)
    tev[k] = {'T_eq2_K': T, 'T_MC_p16_50_84': np.percentile(Ts, [16, 50, 84]).tolist(), 'T_sigma_K': float(np.std(Ts)),
              'P_T_negative': float(np.mean(Ts < 0)), 'Vt_eq3': Vt, 'Vt_table4': T4_VT[k]}
out['temperature_turbulence'] = tev

# ---------------------------------------------------------------- (4) abundances
ab = {}
for k in REG:
    t = T2[k]
    N2 = np.log10(t['N6583'][0] / t['Ha'][0])
    y = np.log10(t['N6583'][0] / (t['S6716'][0] + t['S6731'][0])) + 0.264 * N2
    ab[k] = {'N2': N2, 'PP04_linear_eq6': 8.90 + 0.57 * N2, 'PP04_cubic': 9.37 + 2.03 * N2 + 1.26 * N2 ** 2 + 0.32 * N2 ** 3,
             'Marino13_N2': 8.743 + 0.462 * N2, 'D16_y': y, 'D16_OH': 8.77 + y + 0.45 * (y + 0.3) ** 5,
             'N2_inside_PP04_range(-2.5,-0.3)': bool(-2.5 < N2 < -0.3), 'N2_inside_M13_range(-1.6,-0.2)': bool(-1.6 < N2 < -0.2)}
out['abundances'] = ab

# ---------------------------------------------------------------- (5) energetics (paper inputs)
Vs = EN.vs_dopita(64.0, 1.0); E = EN.E_cioffi(1.0, 119.0, 30.0)
eg = {'Vs_ne64_n01': float(Vs), 'E_n01_Vs119_R30': float(E), 'R_PDS': float(EN.R_pds(E, 1.0)),
      'E_for_RPDS_18': float(1e51 * (18 / 14.0) ** 3.5),
      'age_plus_third': float(EN.t_age(E, 1.0, 30.0, 1 / 3)), 'age_plus_one_printed': float(EN.t_age(E, 1.0, 30.0, 1.0)),
      'age_plus_third_RPDS18': float(9.98e3 * (E / 1e51) ** (3 / 14) * ((30 / 18) ** (10 / 3) + 1 / 3)),
      'age_plus_one_RPDS18': float(9.98e3 * (E / 1e51) ** (3 / 14) * ((30 / 18) ** (10 / 3) + 1))}
g = ven['energetics']['grid']
eg['E_n01_280pc'] = [g['d280_paper_10.74deg_ne64_n01.0']['E'], g['d280_Monogem_12.1deg_ne64_n01.0']['E']]
eg['E_n01_1600pc'] = [g['d1600_paper_10.74deg_ne64_n01.0']['E'], g['d1600_Monogem_12.1deg_ne64_n01.0']['E']]
hb = {k: v for k, v in ven['energetics']['brightness_tied'].items() if k.startswith('Hb')}
eg['Balmer_bound_n0_max'] = [min(v['n0_max'] for v in hb.values()), max(v['n0_max'] for v in hb.values())]
eg['Balmer_bound_E_max'] = [min(v['E_max'] for v in hb.values()), max(v['E_max'] for v in hb.values())]
eg['paper_like_error_MC'] = ven['energetics']['paper_error_mc']
eg['selfconsistent_Vs_ne64_150_p16_50_84'] = pc['energetics']['propagation']['selfconsistent_Vs_p16_50_84']
out['energetics'] = eg

# ---------------------------------------------------------------- (6) [O III] widths and doublet
od = {}
for k in REG[:3]:
    t = T2[k]
    r = t['O5007'] / t['O4959']
    od[k] = {'ratio_5007_4959': r, 'sigma': r * np.hypot(1 / t['SN_O5007'], 1 / t['SN_O4959']),
             'SN4959_paper': t['SN_O4959'], 'SN4959_expected_from_5007': t['SN_O5007'] / 2.98, 'W_O_table2': t['W_O']}
out['oiii_doublet'] = od
fm = pc['line_widths']['fit_mc']
out['oiii_width_mc'] = {k: {'pass_frac': v['frac_pass_cut'], 'pass_p16_50_84': v['Wrec_pass_p16_50_84']} for k, v in fm.items()}

# ---------------------------------------------------------------- (7) radio angular scales
out['radio'] = {'VLA_D_Lband_LAS_full_synthesis_arcsec': 970, 'snapshot_factor': 0.5,
                'NVSS_snapshot_LAS_arcmin': 970 / 2 / 60, 'TGSS_ADR1_extended_limit_arcmin': [10, 20],
                'strand_FWHM_arcmin': pc['sizes']['trace']['main_strand_fwhm_deg'] * 60,
                'spine_length_deg': pc['sizes']['trace']['spine_arc_length_deg']}


def jdef(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, tuple):
        return list(o)
    raise TypeError(type(o))


(DER / 'dispute_other_values.json').write_text(json.dumps(out, indent=1, default=jdef))


# ================================================================ Figure 6
def make_fig():
    from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, SII_C, NII_C, FGDIM, FG, SUN, WARN, OK, ACCENT, DANGER
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5, 'axes.labelsize': 10})
    fig, axs = new_sky_figure(13.0, 10.4, 'Wu+26 red-arm diagnostics, [S II] density, T$_e$ and [O III] width',
                              'LAMOST values from Wu et al. 2026 Table 2 (P4 = nebula-free reference); '
                              'models and literature as labelled', nrows=2, ncols=2)
    cols = {'P1': ACCENT, 'P2': OK, 'P3': NII_C, 'P4': '#ffffff'}
    mk = {'P1': 'o', 'P2': 'D', 'P3': '^', 'P4': 's'}
    LS = 8.3   # annotation size

    def pmark(k, **kw):
        return dict(marker=mk[k], color=cols[k], mfc=cols[k] if k != 'P4' else 'none', mew=1.5, **kw)

    # ---------------- (a) red-line ratios
    ax = axs[0, 0]
    cats = [('SII_sum_Ha', '[S II]/Hα\n(6716+6731)'), ('SII6716_Ha', '[S II]6716/Hα'), ('NII6583_Ha', '[N II]6583/Hα')]
    lit = out['literature_refs']
    bands = {0: lit['M06_Table5_SII_sum_Ha_range_derived_LDL'], 1: lit['M06_Table5_SII6716_Ha_range'],
             2: lit['M06_Table5_NII6583_Ha_range']}
    for i, (q, _) in enumerate(cats):
        lo, hi = bands[i]
        ax.fill_between([i - 0.38, i + 0.38], lo, hi, color=SII_C, alpha=0.16, lw=0, zorder=1)
        for j, k in enumerate(REG):
            v, e = rat[k][q]
            ax.errorbar(i - 0.21 + 0.14 * j, v, yerr=e, ls='', ms=7, capsize=2.5, lw=1.2, zorder=4, **pmark(k))
    ax.plot([0.62, 1.38], [0.11, 0.11], color=FGDIM, ls=':', lw=1.4)
    ax.plot([1.62, 2.38], [0.27, 0.27], color=FGDIM, ls=':', lw=1.4)
    ax.plot([1.62, 2.38], [0.70, 0.70], color=WARN, ls='--', lw=1.2)
    ax.plot([0.62, 1.38], [0.40, 0.40], color=WARN, ls='--', lw=1.2)
    ax.plot([-0.38, 0.38], [0.40, 0.40], color=DANGER, lw=1.4)
    ax.fill_between([-0.38, 0.38], 0.45, 0.50, color=DANGER, alpha=0.32, lw=0)
    label(ax, -0.37, 0.385, 'paper "SNR" cut\n0.4 / 0.45–0.5', color=DANGER, size=LS, va='top')
    ax.set_xticks(range(3)); ax.set_xticklabels([c[1] for c in cats], fontsize=9.3)
    ax.set_xlim(-0.5, 2.5); ax.set_ylim(0, 1.36)
    ax.set_ylabel('line ratio (energy)')
    ax.grid(axis='y', alpha=0.3)
    h = [Line2D([], [], ls='', ms=7, label=k + (' no nebula' if k == 'P4' else ''), **pmark(k)) for k in REG]
    h += [plt.Rectangle((0, 0), 1, 1, color=SII_C, alpha=0.35, label='high-b WIM filaments'),
          Line2D([], [], color=WARN, ls='--', lw=1.2, label='faint WIM lower bound'),
          Line2D([], [], color=FGDIM, ls=':', lw=1.4, label='O-star H II mean')]
    ax.legend(handles=h, loc='upper left', fontsize=8.0, ncol=4, framealpha=0.9, handletextpad=0.35,
              columnspacing=0.9, title='LAMOST (Wu+26)          Madsen+06', title_fontsize=7.6, alignment='left')
    ax.set_title('(a) Red-arm ratios: filament = background = WIM', loc='left', fontsize=10, color=FG)

    # ---------------- (b) [S II] density
    ax = axs[0, 1]
    ne = np.logspace(0, 4.5, 260)
    curves = [('TZ10', 'TZ10', 1e4, ACCENT, '-', '5-level Ω TZ10, A TZ10, 10⁴ K'),
              ('TZ10', 'TZ10', 8000.0, ACCENT, '--', '5-level Ω TZ10, A TZ10, 8000 K'),
              ('TZ10', 'RGJ19', 1e4, OK, ':', '5-level Ω TZ10, A RGJ19, 10⁴ K'),
              ('RBS96', 'MZ82', 8000.0, WARN, '-.', '5-level Ω RBS96, A MZ82, 8000 K')]
    hb = []
    for c, a, T, col, ls, lab in curves:
        hb += ax.plot(ne, [S.ratio(n, T, coll=c, aval=a) for n in ne], color=col, ls=ls, lw=1.5, label=lab)
    rr2 = np.linspace(1.0, 1.435, 400)
    hb += ax.plot(S.proxauf(rr2), rr2, color=SUN, lw=2.0, label='paper Eq. 5 (Proxauf+14), 10⁴ K')
    ax.axvspan(64, 150, color=HA_C, alpha=0.22, lw=0)
    wl = dens['wim_sightline_SII_doublet']['SDSS_WIM_between']
    ax.axhspan(wl[1], wl[2], color=FGDIM, alpha=0.14, lw=0)
    label(ax, 1.12, wl[2] + 0.012, 'WIM on this sightline (SDSS sky fibres, 16–84 %)', color=FGDIM, size=LS, va='bottom')
    for k in REG:
        R = Rtab[k]
        n5 = dens['five_level'][k]['TZ10+TZ10@10000']; nq = dens['eq5_1e4K'][k]
        ax.plot([1, n5], [R, R], color=SII_C, lw=0.6, alpha=0.55)
        ax.plot(n5, R, ls='', ms=7, zorder=5, **pmark(k))
        ax.plot(nq, R, mk[k], color=SUN, mfc='none', mew=1.3, ms=7, zorder=5)
        label(ax, 1.6, R + 0.006, f'{k} R = {R:.2f}', color=cols[k] if k != 'P4' else FG, size=LS, va='bottom')
    sr = dens['R_uncertainty']
    s_sn = [sr[k]['sigma_R_from_SN'] for k in REG]; s_in = [sr[k]['sigma_R_from_intensity_errors'] for k in REG]
    label(ax, 1.6, 0.985, f'σ$_R$ = {min(s_sn):.3f}–{max(s_sn):.3f} (Table 2 S/N 44–65)\n'
          f'or {min(s_in):.2f}–{max(s_in):.2f} (±0.02 intensities)', color=FGDIM, size=LS, va='top')
    hb += [plt.Rectangle((0, 0), 1, 1, color=HA_C, alpha=0.35, label='paper n$_e$ 64–150 cm$^{-3}$'),
           Line2D([], [], marker='o', ls='', color=FG, ms=6.5, label='Table 2 R on 5-level 10⁴ K curve'),
           Line2D([], [], marker='o', ls='', color=SUN, mfc='none', mew=1.3, ms=6.5, label='Table 2 R on Eq. 5')]
    ax.set_xscale('log'); ax.set_xlim(1, 3e4); ax.set_ylim(0.38, 1.58)
    ax.set_xlabel('n$_e$ (cm$^{-3}$)'); ax.set_ylabel('I(6716)/I(6731)')
    ax.legend(handles=hb, loc='lower left', fontsize=7.6, framealpha=0.9)
    ax.grid(alpha=0.3)
    ax.set_title('(b) [S II] doublet → n$_e$: Table 2 gives 170–630, not 64–150', loc='left', fontsize=10, color=FG)

    # ---------------- (c) T_e from Eq. 2
    ax = axs[1, 0]
    ax.axhspan(1e4, 2e4, color=HA_C, alpha=0.18, lw=0)
    ax.axhspan(6000, 9000, color=SII_C, alpha=0.15, lw=0)
    ax.plot([2.72, 3.28], [7500, 7500], color=DANGER, lw=1.8)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=HA_C, alpha=0.4, label='paper: [O III] structure 1–2×10⁴ K'),
                       Line2D([], [], color=DANGER, lw=1.8, label='paper: P4 7.5×10³ K'),
                       plt.Rectangle((0, 0), 1, 1, color=SII_C, alpha=0.35, label='WIM 6000–9000 K (Haffner+99)')],
              loc='upper left', ncol=2, fontsize=8.0, framealpha=0.9, columnspacing=1.0)
    ax.axhline(0, color=FGDIM, lw=0.8)
    for j, k in enumerate(REG):
        p16, p50, p84 = tev[k]['T_MC_p16_50_84']
        ax.errorbar(j, tev[k]['T_eq2_K'], yerr=[[tev[k]['T_eq2_K'] - p16], [p84 - tev[k]['T_eq2_K']]], ls='', ms=7.5,
                    capsize=3, lw=1.3, zorder=4, **pmark(k))
        label(ax, j, p16 - 450, f"P(T<0) = {tev[k]['P_T_negative'] * 100:.0f}%", color=FGDIM, size=LS, ha='center', va='top')
    ax.set_xticks(range(4)); ax.set_xticklabels(['P1', 'P2', 'P3', 'P4\n(no nebula)'])
    ax.set_xlim(-0.6, 3.6); ax.set_ylim(-10200, 2.75e4)
    ax.set_ylabel('T$_e$ from Eq. 2 (K)')
    ax.grid(axis='y', alpha=0.3)
    ax.set_title('(c) Eq. 2 temperature from Table 2 widths (16–84 %)', loc='left', fontsize=10, color=FG)

    # ---------------- (d) [O III] width recovery
    ax = axs[1, 1]
    tw = [5.4, 20, 30, 45]
    ax.axhspan(43.1, 48.9, color=OIII_C, alpha=0.22, lw=0)
    label(ax, 1.0, 49.6, 'paper Table 2: 43.1–48.9 km/s', color=OIII_C, size=LS, va='bottom')
    hd = ax.plot([0, 52], [0, 52], color=FGDIM, lw=1.0, ls='--', label='recovered = true')
    for sn, col, off, lab in [('6.1', ACCENT, -0.6, 'integrated S/N 6.1 (paper P1, P2)'), ('15', FGDIM, 0.6, 'S/N 15'),
                              ('50', '#6f7b86', 1.4, 'S/N 50')]:
        med = [fm[f'SN{sn}_Wint{w}']['Wrec_pass_p16_50_84'] for w in tw]
        y = np.array([m[1] for m in med]); lo = y - np.array([m[0] for m in med]); hi = np.array([m[2] for m in med]) - y
        hd.append(ax.errorbar(np.array(tw) + off, y, yerr=[lo, hi], fmt='o-', ms=6, lw=1.4 if sn == '6.1' else 1.0,
                              color=col, capsize=2.5, label=lab, zorder=4 if sn == '6.1' else 3))
    pf = ', '.join(f"{fm[f'SN6.1_Wint{w}']['frac_pass_cut'] * 100:.0f}" for w in tw)
    label(ax, 51, 10.8, f'S/N 6.1 pass fraction: {pf} %\n(true FWHM 5.4, 20, 30, 45 km/s)', color=ACCENT, size=LS,
          ha='right', va='bottom')
    ax.set_xlim(0, 52); ax.set_ylim(0, 72)
    ax.set_xlabel('true intrinsic [O III] FWHM (km/s)'); ax.set_ylabel('recovered FWHM, fits passing cut (km/s)')
    ax.legend(handles=hd, loc='upper left', fontsize=8.0, framealpha=0.9)
    ax.grid(alpha=0.3)
    ax.set_title('(d) [O III] width at S/N 6: the cut selects ~36–47 km/s', loc='left', fontsize=10, color=FG)

    finish(fig, str(ROOT / 'figures/dispute_other_fig6.png'), dpi=200)


if __name__ == '__main__':
    make_fig()
    o = out
    print(json.dumps({'ratios': {k: {q: [round(x, 3) for x in o['red_ratios'][k][q]] for q in ('SII_sum_Ha', 'SII6716_Ha', 'NII6583_Ha')}
                                 for k in REG},
                      'delta': {k: {q: round(v, 2) for q, v in o['red_ratios'][k]['delta_vs_P4_sigma'].items()} for k in REG[:3]},
                      'dens5': o['sii_density']['five_level_range_all'], 'dens5TZ': o['sii_density']['five_level_range_TZ10_Omega'],
                      'eq5': o['sii_density']['eq5_range'], 'Rfor': o['sii_density']['R_for_ne'],
                      'LDL': o['sii_density']['R_LDL_5level'],
                      'em': {k: {q: round(v, 5) for q, v in o['sii_density']['emission_measure'][k].items()} for k in REG},
                      'T': {k: {q: (round(v, 1) if not isinstance(v, list) else [round(x) for x in v]) for q, v in o['temperature_turbulence'][k].items()} for k in REG},
                      'ab': {k: {q: (round(v, 3) if isinstance(v, float) else v) for q, v in o['abundances'][k].items()} for k in REG},
                      'eg': o['energetics'], 'od': o['oiii_doublet'], 'radio': o['radio']}, indent=1, default=jdef))
