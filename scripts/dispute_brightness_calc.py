"""dispute_brightness_calc.py -- numbers for dispute section 5 (Table 3 [O III] surface brightness).

Re-uses the verified conversion functions of scripts/verify_paper_energetics.py (imported, not modified) and reads
the verified project JSONs. Every number written to data/derived/dispute_brightness_calc.json is either
  (i) copied from a project file (key 'inputs', with its file and JSON path), or
  (ii) computed here from those inputs (keys 'conversion', 'ratios', 'multiplicative', 'additive', 'starlight',
       'fuv', 'scaling').
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_brightness_calc.py
"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
from verify_paper_energetics import mag_to_R, R_to_mag, conversion_cases  # noqa: E402

D = ROOT / 'data/derived'
OUT = D / 'dispute_brightness_calc.json'


def jload(name):
    return json.load(open(D / name))


def r(x, n=3):
    return float(round(float(x), n))


# ----------------------------------------------------------------------------------------------- inputs
T3 = {'P1': 23.60, 'P2': 23.87, 'P3': 23.85}            # Wu+26 Table 3 (p. 8), mag arcsec^-2, "[O III] band"
T3_err = {'P1': 0.028, 'P2': 0.035, 'P3': 0.039}
SQM = 21.5                                               # Wu+26 Sect. 2.1 (p. 4)
FILTER_FWHM = 60.0                                       # Wu+26 Sect. 2.1: CWL 500.7 nm, FWHM 6.0 nm
QHY_MAG = 24.3                                           # QHYCCD note 2024-07-12 (notes/discovery.md §5.4)
FESEN_LIMIT = [0.5, 5.0]                                 # Wu+26 Sect. 4.1 (p. 15), quoting Fesen+2024

nsns = jload('imaging_optical_nsns_measurements.json')
sdss_p = jload('verify_sdss_pregions.json')
sdss_r = jload('verify_sdss_results_base_r0.15.json')
desi = jload('verify_desi_summary.json')
galex = jload('verify_xrayfuv_galex_measure.json')
uvp = jload('verify_xrayfuv_uv_predictions.json')
ener = jload('verify_paper_energetics.json')['brightness_conversion']

inputs = {
    'table3_mag': T3, 'table3_mag_err': T3_err, 'sqm_mag': SQM, 'filter_fwhm_A': FILTER_FWHM,
    'qhyccd_mag': QHY_MAG, 'fesen_limit_R_as_quoted': FESEN_LIMIT,
    'nsns_onridge_excess_R': {k: r(nsns[k]['oiii']['excess_R'], 2) for k in ('P1', 'P2', 'P3')},
    'nsns_peak_excess_R': {k: r(nsns[k]['oiii']['peak_excess_R'], 2) for k in ('P1', 'P2', 'P3')},
    'nsns_P4_excess_R': r(nsns['P4']['oiii']['excess_R'], 2),
    'nsns_source': 'imaging_optical_nsns_measurements.json [Pk][oiii][excess_R / peak_excess_R]; on = |across|<=0.25 deg, 51.5" px',
    'sdss_P_mean_R': {k: r(sdss_p[k]['mean_R'], 2) for k in ('P1', 'P2', 'P3')},
    'sdss_P_err_R': {k: r(sdss_p[k]['err_R'], 2) for k in ('P1', 'P2', 'P3')},
    'sdss_brightest_fibre_R': r(sdss_p['brightest_on_fibres'][0]['I_R'], 2),
    'sdss_source': 'verify_sdss_pregions.json (mean within 0.3 deg, total-sky, f_eff 0.84/0.60; +-0.2 R zero point)',
    'sdss_south_onoff_R': r(sdss_r['regions']['south_dec_lt20']['OIII:OIII5007']['I_R'], 2),
    'sdss_south_onoff_err_stat_R': r(sdss_r['regions']['south_dec_lt20']['OIII:OIII5007']['err_boot'], 2),
    'sdss_south_onoff_err_placebo_R': 0.07,               # notes/verify_sdss.md §1 (placebo term)
    'sdss_sky_continuum_5070_ABmag': {k: r(v, 2) for k, v in sdss_r['sb5070_median_cal'].items()},
    'desi_south_median_R': r(desi['brightness']['south_median']['gfa_chromatic'], 2),
    'desi_south_variants_R': [r(desi['brightness']['south_median'][k], 2) for k in ('gfa_chromatic', 'gfa_halfchrom', 'gfa_achromatic')],
    'desi_sys_frac': 0.08,                                # notes/verify_desi.md (GFA acceptance, chromatic model)
    'desi_nodes_R': {k: v for k, v in desi['brightness']['nodes_gfa'].items()},
    'sdss_imaging_R': {'g_PSP': [2.31, 0.38], 'g_FramesSub': [2.39, 0.46], 'g_Frames': [2.55, 0.53],
                       'all_confounders': [1.65, 0.31], 'source': 'notes/verify_precovery.md §2.1 table'},
    'lamost_mrsn_scaled_R': {'APO_ratio': [2.8, 3.1, 3.2], 'APO_16_84': [1.2, 8.6], 'UVES_ratio': [14.2, 15.5, 16.1],
                             'source': 'notes/spectra_sdss_lamost.md §4.1 (wave 1; indeterminate)'},
    'lamost_lrs_resid_R': [0.88, 0.13],
    'galex_fuv_beta_south_CU_per_NSNS_R': [r(galex['FUV']['beta_OIII_south_dec<17.6']['beta_CU_per_R'], 1),
                                           r(galex['FUV']['beta_OIII_south_dec<17.6']['null_std'], 1)],
    'galex_fuv_beta_all_CU_per_NSNS_R': [r(galex['FUV']['beta_OIII_all']['beta_CU_per_R'], 1),
                                         r(galex['FUV']['beta_OIII_all']['null_std'], 1)],
    'galex_source': 'verify_xrayfuv_galex_measure.json [FUV][beta_OIII_*]',
}

# ----------------------------------------------------------------------------------------------- conversion
conv = {}
conv['R_per_cgs_note'] = '1 R at 5006.84 A = %.4e erg s-1 cm-2 arcsec-2 (verify_paper_energetics.R_5007)' % ener['R_per_cgs']
conv['AB_W60'] = {k: r(mag_to_R(m), 2) for k, m in T3.items()}
conv['by_width_A'] = {str(W): {k: r(mag_to_R(m, W), 2) for k, m in T3.items()} for W in (30, 40, 50, 60, 70)}
conv['AB_Vega_pm0.1'] = {k: [r(mag_to_R(m, ab_offset=0.1), 2), r(mag_to_R(m, ab_offset=-0.1), 2)] for k, m in T3.items()}
conv['4959_fully_transmitted'] = {k: r(mag_to_R(m) / (1 + 1 / 2.98), 2) for k, m in T3.items()}
vals = []
for W in (50, 60, 70):
    for off in (-0.1, 0.0, 0.1):
        for leak in (1.0, 1 / (1 + 1 / 2.98)):
            vals += [mag_to_R(m, W, off) * leak for m in T3.values()]
conv['all_variants_range_R'] = [r(min(vals), 2), r(max(vals), 2)]
conv['mag_for_R_W60'] = {str(I): r(R_to_mag(I), 2) for I in (1.4, 1.62, 2.0, 2.5, 3.0)}
conv['qhyccd_R_W30_W60'] = [r(mag_to_R(QHY_MAG, 30), 2), r(mag_to_R(QHY_MAG, 60), 2)]
conv['sky_R_equiv_W60'] = {str(s): r(mag_to_R(s), 1) for s in (21.0, 21.1, 21.5, 21.75, 22.0, 22.5)}
fc = conversion_cases()


def rng(sub):
    v = [c['true_over_paper'] for k, c in fc.items() if sub(k)]
    return [r(min(v), 2), r(max(v), 2)]


conv['filter_true_over_paper'] = {
    'f0.95_converging_untuned': rng(lambda k: 'f/0.95' in k),
    'f2.8_converging': rng(lambda k: 'f/2.8' in k),
    'collimated_0_5deg': rng(lambda k: 'collimated_field_0deg' in k or 'collimated_field_5deg' in k),
    'collimated_10_15deg': rng(lambda k: any(f'collimated_field_{a}deg' in k for a in (10, 13, 15))),
    'model': '3-cavity 60 A FWHM, n* 1.8-2.05, CWL 5007/5015/5030 A; stars calibrate int(T), line sees T(5007)+T(4959)/2.98',
}

# ----------------------------------------------------------------------------------------------- ratios Table3 / independent
I3 = conv['AB_W60']
P = ('P1', 'P2', 'P3')
ratios = {
    'vs_nsns_onridge': {k: r(I3[k] / inputs['nsns_onridge_excess_R'][k], 2) for k in P},
    'vs_nsns_peak': {k: r(I3[k] / inputs['nsns_peak_excess_R'][k], 2) for k in P},
    'vs_sdss_P': {k: r(I3[k] / inputs['sdss_P_mean_R'][k], 2) for k in P},
    'vs_sdss_south_onoff': [r(min(I3.values()) / inputs['sdss_south_onoff_R'], 2), r(max(I3.values()) / inputs['sdss_south_onoff_R'], 2)],
    'vs_desi_south': [r(min(I3.values()) / inputs['desi_south_median_R'], 2), r(max(I3.values()) / inputs['desi_south_median_R'], 2)],
    'vs_sdss_imaging_2.31_2.55': [r(min(I3.values()) / 2.55, 2), r(max(I3.values()) / 2.31, 2)],
    'vs_qhyccd': [r(min(I3.values()) / conv['qhyccd_R_W30_W60'][1], 2), r(max(I3.values()) / conv['qhyccd_R_W30_W60'][0], 2)],
    'vs_brightest_single_point': {'sdss_fibre_3.65': [r(min(I3.values()) / 3.65, 2), r(max(I3.values()) / 3.65, 2)],
                                  'desi_W@18_3.88': [r(min(I3.values()) / 3.88, 2), r(max(I3.values()) / 3.88, 2)]},
}
allr = list(ratios['vs_nsns_onridge'].values()) + list(ratios['vs_sdss_P'].values()) + ratios['vs_desi_south']
ratios['region_matched_range'] = [r(min(allr), 2), r(max(allr), 2)]
ratios['region_matched_range_excl_P1'] = [r(min(ratios['vs_nsns_onridge']['P2'], ratios['vs_nsns_onridge']['P3'], ratios['vs_sdss_P']['P2'], ratios['vs_sdss_P']['P3']), 2),
                                          r(max(ratios['vs_nsns_onridge']['P2'], ratios['vs_nsns_onridge']['P3'], ratios['vs_sdss_P']['P2'], ratios['vs_sdss_P']['P3']), 2)]
p23 = [ratios['vs_nsns_onridge'][k] for k in ('P2', 'P3')] + [ratios['vs_sdss_P'][k] for k in ('P2', 'P3')] + \
      [I3[k] / inputs['desi_south_median_R'] for k in ('P2', 'P3')]
ratios['region_matched_range_P2P3_incl_desi'] = [r(min(p23), 2), r(max(p23), 2)]
ratios['delta_mag_region_matched'] = [r(2.5 * np.log10(min(allr)), 2), r(2.5 * np.log10(max(allr)), 2)]
ratios['delta_mag_over_quoted_sigma'] = [r(2.5 * np.log10(min(allr)) / max(T3_err.values()), 0),
                                         r(2.5 * np.log10(max(allr)) / min(T3_err.values()), 0)]
req = [1 / x for x in allr]
ratios['required_factor_true_over_T3'] = [r(min(req), 3), r(max(req), 3)]
ratios['fesen_logmid_R'] = r(np.sqrt(FESEN_LIMIT[0] * FESEN_LIMIT[1]), 2)

# ----------------------------------------------------------------------------------------------- multiplicative stack
f_w50 = 50.0 / 60.0                                      # R scales with W: the narrower 50 A band lowers R
f_vega = 10 ** (-0.4 * 0.1)
f_4959 = 1 / (1 + 1 / 2.98)
f_star = 10 ** (-0.4 * 0.3)                              # Wu+26 Sect. 3.1: P1 stellar contamination (verify_paper: ~0.3 mag)
stack = f_w50 * f_vega * f_4959
mult = {
    'factors_true_over_T3': {'AB_Vega_pm0.1mag': [r(f_vega, 3), r(1 / f_vega, 3)], 'width_50_70A': [r(f_w50, 3), r(70 / 60, 3)],
                             '4959_leak_0_to_full': [r(f_4959, 3), 1.0], 'P1_star_0.3mag': r(f_star, 3)},
    'stack_lowering_P2P3': r(stack, 3), 'stack_lowering_P1': r(stack * f_star, 3),
    'stacked_R': {'P1': r(I3['P1'] * stack * f_star, 2), 'P2': r(I3['P2'] * stack, 2), 'P3': r(I3['P3'] * stack, 2)},
    'note': 'generous: 50 A effective width, +0.1 mag Vega/AB, full 4959 transmission and P1 0.3 mag applied simultaneously; '
            'these cannot all hold for one 60-A filter centred on 5007',
}
gap = {}
for k in P:
    st = mult['stacked_R'][k]
    gap[k] = [r(st / inputs['nsns_onridge_excess_R'][k], 2), r(st / inputs['sdss_P_mean_R'][k], 2)]
mult['stacked_over_independent_P_nsns_sdss'] = gap
g = [x for v in gap.values() for x in v]
mult['stacked_over_independent_range'] = [r(min(g), 2), r(max(g), 2)]

# ----------------------------------------------------------------------------------------------- additive requirement
add = {}
dI = {}
for k in P:
    dI[k] = {'vs_nsns': r(I3[k] - inputs['nsns_onridge_excess_R'][k], 2), 'vs_sdss': r(I3[k] - inputs['sdss_P_mean_R'][k], 2)}
add['delta_I_R'] = dI
dv = [x for v in dI.values() for x in v.values()]
add['delta_I_range_R'] = [r(min(dv), 2), r(max(dv), 2)]
for s in (21.1, 21.5, 22.0):
    S = mag_to_R(s)
    add[f'frac_of_sky_{s}'] = [r(min(dv) / S, 3), r(max(dv) / S, 3)]
add['frac_of_sky_range_21.5_22.0'] = [add['frac_of_sky_21.5'][0], add['frac_of_sky_22.0'][1]]
add['frac_of_sky_range_21.1_22.0'] = [add['frac_of_sky_21.1'][0], add['frac_of_sky_22.0'][1]]
# generic map for the figure: required fraction for I_true grid and sky grid (Table 3 P1-P3 each)
add['map_def'] = 'frac = (I_T3 - I_true) / I_sky(W60), I_sky = mag_to_R(sky_mag, 60)'

# ----------------------------------------------------------------------------------------------- starlight (Gaia DR3 cones)
gdir = D / 'verify_paper_raw/gaia_isl'
area = np.pi * (0.25 * 3600) ** 2
star = {}
for k in ('P1', 'P2', 'P3', 'P4'):
    a = np.genfromtxt(gdir / f'isl_{k}_r0.25.csv', delimiter=',', names=True)
    row = {}
    for lab, sel in (('all', np.ones(len(a), bool)), ('G_lt_14', a['gbin'] < 14), ('G_ge_14', a['gbin'] >= 14),
                     ('G_ge_15', a['gbin'] >= 15), ('G_ge_16', a['gbin'] >= 16)):
        f = a['fsum'][sel].sum()
        m = -2.5 * np.log10(f / area)
        row[lab] = {'SB_Gmag': r(m, 2), 'R_eq_W60': r(mag_to_R(m), 2), 'R_eq_W60_colour0.3': r(mag_to_R(m + 0.3), 2)}
    star[k] = row
star['note'] = ('Gaia DR3 G (Vega) integrated starlight in 0.25-deg cones (verify_paper_raw/gaia_isl); R-equivalent treats '
                'G as an AB magnitude in a 60 A band at 5007 A (upper value) or 0.3 mag fainter (colour term; lower value)')
dT3 = {'P1-P2': r(I3['P1'] - I3['P2'], 2), 'P1-P3': r(I3['P1'] - I3['P3'], 2)}
dS = {'P1-P2': r(inputs['sdss_P_mean_R']['P1'] - inputs['sdss_P_mean_R']['P2'], 2),
      'P1-P3': r(inputs['sdss_P_mean_R']['P1'] - inputs['sdss_P_mean_R']['P3'], 2)}
dN = {'P1-P2': r(inputs['nsns_onridge_excess_R']['P1'] - inputs['nsns_onridge_excess_R']['P2'], 2),
      'P1-P3': r(inputs['nsns_onridge_excess_R']['P1'] - inputs['nsns_onridge_excess_R']['P3'], 2)}
p1x = [dT3[k] - dS[k] for k in dT3] + [dT3[k] - dN[k] for k in dT3]
dstar_all = [star['P1']['all']['R_eq_W60'] - star[k]['all']['R_eq_W60'] for k in ('P2', 'P3')]
dstar_all_c = [star['P1']['all']['R_eq_W60_colour0.3'] - star[k]['all']['R_eq_W60_colour0.3'] for k in ('P2', 'P3')]
star['P1_differential'] = {
    'table3_P1_minus_Pk_R': dT3, 'sdss_P1_minus_Pk_R': dS, 'nsns_P1_minus_Pk_R': dN,
    'P1_specific_excess_R': [r(min(p1x), 2), r(max(p1x), 2)],
    'starlight_all_P1_minus_P2P3_R_eq': [r(min(dstar_all_c), 1), r(max(dstar_all), 1)],
    'residual_fraction_needed': [r(min(p1x) / max(dstar_all), 3), r(max(p1x) / min(dstar_all_c), 3)],
}
ge14 = [star[k]['G_ge_14'][c] for k in P for c in ('R_eq_W60', 'R_eq_W60_colour0.3')]
star['unresolved_G_ge_14_R_eq_range_P1P3'] = [r(min(ge14), 2), r(max(ge14), 2)]
star['unresolved_G_ge_14_frac_of_required'] = [r(min(ge14) / add['delta_I_range_R'][1], 2), r(max(ge14) / add['delta_I_range_R'][0], 2)]

# ----------------------------------------------------------------------------------------------- FUV
beta, ebeta = inputs['galex_fuv_beta_south_CU_per_NSNS_R']
fuv = {}
kfac = {k: I3[k] / inputs['nsns_onridge_excess_R'][k] for k in P}
fuv['nsns_scale_factor_if_T3_true'] = {k: r(v, 2) for k, v in kfac.items()}
fuv['FUV_CU_per_trueR_if_T3'] = {k: r(beta / v, 1) for k, v in kfac.items()}
fuv['FUV_CU_per_trueR_if_T3_1sigma_high'] = {k: r((beta + ebeta) / v, 1) for k, v in kfac.items()}
floors = {}
for m, v in uvp['models'].items():
    floors[m] = r(v['per_R']['ebv0.03_civsupp1.0']['FUV']['non_resonant'], 1)
fuv['nonresonant_floor_CU_per_R_ebv0.03'] = floors
fl_min = min(floors.values())
fl_complete = [min(v for m, v in floors.items() if not m.startswith('R88')), max(v for m, v in floors.items() if not m.startswith('R88'))]
fuv['floor_min_all_models'] = fl_min
fuv['floor_complete_range'] = fl_complete
fv = list(fuv['FUV_CU_per_trueR_if_T3'].values())
fv_hi = list(fuv['FUV_CU_per_trueR_if_T3_1sigma_high'].values())
fuv['T3_implied_range'] = [min(fv), max(fv)]
fuv['T3_implied_below_lowest_floor_by'] = [r(fl_min / max(fv), 2), r(fl_min / min(fv), 2)]
fuv['T3_implied_1sigma_high_below_lowest_floor_by'] = r(fl_min / max(fv_hi), 2)
# reddening that would lower the lowest (R88) floor to the Table-3-implied ratio: log-linear in E(B-V) through the
# tabulated E(B-V) = 0, 0.01, 0.03 values (extrapolation, model only)
r88 = uvp['models']['R88_SE']['per_R']
eb = np.array([0.0, 0.01, 0.03])
fl = np.array([r88[f'ebv{e:g}_civsupp1.0']['FUV']['non_resonant'] if e else r88['ebv0.0_civsupp1.0']['FUV']['non_resonant'] for e in eb])
kext = -np.polyfit(eb, np.log(fl), 1)[0]
fuv['R88_floor_vs_ebv'] = {'ebv': eb.tolist(), 'floor': [r(x, 1) for x in fl], 'dlnF_dEBV': r(kext, 2)}
fuv['EBV_needed_to_reach_T3_implied'] = [r(np.log(fl[0] / max(fv)) / kext, 2), r(np.log(fl[0] / min(fv)) / kext, 2)]
fuv['measured_CU_per_R'] = {'south_NSNS_scale': [beta, ebeta],
                            'south_DESI_scale': r(beta * inputs['nsns_onridge_excess_R']['P1'] / inputs['desi_south_median_R'], 0)}
fuv['note'] = ('beta from template regression of GALEX FUV on NSNS [O III] (south, Dec<17.6); if Table 3 were correct the NSNS '
               'scale would be low by I_T3/I_NSNS per region, so FUV per true R = beta / that factor. Floors: non-resonant '
               'FUV lines + two-photon (C IV, Si IV/O IV] excluded), E(B-V)=0.03, from verify_xrayfuv_uv_predictions.json '
               '(R79, CR85, DS96, R88 tables).')

# ----------------------------------------------------------------------------------------------- scaling of derived quantities
lo, hi = ratios['region_matched_range']
scaling = {'EM_and_luminosity_overestimate': [lo, hi], 'ne_from_EM_overestimate_sqrt': [r(np.sqrt(lo), 2), r(np.sqrt(hi), 2)],
           'E_if_n0_tied_to_I_overestimate_pow1.16': [r(lo ** 1.16, 1), r(hi ** 1.16, 1)]}

res = dict(inputs=inputs, conversion=conv, ratios=ratios, multiplicative=mult, additive=add, starlight=star, fuv=fuv,
           scaling=scaling, provenance='scripts/dispute_brightness_calc.py; run 2026-11-17')
json.dump(res, open(OUT, 'w'), indent=1)

if __name__ == '__main__':
    for k in ('conversion', 'ratios', 'multiplicative', 'additive', 'fuv', 'scaling'):
        print('==', k)
        for kk, vv in res[k].items():
            print('  ', kk, vv)
    print('== starlight')
    for k in ('P1', 'P2', 'P3', 'P4'):
        print('  ', k, {kk: (vv['SB_Gmag'], vv['R_eq_W60'], vv['R_eq_W60_colour0.3']) for kk, vv in star[k].items()})
    for kk in ('P1_differential', 'unresolved_G_ge_14_R_eq_range_P1P3', 'unresolved_G_ge_14_frac_of_required'):
        print('  ', kk, star[kk])
    print('wrote', OUT)
