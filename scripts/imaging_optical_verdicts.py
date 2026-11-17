"""Collect the per-survey precovery verdicts of the imaging_optical topic into one JSON
(data/derived/imaging_optical_verdicts.json) from the derived measurement files."""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
from imaging_optical_common import DER, save_json

L = lambda f: json.loads((DER / f).read_text())
sd = L('imaging_optical_sdss_sky.json')['cases']; pl = L('imaging_optical_sdss_placebo_w4.0.json')
ds = L('imaging_optical_dss_summary.json'); ps = L('imaging_optical_ps1_summary.json')
ha = L('imaging_optical_halpha_measurements.json'); ns = L('imaging_optical_nsns_measurements.json')
c = sd['g_vs_OIII_ridge_nostars']['core_dec12_20p5']
out = {
 'units': 'A = broadband response / response expected if the NSNS DR0.2 [OIII] map were exactly calibrated (dimensionless)',
 'SDSS_DR18_sky_g_2002-2007': dict(verdict='detected (statistical, template-correlation)', A=c['A_over_expected'], err=c['err_over_expected'],
     n_runs=c['n_runs'], runs_positive=sd['g_vs_OIII_ridge_nostars']['core_n_positive'],
     window4deg_A=pl['zero_shift'], null_mean=pl['placebo_mean_absshift_ge2'], null_std=pl['placebo_std_absshift_ge2'],
     null_tests={'r_vs_OIII': sd['r_vs_OIII_ridge']['core_dec12_20p5']['A_over_expected'],
                 'r_vs_OIII_err': sd['r_vs_OIII_ridge']['core_dec12_20p5']['err_over_expected'],
                 'g_vs_Halpha': sd['g_vs_Halpha']['core_dec12_20p5']['A_over_expected'],
                 'g_vs_Halpha_err': sd['g_vs_Halpha']['core_dec12_20p5']['err_over_expected']}),
 'PS1_DR2_stacks': dict(verdict='not detected; insensitive by construction (sky subtraction removes >~few-arcmin emission; Medusa PN test)',
     strips={k: (round(v['A_over_expected'], 2), round(v['err_over_expected'], 2)) for k, v in ps['strips'].items()}),
 'DSS2_J_POSSII': dict(verdict='not detected', A=ds['DSS2_Blue']['combined_d13_d19']['A'], err=ds['DSS2_Blue']['combined_d13_d19']['err'],
     upper3sigma=ds['DSS2_Blue']['upper3sigma_trueR_per_NSNS_R_d13_d19'], pct_of_sky_per_R=ds['DSS2_Blue']['pct_of_sky_per_R']),
 'DSS1_O_POSSI': dict(verdict='not detected (5007 outside 103aO band)', A=ds['DSS1_Blue']['combined_d13_d19']['A'], err=ds['DSS1_Blue']['combined_d13_d19']['err']),
 'DSS2_F_and_DSS1_E (controls)': dict(verdict='no [OIII]-shaped signal (as expected)', F=(ds['DSS2_Red']['combined_d13_d19']['A'], ds['DSS2_Red']['combined_d13_d19']['err']),
     E=(ds['DSS1_Red']['combined_d13_d19']['A'], ds['DSS1_Red']['combined_d13_d19']['err'])),
 'SHASSA_Halpha_1997-2000': dict(verdict='not detected (Dec 12.5-15)', excess_R={k: (round(v['Ha_excess_R'], 2), round(max(v['sys_EW_R'], v['strip_scatter_R']), 2)) for k, v in ha['shassa'].items()}),
 'WHAM_SS_1997-98': dict(verdict='not detected along [OIII] south; +0.6-0.9 R (1.5-1.7 sigma) at Dec 17-21 where NSNS shows Halpha; north patch beam +1.5 R (v -10..+30)',
     v_m10_p30=[(r['dec_bin'], round(r['excess_R'], 2), round(r['off_scatter_R'], 2)) for r in ha['wham']['v_m10_p30']]),
 'Finkbeiner2003': dict(verdict='not detected in [OIII]-only south; +0.8-0.9 R at Dec 17-21 (WHAM-resolution, systematics ~0.7-1.6 R)',
     excess_R={k: round(v['fink_excess_R'], 2) for k, v in ha['finkbeiner_vs_ridge'].items()}),
 'VTSS': dict(verdict='no coverage (SkyView: no data for the field)'),
 'NSNS_DR0.2 (reference, 2022-12..2025-04)': dict(verdict='detected (discovery-era data, not precovery)',
     oiii_on_ridge_excess_R={k: round(ns[k]['oiii']['excess_R'], 2) for k in ('P3', 'P1', 'P2', 'south_bright', 'mid_17_19', 'north_arc_21_23', 'north_arc_24_27')},
     oiii_peak_excess_R={k: round(ns[k]['oiii']['peak_excess_R'], 2) for k in ('P3', 'P1', 'P2')},
     Ha_patch_north={ch: round(ns['Ha_patch_north'][ch]['excess_R'], 2) for ch in ('oiii', 'halpha', 'sii')},
     P4={ch: round(ns['P4'][ch]['excess_R'], 2) for ch in ('oiii', 'halpha', 'sii')}),
 'AstroBin_pre-2023-11': dict(verdict='no narrowband image covering the ridge found; 169 wide broadband images cover it (expected contrast 0.3-0.5% of sky, invisible)'),
}
save_json(out, DER / 'imaging_optical_verdicts.json')
print(json.dumps(out, indent=1, default=float)[:3000])
