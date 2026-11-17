"""w3_oii_summary.py -- (wave 3, tag oii) key numbers of the [O II] doublet forensic test in one JSON.

Reads the w3_oii_* products and adds an empirical error-inflation factor for the DESI regression ratios from the
consistency of DISJOINT fibre subsets (independent noise, same filament): pull_k = (R_a - R_b) / sqrt(s_a^2 + s_b^2)
with bootstrap s = (p84 - p16) / 2; inflation = rms(pull).
Output: data/derived/w3_oii_summary.json
"""
import json
import numpy as np

D = '/Users/mickey/astronomy/silk-nebula/data/derived/'


def main():
    F = json.load(open(D + 'w3_oii_desi_fits.json'))
    FR = json.load(open(D + 'w3_oii_desi_fits_robust.json'))
    S = json.load(open(D + 'w3_oii_sdss.json'))
    W = json.load(open(D + 'w3_oii_sdss_wim.json'))
    WC = json.load(open(D + 'w3_oii_sdss_wavecal.json'))['summary']
    SIM = json.load(open(D + 'w3_oii_sim.json'))
    A = json.load(open(D + 'w3_oii_atomic.json'))
    out = {}

    def r(f, sn, model='kin'):
        o = f[f'geo|{sn}|F']['o2'][model]
        return o['ratio'], 0.5 * (o['p84'] - o['p16'])
    pairs = [('qso_hiz', 'lowz_gal'), ('dark', 'bright'), ('south', 'north'), ('Zpos', 'Zneg'), ('w1_only', 'ext_only')]
    pulls = {}
    for a, b in pairs:
        ra, sa = r(F, a)
        rb, sb = r(F, b)
        pulls[f'{a}-{b}'] = dict(Ra=ra, Rb=rb, pull=(ra - rb) / np.hypot(sa, sb))
    pv = np.array([p['pull'] for k, p in pulls.items() if 'ext_only' not in k])   # ext_only bootstrap degenerate
    out['desi_disjoint_pulls'] = pulls
    out['desi_inflation_factor'] = float(np.sqrt(np.mean(pv ** 2)))
    R, s = r(F, 'all')
    infl = max(out['desi_inflation_factor'], 1.5)
    out['desi_all'] = dict(ratio=R, sig_boot=s, sig_inflated=s * infl,
                           sigmas_below_LDL_1p46=(1.46 - R) / (s * infl),
                           models={k: F['geo|all|F']['o2'][k]['ratio'] for k in
                                   ('w2', 'kin', 'kin_s-10', 'kin_s+10', 'kin_v-5', 'kin_v+5', 'vfree', 'vsfree')},
                           robust=FR['geo|all|F']['o2']['kin']['ratio'],
                           profile68=F['geo|all|F']['o2']['profile']['r68'],
                           profile_dchi2_1p46=F['geo|all|F']['o2']['profile']['dchi2_at_1p46'])
    out['desi_subsets_kin'] = {k.split('|')[1]: v['o2']['kin']['ratio'] for k, v in F.items()
                               if k.startswith('geo|') and k.endswith('|F')}
    out['desi_subsets_kin_robust'] = {k.split('|')[1]: v['o2']['kin']['ratio'] for k, v in FR.items()
                                      if k.startswith('geo|') and k.endswith('|F')}
    out['desi_o3_kinematics'] = {k.split('|')[1]: dict(v=v['o3']['v'], sigma_int_kms=v['o3']['sigma_int_kms'])
                                 for k, v in F.items() if k.startswith('geo|') and k.endswith('|F')}
    out['desi_Ha_component'] = {k.split('|')[1]: dict(oii=v['o2']['kin']['ratio'], oii_p16=v['o2']['kin']['p16'],
                                                      oii_p84=v['o2']['kin']['p84'],
                                                      sii=v['s2']['kin']['ratio'] if 's2' in v else None)
                                for k, v in F.items() if k.startswith('ha|') and k.endswith('|F')}
    out['sdss_onoff'] = {rn: {t: dict(ratio=v['o2'][t]['ratio'], p16=v['o2'][t]['p16'], p84=v['o2'][t]['p84'],
                                      v=v['o2'][t]['v'], s_over_wdisp=v['o2'][t].get('s_over_wdisp'))
                              for t in ('fixed', 'vfree', 'vsfree', 'eps_fixed')} | {
        'sii_fixed': dict(ratio=v['s2']['fixed']['ratio'], p16=v['s2']['fixed']['p16'], p84=v['s2']['fixed']['p84'],
                          v_free=v['s2']['vfree']['v'])} for rn, v in S['onoff'].items()}
    out['sdss_wim_control'] = {m: {t: dict(ratio=W[m]['o2'][t]['ratio'], p16=W[m]['o2'][t]['p16'],
                                           p84=W[m]['o2'][t]['p84'], v=W[m]['o2'][t]['v'],
                                           s_over_wdisp=W[m]['o2'][t]['s_over_wdisp'])
                                   for t in ('fixed', 'vfree', 'sfree', 'vsfree')} | {
        'sii_fixed': W[m]['s2']['fixed']['ratio'], 'sii_p16': W[m]['s2']['fixed']['p16'],
        'sii_p84': W[m]['s2']['fixed']['p84'], 'ha_slope_cal': W[m]['ha_slope_cal']['value'],
        'oii_total_per_R_Ha_cal': W[m]['o2']['fixed']['total']} for m in ('between', 'within')}
    out['boss_bluecal'] = WC
    out['simulation_key'] = {inst: dict(sep_over_sigma=SIM[inst]['separation_over_sigma'],
                                        wave2_setup=SIM[inst]['Rtrue_1.46']['wave2_setup'],
                                        v_plus10=SIM[inst]['Rtrue_1.46']['v_fit_minus_true_+10'],
                                        lsf_true_over_fit_1p20=SIM[inst]['Rtrue_1.46']['lsf_true_over_fit_1.20'],
                                        h3_minus0p10=SIM[inst]['Rtrue_1.46']['lsf_h3_-0.10'],
                                        blend_dv20=SIM[inst]['Rtrue_1.46']['blend_two_strands_dv20'],
                                        noise_peakSN17=SIM[inst]['Rtrue_1.46']['noise']['peakSN17'])
                           for inst in ('DESI', 'BOSS')}
    sp = SIM['DESI_specter_replica_Rtrue1.46']
    out['simulation_key']['DESI_extraction_phase_range'] = [min(v['ratio_fit_sigmaR'] for v in sp.values()),
                                                            max(v['ratio_fit_sigmaR'] for v in sp.values())]
    out['atomic'] = dict(low_density_limits=A['low_density_limits'], DESI_all_inversion=A['ne_inversion']['DESI_all'],
                         DESI_all_ne_2sigma_upper=A['DESI_all_ne_2sigma_upper'], sii=A['sii'],
                         desi_7325_refit=A['desi_7325_refit'], T_OII_for_ratio=A['T_OII_for_ratio'],
                         implications=A['implications'])
    json.dump(out, open(D + 'w3_oii_summary.json', 'w'), indent=1, default=float)
    print(json.dumps({k: out[k] for k in ('desi_disjoint_pulls', 'desi_inflation_factor', 'desi_all')}, indent=1,
                     default=float))


if __name__ == '__main__':
    main()
