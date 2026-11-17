"""w3_northsouth_summary.py -- collate the wave-3 northsouth results into one JSON for the synthesis / judges.
Reads data/derived/w3_northsouth_{desi,desi_velcheck,sdss,rm,geom,shock}.json; writes
data/derived/w3_northsouth_summary.json. No new measurements.
"""
import json

DER = '/Users/mickey/astronomy/silk-nebula/data/derived'
L = lambda f: json.load(open(f'{DER}/{f}'))  # noqa: E731
dj, vc, sj, rj, gj, kj = (L('w3_northsouth_desi.json'), L('w3_northsouth_desi_velcheck.json'),
                          L('w3_northsouth_sdss.json'), L('w3_northsouth_rm.json'), L('w3_northsouth_geom.json'),
                          L('w3_northsouth_shock.json'))


def r2(x):
    return None if x is None else round(float(x), 3)


out = {'desi_bins': [], 'sdss': {}, 'desi_zones': {}}
A = dj['A_all_strands']['quantities']
for k, (lo, hi) in enumerate(zip(dj['edges'][:-1], dj['edges'][1:])):
    out['desi_bins'].append(dict(dec=[lo, hi], OIII_R=r2(A['O3R']['slope'][k]), OIII_err=r2(A['O3R']['slope_err_infl'][k]),
                                 OII_OIII=r2(A['OII3727sum']['ratio'][k]),
                                 OII_OIII_err=r2(A['OII3727sum']['ratio_err_infl'][k]),
                                 Hb_OIII=r2(A['Hb']['ratio'][k]), Hb_OIII_err=r2(A['Hb']['ratio_err_infl'][k])))
for name, key in (('halves', 'D_halves'), ('three', 'E_three_zones')):
    q = dj[key]['quantities']
    out['desi_zones'][name] = dict(OIII_R=[r2(x) for x in q['O3R']['slope']],
                                   OIII_err=[r2(x) for x in q['O3R']['slope_err_infl']],
                                   OII_OIII=[r2(x) for x in q['OII3727sum']['ratio']],
                                   OII_OIII_err=[r2(x) for x in q['OII3727sum']['ratio_err_infl']],
                                   Hb_OIII=[r2(x) for x in q['Hb']['ratio']],
                                   Hb_OIII_err=[r2(x) for x in q['Hb']['ratio_err_infl']],
                                   OIII_Hb_2sig_lower=[r2(1 / (h + 2 * e)) for h, e in
                                                       zip(q['Hb']['ratio'], q['Hb']['ratio_err_infl'])])
out['desi_placebo_inflation'] = {k: r2(v['robust']) for k, v in dj['placebo']['inflation'].items()}
for k, v in sj['regions'].items():
    if v is None:
        continue
    rr = v['ratios'].get('OII/OIII5007_energy_BOSS', {})
    vv = v['velocity_OIII']
    out['sdss'][k] = dict(OIII_R=r2(v['OIII:OIII5007']['I_R']), OIII_err=r2(v['OIII:OIII5007']['err_boot']),
                          OII_OIII_energy=r2(rr.get('value')), OII_OIII_p16_p84=[r2(rr.get('p16')), r2(rr.get('p84'))],
                          OII_R_photon=r2(rr.get('OII_R')), Hb_R=r2(v['Hb:Hb']['I_R']), Hb_err=r2(v['Hb:Hb']['err_boot']),
                          OIII_Hb_photon_2sig_lower=r2(v['ratios']['OIII5007/Hb_lower2sig_boot']),
                          v_helio=r2(vv.get('v_corr')), v_err=r2(vv.get('ev_boot')), n_groups=v['n_groups'],
                          n_on=v['n_on_fibres'])
out['desi_velocity'] = dict(
    bins=[dict(bin=b['bin'], v_helio=r2(b['desi_v_helio']), err=r2(b['desi_err_formal_infl']), v_lsr=r2(b['desi_v_lsr']),
               HI_vpeak=r2(b['HI']['vpeak']) if b['HI'] else None) for b in gj['velocity_bins']],
    zones={k: dict(helio_range=[r2(x) for x in z['desi_helio_range']], lsr_range=[r2(x) for x in z['desi_lsr_range']],
                   sdss_helio=r2(z['sdss_helio']), sdss_err=r2(z['sdss_err']), sdss_lsr=r2(z['sdss_lsr']))
           for k, z in gj['velocity_zones'].items()},
    dv_north_minus_south_helio_halves={k: [r2(v['p50']), r2(v['sig'])] for k, v in
                                       vc['halves']['diff_last_minus_first'].items()},
    dv_north_minus_south_helio_zones={k: [r2(v['p50']), r2(v['sig'])] for k, v in
                                      vc['three']['diff_last_minus_first'].items()},
    split_scan=[dict(split=s['split'], south=r2(s['south']['v_helio']), north=r2(s['north']['v_helio']),
                     north_err=r2(s['north']['v_err_formal_infl'])) for s in dj['vel_split_scan']],
    lsr_corr_south_north=[r2(x) for x in gj['lsr_corr_south_north']],
    dv_lsr_range=[r2(x) for x in gj['models']['observed_dv_lsr_range']])
out['rm'] = dict(true={h: dict(S=r2(v['rm']['S']), boot=r2(v['rm']['err_boot']), S_resid=r2(v['resid']['S']),
                               step_loc=r2(v['l1_step_location_deg'])) for h, v in rj['true'].items()},
                 nulls={m: {h: dict(sd=r2(s['null_sd']), p_two=r2(s['p_abs_ge_true']), p_one=r2(s['p_ge_true']),
                                    p_resid=r2(s['resid_p_abs_ge_true']), z=r2(s['z_true']), n=s['n'])
                            for h, s in v.items()} for m, v in rj['nulls'].items()},
                 smooth_gradient_pred=r2(rj['smooth_gradient']['quadratic']['predicted_step']),
                 coord_check=rj['check_vs_wave1'])
c = gj['centres']
out['geometry'] = dict(
    theta_zones={k: {z: [r2(v['mean']), r2(v['rms'])] for z, v in c[k]['zones'].items()} for k in c},
    crossings={k: {kk: r2(vv) for kk, vv in c[k]['dec_where_theta_exceeds'].items()} for k in c},
    circle_fits={k: dict(ra=r2(v['centre_ra']), dec=r2(v['centre_dec']), R=r2(v['R_deg']), R16_84=v['boot_R_p16_p84'],
                         frac_R_lt_15=v['frac_boot_R_lt_15'], chi2_free=r2(v['chi2']),
                         chi2_knies2_fixed=r2(v['knies2_fixed_centre']['chi2']),
                         rms_curv_fixed=r2(v['curvature_centre_fixed']['rms']), n=v['n'],
                         sep_curv=r2(v['sep_from_candidates']['Silk_curvature_centre']),
                         sep_knies2=r2(v['sep_from_candidates']['Monogem_2SNR_Knies24']))
                 for k, v in gj['circle_fits'].items()},
    patch_pgcc=gj['patch_pgcc'], beta_north=r2(gj['models']['T_sky_plane_bend_beta_north_deg']))
out['models'] = dict(S={k: dict(vexp_needed=[r2(x) for x in v['vexp_needed_signed']], G_south=r2(v['G_side_south_median']),
                                G_tangent=r2(v['G_side_tangent']), theta_R=r2(v['theta_R']))
                        for k, v in gj['models']['S_single_sphere'].items()},
                     T=[dict(v_s=t['v_s'], tilt=[r2(x) for x in t['tilt_deg']],
                             depth_pc=[r2(x) for x in t['depth_offset_pc']]) for t in gj['models']['T_tilted_normal']],
                     CH_centroid_range=[r2(min(x['centroid_one_sided'] for x in gj['models']['CH_one_sided_chord'])),
                                        r2(max(x['centroid_one_sided'] for x in gj['models']['CH_one_sided_chord']))])
out['shock'] = dict(L=kj['L_summary'], V=kj['V_summary'],
                    L_end_o2o3=[r2(min(r['o2o3_end'] for r in kj['L'])), r2(max(r['o2o3_end'] for r in kj['L']))],
                    complete_o2o3=[r2(min(c_['ratios']['O2_Hb'] / c_['ratios']['O3_Hb'] for c_ in kj['patch_complete_runs']
                                          if c_['ratios']['O3_Hb'] > 1)),
                                   r2(max(c_['ratios']['O2_Hb'] / c_['ratios']['O3_Hb'] for c_ in kj['patch_complete_runs']
                                          if c_['ratios']['O3_Hb'] > 1))],
                    patch_best=[(p['tag'], {k: r2(v) for k, v in p['ratios'].items()}) for p in kj['patch_best10'][:3]])
json.dump(out, open(f'{DER}/w3_northsouth_summary.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('desi_zones', 'desi_velocity', 'rm', 'models', 'shock')}, indent=0)[:6000])
print(json.dumps(out['geometry']['circle_fits']))
print(json.dumps(out['geometry']['crossings']))
print(json.dumps(out['geometry']['patch_pgcc'])[:900])
