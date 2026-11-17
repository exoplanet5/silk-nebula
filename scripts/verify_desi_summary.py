"""verify_desi_summary.py -- collect the key wave-2 DESI verification numbers into one JSON and print them.
Output: data/derived/verify_desi_summary.json
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D  # noqa: E402


def main():
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    br = json.loads((D / 'verify_desi_brightness.json').read_text())
    ra = json.loads((D / 'verify_desi_ratios.json').read_text())
    ve = json.loads((D / 'verify_desi_velocity.json').read_text())
    rf = json.loads((D / 'verify_desi_regspec_fits.json').read_text())
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    rot = np.load(D / 'verify_desi_detect_nulls.npz')['rot']
    S = {'n_spectra': int(len(lf['sparcl_id'])), 'n_wave1': int(np.sum(lf['src'] == 'w1')),
         'n_wave2_fetch': int(np.sum(lf['src'] == 'ext'))}
    m05 = np.abs(rot[:, 2]) < 0.05
    m02 = np.abs(rot[:, 2]) < 0.02
    S['detection'] = dict(snr=det['true']['snr'], snr_clip8=det['true']['snr_clip8'], n=det['n'],
                          rot_null_std_r005=float(rot[m05, 0].std()), rot_null_max_r005=float(rot[m05, 0].max()),
                          rot_null_std_r002=float(rot[m02, 0].std()), n_r005=int(m05.sum()), n_r002=int(m02.sum()),
                          z_vs_rot_r005=float(det['true']['snr'] / rot[m05, 0].std()),
                          z_vs_rot_r02=det['rigid_rotations']['z_true_vs_clean'],
                          shuffle=det['shuffle_within_petal'], fake_line_snr=det['fake_line_snr'],
                          splits={k: round(v['snr'], 1) for k, v in det['splits'].items()},
                          jackknife_min=min(v['snr'] for v in det['jackknife_nside4'].values()))
    fm = br['forward_model']
    S['brightness'] = dict(conversion_ratio=br['conversion_ratio_gfa_over_wave1'],
                           south_median={k: v['south_median'] for k, v in fm.items()},
                           north_median={k: v['north_median'] for k, v in fm.items()},
                           nodes_gfa={d['node']: [round(d['amp_R'], 2), round(d['err_R'], 2)]
                                      for d in fm['gfa_chromatic']['nodes'] if d['err_R'] < 1},
                           nodes_wave1conv={d['node']: [round(d['amp_R'], 2), round(d['err_R'], 2)]
                                            for d in fm['wave1_psf_to_fiber']['nodes'] if d['err_R'] < 1},
                           wave1_bin=br['wave1_bin_116.27_17.28'], lower_limits_gfa=br['lower_limits']['gfa'],
                           look_elsewhere=br['look_elsewhere_fake_lines'], control=br.get('control_gt5deg'))
    S['ratios'] = {sn: ra['derived'][sn] for sn in ra['derived']}
    keys = ['Hb', 'OIII4959', 'OIII4363_jointHgFe', 'HgI4358_jointFe', 'FeII4359_jointHg', 'HeII4686', 'OII3726',
            'OII3729', 'NeIII3869', 'Ha', 'NII6583', 'SII6716', 'SII6731', 'SIII9069', 'ArIII7136',
            'OII7325_jointCa', 'CaII7291', 'CaII7324_joint', 'FeII7155', 'NI5200', 'NII5755', 'HeI5876', 'sky5577']
    S['line_over_5007'] = {sn: {k: [round(ra[sn][k]['ratio'], 4), round(ra[sn][k]['p16'], 4),
                                    round(ra[sn][k]['p84'], 4), round(ra[sn][k]['p2p5'], 4),
                                    round(ra[sn][k]['p97p5'], 4)] for k in keys if k in ra[sn]}
                           for sn in ('all', 'petalFE', 'south_dec<19.5', 'north_dec>=19.5', 'QSO_or_z>0.5',
                                      'smooth_continuum', 'nbr_instrument', 'dark', 'no_Ha_patch',
                                      'wave1_fibres_only', 'blocks_nside64')}
    S['fake_lines'] = dict(near4363=ra['fake_4363_summary'], all_chi_rms=ra['fake_all_chi']['chi_rms'])
    S['regspec_4363'] = {sn: dict(centre=rf[sn]['scan_feature_near_4363']['best_centre'],
                                  lo=rf[sn]['scan_feature_near_4363']['lo'], hi=rf[sn]['scan_feature_near_4363']['hi'],
                                  amp_4363_Hg=rf[sn]['fixed_4363_Hg']['amp'], err=rf[sn]['fixed_4363_Hg']['err'],
                                  amp_4363_Hg_Fe=rf[sn]['fixed_4363_Hg_Fe']['amp'],
                                  err_Fe=rf[sn]['fixed_4363_Hg_Fe']['err'], expected=rf[sn]['expected'])
                         for sn in rf}
    S['velocity'] = ve
    (D / 'verify_desi_summary.json').write_text(json.dumps(S, indent=1))
    print(json.dumps({k: S[k] for k in ('n_spectra', 'n_wave1', 'n_wave2_fetch', 'detection')}, indent=1))
    print(json.dumps(S['brightness'], indent=1)[:6000])
    print(json.dumps(S['ratios']['all'], indent=1))
    print(json.dumps(S['ratios']['south_dec<19.5'], indent=1)[:3000])
    for sn, d in S['line_over_5007'].items():
        print(sn, {k: v[0] for k, v in d.items()})
    print(json.dumps(S['regspec_4363'], indent=1)[:3000])


if __name__ == '__main__':
    main()
