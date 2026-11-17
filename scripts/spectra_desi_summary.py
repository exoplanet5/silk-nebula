"""spectra_desi_summary.py -- collect the headline DESI numbers (derived) into data/derived/spectra_desi_summary.json.

Derived quantities:
  * [O III]/Hb and [O III]/Ha lower limits (1 / 97.5th percentile of the bootstrap Hb/[O III], Ha/[O III]);
  * R23 = ([O II]3727 + [O III]4959+5007)/Hb lower limit (energy) from the bootstrap of Hb/([O II]+[O III]);
  * n_e from [O II] 3729/3726 (Sanders et al. 2016, 10^4 K: R = a(b+n)/(c+n), a=0.3771, b=2468, c=638.4);
  * n_e from [S II] 6716/6731 for the Halpha-structure component (Sanders+16: a=0.4315, b=2107, c=627.1);
  * ridge-peak [O III] brightness by Dec (forward model) and strict lower limits from positive binned residuals.
Usage: python -I scripts/spectra_desi_summary.py
"""
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
import sys
sys.path.insert(0, str(ROOT / 'scripts'))
from spectra_desi_lines import LINES  # noqa: E402


def ne_sanders(R, a, b, c):
    R = np.asarray(R, float)
    with np.errstate(invalid='ignore', divide='ignore'):
        n = (c * R - a * b) / (a - R)
    return np.where(n > 0, n, 0.0)


def main():
    J = json.load(open(D / 'spectra_desi_ratios.json'))
    res = json.load(open(D / 'spectra_desi_results.json'))
    B = np.load(D / 'spectra_desi_ratios_boot.npz', allow_pickle=False)
    L = LINES
    out = {'n_spectra_analysed': res['n_spectra'], 'n_fibres_OIII_ok': res['OIII_chi_stats']['n']}
    hb = B['Hb_at_vO'][:, 0]
    o3 = B['OIII5007_at_vO'][:, 0]
    ha = B['HaNII_at_vO'][:, 0]
    o2 = B['OII3727_at_vO'][:, 0] + B['OII3727_at_vO'][:, 1]
    # energy conversions (amplitudes are photon units)
    e_hb_o3 = hb / o3 * (L['OIII5007'] / L['Hb'])
    e_ha_o3 = ha / o3 * (L['OIII5007'] / L['Ha'])
    num = o2 * (L['Hb'] / L['OII3726']) + o3 * (1 + 1 / 2.98) * (L['Hb'] / L['OIII5007'])   # in Hb-energy units
    q = hb / num
    out['OIII5007_over_Hb'] = dict(lower_97p5=float(1 / np.percentile(e_hb_o3, 97.5)),
                                   lower_84=float(1 / np.percentile(e_hb_o3, 84)),
                                   median_of_inverse=float(1 / np.median(e_hb_o3)) if np.median(e_hb_o3) > 0 else None)
    out['OIII5007_over_Ha'] = dict(lower_97p5=float(1 / np.percentile(e_ha_o3, 97.5)),
                                   lower_84=float(1 / np.percentile(e_ha_o3, 84)),
                                   median_of_inverse=float(1 / np.median(e_ha_o3)) if np.median(e_ha_o3) > 0 else None)
    out['R23'] = dict(lower_97p5=float(1 / np.percentile(q, 97.5)), lower_84=float(1 / np.percentile(q, 84)),
                      median=float(1 / np.median(q)) if np.median(q) > 0 else None,
                      note='([O II]3726+3729 + [O III]4959+5007)/Hb, energy; Hb at the [O III] velocity')
    r2 = J['model_all']['ratios']['OII3729/3726']
    out['ne_OII'] = dict(ratio=r2['energy'], ratio_p16=r2['p16'], ratio_p84=r2['p84'], ratio_p2p5=r2['p2p5'],
                         ratio_p97p5=r2['p97p5'],
                         ne=float(ne_sanders(r2['energy'], 0.3771, 2468, 638.4)),
                         ne_from_p84_ratio=float(ne_sanders(r2['p84'], 0.3771, 2468, 638.4)),
                         ne_from_p16_ratio=float(ne_sanders(r2['p16'], 0.3771, 2468, 638.4)),
                         ne_upper_from_p2p5=float(ne_sanders(r2['p2p5'], 0.3771, 2468, 638.4)),
                         note='Sanders+2016 at 1e4 K; low-density limit 1.5 (ratio 3729/3726)')
    if 'nbrHa_all' in J:
        rs = J['nbrHa_all']['ratios']['SII6716/6731']
        out['ne_SII_Ha_structure'] = dict(ratio=rs['energy'], p16=rs['p16'], p84=rs['p84'],
                                          ne=float(ne_sanders(rs['energy'], 0.4315, 2107, 627.1)))
    # ridge brightness from the forward model
    fm = res['forward_model']
    nodes = []
    for nm, a, e in zip(fm['basis'], fm['amp_R'], fm['err_R']):
        st, d = nm.split('@')
        nodes.append(dict(strand=st, dec=float(d), amp_R=a, err_R=e))
    good = [n for n in nodes if n['err_R'] < 0.8]
    south = [n['amp_R'] for n in good if n['dec'] < 19.5]
    north = [n['amp_R'] for n in good if n['dec'] >= 21.5]
    out['ridge_peak_R'] = dict(nodes=nodes, south_median=float(np.median(south)) if south else None,
                               south_range=[float(min(south)), float(max(south))] if south else None,
                               north_median=float(np.median(north)) if north else None,
                               north_range=[float(min(north)), float(max(north))] if north else None)
    # strict lower limits: positive binned residuals (fibre minus petal mean <= true SB if SB >= 0 everywhere)
    rows = list(csv.DictReader(open(D / 'spectra_desi_fiber_table.csv')))
    ra = np.array([float(r['ra']) for r in rows])
    dec = np.array([float(r['dec']) for r in rows])
    v = np.array([float(r['OIII_R']) for r in rows])
    e = np.array([float(r['OIII_R_err']) for r in rows])
    ok = np.array([r['ok'] == '1' for r in rows]) & np.isfinite(v) & (e > 0)
    bs = 0.15
    xe = np.arange(113.0, 120.01, bs / np.cos(np.radians(20)))
    ye = np.arange(11.5, 28.81, bs)
    w = 1 / e[ok] ** 2
    sw, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye], weights=w)
    swv, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye], weights=w * v[ok])
    n, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye])
    with np.errstate(invalid='ignore', divide='ignore'):
        m = swv / sw
        em = 1 / np.sqrt(sw)
    lims = []
    for lo, hi in ((12, 15.5), (15.5, 17.5), (17.5, 19.5), (19.5, 22), (22, 25), (25, 28.8)):
        sel = (0.5 * (ye[:-1] + ye[1:]) >= lo) & (0.5 * (ye[:-1] + ye[1:]) < hi)
        mm, ee, nn = m[:, sel], em[:, sel], n[:, sel]
        good_ = (nn >= 4) & np.isfinite(mm) & (mm - 3 * ee > 0)
        if good_.any():
            k = np.argmax(np.where(good_, mm - 2 * ee, -np.inf))
            i, j = np.unravel_index(k, mm.shape)
            lims.append(dict(dec_lo=lo, dec_hi=hi, max_bin_R=float(mm[i, j]), err_R=float(ee[i, j]),
                             lower_limit_2sigma_R=float(mm[i, j] - 2 * ee[i, j]), n_fibres=int(nn[i, j]),
                             ra=float(0.5 * (xe[i] + xe[i + 1])), dec=float(ye[:-1][sel][j] + 0.5 * bs)))
        else:
            lims.append(dict(dec_lo=lo, dec_hi=hi, max_bin_R=None))
    out['strict_lower_limits_0p15deg'] = lims
    out['detection'] = res.get('detection')
    out['OIII_v_vs_vcorr'] = res.get('OIII_v_vs_vcorr')
    out['velocities'] = dict(all=J['model_all']['v_ref_free'], all_boot=J['model_all']['v_ref_boot'],
                             nbr=J['nbr_all']['v_ref_free'], v_by_dec=J.get('v_by_dec'),
                             v_by_subset=J.get('v_by_subset'),
                             Ha_structure=J.get('nbrHa_all', {}).get('v_ref_free'),
                             Ha_structure_boot=J.get('nbrHa_all', {}).get('v_ref_boot'),
                             posP=J.get('model_posP', {}).get('v_ref_free'), negP=J.get('model_negP', {}).get('v_ref_free'))
    out['OIII_width'] = J['model_all']['width']
    out['Te_OIII'] = J['model_all'].get('Te_OIII')
    (D / 'spectra_desi_summary.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: v for k, v in out.items() if k not in ('ridge_peak_R',)}, indent=1, default=float)[:6000])
    print('ridge', json.dumps({k: v for k, v in out['ridge_peak_R'].items() if k != 'nodes'}))


if __name__ == '__main__':
    main()
