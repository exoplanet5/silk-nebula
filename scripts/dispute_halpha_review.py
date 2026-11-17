"""dispute_halpha_review.py -- review checks for Section 4 (tag halpha) of the Comment on Wu et al. 2026
(AJ 172, 298; arXiv:2610.13329v1).  Adds numbers used in notes/dispute/halpha.md that are not in
data/derived/dispute_halpha_results.json:

  1. NIST ASD H I lines at 6564.4-6564.8 A (vacuum), fetched live, to confirm the 2s-3p wavelengths hard-coded in
     scripts/dispute_halpha_analysis.py and the geocoronal fine-structure offset (2:1 weighting).
  2. Rest wavelength that would by itself produce the observed Table 4 Halpha - [N II] offset.
  3. Moon altitude / illumination and target ecliptic latitude at the three exposure mid-times (what kind of
     solar-spectrum continuum could carry solar Halpha absorption).
  4. Zhang+2021 / Wu+2022 sky-ratio relation evaluated per exposure vs the exposure mean.
  5. Paper [S II]/Halpha and [N II]/Halpha (Table 2) after the adopted joint-fit (Fsc) flux factor -- model only.
  6. Summary of the self-consistent fits (Bsc, Fsc, Gsc) next to the earlier non-self-consistent ones (B, F, G).

Reads (not modified): data/derived/dispute_halpha_results.json, data/derived/verify_velocity_lamost_arms.json.
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_halpha_review.py
Output: data/derived/dispute_halpha_review.json
"""
import json
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = 'warn'
from astropy.coordinates import (SkyCoord, EarthLocation, AltAz, get_body, get_sun,   # noqa: E402
                                 GeocentricTrueEcliptic)
from astropy.time import Time                     # noqa: E402
import astropy.units as u                         # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/dispute_halpha_review.json'
C = 299792.458
NIST_URL = ('https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=H+I&limits_type=0&low_w=6564.4&upp_w=6564.8'
            '&unit=0&de=0&format=2&line_out=0&remove_js=on&en_unit=0&output=0&page_size=15&show_obs_wl=1'
            '&show_calc_wl=1&order_out=0&show_av=3&tsb_value=0&A_out=0&intens_out=on&allowed_out=1&forbid_out=1'
            '&conf_out=on&term_out=on&enrg_out=on&J_out=on')


def nist_lines():
    req = urllib.request.Request(NIST_URL, headers={'User-Agent': 'generic-research-client/1.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        txt = r.read().decode('utf-8', 'replace')
    rows = []
    for line in txt.splitlines()[1:]:
        f = [re.sub(r'[="]', '', x) for x in line.split(',')]
        if len(f) < 13:
            continue
        rows.append(dict(obs_vac=f[0] or None, ritz_vac=f[1] or None, Aki=f[3] or None,
                         lower=(f[7] + ' ' + f[9]).strip(), upper=(f[10] + ' ' + f[12]).strip(), type=f[13] or None))
    return rows


def main():
    R = json.load(open(ROOT / 'data/derived/dispute_halpha_results.json'))
    out = {'_retrieved': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}

    # 1. NIST
    try:
        rows = nist_lines()
        out['nist_HI_6564'] = dict(url=NIST_URL, lines=rows)
        l32 = float([r for r in rows if r['lower'] == '2s 1/2' and r['upper'] == '3p 3/2'][0]['ritz_vac'])
        l12 = float([r for r in rows if r['lower'] == '2s 1/2' and r['upper'] == '3p 1/2'][0]['ritz_vac'])
        lg = (2 * l32 + l12) / 3
        out['geocoronal_2to1'] = dict(l_3p32=l32, l_3p12=l12, l_centroid=round(lg, 6),
                                      dv_vs_6564p61=round(C * (lg / 6564.61 - 1), 3),
                                      dv_vs_6564p632=round(C * (lg / 6564.632 - 1), 3),
                                      matches_analysis_script=bool(abs(l32 - 6564.537685) < 1e-6
                                                                   and abs(l12 - 6564.584403) < 1e-6))
    except Exception as e:                       # noqa: BLE001
        out['nist_error'] = str(e)

    # 2. rest wavelength that would give the observed offset by itself
    m, e = R['offsets_table4']['Ha_minus_NII']['wmean'], R['offsets_table4']['Ha_minus_NII']['err']
    lam = 6564.61 * (1 - m / C)
    out['rest_wavelength_needed'] = dict(offset_kms=m, err_kms=e, lambda_rest=round(lam, 3),
                                         err_A=round(6564.61 * e / C, 3),
                                         below_6564p61_A=round(6564.61 - lam, 3), below_6564p632_A=round(6564.632 - lam, 3),
                                         per_region_range_A=[round(6564.61 * x / C, 3) for x in
                                                             (min(v['Ha_minus_NII'] for v in R['offsets_table4']['per_region'].values()),
                                                              max(v['Ha_minus_NII'] for v in R['offsets_table4']['per_region'].values()))])

    # 3. Moon and ecliptic latitude
    lamj = json.load(open(ROOT / 'data/derived/verify_velocity_lamost_arms.json'))['helio']
    times = sorted({x['mid_utc'] for x in lamj['exposures']})
    t = Time(times, scale='utc')
    loc = EarthLocation(lon=117.58 * u.deg, lat=40.39 * u.deg, height=960 * u.m)
    aa = AltAz(obstime=t, location=loc)
    moon = get_body('moon', t, loc)
    sun = get_sun(t)
    elong = sun.separation(moon).rad
    plate = SkyCoord(116.25248 * u.deg, 15.89500 * u.deg)
    out['moon_and_zodi'] = dict(exposure_mid_utc=times,
                                moon_alt_deg=[round(float(x), 1) for x in moon.transform_to(aa).alt.deg],
                                moon_illuminated_frac=[round(float(x), 2) for x in (1 - np.cos(elong)) / 2],
                                plate_alt_deg=[round(float(x), 1) for x in plate.transform_to(aa).alt.deg],
                                plate_ecliptic_lat_deg=round(float(plate.transform_to(GeocentricTrueEcliptic()).lat.deg), 1))

    # 4. per-exposure sky ratio vs mean
    g = R['geocoronal']
    rr = np.array(g['model_ratio_per_exposure'])
    out['sky_ratio_epochs'] = dict(per_exposure=list(rr), mean=g['model_ratio_mean'],
                                   first_over_mean_minus1_pct=round(100 * (rr[0] / g['model_ratio_mean'] - 1), 1),
                                   last_over_mean_minus1_pct=round(100 * (rr[-1] / g['model_ratio_mean'] - 1), 1))

    # 5. [S II]/Halpha after the adopted joint-fit flux factor (model Fsc; Table 2 ratios)
    sii_ha = {'P1': 0.93, 'P2': 0.88, 'P3': 0.83, 'P4': 0.78}          # Table 2 (p. 7)
    nii_ha = {'P1': 1.11, 'P2': 0.88, 'P3': 1.02, 'P4': 0.95}
    fF = {r: R['side_effects_joint_fit_Fsc'][r]['flux_fit_over_true'] for r in sii_ha}
    out['ratios_after_flux_factor_Fsc'] = dict(flux_fit_over_true=fF,
                                             SII_Ha=({r: round(sii_ha[r] * fF[r], 2) for r in sii_ha}),
                                             NII_Ha=({r: round(nii_ha[r] * fF[r], 2) for r in nii_ha}),
                                             note='ratio_true = ratio_paper x (fitted Ha flux / true Ha flux); model only')

    # 6. fit summary
    keys = ('B_refined_inputs', 'F_joint_centroid_width_WHAMfwhm', 'G_F_dm0_only',
            'Bsc_centroid_only_selfcons', 'Fsc_joint_selfcons', 'Gsc_joint_dm0_selfcons')
    out['fit_summary'] = {k: {q: R['fits'][k].get(q) for q in
                              ('oversub_pct', 'oversub_pct_1sigma', 'dm', 'chi2', 'dof', 'pred_shift', 'pred_fwhm_int',
                               'unexplained_offset', 'explained_fraction_of_wmean', 'true_ratio_implied',
                               'g_selfcons', 'flux_factor', 'true_Ha_rel_intensity', 'converged')} for k in keys}
    sF = R['side_effects_joint_fit_Fsc']
    out['Fsc_side_effects_ranges'] = dict(
        dfwhm=[round(min(sF[r]['fwhm_int_fit'] - sF[r]['fwhm_int_true'] for r in sF), 1),
               round(max(sF[r]['fwhm_int_fit'] - sF[r]['fwhm_int_true'] for r in sF), 1)],
        dflux_pct=[round(100 * (min(sF[r]['flux_fit_over_true'] for r in sF) - 1), 0),
                   round(100 * (max(sF[r]['flux_fit_over_true'] for r in sF) - 1), 0)],
        dT_K=[min(sF[r]['dT_eq2_K'] for r in sF), max(sF[r]['dT_eq2_K'] for r in sF)])
    OUT.write_text(json.dumps(out, indent=1))
    for k, v in out.items():
        if k != 'nist_HI_6564':
            print(k, json.dumps(v))


if __name__ == '__main__':
    main()
