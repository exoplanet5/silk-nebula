"""dispute_oiii: numbers for Section 3 of the Comment on Wu et al. 2026 (AJ 172, 298; arXiv:2610.13329v1),
"The [O III] velocity offset is a wavelength-reference artefact".

What this script does (no new spectral fitting except a check refit of the [N I] doublet):
  1. Recomputes the barycentric correction v_corr with astropy for the three MRS exposures of 2024-02-01
     (mid-exposure UTC from data/derived/verify_velocity_lamost_arms.json) at the field centre and at P1-P3
     (Table 1 of the paper), Xinglong (117.58 E, 40.39 N, 960 m).
  2. Recomputes V_LSR - V_helio for the paper's own solar-motion convention (20 km/s toward
     18h03m50.29s +30d00m16.8s J2000; Sect. III.2.3) at P1-P4.
  3. Converts the reference wavelengths ([N I] 5198/5200 airglow, [O III] 5007) to the frames that a plotted
     wavelength scale could be in (air/vacuum x topocentric/heliocentric) and gives the measured offset of the
     paper's plotted airglow from each (Fig. 3a polyline fits from verify_velocity_fig3_measure/solar.json).
  4. Builds the per-region table (published Table 4 [O III], Fig. 3b plotted-scale re-measure, corrected values,
     DESI DR1 within 0.5 deg, red lines, WHAM nearest beam) and the velocity budget.

Inputs (read only): data/derived/verify_velocity_{fig3_measure,fig3_solar,summary,lamost_arms,forensics,
desi_sdss}.json, verify_desi_velocity.json, verify_desi_audit_vel.json, verify_sdss_results_base_r0.15.json,
spectra_ism_wham_results.json, data/derived/verify_velocity_fig3_spectra.npz.
Output: data/derived/dispute_oiii_values.json

Run: ~/.venvs/astro313/bin/python -I scripts/dispute_oiii_values.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
OUT = D / 'dispute_oiii_values.json'
C = 299792.458


def J(name):
    return json.loads((D / name).read_text())


def n_air(lv):  # Morton (2000) / IAU standard, same as the verify_velocity scripts
    s2 = (1e4 / np.asarray(lv, float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


def dv(lam, lam0):
    return C * (lam / lam0 - 1.0)


# ------------------------------------------------------------------ paper values (transcribed)
# Table 1 (p. 6): region centres (deg, J2000)
REG = {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72), 'P4': (115.20, 16.18)}
# Table 4 (p. 13): "RVs" column (stated in Sect. III.2.3, p. 10, to be V_helio) and V_LSR column
T4 = {'P1': {'OIII': (26.8, 4.3), 'SII': (11.8, 2.5), 'NII': (11.8, 1.6), 'Ha': (16.3, 1.2)},
      'P2': {'OIII': (26.8, 4.1), 'SII': (12.0, 1.6), 'NII': (12.0, 1.6), 'Ha': (17.0, 1.2)},
      'P3': {'OIII': (23.8, 3.3), 'SII': (11.7, 1.9), 'NII': (11.7, 1.9), 'Ha': (18.4, 1.4)},
      'P4': {'SII': (12.7, 2.1), 'NII': (13.3, 2.3), 'Ha': (19.7, 1.4)}}
T4_LSR = {'P1': {'OIII': 15.5, 'SII': -1.0, 'NII': -0.9, 'Ha': 4.7},
          'P2': {'OIII': 14.5, 'SII': 0.5, 'NII': -0.5, 'Ha': 4.0},
          'P3': {'OIII': 11.2, 'SII': 1.7, 'NII': 0.4, 'Ha': 6.1},
          'P4': {'SII': 1.1, 'NII': 1.0, 'Ha': 4.5}}


def main():
    res = {}
    # cross-check the transcription against the wave-1 transcription
    w = J('spectra_ism_wham_results.json')
    mism = [(r, l) for r in T4 for l in T4[r] if abs(w['paper_vhelio_table4'][r][l] - T4[r][l][0]) > 1e-6]
    res['table4_transcription_mismatches_vs_spectra_ism'] = mism

    # ------------------------------------------------------------------ 1. v_corr with astropy
    arms = J('verify_velocity_lamost_arms.json')
    mids = sorted({e['mid_utc'] for e in arms['helio']['exposures']})
    loc = EarthLocation(lon=117.58 * u.deg, lat=40.39 * u.deg, height=960 * u.m)
    pos = dict(field_centre=(116.2525, 15.895), **{k: v for k, v in REG.items() if k != 'P4'})
    vc = {}
    for name, (ra, de) in pos.items():
        c = SkyCoord(ra * u.deg, de * u.deg)
        vals = [float(c.radial_velocity_correction('barycentric', obstime=Time(t, scale='utc'),
                                                   location=loc).to_value(u.km / u.s)) for t in mids]
        vc[name] = dict(per_exposure=vals, mean=float(np.mean(vals)))
    allv = np.array([v for d in vc.values() for v in d['per_exposure']])
    res['v_corr'] = dict(mid_utc=mids, by_position=vc, range_all=[float(allv.min()), float(allv.max())],
                         adopted=-8.32, header_HELIO_RV_median=arms['helio']['header_HELIO_RV_median'],
                         header_VACUUM=arms['helio']['header_VACUUM'], header_HELIO=arms['helio']['header_HELIO'],
                         note='V_helio = V_topo + v_corr; LAMOST HELIO_RV = -v_corr')
    VCORR = -8.32

    # ------------------------------------------------------------------ 2. LSR correction (paper convention)
    apex = SkyCoord('18h03m50.29s', '+30d00m16.8s', frame='icrs')
    lsr = {}
    for k, (ra, de) in REG.items():
        th = SkyCoord(ra * u.deg, de * u.deg).separation(apex).to_value(u.rad)
        lsr[k] = float(20.0 * np.cos(th))
    res['vlsr_minus_vhelio'] = lsr
    # internal consistency of the paper's own V_LSR column
    res['paper_vlsr_minus_vhelio_table4'] = {r: {l: round(T4_LSR[r][l] - T4[r][l][0], 1) for l in T4[r]} for r in T4}

    # ------------------------------------------------------------------ 3. frames for the airglow and [O III]
    ni_vac = np.array([1e8 / 19233.177, 1e8 / 19224.464])            # NIST level energies
    ni_air = ni_vac / n_air(ni_vac)
    ni_air_han = np.array([5197.9282, 5200.2856])                      # Hanuschik 2003, table 6 (air)
    o3_air = 5006.843
    o3_vac = float(o3_air * n_air(o3_air * 1.000279))                  # iterate once (good to 1e-4 A)
    for _ in range(4):
        o3_vac = float(o3_air * n_air(o3_vac))
    f_h = 1 + VCORR / C                                                  # topocentric-rest line on a heliocentric grid
    frames = {
        'air_topocentric': ni_air.tolist(),
        'air_heliocentric': (ni_air * f_h).tolist(),
        'vacuum_topocentric': ni_vac.tolist(),
        'vacuum_heliocentric_LAMOST_native': (ni_vac * f_h).tolist(),
    }
    fm = J('verify_velocity_fig3_measure.json')
    fs = J('verify_velocity_fig3_solar.json')
    sk = fs['sky_lines_on_plot_scale']
    mu = np.array([sk['NI5198 (NIST)']['mu'], sk['NI5200 (NIST)']['mu']])
    emu = np.array([sk['NI5198 (NIST)']['e_mu'], sk['NI5200 (NIST)']['e_mu']])
    offs = {}
    for fr, lam in frames.items():
        v = dv(mu, np.array(lam))
        e = C * emu / mu
        wm = float(np.sum(v / e ** 2) / np.sum(1 / e ** 2))
        offs[fr] = dict(per_line=v.round(2).tolist(), wmean=round(wm, 2), err=round(float(1 / np.sqrt(np.sum(1 / e ** 2))), 2))
    v_han = dv(mu, ni_air_han)
    res['airglow_frames'] = dict(ni_vac=ni_vac.tolist(), ni_air_nist=ni_air.tolist(), ni_air_hanuschik=ni_air_han.tolist(),
                                 frame_positions_A=frames, plotted_centres_A=mu.tolist(), plotted_centres_err_A=emu.tolist(),
                                 offset_from_frame_kms=offs, offset_from_hanuschik_kms=v_han.round(2).tolist())
    summ = J('verify_velocity_summary.json')
    zp = summ['sky_zero_point_plot_scale']
    # uncertainty of the sky zero point: fit error of the two [N I] lines (+) half the NIST-Hanuschik difference
    e_fit = offs['air_topocentric']['err']
    e_rest = 0.5 * abs(zp['mean_NIST'] - zp['mean_Hanuschik'])
    e_vc = 0.5 * (allv.max() - allv.min())
    res['sky_zero_point'] = dict(adopted=zp['adopted'], mean_NIST=zp['mean_NIST'], mean_Hanuschik=zp['mean_Hanuschik'],
                                 weak_lines=fs['weak_sky_lines_zero_point'], e_fit=e_fit, e_rest=e_rest,
                                 e_quad=float(np.hypot(e_fit, e_rest)))
    # Frame test: where (relative to the NIST air rest wavelength, in km/s) the airglow must sit on the plotted
    # scale under each hypothesis, against the adopted plotted-scale sky zero point (five sky lines).
    zp_v, zp_e = zp['adopted'], float(np.hypot(e_fit, e_rest))
    vac_minus_air = float(np.mean(dv(ni_vac, ni_air)))
    expect = {'air_topocentric': 0.0, 'air_heliocentric': VCORR, 'air_heliocentric_signflip': -VCORR,
              'vacuum_topocentric': vac_minus_air, 'vacuum_heliocentric_LAMOST_native': vac_minus_air + VCORR}
    res['frame_test'] = {fr: dict(expected=round(ev, 2), measured=round(zp_v, 2), err=round(zp_e, 2),
                                  n_sigma=round(abs(zp_v - ev) / zp_e, 1)) for fr, ev in expect.items()}

    # check refit of the [N I] doublet on the extracted Fig. 3a samples (common sigma, linear continuum)
    sp = np.load(D / 'verify_velocity_fig3_spectra.npz')
    x, y = sp['a_wave'], sp['a_flux']
    m = (x > 5196.4) & (x < 5201.6)

    def dbl(xx, c0, c1, a1, a2, m1, m2, s):
        return (c0 + c1 * (xx - 5199) + a1 * np.exp(-0.5 * ((xx - m1) / s) ** 2) + a2 * np.exp(-0.5 * ((xx - m2) / s) ** 2))
    p, cov = curve_fit(dbl, x[m], y[m], p0=[58, 0, 37, 22, 5197.94, 5200.27, 0.27])
    pe = np.sqrt(np.diag(cov))
    res['NI_doublet_refit'] = dict(c0=p[0], c1=p[1], A1=p[2], A2=p[3], mu1=p[4], mu2=p[5], sig=p[6],
                                   e_mu1=pe[4], e_mu2=pe[5], ratio=p[2] / p[3],
                                   v_mu1_vs_air_nist=float(dv(p[4], ni_air[0])), v_mu2_vs_air_nist=float(dv(p[5], ni_air[1])),
                                   npts=int(m.sum()))

    # ------------------------------------------------------------------ [O III] on the plotted scale
    o = fm['panel_a_fits']['OIII5007']
    sol = fs['OIII_refits']['5007']['with_solar'][1]
    V_plot = float(dv(o['mu'], o3_air))
    budget = dict(V_plot=summ['fig3a']['V_plot'], sky_zero_point=-zp['adopted'], solar_blend=summ['fig3a']['solar_corr'],
                  v_corr=VCORR)
    budget['V_topo'] = budget['V_plot'] + budget['sky_zero_point'] + budget['solar_blend']
    budget['V_helio'] = budget['V_topo'] + VCORR
    budget['err'] = summ['fig3a']['err']
    res['oiii_rest'] = dict(air=o3_air, vac=o3_vac, vac_minus_air_kms=float(dv(o3_vac, o3_air)))
    res['oiii_fig3a'] = dict(mu=o['mu'], e_mu=o['e_mu'], sig=o['sig'], A=o['A'], cont=o['cont'], V_plot=V_plot,
                             e_V_plot=float(C * o['e_mu'] / o['mu']), fig3c_fit_centre=fm['markers_A']['c_right_fitcentre'][0],
                             V_fig3c=float(dv(fm['markers_A']['c_right_fitcentre'][0], o3_air)),
                             V_plot_if_vacuum=fm['velocities']['OIII5007_plot_scale_vac'],
                             solar_refit_mu=sol['mu'], solar_shift=sol['shift_kms'],
                             oiii4959_V_plot=fm['velocities']['OIII4959_plot_scale_air'])
    res['budget_fig3a'] = budget
    OFFSET = -VCORR + zp['adopted']           # 10.09
    E_OFF = 1.0                                # verify_velocity.md sect. 1 (rounded up from e_quad, see sky_zero_point)
    res['correction'] = dict(offset=OFFSET, err=E_OFF, components=dict(missing_v_corr=-VCORR, sky_zero_point=zp['adopted']))
    # positions on the plotted scale of a line with a given true V_helio
    res['plotted_position_for_true_vhelio'] = {str(vh): dict(V_plot=vh + OFFSET, lam=o3_air * (1 + (vh + OFFSET) / C))
                                               for vh in (13.3, 25.0)}

    # ------------------------------------------------------------------ 4. per-region table
    dz = J('verify_desi_velocity.json')['regspec_centroids']
    ism = w['nearest_beams']
    rows = {}
    for r in ('P1', 'P2', 'P3', 'P4'):
        row = dict(lsr=lsr[r])
        if 'OIII' in T4[r]:
            v, e = T4[r]['OIII']
            b = fm['panel_b'][f'b_{r}']
            row['pub_vhelio'] = [v, e]
            row['pub_vlsr_table4'] = T4_LSR[r]['OIII']
            row['fig3b_plotted'] = [b['V_plot_air'], float(C * b['e_mu'] / b['mu'])]
            row['fig3b_minus_table4'] = b['V_plot_air'] - v
            cv, ce = v - OFFSET, float(np.hypot(e, E_OFF))
            row['corr_vhelio'] = [cv, ce]
            row['corr_vlsr'] = cv + lsr[r]
            ds = dz[f'within0.5deg_{r}']['per_fibre_pos_SN4']
            row['desi_0p5'] = [ds['wmean'], ds['err'], ds['n']]
            row['pub_minus_desi'] = [v - ds['wmean'], float(np.hypot(e, ds['err']))]
            row['pub_minus_desi_sigma'] = row['pub_minus_desi'][0] / row['pub_minus_desi'][1]
            row['corr_minus_desi'] = [cv - ds['wmean'], float(np.hypot(ce, ds['err']))]
            row['corr_minus_desi_sigma'] = row['corr_minus_desi'][0] / row['corr_minus_desi'][1]
            nv, ne = T4[r]['NII']
            row['pub_oiii_minus_nii'] = [v - nv, float(np.hypot(e, ne))]
            row['corr_oiii_minus_nii'] = [cv - nv, float(np.hypot(ce, ne))]
            row['pub_vlsr_oiii_minus_nii_table4'] = T4_LSR[r]['OIII'] - T4_LSR[r]['NII']
        row['nii'] = list(T4[r]['NII'])
        row['sii'] = list(T4[r]['SII'])
        row['ha'] = list(T4[r]['Ha'])
        comps = ism[r][0]['components']
        main = max(comps, key=lambda q: q['I_R'])
        row['wham_nearest'] = dict(sep_deg=ism[r][0]['sep_deg'], v_lsr=main['v'], ev=main['ev'], fwhm_int=main['fwhm_int'],
                                   v_helio=main['v'] - lsr[r])
        rows[r] = row
    res['regions'] = rows

    # ------------------------------------------------------------------ independent [O III] velocities
    vd = J('verify_velocity_desi_sdss.json')
    au = J('verify_desi_audit_vel.json')['total_all']['Z']
    sd = J('verify_sdss_results_base_r0.15.json')['regions']['south_dec_lt20']['velocity_OIII']
    comb = summ['combined']
    res['independent'] = dict(
        desi_stack=dict(v=vd['desi']['pred_pos_gt0.5R']['V_helio'], e_stat=vd['desi']['pred_pos_gt0.5R']['err_boot'],
                        e_tot=comb['inputs']['DESI'][1], n=vd['desi']['pred_pos_gt0.5R']['n'], src='verify_velocity_desi_sdss.json'),
        desi_rev2=dict(v=au['v'], p16=au['v_boot_p16_p84'][0], p84=au['v_boot_p16_p84'][1], src='verify_desi_audit_vel.json total_all'),
        sdss_all=dict(v=vd['sdss']['all']['V_helio_5007'], e=vd['sdss']['all']['err_boot'], n_groups=vd['sdss']['all']['ngroups'],
                      src='verify_velocity_desi_sdss.json'),
        sdss_south=dict(v=sd['v_corr'], e=sd['ev_boot'], e_zp=0.8, src='verify_sdss_results_base_r0.15.json south_dec_lt20'),
        adopted=dict(v=comb['V_helio'], e=comb['err'], chi2=comb['chi2'], dof=comb['dof'], v_lsr=comb['V_LSR_standard']),
        adopted_alt_estimators=dict(v=12.4, e=0.9, src='w3_canon.md sect. 5.2 (verify_desi rev-2 + verify_sdss)'),
    )
    nii = np.array([T4[r]['NII'][0] for r in T4])
    sii = np.array([T4[r]['SII'][0] for r in T4])
    wh = np.array([rows[r]['wham_nearest']['v_helio'] for r in rows])
    res['red_lines'] = dict(nii_range=[nii.min(), nii.max()], nii_mean=nii.mean(), sii_range=[sii.min(), sii.max()],
                            wham_vhelio_range=[wh.min(), wh.max()], wham_vhelio_mean=wh.mean(),
                            ha_range=[min(T4[r]['Ha'][0] for r in T4), max(T4[r]['Ha'][0] for r in T4)])
    e_nii_mean = float(np.sqrt(np.sum(np.array([T4[r]['NII'][1] for r in T4]) ** 2)) / len(T4))
    res['nii_mean_err'] = e_nii_mean
    res['adopted_minus_nii_mean'] = [comb['V_helio'] - nii.mean(), float(np.hypot(comb['err'], e_nii_mean))]
    res['fig3a_minus_nii_mean'] = [budget['V_helio'] - nii.mean(), float(np.hypot(budget['err'], e_nii_mean))]

    # ------------------------------------------------------------------ aggregates over P1-P3
    def wmean(pairs):
        v = np.array([p[0] for p in pairs]); e = np.array([p[1] for p in pairs])
        w = 1 / e ** 2
        m = float(np.sum(w * v) / np.sum(w)); em = float(1 / np.sqrt(np.sum(w)))
        return dict(wmean=m, err=em, sigma=m / em, chi2=float(np.sum(w * (v - m) ** 2)), dof=len(v) - 1)
    P = ('P1', 'P2', 'P3')
    res['aggregates_P1P3'] = dict(
        pub_minus_desi=wmean([rows[r]['pub_minus_desi'] for r in P]),
        corr_minus_desi=wmean([rows[r]['corr_minus_desi'] for r in P]),
        pub_oiii_minus_nii=wmean([rows[r]['pub_oiii_minus_nii'] for r in P]),
        corr_oiii_minus_nii=wmean([rows[r]['corr_oiii_minus_nii'] for r in P]),
        fig3b_minus_table4=[rows[r]['fig3b_minus_table4'] for r in P],
        fig3b_minus_table4_mean=float(np.mean([rows[r]['fig3b_minus_table4'] for r in P])),
        table4_vs_fig3b_sigma=[rows[r]['fig3b_minus_table4'] / float(np.hypot(rows[r]['pub_vhelio'][1], rows[r]['fig3b_plotted'][1]))
                               for r in P])
    # [N I] zero point relative to Hanuschik's air wavelengths (the reference the paper cites for the 5198 line)
    eh = C * emu / mu
    res['airglow_frames']['offset_from_hanuschik_wmean'] = [float(np.sum(v_han / eh ** 2) / np.sum(1 / eh ** 2)),
                                                            float(1 / np.sqrt(np.sum(1 / eh ** 2)))]

    # ------------------------------------------------------------------ alternatives (from project files)
    fo = J('verify_velocity_forensics.json')
    air_refs = ('Hanuschik air 5197.928', 'NIST air 5197.902')
    vac_refs = ('vac 5199.349 (correct)', 'Hanuschik air->vac 5199.376')
    av = [q['dV_mult'] for q in fo['scenarios'] if q['measured'].startswith('strong') and
          ((q['reference'] in air_refs and 'vac' in q['o3_rest']) or (q['reference'] in vac_refs and 'air' in q['o3_rest']))]
    res['alternatives'] = dict(
        arms_b0_minus_r0=arms['MN_plate_coadd_b0_minus_r0_all'],
        b0_minus_gaia=arms['MN_plate_vs_gaia']['rv_b0_minus_gaia'], r0_minus_gaia=arms['MN_plate_vs_gaia']['rv_r0_minus_gaia'],
        box_plates_summary=arms['box_plates_summary'],
        air_vac_mixup_abs_range=[float(np.min(np.abs(av))), float(np.max(np.abs(av)))], air_vac_mixup_values=sorted(av),
        solar_blend_shift=sol['shift_kms'], solar_blend_shift_pure_solar=fs['OIII_refits']['5007']['with_solar'][-1]['shift_kms'],
        f_sun=fs['solar_regression']['f_sun'])

    # ------------------------------------------------------------------ review additions (adversarial pass)
    # (i) Which frame is Table 4 itself on?  Weighted mean of (Fig. 3b plotted-scale - Table 4) over P1-P3, against
    #     the value each hypothesis predicts:
    #       Table 4 measured on the plotted (air, topocentric, sky-referenced) scale ............ 0
    #       Table 4 truly heliocentric (plotted - sky offset + v_corr) ........................ +OFFSET (+10.1)
    #       v_corr applied numerically with the wrong sign after measuring on the plotted scale  +zp - |v_corr| (-6.5)
    d3b = [[rows[r]['fig3b_minus_table4'], float(np.hypot(rows[r]['pub_vhelio'][1], rows[r]['fig3b_plotted'][1]))]
           for r in P]
    wm3b = wmean(d3b)
    hyp = {'same_scale_as_plot': 0.0, 'table4_heliocentric': OFFSET, 'vcorr_sign_flipped': zp['adopted'] + VCORR}
    res['table4_frame_test'] = dict(
        fig3b_minus_table4_wmean=wm3b['wmean'], err=wm3b['err'], chi2_same_scale=wm3b['chi2'] + wm3b['sigma'] ** 2,
        hypotheses={k: dict(expected=h, n_sigma=abs(wm3b['wmean'] - h) / float(np.hypot(wm3b['err'], E_OFF)))
                    for k, h in hyp.items()},
        note='errors: Table 4 quoted errors (+) Fig. 3b fit errors; hypothesis sigma adds the +-1.0 correction error')
    # (ii) LAMOST-free check: DESI stack (incl. zero point) + SDSS all-strand, the two non-LAMOST inputs of the
    #      adopted mean, against the calibration-free Fig. 3a value and the P1-P3 published / corrected means.
    ind = wmean([[res['independent']['desi_stack']['v'], res['independent']['desi_stack']['e_tot']],
                 [res['independent']['sdss_all']['v'], res['independent']['sdss_all']['e']]])
    pub_w = wmean([rows[r]['pub_vhelio'] for r in P])
    cor_w = wmean([rows[r]['corr_vhelio'] for r in P])
    res['lamost_free_check'] = dict(
        desi_plus_sdss=dict(v=ind['wmean'], e=ind['err'], chi2=ind['chi2']),
        fig3a_minus_desi_sdss=[budget['V_helio'] - ind['wmean'], float(np.hypot(budget['err'], ind['err']))],
        pub_P1P3_wmean=[pub_w['wmean'], pub_w['err']], corr_P1P3_wmean=[cor_w['wmean'], cor_w['err']],
        pub_minus_desi_sdss=[pub_w['wmean'] - ind['wmean'], float(np.hypot(pub_w['err'], ind['err']))],
        corr_minus_desi_sdss=[cor_w['wmean'] - ind['wmean'], float(np.hypot(cor_w['err'], ind['err']))],
        note='P1-P3 weighted means treat the common +-1.0 correction as independent (errors slightly understated)')
    OUT.write_text(json.dumps(res, indent=1, default=float))

    # ------------------------------------------------------------------ print summary
    print('transcription mismatches:', mism)
    print('v_corr:', {k: round(v['mean'], 3) for k, v in vc.items()}, 'range', res['v_corr']['range_all'])
    print('V_LSR - V_helio:', {k: round(v, 2) for k, v in lsr.items()})
    print('paper V_LSR - V_helio (Table 4):', res['paper_vlsr_minus_vhelio_table4'])
    print('airglow offsets (km/s):', {k: (v['wmean'], v['err']) for k, v in offs.items()}, 'Hanuschik', v_han.round(2))
    print('frame test:', res['frame_test'])
    print('sky zp:', res['sky_zero_point'])
    print('[N I] refit:', {k: round(float(v), 4) for k, v in res['NI_doublet_refit'].items()})
    print('[O III] fig3a:', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res['oiii_fig3a'].items()})
    print('budget:', {k: round(v, 2) for k, v in budget.items()})
    print('correction:', res['correction'], res['plotted_position_for_true_vhelio'])
    for r, row in rows.items():
        print(r, {k: (np.round(v, 2).tolist() if isinstance(v, (list, float)) else v) for k, v in row.items()})
    print('independent:', json.dumps(res['independent'], default=float)[:900])
    print('red lines:', res['red_lines'], 'adopted-[N II]', res['adopted_minus_nii_mean'], 'fig3a-[N II]', res['fig3a_minus_nii_mean'])
    print('alternatives:', json.dumps(res['alternatives'], default=float)[:900])


if __name__ == '__main__':
    main()
