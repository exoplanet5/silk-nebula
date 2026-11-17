"""dispute_halpha_analysis.py -- numbers for Section 4 of the Comment on Wu et al. 2026 (AJ 172, 298;
arXiv:2610.13329v1): "The LAMOST Halpha velocity offset is plausibly a geocoronal-subtraction artefact".

Re-uses (imported, not modified):
  scripts/verify_velocity_halpha_sim.py  -- forward model of WIM Halpha + geocoronal Halpha - MRS-N sky template,
                                            fitted with Gaussian + quadratic (simulate(), gauss(), constants)
  scripts/verify_paper_kd.py             -- Reid+14 URC in the `kd` form, IAU solar motion, traditional solver
Reads (not modified):
  data/derived/verify_hiwham_wham_beams.csv          -- WHAM-SS DR1 beam fits (v_g = V_LSR of main component)
  data/derived/verify_velocity_halpha_sim.json        -- original simulation output (solar-Halpha term, grid)
  data/derived/verify_velocity_lamost_arms.json       -- exposure mid-times, header HELIO_RV
  data/derived/geometry_lamost_halpha_vlsr_fig5_digitized.csv  -- per-fibre Halpha V_LSR digitized from Fig. 5

Paper inputs (Wu+2026, arXiv v1 page numbers): Table 1 (positions), Table 2 p. 7 (Halpha relative intensity and
FWHM), Table 4 p. 13 (V_helio, V_LSR, distances). Atomic data: NIST ASD v5.12 (Kramida et al. 2024), H I lines
6562-6566 A (vacuum), retrieved 2026-11-16: 2s-3p3/2 6564.537685 A, 2s-3p1/2 6564.584403 A, blended 2-3 Ritz
6564.632 A (statistically weighted) and computed blend 6564.608 A.  MRS-N geocoronal model: Zhang+2021 (RAA 21,
280) Sect. 2.4 and Wu+2022 (RAA 22, 075015) eq. I(Ha_sky)/I(6554) = 0.954 - 0.011 (-sunalt).

Review fix (2026-11-17): the sky/nebula ratio g = 0.319 / (Table 2 Halpha) uses the published *fitted* Halpha flux,
which the model itself lowers by a factor f.  The adopted fits (Bsc, Fsc, Gsc) iterate g = 0.319 f / (Table 2 Halpha)
to self-consistency; Fig. 4 and the Comment use Fsc.  The earlier fits A-G are kept for comparison.

Run: ~/.venvs/astro313/bin/python -I scripts/dispute_halpha_analysis.py
Output: data/derived/dispute_halpha_results.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import verify_velocity_halpha_sim as hs          # noqa: E402
import verify_paper_kd as kdm                     # noqa: E402

from astropy.utils import iers                    # noqa: E402
iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = 'warn'
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun   # noqa: E402
from astropy.time import Time                     # noqa: E402
import astropy.units as u                         # noqa: E402

OUT = ROOT / 'data/derived/dispute_halpha_results.json'
C = 299792.458
REGS = ['P1', 'P2', 'P3', 'P4']

# ------------------------------------------------------------------ paper inputs
POS = {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72), 'P4': (115.20, 16.18)}   # Table 1
# Table 4 (p. 13): (V_helio, err, V_LSR printed)
T4 = {
    'P1': {'Ha': (16.3, 1.2, 4.7), 'NII': (11.8, 1.6, -0.9), 'SII': (11.8, 2.5, -1.0)},
    'P2': {'Ha': (17.0, 1.2, 4.0), 'NII': (12.0, 1.6, -0.5), 'SII': (12.0, 1.6, 0.5)},
    'P3': {'Ha': (18.4, 1.4, 6.1), 'NII': (11.7, 1.9, 0.4), 'SII': (11.7, 1.9, 1.7)},
    'P4': {'Ha': (19.7, 1.4, 4.5), 'NII': (13.3, 2.3, 1.0), 'SII': (12.7, 2.1, 1.1)},
}
# Table 4 Halpha distances as printed (MC mode, -err, +err; traditional), footnote unit "100 pc"
T4_HA_DIST = {'P1': (0.47, 0.45, 1.18, 0.61), 'P2': (0.31, 0.31, 1.44, 0.54),
              'P3': (0.34, 0.33, 1.20, 0.58), 'P4': (0.54, 0.52, 1.18, 0.80)}
# Table 2 (p. 7): Halpha relative intensity (units of the 6553 A sky line), Halpha and [N II] FWHM (km/s)
T2_IHA = {'P1': 0.30, 'P2': 0.37, 'P3': 0.26, 'P4': 0.25}
T2_FWHM_HA = {'P1': 31.9, 'P2': 31.0, 'P3': 32.2, 'P4': 32.4}
T2_FWHM_NII = {'P1': 27.0, 'P2': 26.1, 'P3': 30.0, 'P4': 26.9}

# ------------------------------------------------------------------ atomic data (NIST ASD v5.12, vacuum)
L_3P32, L_3P12 = 6564.537685, 6564.584403
L_GEO = (2 * L_3P32 + L_3P12) / 3          # Ly-beta pumping fills 3p by 2J+1 -> 2:1 (A-values equal to 1e-4)
L_REF = {'6564.61 (sim reference)': 6564.61, '6564.608 (NIST calc. blend)': 6564.608,
         '6564.632 (NIST Ritz, stat. weights)': 6564.632}

LOC = EarthLocation(lon=117.58 * u.deg, lat=40.39 * u.deg, height=960 * u.m)   # Xinglong (as verify_velocity)
APEX = SkyCoord('18h03m50.29s', '+30d00m16.8s', frame='icrs')                     # paper Sect. 3.2.3


def corr_lsr(ra, dec):
    """V_LSR - V_helio for 20 km/s toward the paper's apex."""
    return float(20.0 * np.cos(SkyCoord(ra * u.deg, dec * u.deg).separation(APEX).rad))


def wmean(x, s):
    x, s = np.asarray(x, float), np.asarray(s, float)
    w = 1 / s ** 2
    m = np.sum(w * x) / np.sum(w)
    return float(m), float(np.sqrt(1 / np.sum(w))), float(np.sum(((x - m) / s) ** 2)), float(np.sum((x / s) ** 2))


def fit_model(vv, A, mu, s, c0, c1, c2):
    return A * hs.gauss(vv, mu, s) + c0 + c1 * vv / 100 + c2 * (vv / 100) ** 2


def components(g, k, dm, vcorr, vneb, dv_fs, sig_inst=14.6, fwhm_neb=27.0, window_A=8.0):
    """Same construction as hs.simulate(), but returns every component and the full fit."""
    dvpix = 0.1 / hs.L0 * C
    v = np.arange(-window_A / hs.L0 * C, window_A / hs.L0 * C, dvpix)
    s_neb = np.hypot(fwhm_neb / 2.3548, sig_inst)
    s_geo = np.hypot(2.9, sig_inst)
    sep = C * (hs.L_3P12 - hs.L_3P32) / hs.L0
    wim = hs.gauss(v, vneb, s_neb)
    geo = g * (2 / 3 * hs.gauss(v, vcorr + dv_fs - sep / 3, s_geo) + 1 / 3 * hs.gauss(v, vcorr + dv_fs + 2 * sep / 3, s_geo))
    s_mod = np.hypot(sig_inst, 0.02 / hs.L0 * C)
    tmpl = k * g * hs.gauss(v, vcorr + dm, s_mod)
    y = wim + geo - tmpl
    p, _ = curve_fit(fit_model, v, y, p0=[1.0, vneb, s_neb, 0, 0, 0], maxfev=20000)
    fwhm_int = 2.3548 * np.sqrt(max(p[2] ** 2 - sig_inst ** 2, 0))
    return dict(v=v, wim=wim, geo=geo, tmpl=tmpl, resid=y, fit=fit_model(v, *p),
                fitline=p[0] * hs.gauss(v, p[1], abs(p[2])), p=p,
                bias=float(p[1] - vneb), flux_ratio=float(p[0]), fwhm_int_fit=float(fwhm_int),
                geo_centre=float(vcorr + dv_fs), tmpl_centre=float(vcorr + dm))


def run_sim(g, kt, dm, vcorr, vneb, dv_fs, sig_inst=14.6, fwhm_neb=27.0):
    """hs.simulate with module globals set per call; g = model (subtracted) geocoronal / WIM, kt = true/model."""
    hs.VCORR, hs.V_NEB = vcorr, vneb
    b, _ = hs.simulate(g=g * kt, k=1 / kt, dm=dm, sig_inst=sig_inst, fwhm_neb=fwhm_neb, dv_fs=dv_fs)
    return b


def run_sim_w(g, kt, dm, vcorr, vneb, dv_fs, sig_inst=14.6, fwhm_neb=27.0):
    """As run_sim, also returning the fitted intrinsic FWHM (instrument removed in quadrature, as Table 2)."""
    hs.VCORR, hs.V_NEB = vcorr, vneb
    b, s = hs.simulate(g=g * kt, k=1 / kt, dm=dm, sig_inst=sig_inst, fwhm_neb=fwhm_neb, dv_fs=dv_fs)
    return b, float(2.3548 * np.sqrt(max(s ** 2 - sig_inst ** 2, 0.0)))


def fit_grid(obs, sig, gmod, vcorr, vneb, dv_fs, dm_grid, kt_grid, width=None, fwhm_true=None, **kw):
    """chi2 grid over (dm, kt). width = {r: (FWHM, err)} adds the Table 2 Halpha widths to the fit;
    fwhm_true = {r: true intrinsic WIM FWHM} (default: kw fwhm_neb or 27)."""
    chi = np.zeros((len(dm_grid), len(kt_grid)))
    fw = fwhm_true or {r: kw.get('fwhm_neb', 27.0) for r in REGS}
    kw2 = {k: v for k, v in kw.items() if k != 'fwhm_neb'}
    for i, dm in enumerate(dm_grid):
        for j, kt in enumerate(kt_grid):
            c2 = 0.0
            for r in REGS:
                b, w = run_sim_w(gmod[r], kt, dm, vcorr[r], vneb[r], dv_fs, fwhm_neb=fw[r], **kw2)
                c2 += ((b - obs[r]) / sig[r]) ** 2
                if width is not None:
                    c2 += ((w - width[r][0]) / width[r][1]) ** 2
            chi[i, j] = c2
    i, j = np.unravel_index(np.argmin(chi), chi.shape)
    prof = chi.min(axis=0)                     # profile over dm
    ok = kt_grid[prof <= chi.min() + 1.0]
    ndat = 4 * (2 if width is not None else 1)
    best = {r: run_sim_w(gmod[r], kt_grid[j], dm_grid[i], vcorr[r], vneb[r], dv_fs, fwhm_neb=fw[r], **kw2) for r in REGS}
    return dict(dm=float(dm_grid[i]), kt=float(kt_grid[j]), chi2=round(float(chi.min()), 2),
                dof=ndat - (2 if len(dm_grid) > 1 else 1),
                kt_1sigma=[round(float(ok.min()), 2), round(float(ok.max()), 2)],
                kt_1sigma_at_grid_edge=bool(ok.min() <= kt_grid[0] + 1e-9 or ok.max() >= kt_grid[-1] - 1e-9),
                oversub_pct=round(float(100 * (1 / kt_grid[j] - 1)), 1),
                oversub_pct_1sigma=[round(float(100 * (1 / ok.max() - 1)), 1), round(float(100 * (1 / ok.min() - 1)), 1)],
                pred_shift={r: round(best[r][0], 2) for r in REGS},
                pred_fwhm_int={r: round(best[r][1], 1) for r in REGS},
                fwhm_true_used={r: fw[r] for r in REGS})


def flux_factor(g, kt, dm, vcorr, vneb, dv_fs, fwhm):
    """Fitted Halpha flux / true WIM flux (same Gaussian + quadratic fit applied to the WIM line alone)."""
    c = components(g * kt, 1 / kt, dm, vcorr, vneb, dv_fs, fwhm_neb=fwhm)
    t = components(0.0, 1.0, 0.0, vcorr, vneb, dv_fs, fwhm_neb=fwhm)
    return c['flux_ratio'] / t['flux_ratio']


def selfcons_g(gm, kt, dm, vcorr, vneb, dv_fs, fwhm, n_iter=40, tol=1e-3):
    """Self-consistent subtracted-sky / true-WIM ratio g = gm * f(g), where gm = 0.319 / (Table 2 Halpha) uses the
    published (fitted) flux and f = fitted / true flux.  Damped fixed-point iteration; returns (g, f, converged)."""
    f = 1.0
    for _ in range(n_iter):
        fn = flux_factor(gm * f, kt, dm, vcorr, vneb, dv_fs, fwhm)
        if abs(fn - f) < tol:
            return gm * fn, fn, True
        f = 0.5 * (f + fn)
    return gm * f, f, False


def fit_selfcons(obs, sig, gmod, vcorr, vneb, dv_fs, dm_grid, kt_grid, width=None, fwhm_true=None, n_iter=8, tol=3e-3):
    """fit_grid with the sky/nebula ratio made self-consistent: g_r = gmod_r * f_r(best fit), iterated."""
    fw = fwhm_true or {r: 27.0 for r in REGS}
    f = {r: 1.0 for r in REGS}
    for it in range(n_iter):
        g = {r: gmod[r] * f[r] for r in REGS}
        P = fit_grid(obs, sig, g, vcorr, vneb, dv_fs, dm_grid, kt_grid, width=width, fwhm_true=fw)
        fn = {r: flux_factor(g[r], P['kt'], P['dm'], vcorr[r], vneb[r], dv_fs, fw[r]) for r in REGS}
        delta = max(abs(fn[r] - f[r]) for r in REGS)
        f = fn
        if delta < tol:
            break
    P['g_selfcons'] = {r: round(g[r], 3) for r in REGS}
    P['flux_factor'] = {r: round(f[r], 3) for r in REGS}
    P['true_Ha_rel_intensity'] = {r: round(T2_IHA[r] / f[r], 3) for r in REGS}
    P['iterations'] = it + 1
    P['converged'] = bool(delta < tol)
    return P, g


def main():
    out = {'_provenance':'Wu+2026 Tables 1, 2, 4 (arXiv v1 pp. 7, 13); project files listed in the docstring; '
                          'all other numbers computed here'}

    # ---------------- 1. offsets inside the paper's own red arm
    off = {}
    for r in REGS:
        ha, n2, s2 = T4[r]['Ha'], T4[r]['NII'], T4[r]['SII']
        off[r] = dict(Ha_minus_NII=round(ha[0] - n2[0], 2), err_NII=round(np.hypot(ha[1], n2[1]), 2),
                      Ha_minus_SII=round(ha[0] - s2[0], 2), err_SII=round(np.hypot(ha[1], s2[1]), 2),
                      NII_minus_SII=round(n2[0] - s2[0], 2), err_NS=round(np.hypot(n2[1], s2[1]), 2))
    m, e, chi_m, chi0 = wmean([off[r]['Ha_minus_NII'] for r in REGS], [off[r]['err_NII'] for r in REGS])
    m2, e2, chi_m2, chi02 = wmean([off[r]['Ha_minus_SII'] for r in REGS], [off[r]['err_SII'] for r in REGS])
    m3, e3, chi_m3, chi03 = wmean([off[r]['NII_minus_SII'] for r in REGS], [off[r]['err_NS'] for r in REGS])
    out['offsets_table4'] = dict(per_region=off,
                                 Ha_minus_NII=dict(wmean=round(m, 2), err=round(e, 2), z=round(m / e, 1),
                                                   chi2_about_mean=round(chi_m, 2), chi2_zero=round(chi0, 1)),
                                 Ha_minus_SII=dict(wmean=round(m2, 2), err=round(e2, 2), z=round(m2 / e2, 1),
                                                   chi2_zero=round(chi02, 1)),
                                 NII_minus_SII=dict(wmean=round(m3, 2), err=round(e3, 2), chi2_zero=round(chi03, 2)),
                                 note='errors: Table 4 RV errors in quadrature, treated as independent')
    # rest-wavelength convention alone
    out['rest_wavelength_convention_kms'] = round(C * (6564.632 / 6564.61 - 1), 2)

    # ---------------- 2. LSR corrections and WHAM per region
    corr = {r: corr_lsr(*POS[r]) for r in REGS}
    rows = list(csv.DictReader(open(ROOT / 'data/derived/verify_hiwham_wham_beams.csv')))
    bra = np.array([float(x['ra']) for x in rows]); bdec = np.array([float(x['dec']) for x in rows])
    bc = SkyCoord(bra * u.deg, bdec * u.deg)
    vg = np.array([float(x['v_g']) for x in rows]); fg = np.array([float(x['fwhm_g']) for x in rows])
    ig = np.array([float(x['I_g']) for x in rows]); it = np.array([float(x['I_tot']) for x in rows])
    wham = {}
    for r in REGS:
        sep = bc.separation(SkyCoord(*POS[r], unit='deg')).deg
        n = int(np.argmin(sep)); w = sep <= 1.0
        wham[r] = dict(nearest_sep_deg=round(float(sep[n]), 3), v_lsr_nearest=round(float(vg[n]), 2),
                       v_lsr_within1deg_mean=round(float(vg[w].mean()), 2),
                       v_lsr_within1deg_std=round(float(vg[w].std(ddof=1)), 2), n_within1deg=int(w.sum()),
                       fwhm_int_nearest=round(float(np.sqrt(fg[n] ** 2 - 12.0 ** 2)), 1),
                       I_gauss_nearest=round(float(ig[n]), 2), I_tot_nearest=round(float(it[n]), 2),
                       corr_lsr=round(corr[r], 2),
                       v_helio_nearest=round(float(vg[n]) - corr[r], 2),
                       v_helio_within1deg=round(float(vg[w].mean()) - corr[r], 2))
    out['wham_per_region'] = wham
    # recomputed V_LSR of the paper's lines with the paper's own convention
    out['lamost_vlsr_recomputed'] = {r: {ln: round(T4[r][ln][0] + corr[r], 2) for ln in ('Ha', 'NII', 'SII')}
                                     for r in REGS}
    out['Ha_minus_WHAM'] = {r: dict(dv=round(T4[r]['Ha'][0] - wham[r]['v_helio_nearest'], 2),
                                    err=round(np.hypot(T4[r]['Ha'][1], max(wham[r]['v_lsr_within1deg_std'], 0.5)), 2))
                            for r in REGS}
    out['NII_minus_WHAM'] = {r: dict(dv=round(T4[r]['NII'][0] - wham[r]['v_helio_nearest'], 2),
                                     err=round(np.hypot(T4[r]['NII'][1], max(wham[r]['v_lsr_within1deg_std'], 0.5)), 2))
                             for r in REGS}

    for key in ('Ha_minus_WHAM', 'NII_minus_WHAM'):
        mm, ee, cc, _ = wmean([out[key][r]['dv'] for r in REGS], [out[key][r]['err'] for r in REGS])
        out[key]['wmean'] = dict(mean=round(mm, 2), err=round(ee, 2), chi2_about_mean=round(cc, 2),
                                 note='WHAM error = beam-to-beam std within 1 deg (floor 0.5) in quadrature with Table 4')
    # field-wide: digitized Fig. 5 vs WHAM beams inside the plate (HIP 37797, r = 2.5 deg; notes/geometry.md)
    f5 = np.genfromtxt(ROOT / 'data/derived/geometry_lamost_halpha_vlsr_fig5_digitized.csv', delimiter=',',
                       skip_header=1)
    plate = SkyCoord(116.25248 * u.deg, 15.89500 * u.deg)
    inpl = bc.separation(plate).deg <= 2.5
    out['field_wide'] = dict(fig5_n=int(len(f5)), fig5_median_vlsr=round(float(np.median(f5[:, 2])), 2),
                             fig5_layer2_median=round(float(np.median(f5[f5[:, 3] == 1, 2])), 2),
                             fig5_censored_low=int((f5[:, 4] == -1).sum()), fig5_censored_high=int((f5[:, 4] == 1).sum()),
                             wham_in_plate_n=int(inpl.sum()), wham_in_plate_median_vlsr=round(float(np.median(vg[inpl])), 2),
                             wham_in_plate_p16_p84=[round(float(np.percentile(vg[inpl], 16)), 2),
                                                    round(float(np.percentile(vg[inpl], 84)), 2)])

    # ---------------- 3. geocoronal geometry on 2024-02-01
    lam = json.load(open(ROOT / 'data/derived/verify_velocity_lamost_arms.json'))['helio']
    times = sorted({e['mid_utc'] for e in lam['exposures']})
    t = Time(times, scale='utc')
    sun = get_sun(t).transform_to(AltAz(obstime=t, location=LOC))
    sunalt = sun.alt.deg
    ratio = 0.954 - 0.011 * np.abs(sunalt)
    vcorr = {}
    for r in REGS:
        sc = SkyCoord(*POS[r], unit='deg')
        vcorr[r] = float(np.mean(sc.radial_velocity_correction(kind='heliocentric', obstime=t, location=LOC)
                                 .to(u.km / u.s).value))
    dv_fs = {k: round(C * (L_GEO / v - 1), 3) for k, v in L_REF.items()}
    DVFS = C * (L_GEO / hs.L0 - 1)            # on the simulation's 6564.61 A scale
    inst_fwhm_zhang = 2.3548 * 0.32 / 6554.0 * C
    out['geocoronal'] = dict(
        exposure_mid_utc=times, sun_altitude_deg=[round(float(a), 1) for a in sunalt],
        model_ratio_per_exposure=[round(float(x), 3) for x in ratio], model_ratio_mean=round(float(ratio.mean()), 3),
        header_HELIO_RV_median=lam['header_HELIO_RV_median'], header_HELIO_RV_range=lam['header_HELIO_RV_range'],
        v_corr_per_region={r: round(v, 2) for r, v in vcorr.items()},
        lambda_geo_vac=round(L_GEO, 6), dv_fs_kms=dv_fs, dv_fs_sim_original=round(hs.DV_FS, 3),
        v_helio_geo_per_region={r: round(vcorr[r] + DVFS, 2) for r in REGS},
        separation_wim_geo={r: round(T4[r]['NII'][0] - (vcorr[r] + DVFS), 1) for r in REGS},
        RVsep_NII_minus_6554={r: round(T4[r]['NII'][0] - vcorr[r], 1) for r in REGS},
        RVsep_threshold_zhang=25.0, inst_fwhm_kms_zhang_sigma032=round(inst_fwhm_zhang, 1),
        inst_fwhm_kms_R7500=round(C / 7500, 1),
        g_model={r: round(float(ratio.mean()) / T2_IHA[r], 3) for r in REGS})
    gmod = {r: float(ratio.mean()) / T2_IHA[r] for r in REGS}

    # ---------------- 4. forward model fits
    obs = {r: off[r]['Ha_minus_NII'] for r in REGS}
    sig_reg = {r: off[r]['err_NII'] for r in REGS}
    dm_g = np.arange(-12, 8.01, 0.25)
    kt_g = np.arange(0.50, 1.501, 0.01)
    # A: exact reproduction of verify_velocity_halpha_sim (uniform inputs)
    gA = {r: (0.954 - 0.011 * 56.4) / T2_IHA[r] for r in REGS}
    A = fit_grid(obs, {r: 1.9 for r in REGS}, gA, {r: -8.32 for r in REGS}, {r: 12.2 for r in REGS}, hs.DV_FS,
                 np.arange(-12, 4.01, 0.5), np.arange(0.5, 1.81, 0.05))
    vneb = {r: T4[r]['NII'][0] for r in REGS}
    B = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g)
    Cv = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g, fwhm_neb=31.0)
    D = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g, sig_inst=17.0)
    E = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, np.array([0.0]), kt_g)
    t2w = {r: (T2_FWHM_HA[r], {'P1': 4.0, 'P2': 3.0, 'P3': 4.2, 'P4': 3.4}[r]) for r in REGS}
    fw_wham = {r: wham[r]['fwhm_int_nearest'] for r in REGS}
    F = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g, width=t2w, fwhm_true=fw_wham)
    F2 = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g, width=t2w)
    G = fit_grid(obs, sig_reg, gmod, vcorr, vneb, DVFS, np.array([0.0]), kt_g, width=t2w, fwhm_true=fw_wham)
    # self-consistent sky/nebula ratio (review fix): the subtracted template is 0.319 I(6554) while the published
    # Table 2 Halpha is the *fitted* flux, which the model lowers by f; the true WIM flux is Table 2 / f.
    Bsc, gB = fit_selfcons(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g)
    Fsc, gF = fit_selfcons(obs, sig_reg, gmod, vcorr, vneb, DVFS, dm_g, kt_g, width=t2w, fwhm_true=fw_wham)
    Gsc, gG = fit_selfcons(obs, sig_reg, gmod, vcorr, vneb, DVFS, np.array([0.0]), kt_g, width=t2w, fwhm_true=fw_wham)
    # widths predicted by the centroid-only best fit B (not fitted)
    out['fits'] = {'A_reproduce_original': A, 'B_refined_inputs': B, 'C_B_fwhm31': Cv, 'D_B_sig17': D,
                   'E_B_dm0_only': E, 'F_joint_centroid_width_WHAMfwhm': F, 'F2_joint_centroid_width_fwhm27': F2,
                   'G_F_dm0_only': G,
                   'Bsc_centroid_only_selfcons': Bsc, 'Fsc_joint_selfcons': Fsc, 'Gsc_joint_dm0_selfcons': Gsc,
                   'table2_Ha_fwhm': t2w,
                   '_notes': 'A: verify_velocity_halpha_sim settings (sigma 1.9, vcorr -8.32, V_NEB 12.2, DV_FS -2.40). '
                             'B: per-region astropy v_corr, V_NEB=[N II] V_helio, NIST DV_FS, Table 4 errors (Ha+[N II] '
                             'in quadrature), FWHM_neb 27, sig_inst 14.6. C: B with FWHM_neb 31. D: B with sig_inst 17.0 '
                             '(R=7500). E: B with template at RV(6554) (dm=0), kt only. F: B fitted jointly to the Table 4 '
                             'offsets and the Table 2 Halpha FWHM, true WIM FWHM = WHAM nearest-beam intrinsic FWHM. '
                             'F2: as F with true FWHM 27. dm = template offset from topocentric rest (km/s); '
                             'kt = true/model geocoronal intensity; over-subtraction = 1/kt - 1. '
                             'Bsc/Fsc/Gsc: B/F/G with the sky/nebula ratio made self-consistent (g = 0.319 f / '
                             'Table 2 Halpha, f = fitted/true Halpha flux of the best fit, iterated). Fsc is the '
                             'adopted model of the Comment (Sect. 4).'}
    # part of the offset left unexplained by the width-consistent fits, and the implied sky-ratio error
    for key in ('B_refined_inputs', 'F_joint_centroid_width_WHAMfwhm', 'G_F_dm0_only',
                'Bsc_centroid_only_selfcons', 'Fsc_joint_selfcons', 'Gsc_joint_dm0_selfcons'):
        P = out['fits'][key]
        mm, ee, cc, _ = wmean([obs[r] - P['pred_shift'][r] for r in REGS], [sig_reg[r] for r in REGS])
        P['unexplained_offset'] = dict(wmean=round(mm, 2), err=round(ee, 2), z=round(mm / ee, 1))
        P['explained_fraction_of_wmean'] = round(1 - mm / m, 2)
        P['true_ratio_implied'] = round(float(ratio.mean()) * P['kt'], 3)
    # sign test: correct or under-estimated model intensity (kt >= 1), any template offset
    st = {}
    for gl, gset in (('g_published', gmod), ('g_selfcons_Fsc', gF)):
        for kt in (1.0, 1.1, 1.25, 1.5):
            for dm in (-16, -12, -8.24, -4, -2, 0, 2, 4, 8, 12, 16, 20):
                st[f'{gl}_kt{kt}_dm{dm}'] = [round(run_sim(gset[r], kt, dm, vcorr[r], vneb[r], DVFS), 2) for r in REGS]
    allb = [b for v in st.values() for b in v]
    out['sign_test_kt_ge1'] = dict(bias=st, max_bias=max(allb), min_bias=min(allb),
                                   note='kt = true/model >= 1 means no over-subtraction; both sky/nebula normalisations')
    # curves for Fig. 4b/4c: centroid shift and fitted intrinsic FWHM vs over-subtraction at the template position
    # of the adopted joint fit Fsc (true WIM FWHM = WHAM), with the sky/nebula ratio made self-consistent at every
    # point (damped iteration; points that do not converge are stored as null)
    fgrid = np.linspace(-0.30, 0.80, 45)
    cF = {}
    for r in REGS:
        cF[r] = []
        for f in fgrid:
            gs, fs, ok = selfcons_g(gmod[r], 1 / (1 + f), Fsc['dm'], vcorr[r], vneb[r], DVFS, fw_wham[r])
            b_, w_ = run_sim_w(gs, 1 / (1 + f), Fsc['dm'], vcorr[r], vneb[r], DVFS, fwhm_neb=fw_wham[r])
            cF[r].append((b_, w_, gs, ok))
    out['fig4b'] = dict(oversub_frac=[round(float(f), 4) for f in fgrid],
                        dm_F=Fsc['dm'], dm_B=Bsc['dm'], fit_label='Fsc (joint, self-consistent)',
                        shift_F={r: [round(x[0], 3) if x[3] else None for x in cF[r]] for r in REGS},
                        fwhm_F={r: [round(x[1], 2) if x[3] else None for x in cF[r]] for r in REGS},
                        g_F={r: [round(x[2], 3) if x[3] else None for x in cF[r]] for r in REGS},
                        kt1_template_literal_dm_m8p24={r: round(run_sim(gmod[r], 1.0, -8.24, vcorr[r], vneb[r], DVFS), 2)
                                                       for r in REGS},
                        kt1_dm0={r: round(run_sim(gmod[r], 1.0, 0.0, vcorr[r], vneb[r], DVFS), 2) for r in REGS})

    # profiles for Fig. 4a (P1, adopted fit Fsc) and flux / width side effects for each region (fits B, F, Fsc)
    side, sideF, sideFsc = {}, {}, {}
    for r in REGS:
        for lab, P, gg, fwt, dst in (('B', B, gmod[r], 27.0, side), ('F', F, gmod[r], fw_wham[r], sideF),
                                     ('Fsc', Fsc, gF[r], fw_wham[r], sideFsc)):
            cmp_ = components(gg * P['kt'], 1 / P['kt'], P['dm'], vcorr[r], vneb[r], DVFS, fwhm_neb=fwt)
            ref = run_sim(gg, P['kt'], P['dm'], vcorr[r], vneb[r], DVFS, fwhm_neb=fwt)
            assert abs(cmp_['bias'] - ref) < 1e-6, (cmp_['bias'], ref)
            tru = components(0.0, 1.0, 0.0, vcorr[r], vneb[r], DVFS, fwhm_neb=fwt)   # WIM alone, same fit
            w_fit, w_tru = cmp_['fwhm_int_fit'], tru['fwhm_int_fit']
            dst[r] = dict(bias=round(cmp_['bias'], 2), g=round(gg, 3),
                          flux_fit_over_true=round(cmp_['flux_ratio'] / tru['flux_ratio'], 3),
                          fwhm_int_fit=round(w_fit, 2), fwhm_int_true=round(w_tru, 2),
                          dT_eq2_K=round(23.5 * (w_fit ** 2 - w_tru ** 2), -1))
        side[r]['T_eq2_table2_K'] = round(23.5 * (T2_FWHM_HA[r] ** 2 - T2_FWHM_NII[r] ** 2), -1)
        cmp_ = components(gF[r] * Fsc['kt'], 1 / Fsc['kt'], Fsc['dm'], vcorr[r], vneb[r], DVFS, fwhm_neb=fw_wham[r])
        if r == 'P1':
            out['fig4a'] = dict(region='P1', fit_label='Fsc (joint centroid + width, self-consistent)',
                                g_model=round(gmod[r], 3), g_selfcons=round(gF[r], 3),
                                kt=Fsc['kt'], dm=Fsc['dm'], fwhm_true=fw_wham[r],
                                fwhm_int_fit=round(cmp_['fwhm_int_fit'], 2),
                                v=[round(float(x), 3) for x in cmp_['v']],
                                **{k: [round(float(x), 6) for x in cmp_[k]] for k in ('wim', 'geo', 'tmpl', 'resid', 'fit', 'fitline')},
                                geo_centre=round(cmp_['geo_centre'], 2), tmpl_centre=round(cmp_['tmpl_centre'], 2),
                                wim_centre=vneb[r], fit_centre=round(vneb[r] + cmp_['bias'], 2),
                                table4_Ha=T4[r]['Ha'][0], table4_NII=T4[r]['NII'][0])
    out['side_effects_best_fit'] = side
    out['side_effects_joint_fit_F'] = sideF
    out['side_effects_joint_fit_Fsc'] = sideFsc
    # width sensitivity to the template position at fixed over-subtraction (P1, kt of fit B)
    out['width_vs_template_offset_P1'] = {str(dm): [round(x, 2) for x in run_sim_w(gmod['P1'], B['kt'], dm, vcorr['P1'],
                                                                                    vneb['P1'], DVFS)]
                                          for dm in (-4, -2, 0, 2, 4)}
    # solar Halpha absorption term from the original simulation (C_sun in WIM-line flux per A)
    sim0 = json.load(open(ROOT / 'data/derived/verify_velocity_halpha_sim.json'))
    out['solar_halpha_from_verify_sim'] = sim0['solar_Halpha_absorption']

    # ---------------- 5. kinematic distances inherited by the Halpha rows (kd form, Reid+14 URC)
    kdist = {}
    for r in REGS:
        sc = SkyCoord(*POS[r], unit='deg').galactic
        l, b = float(sc.l.deg), float(sc.b.deg)
        v_ha_tab = T4[r]['Ha'][2]
        v_ha_iau = float(kdm.helio_to_lsr(l, b, T4[r]['Ha'][0]))
        v_ha_cor = float(kdm.helio_to_lsr(l, b, T4[r]['Ha'][0] - side[r]['bias']))
        v_ha_corF = float(kdm.helio_to_lsr(l, b, T4[r]['Ha'][0] - sideF[r]['bias']))
        v_ha_corFsc = float(kdm.helio_to_lsr(l, b, T4[r]['Ha'][0] - sideFsc[r]['bias']))
        v_n2 = float(kdm.helio_to_lsr(l, b, T4[r]['NII'][0]))
        v_wh = wham[r]['v_lsr_nearest']
        dd = lambda v: (None if not np.isfinite(kdm.trad_kd(l, b, v)) else round(float(kdm.trad_kd(l, b, v)), 3))
        grad = float((kdm.vlsr_model(l, b, 1.0) - kdm.vlsr_model(l, b, 0.0001)) / 1.0)
        kdist[r] = dict(l=round(l, 3), b=round(b, 3), table4_trad=T4_HA_DIST[r][3], table4_mc=T4_HA_DIST[r][:3],
                        V_Ha_table=v_ha_tab, d_Ha_table=dd(v_ha_tab), V_Ha_iau=round(v_ha_iau, 2), d_Ha_iau=dd(v_ha_iau),
                        V_Ha_corrected_fitB=round(v_ha_cor, 2), d_Ha_corrected_fitB=dd(v_ha_cor),
                        V_Ha_corrected_fitF=round(v_ha_corF, 2), d_Ha_corrected_fitF=dd(v_ha_corF),
                        V_Ha_corrected_fitFsc=round(v_ha_corFsc, 2), d_Ha_corrected_fitFsc=dd(v_ha_corFsc),
                        V_NII_iau=round(v_n2, 2), d_NII=dd(v_n2), V_WHAM=v_wh, d_WHAM=dd(v_wh),
                        grad_kms_per_kpc_0_1=round(grad, 2),
                        distance_shift_per_observed_offset_kpc=round(off[r]['Ha_minus_NII'] / grad, 2))
    out['kinematic_distances_kpc'] = kdist
    OUT.write_text(json.dumps(out, indent=1))

    # console summary
    print(json.dumps(out['offsets_table4'], indent=1))
    print('WHAM', json.dumps(wham, indent=0))
    print('field', out['field_wide'])
    print('geo', json.dumps(out['geocoronal'], indent=0))
    for k in ('A_reproduce_original', 'B_refined_inputs', 'C_B_fwhm31', 'D_B_sig17', 'E_B_dm0_only',
              'F_joint_centroid_width_WHAMfwhm', 'F2_joint_centroid_width_fwhm27', 'G_F_dm0_only',
              'Bsc_centroid_only_selfcons', 'Fsc_joint_selfcons', 'Gsc_joint_dm0_selfcons'):
        print(k, out['fits'][k])
    print('sign', out['sign_test_kt_ge1']['max_bias'], out['sign_test_kt_ge1']['min_bias'])
    print('kt1', out['fig4b']['kt1_dm0'], out['fig4b']['kt1_template_literal_dm_m8p24'])
    print('side', json.dumps(side, indent=0))
    print('sideF', json.dumps(sideF, indent=0))
    print('sideFsc', json.dumps(sideFsc, indent=0))
    print('width vs dm', out['width_vs_template_offset_P1'])
    print('kd', json.dumps(kdist, indent=0))


if __name__ == '__main__':
    main()
