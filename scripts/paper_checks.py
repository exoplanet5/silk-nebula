#!/usr/bin/env python
"""paper_review: independent reproduction of every quantitative step in
Wu et al. 2026, "The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament",
arXiv:2610.13329v1 (AJ, doi:10.3847/1538-3881/aeb11c).

Writes data/derived/paper_checks.json (+ MC sample summaries used by paper_review_figures.py).

Inputs (read-only):
  * Paper Tables 1-4 transcribed below from sources/arxiv/2610.13329v1.pdf (pages 6, 7, 13).
  * WHAM-SS DR1 spectra  data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits (teammate download)
  * Finkbeiner (2003) H-alpha map data/raw/lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits (teammate download)
  * kd source code (T. Wenger, github.com/tvwenger/kd) read as TEXT ONLY in data/raw/github_kd/ to
    emulate the rotation curves used by the paper (code is re-implemented here, never executed).

Run:  ~/.venvs/astro313/bin/python -I scripts/paper_checks.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, brentq
from scipy.stats import gaussian_kde

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import ROOT, REGIONS, gal  # noqa: E402
from paper_review_lib import (C_KMS, U_STD, V_STDc, W_STD, unit_vec, vlsr_kd14, MODELS,  # noqa: E402
                              vs_dopita, e_boumis, r_pds, age_cioffi, ne_proxauf, ne_sanders, R_sanders,
                              SA, SB_, SC, OMEGA_1D, G_3P, T_EXC, BR5007, OH, q_oiii, em_from_oiii,
                              em_from_ha, oiii_ha_photon_ratio, mag_to_R)

import astropy.units as u  # noqa: E402
from astropy.coordinates import SkyCoord, EarthLocation, LSR, AltAz  # noqa: E402
from astropy.time import Time  # noqa: E402

RNG = np.random.default_rng(20261116)
OUT = ROOT / 'data/derived/paper_checks.json'
OUT_MC = ROOT / 'data/derived/paper_review_mc_samples.npz'
R = {}  # results dictionary written to JSON

# ----------------------------------------------------------------------------------------------
# 0. Paper values (transcribed from the PDF; Table 2 p.7, Table 3 p.6-7, Table 4 p.13)
# ----------------------------------------------------------------------------------------------
T2 = {  # relative intensities (to sky line), sigma, S/N ; FWHM intrinsic (km/s) and sigma
    'P1': dict(S6716=(0.14, 0.02), S6731=(0.13, 0.02), N6548=(0.09, 0.02), N6583=(0.24, 0.02),
               Ha=(0.30, 0.02), O5007=1.67, O4959=0.73, SN_O5007=6.1, SN_O4959=1.8,
               W_S=(29.0, 6.1), W_N=(27.0, 4.6), W_H=(31.9, 4.0), W_O=46.9,
               ratio_SII_Ha=0.93, ratio_NII_Ha=1.11, ratio_SII=1.10, SN_red=dict(S=55, N=47, H=50)),
    'P2': dict(S6716=(0.17, 0.02), S6731=(0.15, 0.02), N6548=(0.06, 0.02), N6583=(0.26, 0.03),
               Ha=(0.37, 0.02), O5007=1.83, O4959=0.75, SN_O5007=6.1, SN_O4959=2.0,
               W_S=(31.4, 6.4), W_N=(26.1, 4.7), W_H=(31.0, 3.0), W_O=48.9,
               ratio_SII_Ha=0.88, ratio_NII_Ha=0.88, ratio_SII=1.14, SN_red=dict(S=52, N=48, H=45)),
    'P3': dict(S6716=(0.11, 0.02), S6731=(0.10, 0.02), N6548=(0.07, 0.02), N6583=(0.20, 0.02),
               Ha=(0.26, 0.02), O5007=1.90, O4959=0.75, SN_O5007=6.6, SN_O4959=2.3,
               W_S=(30.1, 6.1), W_N=(30.0, 6.0), W_H=(32.2, 4.2), W_O=43.1,
               ratio_SII_Ha=0.83, ratio_NII_Ha=1.02, ratio_SII=1.19, SN_red=dict(S=65, N=53, H=53)),
    'P4': dict(S6716=(0.10, 0.02), S6731=(0.10, 0.02), N6548=(0.08, 0.02), N6583=(0.20, 0.02),
               Ha=(0.25, 0.02), O5007=None, O4959=None, SN_O5007=None, SN_O4959=None,
               W_S=(32.9, 6.4), W_N=(26.9, 4.8), W_H=(32.4, 3.4), W_O=None,
               ratio_SII_Ha=0.78, ratio_NII_Ha=0.95, ratio_SII=1.03, SN_red=dict(S=44, N=54, H=45)),
}
T3_SB = {'P1': (23.60, 0.028), 'P2': (23.87, 0.035), 'P3': (23.85, 0.039)}  # mag/arcsec^2
# Table 4: line -> (V_helio, sigma, V_LSR, D_MC, D_MC_minus, D_MC_plus, D_trad); V_t per region
T4 = {
    'P1': {'OIII': (26.8, 4.3, 15.5, 1.63, 0.84, 1.01, 1.79), 'SII': (11.8, 2.5, -1.0, 0.06, 0.05, 1.17, 0.16),
           'NII': (11.8, 1.6, -0.9, 0.11, 0.10, 1.12, 0.18), 'Ha': (16.3, 1.2, 4.7, 0.47, 0.45, 1.18, 0.61)},
    'P2': {'OIII': (26.8, 4.1, 14.5, 1.58, 0.83, 1.23, 1.85), 'SII': (12.0, 1.6, 0.5, 0.21, 0.20, 1.06, 0.17),
           'NII': (12.0, 1.6, -0.5, 0.08, 0.07, 1.35, 0.17), 'Ha': (17.0, 1.2, 4.0, 0.31, 0.31, 1.44, 0.54)},
    'P3': {'OIII': (23.8, 3.3, 11.2, 0.98, 0.69, 1.08, 1.20), 'SII': (11.7, 1.9, 1.7, 0.21, 0.21, 1.12, 0.23),
           'NII': (11.7, 1.9, 0.4, 0.09, 0.08, 1.02, 0.14), 'Ha': (18.4, 1.4, 6.1, 0.34, 0.33, 1.20, 0.58)},
    'P4': {'SII': (12.7, 2.1, 1.1, 0.22, 0.21, 1.61, 0.22), 'NII': (13.3, 2.3, 1.0, 0.28, 0.27, 0.93, 0.49),
           'Ha': (19.7, 1.4, 4.5, 0.54, 0.52, 1.18, 0.80)},
}
T4_VT = {'P1': 33.1, 'P2': 32.9, 'P3': 30.6, 'P4': 29.9}
R['inputs'] = {'paper': 'Wu et al. 2026, arXiv:2610.13329v1, AJ doi:10.3847/1538-3881/aeb11c',
               'table2': T2, 'table3_sb_mag_arcsec2': T3_SB, 'table4': T4, 'table4_Vt': T4_VT,
               'regions_radec': REGIONS}

LB = {k: tuple(float(x[0]) for x in gal(*v)) for k, v in REGIONS.items()}
R['regions_lb'] = LB


# ----------------------------------------------------------------------------------------------
# 1. LSR conversion
# ----------------------------------------------------------------------------------------------
apex = SkyCoord('18h03m50.29s', '+30d00m16.8s', frame='icrs').galactic
V_STD = 20.0 * unit_vec(apex.l.deg, apex.b.deg)           # standard solar motion (U,V,W)
V_S10 = np.array([11.1, 12.24, 7.25])                      # Schoenrich, Binney & Dehnen 2010
lsr = {'standard_apex_lb': (apex.l.deg, apex.b.deg), 'standard_UVW': V_STD.tolist(),
       'schoenrich2010_UVW': V_S10.tolist(), 'regions': {}}
for reg, (l, b) in LB.items():
    n = unit_vec(l, b)
    d_std, d_s10 = float(V_STD @ n), float(V_S10 @ n)
    # astropy cross-check (astropy LSR frame default = Schoenrich 2010)
    ra, dec = REGIONS[reg]
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, distance=1 * u.kpc, pm_ra_cosdec=0 * u.mas / u.yr,
                 pm_dec=0 * u.mas / u.yr, radial_velocity=0 * u.km / u.s, frame='icrs')
    d_astropy = float(c.transform_to(LSR()).radial_velocity.to_value(u.km / u.s))
    rows = {}
    for line, (vh, sv, vl, *_r) in T4[reg].items():
        rows[line] = {'V_helio': vh, 'V_LSR_paper': vl, 'paper_offset': round(vl - vh, 2),
                      'V_LSR_std_ours': round(vh + d_std, 2), 'V_LSR_S10_ours': round(vh + d_s10, 2),
                      'paper_minus_ours_std': round(vl - (vh + d_std), 2)}
    lsr['regions'][reg] = {'l': l, 'b': b, 'corr_std': round(d_std, 3), 'corr_S10': round(d_s10, 3),
                           'corr_S10_astropy': round(d_astropy, 3), 'lines': rows}
offs = [r['paper_minus_ours_std'] for g in lsr['regions'].values() for r in g['lines'].values()]
lsr['summary'] = {'expected_corr_std_range': [min(g['corr_std'] for g in lsr['regions'].values()),
                                              max(g['corr_std'] for g in lsr['regions'].values())],
                  'paper_minus_ours_std_range': [min(offs), max(offs)],
                  'paper_minus_ours_std_rms': float(np.sqrt(np.mean(np.square(offs))))}

# Barycentric correction and geocoronal H-alpha position for the LAMOST observation (2024-02-01)
xinglong = EarthLocation(lon=117.5775 * u.deg, lat=40.3933 * u.deg, height=960 * u.m)
tgt = SkyCoord(116.68 * u.deg, 16.42 * u.deg)
times = Time('2024-02-01T12:00:00', scale='utc') + np.linspace(0, 10, 41) * u.hour
alt = tgt.transform_to(AltAz(obstime=times, location=xinglong)).alt.deg
t_transit = times[np.argmax(alt)]
vb = tgt.radial_velocity_correction('barycentric', obstime=t_transit, location=xinglong).to_value(u.km / u.s)
lsr['observation'] = {'site': 'Xinglong (LAMOST) 117.5775E 40.3933N 960 m (assumed)',
                      'transit_utc': t_transit.isot, 'transit_alt_deg': float(alt.max()),
                      'barycentric_correction_kms': round(float(vb), 2),
                      'geocoronal_Ha_V_helio_kms': round(float(vb), 2),
                      'geocoronal_Ha_V_LSR_std_kms': round(float(vb) + lsr['regions']['P1']['corr_std'], 2)}
R['lsr'] = lsr

# ----------------------------------------------------------------------------------------------
# 2. Kinematic distances
# ----------------------------------------------------------------------------------------------
# rotation-curve models (flat, Reid+14/19 linear, kd URC emulations): see paper_review_lib.MODELS


def solve_distance(fn, l, b, v, dmax=10.0):
    """Monotonic solve V_model(d) = v for 0 < d < dmax (outer Galaxy, l~204). NaN if no solution."""
    f0, f1 = fn(l, b, 1e-4) - v, fn(l, b, dmax) - v
    if np.sign(f0) == np.sign(f1):
        return np.nan
    return brentq(lambda d: fn(l, b, d) - v, 1e-4, dmax, xtol=1e-5)


kin = {'models': list(MODELS), 'expected_vlsr': {}, 'trad_distance_kpc': {}, 'gradient_kms_per_kpc': {}}
for reg in ('P1', 'P2', 'P3', 'P4'):
    l, b = LB[reg]
    kin['expected_vlsr'][reg] = {m: {f'{d:g}kpc': round(float(fn(l, b, d)), 2) for d in (0.0001, 0.16, 0.288, 0.30, 1.0, 1.6, 2.0)}
                                 for m, fn in MODELS.items()}
    kin['gradient_kms_per_kpc'][reg] = {m: round(float((fn(l, b, 0.6) - fn(l, b, 0.1)) / 0.5), 2) for m, fn in MODELS.items()}
    kin['trad_distance_kpc'][reg] = {}
    for line, row in T4[reg].items():
        v = row[2]
        kin['trad_distance_kpc'][reg][line] = {'V_LSR': v, 'paper_D_trad': row[6], 'paper_D_MC': row[3],
                                               **{m: (None if np.isnan(dd := solve_distance(fn, l, b, v)) else round(float(dd), 3))
                                                  for m, fn in MODELS.items()}}
# first-order Oort estimate for the reader
A_flat = 236.0 / (2 * 8.15)
kin['oort_first_order'] = {'A_flat_kms_kpc': A_flat,
                           'dV_dd_P1': A_flat * np.sin(np.deg2rad(2 * LB['P1'][0])) * np.cos(np.deg2rad(LB['P1'][1])) ** 2}


def mc_kd14_emulation(l, b, v, verr, n=20000, sig_pec=0.0, rng=RNG):
    """Emulate kd pdf_kd with reid14: resample (a1,a2,a3,R0) and velocity; grid solution with
    0.1 km/s tolerance (no-solution samples dropped, as kd does); report KDE mode and 68.3% HPD."""
    a1 = rng.normal(241, 8, n); a2 = rng.normal(0.90, 0.06, n)
    a3 = rng.normal(1.46, 0.16, n); R0 = rng.normal(8.34, 0.16, n)
    vs = rng.normal(v, verr, n) + (rng.normal(0, sig_pec, n) if sig_pec > 0 else 0)
    grid = np.arange(0.001, 12.0, 0.004)
    out = np.full(n, np.nan)
    for i0 in range(0, n, 2000):
        sl = slice(i0, i0 + 2000)
        vv = vlsr_kd14(l, b, grid[None, :], a1[sl, None], a2[sl, None], a3[sl, None], R0[sl, None])
        dv = np.abs(vv - vs[sl, None])
        j = np.argmin(dv, axis=1)
        ok = dv[np.arange(len(j)), j] < 0.1
        out[sl] = np.where(ok, grid[j], np.nan)
    return out


def mc_general(fn_params, l, b, v, verr, n=20000, sig_pec=0.0, rng=RNG):
    """MC with linear rotation curves (Reid+2014/2019 A5 parameters and errors) + peculiar velocity.
    No-solution samples (v below V(d=0)) are KEPT and counted (assigned d=0)."""
    R0 = rng.normal(fn_params['R0'], fn_params['R0e'], n)
    T0 = rng.normal(fn_params['T0'], fn_params['T0e'], n)
    sl_ = rng.normal(fn_params['dT'], fn_params['dTe'], n)
    U = rng.normal(fn_params['U'], fn_params['Ue'], n)
    V = rng.normal(fn_params['V'], fn_params['Ve'], n)
    W = rng.normal(fn_params['W'], fn_params['We'], n)
    vs = rng.normal(v, verr, n) + (rng.normal(0, sig_pec, n) if sig_pec > 0 else 0)
    grid = np.arange(0.0, 12.0, 0.004)
    out = np.empty(n)
    for i0 in range(0, n, 2000):
        s = slice(i0, i0 + 2000)
        cl, sl = np.cos(np.deg2rad(l)), np.sin(np.deg2rad(l))
        cb, sb = np.cos(np.deg2rad(b)), np.sin(np.deg2rad(b))
        X = grid[None, :] * cb * cl - R0[s, None]
        Y = grid[None, :] * cb * sl
        Rg = np.maximum(np.hypot(X, Y), 1e-6)
        th = T0[s, None] + sl_[s, None] * (Rg - R0[s, None])
        vhel = (th * Y / Rg - U[s, None]) * cb * cl + (-th * X / Rg - T0[s, None] - V[s, None]) * cb * sl - W[s, None] * sb
        vv = vhel + (U_STD * cb * cl + V_STDc * cb * sl + W_STD * sb)
        # monotonic in d: count samples below V(0)
        idx = np.array([np.searchsorted(row, x) for row, x in zip(vv, vs[s])])
        idx = np.clip(idx, 0, len(grid) - 1)
        out[s] = grid[idx]
    return out


def summarize(samples):
    s = samples[np.isfinite(samples)]
    res = {'n_valid': int(s.size), 'frac_nan': float(1 - s.size / samples.size)}
    if s.size < 10:
        return res
    p16, p50, p84 = np.percentile(s, [16, 50, 84])
    kde = gaussian_kde(s)
    g = np.linspace(0, min(8, s.max()), 1601)
    pdf = kde(g)
    mode = g[np.argmax(pdf)]
    # 68.3% HPD from the KDE
    order = np.argsort(pdf)[::-1]
    cum = np.cumsum(pdf[order]) / pdf.sum()
    sel = g[order[:np.searchsorted(cum, 0.683) + 1]]
    res.update({'median': float(p50), 'p16': float(p16), 'p84': float(p84), 'mode': float(mode),
                'hpd68': [float(sel.min()), float(sel.max())],
                'P_d_lt_0.2kpc': float(np.mean(s < 0.2)), 'P_d_lt_0.35kpc': float(np.mean(s < 0.35)),
                'P_d_gt_1kpc': float(np.mean(s > 1.0))})
    return res


REID14 = dict(R0=8.34, R0e=0.16, T0=240, T0e=8, dT=-0.2, dTe=0.4, U=10.7, Ue=1.8, V=15.6, Ve=6.8, W=8.9, We=0.9)
REID19 = dict(R0=8.15, R0e=0.15, T0=236, T0e=7, dT=-0.1, dTe=0.4, U=10.6, Ue=1.2, V=10.7, Ve=6.0, W=7.6, We=0.7)
mc = {'kd14_emulation_no_pec': {}, 'reid14_linear': {}, 'reid19_linear': {}}
mc_store = {}
t0 = time.time()
for reg in ('P1', 'P2', 'P3', 'P4'):
    l, b = LB[reg]
    for line, row in T4[reg].items():
        key = f'{reg}_{line}'
        s = mc_kd14_emulation(l, b, row[2], row[1], n=8000)
        mc['kd14_emulation_no_pec'][key] = {'paper_MC': [row[3], -row[4], row[5]], **summarize(s)}
        if line == 'OIII' or (reg == 'P4' and line == 'Ha'):
            for name, par in (('reid14_linear', REID14), ('reid19_linear', REID19)):
                mc[name][key] = {}
                for sp in (0.0, 5.0, 7.0, 10.0):
                    ss = mc_general(par, l, b, row[2], row[1], n=20000, sig_pec=sp)
                    sm = summarize(ss)
                    sm['frac_below_V0_assigned_d0'] = float(np.mean(ss == 0.0))
                    mc[name][key][f'sig_pec_{sp:g}'] = sm
                    if name == 'reid14_linear' and reg == 'P1':
                        mc_store[f'{key}_pec{sp:g}'] = ss
mc['runtime_s'] = round(time.time() - t0, 1)
kin['monte_carlo'] = mc
# Velocity precision needed to separate 160 pc from 0: dV/dd * 0.16 kpc
grad = kin['gradient_kms_per_kpc']['P1']['reid14_A5_linear']
kin['verdict_numbers'] = {
    'gradient_P1_reid14_kms_per_kpc': grad,
    'V_rot_at_160pc_P1': kin['expected_vlsr']['P1']['reid14_A5_linear']['0.16kpc'],
    'V_rot_at_300pc_P1': kin['expected_vlsr']['P1']['reid14_A5_linear']['0.3kpc'],
    'V_rot_at_1.6kpc_P1': kin['expected_vlsr']['P1']['reid14_A5_linear']['1.6kpc'],
    'distance_shift_per_7kms_pec_kpc': 7.0 / grad,
    'distance_shift_per_10kms_pec_kpc': 10.0 / grad,
}
R['kinematic_distance'] = kin

# ----------------------------------------------------------------------------------------------
# 3. Physical size and height above the plane
# ----------------------------------------------------------------------------------------------
Z_SUN_PC = 20.8   # Bennett & Bovy 2019
b_c = gal(116.68, 16.42)[1][0]
sizes = {}
for dname, dpc in (('160pc_paper', 160.0), ('288pc_PSR_B0656+14', 288.0), ('300pc_Monogem', 300.0), ('1.6kpc_kd_units', 1600.0)):
    sizes[dname] = {'d_pc': dpc,
                    **{f'len_{a:g}deg_pc': round(float(2 * dpc * np.tan(np.deg2rad(a / 2))), 2) for a in (1, 10, 15, 17)},
                    'z_above_sun_pc': round(float(dpc * np.sin(np.deg2rad(b_c))), 1),
                    'z_above_midplane_pc': round(float(dpc * np.sin(np.deg2rad(b_c)) + Z_SUN_PC), 1),
                    'shell_radius_scaled_from_30pc_at_160pc': round(30.0 * dpc / 160.0, 1)}
sizes['angular_radius_of_30pc_at_160pc_deg'] = float(np.rad2deg(np.arctan(30 / 160)))
# Refined geometry from the geometry teammate (data/derived/silk_trace.json), if present
from silk_geometry import load_trace  # noqa: E402
TRACE = load_trace()
if TRACE and not TRACE.get('preliminary', True):
    sp = TRACE['spine']
    cs = SkyCoord([q['ra'] for q in sp] * u.deg, [q['dec'] for q in sp] * u.deg)
    arc = float(np.sum(cs[:-1].separation(cs[1:]).deg))
    fw_main = 2 * TRACE['half_width_deg']
    fw_env = 2 * TRACE['envelope_half_width_deg']['median']
    sizes['trace'] = {'version': TRACE.get('version'), 'spine_arc_length_deg': round(arc, 2),
                      'spine_dec_range': [min(q['dec'] for q in sp), max(q['dec'] for q in sp)],
                      'main_strand_fwhm_deg': fw_main, 'envelope_full_width_deg': fw_env,
                      'curvature_radius_deg': TRACE.get('curvature', {}).get('radius_deg'),
                      'per_distance': {dn: {'spine_length_pc': round(float(np.deg2rad(arc) * v['d_pc']), 1),
                                            'main_fwhm_pc': round(float(np.deg2rad(fw_main) * v['d_pc']), 2),
                                            'envelope_width_pc': round(float(np.deg2rad(fw_env) * v['d_pc']), 2),
                                            'curvature_radius_pc': (round(float(np.deg2rad(TRACE['curvature']['radius_deg']) * v['d_pc']), 1)
                                                                    if TRACE.get('curvature') else None)}
                                       for dn, v in sizes.items() if isinstance(v, dict) and 'd_pc' in v}}
R['sizes'] = sizes

# ----------------------------------------------------------------------------------------------
# 4. Energetics (Dopita 1979; Cioffi+1988 via Boumis+2022)
# (vs_dopita, e_boumis, r_pds, age_cioffi from paper_review_lib)
# ----------------------------------------------------------------------------------------------

en = {}
vs0 = vs_dopita(64.0, 1.0)
e0 = e_boumis(1.0, 119.0, 30.0)
rp0 = r_pds(e0 / 1e51, 1.0)
en['reproduce'] = {'Vs_dopita_ne64_n01': vs0, 'Vs_dopita_ne150_n01': vs_dopita(150.0, 1.0),
                   'E_erg_n01_vs119_R30': e0, 'R_PDS_pc': rp0,
                   'R_PDS_if_E2.0': r_pds(2.0, 1.0), 'E51_needed_for_RPDS_18pc': (18 / 14) ** 3.5,
                   'age_yr_paper_formula_plus1_RPDS_calc': age_cioffi(e0 / 1e51, 1.0, 30.0, rp0, plus=1.0),
                   'age_yr_correct_plus1over3_RPDS_calc': age_cioffi(e0 / 1e51, 1.0, 30.0, rp0, plus=1 / 3),
                   'age_yr_plus1_RPDS18': age_cioffi(e0 / 1e51, 1.0, 30.0, 18.0, plus=1.0),
                   'age_yr_plus1over3_RPDS18': age_cioffi(e0 / 1e51, 1.0, 30.0, 18.0, plus=1 / 3),
                   'age_snowplow_2R_over_7V_yr': (2 / 7) * 30 * 3.0857e13 / 119 / 3.156e7,
                   'age_0.3R_over_V_yr': 0.3 * 30 * 3.0857e13 / 119 / 3.156e7,
                   'paper': {'Vs': '119+-5', 'E': '(2.0+-0.7)e51', 'R_PDS': '18+-2.8', 'age': '(6.8+-1.8)e4'}}
# Error propagation: (a) paper-like (Vs independent 119+-5, n0 1.0+-0.3); (b) self-consistent Vs(ne, n0)
N = 50000
n0s = np.clip(RNG.normal(1.0, 0.3, N), 0.05, None)
vsa = RNG.normal(119, 5, N)
Ea = e_boumis(n0s, vsa, 30.0)
nes = RNG.uniform(64, 150, N)
vsb = vs_dopita(nes, n0s)
Eb = e_boumis(n0s, vsb, 30.0)
en['propagation'] = {'paper_like_E_p16_50_84': np.percentile(Ea, [16, 50, 84]).tolist(),
                     'selfconsistent_Vs_p16_50_84': np.percentile(vsb, [16, 50, 84]).tolist(),
                     'selfconsistent_E_p16_50_84': np.percentile(Eb, [16, 50, 84]).tolist(),
                     'Vs_frac_err_from_n0_only': 0.5 * 0.3}
alt_d = {}
for dname, dpc in (('160pc', 160.0), ('288pc', 288.0), ('1.6kpc', 1600.0)):
    rs = 30.0 * dpc / 160.0
    rows = {}
    for n0 in (1.0, 0.3, 0.1, 0.03, 0.01):
        vs = vs_dopita(64.0, n0)
        E = e_boumis(n0, vs, rs)
        rp = r_pds(E / 1e51, n0)
        rows[f'n0_{n0:g}'] = {'Vs_dopita_ne64': round(float(vs), 1), 'E_erg': float(E), 'E_over_1e51': round(float(E / 1e51), 2),
                              'R_PDS_pc': round(float(rp), 1), 'Rs_over_RPDS': round(float(rs / rp), 2),
                              'age_yr': float(age_cioffi(E / 1e51, n0, rs, rp)),
                              'E_with_paper_Vs119': float(e_boumis(n0, 119.0, rs))}
    alt_d[dname] = {'Rs_pc': rs, 'rows': rows}
en['alternative_distances'] = alt_d
R['energetics'] = en

# ----------------------------------------------------------------------------------------------
# 5. Line widths
# ----------------------------------------------------------------------------------------------
K_B, M_U = 1.380649e-16, 1.66053907e-24
lw = {'thermal_fwhm_1e4K_kms': {}}
for el, A in (('H', 1.008), ('N', 14.007), ('O', 15.999), ('S', 32.06)):
    lw['thermal_fwhm_1e4K_kms'][el] = round(float(np.sqrt(8 * np.log(2) * K_B * 1e4 / (A * M_U)) / 1e5), 2)
wH = lw['thermal_fwhm_1e4K_kms']['H']
lw['eq2_coefficient_derived'] = 1e4 / (wH ** 2 * (1 - 1.008 / 14.007))
lw['eq3_coefficients_derived'] = {'prefactor': float(np.sqrt(14.007 / (14.007 - 1.008))), 'inner': 1.008 / 14.007}
lw['instrumental_fwhm_kms'] = {f'R{Rr}': round(C_KMS / Rr, 1) for Rr in (7000, 7500, 8000)}
regs = {}
for reg, t in T2.items():
    WH, sH = t['W_H']; WN, sN = t['W_N']
    T = 23.5 * (WH ** 2 - WN ** 2)
    Vt = 1.04 * WN * np.sqrt(1 - 0.071 * WH ** 2 / WN ** 2)
    WHs = RNG.normal(WH, sH, 40000); WNs = RNG.normal(WN, sN, 40000)
    Ts = 23.5 * (WHs ** 2 - WNs ** 2)
    regs[reg] = {'T_eq2_K': round(float(T)), 'T_MC_p16_50_84': np.percentile(Ts, [16, 50, 84]).round().tolist(),
                 'T_sigma_K': round(float(np.std(Ts))), 'P_T_negative': round(float(np.mean(Ts < 0)), 3),
                 'Vt_eq3_kms': round(float(Vt), 2), 'Vt_table4': T4_VT[reg],
                 'Ha_obs_fwhm_R7500': round(float(np.hypot(WH, 40.0)), 1),
                 'NII_obs_fwhm_R7500': round(float(np.hypot(WN, 40.0)), 1)}
    if t['W_O']:
        WO = t['W_O']
        regs[reg]['OIII_obs_fwhm_R7500'] = round(float(np.hypot(WO, 40.0)), 1)
        regs[reg]['OIII_nonthermal_fwhm'] = round(float(np.sqrt(WO ** 2 - lw['thermal_fwhm_1e4K_kms']['O'] ** 2)), 1)
        regs[reg]['Ha_nonthermal_fwhm_1e4K'] = round(float(np.sqrt(max(WH ** 2 - wH ** 2, 0))), 1)
lw['regions'] = regs
# H-alpha fine-structure bias on eq. 2 (two-component approximation, 6.4 km/s separation)
fs = {}
for w in (0.5, 0.667):
    var = w * (1 - w) * 6.4 ** 2
    fs[f'w{w:.3f}'] = {'extra_fwhm2': 8 * np.log(2) * var, 'T_bias_K': 23.5 * 8 * np.log(2) * var}
lw['Ha_fine_structure_T_bias'] = fs


# ---- Monte-Carlo of Gaussian line fitting at the paper's S/N, with the paper's FWHM cut ----------
PIX = 6.0  # km/s per pixel (LAMOST MRS blue arm ~0.1 A/pix at 5007 A)
VGRID = np.arange(-240, 240 + PIX, PIX)


def gmodel(v, a, v0, sig, c0, c1):
    return c0 + c1 * v / 100 + a * np.exp(-0.5 * ((v - v0) / sig) ** 2)


def fisher_flux_sn(sig_obs, amp=1.0, noise=1.0):
    p = np.array([amp, 0.0, sig_obs, 0.0, 0.0])
    J = np.empty((VGRID.size, 5)); h = 1e-4
    for i in range(5):
        dp = np.zeros(5); dp[i] = h * max(1, abs(p[i]))
        J[:, i] = (gmodel(VGRID, *(p + dp)) - gmodel(VGRID, *(p - dp))) / (2 * dp[i])
    cov = np.linalg.inv(J.T @ J / noise ** 2)
    F = amp * sig_obs * np.sqrt(2 * np.pi)
    g = np.array([sig_obs * np.sqrt(2 * np.pi), 0, amp * np.sqrt(2 * np.pi), 0, 0])
    return F / np.sqrt(g @ cov @ g)


def sim_fits(w_int, sn, n=1500, w_inst=40.0, w_inst_err=1.5, rng=RNG, v_true=0.0):
    sig_obs = np.hypot(w_int, w_inst) / 2.3548
    noise = fisher_flux_sn(sig_obs) / sn   # noise per pixel giving the requested flux S/N
    res = []
    for _ in range(n):
        y = gmodel(VGRID, 1.0, v_true, sig_obs, 0, 0) + rng.normal(0, noise, VGRID.size)
        try:
            p, cov = curve_fit(gmodel, VGRID, y, p0=[0.8, v_true + rng.normal(0, 5), sig_obs * rng.uniform(0.8, 1.25), 0, 0],
                               sigma=np.full(VGRID.size, noise), absolute_sigma=True, maxfev=4000,
                               bounds=([0, -120, 2.0, -np.inf, -np.inf], [np.inf, 120, 150, np.inf, np.inf]))
            e = np.sqrt(np.diag(cov))
        except Exception:
            continue
        w_obs, ew = 2.3548 * p[2], 2.3548 * e[2]
        winst_meas = rng.normal(w_inst, w_inst_err)
        passed = (w_obs - ew) > (winst_meas + w_inst_err)
        w_rec = np.sqrt(max(w_obs ** 2 - winst_meas ** 2, 0.0))
        res.append((w_obs, ew, passed, w_rec, p[1], e[1], p[0] * p[2] * np.sqrt(2 * np.pi)))
    return np.array(res)


fitmc = {}
fit_store = {}
for sn in (6.1, 15.0, 50.0):
    for w_int in (5.4, 20.0, 30.0, 45.0):
        a = sim_fits(w_int, sn, n=1200 if sn < 10 else 600)
        pas = a[:, 2].astype(bool)
        key = f'SN{sn:g}_Wint{w_int:g}'
        fitmc[key] = {'n': int(len(a)), 'frac_pass_cut': float(pas.mean()),
                      'Wrec_all_p16_50_84': np.percentile(a[:, 3], [16, 50, 84]).round(1).tolist(),
                      'Wrec_pass_p16_50_84': (np.percentile(a[pas, 3], [16, 50, 84]).round(1).tolist() if pas.sum() > 5 else None),
                      'v_centroid_std_kms': round(float(np.std(a[:, 4])), 2),
                      'v_centroid_median_formal_err': round(float(np.median(a[:, 5])), 2)}
        if sn == 6.1:
            fit_store[key] = a
lw['fit_mc'] = fitmc
lw['fit_mc_setup'] = {'pixel_kms': PIX, 'window_kms': [-240, 240], 'W_inst': 40.0, 'W_inst_err': 1.5,
                      'SN_definition': 'integrated-flux S/N of a free Gaussian+linear fit (Fisher-matrix calibrated)',
                      'cut': 'paper: FWHM_obs - err_obs > FWHM_inst + err_inst'}


# [O III] 5007/4959 at the paper's S/N: fit both lines independently
def sim_ratio(sn5007=6.1, w_int=45.0, n=1500, rng=RNG):
    sig_obs = np.hypot(w_int, 40.0) / 2.3548
    noise = fisher_flux_sn(sig_obs) / sn5007
    rat, sn49 = [], []
    for _ in range(n):
        out = []
        for amp in (1.0, 1.0 / 2.98):
            y = gmodel(VGRID, amp, 0, sig_obs, 0, 0) + rng.normal(0, noise, VGRID.size)
            try:
                p, cov = curve_fit(gmodel, VGRID, y, p0=[amp, rng.normal(0, 5), sig_obs, 0, 0],
                                   sigma=np.full(VGRID.size, noise), absolute_sigma=True, maxfev=4000,
                                   bounds=([0, -60, 2.0, -np.inf, -np.inf], [np.inf, 60, 150, np.inf, np.inf]))
                F = p[0] * p[2] * np.sqrt(2 * np.pi)
                g = np.array([p[2] * np.sqrt(2 * np.pi), 0, p[0] * np.sqrt(2 * np.pi), 0, 0])
                out.append((F, F / np.sqrt(g @ cov @ g)))
            except Exception:
                out.append((np.nan, np.nan))
        if out[1][0] > 0:
            rat.append(out[0][0] / out[1][0]); sn49.append(out[1][1])
    rat, sn49 = np.array(rat), np.array(sn49)
    return rat, sn49


rat, sn49 = sim_ratio()
oiii = {'theoretical_5007_4959': 2.98, 'regions': {}}
for reg in ('P1', 'P2', 'P3'):
    t = T2[reg]
    r_ = t['O5007'] / t['O4959']
    sr = r_ * np.sqrt(1 / t['SN_O5007'] ** 2 + 1 / t['SN_O4959'] ** 2)
    oiii['regions'][reg] = {'ratio': round(r_, 3), 'sigma': round(float(sr), 2), 'z_vs_2.98': round(float((r_ - 2.98) / sr), 2),
                            'expected_SN4959_from_SN5007': round(t['SN_O5007'] / 2.98, 2), 'SN4959_paper': t['SN_O4959']}
fin = np.isfinite(rat)
oiii['mc_independent_fits'] = {'n': int(fin.sum()), 'ratio_p16_50_84': np.percentile(rat[fin], [16, 50, 84]).round(2).tolist(),
                               'frac_ratio_below_2.5': float(np.mean(rat[fin] < 2.5)),
                               'SN4959_median': float(np.nanmedian(sn49)),
                               'frac_SN4959_below_3': float(np.nanmean(sn49 < 3))}
R['line_widths'] = lw
R['oiii_doublet'] = oiii

# ----------------------------------------------------------------------------------------------
# 6. [S II] density and physical self-consistency (EM, path length, pressure)
# (ne_proxauf = paper eq. 5; ne_sanders / R_sanders = Sanders+2016; from paper_review_lib)
# ----------------------------------------------------------------------------------------------
dens ={'formula_check': {}, 'regions': {}, 'inverse': {}}
for Rr in (1.03, 1.10, 1.14, 1.19, 1.30, 1.40):
    dens['formula_check'][f'R{Rr:.2f}'] = {'ne_proxauf_eq5': round(float(ne_proxauf(Rr)), 1),
                                            'ne_sanders2016': round(float(ne_sanders(Rr)), 1)}
for ne in (64, 100, 150):
    dens['inverse'][f'ne{ne}'] = {'R_sanders2016': round(float(R_sanders(ne)), 3)}
dens['low_density_limit_sanders'] = SA * SB_ / SC
for reg, t in T2.items():
    s1 = RNG.normal(*t['S6716'], 50000); s2 = RNG.normal(*t['S6731'], 50000)
    ok = (s1 > 0) & (s2 > 0)
    Rs = s1[ok] / s2[ok]
    Rtab = t['ratio_SII']
    Rint = t['S6716'][0] / t['S6731'][0]
    nes = ne_sanders(np.clip(Rs, 0.45, 1.449))
    dens['regions'][reg] = {'R_table': Rtab, 'R_from_intensities': round(Rint, 3),
                            'R_MC_p16_50_84': np.percentile(Rs, [16, 50, 84]).round(3).tolist(),
                            'ne_proxauf_at_Rtable': round(float(ne_proxauf(Rtab)), 0),
                            'ne_sanders_at_Rtable': round(float(ne_sanders(Rtab)), 0),
                            'ne_sanders_at_Rtable_T7500': round(float(ne_sanders(Rtab) * np.sqrt(0.75)), 0),
                            'ne_sanders_MC_p16_50_84': np.percentile(nes, [16, 50, 84]).round(0).tolist(),
                            'frac_MC_R_above_low_density_limit': float(np.mean(Rs > 1.449))}

# WHAM-SS DR1 H-alpha at P1-P4 (nearest 1-deg beam), intensity and centroid velocity
wham = {}
try:
    from astropy.io import fits
    W = fits.open(ROOT / 'data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits')[1].data
    wl, wb = W['GAL-LON'].astype(float), W['GAL-LAT'].astype(float)
    cw = SkyCoord(wl * u.deg, wb * u.deg, frame='galactic')
    for reg, (l, b) in LB.items():
        sep = cw.separation(SkyCoord(l * u.deg, b * u.deg, frame='galactic')).deg
        idx = np.argsort(sep)[:4]
        rows = []
        for i in idx:
            v, f, var = W['VELOCITY'][i].astype(float), W['DATA'][i].astype(float), W['VARIANCE'][i].astype(float)
            m = (v > -60) & (v < 70) & np.isfinite(f)
            try:
                p, cov = curve_fit(lambda x, a, x0, s, c: c + a * np.exp(-0.5 * ((x - x0) / s) ** 2), v[m], f[m],
                                   p0=[f[m].max(), 5, 12, 0], sigma=np.sqrt(np.clip(var[m], 1e-6, None)), maxfev=5000)
                e = np.sqrt(np.diag(cov))
                vc, evc, fw = float(p[1]), float(e[1]), float(2.3548 * abs(p[2]))
            except Exception:
                vc = evc = fw = np.nan
            mm = (v > -40) & (v < 50)
            vmean = float(np.sum(f[mm] * v[mm]) / np.sum(f[mm]))
            rows.append({'l': float(wl[i]), 'b': float(wb[i]), 'sep_deg': round(float(sep[i]), 3),
                         'I_Ha_R': round(float(W['INTEN'][i]), 3), 'err_R': round(float(W['ERROR'][i]), 3),
                         'star_flag': bool(W['STAR'][i]), 'V_LSR_gauss': round(vc, 2), 'V_LSR_gauss_err': round(evc, 2),
                         'FWHM_obs_gauss': round(fw, 1), 'V_LSR_mean_-40_50': round(vmean, 2)})
        wham[reg] = rows
except Exception as ex:  # pragma: no cover
    wham['error'] = repr(ex)
dens['wham_ss_dr1'] = wham

fink = {}
try:
    import healpy as hp
    from astropy.io import fits
    fk = fits.open(ROOT / 'data/raw/lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits')[1].data.field(0).ravel().astype(float)
    for reg, (l, b) in LB.items():
        vals = []
        for dl in np.linspace(-0.2, 0.2, 5):
            for db in np.linspace(-0.2, 0.2, 5):
                vals.append(fk[hp.ang2pix(1024, l + dl / np.cos(np.deg2rad(b)), b + db, lonlat=True, nest=True)])
        fink[reg] = {'I_Ha_R_median_0.4deg_box': round(float(np.median(vals)), 3), 'min': round(float(np.min(vals)), 3),
                     'max': round(float(np.max(vals)), 3)}
except Exception as ex:  # pragma: no cover
    fink['error'] = repr(ex)
dens['finkbeiner2003'] = fink

# [O III] emissivity (collisional excitation 3P -> 1D2, low-density limit): q_oiii, em_from_oiii, em_from_ha (lib)
selfc ={'oiii_emissivity': {f'T{T:g}': {'q_cm3_s': float(q_oiii(T)), 'EM_min_for_10R': float(em_from_oiii(10, T)),
                                         'EM_min_for_13R': float(em_from_oiii(13, T)),
                                         'Ha_R_implied_by_EM_min_10R': float(em_from_oiii(10, T) / (2.75 * (T / 1e4) ** 0.9))}
                             for T in (8e3, 1e4, 1.5e4, 2e4, 3e4, 5e4)},
         'atomic_data': {'Omega_3P_1D_OIII': OMEGA_1D, 'g_3P': G_3P, 'T_exc_K': T_EXC, 'BR_5007': BR5007, 'O/H': OH}}
ha_wham = {reg: rows[0]['I_Ha_R'] for reg, rows in wham.items() if isinstance(rows, list)}
EMs = {}
for reg in ('P1', 'P2', 'P3', 'P4'):
    if reg in ha_wham:
        EMs[reg] = {'EM_Ha_WHAM_T1e4': round(em_from_ha(ha_wham[reg], 1e4), 2)}
selfc['EM_from_WHAM_Ha'] = EMs
path = {}
for dname, dpc in (('160pc', 160.0), ('300pc', 300.0), ('1.6kpc', 1600.0)):
    width = 2 * dpc * np.tan(np.deg2rad(0.5))
    em_ref = em_from_oiii(11.5, 1e4)
    w_main = (np.deg2rad(sizes['trace']['main_strand_fwhm_deg']) * dpc) if 'trace' in sizes else np.nan
    path[dname] = {'filament_width_1deg_pc': round(float(width), 2), 'main_strand_fwhm_pc': round(float(w_main), 2),
                   'ne_if_L_equals_main_fwhm': round(float(np.sqrt(em_ref / w_main)), 2),
                   **{f'L_pc_ne{ne}': float(em_ref / ne ** 2) for ne in (64, 150, 300)},
                   **{f'filling_factor_ne{ne}': float(em_ref / ne ** 2 / width) for ne in (64, 150, 300)},
                   'ne_if_L_equals_width': round(float(np.sqrt(em_ref / width)), 2),
                   'Ha_R_if_ne64_fills_width': float(64 ** 2 * width / 2.75)}
selfc['path_length_vs_width_EM_OIII_11.5R_T1e4'] = path
selfc['pressure'] = {'P_over_k_ne64_T1e4': 2 * 64 * 1e4, 'P_over_k_ne150_T1e4': 2 * 150 * 1e4,
                     'ram_P_over_k_n01_vs119': 1.4 * 1.6726e-24 * 1.0 * (119e5) ** 2 / 1.380649e-16}
# --- H-alpha excess on the filament and the implied [O III]/H-alpha and temperature ---
# Calibrate LAMOST relative H-alpha (sky-line units) against WHAM total H-alpha at P1-P4.
if all(isinstance(wham.get(r), list) for r in ('P1', 'P2', 'P3', 'P4')):
    cal = np.array([wham[r][0]['I_Ha_R'] / T2[r]['Ha'][0] for r in ('P1', 'P2', 'P3', 'P4')])
    k_cal, k_err = float(cal.mean()), float(cal.std(ddof=1))
    excess = {}
    for reg in ('P1', 'P2', 'P3'):
        d = T2[reg]['Ha'][0] - T2['P4']['Ha'][0]
        e = float(np.hypot(T2[reg]['Ha'][1], T2['P4']['Ha'][1]))
        exc, exc_e = d * k_cal, e * k_cal
        ul2 = max(exc, 0) + 2 * exc_e
        excess[reg] = {'Ha_excess_R': round(exc, 2), 'Ha_excess_err_R': round(exc_e, 2), 'Ha_excess_2sigma_UL_R': round(ul2, 2),
                       'WHAM_minus_P4_R': round(wham[reg][0]['I_Ha_R'] - wham['P4'][0]['I_Ha_R'], 2)}
    selfc['Ha_excess_on_filament'] = {'R_per_LAMOST_unit_mean': k_cal, 'R_per_LAMOST_unit_std': k_err,
                                      'per_unit_values': cal.round(2).tolist(), 'regions': excess}

Tg =np.logspace(np.log10(6e3), np.log10(1e5), 400)
selfc['oiii_ha_ratio_vs_T_x1'] = {f'T{T:g}': round(float(oiii_ha_photon_ratio(T)), 2) for T in (8e3, 1e4, 1.5e4, 2e4, 3e4, 5e4)}
R['sii_density'] = dens
R['self_consistency'] = selfc

# ----------------------------------------------------------------------------------------------
# 7 (see oiii_doublet above).  8. Metallicity calibrations
# ----------------------------------------------------------------------------------------------
met = {}
for reg, t in T2.items():
    n6583, ha, s = t['N6583'][0], t['Ha'][0], t['S6716'][0] + t['S6731'][0]
    nsum = n6583 + t['N6548'][0]
    N2 = np.log10(n6583 / ha)
    N2sum = np.log10(nsum / ha)
    y = np.log10(n6583 / s) + 0.264 * np.log10(n6583 / ha)
    met[reg] = {'N2_6583': round(float(N2), 3), 'N2_if_sum': round(float(N2sum), 3),
                'PP04_linear_6583': round(float(8.90 + 0.57 * N2), 3),
                'PP04_linear_sum': round(float(8.90 + 0.57 * N2sum), 3),
                'PP04_cubic_6583': round(float(9.37 + 2.03 * N2 + 1.26 * N2 ** 2 + 0.32 * N2 ** 3), 3),
                'Marino13_N2_6583': round(float(8.743 + 0.462 * N2), 3),
                'D16_y': round(float(y), 3), 'D16_OH': round(float(8.77 + y + 0.45 * (y + 0.3) ** 5), 3),
                'N2_in_PP04_valid_range(-2.5,-0.3)': bool(-2.5 < N2 < -0.3),
                'N2_in_Marino13_range(-1.6,-0.2)': bool(-1.6 < N2 < -0.2)}
R['metallicity'] = {'regions': met, 'paper': {'N2': 8.82, 'N2S2Ha': 8.70, 'solar_quoted': 8.66},
                    'notes': 'PP04 N2 linear: 12+log(O/H)=8.90+0.57 N2 (sigma 0.18 dex, -2.5<N2<-0.3); cubic: 9.37+2.03N2+1.26N2^2+0.32N2^3'}

# ----------------------------------------------------------------------------------------------
# 9. Classification logic: line ratios on vs off nebula; WHAM comparison
# ----------------------------------------------------------------------------------------------
cls = {'regions': {}}
for reg, t in T2.items():
    n = 50000
    s1, s2 = RNG.normal(*t['S6716'], n), RNG.normal(*t['S6731'], n)
    n1, n2 = RNG.normal(*t['N6548'], n), RNG.normal(*t['N6583'], n)
    h = RNG.normal(*t['Ha'], n)
    r_sii, r_nsum, r_n83, r_s16 = (s1 + s2) / h, (n1 + n2) / h, n2 / h, s1 / h
    cls['regions'][reg] = {k: [round(float(np.median(v)), 3), round(float(np.std(v)), 3)] for k, v in
                           (('SII_sum_Ha', r_sii), ('NII_sum_Ha', r_nsum), ('NII6583_Ha', r_n83), ('SII6716_Ha', r_s16),
                            ('NII6583_over_6548', n2 / n1))}
    cls['regions'][reg]['table2_SII_Ha'] = t['ratio_SII_Ha']
    cls['regions'][reg]['table2_NII_Ha'] = t['ratio_NII_Ha']
    cls['regions'][reg]['Ha_minus_P4_skyunits'] = [round(t['Ha'][0] - T2['P4']['Ha'][0], 3), round(float(np.hypot(t['Ha'][1], T2['P4']['Ha'][1])), 3)]
for reg in ('P1', 'P2', 'P3'):
    a, b_ = cls['regions'][reg], cls['regions']['P4']
    cls['regions'][reg]['delta_vs_P4_sigma'] = {k: round((a[k][0] - b_[k][0]) / np.hypot(a[k][1], b_[k][1]), 2)
                                                for k in ('SII_sum_Ha', 'NII_sum_Ha', 'NII6583_Ha', 'SII6716_Ha')}
cls['wham_literature'] = {
    'Madsen2006_WIM_NII6583_Ha': '~0.5 typical, in some cases >1.0 (Madsen, Reynolds & Haffner 2006, ApJ 652, 401, Sect. 4)',
    'Madsen2006_WIM_SII6716_Ha_faint': '>~0.4 in regions of faint Halpha',
    'Madsen2006_WIM_OIII_Ha': '~0.06 (Reynolds 1985), confirmed',
    'Madsen2006_Table5_high_latitude_filaments_l217-236': {
        'A#1 (226.0,+10.2)': [0.49, 0.32, 0.10], 'A#2': [0.44, 0.29, 0.14], 'B#1 (225.2,+22.1)': [0.63, 0.33, 0.06],
        'B#2': [0.58, 0.26, 0.04], 'C (217.1,+28.9)': [0.27, 0.26, 0.01], 'D (235.2,+29.7)': [0.46, 0.29, 0.004],
        'E#1 (226.2,+38.2)': [0.50, 0.40, 0.008], 'E#2': [0.40, 0.36, 0.01], 'F (236.3,+50.9)': [0.46, 0.23, 0.02],
        'Arc (170.8,+56.9)': [0.72, 0.46, 0.06], 'columns': '[NII]6583/Ha, [SII]6716/Ha, [OIII]5007/Ha (energy; upper limits given as value)'},
    'Haffner1999': '[SII]/Ha and [NII]/Ha increase as Ha intensity decreases; [SII]/[NII] ~constant; T 6000-9000 K (ApJ 523, 223)',
}
R['classification'] = cls

# ----------------------------------------------------------------------------------------------
# 10. Surface brightness -> Rayleigh
# (mag_to_R from paper_review_lib)
# ----------------------------------------------------------------------------------------------
sb ={'one_R_at_5007_erg_s_cm2_arcsec2': mag_to_R(23.6)[1], 'regions': {}}
for reg, (m, em) in T3_SB.items():
    sb['regions'][reg] = {'mag': m, 'R_AB_60A': round(float(mag_to_R(m)[0]), 2),
                          'R_AB_60A_vega_minus0.1': round(float(mag_to_R(m - 0.1)[0]), 2),
                          'R_AB_60A_vega_plus0.1': round(float(mag_to_R(m + 0.1)[0]), 2),
                          'R_AB_effwidth_50A': round(float(mag_to_R(m, width=50.0)[0]), 2),
                          'R_AB_effwidth_70A': round(float(mag_to_R(m, width=70.0)[0]), 2)}
# [O III]/H-alpha lower limits on the filament and the minimum temperature they require (x(O++)=1, solar O/H)
if 'Ha_excess_on_filament' in R['self_consistency']:
    lim = {}
    for reg, ex in R['self_consistency']['Ha_excess_on_filament']['regions'].items():
        I_o3 = sb['regions'][reg]['R_AB_60A']
        rat_lo = I_o3 / ex['Ha_excess_2sigma_UL_R']
        ok = oiii_ha_photon_ratio(Tg) >= rat_lo
        lim[reg] = {'I_OIII_R': I_o3, 'OIII_over_Ha_photon_2sigma_LL': round(float(rat_lo), 1),
                    'T_min_K_for_x1': (round(float(Tg[ok][0])) if ok.any() else None)}
    R['self_consistency']['oiii_ha_limits'] = lim
sb['note'] = ('Assumes the tabulated magnitude is an AB-like flux density of pure line emission spread over the '
              '6.0 nm FWHM; the paper does not say whether Table 3 is background-subtracted or Vega/AB.')
R['surface_brightness'] = sb

# ----------------------------------------------------------------------------------------------
# 11. Internal inconsistencies (text locations from the arXiv v1 PDF)
# ----------------------------------------------------------------------------------------------
R['internal_inconsistencies'] = [
    {'id': 'units_distance', 'where': 'Table 4 footnotes 1-2; Sect. 3.3; abstract',
     'issue': "Distance columns labelled 'in units of 100 pc' are kd outputs in kpc (kd documents 'distance (kpc)'); "
              "fiducial 160 pc should read ~1.6 kpc; 28 pc -> ~280 pc; red/blue distance gap '100-150 pc' -> 1-1.5 kpc."},
    {'id': 'vlsr_offsets', 'where': 'Table 4', 'issue': 'V_LSR - V_helio differs between lines of the same region (see lsr.summary)'},
    {'id': 'extent', 'where': 'Fig.1 caption/Sect.1/Sect.3 vs press', 'issue': '1x10 deg in paper vs 1.5x15 deg (press) vs ~17 deg arc (NSNS)'},
    {'id': 'expansion_speed', 'where': 'Sect. 4.2', 'issue': "'expansion speed ~140 km/s, Section 3.4' - no such number in Sect. 3.4; shock speed given is 119 km/s (Sect. 4.1)"},
    {'id': 'temperature', 'where': 'Sect. 3.4 text vs Fig. 7a vs Sect. 3.4 end', 'issue': "median 'about 9000 K' (text) vs 8727 K (Fig 7a) vs 'approximately 10,000 K' (adopted); Fig 7a legend mu=6471, sigma=6798 but caption says sigma~6.47e3"},
    {'id': 'distance_vs_monogem', 'where': 'Sect. 4.1', 'issue': 'adopts 160 pc yet cites Fesen+2026 placing it in the Monogem Ring (~300 pc) without discussion'},
    {'id': 'P1_photometry', 'where': 'Sect. 3.1', 'issue': "P1 'appears the faintest' but has the highest SB (23.60); attributed to stellar contamination -> SBs not reliable to ~0.3 mag"},
    {'id': 'SN_threshold', 'where': 'Sect. 3.2.1 vs Table 2', 'issue': '[OIII]4959 S/N 1.8-2.3 is below the paper\'s own S/N>3 usability threshold yet used as a detection'},
    {'id': 'sii_threshold', 'where': 'abstract/Sect.1/Sect.4.1/Sect.5', 'issue': 'shock threshold quoted as 0.4, 0.45 and 0.45-0.5 in different places'},
    {'id': 'N2_polynomial', 'where': 'Sect. 3.5 eq. 6', 'issue': "text says 'third-order polynomial' but eq. 6 is the linear PP04 N2 relation"},
    {'id': 'solar_OH', 'where': 'Sect. 3.5', 'issue': "solar 12+log(O/H)=8.66 attributed to Allende Prieto+2001 (who give 8.69+-0.05; 8.66 is Asplund+2004)"},
    {'id': 'density_vs_ratio', 'where': 'Sect. 3.4 vs Table 2', 'issue': 'n_e 64-150 cm^-3 requires R=1.25-1.36, but Table 2 lists R=1.03-1.19 (n_e~190-440); background P4 has the LOWEST R (highest n_e), contradicting "background regions show lower electron densities"'},
    {'id': 'Vt_table', 'where': 'Table 4 vs eq. 3 + Table 2', 'issue': 'eq. 3 with Table 2 widths gives V_t 25.8-29.9 km/s, Table 4 lists 29.9-33.1'},
    {'id': 'age_formula', 'where': 'Sect. 4.1', 'issue': "printed age formula has '+1' where Cioffi+1988 inversion gives '+1/3'"},
    {'id': 'fig5_marker', 'where': 'Fig. 5 right', 'issue': 'dashed [OIII] marker at ~+9 km/s, but histogram centroid ~+14 and Table 4 gives +11 to +16'},
    {'id': 'fig4_axis', 'where': 'Fig. 4/5', 'issue': "x-axis labelled 'RA (degrees)' with sexagesimal hour ticks (Fig. 4)"},
    {'id': 'fig2_caption', 'where': 'Fig. 2', 'issue': "caption mentions a 'colored [OIII] overlay' but the figure is greyscale"},
    {'id': 'radius', 'where': 'Sect. 3.3/4.1', 'issue': 'shell radius ~30 pc (10.6 deg at 160 pc) asserted from curvature without a fit; N-S extent only 28 pc'},
    {'id': 'radio', 'where': 'Sect. 4.1', 'issue': 'NVSS/TGSS are interferometric and resolve out degree-scale emission (NVSS LAS ~16 arcmin): non-detection is uninformative'},
]

# ----------------------------------------------------------------------------------------------
R['meta'] = {'script': 'scripts/paper_checks.py', 'run_utc': Time.now().isot, 'access_date': '2026-11-16',
             'python': sys.version.split()[0], 'numpy': np.__version__}


def to_jsonable(o):
    if isinstance(o, dict):
        return {str(k): to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return to_jsonable(o.tolist())
    return o


OUT.write_text(json.dumps(to_jsonable(R), indent=1))
np.savez_compressed(OUT_MC, **{f'kin_{k}': v for k, v in mc_store.items()},
                    **{f'fit_{k}': v for k, v in fit_store.items()}, ratio_mc=rat, sn4959_mc=sn49)
print('wrote', OUT, 'and', OUT_MC)
