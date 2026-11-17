"""verify_velocity: can imperfect geocoronal-Halpha subtraction (and solar Halpha absorption in the sky
continuum) shift the LAMOST MRS-N Halpha centroid redward by 4-7 km/s relative to [N II]/[S II] and WHAM?

Forward model at R ~ 7500 (LSF sigma from the OH 6554 line dispersion 0.32 A quoted by Zhang+2021, RAA 21, 280,
arXiv:2108.08021; R=7500 variant):
  WIM Halpha     : V_helio = +12.2 km/s (V_LSR = 0, WHAM), intrinsic FWHM 27 km/s, intensity 1 (unit)
  geocoronal     : topocentric rest, i.e. V_helio = v_corr = -8.32 km/s on 2024-02-01 (astropy, Xinglong),
                   plus the fine-structure offset of Ly-beta-pumped 3p->2s emission relative to the
                   recombination-weighted Halpha rest (6564.61 A vac): 2:1 components at 6564.523/6564.569 A
                   (NIST H I levels) -> -3.3 km/s; intensity g (relative to WIM).
  MRS-N model    : Gaussian with the RV of OH 6554 (Zhang+2021 Sect. 2.4; Wu+2022 Sect. 3.3), i.e. at v_corr + dm,
                   dm = 0 if the RV is transferred with the recombination rest wavelength, dm = -8.24 km/s if the
                   centroid is literally lam(6554)*6563/6554 as printed; intensity k*g; sigma = sqrt(s6554^2+0.02^2).
  solar Halpha   : optional continuum C_sun (per A, relative to the WIM line flux) with the Kurucz 2005 solar
                   atlas at topocentric rest (zodiacal light / starlight), LSF-convolved.
Fit: single Gaussian + 2nd-order polynomial over +-8 A (as the MRS-N pipeline: Gaussian + 2nd-order polynomial).

Paper Table 2/4 test data: Halpha - [N II] V_helio = +4.5 (P1), +5.0 (P2), +6.7 (P3), +6.4 (P4); relative
Halpha = 0.30, 0.37, 0.26, 0.25 (units of OH 6553); model geocoronal = ratio(sunalt) = 0.954-0.011*|sunalt|
with sunalt -51..-62 deg (exposure mid-times) -> 0.33 in the same units.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_halpha_sim.py
Output: data/derived/verify_velocity_halpha_sim.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/verify_velocity_halpha_sim.json'
C = 299792.458
L0 = 6564.61            # recombination-weighted Halpha, vacuum
VCORR = -8.32
V_NEB = 12.2
# H I fine structure (NIST energy levels, cm^-1): 2s1/2 82258.9544, 3p1/2 97492.2013, 3p3/2 97492.3097
L_3P32 = 1e8 / (97492.3097 - 82258.9544)
L_3P12 = 1e8 / (97492.2013 - 82258.9544)
DV_FS = C * ((2 * L_3P32 + L_3P12) / 3 / L0 - 1)


def gauss(v, mu, s):
    return np.exp(-0.5 * ((v - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def solar_halpha(v):
    d = np.loadtxt(ROOT / 'data/raw/verify_velocity_solar/kurucz2005_flux_650-662nm.asc')
    lam_air = d[:, 0] * 10
    # air -> vacuum (Morton) for consistency with L0
    s2 = (1e4 / lam_air) ** 2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    lv = lam_air * n
    vv = C * (lv / L0 - 1)            # solar lines at topocentric rest -> V_helio = vv + VCORR
    return np.interp(v, vv + VCORR, d[:, 1])


def simulate(g=1.0, k=1.0, dm=0.0, sig_inst=14.6, fwhm_neb=27.0, window_A=8.0, dv_fs=None):
    """g = true geocoronal / WIM; k = subtracted model / true geocoronal; returns (centroid bias, sigma)."""
    dv_fs = DV_FS if dv_fs is None else dv_fs
    dvpix = 0.1 / L0 * C
    v = np.arange(-window_A / L0 * C, window_A / L0 * C, dvpix)
    s_neb = np.hypot(fwhm_neb / 2.3548, sig_inst)
    s_geo = np.hypot(2.9, sig_inst)
    y = gauss(v, V_NEB, s_neb)
    # geocoronal true: 2:1 components separated by 2.13 km/s around VCORR+dv_fs
    sep = C * (L_3P12 - L_3P32) / L0
    y += g * (2 / 3 * gauss(v, VCORR + dv_fs - sep / 3, s_geo) + 1 / 3 * gauss(v, VCORR + dv_fs + 2 * sep / 3, s_geo))
    s_mod = np.hypot(sig_inst, 0.02 / L0 * C)
    y -= k * g * gauss(v, VCORR + dm, s_mod)
    def model(vv, A, mu, s, c0, c1, c2):
        return A * gauss(vv, mu, s) + c0 + c1 * vv / 100 + c2 * (vv / 100) ** 2
    p, _ = curve_fit(model, v, y, p0=[1.0, V_NEB, s_neb, 0, 0, 0], maxfev=20000)
    return float(p[1] - V_NEB), float(abs(p[2]))


def simulate_with_solar(g, k, dm, csun_per_A, sig_inst=14.6, fwhm_neb=27.0, window_A=8.0):
    """Same as simulate() but with an explicit solar continuum: continuum level csun_per_A in units of
    (WIM line flux) per Angstrom."""
    dvpix = 0.1 / L0 * C
    v = np.arange(-window_A / L0 * C, window_A / L0 * C, dvpix)
    s_neb = np.hypot(fwhm_neb / 2.3548, sig_inst)
    s_geo = np.hypot(2.9, sig_inst)
    sep = C * (L_3P12 - L_3P32) / L0
    line = gauss(v, V_NEB, s_neb)
    geo = g * (2 / 3 * gauss(v, VCORR + DV_FS - sep / 3, s_geo) + 1 / 3 * gauss(v, VCORR + DV_FS + 2 * sep / 3, s_geo))
    mod = k * g * gauss(v, VCORR + dm, np.hypot(sig_inst, 0.02 / L0 * C))
    sol = gaussian_filter1d(solar_halpha(v), sig_inst / dvpix)
    cont = csun_per_A * (L0 / C) * sol        # per km/s
    y = line + geo - mod + cont

    def model(vv, A, mu, s, c0, c1, c2):
        return A * gauss(vv, mu, s) + c0 + c1 * vv / 100 + c2 * (vv / 100) ** 2
    p, _ = curve_fit(model, v, y, p0=[1.0, V_NEB, s_neb, cont.mean(), 0, 0], maxfev=40000)
    return float(p[1] - V_NEB), float(abs(p[2]) * 2.3548)


def main():
    out = dict(DV_FS=DV_FS, VCORR=VCORR, V_NEB=V_NEB)
    # grid
    grid = []
    for sig in (14.6, 17.0):
        for dm in (0.0, -8.24):
            for g in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
                for k in (0.6, 0.8, 1.0, 1.2, 1.4):
                    b, s = simulate(g=g, k=k, dm=dm, sig_inst=sig)
                    grid.append(dict(sig_inst=sig, dm=dm, g=g, k=k, bias=round(b, 2)))
    out['grid'] = grid
    # paper regions
    paper = {'P1': (16.3 - 11.8, 0.30), 'P2': (17.0 - 12.0, 0.37), 'P3': (18.4 - 11.7, 0.26), 'P4': (19.7 - 13.3, 0.25)}
    gmod = 0.954 - 0.011 * 56.4
    pred = {}
    for name, (obs, ineb) in paper.items():
        g = gmod / ineb
        row = dict(obs_Ha_minus_NII=obs, I_neb_rel=ineb, g_model=round(g, 2))
        for dm in (0.0, -8.24):
            for kt in (0.7, 1.0, 1.3):      # kt = true / model geocoronal intensity
                row[f'pred_dm{dm}_true/model{kt}'] = round(simulate(g=g * kt, k=1 / kt, dm=dm)[0], 2)
        pred[name] = row
    out['paper_regions'] = pred
    # best (dm, k) by chi2 against observed offsets (sigma ~ 1.9 km/s from Table 4 errors in quadrature)
    best = []
    for dm in np.arange(-12, 4.01, 0.5):
        for kt in np.arange(0.5, 1.81, 0.05):     # kt = true/model geocoronal intensity
            chi = 0
            for name, (obs, ineb) in paper.items():
                gm = gmod / ineb                   # model (subtracted) geocoronal relative to WIM
                b, _ = simulate(g=gm * kt, k=1 / kt, dm=dm)
                chi += ((b - obs) / 1.9) ** 2
            best.append((float(dm), float(kt), float(chi)))
    best = np.array(best)
    i = np.argmin(best[:, 2])
    out['best_fit_dm_kt'] = dict(dm=best[i, 0], k_true_over_model=best[i, 1], chi2=best[i, 2], dof=2)
    # chi2 for the two specific hypotheses with kt=1
    for dm in (0.0, -8.24):
        chi = 0
        bb = []
        for name, (obs, ineb) in paper.items():
            b, _ = simulate(g=gmod / ineb, k=1.0, dm=dm)
            bb.append(round(b, 2))
            chi += ((b - obs) / 1.9) ** 2
        out[f'hyp_dm{dm}_kt1'] = dict(pred=bb, chi2=round(chi, 2))
    # which true/model intensity would be needed with dm=0 (fine-structure-only) to give +5.5 km/s at g~1.1?
    need = []
    for kt in np.arange(0.2, 1.01, 0.05):
        b, _ = simulate(g=1.1 * kt, k=1 / kt, dm=0.0)
        need.append((round(float(kt), 2), round(b, 2)))
    out['dm0_bias_vs_true_over_model'] = need
    # solar Halpha absorption in the continuum
    sol = []
    for cs in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5):
        b, fw = simulate_with_solar(1.1, 1.0, 0.0, cs)
        sol.append(dict(csun_per_A=cs, bias=round(b, 2), fwhm_fit=round(fw, 1)))
    out['solar_Halpha_absorption'] = sol
    OUT.write_text(json.dumps(out, indent=1))
    print('DV_FS', round(DV_FS, 2))
    print(json.dumps(pred, indent=0))
    print(out['best_fit_dm_kt'], out['hyp_dm0.0_kt1'], out['hyp_dm-8.24_kt1'])
    print(need)
    print(sol)
    for r in grid:
        if r['sig_inst'] == 14.6 and r['k'] == 1.0:
            print(r)


if __name__ == '__main__':
    main()
