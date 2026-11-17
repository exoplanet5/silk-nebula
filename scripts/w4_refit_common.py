"""w4_refit_common: corrected (wave-3) observational likelihood for the Silk, shared by the w4_refit_* scripts.

All inputs are canonical wave-3 values (notes/w3_canon.md, w3_oii.md, w3_northsouth.md); energy ratios unless stated.
Every model (shock track, or an altmodel scenario) is scored on the same terms; a term is skipped if the model does not
predict that quantity (the number of terms used is reported).

SOUTH (Dec < 17.6-19.5, single-strand ridge):
  Hb/[O III]   DESI box 0.0122 +- 0.0414 (null-inflated; '[O III]/Hb > 11 at 2 sigma')          [verify_desi_audit_extra.json]
               variant B adds SDSS south [O III]/Hb = 22.8 (energy) +- 0.30 dex (bootstrap + placebo; 2 sigma > 5.3-10)
  [O II]3727/[O III]  1.14 +- 0.08 (DESI 1.12 +- 0.08 at Dec < 17.6; halves 1.17; SDSS 1.17)      [w3_northsouth]
  4363/5007    nominal 0.075 +- 0.030 (0.06-0.09 at 2-3 sigma); variant 'floor': 0.050 +- 0.035
               (one null box gives spurious 0.048, so up to ~0.05 may be spurious)                [verify_desi]
  [O II]7320+30/[O III] 0.076 +- 0.019 (0.071 south - 0.082 all; +-0.015 OH on 7325/3727)      [w3_oii]
  He II/[O III] 0 +- 0.013 (97.5 % < 0.026); He I 5876/[O III] 0 +- 0.0025 (< 0.005)
  [Ne III]/[O III] 0.085 +- 0.02 (0.07-0.10); Ha/[O III] 0.07 +- 0.09 ([O III]/Ha > 4 energy at 2 sigma)
  [N II]6583/[O III] 0.12 x/ 0.30 dex, [S II]/[O III] 0.05 x/ 0.30 dex (systematics-limited x2; wave-2 values)
  [O II]3729/3726 1.31 +- 0.07 (DESI realistic; 2.0 sigma below the LDL) -- reported as a separate, COMMON term
UV (south): FUV 80 +- 15 CU per R (68-92 per NSNS/imaging R; x1.25 higher on the DESI scale); NUV 47 +- 15 (36-57,
            up to 80 from the peak ratio)                                                         [verify_xrayfuv]
NORTH (Dec >= 19.5-21): [O III]/Hb SDSS 7.5 x/ 0.30 dex (16-84 % 2.9-15) and DESI Hb/[O III] 0 +- 0.06 (> 8.4 at 2 sigma);
  [O II]/[O III] 5.5 x/ 0.15 dex (DESI 4.0-5.3, SDSS 7.1-7.7). Brightness N/S 0.36-0.45; FUV/[O III] N/S 0.61 +- 0.15.
Model terms: fractional (linear terms) or dex (log terms) systematics as listed in TERMS.
"""
import numpy as np

# name: (kind, mu, sigma_obs, sigma_model); kind 'lin' (sigma_model fractional of the prediction) or 'log' (dex)
SOUTH_A = {
    'Hb_O3':  ('lin', 0.0122, 0.0414, 0.25),
    'O2_O3':  ('lin', 1.14, 0.08, 0.25),
    'O4363':  ('lin', 0.075, 0.030, 0.25),
    'O2r_O3': ('lin', 0.076, 0.019, 0.30),
    'HeII_O3': ('lin', 0.0, 0.013, 0.25),
    'HeI_O3': ('lin', 0.0, 0.0025, 0.25),
    'Ne3_O3': ('lin', 0.085, 0.020, 0.35),
    'Ha_O3':  ('lin', 0.07, 0.09, 0.25),
    'N2_O3':  ('log', 0.12, 0.30, 0.14),
    'S2_O3':  ('log', 0.05, 0.30, 0.20),
}
SOUTH_B = dict(SOUTH_A, O3_Hb_sdss=('log', 22.8, 0.30, 0.10))
DOUBLET = ('lin', 1.31, 0.07, 0.0)
NORTH = {
    'O3_Hb_sdss': ('log', 7.5, 0.30, 0.10),
    'Hb_O3': ('lin', 0.0, 0.06, 0.25),
    'O2_O3': ('log', 5.5, 0.15, 0.10),
    'HeII_O3': ('lin', 0.0, 0.03, 0.25),
}
UV_OBS = {'FUV': (80.0, 15.0), 'NUV': (47.0, 15.0), 'FUV_NS': (0.61, 0.15), 'I_NS': (0.36, 0.45)}


def derived(r):
    """Fill in derived keys (Hb_O3, Ha_O3, O3_Hb_sdss) from whatever the model provides."""
    r = dict(r)
    if 'Hb_O3' not in r and 'O3_Hb' in r and r['O3_Hb'] is not None:
        r['Hb_O3'] = 1.0 / r['O3_Hb']
    if 'O3_Hb_sdss' not in r and 'O3_Hb' in r and r['O3_Hb'] is not None:
        r['O3_Hb_sdss'] = r['O3_Hb']
    if 'Ha_O3' not in r and 'O3_Ha' in r and r['O3_Ha'] is not None:
        r['Ha_O3'] = 1.0 / r['O3_Ha']
    return r


def term(kind, mu, so, sm, x):
    x = np.asarray(x, float)
    if kind == 'lin':
        return (x - mu)**2 / (so**2 + (sm * x)**2)
    x = np.maximum(x, 1e-30)
    return (np.log10(x / mu))**2 / (so**2 + sm**2)


def chi2(r, obs=SOUTH_A, o4363=None):
    """Return (total chi2, dict of terms, n_terms).  o4363 = (mu, sigma) overrides the 4363 term (e.g. floor)."""
    r = derived(r)
    out = {}
    for k, (kind, mu, so, sm) in obs.items():
        if k == 'O4363' and o4363 is not None:
            mu, so = o4363
        if k not in r or r[k] is None:
            continue
        v = np.asarray(r[k], float)
        out[k] = np.where(np.isfinite(v), term(kind, mu, so, sm, v), 1e6)
    tot = sum(out.values()) if out else np.nan
    return tot, out, len(out)


def chi2_doublet(d):
    return float(term(*DOUBLET, d))


O4363_FLOOR = (0.050, 0.035)
