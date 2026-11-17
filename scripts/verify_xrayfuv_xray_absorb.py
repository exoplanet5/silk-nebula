"""verify_xrayfuv: can foreground absorption make the soft X-ray step at the Silk?

(1) HI column profiles across the spine: HI4PI (16') from the wide SkyView map and GALFA-HI DR2
    (4', wave-1 EQ grid imaging_mw_eq_galfa_nhi_all.fits, -30..+40 km/s) -> Delta N(HI) = W - E
    (x in -2.5..-0.5 vs +0.5..+2.5 deg) and the step c at x=0 (same estimator as the X-ray test).
(2) Implied fractional change of the 0.2-0.6 keV brightness, exp(-sigma f DeltaN) - 1, for
    sigma_eff = 0.5-1.0e-20 cm^2 (approximate band-averaged cross-section for a kT~0.14 keV plasma
    in the eROSITA 0.2-0.6 keV band; literature-order value, not computed from a response) and a
    foreground fraction f = 0.15 (3D-dust E(B-V)<300 pc / SFD total, context note), 0.3, 1.
(3) Re-run the step statistic with realistic absorption corrections sigma*f = 1.5e-21, 3e-21.
Output: data/derived/verify_xrayfuv_xray_absorb.json
"""
import sys, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import *
import verify_xrayfuv_xray_edge as X

ra, dec, _ = load_spine()
E = X.load_erass(); pl = X.Placer(E)
cc = CurveCoords(ra, dec); s_kink = float(np.interp(17.6, cc.dec, cc.s))
segs = {'all': (None, None), 'south': (0.0, s_kink), 'north': (s_kink, float(cc.s[-1]))}
out = {}
# HI4PI on the eRASS pixel table
for seg, (a, b) in segs.items():
    p = pl.profile(ra, dec, E['nhi'].astype(float), np.isfinite(E['nhi']), 0.1, a, b)
    st = X.stats(p[0], p[1])
    xin = (p[0] >= -2.5) & (p[0] <= -0.5); xout = (p[0] >= 0.5) & (p[0] <= 2.5)
    dN = float(np.nanmean(p[1][xin]) - np.nanmean(p[1][xout]))
    out[f'HI4PI_{seg}'] = {'mean_NHI': st['level'], 'dN_W_minus_E': dN, 'step_c': st['c']}
# GALFA (EQ grid from wave-1, read only)
g = fits.open(DER / 'imaging_mw_grids/imaging_mw_eq_galfa_nhi_all.fits')[0]
wg = WCS(g.header)
ny, nx = g.data.shape
yy, xx = np.mgrid[0:ny, 0:nx]
gra, gde = wg.pixel_to_world_values(xx, yy)
GD = {'ra': gra.ravel(), 'dec': gde.ravel(), 'vec': radec_to_vec(gra.ravel(), gde.ravel())}
plg = X.Placer(GD)
img = g.data.astype(float).ravel()
for seg, (a, b) in segs.items():
    p = plg.profile(ra, dec, img, np.isfinite(img) & (img > 0), 0.1, a, b)
    xin = (p[0] >= -2.5) & (p[0] <= -0.5); xout = (p[0] >= 0.5) & (p[0] <= 2.5)
    st = X.stats(p[0], p[1])
    out[f'GALFA_all_{seg}'] = {'mean_NHI': float(np.nanmean(p[1][np.abs(p[0]) <= 2.5])),
                               'dN_W_minus_E': float(np.nanmean(p[1][xin]) - np.nanmean(p[1][xout])),
                               'step_c': st['c'] if st else None}
print(json.dumps(out, indent=1))
imp = {}
for k, v in out.items():
    for sig in (0.5e-20, 1.0e-20):
        for f in (0.15, 0.3, 1.0):
            imp[f'{k}|sigma{sig:.1e}|f{f}'] = float(np.exp(-sig * f * v['dN_W_minus_E']) - 1)
out['implied_fractional_WE_change'] = imp
# realistic absorption-corrected step statistics with nulls
for sf in (1.5e-21, 3e-21):
    img2 = E['soft'] * np.exp(sf * np.nan_to_num(E['nhi'].astype(float)))
    for seg, (a, b) in segs.items():
        r = X.run_placements(pl, img2, E['mask_full'], 0.1, ra, dec, 'soft', 'PSR', a, b)
        out[f'abscorr_{sf:.1e}_{seg}'] = {'real': r['real'], 'c': X.summarize(r, 'c'), 'WE': X.summarize(r, 'WE'),
                                          'cfrac': X.summarize(r, 'cfrac')}
        print(sf, seg, r['real']['WE'], r['real']['cfrac'], out[f'abscorr_{sf:.1e}_{seg}']['c'], out[f'abscorr_{sf:.1e}_{seg}']['WE'], flush=True)
(DER / 'verify_xrayfuv_xray_absorb.json').write_text(json.dumps(out, indent=1, default=float))
