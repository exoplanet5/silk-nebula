"""w4_refit_widths: per-line predictions for the flat-sheet survivors (w4_refit_summary 'flat_plausible_L83' set):
[O III]/[O II] emission-weighted T and thermal FWHM (edge-on flat sheet: no bulk LOS flow), FUV line intensities after
C IV scattering (edge-on S_eff) in CU per R and in LU at I([O III]) = 1.62 R, C IV/Si IV and C IV/O III] photon ratios.
Output: data/derived/w4_refit_widths.json.  Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_widths.py
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import w4_refit_fit as WF
ROOT = '/Users/mickey/astronomy/silk-nebula'
K_B, MO = 1.380649e-16, 15.999 * 1.66054e-24
geo = json.load(open(ROOT + '/data/derived/w4_refit_geom.json'))
fit = json.load(open(ROOT + '/data/derived/w4_refit_fit.json'))
runs = {r['tag']: r for r in fit['runs']}
T = np.load(ROOT + '/data/derived/w4_refit_tracks.npz')
okf = lambda x: abs(x - 80) <= 2 * np.hypot(15, 0.2 * x)
S = [m for m in geo['models'] if (okf(m['flat']['CIV_span0']['FUV']) or okf(m['flat']['CIV_span1']['FUV']))
     and m['flat']['t_star_yr'] <= 1.2e5 and 0.83 <= m['flat']['asym_WE'] <= 1.13 and m['flat']['L_pc'] <= 83]
fw = lambda Tk: 2.3548 * np.sqrt(K_B * Tk / MO) / 1e5
out = []
CU_PER_R = {'CIV1549': 225.0}
for m in S:
    tr = WF.track(T, m['tag']); i = runs[m['tag']]['A_nat']['i']
    To3, To2 = float(tr['wo3_T'][i]), float(tr['wo2_T'][i])
    p = m['uv_parts']; Se = m['flat']['CIV_span0']['S_eff']
    civ = p['CIV1549'] * Se
    out.append(dict(tag=m['tag'], T_O3=To3, T_O2=To2, fwhm_O3=fw(To3), fwhm_O2=fw(To2), ratio_O2_O3_width=fw(To2) / fw(To3),
                    CIV_scat_CU=civ, SiIV_CU=p.get('SiIV1400', 0), OIII1663_CU=p.get('OIII1663', 0), HeII1640_CU=p.get('HeII1640', 0),
                    CIV_over_SiIV_CU=civ / max(p.get('SiIV1400', 1e-9), 1e-9), CIV_thin_over_SiIV=p['CIV1549'] / max(p.get('SiIV1400', 1e-9), 1e-9),
                    CIV_LU_at_1p62R=civ / 225.0 * 1.62 * 79577.47))
rng = lambda k: [float(np.min([o[k] for o in out])), float(np.median([o[k] for o in out])), float(np.max([o[k] for o in out]))]
res = {'n': len(out), 'ranges': {k: rng(k) for k in out[0] if k != 'tag'}, 'models': out,
       'warm_CIE_6e4K': {'fwhm_thermal_O': fw(6.08e4), 'O2_O3_width_ratio': 1.0}}
json.dump(res, open(ROOT + '/data/derived/w4_refit_widths.json', 'w'), indent=1)
print(len(out)); [print(k, [round(x, 3) for x in v]) for k, v in res['ranges'].items()]; print(res['warm_CIE_6e4K'])
