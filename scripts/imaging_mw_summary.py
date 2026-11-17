"""Compile per-band filament contrasts into upper limits / detections with physical conversions.

For each band: on-filament contrast C (on |off|<=h_on minus local reference, whole spine and south segment),
null rms sigma (trace shifted 1.5-4.5 deg), 2-sigma upper limit U = max(C, 0) + 2 sigma (beam-diluted), and an
intrinsic (peak) limit U / D, with D = f / sqrt(f^2 + b^2) the 1-D dilution of a Gaussian filament of FWHM
f = 2 x 0.138 deg (E-strand HWHM, silk_trace.json) by a beam of FWHM b.
Radio: T_b limit -> surface brightness Sigma_nu = 2 k T nu^2 / c^2, scaled to 1 GHz with alpha = 0.5 (S ~ nu^-alpha).
Output: data/derived/imaging_mw_verdicts.json (+ printed table)
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import DER, get_spine

P = json.loads((DER / 'imaging_mw_profiles.json').read_text())
E = json.loads((DER / 'imaging_mw_edges.json').read_text())
T = json.loads((DER / 'imaging_mw_galex_tests.json').read_text())
F = json.loads((DER / 'imaging_mw_fuv_prediction.json').read_text())
sp = get_spine()
f_fwhm = 2 * (sp[2].get('half_width_deg') or 0.138)
k, c = 1.380649e-23, 2.998e8
NU = {'haslam408': 0.408e9, 'stockert1420': 1.42e9, 'drao22': 0.022e9, 'planck030': 28.4e9, 'wmapk': 22.8e9}
TSCALE = {'haslam408': 1.0, 'stockert1420': 1e-3, 'drao22': 1.0, 'planck030': 1.0, 'wmapk': 1e-3}  # -> K


def sigma_1ghz(T_K, nu, alpha=0.5):
    I = 2 * k * T_K * nu ** 2 / c ** 2       # W m-2 Hz-1 sr-1 (Rayleigh-Jeans)
    return I * (1e9 / nu) ** (-alpha)


out = {'filament_fwhm_deg': f_fwhm, 'bands': {}}
print(f'{"band":20s} {"seg":6s} {"C":>11s} {"null_rms":>10s} {"sig":>6s} {"2sig_UL":>10s} {"dil":>5s} {"intrinsic_UL":>12s}  unit')
for name, b in P['bands'].items():
    rec = {'label': b['label'], 'unit': b['unit'], 'beam_deg': b['beam_deg'], 'category': b['category']}
    D = f_fwhm / np.hypot(f_fwhm, b['beam_deg'])
    rec['dilution'] = D
    for seg in ['all', 'south', 'north']:
        sg = b['segments'].get(seg)
        if not sg:
            continue
        C, s = sg['contrast'] - sg['contrast_null_mean'], sg['contrast_null_rms']
        U = max(C, 0) + 2 * s
        r = {'contrast_minus_nullmean': C, 'null_rms': s, 'sig': C / s if s > 0 else None, 'UL2': U, 'UL2_intrinsic': U / D,
             'step_W_minus_E': sg['step_W_minus_E'] - sg['step_null_mean'], 'step_null_rms': sg['step_null_rms'],
             'step_sig': sg['step_sig'], 'on': sg['on'], 'ref': sg['ref']}
        if name in NU:
            r['Sigma1GHz_UL2_W_m2_Hz_sr'] = sigma_1ghz(U / D * TSCALE[name], NU[name])
        rec[seg] = r
        if seg in ('all', 'south'):
            print(f"{name:20s} {seg:6s} {C:+11.4g} {s:10.4g} {r['sig']:+6.1f} {U:10.4g} {D:5.2f} {U / D:12.4g}  {b['unit']}"
                  + (f"  Sigma1GHz<{r['Sigma1GHz_UL2_W_m2_Hz_sr']:.2g}" if 'Sigma1GHz_UL2_W_m2_Hz_sr' in r else ''))
    rec['regions'] = b['regions']
    out['bands'][name] = rec

# edges
edge = {bnd: {kk: v for kk, v in E['along_spine'][bnd].items() if kk != 'segments'} for bnd in E['along_spine']}
wf = [s['wall_over_floor'] for s in E['along_spine']['galfa_nhi_oiii']['segments'] if s['wall_over_floor']]
edge['galfa_nhi_oiii']['wall_over_floor_median'] = float(np.median(wf)); edge['galfa_nhi_oiii']['wall_over_floor_range'] = [float(min(wf)), float(max(wf))]
out['edges'] = edge
ring = {}
for bnd in ['ero_0206_clean', 'rass_r12', 'rass_r45', 'ero_0610_clean']:
    ring[bnd] = [{kk: r[kk] for kk in ('pa0', 'r_peak', 'r_edge', 'peak_over_outside', 'silk_r_in_sector')} for r in E['monogem_ring'][bnd]
                 if r['silk_r_in_sector']]
out['monogem_ring_sectors'] = ring
out['monogem_ring_silk_radius'] = E['monogem_ring']['silk_radius_deg']
out['galex_tile_tests'] = {kk: {x: v[x] for x in ('n_tiles', 'median_contrast_CU', 'err_tile_scatter_CU', 'null_rms_CU', 'null_mean_CU')}
                           for kk, v in T['segments'].items()}
out['fuv_prediction'] = F['predictions_CU']
# FUV in Rayleigh-equivalent for a 1550 A line, W_eff = 300 A
R2PH = 1e6 / (4 * np.pi)
for kk, (v, e) in F['measured_CU_value_nullrms'].items():
    out.setdefault('fuv_measured_R_equiv', {})[kk] = {'CU': v, 'null_rms_CU': e, 'UL2_CU': max(v, 0) + 2 * e,
                                                       'UL2_R_at1549': (max(v, 0) + 2 * e) * 300 / R2PH}
(DER / 'imaging_mw_verdicts.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out['edges'], indent=0, default=float)[:1500])
print(json.dumps(out['fuv_measured_R_equiv'], indent=0))
