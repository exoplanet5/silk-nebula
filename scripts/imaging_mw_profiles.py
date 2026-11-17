"""Cross-filament profiles, on/off contrasts, interior/exterior steps and null tests for every band (EQ grid).

For each band (grid on the 3' EQ grid):
  * perpendicular offset `off` from the E-strand spine (+ = east = convex/exterior side; the spine's
    curvature centre, and PSR B0656+14 / the Monogem Ring centre, lie to the WEST);
  * binned median profiles in segments: all (Dec 12-28.25), south (Dec 12-18), north (Dec 18-28.25), and
    +-0.6 deg of spine around P1, P2, P3;
  * contrast  C = <on> - <ref>, on = |off| <= h_on, ref = h_on + b/2 + 0.15 < |off| < ... + 0.5 deg
    (h_on = max(0.15, b/2), b = beam FWHM);
  * step      S = <-2.5 < off < -0.5> - <0.5 < off < 2.5>  (west/interior minus east/exterior);
  * null tests: same statistics with the trace shifted by +-1.5 ... +-4.5 deg (perpendicular), which samples
    unrelated sky with the same map systematics; significance = (value - null mean) / null rms;
  * P1-P4 disk medians (radius max(0.25, b/2) deg).
Outputs: data/derived/imaging_mw_profiles.json, data/derived/imaging_mw_band_summary.csv
"""
import csv
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import (DER, load_grid, grid_coords, get_spine, trace_coords, spine_s_of, disk_stat, REGIONS)

# name, label, unit, beam FWHM (deg), category
BANDS = [
    ('murthy_fuv', 'GALEX FUV (Murthy14)', 'CU', 0.15, 'UV'),
    ('murthy_nuv', 'GALEX NUV (Murthy14)', 'CU', 0.15, 'UV'),
    ('rass_r12', 'RASS R12 1/4 keV', '1e-6 c/s/arcmin2', 0.4, 'X'),
    ('rass_r45', 'RASS R45 3/4 keV', '1e-6 c/s/arcmin2', 0.4, 'X'),
    ('rass_r67', 'RASS R67 1.5 keV', '1e-6 c/s/arcmin2', 0.4, 'X'),
    ('rass_int_soft_clean', 'RASS3 soft (cleaned)', 'SkyView int', 0.35, 'X'),
    ('ero_0204_clean', 'eRASS1 0.2-0.4 keV', 'cts/s/deg2', 0.35, 'X'),
    ('ero_0406_clean', 'eRASS1 0.4-0.6 keV', 'cts/s/deg2', 0.35, 'X'),
    ('ero_0206_clean', 'eRASS1 0.2-0.6 keV', 'cts/s/deg2', 0.35, 'X'),
    ('ero_0610_clean', 'eRASS1 0.6-1.0 keV', 'cts/s/deg2', 0.35, 'X'),
    ('ero_1023_clean', 'eRASS1 1.0-2.3 keV', 'cts/s/deg2', 0.35, 'X'),
    ('ero_ovii_clean', 'eRASS1 O VII band', 'cts/s/deg2', 0.35, 'X'),
    ('ero_oviii_clean', 'eRASS1 O VIII band', 'cts/s/deg2', 0.35, 'X'),
    ('iris100', 'IRIS 100 um', 'MJy/sr', 0.07, 'IR'),
    ('iris60', 'IRIS 60 um', 'MJy/sr', 0.07, 'IR'),
    ('iris12', 'IRIS 12 um', 'MJy/sr', 0.07, 'IR'),
    ('planck857', 'Planck 857 GHz', 'MJy/sr', 0.08, 'IR'),
    ('planck545', 'Planck 545 GHz', 'MJy/sr', 0.08, 'IR'),
    ('sfd_ebv', 'SFD E(B-V)', 'mag', 0.1, 'IR'),
    ('planck353_pi', 'Planck 353 GHz PI', 'K_CMB', 0.25, 'IR'),
    ('galfa_nhi_all', 'GALFA NHI -30..+40', 'cm-2', 0.07, 'HI'),
    ('galfa_nhi_oiii', 'GALFA NHI +11..+16', 'cm-2', 0.07, 'HI'),
    ('galfa_nhi_red', 'GALFA NHI -1..+6', 'cm-2', 0.07, 'HI'),
    ('galfa_nhi_m10m02', 'GALFA NHI -10..-2', 'cm-2', 0.07, 'HI'),
    ('galfa_nhi_p06p10', 'GALFA NHI +6..+10', 'cm-2', 0.07, 'HI'),
    ('galfa_nhi_p18p26', 'GALFA NHI +18..+26', 'cm-2', 0.07, 'HI'),
    ('galfa_mom1', 'GALFA moment-1', 'km/s', 0.07, 'HI'),
    ('hi4pi_nhi', 'HI4PI NHI', 'cm-2', 0.27, 'HI'),
    ('halpha_f03', 'Halpha (Finkbeiner03)', 'R', 0.1, 'Ha'),
    ('haslam408', 'Haslam 408 MHz', 'K', 0.85, 'radio'),
    ('stockert1420', 'Stockert 1420 MHz', 'mK', 0.6, 'radio'),
    ('drao22', 'DRAO 22 MHz', 'K', 1.5, 'radio'),
    ('planck030', 'Planck 30 GHz I', 'K_CMB', 0.55, 'radio'),
    ('planck030_pi', 'Planck 30 GHz PI', 'K_CMB', 0.55, 'radio'),
    ('wmapk', 'WMAP K 23 GHz', 'mK', 0.88, 'radio'),
    ('gb6_4850', 'GB6 4850 MHz', 'Jy/beam', 0.06, 'radio'),
    ('nvss', 'NVSS 1.4 GHz', 'Jy/beam', 0.05, 'radio'),
    ('tgss', 'TGSS 150 MHz', 'Jy/beam', 0.05, 'radio'),
    ('fermi3', 'Fermi-LAT band 3', 'ph/cm2/s/sr', 1.0, 'gamma'),
    ('fermi4', 'Fermi-LAT band 4', 'ph/cm2/s/sr', 0.5, 'gamma'),
    ('fermi5', 'Fermi-LAT band 5', 'ph/cm2/s/sr', 0.3, 'gamma'),
]
NULL_SHIFTS = [-4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]

lon, lat = grid_coords('eq')
spine = get_spine()
off, s, beyond, L = trace_coords(lon, lat, 'eq', spine)
np.savez_compressed(DER / 'imaging_mw_eq_tracecoords.npz', off=off.astype(np.float32), s=s.astype(np.float32),
                    beyond=beyond, L=L)
s18 = spine_s_of(116.8, 18.0, spine)[0]
SEGS = {'all': (0.0, L), 'south': (0.0, s18), 'north': (s18, L)}
for k in ['P1', 'P2', 'P3']:
    s0, o0 = spine_s_of(*REGIONS[k], spine=spine)
    SEGS[k] = (s0 - 0.6, s0 + 0.6)
REG_OFF = {k: spine_s_of(*REGIONS[k], spine=spine) for k in REGIONS}


def stats(img, offs, sel, b):
    h_on = max(0.15, b / 2)
    r0 = h_on + b / 2 + 0.15
    on = sel & (np.abs(offs) <= h_on)
    ref = sel & (np.abs(offs) > r0) & (np.abs(offs) < r0 + 0.5)
    west = sel & (offs > -2.5) & (offs < -0.5)
    east = sel & (offs > 0.5) & (offs < 2.5)
    f = lambda m: np.nanmedian(img[m]) if np.isfinite(img[m]).sum() > 10 else np.nan
    return f(on) - f(ref), f(west) - f(east), f(on), f(ref)


def profile(img, sel, edges):
    cen = 0.5 * (edges[1:] + edges[:-1]); v = np.full(len(cen), np.nan); e = np.full(len(cen), np.nan)
    for i in range(len(cen)):
        m = sel & (off >= edges[i]) & (off < edges[i + 1]) & np.isfinite(img)
        x = img[m]
        if x.size >= 5:
            v[i] = np.median(x)
            e[i] = 1.2533 * 1.4826 * np.median(np.abs(x - v[i])) / np.sqrt(max(1, x.size / 9))  # 3x3 px correlation
    return cen, v, e


out = {'trace': {'source': spine[2]['source'], 'version': spine[2].get('version'), 'length_deg': L,
                 'regions_s_off': REG_OFF, 'segments': SEGS},
       'method': __doc__, 'bands': {}}
rows = []
for name, lab, unit, b, cat in BANDS:
    try:
        img = load_grid(name, 'eq')
    except FileNotFoundError:
        print('missing', name); continue
    rec = {'label': lab, 'unit': unit, 'beam_deg': b, 'category': cat, 'segments': {}}
    db = 0.1 if b < 0.3 else 0.2
    edges = np.arange(-4.0, 4.0001, db)
    for seg, (lo, hi) in SEGS.items():
        sel = (s >= lo) & (s <= hi) & ~beyond & np.isfinite(img)
        if sel.sum() < 50:
            continue
        C, S, von, vref = stats(img, off, sel, b)
        nulls = []
        for sh in NULL_SHIFTS:
            nulls.append(stats(img, off - sh, sel, b)[:2])
        nulls = np.array(nulls)
        cm, cs = np.nanmean(nulls[:, 0]), np.nanstd(nulls[:, 0])
        sm, ss = np.nanmean(nulls[:, 1]), np.nanstd(nulls[:, 1])
        cen, v, e = profile(img, sel, edges)
        rec['segments'][seg] = {'contrast': float(C), 'contrast_null_mean': float(cm), 'contrast_null_rms': float(cs),
                                'contrast_sig': float((C - cm) / cs) if cs > 0 else None,
                                'step_W_minus_E': float(S), 'step_null_mean': float(sm), 'step_null_rms': float(ss),
                                'step_sig': float((S - sm) / ss) if ss > 0 else None,
                                'on': float(von), 'ref': float(vref),
                                'profile_offsets': cen.round(3).tolist(), 'profile': v.tolist(), 'profile_err': e.tolist()}
        if seg in ('all', 'south', 'north'):
            rows.append([name, lab, unit, cat, b, seg, von, vref, C, cs, rec['segments'][seg]['contrast_sig'], S, ss,
                         rec['segments'][seg]['step_sig']])
    disk = {}
    for k, (r, d) in REGIONS.items():
        med, sig, n = disk_stat(img, lon, lat, r, d, max(0.25, b / 2))
        disk[k] = {'median': float(med), 'robust_sigma': float(sig), 'npix': int(n)}
    for k in ['P1', 'P2', 'P3']:
        disk[k]['minus_P4'] = disk[k]['median'] - disk['P4']['median']
    rec['regions'] = disk
    out['bands'][name] = rec
    a = rec['segments'].get('all', {})
    print(f"{name:20s} C={a.get('contrast', np.nan):+.3g} ({a.get('contrast_sig') or np.nan:+.1f} sig)  "
          f"W-E={a.get('step_W_minus_E', np.nan):+.3g} ({a.get('step_sig') or np.nan:+.1f} sig)  "
          f"P1-P4={disk['P1']['minus_P4']:+.3g} P3-P4={disk['P3']['minus_P4']:+.3g}")

(DER / 'imaging_mw_profiles.json').write_text(json.dumps(out, indent=1, default=float))
with open(DER / 'imaging_mw_band_summary.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['grid', 'label', 'unit', 'category', 'beam_deg', 'segment', 'on_median', 'ref_median', 'contrast_on_minus_ref',
                'contrast_null_rms', 'contrast_sigma', 'step_west_minus_east', 'step_null_rms', 'step_sigma'])
    for r in rows:
        w.writerow([r[0], r[1], r[2], r[3], r[4], r[5]] + [f'{x:.5g}' if isinstance(x, (float, np.floating)) and x is not None else x for x in r[6:]])
print('wrote', len(rows), 'rows')
