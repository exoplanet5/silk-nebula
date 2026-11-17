"""Robustness tests of the GALEX skybg on-filament contrast (reads data/derived/imaging_mw_galex_blocks.npz).

1. Remove the mean radial tile pattern (per band, per survey type) from every block (multiplicative, from blocks
   with |offset| > 0.8 deg so the filament itself does not bias it).
2. Per-tile local contrast: Delta_t = median(blocks |off - shift| <= 0.15) - median(blocks 0.30 < |off - shift| < 0.55)
   within the same tile (tile offsets cancel); combine tiles by median, error from tile-to-tile scatter.
3. Null test: repeat with the trace shifted by shift = +-0.75 ... +-3 deg (perpendicular offset), which samples
   real sky with the same tile pattern; the scatter of null contrasts is the realistic uncertainty.
Output: data/derived/imaging_mw_galex_tests.json
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import DER, get_spine, spine_s_of

B = dict(np.load(DER / 'imaging_mw_galex_blocks.npz'))
spine = get_spine()
s18 = spine_s_of(116.8, 18.0, spine)[0]
tiles = np.unique(B['tile'])
cra = np.zeros(B['tile'].max() + 1); cde = np.zeros_like(cra)
for t in tiles:
    m = B['tile'] == t
    cra[t] = np.mean(B['ra'][m]); cde[t] = np.mean(B['dec'][m])
r = np.hypot((B['ra'] - cra[B['tile']]) * np.cos(np.radians(B['dec'])), B['dec'] - cde[B['tile']])
B['r'] = r
# normalise each tile by its median, derive the mean radial pattern, divide it out
rel = np.zeros_like(B['val'])
for t in tiles:
    m = B['tile'] == t
    rel[m] = B['val'][m] / np.median(B['val'][m])
redges = np.linspace(0, 0.5, 11)
corr = np.ones_like(B['val'])
pattern = {}
for band in ['FUV', 'NUV']:
    for kind in ['AIS', 'MIS']:
        sel = (B['band'] == band) & (np.char.startswith(B['survey'].astype(str), 'AIS') == (kind == 'AIS'))
        far = sel & (np.abs(B['off']) > 0.8)
        rp = []
        for i in range(10):
            mm = far & (r >= redges[i]) & (r < redges[i + 1])
            rp.append(np.median(rel[mm]) if mm.sum() > 50 else 1.0)
        rp = np.array(rp)
        pattern[f'{band}:{kind}'] = rp.round(4).tolist()
        idx = np.clip(np.digitize(r[sel], redges) - 1, 0, 9)
        corr[sel] = rp[idx]
val = B['val'] / corr


def tile_contrasts(band, shift, s_lo=-1, s_hi=99, on_hw=0.15, ref=(0.30, 0.55)):
    off = B['off'] - shift
    sel = (B['band'] == band) & ~B['beyond'] & (B['s'] >= s_lo) & (B['s'] <= s_hi)
    out = []
    for t in np.unique(B['tile'][sel]):
        m = sel & (B['tile'] == t)
        on = m & (np.abs(off) <= on_hw)
        rf = m & (np.abs(off) > ref[0]) & (np.abs(off) < ref[1])
        if on.sum() >= 20 and rf.sum() >= 80:
            out.append((t, np.median(val[on]) - np.median(val[rf]), np.median(val[m])))
    return out


res = {'radial_pattern_rel': pattern, 'radial_edges_deg': redges.tolist(), 'segments': {}}
shifts_null = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.75, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
for band in ['FUV', 'NUV']:
    for seg, (lo, hi) in {'all': (-1, 99), 'south_dec12-18': (-1, s18), 'north_dec18-28': (s18, 99)}.items():
        c = tile_contrasts(band, 0.0, lo, hi)
        d = np.array([x[1] for x in c])
        nulls = []
        for sh in shifts_null:
            cn = tile_contrasts(band, sh, lo, hi)
            if len(cn) >= 3:
                nulls.append(np.median([x[1] for x in cn]))
        nulls = np.array(nulls)
        rec = {'n_tiles': int(len(d)), 'median_contrast_CU': float(np.median(d)) if len(d) else None,
               'mean_contrast_CU': float(np.mean(d)) if len(d) else None,
               'err_tile_scatter_CU': float(1.2533 * 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(len(d))) if len(d) > 2 else None,
               'tile_contrasts_CU': np.round(d, 1).tolist(),
               'null_shifts_deg': shifts_null, 'null_medians_CU': np.round(nulls, 1).tolist(),
               'null_rms_CU': float(np.std(nulls)) if len(nulls) > 2 else None,
               'null_mean_CU': float(np.mean(nulls)) if len(nulls) > 2 else None}
        res['segments'][f'{band}:{seg}'] = rec
        print(f"{band} {seg:15s} Ntile {rec['n_tiles']:3d} median {rec['median_contrast_CU']:+7.1f} "
              f"+-{rec['err_tile_scatter_CU'] or 0:5.1f} (tiles) | null mean {rec['null_mean_CU']:+6.1f} rms {rec['null_rms_CU']:5.1f}")
# FUV/NUV colour of the excess vs the diffuse background, south segment
f = res['segments']['FUV:south_dec12-18']; n = res['segments']['NUV:south_dec12-18']
res['note'] = ('Contrast = on (|off|<=0.15 deg) minus local reference (0.30-0.55 deg) inside the same tile after removing the '
               'mean radial tile pattern. Significance must be judged against null_rms (trace shifted onto blank sky).')
(DER / 'imaging_mw_galex_tests.json').write_text(json.dumps(res, indent=1))
