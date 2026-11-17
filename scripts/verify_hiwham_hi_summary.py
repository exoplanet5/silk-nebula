"""verify_hiwham: summarise the HI profile statistics against the two null families.

Reads data/derived/verify_hiwham_hi_profiles.npz (from verify_hiwham_hi.py) and writes
data/derived/verify_hiwham_hi_results.json.

Significances use ROBUST null widths, sigma_r = (p84 - p16)/2, and the null median; empirical tail fractions
are also given. For the 11-segment wave-1 'wall' statistic the null is built by drawing 11 independent random
segment placements (one per real segment shape) 4000 times, and also from the 120 coherent whole-trace copies.
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import DER

Z = np.load(DER / 'verify_hiwham_hi_profiles.npz', allow_pickle=False)
names = list(Z['names']); groups = list(Z['groups']); STATS = list(Z['stats'])
STr, STc = Z['ST_real'], Z['ST_coh']            # (ng, nm, ns), (ncoh, ng, nm, ns)
STrs, STs, seg_of = Z['ST_real_seg'], Z['ST_seg'], Z['seg_of']     # (11, nm, ns), (nnull, nm, ns)
rng = np.random.default_rng(3)


def rob(x):
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan, np.nan, 0
    return float(np.median(x)), float(0.5 * (np.percentile(x, 84) - np.percentile(x, 16))), int(len(x))


res = {'maps': names, 'groups': groups, 'stats': STATS, 'n_coherent': int(STc.shape[0]), 'n_segment_nulls': int(STs.shape[0]),
       'table': {}}
for j, nm in enumerate(names):
    res['table'][nm] = {}
    for g in ('all', 'south', 'north'):
        gi = groups.index(g)
        row = {}
        for k, st in enumerate(STATS):
            real = float(STr[gi, j, k])
            med, sr, n = rob(STc[:, gi, j, k])
            row[st] = {'real': real, 'coh_median': med, 'coh_sig': sr, 'z_coh': (real - med) / sr if sr else np.nan,
                       'frac_coh_ge': float(np.mean(STc[:, gi, j, k][np.isfinite(STc[:, gi, j, k])] >= real))}
        res['table'][nm][g] = row

# ---------------------------------------------------------------- 11-segment statistics (wave-1 wall test)
seg_tests = {}
for nm in ['old_oiii', 'oiii_helio', 'red', 'all', 'W_helio', 'lsr+11+13', 'lsr+13+15', 'lsr+1+3', 'lsr-1+1']:
    if nm not in names:
        continue
    j = names.index(nm)
    kw, kr = STATS.index('w1_wall'), STATS.index('w1_ratio')
    ks, kd = STATS.index('step_amp'), STATS.index('step_d0')
    real_w = STr[:11, j, kw]                     # from the full-trace grouping (same as wave-1)
    real_w_piece = STrs[:, j, kw]
    real_r = STr[:11, j, kr]
    # independent-segment null: draw one null of each segment's shape
    by = [np.where(seg_of == i)[0] for i in range(11)]
    draws_mean, draws_sem, draws_nin, draws_ratio = [], [], [], []
    for it in range(4000):
        pick = np.array([rng.choice(b) for b in by])
        w = STs[pick, j, kw]; r = STs[pick, j, kr]
        ok = np.isfinite(w)
        draws_mean.append(np.nanmean(w)); draws_sem.append(np.nanstd(w) / np.sqrt(max(ok.sum(), 1)))
        draws_nin.append(np.sum((w > 0.4) & (w < 1.4))); draws_ratio.append(np.nanmedian(r))
    draws_mean, draws_sem, draws_nin, draws_ratio = map(np.array, (draws_mean, draws_sem, draws_nin, draws_ratio))
    coh_w = STc[:, :11, j, kw]
    coh_mean = np.nanmean(coh_w, axis=1)
    coh_nin = np.sum((coh_w > 0.4) & (coh_w < 1.4), axis=1)
    rm = float(np.nanmean(real_w)); rs = float(np.nanstd(real_w) / np.sqrt(np.isfinite(real_w).sum()))
    seg_tests[nm] = {
        'real_wall_offsets': [None if not np.isfinite(x) else round(float(x), 3) for x in real_w],
        'real_wall_offsets_piecewise': [None if not np.isfinite(x) else round(float(x), 3) for x in real_w_piece],
        'real_ratio': [None if not np.isfinite(x) else round(float(x), 3) for x in real_r],
        'real_mean': rm, 'real_sem': rs, 'real_n_in_0.4_1.4': int(np.sum((real_w > 0.4) & (real_w < 1.4))),
        'real_median_ratio': float(np.nanmedian(real_r)),
        'indep_null_mean_of_mean': float(np.mean(draws_mean)), 'indep_null_std_of_mean': float(np.std(draws_mean)),
        'indep_null_p05_p95_mean': [float(np.percentile(draws_mean, 5)), float(np.percentile(draws_mean, 95))],
        'indep_null_median_sem': float(np.median(draws_sem)),
        'indep_null_frac_mean_within_0.15_of_real': float(np.mean(np.abs(draws_mean - rm) < 0.15)),
        'indep_null_frac_n_in_ge_real': float(np.mean(draws_nin >= np.sum((real_w > 0.4) & (real_w < 1.4)))),
        'indep_null_median_ratio': float(np.median(draws_ratio)),
        'indep_null_frac_ratio_ge_real': float(np.mean(draws_ratio >= np.nanmedian(real_r))),
        'single_seg_null_wall_median': float(np.nanmedian(STs[:, j, kw])),
        'single_seg_null_wall_p16_p84': [float(np.nanpercentile(STs[:, j, kw], 16)), float(np.nanpercentile(STs[:, j, kw], 84))],
        'single_seg_null_frac_found': float(np.mean(np.isfinite(STs[:, j, kw]))),
        'coh_null_mean_of_mean': float(np.nanmean(coh_mean)), 'coh_null_std_of_mean': float(np.nanstd(coh_mean)),
        'coh_null_frac_n_in_ge_real': float(np.mean(coh_nin >= np.sum((real_w > 0.4) & (real_w < 1.4)))),
        'real_step_amp_segments': [round(float(x), 3) for x in STr[:11, j, ks] / 1e19],
        'real_step_d0_segments': [round(float(x), 3) for x in STr[:11, j, kd]],
    }
res['wave1_wall_test'] = seg_tests


# ---------------------------------------------------------------- per-segment filament amplitude vs segment nulls
seg_fil = {}
ka = STATS.index('A_fil'); kS = STATS.index('S_loc')
for nm in names:
    j = names.index(nm)
    zs, ss = [], []
    for i in range(11):
        nl = STs[seg_of == i, j, ka]
        med, sr, n = rob(nl)
        zs.append((STrs[i, j, ka] - med) / sr)
        nl2 = STs[seg_of == i, j, kS]; m2, s2, _ = rob(nl2)
        ss.append((STrs[i, j, kS] - m2) / s2)
    seg_fil[nm] = {'z_Afil_segments': [round(float(x), 2) for x in zs], 'z_Sloc_segments': [round(float(x), 2) for x in ss],
                   'chi2_Afil': float(np.nansum(np.array(zs) ** 2)), 'sum_z_Afil_over_sqrt11': float(np.nansum(zs) / np.sqrt(11)),
                   'sum_z_Sloc_over_sqrt11': float(np.nansum(ss) / np.sqrt(11))}
res['segment_z'] = seg_fil
(DER / 'verify_hiwham_hi_results.json').write_text(json.dumps(res, indent=1, default=float))

# ---------------------------------------------------------------- console summary
print(f"coherent nulls {STc.shape[0]}, segment nulls {STs.shape[0]}")
print('\nWhole-spine / south / north, A_fil (1e19 cm^-2) and S_loc, S_wide (1e19) with coherent-null z')
for nm in names:
    line = f'{nm:12s}'
    for g in ('all', 'south', 'north'):
        r = res['table'][nm][g]
        line += f" | {g[:2]} A {r['A_fil']['real'] / 1e19:+5.2f} ({r['A_fil']['z_coh']:+4.1f}) Sl {r['S_loc']['real'] / 1e19:+5.2f} ({r['S_loc']['z_coh']:+4.1f}) Sw {r['S_wide']['real'] / 1e19:+5.2f} ({r['S_wide']['z_coh']:+4.1f})"
    print(line)
print('\nstep fit (whole spine): amp 1e19 / d0 / z_coh(amp)')
for nm in names:
    r = res['table'][nm]['all']
    print(f"{nm:12s} amp {r['step_amp']['real'] / 1e19:+5.2f} d0 {r['step_d0']['real']:+5.2f} z {r['step_amp']['z_coh']:+4.1f} | "
          f"coh d0 median {r['step_d0']['coh_median']:+.2f}  | 1-sig A_fil coh {r['A_fil']['coh_sig'] / 1e19:.2f}")
print('\nwave-1 wall test')
for nm, t in seg_tests.items():
    print(f"{nm:12s} real mean {t['real_mean']:+.2f}+-{t['real_sem']:.2f} n_in {t['real_n_in_0.4_1.4']}  ratio {t['real_median_ratio']:.2f} | "
          f"indep null mean {t['indep_null_mean_of_mean']:+.2f}+-{t['indep_null_std_of_mean']:.2f} sem~{t['indep_null_median_sem']:.2f} "
          f"P(n_in>=real) {t['indep_null_frac_n_in_ge_real']:.2f} ratio {t['indep_null_median_ratio']:.2f} P(ratio>=) {t['indep_null_frac_ratio_ge_real']:.2f} | "
          f"coh mean {t['coh_null_mean_of_mean']:+.2f}+-{t['coh_null_std_of_mean']:.2f} P(n_in>=) {t['coh_null_frac_n_in_ge_real']:.2f}")
    print('      real offsets', t['real_wall_offsets'])
print('\nsegment z (A_fil, S_loc) combined')
for nm, t in seg_fil.items():
    print(f"{nm:12s} A_fil sumz/sqrt11 {t['sum_z_Afil_over_sqrt11']:+.2f} chi2 {t['chi2_Afil']:.1f} | S_loc sumz/sqrt11 {t['sum_z_Sloc_over_sqrt11']:+.2f}  zA {t['z_Afil_segments']}")
