"""verify_precovery step 8: how the PSP sky estimator (Field.sky_<b>) responds along the scan.

For every run/camcol series in the box, fit sky_b[f] = c0 + sum_{k=-K..K} w_k skyFrames_b[f+k] (least squares on the
first differences, which removes slow trends), pooled over all crossing runs. The kernel w_k shows whether PSP is a
temporally smoothed version of the per-frame sky (which would dilute a filament that each camcol crosses in 2-3 fields).
Also reports the high-frequency (second-difference) rms ratio PSP / Frames.
Output: data/derived/verify_precovery_pspkernel.json
"""
import sys, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from verify_precovery_lib import load_fields, DER

a = load_fields()
K = 4
RUNS = [2711, 2740, 2777, 2830, 2864, 2887, 2888, 3462, 3513, 3560, 3605, 3614, 3635, 3644, 4335, 4458, 4469, 4508,
        4517, 4555, 4569, 4576, 5037, 5045, 5046, 6597, 6605, 6728]
res = {}
for b in 'ugriz':
    X, Y, hf = [], [], []
    for r in RUNS:
        for c in range(1, 7):
            k = (a['run'] == r) & (a['camcol'] == c)
            o = np.argsort(a['field'][k])
            p = a[f'sky_{b}'][k][o]; fr = a[f'skyFrames_{b}'][k][o]; fl = a['field'][k][o]
            if len(p) < 2 * K + 10 or np.any(np.diff(fl) != 1):
                continue
            dp = np.diff(p); dfr = np.diff(fr)
            for i in range(K, len(dp) - K):
                X.append(dfr[i - K:i + K + 1]); Y.append(dp[i])
            hf.append((np.std(np.diff(p, 2)), np.std(np.diff(fr, 2))))
    X = np.array(X); Y = np.array(Y)
    ok = np.isfinite(X).all(1) & np.isfinite(Y)
    w, *_ = np.linalg.lstsq(X[ok], Y[ok], rcond=None)
    hf = np.array(hf)
    res[b] = dict(kernel_lags=list(range(-K, K + 1)), kernel=[float(v) for v in w], kernel_sum=float(w.sum()),
                  n=int(ok.sum()), hf_rms_ratio_psp_over_frames=float(np.median(hf[:, 0] / hf[:, 1])))
    print(b, 'kernel', np.round(w, 3), 'sum', round(float(w.sum()), 3), 'HF rms PSP/Frames', round(res[b]['hf_rms_ratio_psp_over_frames'], 3))
(DER / 'verify_precovery_pspkernel.json').write_text(json.dumps(res, indent=1))


# ---------------------------------------------------------------------------------------------------------------
# Part 2: does a time-smoothed template fit the PSP (and Frames) sky better?  T_s[f,c] = sum_k G_sigma(k) T[f+k,c]
# (Gaussian in field number, normalised; template values of neighbouring fields of the same camcol). If PSP were a
# temporally smoothed sky, its amplitude and chi2 would improve with sigma > 0 relative to Frames.
def smoothed_template(idx, T, sig):
    if sig == 0:
        return T.copy()
    import verify_precovery_fit as VF
    Z = VF.Z
    out = np.zeros(len(idx))
    run = Z['run'][idx[0]]
    allk = np.where(Z['run'] == run)[0]
    full = Z['T_gauss028_mean']
    lut = {(int(Z['camcol'][j]), int(Z['field'][j])): full[j] for j in allk}
    ks = np.arange(-int(4 * sig) - 1, int(4 * sig) + 2)
    g = np.exp(-0.5 * (ks / sig) ** 2)
    for n, j in enumerate(idx):
        c, f = int(Z['camcol'][j]), int(Z['field'][j])
        v = np.array([lut.get((c, f + k), np.nan) for k in ks])
        m = np.isfinite(v)
        out[n] = np.sum(g[m] * v[m]) / np.sum(g[m])
    return out


def smoothing_test():
    import verify_precovery_fit as VF
    F = json.loads((DER / 'verify_precovery_fit.json').read_text())
    core = [r for r in F['cases']['g_psp_gauss028_OIIIunits']['runs'] if 12 <= r['dec_cross'] < 20.5]
    eff = {r['run']: max(r['eA_formal'], r['e_null_off'] or r['e_null_in'] or 0) for r in core}
    out = {}
    for est in ('psp', 'frames', 'framessub'):
        y = VF.Z[VF.est_key('g', est)].astype(float)
        row = {}
        rms0, mask0 = {}, {}
        for sig in (0, 0.75, 1.5, 2.5):
            num = den = chi = 0.0
            for r in core:
                idx = VF.select(r['run'])
                T = smoothed_template(idx, VF.Z['T_gauss028_mean'][idx], sig)
                fr = VF.fit(idx, y, [T])
                if sig == 0:
                    rms0[r['run']] = fr['rms']; mask0[r['run']] = fr['mask']
                w = 1 / eff[r['run']] ** 2
                num += w * fr['coef'][0] / VF.CONV_OIII_ONLY; den += w
                m = mask0[r['run']]          # same clipping mask and noise for every sigma -> comparable chi2
                chi += float(np.sum(fr['resid'][m] ** 2)) / rms0[r['run']] ** 2
            row[f'sigma{sig}'] = dict(A=num / den, e=1 / np.sqrt(den), chi2_total=chi)
        c0 = row['sigma0']['chi2_total']
        for v in row.values():
            v['dchi2_vs_sigma0'] = v['chi2_total'] - c0
        out[est] = row
        print('smoothing', est, {k: (round(v['A'], 2), round(v['dchi2_vs_sigma0'], 1)) for k, v in row.items()})
    return out


if __name__ == '__main__':
    res['template_time_smoothing_g'] = smoothing_test()
    (DER / 'verify_precovery_pspkernel.json').write_text(json.dumps(res, indent=1))
