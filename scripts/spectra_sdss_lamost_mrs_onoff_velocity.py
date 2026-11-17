"""LAMOST DR12 MRS: [O III] 5007 velocity from ON-minus-OFF residual stacks per plate (topic spectra_sdss_lamost).

The ON stack alone is (strand at the stars) - (super-sky [O III]); the super-sky term is common to ON and OFF
stars of the same plate, so ON - OFF isolates the strand.  Gaussian with free centre, width fixed to the
instrumental (R = 7500, sigma = 17.0 km/s) convolved with 19 km/s intrinsic, linear baseline, 5005-5012 A.
Also refits the ON and OFF stacks the same way for comparison.

Input : data/derived/spectra_sdss_lamost_mrs_resid_stacks.json (from spectra_sdss_lamost_lamost_resid.py)
Output: data/derived/spectra_sdss_lamost_mrs_onoff_velocity.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
L0 = VAC['OIII5007']
SIG = np.hypot(C_KMS / 7500 / 2.3548, 19.0) / C_KMS * L0      # A


def fit(g, m, e):
    w = (g > 5005.0) & (g < 5012.0) & np.isfinite(m) & (e > 0)
    f = lambda x, a, v, b0, b1: a * np.exp(-0.5 * ((x - L0 * (1 + v / C_KMS)) / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG) + b0 + b1 * (x - L0)
    best = None
    for v0 in (-20.0, 10.0, 40.0):
        for a0 in (-1.0, 1.0):
            try:
                p, c = curve_fit(f, g[w], m[w], p0=[a0, v0, 0, 0], sigma=e[w], absolute_sigma=True, maxfev=20000,
                                 bounds=([-50, -80, -10, -10], [50, 100, 10, 10]))
            except (RuntimeError, ValueError):
                continue
            chi = float(np.sum(((m[w] - f(g[w], *p)) / e[w]) ** 2))
            if best is None or chi < best[2]:
                best = (p, c, chi)
    p, c, chi = best
    dof = int(w.sum()) - 4
    sc = np.sqrt(max(chi / dof, 1))
    return dict(F_R=float(p[0]), F_err=float(np.sqrt(c[0, 0]) * sc), v_helio=float(p[1]), v_err=float(np.sqrt(c[1, 1]) * sc),
                chi2_red=float(chi / dof), snr=float(abs(p[0]) / (np.sqrt(c[0, 0]) * sc)))


def main():
    st = json.loads((ROOT / 'data/derived/spectra_sdss_lamost_mrs_resid_stacks.json').read_text())
    out = {'sigma_fixed_A': float(SIG), 'window_A': [5005.0, 5012.0]}
    tags = [''] + sorted({k.split('_', 2)[2] for k in st if k.startswith('OIII_on_')})
    for t in tags:
        kon, koff = ('OIII_on' + ('_' + t if t else ''), 'OIII_off' + ('_' + t if t else ''))
        if kon not in st or koff not in st:
            continue
        g = np.array(st[kon]['grid'])
        mon, eon = np.array(st[kon]['mean']), np.array(st[kon]['err'])
        mof, eof = np.array(st[koff]['mean']), np.array(st[koff]['err'])
        r = {'n_on': st[kon]['n'], 'n_off': st[koff]['n'],
             'on': fit(g, mon, eon), 'off': fit(g, mof, eof), 'on_minus_off': fit(g, mon - mof, np.hypot(eon, eof))}
        out[t or 'all_plates'] = r
        print(f"{t or 'all':22s} ON {r['on']['F_R']:+.2f}+-{r['on']['F_err']:.2f} v {r['on']['v_helio']:+.1f}+-{r['on']['v_err']:.1f} | "
              f"OFF {r['off']['F_R']:+.2f}+-{r['off']['F_err']:.2f} v {r['off']['v_helio']:+.1f}+-{r['off']['v_err']:.1f} | "
              f"ON-OFF {r['on_minus_off']['F_R']:+.2f}+-{r['on_minus_off']['F_err']:.2f} v {r['on_minus_off']['v_helio']:+.1f}+-{r['on_minus_off']['v_err']:.1f}")
    # the paper's plate MN074209N160100N01 (2024-02-01): 77 public star spectra with RVs; their 5007 residual is
    # I(star position) - <I> over the ~165 nebular fibres per spectrograph -> a dip = mean nebular-fibre [O III]
    pm = ROOT / 'data/derived/spectra_sdss_lamost_mrsmn_resid_stacks.json'
    if pm.exists():
        sm = json.loads(pm.read_text())
        for key in [k for k in sm if k.startswith('OIII_on')]:
            g, m, e = (np.array(sm[key][c]) for c in ('grid', 'mean', 'err'))
            r = fit(g, m, e)
            r['n'] = sm[key]['n']
            out['MN074209N160100N01_' + key] = r
            print(f"MN plate {key:32s} n {r['n']:3d}  F {r['F_R']:+.2f}+-{r['F_err']:.2f} R  v {r['v_helio']:+.1f}+-{r['v_err']:.1f}  snr {r['snr']:.1f}")
    (ROOT / 'data/derived/spectra_sdss_lamost_mrs_onoff_velocity.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
