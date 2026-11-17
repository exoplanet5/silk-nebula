"""spectra_desi_regspec.py -- 'regression spectrum' of the filament component and on/off stacks.

For each fibre i the residual r_i = flux - model is continuum-subtracted window by window (linear, ivar-weighted,
lines and the object's own redshifted lines masked) and converted to in-fibre surface brightness
c_i(lambda) [R/A] = r_i * ffrac_i * k_R(lambda) (ffrac = PSF_TO_FIBER_SPECFLUX, k_R from a 1.52" fibre).
Given an instrument P_i (the expected [OIII]5007 residual of fibre i in R, built WITHOUT fibre i's own data:
 'model' = petal forward-model prediction, 'nbr' = neighbour-averaged DESI [OIII] residual within 0.12 deg),
    S(lambda) = sum_i w_i P_i c_i / sum_i w_i P_i^2,   w_i = ivar_i / (ffrac_i k_R)^2
is the mean spectrum per unit [OIII] residual: integrating a line in S gives (line / [OIII]5007) in energy
units for the component spatially correlated with [OIII]; uncorrelated DIG/sky structure averages out.
Because DESI residuals are differential (fibre minus petal mean), S describes the filament's own spectrum
free of smooth foreground/background emission.

Outputs: data/derived/spectra_desi_regspec.npz (wave, S, eS for each instrument; stacks), and line fits in
         data/derived/spectra_desi_regspec_lines.json
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import LINES, WIN_WAVE, WINDOWS, C_KMS, GAL_LINES  # noqa: E402
import spectra_desi_measure as meas  # noqa: E402
from spectra_desi_analyze import rayleigh_per_unit  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
ALL_LINES = np.array(list(LINES.values()))


def continuum_subtract(R, W, z, excl=6.0, own=6.0):
    """Window-by-window linear continuum removal; returns residual-minus-continuum and modified weights."""
    C = np.zeros_like(R)
    Wout = W.copy()
    N = R.shape[0]
    ownmask = np.zeros(R.shape, bool)
    for gl in GAL_LINES:
        lo = gl * (1 + z)
        ownmask |= np.abs(WIN_WAVE[None, :] - lo[:, None]) < own
    Wout[ownmask] = 0
    for lo, hi in WINDOWS:
        s = (WIN_WAVE >= lo) & (WIN_WAVE <= hi)
        lam = WIN_WAVE[s]
        cm = np.ones(lam.size, bool)
        for l0 in ALL_LINES:
            cm &= np.abs(lam - l0) > excl
        r = R[:, s]
        w = Wout[:, s] * cm[None, :]
        x = (lam - lam.mean()) / 30.0
        for it in range(2):
            S0, S1, S2 = w.sum(1), (w * x).sum(1), (w * x * x).sum(1)
            T0, T1 = (w * r).sum(1), (w * r * x).sum(1)
            det = S0 * S2 - S1 ** 2
            bad = det <= 0
            det[bad] = 1
            a = (S2 * T0 - S1 * T1) / det
            b = (S0 * T1 - S1 * T0) / det
            cont = a[:, None] + b[:, None] * x[None, :]
            if it == 0:
                zz = (r - cont) * np.sqrt(w)
                w = w * (np.abs(zz) < 3)
        C[:, s] = r - cont
        Wout[np.ix_(bad, np.where(s)[0])] = 0
    return C, Wout


def regression_spectrum(C, W, P, okP):
    num = np.zeros(C.shape[1])
    den = np.zeros(C.shape[1])
    for i0 in range(0, C.shape[0], 5000):
        sl = slice(i0, i0 + 5000)
        p = np.where(okP[sl], P[sl], 0.0)[:, None]
        num += np.sum(W[sl] * p * C[sl], 0)
        den += np.sum(W[sl] * p * p, 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        S = num / den
        eS = 1 / np.sqrt(den)
    return S, eS


def stack(C, W, m):
    w = W[m]
    with np.errstate(invalid='ignore', divide='ignore'):
        s = np.sum(w * C[m], 0) / np.sum(w, 0)
        e = 1 / np.sqrt(np.sum(w, 0))
    return s, e


def fit_lines_1d(S, eS, cols, vgrid=np.arange(-150, 150.1, 1.0), sig_int=20.0, half=30.0, free_sigma=False):
    """Fit sum of Gaussians (shared velocity) to a 1-D spectrum; returns amplitudes (integrated), errors, v."""
    lam_c = np.array([c[0][0] for c in cols])
    sel = (WIN_WAVE >= lam_c.min() - half) & (WIN_WAVE <= lam_c.max() + half) & np.isfinite(S) & (eS > 0)
    lam, y, e = WIN_WAVE[sel], S[sel], eS[sel]
    if lam.size < 12:
        return dict(amp=[np.nan] * len(cols), err=[np.nan] * len(cols), v=np.nan, v_lo=np.nan, v_hi=np.nan,
                    chi2=np.nan, dof=0, red_chi2=np.nan)
    sigs = [[np.hypot(meas.lsf_sigma_A(l0), l0 * sig_int / C_KMS) for l0, _ in c] for c in cols]
    best = None
    chis = []
    for v in vgrid:
        G = np.zeros((len(cols) + 1, lam.size))
        for k, c in enumerate(cols):
            for (l0, amp), s in zip(c, sigs[k]):
                mu = l0 * (1 + v / C_KMS)
                G[k] += amp * np.exp(-0.5 * ((lam - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
        G[-1] = 1.0                                   # residual constant offset
        A = (G / e ** 2) @ G.T
        bvec = (G / e ** 2) @ y
        cov = np.linalg.pinv(A)
        p = cov @ bvec
        chi2 = np.sum(((y - p @ G) / e) ** 2)
        chis.append(chi2)
        if best is None or chi2 < best[0]:
            best = (chi2, v, p, cov)
    chis = np.array(chis)
    chi2, v, p, cov = best
    dof = lam.size - len(cols) - 2
    scale = max(chi2 / dof, 1.0)
    within = vgrid[chis <= chi2 + scale]
    return dict(amp=p[:-1].tolist(), err=(np.sqrt(np.diag(cov))[:-1] * np.sqrt(scale)).tolist(), v=float(v),
                v_lo=float(within.min()), v_hi=float(within.max()), chi2=float(chi2), dof=int(dof),
                red_chi2=float(chi2 / dof))


def fit_width(S, eS, l0s, ratios, v0, half=25.0):
    """Fit a single shared Gaussian sigma (km/s, total) to the [OIII] doublet in S (for line-width limits)."""
    sel = (WIN_WAVE >= min(l0s) - half) & (WIN_WAVE <= max(l0s) + half) & np.isfinite(S)
    lam, y, e = WIN_WAVE[sel], S[sel], eS[sel]
    grid_s = np.arange(25, 90.1, 0.5)
    grid_v = np.arange(v0 - 30, v0 + 30.1, 1.0)
    chi = np.full((grid_s.size, grid_v.size), np.inf)
    for i, sk in enumerate(grid_s):
        for j, v in enumerate(grid_v):
            G = np.zeros((2, lam.size))
            for l0, a in zip(l0s, ratios):
                s = l0 * sk / C_KMS
                G[0] += a * np.exp(-0.5 * ((lam - l0 * (1 + v / C_KMS)) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
            G[1] = 1
            A = (G / e ** 2) @ G.T
            p = np.linalg.solve(A, (G / e ** 2) @ y)
            chi[i, j] = np.sum(((y - p @ G) / e) ** 2)
    i, j = np.unravel_index(np.argmin(chi), chi.shape)
    dof = lam.size - 4
    sc = max(chi[i, j] / dof, 1)
    prof = chi.min(1)
    ok = grid_s[prof <= chi[i, j] + sc]
    return dict(sigma_tot_kms=float(grid_s[i]), lo=float(ok.min()), hi=float(ok.max()), v=float(grid_v[j]),
                red_chi2=float(chi[i, j] / dof))


GROUPS_1D = {
    'OIII5007': [[(LINES['OIII5007'], 1.0), (LINES['OIII4959'], 1 / 2.98)]],
    'OIII4959': [[(LINES['OIII4959'], 1.0)]],
    'OIII5007_single': [[(LINES['OIII5007'], 1.0)]],
    'Hb': [[(LINES['Hb'], 1.0)]],
    'Hg': [[(LINES['Hg'], 1.0)]],
    'HaNII': [[(LINES['Ha'], 1.0)], [(LINES['NII6583'], 1.0), (LINES['NII6548'], 1 / 3.05)]],
    'SII': [[(LINES['SII6716'], 1.0)], [(LINES['SII6731'], 1.0)]],
    'OI6300': [[(LINES['OI6300'], 1.0), (LINES['OI6364'], 1 / 3.05)]],
    'HeII4686': [[(LINES['HeII4686'], 1.0)]],
    'HeI5876': [[(LINES['HeI5876'], 1.0)]],
    'NI5200': [[(LINES['NI5198'], 1.0), (LINES['NI5200'], 0.7)]],
    'OII3727': [[(LINES['OII3726'], 1.0)], [(LINES['OII3729'], 1.0)]],
    'NeIII3869': [[(LINES['NeIII3869'], 1.0)]],
    'OIII4363': [[(LINES['OIII4363'], 1.0)]],
    'NII5755': [[(LINES['NII5755'], 1.0)]],
    'ArIII7136': [[(LINES['ArIII7136'], 1.0)]],
    'OII7325': [[(LINES['OII7320'], 1.0), (LINES['OII7330'], 0.8)]],
    'SIII9069': [[(LINES['SIII9069'], 1.0)]],
    'sky5577': [[(LINES['sky5577'], 1.0)]],
}


def main():
    ids, F, I, M, Sg, files = meas.load_windows()
    arr = np.load(D / 'spectra_desi_model_arrays.npz', allow_pickle=False)
    pos = {s: k for k, s in enumerate(arr['sparcl_id'])}
    sel = np.array([i in pos for i in ids])
    ids, F, I, M = ids[sel], F[sel], I[sel], M[sel]
    k = np.array([pos[i] for i in ids])
    cat = {r['sparcl_id']: r for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv'))}
    z = np.array([float(cat[i]['redshift']) for i in ids])
    ffrac = np.clip(np.array([float(cat[i]['mean_psf_to_fiber_specflux']) for i in ids]), 0.3, 1.0)
    C, W = continuum_subtract(F - M, I, z)
    del F, M
    kR = rayleigh_per_unit(WIN_WAVE)
    conv = ffrac[:, None] * kR[None, :]
    C = C * conv
    W = W / conv ** 2
    out = {'wave': WIN_WAVE}
    fits = {}
    for inst in ('model', 'nbr'):
        P = arr['pred'][k] if inst == 'model' else arr['nbO'][k]
        okP = np.isfinite(P)
        S, eS = regression_spectrum(C, W, P, okP)
        out[f'S_{inst}'], out[f'eS_{inst}'] = S, eS
        fits[inst] = {g: fit_lines_1d(S, eS, cols) for g, cols in GROUPS_1D.items()}
        fits[inst]['OIII_width'] = fit_width(S, eS, [LINES['OIII5007'], LINES['OIII4959']], [1, 1 / 2.98],
                                             fits[inst]['OIII5007']['v'])
        fits[inst]['n_fibres'] = int(okP.sum())
    # on/off stacks by model prediction
    pred = arr['pred'][k]
    nb = arr['nbO'][k]
    on = np.isfinite(pred) & (pred > np.nanpercentile(pred, 90))
    neg = np.isfinite(pred) & (pred < np.nanpercentile(pred, 10))
    far = np.array([float(cat[i]['d_any']) > 3 for i in ids])
    for lab, m in (('on', on), ('neg', neg), ('far', far)):
        if m.sum() < 20:
            continue
        s, e = stack(C, W, m)
        out[f'stack_{lab}'], out[f'estack_{lab}'] = s, e
        out[f'n_{lab}'] = int(m.sum())
        fits.setdefault('stacks', {})[lab] = {g: fit_lines_1d(s, e, cols) for g, cols in GROUPS_1D.items()}
        fits['stacks'][lab]['n'] = int(m.sum())
        fits['stacks'][lab]['mean_pred_R'] = float(np.nanmean(pred[m]))
    np.savez_compressed(D / 'spectra_desi_regspec.npz', **out)
    (D / 'spectra_desi_regspec_lines.json').write_text(json.dumps(fits, indent=1))
    for inst in ('model', 'nbr'):
        f = fits[inst]
        print(inst, 'n', f['n_fibres'])
        for g in GROUPS_1D:
            print(f"  {g:12s} amp {np.round(f[g]['amp'], 4)} err {np.round(f[g]['err'], 4)} v {f[g]['v']:.0f} "
                  f"[{f[g]['v_lo']:.0f},{f[g]['v_hi']:.0f}] rchi2 {f[g]['red_chi2']:.2f}")
        print('  width', f['OIII_width'])


if __name__ == '__main__':
    main()
