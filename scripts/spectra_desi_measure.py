"""spectra_desi_measure.py -- Gaussian line fits to DESI sky-subtracted residual spectra (flux - redrock model).

For every saved spectrum and every line group:
  residual r = flux - model (SPARCL 'model' = Redrock best fit, already convolved with the resolution);
  local linear continuum fitted (ivar-weighted, 3-sigma clipped once) to the group window excluding
  +-6 A around every catalogued line and +-6 A around the object's own redshifted lines;
  line model = sum of Gaussians (shared velocity, fixed doublet ratios), sigma^2 = sigma_LSF^2 + sigma_int^2,
  sigma_LSF from DESI 'wave_sigma' (pilot median), sigma_int = 20 km/s (LAMOST [OIII] FWHM 43-49 km/s);
  amplitudes by weighted linear least squares on a velocity grid -150..+150 km/s (2 km/s):
    'fix' = at fixed barycentric velocity V_FIX (+25 km/s for [OIII]/He II, +17 km/s for the rest, LAMOST),
    'free' = at the chi^2 minimum, velocity error from the local chi^2 curvature (dchi2 = 1).
Null groups at line-free wavelengths calibrate the error scale and false-positive rate.

Output: data/derived/spectra_desi_fits.npz (structured arrays per group) and spectra_desi_fits.csv (flat).
Usage : python -I scripts/spectra_desi_measure.py [tags...]   (default: all *_NNNNN.npz in the windows dir)
"""
import csv
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import LINES, WIN_WAVE, C_KMS, GAL_LINES  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
WDIR = ROOT / 'data/derived/spectra_desi_windows'
VGRID = np.arange(-150.0, 150.1, 2.0)
SIG_INT = 20.0
L = LINES

# group: (columns, v_fix, window half-range). Column = list of (lambda_vac, relative amplitude); first entry
# of each column is the reported line.
GROUPS = {
    'OIII5007': ([[(L['OIII5007'], 1.0), (L['OIII4959'], 1 / 2.98)]], 25.0),
    'OIII5007only': ([[(L['OIII5007'], 1.0)]], 25.0),
    'OIII4959only': ([[(L['OIII4959'], 1.0)]], 25.0),
    'Hb': ([[(L['Hb'], 1.0)]], 17.0),
    'HaNII': ([[(L['Ha'], 1.0)], [(L['NII6583'], 1.0), (L['NII6548'], 1 / 3.05)]], 17.0),
    'SII': ([[(L['SII6716'], 1.0)], [(L['SII6731'], 1.0)]], 17.0),
    'OI6300': ([[(L['OI6300'], 1.0), (L['OI6364'], 1 / 3.05)]], 17.0),
    'HeII4686': ([[(L['HeII4686'], 1.0)]], 25.0),
    'NI5200': ([[(L['NI5198'], 1.0), (L['NI5200'], 0.7)]], 17.0),
    'OII3727': ([[(L['OII3726'], 1.0)], [(L['OII3729'], 1.0)]], 17.0),
    'NeIII3869': ([[(L['NeIII3869'], 1.0)]], 25.0),
    'Hg': ([[(L['Hg'], 1.0)]], 17.0),
    'OIII4363': ([[(L['OIII4363'], 1.0)]], 25.0),
    'HeI5876': ([[(L['HeI5876'], 1.0)]], 17.0),
    'NII5755': ([[(L['NII5755'], 1.0)]], 17.0),
    'ArIII7136': ([[(L['ArIII7136'], 1.0)]], 17.0),
    'SIII9069': ([[(L['SIII9069'], 1.0)]], 17.0),
    'sky5577': ([[(L['sky5577'], 1.0)]], 0.0),
    'null4990': ([[(4990.0, 1.0)]], 0.0),
    'null5030': ([[(5030.0, 1.0)]], 0.0),
    'null6620': ([[(6615.0, 1.0)]], 0.0),
    'null4880': ([[(4885.0, 1.0)]], 0.0),
}
ALL_NEB = np.array([v for k, v in L.items()] + [5577.34 * 1.000278])
SKY_OH_MASK = []   # strong OH lines are handled by ivar


def lsf_sigma_A(lam):
    """DESI LSF sigma (A): piecewise-linear in wavelength through the pilot medians (spectra_desi_lsf.json)."""
    p = ROOT / 'data/derived/spectra_desi_lsf.json'
    if p.exists():
        import json
        t = json.loads(p.read_text())
        return float(np.interp(lam, t['wave'], t['sigma_A']))
    return 0.70  # fallback (pilot value near 5008 A)


def load_windows(tags=None):
    files = sorted(WDIR.glob('*_[0-9][0-9][0-9][0-9][0-9].npz'))
    if tags:
        files = [f for f in files if f.name.rsplit('_', 1)[0] in tags]
    ids, F, I, M, S = [], [], [], [], []
    for f in files:
        with np.load(f, allow_pickle=False) as d:
            ids.append(d['sparcl_id'])
            F.append(d['flux'])
            I.append(d['ivar'])
            M.append(d['model'])
            if 'wave_sigma' in d.files:
                S.append(d['wave_sigma'])
    ids = np.concatenate(ids)
    F, I, M = np.concatenate(F), np.concatenate(I), np.concatenate(M)
    u, first = np.unique(ids, return_index=True)        # drop duplicates (re-downloads)
    first = np.sort(first)
    S = np.concatenate(S) if S else None
    return ids[first], F[first], I[first], M[first], S, files


def fit_group(R, W, zobj, cols, v_fix, own_lines_excl=6.0, half=34.0, excl=6.0):
    """R, W: residual and ivar [N, P] on WIN_WAVE. Returns dict of arrays."""
    lam_c = np.array([c[0][0] for c in cols])
    lo, hi = lam_c.min() - half, lam_c.max() + half
    sel = (WIN_WAVE >= lo) & (WIN_WAVE <= hi)
    lam = WIN_WAVE[sel].astype(float)
    r = R[:, sel].astype(float)
    w = W[:, sel].astype(float)
    w[~np.isfinite(r)] = 0
    r[~np.isfinite(r)] = 0
    N, P = r.shape
    # continuum mask: catalogued lines (rest, any velocity within the grid) + object's own lines
    cmask = np.ones(P, bool)
    for l0 in ALL_NEB:
        cmask &= np.abs(lam - l0) > excl
    own = np.zeros((N, P), bool)
    for gl in GAL_LINES:
        lo_obs = gl * (1 + zobj)
        own |= np.abs(lam[None, :] - lo_obs[:, None]) < own_lines_excl
    # object-line flag: own line within 8 A of any of the group's line centres
    flag_own = np.zeros(N, bool)
    for gl in GAL_LINES:
        lo_obs = gl * (1 + zobj)
        for lc in lam_c:
            flag_own |= np.abs(lo_obs - lc) < 8.0
    wc = w * cmask[None, :] * (~own)
    x = (lam - lam.mean()) / 30.0
    # weighted linear continuum, one 3-sigma clip
    for it in range(2):
        S0, S1, S2 = wc.sum(1), (wc * x).sum(1), (wc * x * x).sum(1)
        T0, T1 = (wc * r).sum(1), (wc * r * x).sum(1)
        det = S0 * S2 - S1 ** 2
        det[det == 0] = np.inf
        a = (S2 * T0 - S1 * T1) / det
        b = (S0 * T1 - S1 * T0) / det
        cont = a[:, None] + b[:, None] * x[None, :]
        if it == 0:
            z = (r - cont) * np.sqrt(wc)
            wc = wc * (np.abs(z) < 3)
    rr = r - cont
    sigs = []
    for c in cols:
        sigs.append([np.hypot(lsf_sigma_A(l0), l0 * SIG_INT / C_KMS) for l0, _ in c])

    def design(v):
        G = np.zeros((len(cols), P))
        for k, c in enumerate(cols):
            for (l0, amp), s in zip(c, sigs[k]):
                mu = l0 * (1 + v / C_KMS)
                G[k] += amp * np.exp(-0.5 * ((lam - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
        return G

    def solve(v):
        G = design(v)                           # [K, P]
        A = np.einsum('kp,np,jp->nkj', G, w, G)  # [N, K, K]
        bvec = np.einsum('kp,np,np->nk', G, w, rr)
        K = len(cols)
        good = np.linalg.det(A) > 0
        A[~good] = np.eye(K)
        Ainv = np.linalg.inv(A)
        amp = np.einsum('nkj,nj->nk', Ainv, bvec)
        dchi = np.einsum('nk,nk->n', amp, bvec)  # chi2 reduction
        err = np.sqrt(np.clip(np.einsum('nkk->nk', Ainv), 0, None))
        amp[~good] = np.nan
        err[~good] = np.nan
        dchi[~good] = np.nan
        return amp, err, dchi

    amp_fix, err_fix, dchi_fix = solve(v_fix)
    dchis, amps, errs = [], [], []
    for v in VGRID:
        a_, e_, d_ = solve(v)
        dchis.append(d_)
        amps.append(a_)
        errs.append(e_)
    dchis = np.array(dchis)                     # [V, N]
    dch = np.where(np.isfinite(dchis), dchis, -np.inf)
    ib = np.argmax(dch, 0)
    n_idx = np.arange(N)
    amps = np.array(amps)
    errs = np.array(errs)
    amp_free = amps[ib, n_idx]
    err_free = errs[ib, n_idx]
    v_free = VGRID[ib].astype(float)
    # parabolic refinement and curvature error (chi2 = -dchi)
    ibc = np.clip(ib, 1, len(VGRID) - 2)
    ym, y0, yp = -dchis[ibc - 1, n_idx], -dchis[ibc, n_idx], -dchis[ibc + 1, n_idx]
    h = VGRID[1] - VGRID[0]
    curv = (yp - 2 * y0 + ym) / h ** 2
    with np.errstate(invalid='ignore', divide='ignore'):
        shift = np.where(curv > 0, -0.5 * h * (yp - ym) / (yp - 2 * y0 + ym), 0.0)
        ev = np.where(curv > 0, np.sqrt(2.0 / curv), np.nan)
    edge = (ib == 0) | (ib == len(VGRID) - 1)
    v_free = np.where(edge, v_free, VGRID[ibc] + np.clip(shift, -h, h))
    ev[edge] = np.nan
    npix_good = (w > 0).sum(1)
    return dict(amp_fix=amp_fix, err_fix=err_fix, dchi_fix=dchi_fix, amp_free=amp_free, err_free=err_free,
                dchi_free=dch[ib, n_idx], v_free=v_free, ev_free=ev, v_edge=edge, flag_own=flag_own,
                npix_good=npix_good, cont=a)


def run(tags=None, out_stem='spectra_desi_fits'):
    ids, F, I, M, S, files = load_windows(tags)
    cat = {r['sparcl_id']: r for r in csv.DictReader(open(ROOT / 'data/derived/spectra_desi_catalog.csv'))}
    keep = np.array([i in cat for i in ids])
    ids, F, I, M = ids[keep], F[keep], I[keep], M[keep]
    z = np.array([float(cat[i]['redshift']) for i in ids])
    R = F - M
    print(f'{len(ids)} spectra from {len(files)} files')
    res = {'sparcl_id': ids}
    for g, (cols, vfix) in GROUPS.items():
        out = {}
        for i0 in range(0, len(ids), 4000):
            o = fit_group(R[i0:i0 + 4000], I[i0:i0 + 4000], z[i0:i0 + 4000], cols, vfix)
            for k, v in o.items():
                out.setdefault(k, []).append(v)
        for k, v in out.items():
            res[f'{g}__{k}'] = np.concatenate(v)
        print(g, 'done', flush=True)
    np.savez_compressed(ROOT / f'data/derived/{out_stem}.npz', **res)
    # flat csv (main quantities)
    keys = ['sparcl_id'] + [k for k in res if k != 'sparcl_id' and res[k].ndim == 1] + \
           [f'{g}__amp{j}_{m}' for g, (cols, _) in GROUPS.items() for j in range(len(cols)) for m in ('fix', 'free')]
    with open(ROOT / f'data/derived/{out_stem}.csv', 'w', newline='') as fh:
        wr = csv.writer(fh)
        hdr = ['sparcl_id']
        cols_out = []
        for g, (cols, _) in GROUPS.items():
            for j in range(len(cols)):
                for m in ('fix', 'free'):
                    hdr += [f'{g}_c{j}_amp_{m}', f'{g}_c{j}_err_{m}']
                    cols_out += [(f'{g}__amp_{m}', f'{g}__err_{m}', j)]
            for k in ('dchi_fix', 'dchi_free', 'v_free', 'ev_free', 'v_edge', 'flag_own', 'npix_good'):
                hdr.append(f'{g}_{k}')
                cols_out.append((f'{g}__{k}', None, None))
        wr.writerow(hdr)
        for n in range(len(ids)):
            row = [ids[n]]
            for a, e, j in cols_out:
                if j is None:
                    v = res[a][n]
                    row.append(f'{v:.4g}' if isinstance(v, (float, np.floating)) else int(v))
                else:
                    row += [f'{res[a][n, j]:.4g}', f'{res[e][n, j]:.4g}']
            wr.writerow(row)
    print('wrote', out_stem)


if __name__ == '__main__':
    run(sys.argv[1:] or None)
