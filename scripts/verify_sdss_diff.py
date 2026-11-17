"""verify_sdss stage 2: per-fibre differential line fits against a same-group OFF reference.

For a chosen OFF definition, every fibre's TOTAL sky (R/A, as calibrated) minus the clipped mean of its group's OFF
fibres is fitted with fixed-velocity Gaussians (LSF sigma from wdisp) + linear baseline + eps_c * (group continuum
template, incl. solar/stellar absorption lines) [+ eps_OH * group OH template in the red].  Output is one row per
fibre with line amplitudes (R, as calibrated) and formal errors; the template coefficients guard against ON/OFF
continuum differences filling in or digging out Hbeta/Halpha absorption, and against spatial OH structure.

Usage: python -I verify_sdss_diff.py <tag> <off_lo> <off_hi> <pairing: sp|plate> <exclude_patch:0|1> [v_neb] [use_templates 1|0]
Writes data/derived/verify_sdss_diff_<tag>.npz
"""
import sys, csv
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, VAC, C, load_windows, flam_to_R_per_A, fit_lines, gauss_unit

LINESETS = {   # name: (lo, hi, lines, ties, use_oh)
    'OII': (3700, 3758, ['OII3726', 'OII3729'], {}, False),
    'B4363': (4318, 4398, ['Hgam', 'OIII4363'], {}, False),      # Hg I 4358 added as a topocentric column
    'HeII': (4660, 4715, ['HeII4686'], {}, False),
    'Hb': (4835, 4890, ['Hb'], {}, False),
    'OIII': (4935, 5032, ['OIII4959', 'OIII5007'], {'OIII4959': ('OIII5007', 1 / 3.01)}, False),
    'HaNII': (6522, 6618, ['NII6548', 'Ha', 'NII6583'], {'NII6548': ('NII6583', 1 / 3.06)}, True),
    'SII': (6700, 6752, ['SII6716', 'SII6731'], {}, True),
}


def clipped_mean(X, V, nsig=4.0):
    good = np.isfinite(X) & np.isfinite(V)
    Xm = np.where(good, X, np.nan)
    with np.errstate(all='ignore'):
        med = np.nanmedian(Xm, axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(Xm - med), axis=0)
    ok = good & (np.abs(Xm - med) < nsig * np.maximum(mad, 1e-6))
    n = ok.sum(0)
    mean = np.where(n > 0, np.sum(np.where(ok, X, 0), 0) / np.maximum(n, 1), np.nan)
    var = np.where(n > 0, np.sum(np.where(ok, V, 0), 0) / np.maximum(n, 1) ** 2, np.inf)
    return mean, var


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def main(tag, off_lo, off_hi, pairing, excl_patch, v_neb=13.0, use_tmpl=1):
    D = load_windows()
    lam = D['lam']
    F = read_csv(ROOT / 'data/derived/verify_sdss_fibres.csv')
    G = read_csv(ROOT / 'data/derived/verify_sdss_groups.csv')
    T = np.load(ROOT / 'data/derived/verify_sdss_templates.npz')
    tc, toh = T['cont'], T['oh']
    gid = np.array([int(r['gid']) for r in F])
    dmin = np.array([float(r['dmin']) for r in F])
    inp = np.array([int(r['in_patch30']) for r in F]).astype(bool)
    diam = np.array([float(r['diam']) for r in F])
    plate_key = np.array([f"{r['release']}-{r['plate']}-{r['mjd']}" for r in F])
    TOT = D['tot'].astype(float); IV = D['ivar'].astype(float); WD = D['wdisp'].astype(float)
    helio = D['m_helio_rv'].astype(float)
    nf = len(F)
    out = {ls: np.full((nf, 2), np.nan) for ls in []}
    res = {}
    for ls, (lo, hi, lines, ties, use_oh) in LINESETS.items():
        for L in lines:
            if L not in ties:
                res[f'{ls}:{L}'] = np.full(nf, np.nan); res[f'{ls}:{L}:err'] = np.full(nf, np.nan)
        res[f'{ls}:chi2nu'] = np.full(nf, np.nan)
        res[f'{ls}:eps_c'] = np.full(nf, np.nan)
    res['B4363:Hg4358'] = np.full(nf, np.nan); res['B4363:Hg4358:err'] = np.full(nf, np.nan)
    noff = np.zeros(nf, int)
    # grouping for the OFF reference
    if pairing == 'sp':
        ref_key = gid
    else:
        u = {k: i for i, k in enumerate(sorted(set(plate_key)))}
        ref_key = np.array([u[k] for k in plate_key])
    off_mask = (dmin > off_lo) & (dmin < off_hi)
    if excl_patch:
        off_mask &= ~inp
    for rk in np.unique(ref_key):
        sel = np.where(ref_key == rk)[0]
        offs = sel[off_mask[sel]]
        if len(offs) < 3:
            continue
        X = flam_to_R_per_A(TOT[sel], lam, diam[sel][:, None])
        V = np.where(IV[sel] > 0, flam_to_R_per_A(1.0, lam, diam[sel][:, None]) ** 2 / np.where(IV[sel] > 0, IV[sel], 1), np.inf)
        Xo = X[np.isin(sel, offs)]; Vo = V[np.isin(sel, offs)]
        Mo, MVo = clipped_mean(Xo, Vo)
        noff[sel] = len(offs)
        for j, i in enumerate(sel):
            dX = X[j] - Mo
            dV = V[j] + MVo
            g = gid[i]
            sig = WD[i] * lam * np.log(10) * 1e-4
            for ls, (lo, hi, lines, ties, use_oh) in LINESETS.items():
                m = (lam > lo) & (lam < hi) & np.isfinite(dX) & np.isfinite(dV) & np.isfinite(sig) & np.isfinite(tc[g])
                if m.sum() < 25:
                    continue
                lm = lam[m]
                s_med = float(np.nanmedian(sig[m]))
                ex = {'eps_c': tc[g][m] / max(np.nanmedian(np.abs(tc[g][m])), 1e-6)} if use_tmpl else {}
                if use_tmpl and use_oh and np.any(toh[g][m] != 0):
                    ex['eps_oh'] = toh[g][m] / max(np.max(np.abs(toh[g][m])), 1e-6)
                if ls == 'B4363':
                    # street-lamp Hg I 4358 is topocentric: lam0 (1 + v_corr/c), v_corr = -HELIO_RV
                    mu = VAC['Hg4358'] * (1 - helio[i] / C)
                    ex['Hg4358'] = gauss_unit(lm, mu, s_med)
                    ex['Hg4358_d'] = gauss_unit(lm, mu, s_med) * (lm - mu) / s_med   # small shift / LSF mismatch
                f = fit_lines(lm, dX[m], 1 / dV[m], lines, v_neb, (lambda x, s=s_med: s), npoly=2, ties=ties, extra=ex)
                for L in lines:
                    if L not in ties:
                        res[f'{ls}:{L}'][i] = f[L][0]; res[f'{ls}:{L}:err'][i] = f[L][1]
                res[f'{ls}:chi2nu'][i] = f['_chi2'] / max(f['_dof'], 1)
                if ls == 'B4363':
                    res['B4363:Hg4358'][i] = f['Hg4358'][0]; res['B4363:Hg4358:err'][i] = f['Hg4358'][1]
                res[f'{ls}:eps_c'][i] = f['eps_c'][0] if use_tmpl else np.nan
    np.savez_compressed(ROOT / f'data/derived/verify_sdss_diff_{tag}.npz', noff=noff, off_lo=off_lo, off_hi=off_hi,
                        pairing=pairing, excl_patch=excl_patch, v_neb=v_neb, **{k.replace(':', '__'): v for k, v in res.items()})
    ok = np.isfinite(res['OIII:OIII5007'])
    print(tag, 'fibres fitted', ok.sum(), 'of', nf)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]), a[3], int(a[4]), float(a[5]) if len(a) > 5 else 13.0, int(a[6]) if len(a) > 6 else 1)
