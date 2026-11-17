"""verify_desi_linefit.py -- independent per-fibre line fitter for DESI DR1 sky-subtraction residuals (wave-2, desi).

Differences from wave 1 (spectra_desi_measure.py), on purpose:
  * continuum and lines fitted SIMULTANEOUSLY (quadratic in a +-30 A window + Gaussians), so continuum errors
    propagate; wave 1 fitted a linear continuum first;
  * night-sky lines masked/modelled: Hg I 4047/4358/5461/5770/5791, Na D, [O I]5577 placed at the topocentric
    wavelength shifted into the barycentric frame by the coadd's EFFTIME-weighted v_corr;
  * LSF sigma from my own median of the SPARCL wave_sigma arrays, plus 0.8 A pixel (Delta^2/12) and 12 km/s intrinsic;
  * fixed nebular velocity V_NEB (km/s, barycentric) passed on the command line (default +14), and a free-velocity
    scan for the [O III] doublet only;
  * nulls: fake lines at line-free wavelengths, including six around [O III]4363 (avoiding Hg 4358 and Hgamma).
Group '4363hg' fits [O III]4363 at V_NEB together with the Hg I 4358 sky line at v_corr (topocentric).

Output: data/derived/verify_desi_linefluxes.npz  (amp/err [N,K] per group, ids, plus o3 free-velocity scan)
Usage : ~/.venvs/astro313/bin/python -I scripts/verify_desi_linefit.py [V_NEB]
"""
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (D, C, WAVE, L, SKY_LINES, NEB_LINES, GAL_REST, lsf_sigma, load_spectra,  # noqa
                                load_meta, Exposures, Coadds)

SIG_INT = 12.0
PIX = 0.8


def groups(v_neb):
    n = lambda k: ('neb', L[k])   # noqa: E731
    s = lambda k: ('sky', L[k])   # noqa: E731
    f = lambda lam: ('neb', lam)  # noqa: E731
    G = {
        'o3': [[(n('OIII5007'), 1.0), (n('OIII4959'), 1 / 2.98)]],
        'o3_5007': [[(n('OIII5007'), 1.0)]],
        'o3_4959': [[(n('OIII4959'), 1.0)]],
        'hb': [[(n('Hb'), 1.0)]],
        'hg': [[(n('Hg'), 1.0)]],
        'hd': [[(n('Hd'), 1.0)]],
        '4363': [[(n('OIII4363'), 1.0)]],
        '4363hg': [[(n('OIII4363'), 1.0)], [(s('HgI4358'), 1.0)]],
        'hg4358': [[(s('HgI4358'), 1.0)]],
        'he2': [[(n('HeII4686'), 1.0)]],
        'o2': [[(n('OII3726'), 1.0)], [(n('OII3729'), 1.0)]],
        'o2sum': [[(n('OII3726'), 1.0), (n('OII3729'), 1.3)]],
        'ne3': [[(n('NeIII3869'), 1.0)]],
        'han2': [[(n('Ha'), 1.0)], [(n('NII6583'), 1.0), (n('NII6548'), 1 / 3.05)]],
        's2': [[(n('SII6716'), 1.0)], [(n('SII6731'), 1.0)]],
        'he1': [[(n('HeI5876'), 1.0)]],
        'n2_5755': [[(n('NII5755'), 1.0)]],
        'ar3': [[(n('ArIII7136'), 1.0)]],
        'o2_7325': [[(n('OII7320'), 1.0), (n('OII7330'), 0.8)]],
        'o2ca2': [[(n('OII7320'), 1.0), (n('OII7330'), 0.8)], [(n('CaII7324'), 1.0), (n('CaII7291'), 1.5)]],
        'ca2_7291': [[(n('CaII7291'), 1.0)]],
        'fe2_7155': [[(n('FeII7155'), 1.0)]],
        'fe2_4359': [[(n('FeII4359'), 1.0)]],
        'hg_fe2_4363': [[(n('OIII4363'), 1.0)], [(s('HgI4358'), 1.0)], [(n('FeII4359'), 1.0)]],
        'n1': [[(n('NI5198'), 1.0), (n('NI5200b'), 0.7)]],
        's3': [[(n('SIII9069'), 1.0)]],
        'sky5577': [[(s('sky5577'), 1.0)]],
    }
    for lam in (4885.0, 4990.0, 5030.0, 6615.0, 4700.0, 3760.0, 4354.0, 4369.6, 4374.6, 4385.0, 4395.0, 4330.0,
                7100.0, 7170.0, 7305.0, 7350.0, 5215.0):
        G[f'null{lam:.0f}' if lam not in (4369.6, 4374.6) else f'null{lam:.1f}'] = [[(f(lam), 1.0)]]
    return G


def fit_group(R, W, z, vcorr, cols, v_neb, sigA, half=30.0, vscan=None):
    """Simultaneous quadratic continuum + Gaussians. R, W: [N, P] residual and ivar on WAVE."""
    cent = np.array([lc for c in cols for (kind, lc), _ in c])     # every component line inside the window
    sel = (WAVE >= cent.min() - half) & (WAVE <= cent.max() + half)
    lam = WAVE[sel]
    r = R[:, sel].astype(float)
    w = W[:, sel].astype(float)
    bad = ~np.isfinite(r) | ~np.isfinite(w) | (w <= 0)
    r[bad], w[bad] = 0.0, 0.0
    N, P = r.shape
    own_lams = set(lc for c in cols for (kind, lc), _ in c)
    # masks: other nebular lines (+-3 A), sky lines (+-2.5 A) not part of this group, the object's own lines (+-8 A)
    m = np.ones(P, bool)
    for l0, hw in [(q, 3.0) for q in NEB_LINES] + [(q, 2.5) for q in SKY_LINES]:
        if any(abs(l0 - o) < 0.5 for o in own_lams):
            continue
        m &= np.abs(lam - l0) > hw
    wm = w * m[None, :]
    for gl in GAL_REST:
        wm[np.abs(lam[None, :] - gl * (1 + z)[:, None]) < 8.0] = 0.0
    x = (lam - lam.mean()) / half
    sig = np.sqrt(np.interp(lam, WAVE, sigA) ** 2 + (lam * SIG_INT / C) ** 2 + PIX ** 2 / 12)

    def design(v):
        K = len(cols)
        Gm = np.zeros((N, K, P))
        for k, c in enumerate(cols):
            for (kind, l0), a in c:
                vv = (vcorr[:, None] if kind == 'sky' else v)
                mu = l0 * (1 + np.broadcast_to(vv, (N, 1)) / C)
                s = np.interp(l0, lam, sig)
                Gm[:, k, :] += a * np.exp(-0.5 * ((lam[None, :] - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
        B = np.broadcast_to(np.stack([np.ones(P), x, x * x])[None], (N, 3, P))
        return np.concatenate([B, Gm], 1)            # [N, 3+K, P]

    def solve(v, ww):
        X = design(v)
        A = np.einsum('nkp,np,njp->nkj', X, ww, X)
        b = np.einsum('nkp,np,np->nk', X, ww, r)
        good = np.linalg.cond(A) < 1e12
        A[~good] = np.eye(A.shape[1])
        Ai = np.linalg.inv(A)
        p = np.einsum('nkj,nj->nk', Ai, b)
        model = np.einsum('nk,nkp->np', p, X)
        chi2 = np.sum(ww * (r - model) ** 2, 1)
        p[~good] = np.nan
        return p, Ai, chi2, model, good

    p, Ai, chi2, model, good = solve(v_neb, wm)
    # one 5-sigma clip away from the group's own lines
    near = np.zeros(P, bool)
    for l0 in own_lams:
        near |= np.abs(lam - l0) < 3.0
    zz = (r - model) * np.sqrt(wm)
    wm2 = wm * ((np.abs(zz) < 5) | near[None, :])
    p, Ai, chi2, model, good = solve(v_neb, wm2)
    K = len(cols)
    amp = p[:, 3:]
    err = np.sqrt(np.clip(np.einsum('nkk->nk', Ai)[:, 3:], 0, None))
    err[~good] = np.nan
    dof = (wm2 > 0).sum(1) - 3 - K
    out = dict(amp=amp, err=err, rchi2=chi2 / np.maximum(dof, 1), npix=(wm2 > 0).sum(1))
    if vscan is not None:
        chis, amps = [], []
        for v in vscan:
            pv, Aiv, c2, _, gv = solve(v, wm2)
            chis.append(c2)
            amps.append(pv[:, 3])
        chis, amps = np.array(chis), np.array(amps)
        ib = np.nanargmin(np.where(np.isfinite(chis), chis, np.inf), 0)
        n_ = np.arange(N)
        ibc = np.clip(ib, 1, len(vscan) - 2)
        ym, y0, yp = chis[ibc - 1, n_], chis[ibc, n_], chis[ibc + 1, n_]
        h = vscan[1] - vscan[0]
        curv = (yp - 2 * y0 + ym) / h ** 2
        with np.errstate(invalid='ignore', divide='ignore'):
            sh = np.where(curv > 0, -0.5 * h * (yp - ym) / (yp - 2 * y0 + ym), 0.0)
            ev = np.where(curv > 0, np.sqrt(2.0 / curv), np.nan)
        edge = (ib == 0) | (ib == len(vscan) - 1)
        out.update(v_free=np.where(edge, np.nan, vscan[ibc] + np.clip(sh, -h, h)), ev_free=np.where(edge, np.nan, ev),
                   amp_vfree=amps[ib, n_], dchi_free=chis[np.searchsorted(vscan, v_neb), n_] - chis[ib, n_],
                   chi2_scan=chis.astype(np.float32))
    return out


def main():
    v_neb = float(sys.argv[1]) if len(sys.argv) > 1 else 14.0
    ids, F, I, M, src = load_spectra()
    keep, meta = load_meta(ids)
    ids, F, I, M, src = ids[keep], F[keep], I[keep], M[keep], src[keep]
    ex = Exposures()
    co = Coadds(meta, ex)
    vcorr, vspread = co.vcorr()
    vc = np.where(np.isfinite(vcorr), vcorr, 0.0)
    sigA, nls = lsf_sigma()
    R = F - M
    del F, M
    z = meta['redshift']
    res = dict(sparcl_id=ids, src=src, v_neb=v_neb, vcorr=vcorr, vspread=vspread, lsf_sigma=sigA, n_lsf=nls)
    vscan = np.arange(-120.0, 150.1, 2.0)
    for g, cols in groups(v_neb).items():
        parts = []
        for i0 in range(0, len(ids), 1500):
            sl = slice(i0, i0 + 1500)
            parts.append(fit_group(R[sl], I[sl], z[sl], vc[sl], cols, v_neb, sigA,
                                   vscan=vscan if g == 'o3' else None))
        for k in parts[0]:
            res[f'{g}__{k}'] = np.concatenate([p_[k] for p_ in parts], axis=1 if k == 'chi2_scan' else 0)
        print(g, 'done', flush=True)
    res['vscan'] = vscan
    np.savez_compressed(D / 'verify_desi_linefluxes.npz', **res)
    print('wrote verify_desi_linefluxes.npz', len(ids))


if __name__ == '__main__':
    main()
