"""verify_desi_regspec.py -- pixel-level IV regression spectra of the [O III]-correlated component (wave-2, desi),
used to see WHERE features sit (Hg I 4358 sky line vs [O III]4363 vs [Fe II]4359; [O II]3727; Hb; [Ca II]/[O II]).

For each fibre the residual (flux - Redrock model) in each window is continuum-subtracted (quadratic, all nebular and
sky lines masked, the object's own lines masked) and multiplied by the chromatic acceptance ratio c_i(lambda).
Per pixel, S(lambda) = sum u (Zf - <Zf>)(c - <c>) / sum u (Zf - <Zf>)(O - <O>), u = ivar(lambda),
i.e. flux density per unit [O III]5007 flux of the filament component (energy, per A). A HEALPix nside=32 block
bootstrap (150 resamples) gives the band. Computed in the barycentric frame and, for the Hg test, separately for
fibres with low / high barycentric correction (a topocentric sky line moves by Delta v_corr; a nebular line does not).
Output: data/derived/verify_desi_regspec.npz, verify_desi_regspec_fits.json
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (D, C, WAVE, L, NEB_LINES, SKY_LINES, GAL_REST, load_spectra, load_meta,  # noqa
                                good_fit, lsf_sigma)

RNG = np.random.default_rng(5)
WINS = {'4363': (4301.0, 4405.0), '5007': (4930.0, 5045.0), '3727': (3690.0, 3770.0), 'Hb': (4823.0, 4902.0),
        '7300': (7282.0, 7373.0), '7150': (7100.0, 7178.0), '3869': (3830.0, 3910.0), '4686': (4647.0, 4727.0)}


def contsub(R, W, z, lo, hi):
    s = (WAVE >= lo) & (WAVE <= hi)
    lam = WAVE[s]
    r = R[:, s].astype(float)
    w = W[:, s].astype(float)
    bad = ~np.isfinite(r) | ~np.isfinite(w) | (w <= 0)
    r[bad], w[bad] = 0, 0
    m = np.ones(lam.size, bool)
    for l0 in NEB_LINES:
        m &= np.abs(lam - l0) > 3.0
    for l0 in SKY_LINES:
        m &= np.abs(lam - l0) > 2.5
    own = np.zeros(r.shape, bool)
    for gl in GAL_REST:
        own |= np.abs(lam[None, :] - gl * (1 + z)[:, None]) < 8.0
    w[own] = 0
    x = (lam - lam.mean()) / (hi - lo)
    B = np.stack([np.ones_like(x), x, x * x])
    wc = w * m[None, :]
    for it in range(2):
        A = np.einsum('kp,np,jp->nkj', B, wc, B)
        b = np.einsum('kp,np,np->nk', B, wc, r)
        ok = np.linalg.cond(A) < 1e12
        A[~ok] = np.eye(3)
        p = np.linalg.solve(A, b[..., None])[..., 0]
        cont = p @ B
        if it == 0:
            zz = (r - cont) * np.sqrt(wc)
            wc = wc * (np.abs(zz) < 4)
    w[~ok] = 0
    return lam, r - cont, w


def regspec(c, w, Zf, O, m, blocks, nb=150):
    def one(uu):
        u = w[m] * uu[:, None]
        su = u.sum(0)
        zb = (u * Zf[m, None]).sum(0) / su
        ob = (u * O[m, None]).sum(0) / su
        cb = (u * c[m]).sum(0) / su
        dz = Zf[m, None] - zb
        num = (u * dz * (c[m] - cb)).sum(0)
        den = (u * dz * (O[m, None] - ob)).sum(0)
        return num / den
    S = one(np.ones(m.sum()))
    ub, inv = np.unique(blocks[m], return_inverse=True)
    bs = np.array([one(np.bincount(RNG.integers(0, ub.size, ub.size), minlength=ub.size)[inv].astype(float))
                   for _ in range(nb)])
    return S, bs


def gauss_fit(lam, S, eS, centres, sig, free_idx=None, grid=None):
    """Linear fit of Gaussians (unit integral) + constant; optional centroid scan of one component."""
    def solve(cs):
        G = [np.exp(-0.5 * ((lam - c0) / s) ** 2) / (np.sqrt(2 * np.pi) * s) for c0, s in zip(cs, sig)]
        X = np.vstack(G + [np.ones_like(lam)])
        wv = 1 / eS ** 2
        A = (X * wv) @ X.T
        p = np.linalg.solve(A, (X * wv) @ S)
        chi2 = float(np.sum(wv * (S - p @ X) ** 2))
        return p, np.sqrt(np.diag(np.linalg.inv(A))), chi2
    if free_idx is None:
        p, e, chi2 = solve(centres)
        return dict(amp=p[:-1].tolist(), err=e[:-1].tolist(), chi2=chi2)
    chis = []
    for g in grid:
        cs = list(centres)
        cs[free_idx] = g
        chis.append(solve(cs)[2])
    chis = np.array(chis)
    k = int(np.argmin(chis))
    cs = list(centres)
    cs[free_idx] = grid[k]
    p, e, chi2 = solve(cs)
    within = grid[chis <= chis.min() + 1]
    return dict(best_centre=float(grid[k]), lo=float(within.min()), hi=float(within.max()), amp=p[:-1].tolist(),
                err=e[:-1].tolist(), chi2_min=float(chis.min()), chis=chis.tolist(), grid=grid.tolist())


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids_lf = lf['sparcl_id']
    ids, F, I, M, src = load_spectra()
    pos = {s: k for k, s in enumerate(ids)}
    k = np.array([pos[s] for s in ids_lf])
    R = (F - M)[k]
    I = I[k]
    del F, M
    keep, meta = load_meta(ids_lf)
    import healpy as hp
    blocks = hp.ang2pix(32, meta['ra'], meta['dec'], lonlat=True)
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    O = lf['o3__amp'][:, 0]
    Zf = cm['Z'] / cm['kR_5008']
    base = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & np.isfinite(Zf) & cm['ok_contrib']
    vc = cm['vcorr']
    sigA, _ = lsf_sigma()
    out, fits = {}, {}
    chrom = {'4363': '4363', '5007': '5008', '3727': '3727', 'Hb': '4861', '7300': '7325', '7150': '7136',
             '3869': '3869', '4686': '4686'}
    subsets = {'all': base, 'south': base & (meta['dec'] < 19.5),
               'vlo': base & (vc < np.nanpercentile(vc[base], 33)), 'vhi': base & (vc > np.nanpercentile(vc[base], 67))}
    for wn, (lo, hi) in WINS.items():
        lam, c, w = contsub(R, I, meta['redshift'], lo, hi)
        c = c * (cm[f'f_chrom_{chrom[wn]}'] / cm['f_chrom_5008'])[:, None]
        w = w / (cm[f'f_chrom_{chrom[wn]}'] / cm['f_chrom_5008'])[:, None] ** 2
        out[f'lam_{wn}'] = lam
        for sn, m in subsets.items():
            if wn != '4363' and sn in ('vlo', 'vhi'):
                continue
            S, bs = regspec(c, w, Zf, O, m, blocks)
            out[f'S_{wn}_{sn}'] = S
            out[f'bs_{wn}_{sn}'] = bs.astype(np.float32)
        print(wn, 'done', flush=True)
    out['vcorr_lo_median'] = float(np.nanmedian(vc[subsets['vlo']]))
    out['vcorr_hi_median'] = float(np.nanmedian(vc[subsets['vhi']]))
    # --- fits around 4363 in the regression spectra
    v = float(lf['v_neb'])
    lam = out['lam_4363']
    sel = (lam > 4346) & (lam < 4385)
    sg = lambda l0: float(np.hypot(np.interp(l0, WAVE, sigA), 0.8 / np.sqrt(12)))  # noqa: E731
    for sn in ('all', 'south', 'vlo', 'vhi'):
        S = out[f'S_4363_{sn}'][sel]
        eS = np.nanstd(out[f'bs_4363_{sn}'][:, sel], 0)
        ll = lam[sel]
        c43 = L['OIII4363'] * (1 + v / C)
        vmed = out['vcorr_lo_median'] if sn == 'vlo' else out['vcorr_hi_median'] if sn == 'vhi' else \
            float(np.nanmedian(vc[subsets[sn]]))
        chg = L['HgI4358'] * (1 + vmed / C)
        cfe = L['FeII4359'] * (1 + v / C)
        f = {}
        f['fixed_4363_Hg'] = gauss_fit(ll, S, eS, [c43, chg], [sg(c43), sg(chg)])
        f['fixed_4363_Hg_Fe'] = gauss_fit(ll, S, eS, [c43, chg, cfe], [sg(c43), sg(chg), sg(cfe)])
        f['scan_feature_near_4363'] = gauss_fit(ll, S, eS, [c43, chg], [sg(c43), sg(chg)], free_idx=0,
                                                grid=np.arange(4361.0, 4368.01, 0.05))
        f['scan_feature_near_4359'] = gauss_fit(ll, S, eS, [c43, chg], [sg(c43), sg(chg)], free_idx=1,
                                                grid=np.arange(4357.0, 4363.01, 0.05))
        f['expected'] = dict(OIII4363=c43, HgI4358_at_vcorr=chg, FeII4359=cfe, vcorr_median=vmed)
        fits[sn] = f
        print(sn, 'feature@4363 centre %.2f [%.2f,%.2f] (exp %.2f)' % (
            f['scan_feature_near_4363']['best_centre'], f['scan_feature_near_4363']['lo'],
            f['scan_feature_near_4363']['hi'], c43),
            ' feature@4359 centre %.2f [%.2f,%.2f] (Hg exp %.2f, FeII %.2f)' % (
            f['scan_feature_near_4359']['best_centre'], f['scan_feature_near_4359']['lo'],
            f['scan_feature_near_4359']['hi'], chg, cfe),
            ' amps(4363,Hg)', np.round(f['fixed_4363_Hg']['amp'], 4), np.round(f['fixed_4363_Hg']['err'], 4))
    np.savez_compressed(D / 'verify_desi_regspec.npz', **out)
    for sn in fits:
        for kk in ('scan_feature_near_4363', 'scan_feature_near_4359'):
            fits[sn][kk].pop('chis')
            fits[sn][kk].pop('grid')
    (D / 'verify_desi_regspec_fits.json').write_text(json.dumps(fits, indent=1))


if __name__ == '__main__':
    main()
