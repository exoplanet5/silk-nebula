"""w3_northsouth_desi.py -- wave-3 (tag northsouth): DESI DR1 [O III] brightness, line ratios and velocity in
Dec bins along the Silk spine, with the wave-2-corrected methods.

Method (own code; re-uses only the wave-2 verify_desi data products and helper functions, read-only):
  * Box-integral line fluxes, GFA fibre-acceptance calibration and chromatic correction from
    data/derived/verify_desi_audit_boxfluxes.npz (verify_desi revision 2).
  * NEW: Dec-binned geometry templates.  Every strand sample (E, W, W2 polylines, Gaussian cross-section at the
    trace HWHM) is assigned to the Dec bin of that trace sample; Z_k(fibre) = T_k(fibre) - sum_t w_t <T_k>_(tile,petal)
    (DESI's per-petal sky subtraction, mean over all fiberassign positions, as in verify_desi_common).
    All bins are fitted JOINTLY to all fibres, so the petal-mean contamination of one bin by another is modelled
    (the wave-2 caveat "regression centroids of local subsets are not local" does not apply).
  * Brightness: WLS of O[R] on [1, Z_1..Z_K]; a_k = ridge peak SB (R) of bin k.
  * Ratios: IV-type slope ratio rho_k = b_k(L) / b_k(O), both regressions with the line's inverse variance.
  * Velocity: per-pixel WLS regression of the continuum-subtracted residual around [O III]4960+5008 (prep() of
    verify_desi_audit_vel) on [1, Z_1..Z_K]; each bin's coefficient spectrum is fitted with a pixel-integrated
    Gaussian doublet (4959 tied 1/2.98) on a 0.25 km/s grid.
  * Errors: (i) formal WLS errors (scaled by the reduced chi^2) x an empirical inflation factor from 60 rigid
    placebo placements of the whole trace (non-overlapping, |r| < 0.05), fitted jointly with the true bins;
    (ii) HEALPix nside-32 and nside-64 block bootstraps.  Fine (1.5 deg) bins span only 1-3 nside-32 blocks, so
    their bootstrap is unstable; the placebo-inflated formal error is the primary error for bins.
Variants: (A) all strands per Dec bin [baseline]; (B) E-strand bins with W and W2 as separate nuisance terms;
(C) baseline with Halpha-patch (30 % contour) fibres removed; (D) two halves (trace Dec < / >= 19.5).
Outputs: data/derived/w3_northsouth_desi.json, data/derived/w3_northsouth_desi_spec.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w3_northsouth_desi.py
"""
import json
import sys
import time

import numpy as np
from matplotlib.path import Path as MPath

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (D, WAVE, C, L, load_spectra, lsf_sigma, Exposures, Coadds, Template,  # noqa
                                unitvec, rotmat)
from verify_desi_audit import healpix_ring  # noqa: E402
from verify_desi_audit_vel import prep, pix_gauss  # noqa: E402

RNG = np.random.default_rng(31116)
NBOOT = 400
EDGES = np.array([12.0, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0, 25.5, 27.0, 28.6])
OMEGA = np.pi * 0.76 ** 2
V_EDGES = np.array([12.0, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 25.5, 28.6])
NPLACE = 60


def bin_values(T, ra, dec, edges, strands=None, R=None):
    """(nbins, n) template values; each strand sample contributes to the bin of its own (trace) Dec.
    If the template was built rotated (Template(R=R)), bins use the unrotated Dec of each sample (P @ R)."""
    X = unitvec(ra, dec)
    nb = len(edges) - 1
    out = np.zeros((nb, len(ra)))
    for s, tree, sig, a, P in T.parts:
        if strands is not None and s not in strands:
            continue
        dist, idx = tree.query(X, distance_upper_bound=np.radians(2.5))
        ok = np.isfinite(dist)
        ii = np.where(ok)[0]
        d = np.degrees(dist[ok])
        val = a * np.exp(-0.5 * (d / sig[idx[ok]]) ** 2)
        P0 = P[idx[ok]] if R is None else P[idx[ok]] @ R
        sdec = np.degrees(np.arcsin(np.clip(P0[:, 2], -1, 1)))
        k = np.digitize(sdec, edges) - 1
        g = (k >= 0) & (k < nb)
        np.add.at(out, (k[g], ii[g]), val[g])
    return out


def petal_sub(T, co, ra, dec, edges, strands=None, R=None):
    keys = co.pet_keys
    allra = np.concatenate([co.pet_ra[k] for k in keys])
    alldec = np.concatenate([co.pet_dec[k] for k in keys])
    lens = np.array([co.pet_ra[k].size for k in keys])
    ends = np.cumsum(lens)
    starts = ends - lens
    va = bin_values(T, allra, alldec, edges, strands, R)
    vf = bin_values(T, ra, dec, edges, strands, R)
    Z = np.zeros_like(vf)
    for kb in range(va.shape[0]):
        cs = np.concatenate([[0.0], np.cumsum(va[kb])])
        pm = {k: (cs[e] - cs[s]) / (e - s) for k, s, e in zip(keys, starts, ends)}
        sub = np.array([sum(w * pm[(t, p)] for t, p, w in c) if c else np.nan for c in co.contrib])
        Z[kb] = vf[kb] - sub
    return Z, vf


def wls(X, y, w):
    XtW = X.T * w
    return np.linalg.solve(XtW @ X, XtW @ y)


class VelFitter:
    """chi^2(v) grid fit of a pixel-integrated Gaussian doublet + per-segment linear continuum (batched)."""

    def __init__(self, lam, seg, lam0s, ratios, sig, vgrid):
        self.vgrid = vgrid
        nseg = len(lam0s)
        base = []
        for k in range(nseg):
            m = (seg == k).astype(float)
            base += [m, m * (lam - lam0s[k])]
        self.base = np.array(base).T                       # (npix, 2 nseg)
        G = np.zeros((vgrid.size, lam.size))
        for j, v in enumerate(vgrid):
            for k in range(nseg):
                m = seg == k
                G[j, m] = ratios[k] * pix_gauss(lam[m], lam0s[k] * (1 + v / C), sig[k])
        self.G = G
        nv, npix = G.shape
        self.A = np.concatenate([np.broadcast_to(self.base, (nv, npix, self.base.shape[1])), G[:, :, None]], 2)

    def fit(self, b, sb):
        wt = 1 / sb ** 2
        AtW = np.transpose(self.A, (0, 2, 1)) * wt[None, None, :]
        M = AtW @ self.A
        r = AtW @ b
        p = np.linalg.solve(M, r[..., None])[..., 0]
        res = b[None, :] - np.einsum('vpk,vk->vp', self.A, p)
        chi = np.sum(wt * res ** 2, 1)
        j = int(np.argmin(chi))
        vb = self.vgrid[j]
        if 0 < j < self.vgrid.size - 1:
            y0, y1, y2 = chi[j - 1:j + 2]
            den = y0 - 2 * y1 + y2
            if den > 0:
                vb = self.vgrid[j] + 0.5 * (y0 - y2) / den * (self.vgrid[1] - self.vgrid[0])
        okc = chi <= chi.min() + 1
        return float(vb), float(p[j, -1]), (float(self.vgrid[okc].min()), float(self.vgrid[okc].max())), chi


def summ(x, est=None):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return None
    p = np.percentile(x, [2.5, 16, 50, 84, 97.5])
    d = dict(p2p5=float(p[0]), p16=float(p[1]), p50=float(p[2]), p84=float(p[3]), p97p5=float(p[4]),
             sig=float(0.5 * (p[3] - p[1])), sd=float(np.std(x)))
    if est is not None:
        d['est'] = float(est)
    return d


def robust_sd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else np.nan


def placement_rotations(T, ra, dec, Ttrue, nwant, rng):
    """Rigid placebo placements of the whole trace: spin about the trace centre (+-20 deg) then a 1-4 deg
    translation in a random direction; kept if the template correlation with the true one is |r| < 0.05."""
    allpts = np.concatenate([p[4] for p in T.parts])
    c = allpts.mean(0)
    c /= np.linalg.norm(c)
    out = []
    tries = 0
    while len(out) < nwant and tries < 20 * nwant:
        tries += 1
        Rs = rotmat(c, rng.uniform(-20, 20))
        a = np.cross(c, rng.normal(size=3))
        a /= np.linalg.norm(a)
        Rt = rotmat(a, rng.uniform(1.0, 4.0))
        R = Rt @ Rs
        Tn = Template(R=R)
        tv = Tn(ra, dec)
        if np.std(tv) == 0:
            continue
        r = np.corrcoef(tv, Ttrue)[0, 1]
        if abs(r) < 0.05:
            out.append((R, Tn, float(r)))
    return out


def main():
    t0 = time.time()
    meta = np.load(D / 'verify_desi_coadd_meta.npz')
    box = np.load(D / 'verify_desi_audit_boxfluxes.npz')
    assert (meta['sparcl_id'] == box['sparcl_id']).all()
    ra, dec, z = meta['ra'], meta['dec'], meta['redshift']
    n = ra.size
    ok = meta['ok_contrib'] & np.isfinite(meta['Z']) & np.isfinite(meta['f_chrom_5008'])
    K = len(EDGES) - 1
    # ---------------------------------------------------------------- Dec-bin templates
    mdict = {k: meta[k] for k in ('targetid', 'survey', 'program', 'ra', 'dec')}
    ex = Exposures()
    co = Coadds(mdict, ex)
    T = Template()
    Zb, Tb = petal_sub(T, co, ra, dec, EDGES)                       # variant A: all strands
    chk = np.nanmax(np.abs(np.nansum(Zb, 0)[ok] - meta['Z'][ok]))
    print('bins', K, 'check sum_k Z_k vs meta Z: max |diff| =', chk, f'({time.time()-t0:.0f}s)', flush=True)
    ZbE, _ = petal_sub(T, co, ra, dec, EDGES, strands=('E',))      # variant B: E-strand bins
    ZW, _ = petal_sub(T, co, ra, dec, np.array([12.0, 17.5, 20.5]), strands=('W',))
    ZW2, _ = petal_sub(T, co, ra, dec, np.array([20.0, 25.0]), strands=('W2',))
    Zh, _ = petal_sub(T, co, ra, dec, np.array([0.0, 19.5, 30.0]))  # variant D: halves
    Z3, _ = petal_sub(T, co, ra, dec, np.array([12.0, 17.6, 21.0, 28.6]))
    Zv, Tv = petal_sub(T, co, ra, dec, V_EDGES)
    ZvE, _ = petal_sub(T, co, ra, dec, V_EDGES, strands=('E',))
    tr = json.loads((D / 'silk_trace.json').read_text())
    patch30 = np.array(tr['halpha_patch_info']['outer_polygon_30pct'], float)
    inpatch = MPath(patch30).contains_points(np.c_[ra, dec])
    print('fibres in Halpha-patch 30% polygon:', int((inpatch & ok).sum()), flush=True)

    # ---------------------------------------------------------------- box fluxes
    kR = 1e-17 / OMEGA / float(meta['one_R_5008'])
    E = {k[2:]: box[k] for k in box.files if k.startswith('E_')}
    U = {k[2:]: box[k] for k in box.files if k.startswith('U_')}
    E['SII'] = E['SII6716'] + E['SII6731']
    U['SII'] = 1 / (1 / U['SII6716'] + 1 / U['SII6731'])
    E['O3R'] = E['OIII5007'] * kR
    U['O3R'] = U['OIII5007'] / kR ** 2
    good = {}
    for nm in E:
        e = 1 / np.sqrt(U[nm])
        m = ok & np.isfinite(E[nm]) & np.isfinite(e) & (e > 0)
        m &= e < 4 * np.nanmedian(e[m])
        m &= np.abs(E[nm] / e) < 30
        good[nm] = m
    good['O3R'] = good['OIII5007']
    lines = ['OII3727sum', 'Hb', 'Ha', 'NII6583', 'SII', 'OIII4959', 'NeIII3869', 'OIII4363', 'HeII4686',
             'SIII9069', 'null_OII3726_-12', 'null_Hb_+11', 'null_OIII5007_+14', 'null_OIII5007_-14', 'sky5577']

    def boot_sets(nside):
        hp = healpix_ring(ra, dec, nside)
        ub, inv = np.unique(hp, return_inverse=True)
        mem = [np.where(inv == b)[0] for b in range(ub.size)]
        return [np.concatenate([mem[b] for b in RNG.integers(0, ub.size, ub.size)]) for _ in range(NBOOT)]
    BS = {'nside32': boot_sets(32), 'nside64': boot_sets(64)}

    def wls_cov(X, y, w):
        XtW = X.T * w
        Minv = np.linalg.inv(XtW @ X)
        b = Minv @ (XtW @ y)
        r = y - X @ b
        s2 = float(np.sum(w * r ** 2) / max(len(y) - X.shape[1], 1))
        return b, Minv, s2

    # ---------------------------------------------------------------- placebo inflation of formal errors
    Ttrue = meta['T']
    plc = placement_rotations(T, ra[ok], dec[ok], Ttrue[ok], NPLACE, RNG)
    print('placebo placements', len(plc), f'({time.time()-t0:.0f}s)', flush=True)
    chis = {nm: [] for nm in ['O3R'] + lines}
    for R, Tn, r in plc:
        Zn, _ = petal_sub(Tn, co, ra, dec, EDGES, R=R)
        for nm in chis:
            sel = good[nm] & good['OIII5007']
            idx = np.where(sel)[0]
            wz = np.array([np.sum(U[nm][idx] * Zn[k][idx] ** 2) for k in range(K)])
            wt = np.array([np.sum(U[nm][idx] * Zb[k][idx] ** 2) for k in range(K)])
            use = [k for k in range(K) if wz[k] > 0.2 * wt[k]]   # placed bin must sample comparable power
            if not use:
                continue
            X = np.c_[np.ones(n), *Zb, *[Zn[k] for k in use]]
            try:
                b, Minv, s2 = wls_cov(X[idx], E[nm][idx], U[nm][idx])
            except np.linalg.LinAlgError:
                continue
            for jj, k in enumerate(use):
                chis[nm].append(b[1 + K + jj] / np.sqrt(Minv[1 + K + jj, 1 + K + jj] * s2))
    infl = {nm: dict(robust=robust_sd(v), sd=float(np.std(v)), n=len(v), mean=float(np.mean(v)) if v else None)
            for nm, v in chis.items()}
    print('inflation (robust)', {k: round(v['robust'], 2) for k, v in infl.items()}, f'({time.time()-t0:.0f}s)',
          flush=True)

    def fit_block(Zs, base_sel, tag):
        Kz = len(Zs)
        X = np.c_[np.ones(n), *Zs]
        res = {'n_templates': Kz, 'quantities': {}}
        for nm in ['O3R'] + lines:
            sel = base_sel & good[nm] & good['OIII5007']
            idx = np.where(sel)[0]
            msk = np.zeros(n, bool)
            msk[idx] = True
            w = U[nm]
            bl, Ml, s2l = wls_cov(X[idx], E[nm][idx], w[idx])
            bo, Mo, s2o = wls_cov(X[idx], E['OIII5007'][idx], w[idx])
            F_l = max(infl[nm]['robust'], 1.0)
            F_o = max(infl['O3R']['robust'], 1.0)
            sl = np.sqrt(np.diag(Ml) * s2l)[1:] * F_l
            so = np.sqrt(np.diag(Mo) * s2o)[1:] * F_o
            q = {'slope': bl[1:].tolist(), 'slope_err_infl': sl.tolist(), 'n': int(idx.size),
                 'sum_wZ2': [float(np.sum(w[idx] * Zs[k][idx] ** 2)) for k in range(Kz)]}
            if nm != 'O3R':
                rho = bl[1:] / bo[1:]
                q['ratio'] = rho.tolist()
                q['ratio_err_infl'] = (np.sqrt(sl ** 2 + (rho * so) ** 2) / np.abs(bo[1:])).tolist()
            for bsn, sets in BS.items():
                bt = []
                for bi in sets:
                    ii = bi[msk[bi]]
                    try:
                        a1 = wls(X[ii], E[nm][ii], w[ii])
                        if nm == 'O3R':
                            bt.append(a1[1:])
                        else:
                            bt.append(a1[1:] / wls(X[ii], E['OIII5007'][ii], w[ii])[1:])
                    except np.linalg.LinAlgError:
                        pass
                bt = np.array(bt)
                q[f'boot_{bsn}'] = [summ(bt[:, k]) for k in range(Kz)]
            res['quantities'][nm] = q
        b = res['quantities']['O3R']
        print(tag, 'SB(R)', [f"{v:.2f}+-{e:.2f}" for v, e in zip(b['slope'], b['slope_err_infl'])], flush=True)
        o2 = res['quantities']['OII3727sum']
        print(tag, 'OII/OIII', [f"{v:.2f}+-{e:.2f}" for v, e in zip(o2['ratio'], o2['ratio_err_infl'])], flush=True)
        hb = res['quantities']['Hb']
        print(tag, 'Hb/OIII', [f"{v:.3f}+-{e:.3f}" for v, e in zip(hb['ratio'], hb['ratio_err_infl'])], flush=True)
        return res

    out = dict(edges=EDGES.tolist(), v_edges=V_EDGES.tolist(), nboot=NBOOT, n_usable=int(ok.sum()),
               check_sumZ_vs_metaZ=float(chk), placebo=dict(n=len(plc), r=[p[2] for p in plc], inflation=infl),
               fibres_per_bin_T_gt_0p3=[int((ok & (Tb[k] > 0.3)).sum()) for k in range(K)],
               fibres_per_vbin_T_gt_0p3=[int((ok & (Tv[k] > 0.3)).sum()) for k in range(len(V_EDGES) - 1)])
    out['A_all_strands'] = fit_block(list(Zb), ok, 'A')
    out['B_Estrand_bins'] = fit_block(list(ZbE) + list(ZW) + list(ZW2), ok, 'B')
    out['B_labels'] = ([f'E{lo:g}-{hi:g}' for lo, hi in zip(EDGES[:-1], EDGES[1:])] +
                       ['W12-17.5', 'W17.5-20.5', 'W2'])
    out['C_no_patch'] = fit_block(list(Zb), ok & ~inpatch, 'C')
    out['D_halves'] = fit_block(list(Zh), ok, 'D')
    out['E_three_zones'] = fit_block(list(Z3), ok, 'E3')
    print(f'ratios done ({time.time()-t0:.0f}s)', flush=True)

    # ---------------------------------------------------------------- velocity
    ids, F, I, M, src = load_spectra()
    pos = {s: i for i, s in enumerate(ids)}
    take = np.array([pos[s] for s in meta['sparcl_id']])
    R = F[take] - M[take]
    I = I[take]
    del F, M
    lsf, _ = lsf_sigma()
    lam0s = [L['OIII4959'], L['OIII5007']]
    ratios_d = [1 / 2.98, 1.0]
    sig = [float(np.sqrt(np.interp(l0, WAVE, lsf) ** 2 + (12.0 / C * l0) ** 2)) for l0 in lam0s]
    keep, Y, w, pix, seg = prep(R, I, z, ok, lam0s, meta=meta)
    lam = WAVE[pix]
    print('velocity sample', keep.size, f'({time.time()-t0:.0f}s)', flush=True)
    vf = VelFitter(lam, seg, lam0s, ratios_d, sig, np.arange(-40, 70.01, 0.25))
    hpk = healpix_ring(ra[keep], dec[keep], 32)
    ub2, inv2 = np.unique(hpk, return_inverse=True)
    mem2 = [np.where(inv2 == b)[0] for b in range(ub2.size)]
    bidx2 = [np.concatenate([mem2[b] for b in RNG.integers(0, ub2.size, ub2.size)]) for _ in range(NBOOT)]
    F_v = max(infl['O3R']['robust'], 1.0)
    specs = {}

    def vel_block(Zs, rows, tag, labels):
        X = np.c_[np.ones(rows.size), *[Zz[keep][rows] for Zz in Zs]]
        Yr, wr = Y[rows], w[rows]
        XtW = X.T * wr
        Minv = np.linalg.inv(XtW @ X)
        B = Minv @ (XtW @ Yr)
        resid = Yr - X @ B
        s2 = np.sum(wr[:, None] * resid ** 2, 0) / max(rows.size - X.shape[1], 1)
        inrow = np.zeros(keep.size, bool)
        inrow[rows] = True
        pos_r = -np.ones(keep.size, int)
        pos_r[rows] = np.arange(rows.size)
        boots = []
        for bi in bidx2:
            ii = pos_r[bi[inrow[bi]]]
            try:
                boots.append(wls(X[ii], Yr[ii], wr[ii]))
            except np.linalg.LinAlgError:
                pass
        boots = np.array(boots)
        res = []
        for k in range(len(Zs)):
            b = B[k + 1]
            sbf = np.sqrt(Minv[k + 1, k + 1] * s2) * F_v                 # formal x placebo inflation
            v, amp, frange, chi = vf.fit(b, sbf)
            vb, ab = [], []
            for bb in boots:
                vv, aa, _, _ = vf.fit(bb[k + 1], sbf)
                vb.append(vv)
                ab.append(aa)
            vb, ab = np.array(vb), np.array(ab)
            # formal amplitude error at the best velocity
            j = int(np.argmin(np.abs(vf.vgrid - v)))
            A = vf.A[j]
            Wt = 1 / sbf ** 2
            cov = np.linalg.inv((A.T * Wt) @ A)
            amp_err = float(np.sqrt(cov[-1, -1]))
            edge = (frange[0] <= vf.vgrid[0] + 0.1) or (frange[1] >= vf.vgrid[-1] - 0.1)
            snr = amp / amp_err
            res.append(dict(label=labels[k], v_helio=v, v_err_formal_infl=0.5 * (frange[1] - frange[0]),
                            v_formal_range=frange, boot=summ(vb), amp_R=float(amp * kR), amp_err_R=float(amp_err * kR),
                            amp_snr_formal_infl=float(snr), amp_boot_robust_sd_R=robust_sd(ab) * kR,
                            reliable=bool(snr >= 3 and not edge), sum_wZ2=float(np.sum(wr * X[:, k + 1] ** 2))))
            specs[f'{tag}_{k}'] = np.vstack([lam, b, sbf])
        print(tag, [(r['label'], round(r['v_helio'], 1), round(r['v_err_formal_infl'], 1),
                     round(r['boot']['sig'], 1) if r['boot'] else None, round(r['amp_snr_formal_infl'], 1))
                    for r in res], flush=True)
        return res

    allr = np.arange(keep.size)
    labV = [f'{lo:g}-{hi:g}' for lo, hi in zip(V_EDGES[:-1], V_EDGES[1:])]
    out['vel_A_bins'] = vel_block(list(Zv), allr, 'velA', labV)
    out['vel_B_Ebins'] = vel_block(list(ZvE) + list(ZW) + list(ZW2), allr, 'velB',
                                   ['E' + x for x in labV] + ['W12-17.5', 'W17.5-20.5', 'W2'])
    out['vel_D_halves'] = vel_block(list(Zh), allr, 'velD', ['dec<19.5', 'dec>=19.5'])
    out['vel_D_halves_no_patch'] = vel_block(list(Zh), np.where(~inpatch[keep])[0], 'velDnp', ['dec<19.5', 'dec>=19.5'])
    out['vel_E_three_zones'] = vel_block(list(Z3), allr, 'velE', ['<17.6', '17.6-21', '>=21'])
    # split scan: two-zone model with a moving split
    scan = []
    for sdec in (16.5, 17.6, 18.5, 19.5, 20.5, 21.5, 23.0):
        Zs2, _ = petal_sub(T, co, ra, dec, np.array([0.0, sdec, 30.0]))
        r2 = vel_block(list(Zs2), allr, f'split{sdec:g}', [f'<{sdec:g}', f'>={sdec:g}'])
        scan.append(dict(split=sdec, south=r2[0], north=r2[1]))
    out['vel_split_scan'] = scan
    # audit-style reproduction (fibres selected by fibre Dec, global template)
    dk = dec[keep]
    for nm, rows in (('audit_style_dec<19.5', np.where(dk < 19.5)[0]),
                     ('audit_style_dec>=19.5', np.where(dk >= 19.5)[0])):
        out[f'vel_{nm}'] = vel_block([meta['Z']], rows, nm, [nm])
    out['lsf_sigma_A'] = sig
    out['velocity_zero_point_note'] = ('DESI barycentric frame confirmed (verify_desi); absolute zero point +-1.5 km/s '
                                       '(stellar RVs vs Gaia +1..+3). V_LSR = V_helio - 12.1 km/s.')
    json.dump(out, open(D / 'w3_northsouth_desi.json', 'w'), indent=1)
    np.savez_compressed(D / 'w3_northsouth_desi_spec.npz', **specs)
    print(f'done ({time.time()-t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
