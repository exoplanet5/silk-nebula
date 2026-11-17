"""verify_desi_audit_vel.py -- second-pass (resume) audit of DESI [O III] velocities and the wavelength frame.

Own code, independent of verify_desi_velocity.py:
  * per-pixel regression of continuum-subtracted residuals around BOTH [O III]4960 and 5008 on the strand
    templates (joint Z_E, Z_W, Z_W2 with intercept; or a single strand's Z for Z>0 fibres), weights = median ivar;
  * velocity from a chi^2 grid (0.25 km/s) of a pixel-integrated Gaussian doublet model (4959 tied at 1/2.98),
    linear in the amplitudes/continua, so no non-linear optimiser can stall; HEALPix nside-32 block bootstrap;
  * frame check: Gaussian (pixel-integrated) centroid of the [O I]5577 and Hg I 4358 VARIANCE peaks (1/ivar) of
    coadds with a small v_corr spread, regressed on v_corr (barycentric spectra -> slope +1).
Output: data/derived/verify_desi_audit_vel.json
"""
import json
import sys

import numpy as np
from scipy.special import erf

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, WAVE, C, L, GAL_REST, load_spectra, lsf_sigma  # noqa: E402
from verify_desi_audit import healpix_ring  # noqa: E402

RNG = np.random.default_rng(1116)
DL = 0.8


def pix_gauss(lam, mu, sig):
    """Pixel-integrated unit-area Gaussian on 0.8 A pixels."""
    a = (lam - DL / 2 - mu) / (np.sqrt(2) * sig)
    b = (lam + DL / 2 - mu) / (np.sqrt(2) * sig)
    return 0.5 * (erf(b) - erf(a)) / DL


def prep(R, I, z, sel, lam0s, half=23.0, meta=None):
    """Continuum-subtracted residual (x f(5008)) and weights in windows around lam0s; fibres with any bad pixel
    or own-redshift line inside a window are dropped."""
    pix = np.concatenate([np.where(np.abs(WAVE - l0) < half)[0] for l0 in lam0s])
    seg = np.concatenate([np.full(np.sum(np.abs(WAVE - l0) < half), k) for k, l0 in enumerate(lam0s)])
    idx = np.where(sel)[0]
    rows, wts, keep = [], [], []
    for i in idx:
        if np.any(I[i, pix] <= 0):
            continue
        zl = GAL_REST * (1 + z[i])
        if any(np.any(np.abs(zl - l0) < half + 8) for l0 in lam0s):
            continue
        rc = np.empty(pix.size)
        for k, l0 in enumerate(lam0s):
            m = seg == k
            x = WAVE[pix][m] - l0
            cm = np.abs(x) > 3.5
            w = I[i, pix][m]
            A = np.c_[np.ones(cm.sum()), x[cm]]
            p = np.linalg.lstsq(A * np.sqrt(w[cm])[:, None], R[i, pix][m][cm] * np.sqrt(w[cm]), rcond=None)[0]
            rc[m] = R[i, pix][m] - (p[0] + p[1] * x)
        f = meta['f_chrom_5008'][i]
        rows.append(rc * f)
        wts.append(np.median(I[i, pix]) / f ** 2)
        keep.append(i)
    return np.array(keep), np.array(rows), np.array(wts), pix, seg


def regress(X, Y, w):
    XtW = X.T * w
    return np.linalg.solve(XtW @ X, XtW @ Y)


def vel_grid(lam, seg, b, sb, lam0s, ratios, sig, vgrid):
    """chi^2(v) of b = sum_k (a_k + s_k x) [segment k] + A sum_k r_k G(lam; lam0_k (1+v/c), sig_k)."""
    nseg = len(lam0s)
    chi = np.empty(vgrid.size)
    amps = np.empty(vgrid.size)
    base = []
    for k in range(nseg):
        m = (seg == k).astype(float)
        base += [m, m * (lam - lam0s[k])]
    base = np.array(base).T
    wt = 1 / sb ** 2
    for j, v in enumerate(vgrid):
        g = np.zeros(lam.size)
        for k in range(nseg):
            m = seg == k
            g[m] = ratios[k] * pix_gauss(lam[m], lam0s[k] * (1 + v / C), sig[k])
        X = np.c_[base, g]
        p = regress(X, b, wt)
        chi[j] = np.sum(wt * (b - X @ p) ** 2)
        amps[j] = p[-1]
    j = int(np.argmin(chi))
    # parabolic refinement
    if 0 < j < vgrid.size - 1:
        y0, y1, y2 = chi[j - 1:j + 2]
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) if (y0 - 2 * y1 + y2) > 0 else 0.0
        vbest = vgrid[j] + d * (vgrid[1] - vgrid[0])
    else:
        vbest = vgrid[j]
    return vbest, amps[j], chi


def main():
    meta = np.load(D / 'verify_desi_coadd_meta.npz')
    ids, F, I, M, src = load_spectra()
    pos = {s: i for i, s in enumerate(ids)}
    take = np.array([pos[s] for s in meta['sparcl_id']])
    F, I, M = F[take], I[take], M[take]
    R = F - M
    z, dec, ra = meta['redshift'], meta['dec'], meta['ra']
    ok = meta['ok_contrib'] & np.isfinite(meta['Z'])
    lsf, _ = lsf_sigma()
    lam0s = [L['OIII4959'], L['OIII5007']]
    ratios = [1 / 2.98, 1.0]
    sig = [float(np.sqrt(np.interp(l0, WAVE, lsf) ** 2 + (12.0 / C * l0) ** 2)) for l0 in lam0s]
    keep, Y, w, pix, seg = prep(R, I, z, ok, lam0s, meta=meta)
    lam = WAVE[pix]
    hp = healpix_ring(ra[keep], dec[keep], 32)
    ub, inv = np.unique(hp, return_inverse=True)
    members = [np.where(inv == b)[0] for b in range(ub.size)]
    print('fibres in velocity sample', keep.size)
    vgrid = np.arange(-40, 70.01, 0.25)
    out = dict(n=int(keep.size), lsf_sigma=sig, v_grid=[float(vgrid[0]), float(vgrid[-1])])

    def run(name, rows, instr_names, nboot=300, intercept=True):
        Xs = [meta[nm][keep][rows] for nm in instr_names]
        keepc = [j for j, x in enumerate(Xs) if np.std(x) > 1e-4]
        X = np.c_[np.ones(rows.size), *[Xs[j] for j in keepc]] if intercept else np.c_[[Xs[j] for j in keepc]].T
        B = regress(X, Y[rows], w[rows])
        sub_hp = inv[rows]
        ublk = np.unique(sub_hp)
        mem = [rows_i for rows_i in (np.where(sub_hp == b)[0] for b in ublk)]
        boots = []
        for _ in range(nboot):
            pick = RNG.integers(0, ublk.size, ublk.size)
            ii = np.concatenate([mem[b] for b in pick])
            try:
                boots.append(regress(X[ii], Y[rows][ii], w[rows][ii]))
            except np.linalg.LinAlgError:
                pass
        boots = np.array(boots)
        res = {}
        for jj, j in enumerate(keepc):
            k = jj + (1 if intercept else 0)
            b = B[k]
            q = np.percentile(boots[:, k, :], [16, 84], axis=0)
            sb = 0.5 * (q[1] - q[0])
            sb = np.maximum(sb, 0.1 * np.median(sb[sb > 0]) if np.any(sb > 0) else 1e-12)
            v, a, chi = vel_grid(lam, seg, b, sb, lam0s, ratios, sig, vgrid)
            vb = []
            for bb in boots:
                try:
                    vv, aa, _ = vel_grid(lam, seg, bb[k], sb, lam0s, ratios, sig, vgrid)
                    vb.append(vv)
                except np.linalg.LinAlgError:
                    pass
            vb = np.array(vb)
            # formal error from delta chi2 = 1
            okc = chi <= chi.min() + 1
            res[instr_names[j]] = dict(v=float(v), v_boot_p16_p84=np.percentile(vb, [16, 84]).tolist(),
                                       v_boot_sd=float(np.std(vb)), v_boot_median=float(np.median(vb)),
                                       v_formal_range=[float(vgrid[okc].min()), float(vgrid[okc].max())],
                                       amp=float(a))
        res['n_fibres'] = int(rows.size)
        res['n_blocks'] = int(ublk.size)
        out[name] = res
        print(name, rows.size, ublk.size, {k: (round(v['v'], 1), [round(x, 1) for x in v['v_boot_p16_p84']],
                                            round(v['v_formal_range'][0], 1), round(v['v_formal_range'][1], 1))
                                        for k, v in res.items() if isinstance(v, dict)})

    dk = dec[keep]
    allr = np.arange(keep.size)
    run('joint_all', allr, ['Z_E', 'Z_W', 'Z_W2'])
    run('total_all', allr, ['Z'])
    run('joint_dec<15.5', np.where(dk < 15.5)[0], ['Z_E', 'Z_W'])
    run('joint_dec15.5-17.5', np.where((dk >= 15.5) & (dk < 17.5))[0], ['Z_E', 'Z_W'])
    run('joint_dec17.5-19.5', np.where((dk >= 17.5) & (dk < 19.5))[0], ['Z_E', 'Z_W'])
    run('joint_dec>=19.5', np.where(dk >= 19.5)[0], ['Z_E', 'Z_W', 'Z_W2'])
    run('total_dec<19.5', np.where(dk < 19.5)[0], ['Z'])
    run('total_dec>=19.5', np.where(dk >= 19.5)[0], ['Z'])
    # on-ridge fibres of one strand only, regressed on that strand's template (more local)
    for s in ('E', 'W'):
        T = meta[f'T_{s}'][keep]
        other = meta['T'][keep] - T
        rows = np.where((T > 0.3) & (other < 0.1))[0]
        run(f'onridge_{s}_only', rows, [f'Z_{s}'])
        rows = np.where((T > 0.3) & (other < 0.1) & (dk >= 15.5) & (dk < 17.5))[0]
        if rows.size > 30:
            run(f'onridge_{s}_only_dec15.5-17.5', rows, [f'Z_{s}'])
    # LAMOST regions (Wu+2026 P1-P3 centres from the trace file), on-ridge fibres within 0.5 deg, total Z
    tr = json.loads((D / 'silk_trace.json').read_text())
    preg = {}
    for k in ('P1', 'P2', 'P3'):
        for container in (tr.get('paper_regions', {}), tr.get('regions', {}), tr):
            if isinstance(container, dict) and k in container and isinstance(container[k], dict):
                preg[k] = (container[k].get('ra'), container[k].get('dec'))
                break
    out['paper_regions_used'] = preg
    for k, (pra, pde) in preg.items():
        if pra is None:
            continue
        dd = np.hypot((ra[keep] - pra) * np.cos(np.radians(pde)), dec[keep] - pde)
        rows = np.where((dd < 0.5) & (meta['T'][keep] > 0.3))[0]
        if rows.size > 20:
            run(f'onridge_near_{k}', rows, ['Z'])

    # -------------------------------------------------------------- frame check from variance peaks
    frame = {}
    small = ok & (meta['vspread'] < 3.0) & np.isfinite(meta['vcorr'])
    for nm, l0 in (('sky5577', L['sky5577']), ('HgI4358', L['HgI4358'])):
        pw = np.where(np.abs(WAVE - l0) < 9)[0]
        lamw = WAVE[pw]
        sg = float(np.interp(l0, WAVE, lsf))
        mus = np.arange(l0 - 1.5, l0 + 1.5, 0.005)
        G = np.array([pix_gauss(lamw, mu, sg) for mu in mus])
        vs, vc = [], []
        for i in np.where(small)[0]:
            if np.any(I[i, pw] <= 0):
                continue
            V = 1 / I[i, pw]
            best, bchi = None, np.inf
            for jm in range(mus.size):
                X = np.c_[np.ones(lamw.size), lamw - l0, G[jm]]
                p = np.linalg.lstsq(X, V, rcond=None)[0]
                r = np.sum((V - X @ p) ** 2)
                if r < bchi and p[2] > 0:
                    bchi, best = r, mus[jm]
            if best is None:
                continue
            vs.append((best / l0 - 1) * C)
            vc.append(meta['vcorr'][i])
        vs, vc = np.array(vs), np.array(vc)
        A = np.c_[np.ones(vc.size), vc]
        keepm = np.ones(vc.size, bool)
        for _ in range(5):
            p = np.linalg.lstsq(A[keepm], vs[keepm], rcond=None)[0]
            r = vs - A @ p
            s = 1.4826 * np.median(np.abs(r[keepm]))
            keepm = np.abs(r) < 4 * s
        cov = np.linalg.inv(A[keepm].T @ A[keepm]) * np.var(r[keepm])
        frame[nm] = dict(n=int(keepm.sum()), slope=float(p[1]), slope_err=float(np.sqrt(cov[1, 1])),
                         intercept_kms=float(p[0]), rms=float(np.std(r[keepm])),
                         vcorr_range=[float(vc.min()), float(vc.max())])
        print('frame', nm, frame[nm])
    out['frame_variance_peaks'] = frame
    json.dump(out, open(D / 'verify_desi_audit_vel.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
