"""verify_desi_velocity.py -- (d) [O III] velocity of the Silk filament in DESI DR1 and the wavelength frame.

1. Centroid of the [O III]4959+5007 doublet in the pixel-level IV regression spectrum (own code, as in
   verify_desi_regspec.py) for all fibres, E-strand fibres, W-strand fibres and Dec bands; velocity by chi^2 scan of a
   tied Gaussian doublet (LSF from SPARCL wave_sigma + pixel + 12 km/s); errors from a HEALPix nside=32 block
   bootstrap (each bootstrap spectrum re-fitted).
   E / W membership: unsubtracted template T_E > 0.3 and T_W < 0.1 (E), and vice versa (W).
2. Per-fibre free-velocity fits (my o3 scan) for |S/N| > 5, inverse-variance means, E/W, and slope vs v_corr.
3. Frame test with the night sky: the [O I]5577 sky line raises the pixel variance; the centroid of 1/ivar excess in
   every coadd (single night, v_corr spread < 1 km/s) is regressed on the coadd's barycentric correction. Barycentric
   wavelengths => slope +1 (sky line appears at +v_corr) and intercept 0 (if the vacuum 5578.888 A is right).
   The same for Hg I 4358 (weaker). Tiles 41485/41487 (wave-1 -19 km/s stellar offset) are checked separately.
4. List of [O III]-bearing fibres (|S/N| > 3 or within 0.3 deg of a strand) that include tiles 41485/41487.
Output: data/derived/verify_desi_velocity.json
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, C, WAVE, L, load_spectra, load_meta, good_fit, lsf_sigma  # noqa: E402
from verify_desi_regspec import contsub, regspec  # noqa: E402

RNG = np.random.default_rng(9)


def doublet_v(lam, S, eS, sigA, vgrid=np.arange(-40.0, 70.01, 0.25)):
    l1, l2 = L['OIII5007'], L['OIII4959']
    s1 = np.sqrt(np.interp(l1, WAVE, sigA) ** 2 + (l1 * 12 / C) ** 2 + 0.8 ** 2 / 12)
    s2 = np.sqrt(np.interp(l2, WAVE, sigA) ** 2 + (l2 * 12 / C) ** 2 + 0.8 ** 2 / 12)
    w = 1 / eS ** 2
    chis = []
    for v in vgrid:
        g = (np.exp(-0.5 * ((lam - l1 * (1 + v / C)) / s1) ** 2) / s1 +
             np.exp(-0.5 * ((lam - l2 * (1 + v / C)) / s2) ** 2) / s2 / 2.98) / np.sqrt(2 * np.pi)
        X = np.vstack([g, np.ones_like(lam)])
        A = (X * w) @ X.T
        p = np.linalg.solve(A, (X * w) @ S)
        chis.append(np.sum(w * (S - p @ X) ** 2))
    chis = np.array(chis)
    k = int(np.argmin(chis))
    return float(vgrid[k]), chis


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids_lf = lf['sparcl_id']
    keep, meta = load_meta(ids_lf)
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    sc = det['err_scale_from_nulls']
    O, EO = lf['o3__amp'][:, 0], lf['o3__err'][:, 0] * sc
    Zf = cm['Z'] / cm['kR_5008']
    base = good_fit(lf, 'o3', scale=sc) & np.isfinite(Zf) & cm['ok_contrib']
    import healpy as hp
    blocks = hp.ang2pix(32, meta['ra'], meta['dec'], lonlat=True)
    sigA, _ = lsf_sigma()
    ids, F, I, M, _ = load_spectra()
    pos = {s: k for k, s in enumerate(ids)}
    kk = np.array([pos[s] for s in ids_lf])
    R, I = (F - M)[kk], I[kk]
    Fsk = F[kk]
    del F, M
    out = {}
    # ---------------- 1. regression-spectrum centroids
    lam, c, w = contsub(R, I, meta['redshift'], 4930.0, 5045.0)
    sel = (lam > 4945) & (lam < 5025)
    E_ = (cm['T_E'] > 0.3) & (cm['T_W'] < 0.1)
    W_ = (cm['T_W'] > 0.3) & (cm['T_E'] < 0.1)
    dec = meta['dec']
    subsets = {'all': base, 'E_strand': base & E_, 'W_strand': base & W_,
               'E_strand_dec15.5-17.5': base & E_ & (dec >= 15.5) & (dec < 17.5),
               'W_strand_dec15.5-17.5': base & W_ & (dec >= 15.5) & (dec < 17.5),
               'dec12-15.5': base & (dec < 15.5), 'dec15.5-17.5': base & (dec >= 15.5) & (dec < 17.5),
               'dec17.5-19.5': base & (dec >= 17.5) & (dec < 19.5), 'dec19.5-22': base & (dec >= 19.5) & (dec < 22),
               'dec22-25': base & (dec >= 22) & (dec < 25), 'dec25-28.6': base & (dec >= 25),
               'dark': base & (meta['program'] == 'dark'), 'bright': base & (meta['program'] == 'bright'),
               'main': base & (meta['survey'] == 'main'), 'sv': base & (meta['survey'] != 'main')}
    from silk_geometry import REGIONS
    for nm in ('P1', 'P2', 'P3'):
        x0, y0 = REGIONS[nm]
        dd = np.hypot((meta['ra'] - x0) * np.cos(np.radians(y0)), dec - y0)
        subsets[f'within0.5deg_{nm}'] = base & (dd < 0.5)
    rs = {}
    vf, evf = lf['o3__v_free'], lf['o3__ev_free']
    snr_ = O / EO
    for sn, m in subsets.items():
        if m.sum() < 40:
            continue
        # per-fibre companion estimate: positive-residual fibres with S/N > 4, error from the scatter
        mp = m & (snr_ > 4) & np.isfinite(vf) & np.isfinite(evf) & (evf < 30)
        pf_ = None
        if mp.sum() >= 3:
            ww_ = 1 / (evf[mp] ** 2 + 5.0 ** 2)
            pf_ = dict(n=int(mp.sum()), wmean=float(np.sum(ww_ * vf[mp]) / ww_.sum()),
                       err=float(max(ww_.sum() ** -0.5, np.std(vf[mp]) / np.sqrt(mp.sum()))),
                       n_blocks=int(np.unique(blocks[m]).size))
        S, bs = regspec(c, w, Zf, O, m, blocks, nb=150)
        eS = np.nanstd(bs[:, sel], 0)
        okp = np.isfinite(S[sel]) & (eS > 0)
        v0, chis = doublet_v(lam[sel][okp], S[sel][okp], eS[okp], sigA)
        vb = []
        for b in range(bs.shape[0]):
            if np.all(np.isfinite(bs[b, sel][okp])):
                vb.append(doublet_v(lam[sel][okp], bs[b, sel][okp], eS[okp], sigA)[0])
        vb = np.array(vb)
        rs[sn] = dict(n=int(m.sum()), v=v0, p16=float(np.percentile(vb, 16)), p84=float(np.percentile(vb, 84)),
                      sig=float(0.5 * (np.percentile(vb, 84) - np.percentile(vb, 16))),
                      vcorr_median=float(np.nanmedian(cm['vcorr'][m])), n_blocks=int(np.unique(blocks[m]).size),
                      per_fibre_pos_SN4=pf_)
        print(sn, rs[sn], flush=True)
    out['regspec_centroids'] = rs
    # ---------------- 2. per-fibre velocities
    v, ev = lf['o3__v_free'], lf['o3__ev_free']
    snr = O / EO
    pf = {}
    for lab, m in (('pos_SN>5', base & (snr > 5)), ('neg_SN>5', base & (snr < -5)),
                   ('pos_SN>5_E', base & (snr > 5) & E_), ('pos_SN>5_W', base & (snr > 5) & W_),
                   ('neg_SN>5_E', base & (snr < -5) & E_), ('neg_SN>5_W', base & (snr < -5) & W_)):
        mm = m & np.isfinite(v) & np.isfinite(ev) & (ev < 30)
        if mm.sum() < 3:
            continue
        ww = 1 / (ev[mm] ** 2 + 5.0 ** 2)
        vm = float(np.sum(ww * v[mm]) / ww.sum())
        pf[lab] = dict(n=int(mm.sum()), wmean=vm, err=float(max(ww.sum() ** -0.5, np.std(v[mm]) / np.sqrt(mm.sum()))),
                       median=float(np.median(v[mm])), std=float(np.std(v[mm])))
    mm = base & (np.abs(snr) > 5) & np.isfinite(v) & (ev < 30) & np.isfinite(cm['vcorr'])
    X = np.c_[np.ones(mm.sum()), cm['vcorr'][mm]]
    ww = 1 / (ev[mm] ** 2 + 5.0 ** 2)
    A = X.T @ (X * ww[:, None])
    b = np.linalg.solve(A, X.T @ (ww * v[mm]))
    res = v[mm] - X @ b
    s2 = np.sum(ww * res ** 2) / (mm.sum() - 2)
    cov = np.linalg.inv(A) * max(s2, 1)
    pf['v_vs_vcorr'] = dict(n=int(mm.sum()), intercept=float(b[0]), slope=float(b[1]),
                            err_slope=float(np.sqrt(cov[1, 1])), vcorr_range=[float(cm['vcorr'][mm].min()),
                                                                             float(cm['vcorr'][mm].max())])
    # same with Dec-band fixed effects (v_corr is confounded with position through the observing dates)
    band = np.digitize(dec[mm], [15.5, 17.5, 19.5, 22, 25])
    Xf = np.c_[np.eye(6)[band], cm['vcorr'][mm]]
    keepc = Xf[:, :6].sum(0) > 0
    Xf = np.c_[Xf[:, :6][:, keepc], Xf[:, 6]]
    Af = Xf.T @ (Xf * ww[:, None])
    bf = np.linalg.solve(Af, Xf.T @ (ww * v[mm]))
    rf = v[mm] - Xf @ bf
    covf = np.linalg.inv(Af) * max(np.sum(ww * rf ** 2) / (mm.sum() - Xf.shape[1]), 1)
    pf['v_vs_vcorr_decFE'] = dict(slope=float(bf[-1]), err_slope=float(np.sqrt(covf[-1, -1])))
    pf['corr_vcorr_dec'] = float(np.corrcoef(cm['vcorr'][mm], dec[mm])[0, 1])
    out['per_fibre'] = pf
    print(json.dumps(pf, indent=1))
    # ---------------- 3. sky-line frame test (variance centroid)
    sky = {}
    single = (meta['max_mjd'] - meta['min_mjd'] < 0.5) & (cm['vspread'] < 1.0) & np.isfinite(cm['vcorr'])
    for nm, l0, half in (('OI5577', L['sky5577'], 2.6), ('HgI4358', L['HgI4358'], 2.0)):
        s = (WAVE > l0 - 12) & (WAVE < l0 + 12)
        lam_s = WAVE[s]
        with np.errstate(divide='ignore', invalid='ignore'):
            V = np.where(I[:, s] > 0, 1 / I[:, s], np.nan)
        core = np.abs(lam_s - l0) < half
        side = (np.abs(lam_s - l0) > 5) & (np.abs(lam_s - l0) < 11)
        bg = np.nanmedian(V[:, side], 1)
        ex = V[:, core] - bg[:, None]
        ex = np.where(ex > 0, ex, 0)
        cen = np.nansum(ex * lam_s[core][None, :], 1) / np.nansum(ex, 1)
        amp = np.nansum(ex, 1) / bg
        vcen = (cen / l0 - 1) * C
        good = single & np.isfinite(vcen) & (amp > (3.0 if nm == 'OI5577' else 0.3))
        X = np.c_[np.ones(good.sum()), cm['vcorr'][good]]
        bb, *_ = np.linalg.lstsq(X, vcen[good], rcond=None)
        r = vcen[good] - X @ bb
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        g2 = np.abs(r) < 5 * mad
        bb, *_ = np.linalg.lstsq(X[g2], vcen[good][g2], rcond=None)
        r2 = vcen[good][g2] - X[g2] @ bb
        cov = np.linalg.inv(X[g2].T @ X[g2]) * np.var(r2)
        t41 = cm['has_41485_41487'][good][g2]
        # stacked (median) normalised variance profiles in v_corr bins, Gaussian+constant with free width
        from scipy.optimize import curve_fit
        Vn = V / bg[:, None] - 1
        stk = []
        for lo_, hi_ in ((-25, -15), (-15, -5), (-5, 5), (5, 15), (15, 25), (25, 32)):
            mb = single & (cm['vcorr'] >= lo_) & (cm['vcorr'] < hi_)
            if mb.sum() < 20:
                continue
            prof = np.nanmedian(Vn[mb], 0)
            fg = lambda x, a, c0, sg, b0: a * np.exp(-0.5 * ((x - c0) / sg) ** 2) + b0  # noqa: E731
            pp, pc = curve_fit(fg, lam_s, prof, p0=[np.nanmax(prof), l0, 1.0, 0.0])
            stk.append((float(np.median(cm['vcorr'][mb])), float((pp[1] / l0 - 1) * C),
                        float(np.sqrt(pc[1, 1]) / l0 * C), float(abs(pp[2])), int(mb.sum())))
        stk = np.array(stk)
        ps = np.polyfit(stk[:, 0], stk[:, 1], 1)
        sky[nm + '_stacked'] = dict(bins=stk.tolist(), slope=float(ps[0]), intercept=float(ps[1]),
                                    note='[vcorr_median, centre_kms, err, sigma_A, n] per v_corr bin')
        sky[nm] = dict(n=int(g2.sum()), slope=float(bb[1]), err_slope=float(np.sqrt(cov[1, 1])),
                       intercept_kms=float(bb[0]), err_intercept=float(np.sqrt(cov[0, 0])), rms_kms=float(np.std(r2)),
                       vcorr_range=[float(cm['vcorr'][good].min()), float(cm['vcorr'][good].max())],
                       n_41485_41487=int(t41.sum()),
                       resid_41485_41487_mean=float(np.mean(r2[t41])) if t41.any() else None,
                       resid_41485_41487_err=float(np.std(r2[t41]) / np.sqrt(max(t41.sum(), 1))) if t41.any() else None,
                       note='variance-weighted centroid of the sky-line 1/ivar excess, single-night coadds')
        # no-correction hypothesis: slope 0 ; barycentric: slope 1
        out.setdefault('sky_arrays', {})[nm] = dict(vcorr=cm['vcorr'][good][g2].tolist(),
                                                    vcen=vcen[good][g2].tolist(), t41=t41.tolist())
        print(nm, sky[nm], flush=True)
    # simulate the moment-centroid estimator (Gaussian line on the 0.8 A grid, fixed window, side-band background)
    for nm, l0, half in (('OI5577', L['sky5577'], 2.6), ('HgI4358', L['HgI4358'], 2.0)):
        grid = 3600.0 + 0.8 * np.arange(7781)
        for lab, s_ in (('lsf', float(np.hypot(np.interp(l0, WAVE, sigA), 0.8 / np.sqrt(12)))),
                        ('stacked_width', float(np.median(np.array(sky[nm + '_stacked']['bins'])[:, 3])))):
            vs = np.linspace(-25, 32, 58)
            cs = []
            for vv in vs:
                mu = l0 * (1 + vv / C)
                g = grid[np.abs(grid - l0) < half]
                prof = np.exp(-0.5 * ((g - mu) / s_) ** 2)
                cs.append((np.sum(prof * g) / np.sum(prof) / l0 - 1) * C)
            sl = np.polyfit(vs, cs, 1)
            sky[nm][f'sim_slope_{lab}'] = float(sl[0])
            sky[nm][f'sim_sigma_{lab}'] = s_
        print(nm, {k: round(v, 3) for k, v in sky[nm].items() if k.startswith('sim')},
              'stacked slope', round(sky[nm + '_stacked']['slope'], 3))
    out['sky_frame_test'] = sky
    # ---------------- 4. [O III] fibres on tiles 41485/41487
    t41 = cm['has_41485_41487']
    near = cm['T'] > 0.1
    sig = base & (np.abs(snr) > 3)
    out['tiles_41485_41487'] = {'n_coadds_total': int(t41.sum()), 'n_with_good_o3': int((t41 & base).sum()),
                                'n_with_abs_snr_gt3': int((t41 & sig).sum()),
                                'n_near_trace_T_gt0.1': int((t41 & near).sum()),
                                'n_in_velocity_fits_abs_snr_gt5': int((t41 & base & (np.abs(snr) > 5) &
                                                                       np.isfinite(v)).sum()),
                                'ids': [str(x) for x in ids_lf[t41 & (sig | near)]],
                                'decs': [float(x) for x in dec[t41 & (sig | near)]]}
    # weight of those coadds in the all-fibre regression (sum u Z^2 share)
    u = 1 / EO ** 2
    share = np.nansum((u * Zf ** 2)[base & t41]) / np.nansum((u * Zf ** 2)[base])
    out['tiles_41485_41487']['share_of_instrument_weight'] = float(share)
    print(out['tiles_41485_41487'])
    sa = out.pop('sky_arrays')
    (D / 'verify_desi_velocity.json').write_text(json.dumps(out, indent=1))
    np.savez_compressed(D / 'verify_desi_velocity_sky.npz',
                        **{f'{k}_{q}': np.array(v_) for k, d in sa.items() for q, v_ in d.items()})


if __name__ == '__main__':
    main()
