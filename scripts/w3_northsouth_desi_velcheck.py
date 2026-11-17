"""w3_northsouth_desi_velcheck.py -- wave-3 (tag northsouth): robustness of the DESI north/south [O III]
velocity difference to the pixel weighting and line choice.

For the two halves (trace Dec < / >= 19.5) and the three zones (< 17.6, 17.6-21, >= 21), joint per-pixel regression
of the continuum-subtracted residual (prep() of verify_desi_audit_vel) on the Dec-zone templates (same templates as
w3_northsouth_desi.py), then the zone coefficient spectrum is fitted with:
  W1  both lines, weights from the nside-32 block-bootstrap pixel sigma (as verify_desi_audit_vel)
  W2  both lines, weights from the formal per-pixel WLS sigma
  W3  both lines, uniform weights
  L5  [O III]5008 only (W1 weights);  L4  [O III]4960 only (W1 weights)
  width test: chi^2 of the doublet fit for intrinsic FWHM 28 (baseline), 40, 55, 70 km/s (W1 weights)
  MOM model-free first moment of the coefficient spectrum within +-1.4 A of each line (continuum from 4-12 A
      side bands), photon-weighted, both lines combined
Each estimator is bootstrapped over the same nside-32 resamples (velocity differences north - south are computed
per resample, so their errors include the covariance).
Output: data/derived/w3_northsouth_desi_velcheck.json
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, WAVE, C, L, load_spectra, lsf_sigma, Exposures, Coadds, Template  # noqa: E402
from verify_desi_audit import healpix_ring  # noqa: E402
from verify_desi_audit_vel import prep  # noqa: E402
from w3_northsouth_desi import petal_sub, wls, VelFitter, summ  # noqa: E402

RNG = np.random.default_rng(777)
NBOOT = 300


def moment_v(lam, b, lam0s, half=1.4):
    num, den = 0.0, 0.0
    for l0 in lam0s:
        c = l0 * (1 + 12.0 / C)
        side = (np.abs(lam - c) > 4) & (np.abs(lam - c) < 12)
        p = np.polyfit(lam[side] - c, b[side], 1)
        inn = np.abs(lam - c) < half
        f = b[inn] - np.polyval(p, lam[inn] - c)
        num += np.sum(f * (lam[inn] / l0 - 1) * C)
        den += np.sum(f)
    return num / den if den > 0 else np.nan


def main():
    meta = np.load(D / 'verify_desi_coadd_meta.npz')
    ra, dec, z = meta['ra'], meta['dec'], meta['redshift']
    ok = meta['ok_contrib'] & np.isfinite(meta['Z']) & np.isfinite(meta['f_chrom_5008'])
    co = Coadds({k: meta[k] for k in ('targetid', 'survey', 'program', 'ra', 'dec')}, Exposures())
    T = Template()
    zones = {'halves': np.array([0.0, 19.5, 30.0]), 'three': np.array([12.0, 17.6, 21.0, 28.6])}
    Zz = {k: petal_sub(T, co, ra, dec, e)[0] for k, e in zones.items()}
    ids, F, I, M, src = load_spectra()
    pos = {s: i for i, s in enumerate(ids)}
    take = np.array([pos[s] for s in meta['sparcl_id']])
    R = F[take] - M[take]
    I = I[take]
    del F, M
    lsf, _ = lsf_sigma()
    lam0s = [L['OIII4959'], L['OIII5007']]
    sig = [float(np.sqrt(np.interp(l0, WAVE, lsf) ** 2 + (12.0 / C * l0) ** 2)) for l0 in lam0s]
    keep, Y, w, pix, seg = prep(R, I, z, ok, lam0s, meta=meta)
    lam = WAVE[pix]
    vg = np.arange(-40, 70.01, 0.25)
    fit_both = VelFitter(lam, seg, lam0s, [1 / 2.98, 1.0], sig, vg)
    m5, m4 = seg == 1, seg == 0
    fit_5 = VelFitter(lam[m5], np.zeros(m5.sum(), int), [lam0s[1]], [1.0], [sig[1]], vg)
    fit_4 = VelFitter(lam[m4], np.zeros(m4.sum(), int), [lam0s[0]], [1.0], [sig[0]], vg)
    wide = [VelFitter(lam, seg, lam0s, [1 / 2.98, 1.0],
                      [float(np.sqrt(s ** 2 + ((fw / 2.3548) ** 2 - 12.0 ** 2) / C ** 2 * l0 ** 2))
                       for s, l0 in zip(sig, lam0s)], vg) for fw in (40.0, 55.0, 70.0)]
    hp = healpix_ring(ra[keep], dec[keep], 32)
    ub, inv = np.unique(hp, return_inverse=True)
    mem = [np.where(inv == b)[0] for b in range(ub.size)]
    bidx = [np.concatenate([mem[b] for b in RNG.integers(0, ub.size, ub.size)]) for _ in range(NBOOT)]
    out = {}
    for zn, Zs in Zz.items():
        X = np.c_[np.ones(keep.size), *[Zk[keep] for Zk in Zs]]
        XtW = X.T * w
        Minv = np.linalg.inv(XtW @ X)
        B = Minv @ (XtW @ Y)
        resid = Y - X @ B
        s2 = np.sum(w[:, None] * resid ** 2, 0) / (keep.size - X.shape[1])
        boots = np.array([wls(X[ii], Y[ii], w[ii]) for ii in bidx])
        res = {}
        est = {}
        for k in range(len(Zs)):
            b = B[k + 1]
            q = np.percentile(boots[:, k + 1, :], [16, 84], axis=0)
            sb = np.maximum(0.5 * (q[1] - q[0]), 1e-30)
            sf = np.sqrt(Minv[k + 1, k + 1] * s2)
            su = np.full(b.size, np.median(sb))
            ests = {
                'W1_boot_weights': lambda bb: fit_both.fit(bb, sb)[0],
                'W2_formal_weights': lambda bb: fit_both.fit(bb, sf)[0],
                'W3_uniform': lambda bb: fit_both.fit(bb, su)[0],
                'L5_5008_only': lambda bb: fit_5.fit(bb[m5], sb[m5])[0],
                'L4_4960_only': lambda bb: fit_4.fit(bb[m4], sb[m4])[0],
                'MOM_moment': lambda bb: moment_v(lam, bb, lam0s),
            }
            r = {}
            for nm, fn in ests.items():
                v0 = fn(b)
                vb = np.array([fn(bb[k + 1]) for bb in boots])
                r[nm] = dict(v=float(v0), boot=summ(vb))
                est.setdefault(nm, []).append(vb)
            chis = [wf.fit(b, sb)[3].min() for wf in wide]
            chi0 = fit_both.fit(b, sb)[3].min()
            r['width_test_chi2'] = dict(FWHM_kms=[28.3, 40.0, 55.0, 70.0],
                                        chi2=[float(chi0)] + [float(c) for c in chis])
            r['ratio_sb_over_sf_median'] = float(np.median(sb / sf))
            res[f'zone{k}'] = r
        # north - south differences per resample (last zone minus first zone)
        diffs = {nm: summ(v[-1] - v[0]) for nm, v in est.items()}
        res['diff_last_minus_first'] = diffs
        out[zn] = res
        print(zn)
        for k in range(len(Zs)):
            print('  zone', k, {nm: (round(v['v'], 1), round(v['boot']['sig'], 1)) for nm, v in res[f'zone{k}'].items()
                                if isinstance(v, dict) and 'v' in v}, 'sb/sf', round(res[f'zone{k}']['ratio_sb_over_sf_median'], 2))
        print('  diff', {nm: (round(v['p50'], 1), round(v['sig'], 1)) for nm, v in diffs.items()})
    json.dump(out, open(D / 'w3_northsouth_desi_velcheck.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
