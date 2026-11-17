"""w3_oii_desi_fit.py -- (wave 3, tag oii) doublet fits to the DESI regression spectra of w3_oii_desi_regress.py.

For each subset / weighting (geometry instrument 'geo', neighbour-Halpha instrument 'ha'):
  [O III] 4959+5007 (ratio 1/2.98 fixed): free v and sigma (point-sampled Gaussians, grid) -> v_O3, sigma_O3,
      sigma_int = sqrt(sigma_O3^2 - sigma_R(5008)^2)  (ha instrument: v = +12, sigma_int = 12 km/s assumed).
  Doublet models (point-sampled unit Gaussians + linear baseline; [O II] fit 3712-3745 A, [S II] 6700-6752 A):
      w2        v = +14 km/s, sigma^2 = sigma_R^2 + 0.8^2/12 + (12 km/s)^2      (wave-2 per-fibre fitter replica)
      kin       v = v_O3(subset), sigma^2 = sigma_R^2 + sigma_int^2               (same kinematics as [O III])
      kin_s+-5/10, kin_v+-5   kin with sigma x0.90..1.10 or v +-5 km/s
      vfree     v free, sigma as kin;      vsfree   v and sigma free
  Ratio intervals: every bootstrap replicate refitted with the same model (free parameters re-optimised).
  Profile likelihood: chi^2(ratio) minimised over v, sigma, amplitude and baseline with the shrunk bootstrap
      covariance (and the same with kin kinematics fixed).
  [O II] 7320+7330 / 3726+3729 from Gaussian fits with kin kinematics.
Output: data/derived/w3_oii_desi_fits.json  (argument 'robust': fits of w3_oii_desi_regspec_robust.npz ->
        data/derived/w3_oii_desi_fits_robust.json)
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_oii_common import L, C, design, shrink_cov, grid_fit, chi2_curve, profile  # noqa: E402

D = '/Users/mickey/astronomy/silk-nebula/data/derived/'
TAG = sys.argv[1] if len(sys.argv) > 1 else ''          # '' or 'robust'
Z = np.load(D + ('w3_oii_desi_regspec_robust.npz' if TAG == 'robust' else 'w3_oii_desi_regspec.npz'))
LSF = json.load(open(D + 'w3_oii_desi_lsf.json'))
SR = {'3727': LSF['3727.1']['median'], '5008': LSF['5008.2']['median'],
      '6724': 0.5 * (LSF['6718.3']['median'] + LSF['6732.7']['median']),
      '7325': 0.5 * (LSF['7321.0']['median'] + LSF['7331.4']['median'])}
O2 = [L['OII3726'], L['OII3729']]
O3 = [L['OIII4959'], L['OIII5007']]
S2 = [L['SII6731'], L['SII6716']]   # order chosen so that ratio = line2/line1 = 6716/6731
O7 = [L['OII7319'], L['OII7320'], L['OII7330'], L['OII7331']]
VGRID = np.arange(-10.0, 40.01, 0.5)
FR = np.arange(0.80, 1.40001, 0.01)     # sigma / sigma_R
RGRID = np.arange(0.80, 1.901, 0.01)


def get(inst, sn, wn, wt):
    k = f'S_{inst}_{sn}_{wn}_{wt}'
    if k not in Z.files:
        return None
    lam, S, bs = Z[f'lam_{wn}'], Z[k], Z[f'bs_{inst}_{sn}_{wn}_{wt}'].astype(float)
    return lam, S, bs


def prep(lam, S, bs, lo, hi):
    """Diagonal weights from the bootstrap variance after removing the common line-strength (scale) mode, smoothed
    over 5 pixels. (A GLS with the bootstrap covariance would be fitted to the very replicates used to estimate the
    intervals and shrink them artificially; the block bootstrap also has rank <= n_blocks.)"""
    sel = (lam >= lo) & (lam <= hi)
    b = bs[:, sel]
    b = b[np.all(np.isfinite(b), 1)]
    y = S[sel]
    a = (b @ y) / (y @ y)
    var = np.var(b - a[:, None] * y[None, :], 0)
    k = 2
    vs = np.array([np.median(var[max(0, i - k):i + k + 1]) for i in range(var.size)])
    Wi = np.diag(1.0 / np.maximum(vs, 1e-12))
    return lam[sel], np.vstack([y[None, :], b]), Wi


def pair_design(x, lines, v, s, ratio=None, npoly=2):
    X = design(x, lines, v, [s] * len(lines), npoly=npoly)
    if ratio is None:
        return X
    return np.vstack([X[0] + ratio * X[1], X[2:]])


def o3_kin(lam, S, bs):
    x, Ys, Wi = prep(lam, S, bs, 4945.0, 5025.0)

    def builders(vg, fg):
        for v in vg:
            for f in fg:
                s = f * SR['5008']
                X = design(x, O3, v, [s * O3[0] / O3[1], s], npoly=2)
                yield (v, s), np.vstack([X[0] / 2.98 + X[1], X[2], X[3]])
    c2, par, p = grid_fit(Ys, Wi, builders(VGRID, FR))
    v, s = par[0]
    vs = np.array([q[0] for q in par[1:]])
    ss = np.array([q[1] for q in par[1:]])
    sint2 = s ** 2 - SR['5008'] ** 2
    return dict(v=float(v), v_p16=float(np.percentile(vs, 16)), v_p84=float(np.percentile(vs, 84)),
                sigma=float(s), sigma_p16=float(np.percentile(ss, 16)), sigma_p84=float(np.percentile(ss, 84)),
                sigma_int_A=float(np.sqrt(max(sint2, 0.0))), sigma_int_kms=float(np.sqrt(max(sint2, 0.0)) / 5008.2 * C),
                amp5007=float(p[0][0]), chi2=float(c2[0]), npix=int(x.size))


def summarise(ratios):
    q = np.nanpercentile(ratios[1:], [2.5, 16, 50, 84, 97.5])
    return dict(ratio=float(ratios[0]), p2p5=float(q[0]), p16=float(q[1]), p50=float(q[2]), p84=float(q[3]),
                p97p5=float(q[4]))


def doublet_models(lam, S, bs, lines, sR, v_kin, sint_A, lo, hi):
    x, Ys, Wi = prep(lam, S, bs, lo, hi)
    lc = 0.5 * (lines[0] + lines[1])
    s_kin = float(np.sqrt(sR ** 2 + (sint_A * lc / 5008.2) ** 2))
    s_w2 = float(np.sqrt(sR ** 2 + 0.8 ** 2 / 12 + (lc * 12 / C) ** 2))
    specs = {'w2': ([14.0], [s_w2]), 'kin': ([v_kin], [s_kin])}
    for tag, f in (('kin_s-10', 0.90), ('kin_s-5', 0.95), ('kin_s+5', 1.05), ('kin_s+10', 1.10)):
        specs[tag] = ([v_kin], [s_kin * f])
    for tag, dv in (('kin_v-5', -5.0), ('kin_v+5', 5.0)):
        specs[tag] = ([v_kin + dv], [s_kin])
    specs['vfree'] = (VGRID, [s_kin])
    specs['vsfree'] = (VGRID, list(FR * sR))
    out = {}
    for tag, (vg, sg) in specs.items():
        c2, par, p = grid_fit(Ys, Wi, (((v, s), pair_design(x, lines, v, s)) for v in vg for s in sg))
        ratios = np.array([q[1] / q[0] for q in p])
        d = summarise(ratios)
        d.update(v=float(par[0][0]), sigma=float(par[0][1]), sigma_over_sR=float(par[0][1] / sR), chi2=float(c2[0]),
                 npix=int(x.size), sum_amp=float(p[0][0] + p[0][1]),
                 v_p16=float(np.percentile([q[0] for q in par[1:]], 16)),
                 v_p84=float(np.percentile([q[0] for q in par[1:]], 84)),
                 s_p16=float(np.percentile([q[1] for q in par[1:]], 16) / sR),
                 s_p84=float(np.percentile([q[1] for q in par[1:]], 84) / sR))
        out[tag] = d
    # white-noise Monte Carlo around the kin best fit (noise = per-pixel weights) -> independent interval
    rng = np.random.default_rng(7)
    Xk = pair_design(x, lines, v_kin, s_kin)
    pk = np.linalg.solve(Xk @ Wi @ Xk.T, Xk @ Wi @ Ys[0])
    sd = 1 / np.sqrt(np.diag(Wi))
    mc = rng.normal(0, 1, (400, x.size)) * sd[None, :] + (pk @ Xk)[None, :]
    c2m, parm, pm = grid_fit(mc, Wi, [((v_kin, s_kin), Xk)])
    rm = np.array([q[1] / q[0] for q in pm])
    out['kin']['mc_p16'], out['kin']['mc_p84'] = float(np.percentile(rm, 16)), float(np.percentile(rm, 84))
    out['kin']['mc_p2p5'], out['kin']['mc_p97p5'] = float(np.percentile(rm, 2.5)), float(np.percentile(rm, 97.5))
    out['kin']['chi2_dof'] = float(out['kin']['chi2'] / max(x.size - 4, 1))
    # profile likelihood in the ratio (data row only)
    y0 = Ys[:1]
    prof = []
    for r in RGRID:
        Xs = [pair_design(x, lines, v, f * sR, ratio=r) for v in VGRID[::2] for f in FR[::3]]
        prof.append(chi2_curve(y0, Wi, Xs).min())
    prof = np.array(prof)
    scale = max(prof.min() / max(x.size - 5, 1), 1.0)       # inflate if chi2/dof > 1
    d = (prof - prof.min()) / scale
    prof_k = chi2_curve(y0, Wi, [pair_design(x, lines, v_kin, s_kin, ratio=r) for r in RGRID])[:, 0]
    dk = (prof_k - prof_k.min()) / max(prof_k.min() / max(x.size - 4, 1), 1.0)

    def summ(dd):
        return dict(r_best=float(RGRID[int(np.argmin(dd))]), r68=[float(RGRID[dd <= 1].min()), float(RGRID[dd <= 1].max())],
                    r95=[float(RGRID[dd <= 4].min()), float(RGRID[dd <= 4].max())],
                    dchi2_at_1p46=float(np.interp(1.46, RGRID, dd)), dchi2_at_1p30=float(np.interp(1.30, RGRID, dd)),
                    dchi2=[float(t) for t in dd])
    out['profile'] = summ(d)
    out['profile_kin'] = summ(dk)
    out['rgrid'] = [float(t) for t in RGRID]
    out['s_kin'] = s_kin
    out['s_w2'] = s_w2
    return out


def total_flux(lam, S, bs, lines, sR, v, sint_A, lo, hi):
    x, Ys, Wi = prep(lam, S, bs, lo, hi)
    s = float(np.sqrt(sR ** 2 + (sint_A * np.mean(lines) / 5008.2) ** 2))
    X = design(x, lines, v, [s] * len(lines), npoly=2)
    c2, par, p = grid_fit(Ys, Wi, [((v, s), X)])
    return np.array([q[:len(lines)].sum() for q in p])


def main():
    res = {'sigma_R': SR, 'notes': 'ratios 3729/3726 and 6716/6731 (energy); point-sampled Gaussians; '
                                   'bootstrap = HEALPix nside-32 blocks'}
    subs_geo = sorted({k[len('S_geo_'):].rsplit('_', 2)[0] for k in Z.files if k.startswith('S_geo_')})
    for wt in ('F', 'P'):
        for sn in subs_geo:
            g5, g3 = get('geo', sn, '5007', wt), get('geo', sn, '3727', wt)
            if g5 is None or g3 is None or not np.all(np.isfinite(g5[1])) or not np.all(np.isfinite(g3[1])):
                continue
            key = f'geo|{sn}|{wt}'
            kin = o3_kin(*g5)
            r = {'o3': kin, 'n': int(Z[f'n_geo_{sn}']), 'nblk': int(Z[f'nblk_geo_{sn}'])}
            r['o2'] = doublet_models(*g3, O2, SR['3727'], kin['v'], kin['sigma_int_A'], 3712.0, 3745.0)
            gs = get('geo', sn, 'SII', wt)
            if gs is not None and np.all(np.isfinite(gs[1])):
                r['s2'] = doublet_models(*gs, S2, SR['6724'], kin['v'], kin['sigma_int_A'], 6700.0, 6752.0)
            g7 = get('geo', sn, '7325', wt)
            if g7 is not None and np.all(np.isfinite(g7[1])):
                t7 = total_flux(*g7, O7, SR['7325'], kin['v'], kin['sigma_int_A'], 7300.0, 7352.0)
                t3 = total_flux(*g3, O2, SR['3727'], kin['v'], kin['sigma_int_A'], 3712.0, 3745.0)
                n = min(len(t7), len(t3))
                rr = t7[:n] / t3[:n]
                d = summarise(rr)
                d.update(o2_3727_over_o3=float(t3[0]), o2_7325_over_o3=float(t7[0]))
                r['o2_7325_over_3727'] = d
            res[key] = r
            o = r['o2']
            print(f"{key:24s} n={r['n']:5d} blk={r['nblk']:3d} vO3={kin['v']:+5.1f} sint={kin['sigma_int_kms']:4.1f} | "
                  f"w2 {o['w2']['ratio']:.3f} [{o['w2']['p16']:.3f},{o['w2']['p84']:.3f}] "
                  f"kin {o['kin']['ratio']:.3f} [{o['kin']['p16']:.3f},{o['kin']['p84']:.3f}] "
                  f"vfree {o['vfree']['ratio']:.3f}(v{o['vfree']['v']:+.0f}) "
                  f"vsfree {o['vsfree']['ratio']:.3f}(v{o['vsfree']['v']:+.0f},s{o['vsfree']['sigma_over_sR']:.2f}) "
                  f"prof {o['profile']['r_best']:.2f} {o['profile']['r68']} d146={o['profile']['dchi2_at_1p46']:.1f}"
                  + (f" | SII {r['s2']['kin']['ratio']:.2f} [{r['s2']['kin']['p16']:.2f},{r['s2']['kin']['p84']:.2f}]"
                     if 's2' in r else '')
                  + (f" | 7325/3727 {r['o2_7325_over_3727']['ratio']:.4f}" if 'o2_7325_over_3727' in r else ''),
                  flush=True)
    sint = 12.0 / C * 5008.2
    for wt in ('F', 'P'):
        for sn in ('all', 'no_patch', 'south', 'north', 'patch'):
            g3, gs = get('ha', sn, '3727', wt), get('ha', sn, 'SII', wt)
            if g3 is None:
                continue
            key = f'ha|{sn}|{wt}'
            r = {'n': int(Z[f'n_ha_{sn}']), 'nblk': int(Z[f'nblk_ha_{sn}'])}
            if np.all(np.isfinite(g3[1])):
                r['o2'] = doublet_models(*g3, O2, SR['3727'], 12.0, sint, 3712.0, 3745.0)
            if gs is not None and np.all(np.isfinite(gs[1])):
                r['s2'] = doublet_models(*gs, S2, SR['6724'], 12.0, sint, 6700.0, 6752.0)
            res[key] = r
            s1 = ('O2 kin %.3f [%.3f,%.3f] vfree %.3f(v%+.0f) prof %s' % (
                r['o2']['kin']['ratio'], r['o2']['kin']['p16'], r['o2']['kin']['p84'], r['o2']['vfree']['ratio'],
                r['o2']['vfree']['v'], r['o2']['profile']['r68'])) if 'o2' in r else ''
            s2 = ('SII kin %.3f [%.3f,%.3f] prof %s' % (r['s2']['kin']['ratio'], r['s2']['kin']['p16'],
                                                       r['s2']['kin']['p84'], r['s2']['profile']['r68'])) if 's2' in r else ''
            print(f'{key:24s} n={r["n"]} blk={r["nblk"]}', s1, '|', s2, flush=True)
    json.dump(res, open(D + ('w3_oii_desi_fits_robust.json' if TAG == 'robust' else 'w3_oii_desi_fits.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
