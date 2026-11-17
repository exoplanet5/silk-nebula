"""w4_geostats task 2: ionisation order across the filament. Ridge offsets of GALEX FUV, NUV, DESI [O III], DESI [O II]
and eRASS1 0.2-0.6 keV relative to the NSNS [O III] ridge, with one estimator and placement + injection nulls.

Estimator (all bands): block-stacked cross profile (w4_geostats_common.stack_profile; 0.025-deg bins for maps,
0.1-deg bins and a single block for DESI fibres); reference = split-normal fit to the NSNS [O III] stack of the same
segment, giving the template T_N (unit peak, centred at the NSNS ridge x0_N). Offset delta_B = template shift of band
B against T_N(x - x0_N) (y = a T + b0 + b1 x, |x| <= 0.9 deg; DESI 1.2 deg). Variants: free-centre Gaussian
(minus the NSNS Gaussian centre), fit windows 0.6 / 0.9 / 1.2 deg. DESI [O II] is also measured against the DESI
[O III] template of the same fibres (the judge_skeptic definition). + = east (outside, away from Monogem).

Nulls:
  GALEX (placement + injection; three injections: co-shaped NSNS template [judge_obs definition], the band's own
    shape with its mode on the [O III] ridge [mode displacement], NSNS shape displaced by the prediction [power]): the trace is shifted rigidly to off-Silk sky (dRA +-0.8..3.0 deg in 0.2 steps x
    dDec -1.5/0/+1.5; pixels within 0.6 deg of the real spine excluded; >= 85 % finite bins required). At each
    placement a CO-SPATIAL copy of the observed band ridge (its own fitted split-normal shape and amplitude, centred
    on the placed [O III] ridge) is injected at pixel level and the estimator rerun -> null distribution of delta
    for a co-spatial truth with the real GALEX sky structure. A second injection at the shock-predicted offset
    (+0.065 deg, 50 % C IV mix) gives the power.
  NSNS ridge (placement): the NSNS stack + injected copy of itself at off-Silk placements (dRA +-1.5..5, dDec
    -2/0/+2; pixels within 1 deg excluded) -> scatter of x0_N itself, which enters every offset.
  DESI (fibres are confined to |x| < ~2 deg, so no off-Silk placements): (i) petal-group bootstrap (resample
    group_petal with replacement); (ii) petal-orientation injection null: fit a co-spatial model, keep each fibre's
    residual, mirror x for a random half of the petal groups, inject the co-spatial model at the new positions and
    re-measure (tests whether petal-coherent baselines can fake an offset); (iii) along-spine jackknife.
  eRASS1: no X-ray ridge. Ridge amplitude at the [O III] position (Gaussian FWHM 0.28 deg, linear baseline) and the
    best step location x_step (scan -2..+2 deg) with rotations of the trace about PSR B0656+14 (same radius;
    |theta| 20-180 deg in 10-deg steps) as the null.
Output: data/derived/w4_geostats_order.json, w4_geostats_order_profiles.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w4_geostats_order.py
"""
import sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w4_geostats_common import (CC0, Placed, stack_profile, fit_split, fit_gauss_centre, shift_fit, template_from_fit,  # noqa
                                splitnorm, shifted_cc, rotated_cc, load_nsns, load_galex, load_desi, load_erass,
                                empirical_p, robust_sd, DER, json_dump, radec_to_vec, s_at, SEGMENTS)

RNG = np.random.default_rng(41172)
OUT = DER / 'w4_geostats_order.json'
E25 = np.arange(-1.2, 1.2001, 0.025); X25 = 0.5 * (E25[1:] + E25[:-1])
E10 = np.arange(-1.5, 1.5001, 0.1); X10 = 0.5 * (E10[1:] + E10[:-1])
E01 = np.arange(-3.0, 3.0001, 0.1); X01 = 0.5 * (E01[1:] + E01[:-1])
PRED = {'FUV': 0.065, 'NUV': 0.035, 'O2': -0.045}
GSHIFTS = [(a, b) for a in np.r_[-3.0:-0.79:0.2, 0.8:3.01:0.2] for b in (-1.5, 0.0, 1.5)]
NSHIFTS = [(a, b) for a in np.r_[-5:-1.4:0.5, 1.5:5.01:0.5] for b in (-2.0, 0.0, 2.0)]


def tmpl_at(T, x0):
    return lambda xx: T(xx - x0)


def ref_nsns(N, seg):
    P = Placed(CC0, N, seg)
    st, B = stack_profile(P, N['v'][P.idx], E25, 'median')
    f = fit_split(X25, st)
    gc = fit_gauss_centre(X25, st)
    return P, st, B, f, gc


def measure(x, y, TN, x0N, gcN, xmax=0.9):
    d, a = shift_fit(x, y, tmpl_at(TN, x0N), xmax=xmax)
    gc = fit_gauss_centre(x, y, xmax=xmax)
    return dict(delta=d, amp=a, delta_gauss=(gc - gcN) if np.isfinite(gc) else np.nan)


# ------------------------------------------------------------------------------------------------ NSNS ridge scatter
def nsns_ridge_null(N, seg, fN):
    TN = template_from_fit(fN)
    out = []
    for dra, ddec in NSHIFTS:
        cc = shifted_cc(dra, ddec)
        P = Placed(cc, N, seg, excl_real=1.0)
        v = N['v'][P.idx] + fN['A'] * TN(P.x - fN['x0'])
        st, _ = stack_profile(P, v, E25, 'median')
        if np.isfinite(st[np.abs(X25) <= 0.9]).mean() < 0.85:
            continue
        f = fit_split(X25, st)
        if f is not None:
            out.append((f['x0'] - fN['x0'], f['ratio_we']))
    a = np.array(out)
    return dict(n=int(a.shape[0]), x0_sd=float(np.std(a[:, 0])), x0_rsd=robust_sd(a[:, 0]),
                x0_median=float(np.median(a[:, 0])))


# ------------------------------------------------------------------------------------------------ GALEX
def galex_band(G, band, seg, fN, gcN, stN):
    TN = template_from_fit(fN)
    P = Placed(CC0, G, seg)
    st, B = stack_profile(P, G['v'][P.idx], E25, 'mean', tile_demean=True, min_n=10)
    obs = {}
    for xm in (0.6, 0.9, 1.2):
        obs[f'xmax{xm}'] = measure(X25, st, TN, fN['x0'], gcN, xmax=xm)
    fB = fit_split(X25, st)
    # along-spine jackknife (10 chunks) of the primary statistic
    K = 10
    edges = np.linspace(0, P.nblock, K + 1).astype(int)
    jk = []
    for k in range(K):
        keep = np.ones(P.nblock, bool); keep[edges[k]:edges[k + 1]] = False
        with np.errstate(all='ignore'):
            y = np.nanmean(B[keep], axis=0)
        jk.append(shift_fit(X25, y, tmpl_at(TN, fN['x0']))[0])
    jk = np.array(jk)
    jk_err = float(np.sqrt((K - 1) / K * np.sum((jk - jk.mean()) ** 2)))
    # Three injections at each off-Silk placement (pixel level, observed band amplitude):
    #  'coshape'  : NSNS [O III] template shape, centred on the [O III] ridge  (judge_obs definition; tests whether the
    #               band is distributed like [O III] at all -> delta = shape + displacement)
    #  'ownshape' : the band's own fitted split-normal, its mode placed on the [O III] ridge (tests mode displacement)
    #  'pred'     : NSNS shape displaced by the shock prediction PRED[band]
    a_obs = obs['xmax0.9']['amp']
    if fB is not None:
        Tb = lambda xx: np.exp(-0.5 * (xx / np.where(xx < 0, fB['sig_w'], fB['sig_e'])) ** 2)
        A_own = fB['A']
    else:
        Tb = TN; A_own = a_obs
    inj = {'coshape': lambda xx: a_obs * TN(xx - fN['x0']),
           'ownshape': lambda xx: A_own * Tb(xx - fN['x0']),
           'pred': lambda xx: a_obs * TN(xx - fN['x0'] - PRED[band])}
    nul = {k: dict(delta=[], gauss=[], mode=[], ratio=[]) for k in inj}
    used = []
    for dra, ddec in GSHIFTS:
        cc = shifted_cc(dra, ddec)
        Pn = Placed(cc, G, seg, excl_real=0.6)
        if Pn.idx.size < 2000:
            continue
        base = G['v'][Pn.idx]
        st0, _ = stack_profile(Pn, base + inj['coshape'](Pn.x), E25, 'mean', tile_demean=True, min_n=10)
        if np.isfinite(st0[np.abs(X25) <= 0.9]).mean() < 0.85:
            continue
        used.append((dra, ddec))
        for k, f in inj.items():
            stk = st0 if k == 'coshape' else stack_profile(Pn, base + f(Pn.x), E25, 'mean', tile_demean=True, min_n=10)[0]
            m = measure(X25, stk, TN, fN['x0'], gcN)
            fs = fit_split(X25, stk)
            nul[k]['delta'].append(m['delta']); nul[k]['gauss'].append(m['delta_gauss'])
            nul[k]['mode'].append(fs['x0'] - fN['x0'] if fs else np.nan)
            nul[k]['ratio'].append(fs['ratio_we'] if fs else np.nan)
    nul = {k: {q: np.array(v) for q, v in d.items()} for k, d in nul.items()}
    d_obs = obs['xmax0.9']['delta']
    mode_obs = (fB['x0'] - fN['x0']) if fB else np.nan
    ratio_obs = fB['ratio_we'] if fB else np.nan

    def zrow(null, o, side):
        nn = null[np.isfinite(null)]
        med, sd = float(np.median(nn)), float(np.std(nn))
        return dict(obs=float(o), median=med, sd=sd, rsd=robust_sd(nn), z=float((o - med) / sd),
                    p_one=empirical_p(nn - med, o - med, side), p_two=empirical_p(nn - med, o - med, 'two'), n=int(nn.size))
    res = dict(obs=obs, fit=fB, jk_err=jk_err, n_null=len(used), placements=used,
               test_coshape_delta=zrow(nul['coshape']['delta'], d_obs, 'greater'),
               test_coshape_gauss=zrow(nul['coshape']['gauss'], obs['xmax0.9']['delta_gauss'], 'greater'),
               test_mode_displacement=zrow(nul['ownshape']['mode'], mode_obs, 'greater'),
               test_wing_ratio_we=zrow(nul['coshape']['ratio'], ratio_obs, 'less'),
               ownshape_delta=zrow(nul['ownshape']['delta'], d_obs, 'greater'),
               inj_pred=dict(delta_true=PRED[band], median=float(np.median(nul['pred']['delta'])),
                             sd=float(np.std(nul['pred']['delta'])),
                             frac_ge_obs=float(np.mean(nul['pred']['delta'] >= d_obs)),
                             mode_median=float(np.nanmedian(nul['pred']['mode'])), mode_sd=float(np.nanstd(nul['pred']['mode']))),
               null_values={k: {q: v.tolist() for q, v in d.items()} for k, d in nul.items()})
    s0 = res['test_coshape_delta']['sd']; s1 = res['inj_pred']['sd']; sp = 0.5 * (s0 + s1)
    res['LR_pred_vs_cospatial'] = float(np.exp(-0.5 * ((d_obs - res['inj_pred']['median']) / sp) ** 2
                                               + 0.5 * ((d_obs - res['test_coshape_delta']['median']) / sp) ** 2))
    ms = res['test_mode_displacement']['sd']
    res['LR_mode_pred_vs_cospatial'] = float(np.exp(-0.5 * ((mode_obs - PRED[band]) / ms) ** 2 + 0.5 * (mode_obs / ms) ** 2))
    return res, st


# ------------------------------------------------------------------------------------------------ DESI
def desi_offsets(Dd, seg, fN, gcN):
    TN = template_from_fit(fN)
    P = Placed(CC0, Dd, seg, maxd=1.8, block=100.0)
    x = P.x
    o3 = Dd['o3'][P.idx]; o2 = Dd['o2'][P.idx]; pet = Dd['petal'][P.idx]
    s = P.s

    def prof(xx, vv):
        st, _ = stack_profile(_FakeP(xx), vv, E10, 'mean', min_n=3)
        return st

    def stats(xx, v3, v2):
        p3 = prof(xx, v3); p2 = prof(xx, v2)
        f3 = fit_split(X10, p3, xmax=1.2)
        d3N = shift_fit(X10, p3, tmpl_at(TN, fN['x0']), xmax=1.2)
        d2N = shift_fit(X10, p2, tmpl_at(TN, fN['x0']), xmax=1.2)
        if f3 is not None:
            T3 = template_from_fit(f3)
            d23 = shift_fit(X10, p2, tmpl_at(T3, f3['x0']), xmax=1.2)
        else:
            d23 = (np.nan, np.nan)
        return dict(d_O3_vs_NSNS=d3N[0], a3=d3N[1], d_O2_vs_NSNS=d2N[0], a2=d2N[1], d_O2_vs_O3=d23[0],
                    O3_fit=f3, p3=p3, p2=p2)
    obs = stats(x, o3, o2)
    # (i) petal-group bootstrap
    ug, inv = np.unique(pet, return_inverse=True)
    bs = []
    for _ in range(300):
        pick = RNG.integers(0, ug.size, ug.size)
        idx = np.concatenate([np.nonzero(inv == g)[0] for g in pick])
        r = stats(x[idx], o3[idx], o2[idx])
        bs.append((r['d_O3_vs_NSNS'], r['d_O2_vs_NSNS'], r['d_O2_vs_O3']))
    bs = np.array(bs)
    # (ii) petal-orientation injection null (co-spatial truth)
    a3, a2 = obs['a3'], obs['a2']
    m3 = a3 * TN(x - fN['x0']); m2 = a2 * TN(x - fN['x0'])
    r3 = o3 - m3; r2 = o2 - m2
    nul = []; nulp = []
    for _ in range(300):
        sg = np.where(RNG.random(ug.size) < 0.5, -1.0, 1.0)[inv]
        xs = sg * x
        v3 = a3 * TN(xs - fN['x0']) + r3
        v2 = a2 * TN(xs - fN['x0']) + r2
        r = stats(xs, v3, v2)
        nul.append((r['d_O3_vs_NSNS'], r['d_O2_vs_NSNS'], r['d_O2_vs_O3']))
        v2p = a2 * TN(xs - fN['x0'] - PRED['O2']) + r2
        rp = stats(xs, v3, v2p)
        nulp.append(rp['d_O2_vs_O3'])
    nul = np.array(nul); nulp = np.array(nulp)
    # (iii) along-spine jackknife
    K = 8
    se = np.linspace(s.min(), s.max() + 1e-6, K + 1)
    jk = []
    for k in range(K):
        keep = ~((s >= se[k]) & (s < se[k + 1]))
        r = stats(x[keep], o3[keep], o2[keep])
        jk.append((r['d_O3_vs_NSNS'], r['d_O2_vs_NSNS'], r['d_O2_vs_O3']))
    jk = np.array(jk)
    jke = np.sqrt((K - 1) / K * np.nansum((jk - np.nanmean(jk, 0)) ** 2, 0))
    names = ('d_O3_vs_NSNS', 'd_O2_vs_NSNS', 'd_O2_vs_O3')
    res = dict(n_fibres=int(x.size), n_petal_groups=int(ug.size), O3_fit=obs['O3_fit'], a3=a3, a2=a2)
    for j, nm in enumerate(names):
        o = obs[nm]
        b = bs[:, j]; nn = nul[:, j]
        res[nm] = dict(obs=o, boot_sd=float(np.nanstd(b)), boot_p16_p84=(float(np.nanpercentile(b, 16)), float(np.nanpercentile(b, 84))),
                       orient_null_median=float(np.nanmedian(nn)), orient_null_sd=float(np.nanstd(nn)),
                       orient_p_two=empirical_p(nn - np.nanmedian(nn), o - np.nanmedian(nn), 'two'),
                       z_orient=float((o - np.nanmedian(nn)) / np.nanstd(nn)), jk_err=float(jke[j]),
                       z_boot=float(o / np.nanstd(b)))
    res['d_O2_vs_O3']['inj_pred'] = dict(delta_true=PRED['O2'], median=float(np.nanmedian(nulp)), sd=float(np.nanstd(nulp)),
                                         frac_ge_obs=float(np.mean(nulp >= obs['d_O2_vs_O3'])))
    return res, obs['p3'], obs['p2']


class _FakeP:
    """Single-block 'Placed' for fibre lists."""
    def __init__(self, x):
        self.x = x; self.nblock = 1; self.blk = np.zeros(x.size, int); self.tile = None


# ------------------------------------------------------------------------------------------------ eRASS1
def xray(E, seg, fN):
    def prof(cc):
        P = Placed(cc, E, seg, maxd=3.2, block=100.0)
        st, _ = stack_profile(P, E['v'][P.idx], E01, 'mean', min_n=3)
        return st

    def stats(st):
        m = np.isfinite(st) & (np.abs(X01) <= 2.5)
        if m.sum() < 0.7 * (np.abs(X01) <= 2.5).sum():
            return None
        x, y = X01[m], st[m]
        g = np.exp(-0.5 * ((x - fN['x0']) / (0.28 / 2.3548)) ** 2)
        A = np.linalg.lstsq(np.vstack([g, np.ones_like(x), x]).T, y, rcond=None)[0][0]
        best = (-np.inf, np.nan)
        for xs in np.arange(-2.0, 2.01, 0.1):
            c = np.linalg.lstsq(np.vstack([np.ones_like(x), x, (x < xs).astype(float)]).T, y, rcond=None)[0]
            if c[2] > best[0]:
                best = (c[2], xs)
        return dict(ridge_amp=float(A), ridge_frac=float(A / np.mean(y)), x_step=float(best[1]),
                    step=float(best[0]), level=float(np.mean(y)))
    st = prof(CC0)
    obs = stats(st)
    nul = []
    for th in np.r_[-180:-19:10, 20:181:10]:
        r = stats(prof(rotated_cc(th, 'PSR')))
        if r is not None:
            nul.append((r['ridge_amp'], r['x_step'], r['ridge_frac'], r['step']))
    nul = np.array(nul)
    return dict(obs=obs, n_null=int(nul.shape[0]),
                ridge_amp_null=dict(median=float(np.median(nul[:, 0])), sd=float(np.std(nul[:, 0])),
                                    z=float((obs['ridge_amp'] - np.median(nul[:, 0])) / np.std(nul[:, 0]))),
                ridge_frac_null_sd=float(np.std(nul[:, 2])),
                step_null=dict(median=float(np.median(nul[:, 3])), sd=float(np.std(nul[:, 3])),
                               p_ge_obs=empirical_p(nul[:, 3], obs['step'], 'greater')),
                x_step_null=dict(median=float(np.median(nul[:, 1])), sd=float(np.std(nul[:, 1])),
                                 frac_within_0p3_of_obs=float(np.mean(np.abs(nul[:, 1] - obs['x_step']) <= 0.3))),
                x_step_null_values=nul[:, 1].tolist()), st


def main():
    N = load_nsns()
    G = {b: load_galex(b) for b in ('FUV', 'NUV')}
    Dd = load_desi()
    E = load_erass()
    res = dict(meta=dict(desi_fibres_ok=Dd['n_ok'], desi_fibres_all=Dd['n_all'], pred=PRED,
                         n_galex_placements_tried=len(GSHIFTS)))
    prof_store = {}
    for seg in ('S1', 'S12', 'N'):
        Pn, stN, BN, fN, gcN = ref_nsns(N, seg)
        r = dict(nsns_fit=fN, nsns_gauss_centre=gcN)
        prof_store[f'{seg}_nsns'] = stN
        r['nsns_ridge_null'] = nsns_ridge_null(N, seg, fN)
        print(seg, 'NSNS', fN, r['nsns_ridge_null'], flush=True)
        if seg != 'N':
            for b in ('FUV', 'NUV'):
                rb, stb = galex_band(G[b], b, seg, fN, gcN, stN)
                r[b] = rb
                prof_store[f'{seg}_{b}'] = stb
                print(seg, b, rb['obs']['xmax0.9'], 'coshape', rb['test_coshape_delta'], 'mode', rb['test_mode_displacement'],
                      'wing', rb['test_wing_ratio_we'], 'pred', rb['inj_pred'], 'jk', rb['jk_err'], flush=True)
            rx, stx = xray(E, seg, fN)
            r['xray'] = rx; prof_store[f'{seg}_xray'] = stx
            print(seg, 'xray', rx['obs'], rx['ridge_amp_null'], rx['x_step_null'], flush=True)
        rd, p3, p2 = desi_offsets(Dd, seg, fN, gcN)
        r['desi'] = rd
        prof_store[f'{seg}_desi_o3'] = p3; prof_store[f'{seg}_desi_o2'] = p2
        print(seg, 'DESI', {k: rd[k] for k in ('n_fibres', 'd_O3_vs_NSNS', 'd_O2_vs_O3')}, flush=True)
        res[seg] = r
    np.savez_compressed(DER / 'w4_geostats_order_profiles.npz', x25=X25, x10=X10, x01=X01, **prof_store)
    json_dump(res, OUT)
    print('saved', OUT)


if __name__ == '__main__':
    main()
