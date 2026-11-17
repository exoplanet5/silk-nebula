"""w3_northsouth_rm.py -- wave-3 (tag northsouth): is the NVSS rotation-measure step across the northern Silk
(spectra_ism: E(2-6 deg) - W(2-6 deg) median RM = +36.8 +- 4.6 rad m^-2 for spine Dec 20-28) tied to the filament?

Input: data/derived/spectra_ism_rm_sources.csv (Taylor+2009 NVSS RMs in a 20-deg cone on (116, +20); wave-1 table,
read only).  Own geometry code (no astropy matching):
  * the spine (silk_trace.json) is densified; each source gets a signed perpendicular offset (+ = the side the
    spine's local east normal points to) and the along-spine Dec of its nearest spine sample;
  * statistic S = median RM(+2..+6 deg) - median RM(-6..-2 deg) for sources whose nearest sample lies at
    original-spine Dec 20-28.2 (north) or 12-20 (south), not beyond the ends;  also on the wave-1 residual RM
    (quadratic foreground removed) and an L1 step location.
Placement nulls (whole trace moved rigidly, same along-spine selection):
  (i) translations only (orientation kept): 1-10 deg in random directions;
  (ii) random spin about the cone centre (0-360 deg) plus 0-8 deg translation;
  kept if both side bands hold >= 15 sources and the placed north-half side bands lie inside the 20-deg cone.
  A placement 'overlaps' the true one if > 30 % of its side-band sources are also true side-band sources.
Also: the E-W difference predicted by a smooth (planar / quadratic) RM gradient fitted to sources >= 7 deg from
the northern spine.
Outputs: data/derived/w3_northsouth_rm.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_northsouth_rm.py
"""
import csv
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import unitvec, rotmat  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula'
RNG = np.random.default_rng(2026)
CONE = (116.0, 20.0, 20.0)


def load_sources():
    rows = list(csv.DictReader(open(f'{ROOT}/data/derived/spectra_ism_rm_sources.csv')))
    g = lambda k: np.array([float(r[k]) for r in rows])  # noqa: E731
    return dict(ra=g('ra'), dec=g('dec'), rm=g('RM'), erm=g('eRM'), resid=g('resid_RM'), dperp_w1=g('dperp_deg'),
                dsp_w1=g('dec_spine'), beyond_w1=g('beyond').astype(bool))


def dense_spine(step=0.02):
    tr = json.load(open(f'{ROOT}/data/derived/silk_trace.json'))
    ra = np.array([p['ra'] for p in tr['spine']])
    de = np.array([p['dec'] for p in tr['spine']])
    v = unitvec(ra, de)
    pts = []
    for k in range(len(ra) - 1):
        seg = np.degrees(np.arccos(np.clip(v[k] @ v[k + 1], -1, 1)))
        nn = max(int(np.ceil(seg / step)), 1)
        f = np.arange(nn) / nn
        pts.append(v[k][None] * (1 - f[:, None]) + v[k + 1][None] * f[:, None])
    pts.append(v[-1:])
    P = np.concatenate(pts)
    P /= np.linalg.norm(P, axis=1)[:, None]
    # local east normal: tangent x position gives a vector perpendicular to the spine on the sphere
    tg = np.gradient(P, axis=0)
    nrm = np.cross(P, tg)                 # for a S->N spine this points to smaller RA (west); flip below
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    east = np.c_[-np.sin(np.radians(ra.mean())), np.cos(np.radians(ra.mean())), 0.0][0]
    if np.mean(nrm @ east) < 0:
        nrm = -nrm
    sdec = np.degrees(np.arcsin(P[:, 2]))
    return P, nrm, sdec


def coords(S, P, nrm, sdec, R=None):
    """Signed offset (deg) and along-spine (original) Dec for the trace rotated by R."""
    X = unitvec(S['ra'], S['dec'])
    Pr = P if R is None else P @ R.T
    Nr = nrm if R is None else nrm @ R.T
    dots = np.clip(X @ Pr.T, -1, 1)
    j = np.argmax(dots, 1)
    d = np.degrees(np.arccos(dots[np.arange(len(j)), j]))
    sgn = np.sign(np.sum((X - Pr[j]) * Nr[j], 1))
    sgn[sgn == 0] = 1
    beyond = (j == 0) | (j == len(P) - 1)
    return sgn * d, sdec[j], beyond, Pr


def step_stat(val, dp, ds, beyond, dlo, dhi, inner=2.0, outer=6.0):
    m = (~beyond) & (ds >= dlo) & (ds < dhi)
    e = m & (dp > inner) & (dp < outer)
    w = m & (dp < -inner) & (dp > -outer)
    if e.sum() < 15 or w.sum() < 15:
        return None, e, w
    return float(np.median(val[e]) - np.median(val[w])), e, w


def l1_step(dd, rr):
    grid = np.arange(-4.0, 4.01, 0.1)
    cost = []
    for d0 in grid:
        a, b = rr[dd < d0], rr[dd >= d0]
        if len(a) < 5 or len(b) < 5:
            cost.append(np.inf)
            continue
        cost.append(np.sum(np.abs(a - np.median(a))) + np.sum(np.abs(b - np.median(b))))
    k = int(np.argmin(cost))
    return float(grid[k]), float(np.median(rr[dd >= grid[k]]) - np.median(rr[dd < grid[k]]))


def in_cone(Pr_sel):
    c = unitvec(np.array([CONE[0]]), np.array([CONE[1]]))[0]
    return np.degrees(np.arccos(np.clip(Pr_sel @ c, -1, 1)))


def main():
    S = load_sources()
    P, nrm, sdec = dense_spine()
    dp, ds, by, _ = coords(S, P, nrm, sdec)
    out = {'n_sources': int(S['ra'].size)}
    # consistency with wave-1 coordinates
    out['check_vs_wave1'] = dict(median_abs_dperp_diff=float(np.median(np.abs(dp - S['dperp_w1']))),
                                 frac_sign_agree=float(np.mean(np.sign(dp) == np.sign(S['dperp_w1']))))
    halves = {'north': (20.0, 28.2), 'south': (12.0, 20.0)}
    true = {}
    for h, (lo, hi) in halves.items():
        r = {}
        for key in ('rm', 'resid'):
            s0, e, w = step_stat(S[key], dp, ds, by, lo, hi)
            bs = []
            ie, iw = np.where(e)[0], np.where(w)[0]
            for _ in range(4000):
                bs.append(np.median(S[key][RNG.choice(ie, ie.size)]) - np.median(S[key][RNG.choice(iw, iw.size)]))
            r[key] = dict(S=s0, err_boot=float(np.std(bs)), n_E=int(e.sum()), n_W=int(w.sum()),
                          frac_neg_E=float(np.mean(S['rm'][e] < 0)), frac_neg_W=float(np.mean(S['rm'][w] < 0)))
        m = (~by) & (ds >= lo) & (ds < hi) & (np.abs(dp) < 6)
        r['l1_step_location_deg'], r['l1_step_size'] = l1_step(dp[m], S['rm'][m])
        true[h] = r
        print('true', h, r)
    out['true'] = true
    sel_true = {h: step_stat(S['rm'], dp, ds, by, *halves[h])[1] | step_stat(S['rm'], dp, ds, by, *halves[h])[2]
                for h in halves}

    # ---------------------------------------------------------------- placement nulls
    allp = P
    cen = allp.mean(0)
    cen /= np.linalg.norm(cen)
    cone_c = unitvec(np.array([CONE[0]]), np.array([CONE[1]]))[0]
    nulls = {}
    for mode in ('translate', 'spin_translate'):
        res = {h: dict(rm=[], resid=[], l1=[], overlap=[]) for h in halves}
        tries = 0
        while min(len(res[h]['rm']) for h in halves) < 600 and tries < 6000:
            tries += 1
            if mode == 'translate':
                Rs = np.eye(3)
                ang = RNG.uniform(1.0, 10.0)
            else:
                Rs = rotmat(cone_c, RNG.uniform(0, 360))
                ang = RNG.uniform(0.0, 8.0)
            a = np.cross(cen, RNG.normal(size=3))
            a /= np.linalg.norm(a)
            R = rotmat(a, ang) @ Rs
            dpn, dsn, byn, Pr = coords(S, P, nrm, sdec, R)
            for h, (lo, hi) in halves.items():
                seg = (sdec >= lo) & (sdec < hi)
                # side bands (+-6 deg) of the placed half must lie inside the cone
                if in_cone(Pr[seg]).max() > CONE[2] - 6.5:
                    continue
                s1, e, w = step_stat(S['rm'], dpn, dsn, byn, lo, hi)
                if s1 is None:
                    continue
                s2, _, _ = step_stat(S['resid'], dpn, dsn, byn, lo, hi)
                m = (~byn) & (dsn >= lo) & (dsn < hi) & (np.abs(dpn) < 6)
                loc, size = l1_step(dpn[m], S['rm'][m])
                ov = float(np.mean(sel_true[h][e | w]))
                res[h]['rm'].append(s1)
                res[h]['resid'].append(s2)
                res[h]['l1'].append((loc, size))
                res[h]['overlap'].append(ov)
        summ = {}
        for h in halves:
            x = np.array(res[h]['rm'])
            xr = np.array(res[h]['resid'])
            ov = np.array(res[h]['overlap'])
            l1 = np.array(res[h]['l1'])
            t = true[h]['rm']['S']
            tr_ = true[h]['resid']['S']
            lo_ov = ov < 0.3
            summ[h] = dict(
                n=int(x.size), n_nonoverlap=int(lo_ov.sum()), tries=tries,
                null_sd=float(np.std(x)), null_mean=float(np.mean(x)),
                p_abs_ge_true=float(np.mean(np.abs(x) >= abs(t))),
                p_ge_true=float(np.mean(x >= t)),
                nonoverlap_null_sd=float(np.std(x[lo_ov])) if lo_ov.sum() > 10 else None,
                nonoverlap_p_abs_ge_true=float(np.mean(np.abs(x[lo_ov]) >= abs(t))) if lo_ov.sum() > 10 else None,
                resid_null_sd=float(np.std(xr)), resid_p_abs_ge_true=float(np.mean(np.abs(xr) >= abs(tr_))),
                z_true=float((t - np.mean(x)) / np.std(x)),
                l1_size_null_sd=float(np.std(l1[:, 1])),
                p_l1_abs_ge_true=float(np.mean(np.abs(l1[:, 1]) >= abs(true[h]['l1_step_size']))),
                pcts_abs=np.percentile(np.abs(x), [50, 84, 95, 99]).tolist())
            print(mode, h, summ[h])
        nulls[mode] = dict(summary=summ, samples={h: dict(rm=res[h]['rm'], overlap=res[h]['overlap'],
                                                         resid=res[h]['resid']) for h in halves})
    out['nulls'] = {k: v['summary'] for k, v in nulls.items()}

    # ---------------------------------------------------------------- smooth-gradient prediction
    X = unitvec(S['ra'], S['dec'])
    x = (S['ra'] - CONE[0]) * np.cos(np.radians(S['dec']))
    y = S['dec'] - CONE[1]
    north_band = (~by) & (ds >= 20) & (ds < 28.2) & (np.abs(dp) < 7)
    fit = ~north_band & (np.hypot(x, y) < 18)
    pred = {}
    for name, A in (('plane', np.c_[np.ones_like(x), x, y]),
                    ('quadratic', np.c_[np.ones_like(x), x, y, x * x, x * y, y * y])):
        use = fit.copy()
        for _ in range(5):
            cf = np.linalg.lstsq(A[use], S['rm'][use], rcond=None)[0]
            r = S['rm'] - A @ cf
            sd = 1.4826 * np.median(np.abs(r[use] - np.median(r[use])))
            use = fit & (np.abs(r) < 3 * sd)
        model = A @ cf
        s_m, e, w = step_stat(model, dp, ds, by, 20.0, 28.2)
        s_r, _, _ = step_stat(S['rm'] - model, dp, ds, by, 20.0, 28.2)
        pred[name] = dict(coef=cf.tolist(), predicted_step=s_m, residual_step=s_r, robust_sd=float(sd),
                          n_fit=int(use.sum()))
        print('gradient', name, pred[name])
    out['smooth_gradient'] = pred
    # cross-spine profiles (north / south) for the figure
    prof = {}
    for h, (lo, hi) in halves.items():
        m = (~by) & (ds >= lo) & (ds < hi)
        edges = np.arange(-8, 8.01, 1.0)
        med, err, nn = [], [], []
        for b1, b2 in zip(edges[:-1], edges[1:]):
            mm = m & (dp >= b1) & (dp < b2)
            if mm.sum() >= 3:
                ii = np.where(mm)[0]
                med.append(float(np.median(S['rm'][mm])))
                err.append(float(np.std([np.median(S['rm'][RNG.choice(ii, ii.size)]) for _ in range(1000)])))
            else:
                med.append(None)
                err.append(None)
            nn.append(int(mm.sum()))
        prof[h] = dict(edges=edges.tolist(), median_rm=med, err=err, n=nn)
    out['profiles'] = prof
    json.dump(out, open(f'{ROOT}/data/derived/w3_northsouth_rm.json', 'w'), indent=1)
    np.savez_compressed(f'{ROOT}/data/derived/w3_northsouth_rm_nulls.npz',
                        **{f'{m}_{h}_{k}': np.array(nulls[m]['samples'][h][k]) for m in nulls for h in halves
                           for k in ('rm', 'overlap', 'resid')})


if __name__ == '__main__':
    main()
