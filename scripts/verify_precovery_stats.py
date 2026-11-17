"""verify_precovery step 5: jackknives, epoch splits, Dec profile, joint g+u significance, band ratios.

Reads data/derived/verify_precovery_fit.json (per-run amplitudes with off-window and in-window nulls).
Per-run error = max(formal, std of the off-window null amplitudes) ('off' mode; validated by chi2/dof ~ 1.1),
falling back to max(formal, robust in-window null scatter) when a run has < 5 off-window positions.
Output: data/derived/verify_precovery_stats.json
"""
import sys, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from astropy.time import Time
from verify_precovery_lib import DER

F = json.loads((DER / 'verify_precovery_fit.json').read_text())
rng = np.random.default_rng(7)


def eff(r):
    eo = r.get('e_null_off'); en = r.get('e_null_in') or 0.0
    return max(r['eA_formal'], eo) if eo else max(r['eA_formal'], en)


def wmean(A, e):
    w = 1 / np.asarray(e) ** 2
    m = np.sum(w * A) / w.sum(); em = 1 / np.sqrt(w.sum())
    chi = np.sum(w * (np.asarray(A) - m) ** 2) / max(len(A) - 1, 1)
    return float(m), float(em * max(1, np.sqrt(chi))), float(chi)


def runs_of(case, sel):
    return [r for r in F['cases'][case]['runs'] if sel(r)]


CORE = lambda r: 12.0 <= r['dec_cross'] < 20.5
NORTH = lambda r: r['dec_cross'] >= 20.5


def season(mjd):
    t = Time(mjd, format='mjd').datetime
    y = t.year if t.month >= 8 else t.year - 1
    return f'{y}-{(y + 1) % 100:02d}'


def jackknife(case, sel=CORE):
    rr = runs_of(case, sel)
    A = np.array([r['A'] for r in rr]); e = np.array([eff(r) for r in rr])
    full = wmean(A, e)
    loo = []
    for i in range(len(rr)):
        k = np.ones(len(rr), bool); k[i] = False
        loo.append(dict(run=rr[i]['run'], A=wmean(A[k], e[k])[0]))
    v = np.array([x['A'] for x in loo]); n = len(v)
    jk_err = float(np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2)))
    seas = {}
    for s in sorted(set(season(r['mjd']) for r in rr)):
        k = np.array([season(r['mjd']) == s for r in rr])
        if k.sum() >= 1:
            m, em, chi = wmean(A[k], e[k]) if k.sum() > 1 else (float(A[k][0]), float(e[k][0]), np.nan)
            mo, emo, _ = wmean(A[~k], e[~k])
            seas[s] = dict(n=int(k.sum()), runs=[r['run'] for r, kk in zip(rr, k) if kk], A=m, e=em,
                           A_without=mo, e_without=emo)
    return dict(full=full, loo=loo, jk_err=jk_err, loo_min=float(v.min()), loo_max=float(v.max()),
                loo_min_run=loo[int(np.argmin(v))]['run'], loo_max_run=loo[int(np.argmax(v))]['run'], seasons=seas)


def null_sig(case, sel=CORE, ndraw=20000):
    rr = runs_of(case, sel)
    A = np.array([r['A'] for r in rr]); e = np.array([eff(r) for r in rr]); w = 1 / e ** 2
    m0 = float(np.sum(w * A) / w.sum())
    keys = sorted({k for r in rr for k in r['null_off']}, key=float)
    same = []
    for k in keys:
        have = [(i, r['null_off'][k]['A']) for i, r in enumerate(rr) if k in r['null_off']]
        if len(have) >= max(3, len(rr) // 2):
            ii = np.array([h[0] for h in have]); v = np.array([h[1] for h in have])
            same.append(float(np.sum(w[ii] * v) / np.sum(w[ii])))
    pools = []
    for r in rr:
        p = np.array([v['A'] for v in r['null_off'].values()])
        q = np.array(r['null_in'], float); q = q[np.isfinite(q)]
        pools.append(p if len(p) >= 3 else (q if len(q) else np.array([0.0])))
    draws = np.array([[p[rng.integers(len(p))] for p in pools] for _ in range(ndraw)])
    comb = (draws * w).sum(1) / w.sum()
    return dict(A=m0, same_shift=same, same_std=float(np.std(same)) if len(same) > 2 else None,
                sig_same=float(m0 / np.std(same)) if len(same) > 2 else None,
                rand_std=float(np.std(comb)), sig_rand=float(m0 / np.std(comb)),
                p_rand=float(np.mean(comb >= m0)), n_runs=len(rr))


def joint_gu(gcase, ucase, sel=CORE, ndraw=20000):
    """Joint significance of g ([OIII]) and u ([OII]) at the known ridge position against nulls drawn at the SAME
    sky position in both bands (keeps any common sky structure)."""
    G = {r['run']: r for r in runs_of(gcase, sel)}; U = {r['run']: r for r in runs_of(ucase, sel)}
    runs = sorted(set(G) & set(U))
    eg = np.array([eff(G[k]) for k in runs]); eu = np.array([eff(U[k]) for k in runs])
    wg = 1 / eg ** 2; wu = 1 / eu ** 2
    Ag = np.sum(wg * [G[k]['A'] for k in runs]) / wg.sum(); Au = np.sum(wu * [U[k]['A'] for k in runs]) / wu.sum()
    keys = [sorted(set(G[k]['null_off']) & set(U[k]['null_off'])) for k in runs]
    sg, su = [], []
    for _ in range(ndraw):
        a, b = [], []
        for k, ks in zip(runs, keys):
            if len(ks) >= 3:
                s = ks[rng.integers(len(ks))]
                a.append(G[k]['null_off'][s]['A']); b.append(U[k]['null_off'][s]['A'])
            else:
                qg = np.array(G[k]['null_in'], float); qu = np.array(U[k]['null_in'], float)
                ok = np.isfinite(qg) & np.isfinite(qu)
                j = rng.choice(np.where(ok)[0]) if ok.any() else None
                a.append(qg[j] if j is not None else 0.0); b.append(qu[j] if j is not None else 0.0)
        sg.append(np.sum(wg * a) / wg.sum()); su.append(np.sum(wu * b) / wu.sum())
    sg = np.array(sg); su = np.array(su)
    zg = Ag / sg.std(); zu = Au / su.std()
    Zn = sg / sg.std() + su / su.std()
    Z0 = zg + zu
    # same-shift (coherent sky position for all runs) joint statistic
    allk = sorted({s for ks in keys for s in ks}, key=float)
    ss = []
    for s in allk:
        ii = [i for i, ks in enumerate(keys) if s in ks]
        if len(ii) >= len(runs) // 2:
            a = np.array([G[runs[i]]['null_off'][s]['A'] for i in ii]); b = np.array([U[runs[i]]['null_off'][s]['A'] for i in ii])
            ss.append((np.sum(wg[ii] * a) / wg[ii].sum(), np.sum(wu[ii] * b) / wu[ii].sum()))
    ss = np.array(ss)
    out = dict(n_runs=len(runs), A_g=float(Ag), A_u=float(Au), z_g_rand=float(zg), z_u_rand=float(zu),
               corr_null_gu=float(np.corrcoef(sg, su)[0, 1]), joint_sig_rand=float(Z0 / Zn.std()),
               p_joint_rand=float(np.mean(Zn >= Z0)))
    if len(ss) > 3:
        zs = ss[:, 0] / ss[:, 0].std() + ss[:, 1] / ss[:, 1].std()
        out.update(same_shift_pairs=ss.tolist(), joint_sig_same=float((Ag / ss[:, 0].std() + Au / ss[:, 1].std()) / zs.std()),
                   corr_same_gu=float(np.corrcoef(ss[:, 0], ss[:, 1])[0, 1]))
    # band ratio u/g per run and combined (photon units of the cases)
    ratio = Au / Ag
    # bootstrap over runs for the ratio
    bs = []
    Agr = np.array([G[k]['A'] for k in runs]); Aur = np.array([U[k]['A'] for k in runs])
    for _ in range(4000):
        i = rng.integers(len(runs), size=len(runs))
        bs.append(np.sum(wu[i] * Aur[i]) / wu[i].sum() / (np.sum(wg[i] * Agr[i]) / wg[i].sum()))
    out.update(ratio_u_over_g=float(ratio), ratio_boot_16_50_84=[float(x) for x in np.percentile(bs, [16, 50, 84])])
    return out


def dec_profile(case):
    rr = F['cases'][case]['runs']
    bins = [(12, 14), (14, 15.5), (15.5, 17.0), (17.0, 18.0), (18.0, 19.5), (19.5, 21.0), (21.0, 23.5), (23.5, 27)]
    out = []
    for lo, hi in bins:
        k = [r for r in rr if lo <= r['dec_cross'] < hi]
        if not k:
            continue
        A = np.array([r['A'] for r in k]); e = np.array([eff(r) for r in k])
        m, em, chi = wmean(A, e) if len(k) > 1 else (float(A[0]), float(e[0]), np.nan)
        out.append(dict(dec_lo=lo, dec_hi=hi, n=len(k), A=m, e=em, chi2=chi, runs=[r['run'] for r in k]))
    return out


def per_run_table(case):
    rows = []
    for r in sorted(F['cases'][case]['runs'], key=lambda r: r['mjd']):
        rows.append(dict(run=r['run'], date=Time(r['mjd'], format='mjd').iso[:10], dec_cross=round(r['dec_cross'], 2),
                         A=r['A'], e=eff(r), snr=r['A'] / eff(r), e_formal=r['eA_formal'], e_off=r.get('e_null_off'),
                         e_in=r.get('e_null_in')))
    return rows


def main():
    res = {'jackknife': {}, 'nulls': {}, 'joint': {}, 'dec_profile': {}, 'per_run': {}}
    for c in ['g_psp_gauss028_OIIIunits', 'g_framessub_gauss028_OIIIunits', 'g_frames_gauss028_OIIIunits',
              'u_psp_gauss028_OIIunits', 'u_framessub_gauss028_OIIunits', 'u_frames_gauss028_OIIunits',
              'g_psp_desi', 'g_framessub_desi', 'g_psp_fibre', 'g_framessub_fibre', 'r_psp_gauss028', 'i_psp_gauss028',
              'z_psp_gauss028', 'r_framessub_gauss028', 'i_framessub_gauss028']:
        res['jackknife'][c] = jackknife(c)
        res['nulls'][c] = null_sig(c)
        res['nulls'][c + '_north'] = null_sig(c, NORTH)
        res['dec_profile'][c] = dec_profile(c)
        res['per_run'][c] = per_run_table(c)
        j = res['jackknife'][c]; n = res['nulls'][c]; nn = res['nulls'][c + '_north']
        print(f"{c:34s} A={j['full'][0]:+.2f}+-{j['full'][1]:.2f} chi2 {j['full'][2]:.2f} jk {j['jk_err']:.2f} "
              f"LOO [{j['loo_min']:+.2f} (-{j['loo_min_run']}), {j['loo_max']:+.2f} (-{j['loo_max_run']})] "
              f"sig same {n['sig_same'] if n['sig_same'] is None else round(n['sig_same'], 1)} rand {n['sig_rand']:.1f} p {n['p_rand']:.4f}"
              f" | north A={nn['A']:+.2f} sig_rand {nn['sig_rand']:.1f} same {nn['sig_same'] if nn['sig_same'] is None else round(nn['sig_same'], 1)}")
        for s, v in j['seasons'].items():
            print(f"      season {s}: n={v['n']} A={v['A']:+.2f}+-{v['e']:.2f}   without: {v['A_without']:+.2f}+-{v['e_without']:.2f}")
    for est in ('psp', 'framessub', 'frames'):
        for sel, nm in ((CORE, 'core'), (NORTH, 'north')):
            k = f'{est}_{nm}'
            res['joint'][k] = joint_gu(f'g_{est}_gauss028_OIIIunits', f'u_{est}_gauss028_OIIunits', sel)
            j = res['joint'][k]
            print('joint', k, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in j.items() if kk != 'same_shift_pairs'})
    (DER / 'verify_precovery_stats.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
