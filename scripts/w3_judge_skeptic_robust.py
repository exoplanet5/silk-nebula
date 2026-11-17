"""w3_judge_skeptic: robustness of the cross-profile results of w3_judge_skeptic_tests.py, plus a flat-sheet
model variant.

Observed variants (Dec 12-15 single strand, and Dec 12-17.6):
  NSNS [O III] asymmetry: median/mean, 0.025/0.05 deg bins, fit range 0.6/0.9/1.2 deg; non-parametric
  ratio R_asym = <I>(-0.45..-0.15) / <I>(+0.15..+0.45) after a linear baseline fitted at 0.6 < |x| < 1.1.
  GALEX FUV/NUV shift vs the NSNS template: per-tile reference |x| > 0.3 or > 0.5; fit range 0.6/0.9/1.2.
  DESI [O II] vs [O III] shift: fit range 0.8/1.2/1.5, bins 0.1/0.15; and a first-moment difference.
Model variant: the 25 consistent shock models seen as a FLAT sheet edge-on (profile = emissivity j(s) across the
layer, no curvature tail), convolved to the observed FWHM; R_asym computed for both curved and flat geometries.
Output: data/derived/w3_judge_skeptic_robust.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_skeptic_robust.py
"""
import sys
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import DER
import w3_judge_skeptic_tests as T
import w3_judge_skeptic_modelprof as MP


def r_asym(x, y):
    m = np.isfinite(y)
    b = m & (np.abs(x) > 0.6) & (np.abs(x) < 1.1)
    if b.sum() < 4:
        return np.nan
    c = np.polyfit(x[b], y[b], 1)
    yy = y - np.polyval(c, x)
    w = m & (x >= -0.45) & (x <= -0.15); e = m & (x >= 0.15) & (x <= 0.45)
    return float(np.mean(yy[w]) / np.mean(yy[e])) if np.mean(yy[e]) > 0 else np.nan


def jk(vals):
    v = np.array([q for q in vals if np.isfinite(q)]); K = len(v)
    return float(np.sqrt((K - 1) / K * np.sum((v - v.mean()) ** 2))) if K > 2 else np.nan


def observed():
    import verify_xrayfuv_galex_measure as GM
    from verify_desi_common import good_fit
    from verify_xrayfuv_common import radec_to_vec
    ra_n, de_n, v_n = T.nsns_data()
    xn, sn = T.CC0.coords(radec_to_vec(ra_n, de_n), maxdist_deg=1.6)
    G = {b: GM.prep_band(b) for b in ('FUV', 'NUV')}
    Gs = {b: T.CC0.coords(radec_to_vec(G[b]['ra'], G[b]['dec']), maxdist_deg=1.6)[1] for b in G}
    lf = np.load(DER / 'verify_desi_linefluxes.npz'); cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    det = json.loads((DER / 'verify_desi_detect.json').read_text())
    okd = (good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & good_fit(lf, 'o2sum', scale=det['err_scale_from_nulls'])
           & cm['ok_contrib'] & np.isfinite(cm['kR_5008']) & np.isfinite(cm['kR_3727']))
    d3 = np.where(okd, lf['o3__amp'][:, 0] * cm['kR_5008'], np.nan)
    d2 = np.where(okd, lf['o2sum__amp'][:, 0] * cm['kR_3727'], np.nan)
    ds = T.CC0.coords(radec_to_vec(cm['ra'], cm['dec']), maxdist_deg=1.6)[1]
    out = {}
    for seg, (smin, smax) in (('dec12_15', (0.0, T.S15)), ('dec12_17p6', (0.0, T.S176))):
        chunks = np.linspace(smin, smax, 11)
        res = {'nsns': [], 'galex': [], 'desi': []}

        def drop(vals, s, k):
            return vals if k is None else np.where((s >= chunks[k]) & (s < chunks[k + 1]), np.nan, vals)
        # ---- NSNS asymmetry variants
        for stat in ('median', 'mean'):
            for bw in (0.025, 0.05):
                e = np.arange(-1.2, 1.2001, bw); xc = 0.5 * (e[1:] + e[:-1])
                profs = [T.prof_generic(ra_n, de_n, drop(v_n, sn, k), e, smin, smax, stat=stat) for k in [None] + list(range(10))]
                ra_full = r_asym(xc, profs[0]); ra_jk = jk([r_asym(xc, p) for p in profs[1:]])
                for xmax in (0.6, 0.9, 1.2):
                    f = T.fit_split(xc, profs[0], xmax=xmax)
                    fj = [T.fit_split(xc, p, xmax=xmax) for p in profs[1:]]
                    res['nsns'].append(dict(stat=stat, bin=bw, xmax=xmax, asym=f['asym_w_over_e'] if f else np.nan,
                                            asym_jk=jk([q['asym_w_over_e'] for q in fj if q]),
                                            x0=f['x0'] if f else np.nan, R_asym=ra_full, R_asym_jk=ra_jk))
                if stat == 'median' and bw == 0.025:
                    tmpl_prof = profs
        # NSNS template (median, 0.025, xmax 0.9) for the GALEX shifts
        e = np.arange(-1.2, 1.2001, 0.025); xc = 0.5 * (e[1:] + e[:-1])
        fN = T.fit_split(xc, tmpl_prof[0], xmax=0.9)
        tx = np.linspace(-1.2, 1.2, 961)
        tmpl = T.splitnorm([fN['A'], fN['x0'], fN['sig_w'], fN['sig_e'], 0, 0], tx)
        for band in ('FUV', 'NUV'):
            D = G[band]
            for dout in (0.3, 0.5):
                profs = [T.prof_generic(D['ra'], D['dec'], drop(D['v'].astype(float), Gs[band], k), e, smin, smax,
                                        stat='mean', tile=D['tile'], demean_out=dout) for k in [None] + list(range(10))]
                for xmax in (0.6, 0.9, 1.2):
                    d0 = T.shift_fit(xc, profs[0], tx, tmpl, xmax=xmax)[0]
                    dj = jk([T.shift_fit(xc, p, tx, tmpl, xmax=xmax)[0] for p in profs[1:]])
                    f = T.fit_split(xc, profs[0], xmax=xmax)
                    res['galex'].append(dict(band=band, demean_out=dout, xmax=xmax, shift=d0, shift_jk=dj,
                                             asym=f['asym_w_over_e'] if f else np.nan, R_asym=r_asym(xc, profs[0])))
        # ---- DESI [O II] vs [O III]
        for bw in (0.1, 0.15):
            e2 = np.arange(-1.5, 1.5001, bw); x2 = 0.5 * (e2[1:] + e2[:-1])
            p3 = [T.prof_generic(cm['ra'], cm['dec'], drop(d3, ds, k), e2, smin, smax) for k in [None] + list(range(10))]
            p2 = [T.prof_generic(cm['ra'], cm['dec'], drop(d2, ds, k), e2, smin, smax) for k in [None] + list(range(10))]

            def one(a3, a2, xmax):
                f3 = T.fit_split(x2, a3, xmax=xmax)
                if f3 is None:
                    return np.nan, np.nan
                t3 = T.splitnorm([f3['A'], f3['x0'], f3['sig_w'], f3['sig_e'], 0, 0], tx)
                dl = T.shift_fit(x2, a2, tx, t3, xmax=xmax)[0]
                f2 = T.fit_split(x2, a2, xmax=xmax)
                mom = (f2['moment'] - f3['moment']) if (f2 and np.isfinite(f2['moment']) and np.isfinite(f3['moment'])) else np.nan
                return dl, mom
            for xmax in (0.8, 1.2, 1.5):
                d0, m0 = one(p3[0], p2[0], xmax)
                rr = [one(a, b, xmax) for a, b in zip(p3[1:], p2[1:])]
                res['desi'].append(dict(bin=bw, xmax=xmax, shift=d0, shift_jk=jk([r[0] for r in rr]),
                                        dmoment=m0, dmoment_jk=jk([r[1] for r in rr]),
                                        R_asym_O3=r_asym(x2, p3[0]), R_asym_O2=r_asym(x2, p2[0])))
        out[seg] = res
        print(seg, json.dumps(res, default=float)[:3000], flush=True)
    return out


def model_variants():
    summ = json.loads((DER / 'verify_shockmodel_summary.json').read_text())
    TR = np.load(DER / 'verify_shockmodel_tracks.npz')
    xg = np.arange(-1.5, 1.5001, 0.005)
    rows = []
    for m in summ['consistent']:
        tag = f"v{m['v']}_x{m['xn']}_b{m['b']}_asplund09_trn"
        N = TR[tag + '__N']; s_ref = TR[tag + '__s']
        istar = int(np.argmin(np.abs(np.log10(np.maximum(N, 1)) - m['logN'])))
        rec = dict(tag=tag)
        for geom in ('curved', 'flat'):
            sky = {}
            for key in ('F_O3', 'F_O2', 'F_CIV'):
                if geom == 'curved':
                    hs, I, _ = MP.proj(s_ref, TR[tag + '__' + key], istar, m['n0'])
                else:
                    s = s_ref[:istar + 1] * MP.N0REF / m['n0']; Fc = TR[tag + '__' + key][:istar + 1]
                    hs = 0.5 * (s[1:] + s[:-1]); I = np.diff(Fc) / np.maximum(np.diff(s), 1e-30)
                hd = np.degrees(hs / (MP.D_PC * MP.PC))
                sky[key] = (hd, I)
            hpk = sky['F_O3'][0][np.argmax(sky['F_O3'][1])]
            arr = {k: np.interp(hpk - xg, hd, I, left=0.0, right=0.0) for k, (hd, I) in sky.items()}
            f0 = T.fit_split(xg, arr['F_O3'], xmax=1.2)
            for seg, fw in MP.FWHM_OBS_STACK.items():
                sig_k = np.sqrt(max(fw ** 2 - f0['fwhm'] ** 2, 1e-6)) / 2.3548
                conv = {k: gaussian_filter1d(v, sig_k / 0.005, mode='constant') for k, v in arr.items()}
                fO3 = T.fit_split(xg, conv['F_O3'], xmax=0.9)
                tmpl = T.splitnorm([fO3['A'], fO3['x0'], fO3['sig_w'], fO3['sig_e'], 0, 0], xg)
                rec[f'{geom}_{seg}'] = dict(asym=fO3['asym_w_over_e'], R_asym=r_asym(xg, conv['F_O3']),
                                            shift_CIV=T.shift_fit(xg, conv['F_CIV'] / conv['F_CIV'].max(), xg, tmpl / tmpl.max())[0],
                                            shift_O2=T.shift_fit(xg, conv['F_O2'] / conv['F_O2'].max(), xg, tmpl / tmpl.max())[0],
                                            intrinsic_fwhm=f0['fwhm'])
        rows.append(rec)
    summary = {}
    for key in rows[0]:
        if key == 'tag':
            continue
        summary[key] = {q: [float(np.min([r[key][q] for r in rows])), float(np.median([r[key][q] for r in rows])),
                            float(np.max([r[key][q] for r in rows]))] for q in rows[0][key]}
    print(json.dumps(summary, indent=1))
    return dict(rows=rows, summary_min_median_max=summary)


def main():
    mod = model_variants()
    obs = observed()
    (DER / 'w3_judge_skeptic_robust.json').write_text(json.dumps(dict(observed=obs, model=mod), indent=1, default=float))
    print('saved')


if __name__ == '__main__':
    main()
