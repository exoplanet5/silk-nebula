"""verify_precovery step 7: cross-calibration, band SED, null calibration, epoch and confounder tests.

All inputs are my own wave-2 products (verify_precovery_fields.npz, verify_precovery_fit.json) plus read-only wave-1
tables used as independent comparison data (spectra_sdss_lamost_sdss_fibres.csv, spectra_desi_results.json).

(1) Template cross-calibration: noise-free forward model of the DESI (forward-model amplitudes, trace widths) and
    SDSS sky-fibre (strand amplitudes, FWHM 0.28 deg) ridge models through exactly the same Field-level fit as the
    data (same runs, fields, nuisance model, gauss028 unit template, footprint averaging). Gives the gauss028 amplitude
    the imaging WOULD measure if the sky contained the DESI or the fibre model -> direct comparison with A_obs.
    Also a template-mismatch injection: the DESI model injected at the off-window (shifted) positions and recovered
    with the shifted gauss028 template (bias of the fixed-width template when the true widths / amplitudes vary).
(2) Field-level band SED of the core ridge amplitude (u, g, r, i, z in nmgy arcsec^-2) vs an emission-line spectrum
    (DESI ratios) and a flat-f_nu continuum.
(3) Null calibration: distribution of per-run off-window null z = A_null / e_eff (all bands); tail fractions.
(4) Temporal test at fixed position: per-run scale relative to the DESI-model prediction vs date (core runs),
    and the Dec 17.2-17.7 run group (2003-12 ... 2007-02).
(5) Confounder correlations: per-run g amplitude vs Moon altitude / illumination, ecliptic latitude, sky level.
(6) SDSS sky fibres, gauss028 template fitted to Dec 12-20.5 only (same Dec range as the imaging core).
Output: data/derived/verify_precovery_extra.json
"""
import sys, json, csv
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from scipy.stats import spearmanr
from astropy.time import Time
import verify_precovery_fit as VF
from verify_precovery_lib import DER, load_strands, make_model, band_response_per_R5007

Z = VF.Z
F = json.loads((DER / 'verify_precovery_fit.json').read_text())
CORE = lambda r: 12.0 <= r['dec_cross'] < 20.5
NORTH = lambda r: r['dec_cross'] >= 20.5


def eff(r):
    eo = r.get('e_null_off'); en = r.get('e_null_in') or 0.0
    return max(r['eA_formal'], eo) if eo else max(r['eA_formal'], en)


def wmean(A, e):
    A = np.asarray(A, float); w = 1 / np.asarray(e, float) ** 2
    m = np.sum(w * A) / w.sum(); em = 1 / np.sqrt(w.sum())
    chi = np.sum(w * (A - m) ** 2) / max(len(A) - 1, 1)
    return float(m), float(em * max(1.0, np.sqrt(chi))), float(chi)


# ------------------------------------------------------------------------------------------------ (1)
def cross_calibration():
    out = {}
    runs = [r for r in F['cases']['g_psp_gauss028_OIIIunits']['runs'] if CORE(r)]
    w = {r['run']: 1 / eff(r) ** 2 for r in runs}
    pred = {'desi': [], 'fibre': []}
    for r in runs:
        idx = VF.select(r['run'])
        T0 = Z['T_gauss028_mean'][idx]
        for k in pred:
            yy = np.zeros(len(Z['ra'])); yy[idx] = Z[f'T_{k}_mean'][idx]          # model sky in R ([OIII]5007)
            fr = VF.fit(idx, yy, [T0], clip=99)
            pred[k].append((r['run'], r['dec_cross'], float(fr['coef'][0]), r['A'], eff(r)))
    for k, v in pred.items():
        v = np.array(v)
        ww = 1 / v[:, 4] ** 2
        Ap = float(np.sum(ww * v[:, 2]) / ww.sum())
        out[f'A_pred_{k}'] = Ap
        out[f'per_run_pred_{k}'] = [dict(run=int(a), dec=float(b), A_pred=float(c), A_obs=float(d), e=float(e))
                                    for a, b, c, d, e in v]
    A_obs, e_obs, chi = wmean([r['A'] for r in runs], [eff(r) for r in runs])
    out['A_obs_g_psp_OIIIunits'] = A_obs; out['e_obs'] = e_obs
    out['ratio_obs_over_desi'] = A_obs / out['A_pred_desi']; out['ratio_obs_over_desi_err'] = e_obs / out['A_pred_desi']
    out['ratio_obs_over_fibre'] = A_obs / out['A_pred_fibre']; out['ratio_obs_over_fibre_err'] = e_obs / out['A_pred_fibre']
    # template-mismatch injection at off-window positions (DESI model injected, gauss028 recovered)
    SH = Z['shifts']; rec = []
    for r in runs:
        for s in VF.OFF_SHIFTS:
            i = int(np.where(np.isclose(SH, s))[0][0])
            jdx = VF.select(r['run'], i)
            if jdx is None:
                continue
            yy = np.zeros(len(Z['ra'])); yy[jdx] = Z['TS_desi'][i][jdx]
            Tg = Z['TS_gauss028'][i][jdx]
            fr = VF.fit(jdx, yy, [Tg], clip=99)
            rec.append(float(fr['coef'][0]))
    rec = np.array(rec)
    out['inject_desi_recover_gauss028_offwindow'] = dict(n=len(rec), mean=float(rec.mean()), std=float(rec.std()),
                                                         p16_50_84=[float(x) for x in np.percentile(rec, [16, 50, 84])])
    return out


# ------------------------------------------------------------------------------------------------ (2)
def band_sed():
    res = {}
    conv = F['conv_nmgy_arcsec2_per_R5007']
    for est in ('psp', 'framessub'):
        for sel, nm in ((CORE, 'core'), (NORTH, 'north')):
            d = {}
            for b in 'ugriz':
                case = f'{b}_{est}_gauss028' if not (b == 'g' and est == 'psp') else 'g_psp_gauss028'
                rr = [r for r in F['cases'][case]['runs'] if sel(r)]
                m, e, chi = wmean([r['A'] for r in rr], [eff(r) for r in rr])
                d[b] = dict(nmgy=m * conv[b] * 1e3, err=e * conv[b] * 1e3, chi2=chi, n=len(rr))   # 1e-3 nmgy arcsec^-2
            g = d['g']['nmgy']
            line = {b: band_response_per_R5007(b)[0] / band_response_per_R5007('g')[0] for b in 'ugriz'}
            # best-fit scale and chi2 for each SED hypothesis (u, g, r, i; z excluded: jackknife-unstable)
            hyp = {}
            for hn, shape in (('line_DESI', line), ('flat_fnu', {b: 1.0 for b in 'ugriz'})):
                x = np.array([shape[b] for b in 'ugri']); y = np.array([d[b]['nmgy'] for b in 'ugri'])
                s = np.array([d[b]['err'] for b in 'ugri'])
                k = np.sum(x * y / s ** 2) / np.sum(x ** 2 / s ** 2)
                hyp[hn] = dict(scale=float(k), chi2=float(np.sum((y - k * x) ** 2 / s ** 2)), dof=3)
            # g amplitude fixed: r and i prediction vs observation
            ri = {}
            for hn, shape in (('line_DESI', line), ('flat_fnu', {b: 1.0 for b in 'ugriz'})):
                ri[hn] = {b: dict(pred=g * shape[b], obs=d[b]['nmgy'], err=d[b]['err'],
                                  z=(d[b]['nmgy'] - g * shape[b]) / np.hypot(d[b]['err'], d['g']['err'] * shape[b]))
                          for b in 'uriz'}
            res[f'{est}_{nm}'] = dict(bands=d, ratio_to_g={b: (d[b]['nmgy'] / g if g else None) for b in 'ugriz'},
                                      expected_line_ratio=line, hypotheses=hyp, given_g=ri)
    return res


# ------------------------------------------------------------------------------------------------ (3)
def null_calibration():
    out = {}
    for case in ('g_psp_gauss028', 'g_framessub_gauss028', 'u_psp_gauss028', 'r_psp_gauss028', 'i_psp_gauss028'):
        z = []
        for r in F['cases'][case]['runs']:
            e = eff(r)
            z += [v['A'] / e for v in r['null_off'].values()]
        z = np.array(z)
        out[case] = dict(n=len(z), std=float(np.std(z)), mean=float(np.mean(z)),
                         frac_ge3=float(np.mean(z >= 3)), frac_le_m3=float(np.mean(z <= -3)),
                         frac_abs_ge3=float(np.mean(np.abs(z) >= 3)))
    return out


# ------------------------------------------------------------------------------------------------ (4)
def temporal():
    out = {}
    rr = [r for r in F['cases']['g_psp_desi']['runs'] if CORE(r)]
    t = np.array([Time(r['mjd'], format='mjd').decimalyear for r in rr])
    A = np.array([r['A'] for r in rr]); e = np.array([eff(r) for r in rr])
    w = 1 / e ** 2
    X = np.c_[np.ones(len(t)), t - 2004.0]
    cov = np.linalg.inv(X.T @ (X * w[:, None])); c = cov @ (X.T @ (w * A))
    chi = np.sum(w * (A - X @ c) ** 2) / (len(A) - 2)
    out['desi_scale_vs_year'] = dict(intercept_2004=float(c[0]), slope_per_yr=float(c[1]),
                                     slope_err=float(np.sqrt(cov[1, 1] * max(chi, 1))), chi2_dof=float(chi), n=len(A))
    for case in ('g_psp_gauss028_OIIIunits', 'g_framessub_gauss028_OIIIunits', 'u_psp_gauss028_OIIunits'):
        g = [r for r in F['cases'][case]['runs'] if 17.2 <= r['dec_cross'] <= 17.7]
        early = [r for r in g if r['mjd'] < 54000]; late = [r for r in g if r['mjd'] >= 54000]
        me, ee, _ = wmean([r['A'] for r in early], [eff(r) for r in early])
        ml, el, _ = wmean([r['A'] for r in late], [eff(r) for r in late])
        out[f'dec17p5_{case}'] = dict(early_runs=[r['run'] for r in early], early=(me, ee),
                                      late_runs=[r['run'] for r in late], late=(ml, el), diff=ml - me,
                                      diff_err=float(np.hypot(ee, el)))
    return out


# ------------------------------------------------------------------------------------------------ (5)
def confounders():
    out = {}
    for case in ('g_psp_gauss028_OIIIunits', 'g_framessub_gauss028_OIIIunits'):
        rr = [r for r in F['cases'][case]['runs'] if CORE(r)]
        snr = np.array([r['A'] / eff(r) for r in rr])
        d = {}
        for k in ('moon_alt', 'moon_illum', 'moon_sep', 'ecl_beta', 'sky'):
            v = np.array([r[k] for r in rr])
            rho, p = spearmanr(v, snr)
            d[k] = dict(spearman=float(rho), p=float(p))
        d['runs'] = [dict(run=r['run'], moon_alt=round(r['moon_alt'], 1), moon_illum=round(r['moon_illum'], 2),
                          ecl_beta=round(r['ecl_beta'], 1), sky=round(r['sky'], 3), snr=round(r['A'] / eff(r), 2)) for r in rr]
        out[case] = d
    return out


# ------------------------------------------------------------------------------------------------ (6)
def fibres_core():
    strands, _ = load_strands()
    r = list(csv.DictReader(open(DER / 'spectra_sdss_lamost_sdss_fibres.csv')))
    ra = np.array([float(x['ra']) for x in r]); dec = np.array([float(x['dec']) for x in r])
    I = np.array([float(x['I_OIII_cal']) for x in r]); e = np.array([float(x['I_OIII_cal_err']) for x in r])
    grp = np.array([f"{x['plate']}-{x['mjd']}-{x['sp']}" for x in r])
    T = make_model('gauss028', strands)(ra, dec)
    good = np.isfinite(I) & np.isfinite(e) & (e > 0)
    out = {}
    for nm, sel in (('dec12_20p5', (dec >= 12) & (dec < 20.5)), ('dec20p5_28', dec >= 20.5), ('all', np.ones(len(dec), bool))):
        idx = np.where(good & sel)[0]
        g, gi = np.unique(grp[idx], return_inverse=True)
        X = np.zeros((len(idx), len(g))); X[np.arange(len(idx)), gi] = 1; X = np.hstack([X, T[idx, None]])
        y = I[idx]; w = 1 / e[idx]; ok = np.ones(len(idx), bool)
        for _ in range(6):
            c, *_ = np.linalg.lstsq(X[ok] * w[ok, None], y[ok] * w[ok], rcond=None)
            res = (y - X @ c) * w; new = np.abs(res) < 5
            if (new == ok).all():
                break
            ok = new
        chi2 = np.sum(res[ok] ** 2) / max(ok.sum() - X.shape[1], 1)
        cov = np.linalg.pinv((X[ok] * w[ok, None]).T @ (X[ok] * w[ok, None])) * max(chi2, 1)
        # bootstrap over plate groups for a systematic-aware error
        rng = np.random.default_rng(3); bs = []
        ug = np.unique(gi)
        for _ in range(300):
            pick = rng.choice(ug, len(ug))
            rows = np.concatenate([np.where((gi == p) & ok)[0] for p in pick])
            Xb = X[rows]; keep = np.abs(Xb).sum(0) > 0
            cb, *_ = np.linalg.lstsq(Xb[:, keep] * w[rows, None], y[rows] * w[rows], rcond=None)
            bs.append(cb[-1])
        out[nm] = dict(A_R=float(c[-1]), e_formal=float(np.sqrt(cov[-1, -1])), e_boot=float(np.std(bs)),
                       chi2nu=float(chi2), n_fib=int(ok.sum()), n_groups=len(g), n_on=int((T[idx] > 0.5).sum()))
    return out


# ------------------------------------------------------------------------------------------------ (7)
def imaging_epochs():
    """Per-night combination and date-ordered cumulative significance of the CORE (Dec 12-20.5) g crossings,
    with the simultaneous u, r, i amplitudes of each night (line vs broadband diagnostic)."""
    out = {}
    cases = {'g': 'g_psp_gauss028_OIIIunits', 'gF': 'g_framessub_gauss028_OIIIunits', 'u': 'u_psp_gauss028_OIIunits',
             'r': 'r_psp_gauss028', 'i': 'i_psp_gauss028'}
    R = {k: {r['run']: r for r in F['cases'][c]['runs']} for k, c in cases.items()}
    nights = {}
    for run, r in R['g'].items():
        nights.setdefault(Time(r['mjd'], format='mjd').iso[:10], []).append(run)
    rows = []
    for d in sorted(nights):
        rr = nights[d]
        row = dict(date=d, runs=rr, dec=[round(R['g'][k]['dec_cross'], 2) for k in rr],
                   core=bool(all(12 <= R['g'][k]['dec_cross'] < 20.5 for k in rr)))
        for b in cases:
            ks = [k for k in rr if k in R[b]]
            m, e, _ = wmean([R[b][k]['A'] for k in ks], [eff(R[b][k]) for k in ks]) if len(ks) > 1 else \
                (R[b][ks[0]]['A'], eff(R[b][ks[0]]), np.nan)
            row[b] = (float(m), float(e), float(m / e))
        rows.append(row)
    # cumulative core g significance in date order
    A, E = [], []
    for row in rows:
        if row['core']:
            A.append(row['g'][0]); E.append(row['g'][1])
            w = 1 / np.array(E) ** 2
            row['cum_core_g_snr'] = float(np.sum(w * np.array(A)) / w.sum() * np.sqrt(w.sum()))
    out['nights'] = rows
    return out


def main():
    res = dict(cross_calibration=cross_calibration(), band_sed=band_sed(), null_calibration=null_calibration(),
               temporal=temporal(), confounders=confounders(), fibres_core=fibres_core(),
               imaging_epochs=imaging_epochs())
    for row in res['imaging_epochs']['nights']:
        print('night', row['date'], row['runs'], row['dec'], 'core' if row['core'] else 'north',
              ' '.join(f"{b}:{row[b][0]:+.2f}+-{row[b][1]:.2f}({row[b][2]:+.1f})" for b in ('g', 'gF', 'u', 'r', 'i')),
              'cum', round(row.get('cum_core_g_snr', np.nan), 2))
    cc = res['cross_calibration']
    print('cross-cal: A_obs %.2f+-%.2f  pred DESI %.2f  pred fibre %.2f  obs/DESI %.2f+-%.2f  obs/fibre %.2f+-%.2f' % (
        cc['A_obs_g_psp_OIIIunits'], cc['e_obs'], cc['A_pred_desi'], cc['A_pred_fibre'], cc['ratio_obs_over_desi'],
        cc['ratio_obs_over_desi_err'], cc['ratio_obs_over_fibre'], cc['ratio_obs_over_fibre_err']))
    print('inject DESI -> gauss028 (off windows):', cc['inject_desi_recover_gauss028_offwindow'])
    for k, v in res['band_sed'].items():
        print('SED', k, {b: (round(x['nmgy'], 2), round(x['err'], 2)) for b, x in v['bands'].items()},
              {h: (round(x['chi2'], 1)) for h, x in v['hypotheses'].items()},
              {h: {b: round(y['z'], 1) for b, y in x.items()} for h, x in v['given_g'].items()})
    for k, v in res['null_calibration'].items():
        print('nullcal', k, v)
    print('temporal', json.dumps(res['temporal'], indent=None))
    for k, v in res['confounders'].items():
        print('conf', k, {kk: vv for kk, vv in v.items() if kk != 'runs'})
    print('fibres', res['fibres_core'])
    (DER / 'verify_precovery_extra.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
