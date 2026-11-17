"""verify_precovery step 4: dated precovery table across all datasets known to contain the Silk [OIII] filament.

Per-epoch detection significance, computed here with my own fits wherever the per-epoch data are on disk:
  SDSS imaging sky    : per-run amplitudes from verify_precovery_fit.json (g, PSP estimator, gauss028 template);
                        significance = A / e_eff (e_eff = max(formal, in-window shifted-template null scatter)).
  SDSS sky fibres     : spectra_sdss_lamost_sdss_fibres.csv (wave-1 per-fibre absolute [OIII], diffuse-corrected).
                        Per plate-MJD: I_f = b_(plate,mjd,sp) + A T(ra,dec) with T = unit E+W Gaussian ridges
                        (FWHM 0.28 deg); 5-sigma clipping; error = max(formal x sqrt(chi2_nu), null scatter of the
                        template shifted by +-0.8..2.0 deg, fitted jointly).
  LAMOST DR12 LRS     : spectra_sdss_lamost_lrs_resid_stars.csv (per-star [OIII]5007 residual, on |d|<0.2 / off 0.6-2.5
                        deg) joined to obsdate (data/raw/lamost_dr12/lrs_catalogue_strandbox.csv); ON-OFF per plan.
  DESI DR1            : spectra_desi_fiber_table.csv + spectra_desi_catalog.csv (mean_mjd); per calendar month of
                        mean_mjd, slope k of OIII_R on the wave-1 forward-model residual prediction (model_resid_R),
                        significance k / sigma_k  (the model was fitted to all epochs: mild circularity, noted).
  Literature/archival : NSNS DR0.2 HiPS t_min, SWONS/AstroBin dates (notes/discovery.md, discovery_timeline.json),
                        DSS non-detections (notes/imaging_optical.md).
Output: data/derived/verify_precovery_table.json, verify_precovery_table.csv
"""
import sys, json, csv
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from astropy.time import Time
from verify_precovery_lib import DER, RAW, load_strands, make_model


def iso(mjd):
    return Time(mjd, format='mjd').iso[:10]


def eff(r):
    """per-run error: max(formal, off-window null std) if >= 5 off windows, else max(formal, in-window null MAD)."""
    eo = r.get('e_null_off'); en = r.get('e_null_in') or 0.0
    return max(r['eA_formal'], eo) if eo else max(r['eA_formal'], en)


IMG_CASES = ('g_psp_gauss028_OIIIunits', 'g_framessub_gauss028_OIIIunits', 'g_frames_gauss028_OIIIunits',
             'u_psp_gauss028_OIIunits', 'u_framessub_gauss028_OIIunits', 'u_frames_gauss028_OIIunits')


def sdss_imaging():
    F = json.loads((DER / 'verify_precovery_fit.json').read_text())
    out = {}
    for case in IMG_CASES:
        runs = F['cases'][case]['runs']
        rows = []
        for r in sorted(runs, key=lambda r: r['mjd']):
            e = eff(r)
            rows.append(dict(run=r['run'], date=iso(r['mjd']), mjd=r['mjd'], dec_cross=r['dec_cross'], A=r['A'],
                             e=e, snr=r['A'] / e))
        # cumulative significance in date order (inverse-variance)
        A = np.array([x['A'] for x in rows]); e = np.array([x['e'] for x in rows])
        core = np.ones(len(rows), bool)       # cumulative over all crossings (south and north)
        cum = []
        for i in range(len(rows)):
            m = core[:i + 1]
            if m.sum() == 0:
                cum.append(None); continue
            w = 1 / e[:i + 1][m] ** 2
            mm = np.sum(w * A[:i + 1][m]) / w.sum(); em = 1 / np.sqrt(w.sum())
            cum.append(float(mm / em))
        for x, c in zip(rows, cum):
            x['cum_core_snr'] = c
        out[case] = rows
    return out


def sdss_fibres():
    strands, _ = load_strands()
    r = list(csv.DictReader(open(DER / 'spectra_sdss_lamost_sdss_fibres.csv')))
    ra = np.array([float(x['ra']) for x in r]); dec = np.array([float(x['dec']) for x in r])
    I = np.array([float(x['I_OIII_cal']) for x in r]); e = np.array([float(x['I_OIII_cal_err']) for x in r])
    grp = np.array([f"{x['plate']}-{x['mjd']}-{x['sp']}" for x in r])
    pm = np.array([f"{x['plate']}-{x['mjd']}" for x in r]); mjd = np.array([float(x['mjd']) for x in r])
    inst = np.array([x['instrument'] for x in r])
    T = make_model('gauss028', strands)(ra, dec)
    Td = make_model('desi', strands)(ra, dec)
    shifts = [-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    TS = [make_model('gauss028', strands, ra_shift=s)(ra, dec) for s in shifts]
    good = np.isfinite(I) & np.isfinite(e) & (e > 0)

    def fit(idx, cols):
        g, gi = np.unique(grp[idx], return_inverse=True)
        X = np.zeros((len(idx), len(g))); X[np.arange(len(idx)), gi] = 1
        X = np.hstack([X, np.column_stack(cols)])
        y = I[idx]; w = 1 / e[idx]
        ok = np.ones(len(idx), bool)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(X[ok] * w[ok, None], y[ok] * w[ok], rcond=None)
            res = (y - X @ c) * w
            new = np.abs(res) < 5
            if (new == ok).all():
                break
            ok = new
        dof = max(ok.sum() - X.shape[1], 1)
        chi2 = float(np.sum(res[ok] ** 2) / dof)
        cov = np.linalg.pinv((X[ok] * w[ok, None]).T @ (X[ok] * w[ok, None])) * max(chi2, 1.0)
        k = len(cols)
        return c[-k:], np.sqrt(np.diag(cov)[-k:]), chi2, int(ok.sum())

    rows = []
    for p in sorted(set(pm), key=lambda s: float(s.split('-')[1])):
        idx = np.where((pm == p) & good)[0]
        if len(idx) < 15 or np.ptp(T[idx]) < 0.3 or (T[idx] > 0.3).sum() < 2:
            continue
        c, ce, chi2, n = fit(idx, [T[idx]])
        nulls = []
        for Ts in TS:
            if np.ptp(Ts[idx]) < 0.3:
                continue
            cc, _, _, _ = fit(idx, [T[idx], Ts[idx]])
            nulls.append(cc[1])
        en = float(np.std(nulls)) if len(nulls) >= 5 else np.nan
        ee = float(max(ce[0], en if np.isfinite(en) else 0))
        rows.append(dict(plate_mjd=p, date=iso(float(p.split('-')[1])), mjd=float(p.split('-')[1]),
                         instrument=str(inst[idx][0]), ra=float(np.median(ra[idx])), dec=float(np.median(dec[idx])),
                         n_fib=n, n_on=int((T[idx] > 0.5).sum()), A=float(c[0]), e_formal=float(ce[0]), e_null=en,
                         e=ee, snr=float(c[0] / ee), chi2nu=chi2))
    # all-plate fits for the cross-calibration
    idx = np.where(good)[0]
    cg, eg, chig, ng = fit(idx, [T[idx]])
    cd, ed, chid, nd = fit(idx, [Td[idx]])
    south = np.where(good & (dec < 20.5))[0]
    cds, eds, _, _ = fit(south, [Td[south]])
    calib = dict(all_gauss028_A_R=float(cg[0]), all_gauss028_e=float(eg[0]),
                 all_desi_scale=float(cd[0]), all_desi_scale_e=float(ed[0]),
                 south_desi_scale=float(cds[0]), south_desi_scale_e=float(eds[0]), n=ng)
    return rows, calib


def lamost():
    l = list(csv.DictReader(open(DER / 'spectra_sdss_lamost_lrs_resid_stars.csv')))
    cat = {}
    with open(RAW / 'lamost_dr12/lrs_catalogue_strandbox.csv', encoding='utf-8-sig') as fh:
        for x in csv.DictReader(fh):
            cat[x['obsid']] = (x['obsdate'], x['lmjd'])
    rows = []
    by = {}
    for x in l:
        d = cat.get(x['obsid'])
        if d is None:
            continue
        try:
            v = float(x['I_OIII5007']); e = float(x['I_OIII5007_err'])
        except ValueError:
            continue
        if not (np.isfinite(v) and np.isfinite(e) and e > 0):
            continue
        by.setdefault((x['planid'], d[0]), []).append((x['sample'], v, e))
    for (plan, date), v in sorted(by.items(), key=lambda kv: kv[0][1]):
        on = np.array([(a, b) for s, a, b in v if s == 'on']); off = np.array([(a, b) for s, a, b in v if s == 'off'])
        if len(on) < 5 or len(off) < 5:
            rows.append(dict(plan=plan, date=date, n_on=len(on), n_off=len(off), A=None, e=None, snr=None)); continue

        def wmean(a):
            # robust: median with 1.2533 sigma/sqrt(n)
            return float(np.median(a[:, 0])), float(1.2533 * 1.4826 * np.median(np.abs(a[:, 0] - np.median(a[:, 0]))) / np.sqrt(len(a)))
        mo, eo = wmean(on); mf, ef = wmean(off)
        rows.append(dict(plan=plan, date=date, n_on=len(on), n_off=len(off), A=mo - mf, e=float(np.hypot(eo, ef)),
                         snr=float((mo - mf) / np.hypot(eo, ef))))
    return rows


def desi():
    cat = {}
    with open(DER / 'spectra_desi_catalog.csv') as fh:
        for x in csv.DictReader(fh):
            cat[x['sparcl_id']] = (float(x['mean_mjd']), float(x['min_mjd']), float(x['max_mjd']))
    rows = []
    with open(DER / 'spectra_desi_fiber_table.csv') as fh:
        for x in csv.DictReader(fh):
            if x['ok'] not in ('1', 'True', 'true'):
                continue
            m = cat.get(x['sparcl_id'])
            if m is None:
                continue
            try:
                rows.append((m[0], m[1], m[2], float(x['OIII_R']), float(x['OIII_R_err']), float(x['model_resid_R'])))
            except ValueError:
                pass
    a = np.array(rows)
    a = a[np.isfinite(a).all(1)]
    month = np.array([iso(v)[:7] for v in a[:, 0]])
    out = []
    for mo in sorted(set(month)):
        k = month == mo
        if k.sum() < 30:
            continue
        y = a[k, 3]; e = a[k, 4]; x = a[k, 5]
        if np.ptp(x) < 0.3:
            continue
        ok = np.ones(k.sum(), bool)
        for _ in range(5):
            X = np.c_[np.ones(ok.sum()), x[ok]]
            c, *_ = np.linalg.lstsq(X / e[ok, None], y[ok] / e[ok], rcond=None)
            r = (y - c[0] - c[1] * x) / e
            new = np.abs(r) < 6
            if (new == ok).all():
                break
            ok = new
        X = np.c_[np.ones(ok.sum()), x[ok]] / e[ok, None]
        chi2 = np.sum(r[ok] ** 2) / max(ok.sum() - 2, 1)
        cov = np.linalg.pinv(X.T @ X) * max(chi2, 1)
        out.append(dict(month=mo, first_date=iso(a[k, 1].min()), n=int(k.sum()), slope=float(c[1]),
                        e=float(np.sqrt(cov[1, 1])), snr=float(c[1] / np.sqrt(cov[1, 1])),
                        n_pred_gt1R=int((np.abs(x) > 1).sum())))
    return out


def main():
    res = {'sdss_imaging': sdss_imaging()}
    res['sdss_fibres'], res['sdss_fibres_calib'] = sdss_fibres()
    res['lamost_lrs'] = lamost()
    res['desi'] = desi()
    # earliest >= 3 sigma per dataset
    first = {}
    for k, rows in [('sdss_imaging_' + c, res['sdss_imaging'][c]) for c in IMG_CASES] + \
                   [('sdss_fibres', res['sdss_fibres']), ('lamost_lrs', res['lamost_lrs'])]:
        ok = [r for r in rows if r.get('snr') is not None and r['snr'] >= 3]
        first[k] = (ok[0]['date'], ok[0]['snr'], ok[0].get('run', ok[0].get('plate_mjd', ok[0].get('plan')))) if ok else None
        cum = [r for r in rows if r.get('cum_core_snr') is not None and r['cum_core_snr'] >= 3]
        if cum:
            first[k + '_cumulative'] = (cum[0]['date'], cum[0]['cum_core_snr'], cum[0]['run'])
    okd = [r for r in res['desi'] if r['snr'] >= 3]
    first['desi'] = (okd[0]['month'], okd[0]['snr']) if okd else None
    res['first_3sigma'] = first
    (DER / 'verify_precovery_table.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(first, indent=1))
    print('fibre calib', res['sdss_fibres_calib'])
    for r in res['sdss_fibres'][:30]:
        print(f"{r['date']} {r['plate_mjd']:12s} {r['instrument']:4s} dec {r['dec']:5.1f} n {r['n_fib']:3d} on {r['n_on']:2d} "
              f"A {r['A']:+.2f} +- {r['e']:.2f} (formal {r['e_formal']:.2f} null {r['e_null']:.2f}) snr {r['snr']:+.1f}")
    for r in res['lamost_lrs']:
        if r['A'] is not None:
            print('LAMOST', r['date'], r['plan'], r['n_on'], r['n_off'], f"{r['A']:+.2f} +- {r['e']:.2f} snr {r['snr']:+.1f}")
    for r in res['desi']:
        print('DESI', r['month'], r['first_date'], r['n'], f"{r['slope']:.2f}+-{r['e']:.2f} snr {r['snr']:.1f}")


if __name__ == '__main__':
    main()
