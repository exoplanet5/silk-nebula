"""SDSS drift-scan sky levels as a broadband [OIII] photometer across the Silk Nebula.

Idea: the SDSS photo pipeline records, for every 13.5'x9.8' field, the median sky surface brightness
(Field.sky_g, Field.sky_r in nanomaggies/arcsec^2, DR18 SkyServer).  Diffuse emission on degree scales is part of that
'sky'.  SDSS Legacy/SEGUE stripes here run at PA ~47 deg, i.e. obliquely across the N-S filament, and the six camera
columns (observed simultaneously, offset 0.42 deg perpendicular to the scan) cross the [OIII] ridge at different times.
Per run we therefore fit
      sky[f, c] = S[f] + g[c] + h[c] (f - f0) + A * T[f, c]
  S[f]  free common-mode term per field number (same time in all camcols; absorbs airglow / moon / twilight changes),
  g, h  per-camcol offset and linear drift (flat-field / scattered light),
  T     NSNS DR0.2 [OIII] (or Halpha) map, smoothed to the field size (Gaussian sigma 4'), at the field centre, in R.
A (nmgy arcsec^-2 per NSNS-R) is the broadband response to the NSNS template; 4-sigma iterative clipping.
Errors: (i) formal from residual scatter, inflated by sqrt(2) for the pairwise-repeated sky_g values;
        (ii) placebo: refit with the template sampled at RA offsets of +-0.8..2.0 deg / cos(dec) (off-filament).
Expected A for an exactly calibrated NSNS map: [OIII]5007 1 R (+4959 0.33 R) -> 1.80e-3 nmgy arcsec^-2 in g;
Halpha 1 R -> 1.88e-3 in r (imaging_optical_expected.py, SDSS 2010 response curves).
Input : data/raw/sdss_dr18/fields_box.csv (SkyServer DR18 Field table, see fields_box.json for the query URL)
        data/derived/imaging_optical_nsns_maps.npz, imaging_optical_nsns_ridge.json
Output: data/derived/imaging_optical_sdss_sky.json, imaging_optical_sdss_sky_fields.csv
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from imaging_optical_common import RAW, DER, save_json, along_across
from imaging_optical_nsns import ridge_spine

EXPECT = {'g': 1.80e-3, 'r': 1.88e-3}   # nmgy arcsec^-2 per R ([OIII]5007+4959 for g, Halpha for r)
PLACEBO_SHIFTS = [s for s in np.arange(-2.0, 2.01, 0.2) if abs(s) >= 0.8]


def load_fields():
    a = np.genfromtxt(RAW / 'sdss_dr18/fields_box.csv', delimiter=',', names=True, skip_header=1)
    return a


_TCACHE = {}


def template_sampler(chan, sigma_pix=4.0, variant='raw'):
    """variant: 'raw' = NSNS map Gaussian-smoothed (sigma 4');
    'medfilt' = 11'x11' median filter first (removes star residuals / point sources, keeps the ~20' wide filament);
    'ridge' = medfilt, tapered to zero beyond 1.0 deg from the [OIII] ridge (pure-filament template)."""
    key = (chan, sigma_pix, variant)
    if key in _TCACHE:
        return _TCACHE[key]
    from scipy.ndimage import median_filter
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    w = WCS(str(z['wcs_header']))
    base = np.nan_to_num(z[chan])
    if variant in ('medfilt', 'ridge'):
        base = median_filter(base, size=11)
    img = gaussian_filter(base, sigma_pix)
    if variant == 'ridge':
        sra, sdec = ridge_spine()
        acr, alo, bey = along_across(z['ra'], z['dec'], sra, sdec)
        taper = np.clip((1.2 - np.abs(acr)) / 0.4, 0, 1)
        img = img - np.nanmedian(img[(np.abs(acr) > 1.0) & (np.abs(acr) < 1.6)])
        img = img * taper
    def sample(ra, dec):
        x, y = w.wcs_world2pix(ra, dec, 0)
        ix = np.round(x).astype(int); iy = np.round(y).astype(int)
        ok = (ix >= 0) & (ix < img.shape[1]) & (iy >= 0) & (iy < img.shape[0])
        out = np.full(np.shape(ra), np.nan)
        out[ok] = img[iy[ok], ix[ok]]
        return out
    _TCACHE[key] = sample
    return sample


def bright_star_fields(a, gmax=7.0):
    """Boolean: field centre within 0.25 deg of a G<5 star or 0.12 deg of a 5<=G<7 star (Gaia DR3, Data Lab)."""
    import glob
    from scipy.spatial import cKDTree
    rows = []
    for f in sorted(glob.glob(str(RAW / 'gaia_dr3_datalab/gaia_g175_dec*.csv'))):
        g = np.genfromtxt(f, delimiter=',', names=True)
        m = g['phot_g_mean_mag'] < gmax
        rows.append(np.c_[g['ra'][m], g['dec'][m], g['phot_g_mean_mag'][m]])
    st = np.vstack(rows)
    def xyz(ra, de):
        ra = np.radians(ra); de = np.radians(de)
        return np.c_[np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)]
    t = cKDTree(xyz(a['ra'], a['dec']))
    bad = np.zeros(len(a), bool)
    for ra, de, G in st:
        r = 0.25 if G < 5 else 0.12
        for i in t.query_ball_point(xyz(np.array([ra]), np.array([de]))[0], np.radians(r)):
            bad[i] = True
    return bad, len(st)


def fit_run(sub, T, clip=4.0, niter=5):
    """sub: structured rows of one run; T: template values. Returns A, formal err, n_used, residual rms."""
    f = sub['field'].astype(int); c = sub['camcol'].astype(int); y = sub['sky']
    ok = np.isfinite(y) & np.isfinite(T) & (y > 0)
    for _ in range(niter):
        fu = np.unique(f[ok]); cu = np.unique(c[ok])
        if len(cu) < 3 or ok.sum() < 3 * len(fu) or len(fu) < 8:
            return None
        fi = {v: i for i, v in enumerate(fu)}; ci = {v: i for i, v in enumerate(cu)}
        n = ok.sum(); nf = len(fu); nc = len(cu)
        X = np.zeros((n, nf + 2 * nc))     # S_f, g_c (first camcol fixed by S), h_c, + A
        rows = np.where(ok)[0]
        f0 = np.median(f[ok])
        for k, i in enumerate(rows):
            X[k, fi[f[i]]] = 1.0
            if ci[c[i]] > 0:
                X[k, nf + ci[c[i]]] = 1.0
            X[k, nf + nc + ci[c[i]]] = (f[i] - f0) / 100.0
        X = np.delete(X, nf, axis=1)        # drop g_c of the reference camcol
        X = np.c_[X, T[rows]]
        coef, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
        res = y[rows] - X @ coef
        s = 1.4826 * np.median(np.abs(res - np.median(res)))
        bad = np.abs(res) > clip * max(s, 1e-6)
        if not bad.any():
            break
        ok[rows[bad]] = False
    dof = max(len(rows) - X.shape[1], 1)
    s2 = np.sum(res ** 2) / dof
    try:
        cov = s2 * np.linalg.pinv(X.T @ X)
    except Exception:
        return None
    # leverage of the template after removing the nuisance terms
    eA = float(np.sqrt(max(cov[-1, -1], 0)) * np.sqrt(2.0))
    return dict(A=float(coef[-1]), eA_formal=eA, n=int(len(rows)), n_fields=int(nf), rms=float(np.sqrt(s2)),
                T_range=float(np.nanmax(T[rows]) - np.nanmin(T[rows])))


def run_band(a, band, chan, sra, sdec, ridge, variant='raw', exclude=None, mjd_range=None):
    samp = template_sampler(chan, variant=variant)
    sky = a[f'sky_{band}']
    out_runs = []
    rows_csv = []
    across, along, beyond = along_across(a['ra'], a['dec'], sra, sdec)
    for run in np.unique(a['run']):
        s = (a['run'] == run) & (np.abs(across) <= 2.5) & ~beyond & (a['quality'] >= 2)
        if exclude is not None:
            s &= ~exclude
        if mjd_range is not None and not (mjd_range[0] <= np.nanmedian(a['mjd_g'][a['run'] == run]) < mjd_range[1]):
            continue
        if s.sum() < 40:
            continue
        sub = np.zeros(s.sum(), dtype=[('field', 'i4'), ('camcol', 'i4'), ('sky', 'f8')])
        sub['field'] = a['field'][s]; sub['camcol'] = a['camcol'][s]; sub['sky'] = sky[s]
        T = samp(a['ra'][s], a['dec'][s])
        # require that the run actually crosses the ridge with fields on both sides
        acr = across[s]
        if not ((acr < -0.7).sum() >= 6 and (acr > 0.7).sum() >= 6 and (np.abs(acr) < 0.3).sum() >= 4):
            continue
        r = fit_run(sub, T)
        if r is None:
            continue
        pl = []
        for sh in PLACEBO_SHIFTS:
            Ts = samp(a['ra'][s] + sh / np.cos(np.radians(a['dec'][s])), a['dec'][s])
            rp = fit_run(sub, Ts)
            if rp is not None:
                pl.append(rp['A'])
        r.update(run=int(run), mjd=float(np.nanmedian(a['mjd_g'][s])), dec_cross=float(np.median(a['dec'][s][np.abs(acr) < 0.3])),
                 A_placebo_std=float(np.std(pl)) if len(pl) >= 5 else None, n_placebo=len(pl),
                 sky_median=float(np.median(sky[s])))
        out_runs.append(r)
        for i in np.where(s)[0][:0]:
            pass
    return out_runs


def combine(runs, key_err='eA'):
    A = np.array([r['A'] for r in runs]); e = np.array([r[key_err] for r in runs])
    w = 1 / e ** 2
    m = np.sum(w * A) / np.sum(w); em = 1 / np.sqrt(np.sum(w))
    chi2 = float(np.sum(w * (A - m) ** 2)); dof = len(A) - 1
    # scale error by sqrt(chi2/dof) if > 1
    sc = max(1.0, np.sqrt(chi2 / max(dof, 1)))
    return dict(A=float(m), eA=float(em * sc), chi2_dof=float(chi2 / max(dof, 1)), n_runs=len(A))


def summarize(runs, band):
    for r in runs:
        r['eA'] = max(r['eA_formal'], r['A_placebo_std'] or 0.0)
    comb = combine(runs) if runs else None
    bins = {}
    for d0, d1 in [(9, 14), (14, 17), (17, 20.5), (20.5, 24), (24, 30)]:
        rr = [r for r in runs if d0 <= r['dec_cross'] < d1]
        if len(rr) >= 2:
            bins[f'{d0}-{d1}'] = combine(rr)
    core = [r for r in runs if 12 <= r['dec_cross'] < 20.5]
    cc = combine(core) if len(core) >= 2 else None
    out = dict(combined=comb, by_dec=bins, core_dec12_20p5=cc)
    for k in ('combined', 'core_dec12_20p5'):
        if out[k]:
            out[k]['A_over_expected'] = out[k]['A'] / EXPECT[band]
            out[k]['err_over_expected'] = out[k]['eA'] / EXPECT[band]
    for k, v in bins.items():
        v['A_over_expected'] = v['A'] / EXPECT[band]; v['err_over_expected'] = v['eA'] / EXPECT[band]
    if core:
        x = np.array([r['A'] / r['eA'] for r in core])
        out['core_n_positive'] = int((x > 0).sum()); out['core_n'] = len(core)
        out['core_median_A_over_expected'] = float(np.median([r['A'] for r in core]) / EXPECT[band])
    return out


def main():
    a = load_fields()
    sra, sdec = ridge_spine()
    ridge = json.loads((DER / 'imaging_optical_nsns_ridge.json').read_text())
    bad, nst = bright_star_fields(a)
    res = {'method': __doc__.strip().split('\n')[0], 'expected_A_if_NSNS_calibrated': EXPECT,
           'bright_star_exclusion': f'{int(bad.sum())} of {len(a)} fields within 0.25 deg (G<5) / 0.12 deg (G<7) of {nst} Gaia stars',
           'cases': {}}
    tests = [('g_vs_OIII', 'g', 'oiii', 'raw', None, None),
             ('g_vs_OIII_medfilt', 'g', 'oiii', 'medfilt', None, None),
             ('g_vs_OIII_ridge', 'g', 'oiii', 'ridge', None, None),
             ('g_vs_OIII_ridge_nostars', 'g', 'oiii', 'ridge', bad, None),
             ('g_vs_OIII_ridge_mjd_lt_53200', 'g', 'oiii', 'ridge', bad, (0, 53200)),
             ('g_vs_OIII_ridge_mjd_ge_53200', 'g', 'oiii', 'ridge', bad, (53200, 1e9)),
             ('r_vs_Halpha', 'r', 'halpha', 'medfilt', bad, None),
             ('r_vs_OIII_ridge', 'r', 'oiii', 'ridge', bad, None),
             ('g_vs_Halpha', 'g', 'halpha', 'medfilt', bad, None)]
    for label, band, chan, var, exc, mr in tests:
        runs = run_band(a, band, chan, sra, sdec, ridge, variant=var, exclude=exc, mjd_range=mr)
        summ = summarize(runs, band)
        res['cases'][label] = dict(band=band, template=chan, variant=var, star_excluded=exc is not None, mjd_range=mr,
                                   runs=runs, **summ)
        c = summ['core_dec12_20p5']; t = summ['combined']
        print(f"{label:32s} all: {t['A_over_expected']:.2f}+-{t['err_over_expected']:.2f} (n={t['n_runs']}, chi2/dof {t['chi2_dof']:.1f})"
              f"  core Dec12-20.5: {c['A_over_expected']:.2f}+-{c['err_over_expected']:.2f} (n={c['n_runs']}, chi2/dof {c['chi2_dof']:.1f}, "
              f"pos {summ.get('core_n_positive')}/{summ.get('core_n')}, median {summ.get('core_median_A_over_expected', np.nan):.2f})")
        for k, v in summ['by_dec'].items():
            print('      ', k, f"{v['A_over_expected']:.2f} +- {v['err_over_expected']:.2f}  n={v['n_runs']}")
    save_json(res, DER / 'imaging_optical_sdss_sky.json')


if __name__ == '__main__':
    main()
