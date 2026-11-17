"""verify_desi_brightness.py -- (b) absolute [O III]5007 surface brightness of the Silk ridge from DESI residuals.

Conversion (derived from desispec fluxcalibration.py + fiberfluxcorr.py, see notes/verify_desi.md):
  DESI calibrates every fibre for POINT sources: F_cal = counts_flat * FLAT_TO_PSF / C, with C fitted on standard
  stars, so for a uniform surface brightness SB, F_cal = SB * Omega_fibre / f_ps, where f_ps is the TRUE
  point-source fibre acceptance of that exposure (seeing, fibre angular size, offsets). Hence
      SB = F_cal * f_ps / Omega.
  PSF_TO_FIBER_SPECFLUX (used by wave 1 as f) is NOT f_ps: for PSF-morphology targets it is the constant 0.789
  (Legacy-Surveys 1.5"-fibre / 1"-Gaussian convention); for galaxies 0.789 * f_ps/f_gal * f_gal,nom/f_ps,nom.
  I use f_ps from the GFA (desi_dr1.exposure.fiber_fracflux_gfa), EFFTIME-weighted harmonic mean over the
  coadd's exposures, scaled from the GFA band to 5008 A with Moffat(3.5), FWHM ~ lambda^-0.2.
Forward model: own implementation (unit-peak Gaussian ridges on E/W/W2, triangular 1-deg Dec nodes, minus the
  EFFTIME-weighted tile-petal mean over all fiberassign positions), WLS with 4-sigma clipping, fitted twice:
  with the GFA conversion (mine) and with the wave-1 conversion (to separate model from calibration effects).
Model-free lower limits: wave-1 0.15-deg binning on my fluxes, with a look-elsewhere null from 12 fake lines.
Output: data/derived/verify_desi_brightness.json
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, load_meta, Exposures, Coadds, Template, trace_strands, good_fit  # noqa: E402


def node_list():
    out = []
    for s, arr in trace_strands().items():
        dmin, dmax = arr[:, 1].min(), arr[:, 1].max()
        nodes = np.unique(np.concatenate([[dmin], np.arange(np.ceil(dmin), dmax + 1e-6, 1.0), [dmax]]))
        out += [(s, j, nodes) for j in range(len(nodes))]
    return out


def tri(dec, j, nodes):
    w = np.zeros_like(dec)
    d0 = nodes[j]
    if j > 0:
        lo = nodes[j - 1]
        m = (dec >= lo) & (dec <= d0)
        w[m] = (dec[m] - lo) / (d0 - lo)
    else:
        w[dec <= d0] = 1.0
    if j < len(nodes) - 1:
        hi = nodes[j + 1]
        m = (dec > d0) & (dec <= hi)
        w[m] = (hi - dec[m]) / (hi - d0)
    else:
        w[dec > d0] = 1.0
    return w


def design(meta, co):
    nl = node_list()
    keys = co.pet_keys
    ra = np.concatenate([co.pet_ra[k] for k in keys])
    de = np.concatenate([co.pet_dec[k] for k in keys])
    lens = np.array([co.pet_ra[k].size for k in keys])
    ends = np.cumsum(lens)
    starts = ends - lens
    Tall = {s: Template(strands=(s,))(ra, de) for s in ('E', 'W', 'W2')}
    Tfib = {s: Template(strands=(s,))(meta['ra'], meta['dec']) for s in ('E', 'W', 'W2')}
    X = np.zeros((len(meta['ra']), len(nl)))
    names = []
    for k, (s, j, nodes) in enumerate(nl):
        va = Tall[s] * tri(de, j, nodes)
        cs = np.concatenate([[0.0], np.cumsum(va)])
        pm = {kk: (cs[e] - cs[b]) / (e - b) for kk, b, e in zip(keys, starts, ends)}
        sub = np.array([sum(w * pm[(t, p)] for t, p, w in c) if c else np.nan for c in co.contrib])
        X[:, k] = Tfib[s] * tri(meta['dec'], j, nodes) - sub
        names.append(f'{s}@{nodes[j]:.2f}')
    return X, names


def wls_clip(y, e, X, nsig=4.0, iters=5):
    keep = np.ones(len(y), bool)
    for _ in range(iters):
        w = 1 / e[keep] ** 2
        A = X[keep].T @ (X[keep] * w[:, None])
        act = np.diag(A) > 0
        b = np.zeros(X.shape[1])
        cov = np.full((X.shape[1], X.shape[1]), np.nan)
        Ai = np.linalg.pinv(A[np.ix_(act, act)])
        b[act] = Ai @ (X[keep][:, act].T @ (w * y[keep]))
        cov[np.ix_(act, act)] = Ai
        r = (y - X @ b) / e
        keep = np.abs(r) < nsig
    chi2 = float(np.sum(((y[keep] - X[keep] @ b) / e[keep]) ** 2))
    dof = int(keep.sum() - act.sum())
    eb = np.sqrt(np.diag(cov) * max(chi2 / dof, 1))
    return b, eb, chi2, dof, keep, act


def binned_max(ra, dec, v, e, ok, lo, hi, bs=0.15):
    xe = np.arange(113.0, 120.01, bs / np.cos(np.radians(20)))
    ye = np.arange(11.5, 28.81, bs)
    w = 1 / e[ok] ** 2
    sw, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye], weights=w)
    swv, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye], weights=w * v[ok])
    n, _, _ = np.histogram2d(ra[ok], dec[ok], [xe, ye])
    with np.errstate(invalid='ignore', divide='ignore'):
        m, em = swv / sw, 1 / np.sqrt(sw)
    yc = 0.5 * (ye[:-1] + ye[1:])
    sel = (yc >= lo) & (yc < hi)
    mm, ee, nn = m[:, sel], em[:, sel], n[:, sel]
    good = (nn >= 4) & np.isfinite(mm)
    if not good.any():
        return None
    snr = np.where(good, mm / ee, -np.inf)
    lim = np.where(good, mm - 2 * ee, -np.inf)
    k = np.argmax(lim)
    i, j = np.unravel_index(k, mm.shape)
    return dict(max_bin=float(mm[i, j]), err=float(ee[i, j]), lower2s=float(lim[i, j]), n=int(nn[i, j]),
                ra=float(0.5 * (xe[i] + xe[i + 1])), dec=float(yc[sel][j]), snr=float(mm[i, j] / ee[i, j]),
                max_snr=float(snr.max()), n_bins=int(good.sum()))


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids = lf['sparcl_id']
    keep, meta = load_meta(ids)
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    scale = det['err_scale_from_nulls']
    O, E = lf['o3__amp'][:, 0], lf['o3__err'][:, 0] * scale
    ok = np.isfinite(O) & np.isfinite(E) & (E > 0) & cm['ok_contrib'] & np.isfinite(cm['kR_5008'])
    ok &= good_fit(lf, 'o3')
    ex = Exposures()
    co = Coadds(meta, ex)
    X, names = design(meta, co)
    out = {'n_fibres': int(ok.sum()), 'nodes': names}
    conv = {'gfa_chromatic': cm['kR_5008'], 'wave1_psf_to_fiber': cm['kR_w1_5008'],
            'gfa_achromatic': cm['kR_5008'] * cm['f_achrom'] / cm['f_chrom_5008'],
            'gfa_halfchrom': cm['kR_5008'] * cm['f_half_5008'] / cm['f_chrom_5008'],
            'gfa_lamref6100': cm['kR_5008'] * cm['f_chrom_5008_lamref6100'] / cm['f_chrom_5008'],
            'gfa_lamref6700': cm['kR_5008'] * cm['f_chrom_5008_lamref6700'] / cm['f_chrom_5008']}
    fits = {}
    for lab, k in conv.items():
        y, e = O * k, E * k
        b, eb, chi2, dof, kp, act = wls_clip(y[ok], e[ok], X[ok])
        nodes = [dict(node=n, amp_R=float(a), err_R=float(s)) for n, a, s, ac in zip(names, b, eb, act) if ac]
        good_ = [d for d in nodes if d['err_R'] < 1.0]          # drop unconstrained end nodes
        south = [d['amp_R'] for d in good_ if float(d['node'].split('@')[1]) < 19.5 and
                 not d['node'].startswith('W2')]
        north = [d['amp_R'] for d in good_ if float(d['node'].split('@')[1]) >= 21.5]   # E and W2, as wave 1
        fits[lab] = dict(nodes=nodes, chi2=chi2, dof=dof, n_clipped=int((~kp).sum()),
                         south_median=float(np.median(south)), north_median=float(np.median(north)))
    out['forward_model'] = fits
    # model sky SB (before petal subtraction) and measured residuals at the paper's P1-P4
    from silk_geometry import REGIONS
    nl = node_list()
    pra = np.array([v[0] for v in REGIONS.values()])
    pde = np.array([v[1] for v in REGIONS.values()])
    Tp = {s: Template(strands=(s,))(pra, pde) for s in ('E', 'W', 'W2')}
    Bp = np.array([Tp[s] * tri(pde, j, nodes) for (s, j, nodes) in nl]).T
    reg = {}
    for lab in ('gfa_chromatic', 'wave1_psf_to_fiber'):
        b = np.array([0.0 if not np.isfinite(d) else d for d in
                      [next((q['amp_R'] for q in fits[lab]['nodes'] if q['node'] == n), np.nan) for n in names]])
        for k, nm in enumerate(REGIONS):
            reg.setdefault(nm, {})[f'model_SB_{lab}'] = float(Bp[k] @ b)
    for k, nm in enumerate(REGIONS):
        dd = np.hypot((meta['ra'] - pra[k]) * np.cos(np.radians(pde[k])), meta['dec'] - pde[k])
        for rad in (0.1, 0.2):
            m = ok & (dd < rad)
            if m.sum() >= 3:
                v_, e_ = O[m] * cm['kR_5008'][m], E[m] * cm['kR_5008'][m]
                ww = 1 / e_ ** 2
                reg[nm][f'resid_R_r{rad}'] = [float(np.sum(ww * v_) / ww.sum()), float(ww.sum() ** -0.5), int(m.sum())]
    out['regions'] = reg
    # ratio of my calibration to wave 1 for the fibres that carry the ridge signal (|Z| weighted)
    Z = cm['Z']
    w = np.abs(Z) * ok
    r = cm['kR_5008'] / cm['kR_w1_5008']
    out['conversion_ratio_gfa_over_wave1'] = dict(
        median_all=float(np.nanmedian(r[ok])), p16=float(np.nanpercentile(r[ok], 16)),
        p84=float(np.nanpercentile(r[ok], 84)), Zweighted_mean=float(np.nansum(w * r) / np.nansum(w)),
        positive_ridge_fibres_median=float(np.nanmedian(r[ok & (Z > 0.3)])))
    # ---- model-free lower limits: wave-1 binning, my fluxes and conversion
    bands = ((12, 15.5), (15.5, 17.5), (17.5, 19.5), (19.5, 22), (22, 25), (25, 28.8))
    lims = {}
    for lab, k in (('gfa', cm['kR_5008']), ('wave1conv', cm['kR_w1_5008'])):
        v, e = O * k, E * k
        lims[lab] = {f'{lo}-{hi}': binned_max(meta['ra'], meta['dec'], v, e, ok, lo, hi) for lo, hi in bands}
    # look-elsewhere: max bin S/N in each band for the fake lines (same fibres, same conversion)
    fake = [g[:-5] for g in lf.files if g.startswith('null') and g.endswith('__amp')]
    le = {}
    for lo, hi in bands:
        mx = []
        for g in fake:
            a, e = lf[g + '__amp'][:, 0], lf[g + '__err'][:, 0] * scale
            m = ok & good_fit(lf, g)
            bm = binned_max(meta['ra'], meta['dec'], a * cm['kR_5008'], e * cm['kR_5008'], m, lo, hi)
            if bm:
                mx.append(bm['max_snr'])
        le[f'{lo}-{hi}'] = dict(fake_max_snr_median=float(np.median(mx)) if mx else None,
                                fake_max_snr_max=float(np.max(mx)) if mx else None, n_fake=len(mx))
    out['lower_limits'] = lims
    out['look_elsewhere_fake_lines'] = le
    # reproduce the wave-1 bin at (116.27, +17.28): fibres within its 0.15-deg cell
    cx, cy = 116.2723, 17.275
    hw = 0.15 / np.cos(np.radians(20)) / 2
    cell = ok & (np.abs(meta['ra'] - cx) <= hw) & (np.abs(meta['dec'] - cy) <= 0.075)
    for lab, k in (('gfa', cm['kR_5008']), ('wave1conv', cm['kR_w1_5008'])):
        v, e = O[cell] * k[cell], E[cell] * k[cell]
        ww = 1 / e ** 2
        out.setdefault('wave1_bin_116.27_17.28', {})[lab] = dict(
            n=int(cell.sum()), mean_R=float(np.sum(ww * v) / ww.sum()), err_R=float(ww.sum() ** -0.5),
            unweighted_median_R=float(np.median(v)), min_R=float(v.min()), max_R=float(v.max()))
    out['wave1_bin_116.27_17.28']['f_true_median'] = float(np.median(cm['f_chrom_5008'][cell]))
    out['wave1_bin_116.27_17.28']['psf_to_fiber_median'] = float(np.median(meta['mean_psf_to_fiber_specflux'][cell]))
    out['wave1_bin_116.27_17.28']['programs'] = sorted(set(meta['program'][cell]))
    # far-field zero point: [O III] residual of control fibres > 5 deg from the trace
    far = ok & (meta['d_any'] > 5)
    if far.sum() > 10:
        v, e = O[far] * cm['kR_5008'][far], E[far] * cm['kR_5008'][far]
        ww = 1 / e ** 2
        out['control_gt5deg'] = dict(n=int(far.sum()), wmean_R=float(np.sum(ww * v) / ww.sum()),
                                     err_R=float(ww.sum() ** -0.5), median_R=float(np.median(v)),
                                     chi_std=float(np.std(v / e)), frac_gt3sig=float(np.mean(v / e > 3)),
                                     frac_lt_m3sig=float(np.mean(v / e < -3)))
    (D / 'verify_desi_brightness.json').write_text(json.dumps(out, indent=1))
    for lab, f in fits.items():
        print(lab, 'south median %.2f north median %.2f chi2/dof %.0f/%d' % (f['south_median'], f['north_median'],
                                                                          f['chi2'], f['dof']))
        print('   ', ' '.join(f"{d['node']}:{d['amp_R']:.2f}+-{d['err_R']:.2f}" for d in f['nodes']))
    print(json.dumps(out['conversion_ratio_gfa_over_wave1']))
    print(json.dumps(lims, indent=0)[:3000])
    print(json.dumps(le))
    print(json.dumps(out['wave1_bin_116.27_17.28']))
    print(json.dumps(out.get('control_gt5deg')))


if __name__ == '__main__':
    main()
