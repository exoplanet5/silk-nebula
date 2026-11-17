"""Fit the [OIII] strands of the Silk Nebula in every registered image.

For each Dec bin (0.1 deg) and each source map, an RA profile (background = running median over 1.6 deg RA)
is fitted with a sum of Gaussians, one per strand expected at that Dec (east/main strand E, west strand W),
initialised from a smooth prior track and constrained to +-0.25 deg around it. Two passes: pass 2 uses
smoothing splines of the pass-1 centres as priors.
Output: data/derived/geometry_strand_fits.csv
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_trace_fit.py
"""
import csv
import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares
from scipy.interpolate import UnivariateSpline

ROOT = '/Users/mickey/astronomy/silk-nebula'
z = np.load(f'{ROOT}/data/derived/geometry_maps.npz')
ra, de = z['ra_grid'], z['dec_grid']
STEP = abs(ra[1] - ra[0])
SRC = {'swons': z['swons_B'] - z['swons_R'], 'fig4': -z['fig4_grey'], 'nsnsN': z['nsnsNs_OIII'], 'nsnsW': z['nsnsW_OIII']}

# visual priors (RA, Dec) read off the registered NSNS-W [OIII] map; refined in pass 2
PRIOR = {
    # E: main (long) strand; north of Dec 22 and south of Dec 13.5 positions come from the NSNS-wide ridge peaks
    'E': np.array([(116.30, 12.0), (116.34, 12.7), (116.37, 13.3), (116.45, 14.0), (116.6, 15), (116.72, 16),
                   (116.82, 17), (116.92, 18), (116.95, 19), (116.9, 20), (116.72, 21.2), (116.62, 21.7),
                   (116.5, 22.1), (116.34, 22.7), (116.12, 23.3), (115.98, 23.9), (115.6, 24.7), (115.25, 25.5),
                   (114.95, 26.2), (114.67, 26.9), (114.47, 27.6), (114.3, 28.2), (114.2, 28.9)]),
    # W: bright western strand that bends west and ends near the Halpha patch
    'W': np.array([(116.42, 15.0), (116.4, 16), (116.35, 17), (116.25, 18), (116.1, 19), (116.0, 20), (115.9, 20.6)]),
    # W2: faint far-western strand seen in NSNS-wide [OIII] (Dec ~21.4-24)
    'W2': np.array([(115.18, 21.4), (115.22, 22.5), (115.25, 23.5), (115.3, 24.0)]),
}
DRANGE = {'E': (12.0, 28.9), 'W': (15.0, 20.6), 'W2': (21.4, 24.0)}


def prior_ra(strand, d, tracks=None):
    if tracks is not None and strand in tracks:
        return float(tracks[strand](d))
    p = PRIOR[strand]
    return float(np.interp(d, p[:, 1], p[:, 0]))


def gmodel(x, pars, n):
    y = np.full_like(x, pars[-2]) + pars[-1] * (x - x.mean())
    for k in range(n):
        a, mu, s = pars[3 * k:3 * k + 3]
        y += a * np.exp(-0.5 * ((x - mu) / s) ** 2)
    return y


def fit_bin(name, m, d0, dbin, tracks=None):
    j = (de >= d0 - dbin / 2) & (de < d0 + dbin / 2)
    with np.errstate(all='ignore'):
        p = np.nanmean(m[j], axis=0)
    strands = [s for s in ('E', 'W', 'W2') if DRANGE[s][0] <= d0 <= DRANGE[s][1]]
    mus = [prior_ra(s, d0, tracks) for s in strands]
    lo, hi = min(mus) - 0.7, max(mus) + 0.7
    sel = np.isfinite(p) & (ra >= lo) & (ra <= hi)
    if sel.sum() < 0.8 * ((hi - lo) / STEP):
        return []
    good = np.isfinite(p)
    pp = np.interp(np.arange(len(p)), np.where(good)[0], p[good])
    bg = ndimage.median_filter(pp, size=int(1.6 / STEP))
    x = ra[sel]; y = (pp - bg)[sel]
    cosd = np.cos(np.radians(d0))
    noise = 1.4826 * np.median(np.abs(np.diff(y))) / np.sqrt(2)
    p0, lb, ub = [], [], []
    for mu in mus:
        amp0 = max(np.interp(mu, x[::-1], y[::-1]), noise)
        p0 += [amp0, mu, 0.08 / cosd]; lb += [0, mu - 0.25, 0.015 / cosd]; ub += [np.inf, mu + 0.25, 0.35 / cosd]
    p0 += [0.0, 0.0]; lb += [-np.inf, -np.inf]; ub += [np.inf, np.inf]
    r = least_squares(lambda q: (gmodel(x, q, len(mus)) - y) / noise, p0, bounds=(lb, ub))
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, 2 * r.cost / max(1, len(x) - len(p0)))
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(len(p0), np.nan)
    out = []
    for k, s in enumerate(strands):
        a, mu, sg = r.x[3 * k:3 * k + 3]
        out.append(dict(source=name, strand=s, dec=round(d0, 3), ra=mu, ra_err=err[3 * k + 1], amp_snr=a / noise,
                        fwhm_sky_deg=2.3548 * sg * cosd, at_bound=bool(abs(mu - mus[k]) > 0.245)))
    return out


def run(tracks=None, dbin=0.1):
    res = []
    for name, m in SRC.items():
        for d0 in np.arange(12.0, 28.91, dbin):
            res += fit_bin(name, m, d0, dbin, tracks)
    return res


def smooth_tracks(res, min_snr=5):
    tr = {}
    for s in ('E', 'W', 'W2'):
        r = [q for q in res if q['strand'] == s and q['amp_snr'] > min_snr and not q['at_bound'] and np.isfinite(q['ra_err'])]
        d = np.array([q['dec'] for q in r]); x = np.array([q['ra'] for q in r])
        o = np.argsort(d); d, x = d[o], x[o]
        # robust binned median then spline
        bins = np.arange(DRANGE[s][0], DRANGE[s][1] + 0.26, 0.25)
        bc, bm = [], []
        for b0, b1 in zip(bins[:-1], bins[1:]):
            k = (d >= b0) & (d < b1)
            if k.sum() >= 2:
                bc.append(np.median(d[k])); bm.append(np.median(x[k]))
        spl = UnivariateSpline(bc, bm, k=3, s=len(bc) * 0.03 ** 2)
        tr[s] = spl
    return tr


if __name__ == '__main__':
    r1 = run()
    tr = smooth_tracks(r1)
    r2 = run(tr)
    with open(f'{ROOT}/data/derived/geometry_strand_fits.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(r2[0].keys()))
        w.writeheader()
        for q in r2:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in q.items()})
    for s in ('E', 'W', 'W2'):
        q = [x for x in r2 if x['strand'] == s and x['amp_snr'] > 5 and not x['at_bound']]
        print(s, len(q), 'fits; median fwhm %.3f deg' % np.median([x['fwhm_sky_deg'] for x in q]))
