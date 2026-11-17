"""verify_sdss: absolute per-fibre [O III]5007 near P1-P3 and the brightest fibres (diffuse-corrected).

I_f = f * ( I_abs(group mean total sky) + Delta_f - <Delta>_group ), with Delta_f the fibre's [O III] relative to the
group OFF mean (verify_sdss_diff_<tag>.npz, from flux+sky) and I_abs from the group-mean absolute fit (solar template
+ lines; template systematic of order 0.2 R).  Writes data/derived/verify_sdss_pregions.json
Usage: python -I verify_sdss_pregions.py <tag> f_sdss f_boss
"""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT
from verify_sdss_results import load, good_groups
from silk_geometry import REGIONS
from astropy.coordinates import SkyCoord
import astropy.units as u


def main(tag, f_sdss, f_boss):
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    a = amp['OIII:OIII5007']; e = amp['OIII:OIII5007:err']
    g = fib['gid']
    absg = grp['I_OIII5007_abs_cal']
    mean_d = np.full(len(absg), np.nan)
    for k in np.where(gok)[0]:
        s = (g == k) & np.isfinite(a)
        if s.sum() >= 3:
            mean_d[k] = np.mean(a[s])
    I = (absg[g] + a - mean_d[g]) * fac[g]
    E = e * fac[g]
    ok = gok[g] & np.isfinite(I)
    c = SkyCoord(fib['ra'] * u.deg, fib['dec'] * u.deg)
    out = {}
    for p, (ra, dec) in REGIONS.items():
        sep = c.separation(SkyCoord(ra * u.deg, dec * u.deg)).deg
        s = ok & (sep < 0.3)
        if s.sum():
            w = 1 / E[s] ** 2
            j = np.argmax(I[s])
            out[p] = {'n': int(s.sum()), 'mean_R': float(np.sum(w * I[s]) / np.sum(w)), 'err_R': float(1 / np.sqrt(np.sum(w))),
                      'median_R': float(np.median(I[s])), 'max_R': float(I[s][j]), 'max_err_R': float(E[s][j])}
    on = ok & (fib['dmin'] < 0.15)
    top = np.argsort(-np.where(on, I, -np.inf))[:10]
    out['brightest_on_fibres'] = [{'ra': float(fib['ra'][i]), 'dec': float(fib['dec'][i]), 'I_R': float(I[i]), 'err_R': float(E[i]),
                                   'instrument': str(fib['instrument'][i]), 'mjd': int(fib['mjd'][i])} for i in top]
    off = ok & (fib['dmin'] > 1.0) & (fib['dmin'] < 3.0)
    out['offstrand_1to3deg_median_R'] = float(np.median(I[off]))
    out['offstrand_1to3deg_p16_p84'] = [float(np.percentile(I[off], 16)), float(np.percentile(I[off], 84))]
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_pregions.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
