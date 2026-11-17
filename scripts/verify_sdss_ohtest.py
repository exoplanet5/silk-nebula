"""verify_sdss: are OH residuals or geocoronal Halpha responsible for the high north-strand [N II]/Halpha, [S II]/Halpha?

Per group, ON-OFF total-sky spectra (R/A, diffuse-corrected) are shifted to the TOPOCENTRIC frame (lam / (1 + v_corr/c),
v_corr = -HELIO_RV) and stacked; OH lines are then sharp and nebular lines are smeared by the spread of v_corr (fitted
with a broadened Gaussian).  Fit: Halpha, [N II]6548/6583, [S II], OH 6555.4, 6579.1, 6535.0 (free), linear baseline.
The ON-OFF OH amplitudes are compared with the absolute OH intensities (fraction = spatial OH structure leaking into
the differential).  Writes data/derived/verify_sdss_ohtest.json
Usage: python -I verify_sdss_ohtest.py <tag> f_sdss f_boss
"""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, VAC, C, load_windows, flam_to_R_per_A, fit_lines
from verify_sdss_results import load, good_groups

OHV = {'OH6535': 6535.0, 'OH6555': 6555.4, 'OH6579': 6579.1, 'OH6606': 6606.2}


def main(tag, f_sdss, f_boss):
    D = load_windows(); lam = D['lam']
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    dec = fib['dec']
    ON = fib['dmin'] < 0.15
    OFF = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    grid = np.arange(6520.0, 6620.0, 0.5)
    out = {}
    for nm, sel in (('south', ON & (dec < 20)), ('north', ON & (dec >= 20)), ('patch', fib['inpatch30'])):
        num = np.zeros(grid.size); den = np.zeros(grid.size); tot_abs = np.zeros(grid.size); vc = []; ws = []
        for g in np.where(gok)[0]:
            s = np.where(fib['gid'] == g)[0]
            so = s[sel[s]]; sf = s[(OFF & ~sel)[s]]
            if len(so) < 1 or len(sf) < 3:
                continue
            m = (lam > 6510) & (lam < 6630)
            X = flam_to_R_per_A(D['tot'][s][:, m].astype(float), lam[m], fib['diam'][s[0]]) * fac[g]
            X[D['ivar'][s][:, m] <= 0] = np.nan
            d = np.nanmean(X[np.isin(s, so)], 0) - np.nanmean(X[np.isin(s, sf)], 0)
            a = np.nanmean(X, 0)
            v = -float(grp['helio_rv'][g])
            lt = lam[m] / (1 + v / C)
            ok = np.isfinite(d)
            w = len(so) * len(sf) / (len(so) + len(sf))
            num += w * np.interp(grid, lt[ok], d[ok]); den += w
            tot_abs += w * np.interp(grid, lt[ok], a[ok])
            vc.append(v); ws.append(w)
        sp = num / den; ab = tot_abs / den
        vc = np.array(vc); ws = np.array(ws)
        vrms = float(np.sqrt(np.sum(ws * vc ** 2) / np.sum(ws) - (np.sum(ws * vc) / np.sum(ws)) ** 2))
        vmean = float(np.sum(ws * vc) / np.sum(ws))
        # nebular lines in the topocentric frame sit at V_helio(+13) - v_corr -> mean offset, smeared by vrms
        s_lsf = 1.4
        s_neb = np.sqrt(s_lsf ** 2 + (6564.6 * vrms / C) ** 2)
        lines_neb = [('Ha', VAC['Ha'] * (1 + (13 - vmean) / C)), ('NII6548', VAC['NII6548'] * (1 + (13 - vmean) / C)),
                     ('NII6583', VAC['NII6583'] * (1 + (13 - vmean) / C))]
        cols = [np.ones_like(grid), (grid - 6570) / 50]
        names = ['c0', 'c1']
        for n, l0 in lines_neb:
            if n == 'NII6548':
                continue
            g_ = np.exp(-0.5 * ((grid - l0) / s_neb) ** 2) / (np.sqrt(2 * np.pi) * s_neb)
            if n == 'NII6583':
                l2 = lines_neb[1][1]
                g_ = g_ + np.exp(-0.5 * ((grid - l2) / s_neb) ** 2) / (np.sqrt(2 * np.pi) * s_neb) / 3.06
            cols.append(g_); names.append(n)
        for n, l0 in OHV.items():
            cols.append(np.exp(-0.5 * ((grid - l0) / s_lsf) ** 2) / (np.sqrt(2 * np.pi) * s_lsf)); names.append(n)
        A = np.array(cols).T
        ok = np.isfinite(sp)
        c, *_ = np.linalg.lstsq(A[ok], sp[ok], rcond=None)
        r = sp[ok] - A[ok] @ c
        cov = np.linalg.pinv(A[ok].T @ A[ok]) * np.var(r) * 1.0
        # absolute OH amplitudes from the mean total sky (same model, no nebular lines needed for the ratio)
        ca, *_ = np.linalg.lstsq(A[ok], ab[ok], rcond=None)
        res = {n: (float(c[i]), float(np.sqrt(cov[i, i]))) for i, n in enumerate(names)}
        res['abs'] = {n: float(ca[i]) for i, n in enumerate(names) if n.startswith('OH') or n == 'Ha'}
        res['OH_fraction'] = {n: float(c[i] / ca[i]) for i, n in enumerate(names) if n.startswith('OH')}
        # [N II] fit without OH terms (to see how much OH matters)
        keep = [i for i, n in enumerate(names) if not n.startswith('OH')]
        c2, *_ = np.linalg.lstsq(A[ok][:, keep], sp[ok], rcond=None)
        res['without_OH_terms'] = {names[i]: float(c2[j]) for j, i in enumerate(keep) if names[i] in ('Ha', 'NII6583')}
        res['vcorr_mean'] = vmean; res['vcorr_rms'] = vrms; res['n_groups'] = int(len(vc))
        out[nm] = res
        np.save(ROOT / f'data/derived/verify_sdss_ohtest_{nm}.npy', np.vstack([grid, sp, ab]))
        print(nm, {k: (round(v[0], 3), round(v[1], 3)) for k, v in res.items() if isinstance(v, tuple)}, 'OH frac', {k: round(v, 4) for k, v in res['OH_fraction'].items()},
              'noOH', {k: round(v, 3) for k, v in res['without_OH_terms'].items()}, 'vcorr', round(vmean, 1), round(vrms, 1))
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_ohtest.json', 'w'), indent=1)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
