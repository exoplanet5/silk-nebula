"""w3_oii_sdss_wavecal.py -- (wave 3, tag oii) BOSS wavelength zero point and LSF width at the blue end (3650-4360 A),
where the [O II] 3726/3729 doublet sits, from street-lamp Hg I lines in the total sky of BOSS sky fibres.

For a 2.8-sigma-separated doublet (BOSS at 3727 A) the fitted 3729/3726 ratio moves by ~0.01 per km/s of centroid
error and by ~0.03 per 10 % of LSF-width error (w3_oii_sim.py), so the blue-end calibration has to be known to
a few km/s and a few per cent. Lines (air -> vacuum): Hg I 3650.153 (+3654.836, 3662.880/3663.276 fitted jointly),
4046.563, 4358.328. Topocentric lines sit at lambda0 (1 - HELIO_RV/c) on the heliocentric pipeline grid (wave-2 sign
test). Per BOSS dark group: clipped-mean total sky (flux + sky), Gaussian fits with a local quadratic continuum and
(i) sigma = wdisp, free centroid; (ii) free sigma. [O I]5577 (airglow, isolated) is added as a mid-band LSF check.
Reads spec-lite FITS tables only (python -I).
Output: data/derived/w3_oii_sdss_wavecal.json
"""
import csv
import json
import sys
import numpy as np
from astropy.io import fits

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_oii_common import air2vac, C, gauss  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula/'
HG = {'Hg3650': [3650.153, 3654.836], 'Hg3663': [3662.880, 3663.276], 'Hg4047': [4046.563], 'Hg4358': [4358.328],
      'OI5577': [5577.339]}
REL = {'Hg3650': [1.0, 0.12], 'Hg3663': [0.25, 1.0], 'Hg4047': [1.0], 'Hg4358': [1.0], 'OI5577': [1.0]}   # blend weights
WIN = {'Hg3650': (3635.0, 3672.0), 'Hg3663': (3655.0, 3680.0), 'Hg4047': (4030.0, 4065.0), 'Hg4358': (4340.0, 4378.0),
       'OI5577': (5560.0, 5595.0)}


def fcol(rows, k):
    out = []
    for r in rows:
        try:
            out.append(float(r[k]))
        except (ValueError, TypeError):
            out.append(np.nan)
    return np.array(out)


def fit_line(lam, y, w, lines, rel, v, s, npoly=3):
    xc = (lam - lam.mean()) / 20.0
    g = sum(r * gauss(lam, l0 * (1 + v / C), s) for l0, r in zip(lines, rel))
    X = np.vstack([g] + [xc ** p for p in range(npoly)])
    A = (X * w) @ X.T
    p = np.linalg.solve(A, (X * w) @ y)
    res = y - p @ X
    return p[0], float(np.sum(w * res * res)), float(np.sqrt(np.linalg.inv(A)[0, 0]))


def main():
    Z = np.load(ROOT + 'data/derived/verify_sdss_windows.npz')
    files, ins = Z['m_file'], Z['m_instrument']
    F = list(csv.DictReader(open(ROOT + 'data/derived/verify_sdss_fibres.csv')))
    G = list(csv.DictReader(open(ROOT + 'data/derived/verify_sdss_groups.csv')))
    gid = np.array([int(r['gid']) for r in F])
    gins = np.array([r['instrument'] for r in G])
    sb, nf = fcol(G, 'sb5070_AB_cal'), fcol(G, 'n_fib')
    ok = (nf >= 6) & np.isfinite(sb)
    med = np.nanmedian(sb[ok & (gins == 'BOSS')])
    boss = np.where(ok & (gins == 'BOSS') & (sb > med - 0.75))[0]
    helio = fcol(G, 'helio_rv')
    out = {'lines_air': HG, 'groups': []}
    vg = np.arange(-60.0, 60.01, 0.5)
    for g in boss:
        sel = np.where(gid == g)[0][:40]
        specs, ivs, wds, lam0 = [], [], [], None
        for i in sel:
            with fits.open(ROOT + str(files[i]), memmap=False) as h:
                t = h[1].data
                lam = 10 ** np.asarray(t['loglam'], float)
                m = (lam > 3620) & (lam < 5600)
                lam = lam[m]
                tot = np.asarray(t['flux'], float)[m] + np.asarray(t['sky'], float)[m]
                iv = np.asarray(t['ivar'], float)[m]
                wd = np.asarray(t['wdisp'], float)[m]
            if lam0 is None:
                lam0 = lam
            if lam.size != lam0.size or np.max(np.abs(lam - lam0)) > 1e-3:
                continue
            specs.append(tot)
            ivs.append(iv)
            wds.append(wd)
        if len(specs) < 6:
            continue
        S = np.array(specs)
        IV = np.array(ivs)
        S[IV <= 0] = np.nan
        med_s = np.nanmedian(S, 0)
        mad = 1.4826 * np.nanmedian(np.abs(S - med_s), 0)
        S[np.abs(S - med_s) > 4 * np.maximum(mad, 1e-6)] = np.nan
        M = np.nanmean(S, 0)
        V = np.nanvar(S, 0) / np.maximum(np.sum(np.isfinite(S), 0), 1)
        sig = np.nanmedian(np.array(wds), 0) * lam0 * np.log(10) * 1e-4
        vtop = -helio[g]
        row = dict(gid=int(g), plate=G[g]['plate'], mjd=G[g]['mjd'], helio_rv=float(helio[g]), n=len(specs))
        for nm, lines_air in HG.items():
            lines = list(air2vac(np.array(lines_air)))
            lo, hi = WIN[nm]
            m = (lam0 > lo) & (lam0 < hi) & np.isfinite(M) & np.isfinite(V) & (V > 0)
            if m.sum() < 20:
                continue
            x, y, w = lam0[m], M[m], 1 / V[m]
            s0 = float(np.interp(lines[np.argmax(REL[nm])], lam0, sig))
            chi = [fit_line(x, y, w, lines, REL[nm], vtop + dv, s0)[1] for dv in vg]
            k = int(np.argmin(chi))
            kk = np.clip(k, 1, len(vg) - 2)
            c3 = np.polyfit(vg[kk - 1:kk + 2], np.array(chi)[kk - 1:kk + 2], 2)
            dv = -c3[1] / (2 * c3[0]) if c3[0] > 0 else vg[k]
            amp, _, eamp = fit_line(x, y, w, lines, REL[nm], vtop + dv, s0)
            # free sigma at that centroid
            sg = s0 * np.arange(0.7, 2.2, 0.01)
            chs = [fit_line(x, y, w, lines, REL[nm], vtop + dv, ss)[1] for ss in sg]
            row[nm] = dict(dv_kms=float(dv), amp=float(amp), snr=float(amp / eamp), sigma_wdisp=s0,
                           sigma_fit_over_wdisp=float(sg[int(np.argmin(chs))] / s0))
        out['groups'].append(row)
    summ = {}
    for nm in HG:
        rows = [r[nm] for r in out['groups'] if nm in r and r[nm]['snr'] > 5]
        if not rows:
            summ[nm] = None
            continue
        dv = np.array([r['dv_kms'] for r in rows])
        sr = np.array([r['sigma_fit_over_wdisp'] for r in rows])
        summ[nm] = dict(n=len(rows), dv_median=float(np.median(dv)), dv_err=float(1.2533 * np.std(dv) / np.sqrt(len(dv))),
                        dv_p16=float(np.percentile(dv, 16)), dv_p84=float(np.percentile(dv, 84)),
                        sigma_over_wdisp_median=float(np.median(sr)), sigma_over_wdisp_p16=float(np.percentile(sr, 16)),
                        sigma_over_wdisp_p84=float(np.percentile(sr, 84)),
                        snr_median=float(np.median([r['snr'] for r in rows])))
        print(nm, summ[nm], flush=True)
    out['summary'] = summ
    json.dump(out, open(ROOT + 'data/derived/w3_oii_sdss_wavecal.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
