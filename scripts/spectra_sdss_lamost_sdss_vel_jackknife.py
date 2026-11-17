"""Robustness of the SDSS strand [O III] velocity and Hbeta limit (topic spectra_sdss_lamost):
instrument split (SDSS-I/II 3" fibres vs BOSS 2" fibres) and delete-one-group jackknife of the
south (Dec < 20) ON-minus-OFF stack built by spectra_sdss_lamost_sdss_substacks.py.

Output: data/derived/spectra_sdss_lamost_sdss_vel_jackknife.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import spectra_sdss_lamost_sdss_substacks as SS
import spectra_sdss_lamost_trace as TR

ROOT = Path('/Users/mickey/astronomy/silk-nebula')


def run(d, lam, groups, keys_arr, onm, offm, fpsf, zb, wd):
    S = SS.stack(d, lam, groups, keys_arr, onm, offm, fpsf)
    if S is None:
        return None
    w_on = dict(S['used'])
    zp = np.nansum([zb.get(k, np.nan) * n for k, n in w_on.items() if np.isfinite(zb.get(k, np.nan))]) / max(
        sum(n for k, n in w_on.items() if np.isfinite(zb.get(k, np.nan))), 1)
    out = dict(n_on=S['n_on'], n_groups=S['n_groups'], zp_blue=float(zp))
    for nm, lines, lo, hi, cam in SS.SETS[:2]:
        sig_pix = np.nanmedian(wd[(lam > lo) & (lam < hi)]) * 1e-4 * np.log(10)
        res, _ = SS.fit_set(lam, S['stack'], S['err'], lines, lo, hi, sig_pix, free_v=(nm == 'OIII'))
        if nm == 'OIII':
            out['OIII5007'] = res['OIII5007']['I']; out['OIII5007_err'] = res['OIII5007']['err']
            out['v_helio'] = res.get('v_raw', np.nan) - zp; out['v_err'] = res.get('v_raw_err', np.nan)
        else:
            out['Hb'] = res['Hb']['I']; out['Hb_err'] = res['Hb']['err']
    return out, [k for k, _ in S['used']]


def main():
    d, G = SS.load()
    lam = 10 ** d['loglam']
    keys_arr = np.array([f'{p}-{m}-{s}' for p, m, s in zip(d['plate'], d['mjd'], d['spectrograph'])])
    gkey = lambda r: f"{int(r['plate'])}-{int(r['mjd'])}-{int(r['sp'])}"
    bad = {gkey(r) for r in G if (r['sb5040_ABmag'] < 20.4) or (r['I_OIII_cal_err'] > 0.6)}
    fpsf = {gkey(r): float(r['f_psf']) for r in G}
    inst = {gkey(r): str(r['instrument']) for r in G}
    good = [gkey(r) for r in G if gkey(r) not in bad]
    Z = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_zeropoints.csv')
    zb = {str(r['group']): float(r['dv_OI5577']) for r in Z}
    wd = np.nanmedian(np.where(d['wdisp'] > 0, d['wdisp'], np.nan), 0)
    geo = TR.load()
    dist, inside, _ = TR.nearest_strand_distance(d['ra'], d['dec'], geo)
    ad = np.abs(dist)
    onm = inside & (ad < SS.ON) & (d['dec'] < 20)
    offm = inside & (ad > SS.OFF[0]) & (ad < SS.OFF[1]) & (d['dec'] < 20)
    res = {}
    full, used = run(d, lam, good, keys_arr, onm, offm, fpsf, zb, wd)
    res['south_full'] = full
    for name in ('SDSS', 'BOSS'):
        r = run(d, lam, [g for g in good if inst[g] == name], keys_arr, onm, offm, fpsf, zb, wd)
        res[f'south_{name}_only'] = r[0] if r else None
    jk = []
    for g in used:
        r = run(d, lam, [x for x in good if x != g], keys_arr, onm, offm, fpsf, zb, wd)
        if r:
            jk.append(dict(dropped=g, **r[0]))
    v = np.array([j['v_helio'] for j in jk])
    n = len(v)
    res['jackknife'] = dict(n=n, v_mean=float(v.mean()), v_jk_err=float(np.sqrt((n - 1) / n * np.sum((v - v.mean()) ** 2))),
                            v_min=float(v.min()), v_max=float(v.max()),
                            hb_max=float(max(j['Hb'] for j in jk)), oiii_min=float(min(j['OIII5007'] for j in jk)),
                            dropped_extremes=[jk[int(np.argmin(v))]['dropped'], jk[int(np.argmax(v))]['dropped']])
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_vel_jackknife.json').write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
