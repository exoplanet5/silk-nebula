"""verify_sdss: wavelength zero point near [O III]5007 from topocentric sky lines in each group's mean total sky.

Lines: Hg I 4358.3 (street lamps; vac 4359.56), [N I] 5198/5200 blend (vac 5199.35/5201.71), [O I] 5577 (vac 5578.89),
[O I] 6300 (red camera).  For each, the offset dv = v_measured - v_expected, with v_expected = -HELIO_RV (topocentric
emission on the pipeline heliocentric grid; sign checked against astropy in verify_sdss_groups).  The zero point at
5008 A is interpolated linearly in wavelength between Hg4358 and [O I]5577 where both are measured.
Writes data/derived/verify_sdss_zeropoints.csv
"""
import sys, csv
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, VAC, C, load_windows, flam_to_R_per_A, fit_lines
from verify_sdss_groups import clipped_mean, sun_on


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def scan(lm, y, w, lines, sig, npoly, ties=None, extra=None, vr=(-80, 80)):
    vs = np.arange(vr[0], vr[1] + 0.01, 1.0)
    ch = np.array([fit_lines(lm, y, w, lines, v, sig, npoly=npoly, ties=ties, extra=extra)['_chi2'] for v in vs])
    j = int(np.argmin(ch))
    if j < 3 or j > len(vs) - 4:
        return np.nan, np.nan, np.nan
    v2 = np.arange(vs[j] - 3, vs[j] + 3.01, 0.2)
    c2 = np.array([fit_lines(lm, y, w, lines, v, sig, npoly=npoly, ties=ties, extra=extra)['_chi2'] for v in v2])
    pc = np.polyfit(v2, c2, 2)
    if pc[0] <= 0:
        return np.nan, np.nan, np.nan
    vb = -pc[1] / (2 * pc[0])
    f = fit_lines(lm, y, w, lines, vb, sig, npoly=npoly, ties=ties, extra=extra)
    chi2nu = f['_chi2'] / max(f['_dof'], 1)
    free = [l for l in lines if not (ties and l in ties)]
    return float(vb), float(np.sqrt(2 / pc[0]) * np.sqrt(max(chi2nu, 1))), f[free[0]][0]


def main():
    D = load_windows()
    lam = D['lam']
    F = read_csv(ROOT / 'data/derived/verify_sdss_fibres.csv')
    G = read_csv(ROOT / 'data/derived/verify_sdss_groups.csv')
    gid = np.array([int(r['gid']) for r in F])
    rows = []
    for r in G:
        g = int(r['gid'])
        sel = np.where(gid == g)[0]
        diam = float(r['diam'])
        X = flam_to_R_per_A(D['tot'][sel].astype(float), lam, diam)
        iv = D['ivar'][sel].astype(float)
        V = np.where(iv > 0, flam_to_R_per_A(1.0, lam, diam) ** 2 / np.where(iv > 0, iv, 1), np.inf)
        M, MV = clipped_mean(X, V)[:2]
        sig = np.nanmedian(D['wdisp'][sel], 0) * lam * np.log(10) * 1e-4
        vexp = -float(r['helio_rv'])
        o = dict(gid=g, plate=r['plate'], mjd=r['mjd'], sp=r['sp'], instrument=r['instrument'], helio_rv=float(r['helio_rv']))
        for nm, lo, hi, lines, ties, usesun in (
                ('Hg4358', 4345, 4375, ['Hg4358', 'OIII4363'], None, True),
                ('NI', 5185, 5216, ['NI5198', 'NI5200'], {'NI5198': ('NI5200', 1.6)}, True),
                ('OI5577', 5560, 5598, ['OI5577'], None, False),
                ('OI6300', 6288, 6316, ['OI6300'], None, False),
                ('OIII5007abs', 4990, 5025, ['OIII5007'], None, True)):
            m = (lam > lo) & (lam < hi) & np.isfinite(M) & np.isfinite(MV)
            if m.sum() < 15:
                o[f'dv_{nm}'] = np.nan; o[f'edv_{nm}'] = np.nan; o[f'I_{nm}'] = np.nan
                continue
            s_ = float(np.nanmedian(sig[m]))
            ex = {'sun': sun_on(lam[m], sig[m])} if usesun else None
            try:
                vb, ev, amp = scan(lam[m], M[m], 1 / MV[m], lines, s_, 2, ties=ties, extra=ex)
            except Exception:
                vb, ev, amp = np.nan, np.nan, np.nan
            o[f'dv_{nm}'] = vb - vexp if nm != 'OIII5007abs' else vb       # [O III] absolute: helio velocity of total sky line
            o[f'edv_{nm}'] = ev; o[f'I_{nm}'] = amp
        # [N I] ratio sensitivity: ratio 5198/5200 = 1.3 and 2.0
        m = (lam > 5185) & (lam < 5216) & np.isfinite(M) & np.isfinite(MV)
        ex = {'sun': sun_on(lam[m], sig[m])} if m.sum() >= 15 else None
        for rr in (1.3, 2.0):
            if ex is None:
                o[f'dv_NI_r{rr}'] = np.nan
                continue
            vb, ev, amp = scan(lam[m], M[m], 1 / MV[m], ['NI5198', 'NI5200'], float(np.nanmedian(sig[m])), 2,
                               ties={'NI5198': ('NI5200', rr)}, extra=ex)
            o[f'dv_NI_r{rr}'] = vb - vexp
        # interpolate zero point at 5008.24 A between Hg4358 and OI5577
        if np.isfinite(o['dv_Hg4358']) and np.isfinite(o['dv_OI5577']):
            t = (5008.24 - VAC['Hg4358']) / (VAC['OI5577'] - VAC['Hg4358'])
            o['zp5007_interp'] = (1 - t) * o['dv_Hg4358'] + t * o['dv_OI5577']
        else:
            o['zp5007_interp'] = np.nan
        rows.append(o)
        if g % 40 == 0:
            print(g, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in o.items()}, flush=True)
    keys = list(rows[0].keys())
    with open(ROOT / 'data/derived/verify_sdss_zeropoints.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for o in rows:
            w.writerow(o)
    print('done', len(rows))


if __name__ == '__main__':
    main()
