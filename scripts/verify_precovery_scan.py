"""verify_precovery step 6: RA-shift scan of the combined core amplitude (template and window move together).

For each shift s (deg of RA at fixed Dec) the unit E+W ridge template (FWHM 0.28 deg) and the +-2.5 deg window are both
centred on the shifted spine; no other masking. Per-run weights are fixed to the s = 0 'off' errors (verify_precovery_fit
per-run e = max(formal, off-window null std)). A real ridge appears as a peak at s = 0 with the template's
nuisance-filtered autocorrelation shape; injected copies of the template (x the observed s=0 amplitude) give the
expected shape.  Bands/estimators: g, u (psp, framessub), r, i (framessub).
Output: data/derived/verify_precovery_scan.json
"""
import sys, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
import verify_precovery_fit as VF
from verify_precovery_lib import DER

Z = VF.Z; SH = VF.SH
F = json.loads((DER / 'verify_precovery_fit.json').read_text())


def eff(r):
    eo = r.get('e_null_off'); en = r.get('e_null_in') or 0.0
    return max(r['eA_formal'], eo) if eo else max(r['eA_formal'], en)


def select_shift(run, s):
    m = (Z['run'] == run) & VF.base_mask(2)
    if s == 0:
        d = Z['dE']; o = Z['overE']
    else:
        i = int(np.where(np.isclose(SH, s))[0][0]); d = Z['dE_shift'][i]; o = Z['overE_shift'][i]
    idx = np.where(m & (np.abs(d) <= 2.5) & (o == 0))[0]
    dd = d[idx]
    ok = (dd < -0.7).sum() >= 6 and (dd > 0.7).sum() >= 6 and (np.abs(dd) < 0.3).sum() >= 4
    return idx if ok and len(np.unique(Z['camcol'][idx])) >= 3 else None


def scan(case, band, est, conv):
    runs = {r['run']: r for r in F['cases'][case]['runs'] if 12 <= r['dec_cross'] < 20.5}
    y = Z[VF.est_key(band, est)].astype(float)
    shifts = [0.0] + [float(s) for s in SH if -3.01 <= s <= 6.01]
    shifts = sorted(shifts)
    out = []
    A0 = np.sum([runs[k]['A'] / eff(runs[k]) ** 2 for k in runs]) / np.sum([1 / eff(runs[k]) ** 2 for k in runs])
    for s in shifts:
        num = den = 0.0; inj_num = 0.0; n = 0
        for run, r in runs.items():
            idx = select_shift(run, s)
            if idx is None:
                continue
            Ts = Z['T_gauss028_mean'][idx] if s == 0 else Z['TS_gauss028'][int(np.where(np.isclose(SH, s))[0][0])][idx]
            if np.ptp(Ts) < 0.3:
                continue
            fr = VF.fit(idx, y, [Ts])
            # expected response at this shift to the REAL ridge (template at s=0) of amplitude A0
            T0 = Z['T_gauss028_mean'][idx]
            yy = np.zeros(len(y)); yy[idx] = A0 * conv * T0
            fi = VF.fit(idx, yy, [Ts], clip=99)
            w = 1 / eff(r) ** 2
            num += w * fr['coef'][0] / conv; inj_num += w * fi['coef'][0] / conv; den += w; n += 1
        if n >= 5:
            out.append(dict(shift=s, A=num / den, A_expected_from_ridge=inj_num / den, n=n, err=1 / np.sqrt(den)))
    return dict(case=case, A0=float(A0), scan=out)


def main():
    res = {}
    for case, band, est, conv in [('g_psp_gauss028_OIIIunits', 'g', 'psp', VF.CONV_OIII_ONLY),
                                  ('g_framessub_gauss028_OIIIunits', 'g', 'framessub', VF.CONV_OIII_ONLY),
                                  ('u_psp_gauss028_OIIunits', 'u', 'psp', VF.CONV_OII_U),
                                  ('u_framessub_gauss028_OIIunits', 'u', 'framessub', VF.CONV_OII_U),
                                  ('r_framessub_gauss028', 'r', 'framessub', VF.CONV['r']),
                                  ('i_framessub_gauss028', 'i', 'framessub', VF.CONV['i'])]:
        res[case] = scan(case, band, est, conv)
        sc = res[case]['scan']
        off = np.array([x['A'] for x in sc if abs(x['shift']) >= 2.0])
        z = [x for x in sc if x['shift'] == 0][0]
        print(f"{case:34s} A(0)={z['A']:+.2f}  off(|s|>=2) mean {off.mean():+.2f} std {off.std():.2f}  -> {z['A'] / off.std():.1f} sigma; n_off={len(off)}")
        res[case]['offpeak_mean'] = float(off.mean()); res[case]['offpeak_std'] = float(off.std())
    (DER / 'verify_precovery_scan.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
