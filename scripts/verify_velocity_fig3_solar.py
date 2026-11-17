"""verify_velocity: zero point and continuum analysis of the paper's Fig. 3a blue-arm spectrum (vector-extracted
by verify_velocity_fig3_vectors.py).

1. Sky-line zero point of the plotted wavelength scale from all UVES-atlas sky lines that are detected
   (Hanuschik 2003, table 6, air wavelengths; [N I] also vs NIST levels).
2. Solar (Fraunhofer) content of the continuum: the continuum-normalised spectrum is regressed on the
   continuum-normalised Kurucz 2005 solar flux atlas (air, data/raw/verify_velocity_solar/) convolved to the
   measured sky-line LSF, on a grid of Doppler shifts.  The regression slope f_sun is the fraction of the
   continuum that carries solar-type absorption lines (zodiacal light, scattered sunlight, G/K starlight).
3. Refit of [O III] 5007 and 4959 with and without that solar component, to see how much a plain
   Gaussian + polynomial centroid is pulled by the Fraunhofer lines (Fe I 5005.7/5006.1, 5007.2 blend;
   Fe I 4957.3/4957.6).

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_fig3_solar.py
Output: data/derived/verify_velocity_fig3_solar.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter, gaussian_filter1d
from scipy.optimize import curve_fit

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
C = 299792.458
VCORR = -8.32
OUT = ROOT / 'data/derived/verify_velocity_fig3_solar.json'


def n_air(lv):
    s2 = (1e4 / np.asarray(lv, float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


NI_NIST_AIR = np.array([1e8 / 19233.177, 1e8 / 19224.464]) / n_air(np.array([1e8 / 19233.177, 1e8 / 19224.464]))
O3_AIR = {'5007': 5006.843, '4959': 4958.911}


def load():
    z = np.load(ROOT / 'data/derived/verify_velocity_fig3_spectra.npz')
    w, f = z['a_wave'], z['a_flux']
    o = np.argsort(w)
    return w[o], f[o]


def solar():
    d = np.loadtxt(ROOT / 'data/raw/verify_velocity_solar/kurucz2005_flux_494-526nm.asc')
    return d[:, 0] * 10.0, d[:, 1]


def g1(x, A, mu, s, c0, c1, c2=0.0):
    t = x - 5100.0
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2) + c0 + c1 * t + c2 * t * t


def fit_line(w, f, lo, hi, mu0, s0=0.3, quad=False):
    m = (w > lo) & (w < hi)
    x, y = w[m], f[m]
    if quad:
        fn = lambda xx, A, mu, s, c0, c1, c2: g1(xx, A, mu, s, c0, c1, c2)
        p0 = [y.max() - np.median(y), mu0, s0, np.median(y), 0, 0]
    else:
        fn = lambda xx, A, mu, s, c0, c1: g1(xx, A, mu, s, c0, c1)
        p0 = [y.max() - np.median(y), mu0, s0, np.median(y), 0]
    p, cov = curve_fit(fn, x, y, p0=p0, maxfev=20000)
    return p, np.sqrt(np.diag(cov)), float(np.std(y - fn(x, *p)))


def main():
    w, f = load()
    out = {}
    # ------------------------------------------------------------------ 1. sky lines
    sky = {'NI5198 (Hanuschik)': 5197.9282, 'NI5200 (Hanuschik)': 5200.2856,
           'NI5198 (NIST)': NI_NIST_AIR[0], 'NI5200 (NIST)': NI_NIST_AIR[1],
           's5202.98': 5202.9790, 's5217.39': 5217.3906, 's5224.15': 5224.1450, 's5232.15': 5232.1528,
           's5238.75': 5238.7520, 's5276.14': 5276.1431}
    skyfit = {}
    for k, l0 in sky.items():
        try:
            p, e, rms = fit_line(w, f, l0 - 1.0, l0 + 1.0, l0, 0.28)
            if p[0] <= 0 or abs(p[1] - l0) > 0.4 or not (0.15 < abs(p[2]) < 0.6):
                skyfit[k] = dict(rejected=True, mu=float(p[1]), A=float(p[0]), sig=float(abs(p[2])))
                continue
            skyfit[k] = dict(mu=float(p[1]), e_mu=float(e[1]), A=float(p[0]), sig=float(abs(p[2])),
                             snr=float(p[0] / rms), v_plot=float(C * (p[1] / l0 - 1)), e_v=float(C * e[1] / l0))
        except Exception as ex:
            skyfit[k] = dict(error=str(ex))
    out['sky_lines_on_plot_scale'] = skyfit
    others = [v for k, v in skyfit.items() if k.startswith('s') and 'v_plot' in v]
    if others:
        vv = np.array([o['v_plot'] for o in others])
        ee = np.array([max(o['e_v'], 1.0) for o in others])
        wgt = 1 / ee ** 2
        out['weak_sky_lines_zero_point'] = dict(n=len(vv), wmean=float((vv * wgt).sum() / wgt.sum()),
                                                err=float(1 / np.sqrt(wgt.sum())), median=float(np.median(vv)),
                                                values=vv.round(2).tolist())

    # LSF from [N I] 5198
    sig_lsf = skyfit['NI5198 (Hanuschik)']['sig']
    out['lsf_sigma_A'] = sig_lsf
    out['lsf_R'] = float(5198 / (2.3548 * sig_lsf))

    # ------------------------------------------------------------------ 2. solar content
    ws, fs = solar()
    dws = ws[1] - ws[0]
    fs_c = gaussian_filter1d(fs, sig_lsf / dws)

    def norm(x, y, width_A=12.0):
        step = np.median(np.diff(x))
        k = int(width_A / step) | 1
        cont = median_filter(y, size=k, mode='nearest')
        # upper-envelope style: iterate once with clipping of absorption
        return y / cont

    mask = np.ones_like(w, bool)
    for l0 in (4958.9, 5006.8, 5197.9, 5200.3, 5203.0, 5217.4, 5224.1, 5232.2, 5238.8, 5256.1, 5276.2):
        mask &= np.abs(w - l0 * (1 + 25 / C)) > 1.6
    fit_reg = (w > 4945) & (w < 5360) & mask
    yn = norm(w, f)
    results = []
    for v in np.arange(-60, 61, 2.0):
        sm = np.interp(w, ws * (1 + v / C), fs_c)
        sn = norm(w, sm)
        x = sn[fit_reg] - 1
        y = yn[fit_reg] - 1
        fsun = float((x * y).sum() / (x * x).sum())
        r = y - fsun * x
        results.append((v, fsun, float((r * r).sum()), float(np.corrcoef(x, y)[0, 1])))
    results = np.array(results)
    ib = np.argmin(results[:, 2])
    vbest, fbest = results[ib, 0], results[ib, 1]
    # bootstrap in 10 A blocks for f_sun uncertainty at vbest
    sm = np.interp(w, ws * (1 + vbest / C), fs_c)
    sn = norm(w, sm)
    x_all, y_all = sn - 1, yn - 1
    blocks = np.floor((w - 4940) / 10).astype(int)
    ub = np.unique(blocks[fit_reg])
    rng = np.random.default_rng(5)
    fb = []
    for _ in range(500):
        pick = rng.choice(ub, ub.size)
        idx = np.concatenate([np.where(fit_reg & (blocks == b))[0] for b in pick])
        fb.append((x_all[idx] * y_all[idx]).sum() / (x_all[idx] ** 2).sum())
    out['solar_regression'] = dict(v_best_plot_scale=float(vbest), f_sun=float(fbest), f_sun_boot_std=float(np.std(fb)),
                                   corr_at_best=float(results[ib, 3]),
                                   f_sun_vs_v=[(float(a), round(float(b), 3)) for a, b, *_ in results[::5]])
    # Mg b check (5167.3, 5172.7, 5183.6 air)
    mg = {}
    for l0 in (5167.32, 5172.68, 5183.60):
        m = np.abs(w - l0 * (1 + vbest / C)) < 0.35
        ms = np.abs(w - l0 * (1 + vbest / C)) < 0.35
        mg[str(l0)] = dict(data_depth=float(1 - np.median(yn[m])), solar_depth_conv=float(1 - np.median(sn[ms])))
    out['MgB_depths'] = mg

    # ------------------------------------------------------------------ 3. [O III] refits
    cont_sun = sm  # solar spectrum (unnormalised residual flux) at LSF and best shift
    ref = {}
    for line, lo, hi, mu0 in (('5007', 5002.0, 5013.0, 5007.25), ('4959', 4954.0, 4964.0, 4959.4)):
        m = (w > lo) & (w < hi)
        x, y = w[m], f[m]
        S = cont_sun[m]
        # (a) plain Gaussian + quadratic (as the paper)
        pa, ea, _ = fit_line(w, f, lo, hi, mu0, 0.45, quad=True)
        # (b) Gaussian + quadratic * [(1-f) + f*S]  with f = f_sun (fixed) and +-1 sigma
        rows = []
        for fsun in (fbest - np.std(fb), fbest, fbest + np.std(fb), 0.5, 1.0):
            fsun = float(np.clip(fsun, 0, 1))

            def model(xx, A, mu, s, c0, c1, c2):
                t = xx - 5100.0
                Sx = np.interp(xx, x, S)
                return A * np.exp(-0.5 * ((xx - mu) / s) ** 2) + (c0 + c1 * t + c2 * t * t) * ((1 - fsun) + fsun * Sx)
            p, cov = curve_fit(model, x, y, p0=[pa[0], pa[1], abs(pa[2]), pa[3], pa[4], pa[5]], maxfev=40000)
            rows.append(dict(f_sun=fsun, mu=float(p[1]), e_mu=float(np.sqrt(cov[1, 1])), sig=float(abs(p[2])),
                             shift_kms=float(C * (p[1] - pa[1]) / pa[1]), A=float(p[0])))
        ref[line] = dict(plain_mu=float(pa[1]), plain_e_mu=float(ea[1]), plain_sig=float(abs(pa[2])), with_solar=rows)
    out['OIII_refits'] = ref

    # ------------------------------------------------------------------ 4. calibration-free velocities
    vel = {}
    for line in ('5007', '4959'):
        for tag, mu in (('plain', ref[line]['plain_mu']), ('solar_fsun', ref[line]['with_solar'][1]['mu'])):
            vplot = C * (mu / O3_AIR[line] - 1)
            for sk in ('NI5198 (NIST)', 'NI5198 (Hanuschik)', 'NI5200 (NIST)', 'NI5200 (Hanuschik)'):
                if 'v_plot' in skyfit.get(sk, {}):
                    vel[f'{line}_{tag}_vs_{sk}'] = float(vplot - skyfit[sk]['v_plot'] + VCORR)
            if 'weak_sky_lines_zero_point' in out:
                vel[f'{line}_{tag}_vs_weaksky'] = float(vplot - out['weak_sky_lines_zero_point']['wmean'] + VCORR)
            vel[f'{line}_{tag}_plot_scale'] = float(vplot)
    out['V_helio_calibration_free'] = vel
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:6000])


if __name__ == '__main__':
    main()
