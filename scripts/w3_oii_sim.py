"""w3_oii_sim.py -- (wave 3, tag oii) simulations of the [O II] 3729/3726 doublet fits used by waves 1-2.

Instruments (LSF from the data on disk):
  DESI : linear 0.8 A grid; sigma_R(3727) = median SPARCL wave_sigma (w3_oii_desi_lsf.json, 430 pilot coadds);
         wave-2 per-fibre fitter: v = +14 km/s, sigma^2 = sigma_R^2 + 0.8^2/12 + (12 km/s)^2, quadratic baseline.
  BOSS : log-lambda grid 1e-4 dex (0.858 A at 3727); sigma = median spec-lite wdisp of BOSS sky fibres at 3728 A;
         wave-2 differential fitter: v = +13 km/s, sigma = wdisp, linear baseline (+ continuum template).
Truth: two Gaussians at true ratio R_true (1.46 low-density limit, or 1.30), v_true = +11 km/s, sigma_int 14 km/s.
Systematics simulated (noiseless unless stated; all linear, so the bias is deterministic):
  LSF scale mismatch; fixed vs free width; wavelength-dependent LSF (sigma_3729 != sigma_3726); velocity offset
  (fit v - true v); two-strand velocity blending; non-Gaussian LSF (Gauss-Hermite h3 skew, h4 / broad wings);
  ON-OFF differencing with a WIM [O II] (low-density) in both ON and OFF but a different LSF width / wavelength
  zero point / throughput in the OFF; a solar Fe I 3727.6 (air) absorption between the lines from an ON-OFF continuum
  mismatch; DESI spectroperfectionism extraction of an off-grid line (1D specter replica: delta-function basis,
  resolution matrix from sqrt(iCov)); ratio-estimator noise bias at the observed S/N.
Output: data/derived/w3_oii_sim.json
"""
import json
import sys
import numpy as np
import scipy.linalg
from scipy.special import erf

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_oii_common import L, C  # noqa: E402

D = '/Users/mickey/astronomy/silk-nebula/data/derived/'
L26, L29 = L['OII3726'], L['OII3729']
FE3728 = 3728.68      # Fe I 3727.62 air -> vacuum (solar line between the [O II] components)
LSF = json.load(open(D + 'w3_oii_desi_lsf.json'))
SR_DESI = LSF['3727.1']['median']
try:
    SR_BOSS = json.load(open(D + 'w3_oii_sdss_lsf.json'))['boss_wdisp_sigma_3728_median']
except (FileNotFoundError, KeyError):
    SR_BOSS = 1.012
GRID = {'DESI': 3600.0 + 0.8 * np.arange(int((3800 - 3600) / 0.8)),
        'BOSS': 10 ** (3.5500 + 1e-4 * np.arange(0, 400))}
for k in GRID:
    g = GRID[k]
    GRID[k] = g[(g > 3700) & (g < 3760)]
FITWIN = (3712.0, 3745.0)
V_TRUE, SINT_KMS = 11.0, 14.0


def gh_profile(x, mu, s, h3=0.0, h4=0.0):
    """Unit-integral Gauss-Hermite profile (h3 skew, h4 kurtosis), normalised numerically on a fine grid."""
    y = (x - mu) / s
    H3 = (2 * np.sqrt(2) * y ** 3 - 3 * np.sqrt(2) * y) / np.sqrt(6)
    H4 = (4 * y ** 4 - 12 * y ** 2 + 3) / np.sqrt(24)
    g = np.exp(-0.5 * y ** 2) * (1 + h3 * H3 + h4 * H4)
    xf = np.linspace(mu - 12 * s, mu + 12 * s, 6001)
    yf = (xf - mu) / s
    gf = np.exp(-0.5 * yf ** 2) * (1 + h3 * (2 * np.sqrt(2) * yf ** 3 - 3 * np.sqrt(2) * yf) / np.sqrt(6)
                                   + h4 * (4 * yf ** 4 - 12 * yf ** 2 + 3) / np.sqrt(24))
    return g / np.trapezoid(gf, xf)


def G(x, mu, s):
    return np.exp(-0.5 * ((x - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def doublet(x, R, v, s26, s29=None, h3=0.0, h4=0.0, wing=0.0, wing_fac=2.0):
    s29 = s26 if s29 is None else s29
    out = 0
    for amp, l0, s in ((1.0, L26, s26), (R, L29, s29)):
        mu = l0 * (1 + v / C)
        prof = gh_profile(x, mu, s, h3, h4) if (h3 or h4) else G(x, mu, s)
        if wing:
            prof = (1 - wing) * prof + wing * G(x, mu, wing_fac * s)
        out = out + amp * prof
    return out


def fit(x, y, v, s, npoly=2, w=None, extra=None, s29=None):
    m = (x > FITWIN[0]) & (x < FITWIN[1])
    xx = x[m]
    s29 = s if s29 is None else s29
    cols = [G(xx, L26 * (1 + v / C), s), G(xx, L29 * (1 + v / C), s29)]
    cols += [((xx - 3728.0) / 30.0) ** p for p in range(npoly)]
    if extra is not None:
        cols += [e[m] for e in extra]
    X = np.vstack(cols)
    ww = np.ones(m.sum()) if w is None else w[m]
    p = np.linalg.solve((X * ww) @ X.T, (X * ww) @ y[m])
    r = y[m] - p @ X
    return p[1] / p[0], p, float(np.sum(ww * r * r))


def fit_free(x, y, s0, vgrid, sgrid, npoly=2):
    best = (np.inf, None)
    for v in vgrid:
        for s in sgrid:
            r, p, c2 = fit(x, y, v, s, npoly)
            if c2 < best[0]:
                best = (c2, (r, v, s))
    return best[1]


def specter_1d(R, v, s_ccd=0.60, ccd_dl=0.57, lo=3680.0, hi=3780.0):
    """1D replica of specter ex2d: delta-function wavelength basis on the 0.8 A grid, CCD pixel-integrated PSF,
    least-squares deconvolution, resolution matrix R = sqrt(iCov) row-normalised, flux = R xflux (per A)."""
    grid = np.arange(3690.0, 3770.01, 0.8)
    edges = np.arange(lo, hi + ccd_dl, ccd_dl)

    def pixint(mu, s):
        a = (edges[:-1] - mu) / (np.sqrt(2) * s)
        b = (edges[1:] - mu) / (np.sqrt(2) * s)
        return 0.5 * (erf(b) - erf(a))
    A = np.array([pixint(g, s_ccd) for g in grid]).T
    img = pixint(L26 * (1 + v / C), s_ccd) + R * pixint(L29 * (1 + v / C), s_ccd)
    icov = A.T @ A
    xflux = np.linalg.solve(icov, A.T @ img)
    w, vv = scipy.linalg.eigh(icov)
    sq = (vv * np.sqrt(np.clip(w, 0, None))) @ vv.T
    Rm = sq / sq.sum(1)[:, None]
    k = grid.size // 2
    d = np.arange(grid.size) - k
    sR = float(np.sqrt(np.sum(Rm[k] * d ** 2)) * 0.8)
    return grid, Rm @ xflux / 0.8, sR


def main():
    out = {'sigma_R': {'DESI': SR_DESI, 'BOSS': SR_BOSS}, 'v_true': V_TRUE, 'sigma_int_kms': SINT_KMS}
    for inst in ('DESI', 'BOSS'):
        x = GRID[inst]
        sR = SR_DESI if inst == 'DESI' else SR_BOSS
        st = np.sqrt(sR ** 2 + (3728.0 * SINT_KMS / C) ** 2)           # true total sigma
        s_w2 = np.sqrt(sR ** 2 + 0.8 ** 2 / 12 + (3728.0 * 12 / C) ** 2) if inst == 'DESI' else sR
        v_w2 = 14.0 if inst == 'DESI' else 13.0
        res = {'sigma_true_A': st, 'sigma_wave2_A': s_w2, 'v_wave2': v_w2, 'separation_over_sigma': (L29 - L26) / st}
        for Rt in (1.46, 1.30):
            r = {}
            y = doublet(x, Rt, V_TRUE, st)
            r['exact_model'] = fit(x, y, V_TRUE, st)[0]
            r['wave2_setup'] = fit(x, y, v_w2, s_w2)[0]
            for f in (0.90, 0.95, 1.05, 1.10, 1.20):        # true sigma = f x fitted sigma
                r[f'lsf_true_over_fit_{f:.2f}'] = fit(x, doublet(x, Rt, V_TRUE, f * st), V_TRUE, st)[0]
            for dv in (-10, -5, -3, 3, 5, 10):                # fitted v - true v
                r[f'v_fit_minus_true_{dv:+d}'] = fit(x, y, V_TRUE + dv, st)[0]
            for q in (0.95, 0.98, 1.02, 1.05):                # sigma(3729)/sigma(3726) true, fitted common
                r[f'sig29_over_sig26_{q:.2f}'] = fit(x, doublet(x, Rt, V_TRUE, st, q * st), V_TRUE, st)[0]
            for dvb in (5, 10, 20):                           # two strands at v +- dvb/2, equal [O II]
                yb = 0.5 * doublet(x, Rt, V_TRUE - dvb / 2, st) + 0.5 * doublet(x, Rt, V_TRUE + dvb / 2, st)
                r[f'blend_two_strands_dv{dvb}'] = fit(x, yb, V_TRUE, st)[0]
            for h3 in (-0.10, -0.05, 0.05, 0.10):
                r[f'lsf_h3_{h3:+.2f}'] = fit(x, doublet(x, Rt, V_TRUE, st, h3=h3), V_TRUE, st)[0]
            for h4 in (0.05, 0.10):
                r[f'lsf_h4_{h4:+.2f}'] = fit(x, doublet(x, Rt, V_TRUE, st, h4=h4), V_TRUE, st)[0]
            for wg in (0.05, 0.10, 0.20):
                r[f'lsf_wing{wg:.2f}_x2'] = fit(x, doublet(x, Rt, V_TRUE, st, wing=wg), V_TRUE, st)[0]
            # free width / free velocity on the exact data
            r['free_sigma'] = fit_free(x, y, st, [V_TRUE], st * np.arange(0.8, 1.31, 0.005))[0]
            r['free_v_sigma'] = fit_free(x, y, st, np.arange(-10, 35.1, 0.5), st * np.arange(0.8, 1.31, 0.01))[0]
            # ON-OFF: ON = filament (1) + WIM (W x filament [O II], ratio 1.46); OFF = WIM with a different LSF width,
            # wavelength offset or throughput (and a solar continuum dip from a continuum level mismatch)
            onoff = {}
            for W in (1.0, 3.0):
                for dsig in (-0.05, 0.05):
                    yy = y + W * (1 + Rt * 0) * 0
                    on = y + W * doublet(x, 1.46, 12.0, st) / 2.46 * (1 + Rt)
                    offs = W * doublet(x, 1.46, 12.0, st * (1 + dsig)) / 2.46 * (1 + Rt)
                    onoff[f'WIM{W:.0f}x_offLSF{dsig:+.2f}'] = fit(x, on - offs, V_TRUE, st)[0]
                for dl in (-0.05, 0.05):                            # OFF wavelength offset (A)
                    on = y + W * doublet(x, 1.46, 12.0, st) / 2.46 * (1 + Rt)
                    offs = W * doublet(x, 1.46, 12.0 + dl / 3728.0 * C, st) / 2.46 * (1 + Rt)
                    onoff[f'WIM{W:.0f}x_offdlam{dl:+.2f}A'] = fit(x, on - offs, V_TRUE, st)[0]
                for eps in (-0.05, 0.05):                           # OFF throughput mismatch
                    on = y + W * doublet(x, 1.46, 12.0, st) / 2.46 * (1 + Rt)
                    offs = (1 + eps) * W * doublet(x, 1.46, 12.0, st) / 2.46 * (1 + Rt)
                    onoff[f'WIM{W:.0f}x_offscale{eps:+.2f}'] = fit(x, on - offs, V_TRUE, st)[0]
            # solar Fe I 3728.7 dip of depth d (in units of the 3726 line peak) from a continuum mismatch
            pk = G(0, 0, st)
            for dd in (0.02, 0.05, 0.10):
                dip = -dd * pk * np.exp(-0.5 * ((x - FE3728) / st) ** 2)
                onoff[f'FeI3728_dip_{dd:.2f}_of_peak'] = fit(x, y + dip, V_TRUE, st)[0]
            r['onoff'] = onoff
            # noise: ratio-estimator bias and scatter at the observed per-pixel S/N (peak of 3726 / pixel noise)
            rng = np.random.default_rng(3)
            noise = {}
            for snr in (8.0, 12.0, 17.0, 25.0):
                sd = G(0, 0, st) / snr
                rr = []
                for _ in range(3000):
                    rr.append(fit(x, y + rng.normal(0, sd, x.size), V_TRUE, st)[0])
                rr = np.array(rr)
                noise[f'peakSN{snr:.0f}'] = dict(mean=float(np.mean(rr)), median=float(np.median(rr)),
                                                 p16=float(np.percentile(rr, 16)), p84=float(np.percentile(rr, 84)))
            r['noise'] = noise
            res[f'Rtrue_{Rt:.2f}'] = r
        # velocity / width offsets that turn 1.46 into 1.30 (wave-2 width, scan v)
        y = doublet(x, 1.46, V_TRUE, st)
        scan = {}
        for dv in np.arange(-20, 20.1, 1.0):
            scan[f'{dv:+.0f}'] = fit(x, y, V_TRUE + dv, st)[0]
        res['ratio_vs_v_offset_Rtrue1.46'] = scan
        sc2 = {}
        for f in np.arange(0.80, 1.31, 0.05):
            sc2[f'{f:.2f}'] = fit(x, doublet(x, 1.46, V_TRUE, f * st), V_TRUE, st)[0]
        res['ratio_vs_lsf_true_over_fit_Rtrue1.46'] = sc2
        out[inst] = res
        print(inst, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in res['Rtrue_1.46'].items()
                                if not isinstance(v, dict)}, indent=0))
        print(inst, 'onoff', {k: round(v, 3) for k, v in res['Rtrue_1.46']['onoff'].items()})
        print(inst, 'noise', res['Rtrue_1.46']['noise'])
        print(inst, 'v-scan', {k: round(v, 3) for k, v in list(res['ratio_vs_v_offset_Rtrue1.46'].items())[::4]})
    # DESI extraction replica: phase dependence
    sp = {}
    for s_ccd in (0.55, 0.62):
        for v in (0.0, 5.0, 10.6, 14.0, 17.0, 20.9, 25.0, 30.0):
            grid, F, sR = specter_1d(1.46, v, s_ccd)
            ph26 = ((L26 * (1 + v / C) - 3600) / 0.8) % 1
            ph29 = ((L29 * (1 + v / C) - 3600) / 0.8) % 1
            r1 = fit(grid, F, v, sR)[0]
            r2 = fit(grid, F, v, np.sqrt(sR ** 2 + 0.8 ** 2 / 12))[0]
            sp[f's{s_ccd:.2f}_v{v:+.1f}'] = dict(phase3726=float(ph26), phase3729=float(ph29), sigma_R=sR,
                                                  ratio_fit_sigmaR=float(r1), ratio_fit_sigmaR_plus_pix=float(r2))
    out['DESI_specter_replica_Rtrue1.46'] = sp
    rs = [d['ratio_fit_sigmaR'] for d in sp.values()]
    print('specter replica ratio range', min(rs), max(rs))
    json.dump(out, open(D + 'w3_oii_sim.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
