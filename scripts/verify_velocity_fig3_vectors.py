"""verify_velocity: recover the paper's own MRS-N blue-arm spectra from the vector graphics of Fig. 3
(arXiv:2610.13329v1 PDF, page 7) and measure [O III] relative to the simultaneously recorded [N I] 5198/5200
airglow doublet.  This is calibration-free: the airglow is at topocentric rest, so
    V_topo([O III]) = c * [ (lam_O / lam_N) * (lam_N0 / lam_O0) - 1 ]
and V_helio = V_topo + v_corr (v_corr = -8.32 km/s for this plate; verify_velocity_lamost_arms.json).

Axis calibration uses the matplotlib tick marks drawn as vector lines (listed in the output JSON).
The plotted curves are the paper's data samples (polyline vertices).  The PDF is the arXiv v1 file already on
disk (sources/arxiv/2610.13329v1.pdf); it is only parsed, never executed.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_fig3_vectors.py
Output: data/derived/verify_velocity_fig3_spectra.npz, data/derived/verify_velocity_fig3_measure.json
"""
import json
from pathlib import Path

import fitz
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
PDF = ROOT / 'sources/arxiv/2610.13329v1.pdf'
OUTN = ROOT / 'data/derived/verify_velocity_fig3_spectra.npz'
OUTJ = ROOT / 'data/derived/verify_velocity_fig3_measure.json'
C = 299792.458
VCORR = -8.32


def n_air(lv):
    s2 = (1e4 / np.asarray(lv, float)) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


def vac2air(lv):
    return lv / n_air(lv)


NI_VAC = np.array([1e8 / 19233.177, 1e8 / 19224.464])
NI_AIR = vac2air(NI_VAC)
O3_AIR = np.array([4958.911, 5006.843])


def poly_pts(d):
    pts = []
    for it in d['items']:
        if it[0] == 'l':
            if not pts:
                pts.append((it[1].x, it[1].y))
            pts.append((it[2].x, it[2].y))
        elif it[0] == 'c':
            if not pts:
                pts.append((it[1].x, it[1].y))
            pts.append((it[4].x, it[4].y))
    return np.array(pts)


def lin_cal(pix, val):
    A = np.vstack([pix, np.ones_like(pix)]).T
    (a, b), res, *_ = np.linalg.lstsq(A, val, rcond=None)
    rms = np.sqrt(np.mean((a * pix + b - val) ** 2))
    return a, b, rms


def g1(x, A, mu, s, c0, c1):
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2) + c0 + c1 * (x - mu)


def g2(x, A1, A2, mu1, s, c0, c1, sep):
    return (A1 * np.exp(-0.5 * ((x - mu1) / s) ** 2) + A2 * np.exp(-0.5 * ((x - mu1 - sep) / s) ** 2)
            + c0 + c1 * (x - mu1))


def fit1(w, f, lo, hi, mu0, s0=0.35):
    m = (w > lo) & (w < hi)
    x, y = w[m], f[m]
    p0 = [y.max() - np.median(y), mu0, s0, np.median(y), 0]
    p, cov = curve_fit(g1, x, y, p0=p0, maxfev=20000)
    r = y - g1(x, *p)
    # rescale covariance by residual scatter (curve_fit default absolute_sigma=False already does)
    e = np.sqrt(np.diag(cov))
    return dict(A=p[0], mu=p[1], sig=abs(p[2]), cont=p[3], e_mu=e[1], e_sig=e[2], rms=float(r.std()), npts=int(m.sum()))


def main():
    doc = fitz.open(PDF)
    page = doc[6]
    dr = page.get_drawings()
    # ---- axis calibrations from tick marks (PDF points); values read from the rendered tick labels
    xt_a = np.array([89.67, 144.41, 199.15, 253.9, 308.64, 363.38, 418.12, 472.87, 527.61])
    xv_a = np.arange(4950, 5351, 50.0)
    yt_a = np.array([209.44, 191.17, 172.9, 154.63, 136.36, 118.1, 99.83, 81.56])
    yv_a = np.arange(0, 141, 20.0)
    xt_b = np.array([92.99, 124.32, 155.65, 186.99, 218.32, 249.65, 280.98])
    xv_b = np.arange(4960, 5021, 10.0)
    yt_b = np.array([360.16, 338.88, 317.61, 296.33, 275.06])
    yv_b = np.arange(0, 201, 50.0)
    xt_cr = np.array([464.56, 477.32, 490.08, 502.84, 515.6, 528.36, 541.12, 553.88])
    xv_cr = np.arange(5000, 5015, 2.0)
    xt_cl = np.array([343.96, 368.77, 393.59, 418.4])
    xv_cl = np.arange(4950, 4966, 5.0)
    cal = {}
    for name, (xp, xv) in dict(a=(xt_a, xv_a), b=(xt_b, xv_b), c_right=(xt_cr, xv_cr), c_left=(xt_cl, xv_cl)).items():
        a, b, rms = lin_cal(xp, xv)
        cal[name] = dict(a=a, b=b, rms_A=rms)
    ya = lin_cal(yt_a, yv_a)
    yb = lin_cal(yt_b, yv_b)

    big = sorted(dr, key=lambda d: -len(d['items']))
    spec = {}
    # panel a: black, ~3589 segments, within x 71..554, y 72..209
    for d in big[:12]:
        r = d['rect']
        col = d.get('color')
        n = len(d['items'])
        if n > 3000 and r.y1 < 215:
            p = poly_pts(d)
            spec['a'] = (cal['a']['a'] * p[:, 0] + cal['a']['b'], ya[0] * p[:, 1] + ya[1])
        elif 590 < n < 700 and 295 < r.y0 < 340:
            key = {(0.0, 0.0, 0.0): 'b_P4', (0.733, 0.345, 0.424): 'b_P2', (0.541, 0.247, 0.502): 'b_P1',
                   (0.957, 0.725, 0.176): 'b_P3'}[tuple(round(c, 3) for c in col)]
            p = poly_pts(d)
            spec[key] = (cal['b']['a'] * p[:, 0] + cal['b']['b'], yb[0] * p[:, 1] + yb[1])
    # vertical marker / fit-centre lines
    marks = {}
    for d in dr:
        if len(d['items']) == 1 and d['items'][0][0] == 'l':
            p1, p2 = d['items'][0][1], d['items'][0][2]
            if abs(p1.x - p2.x) < 0.01 and abs(p1.y - p2.y) > 80 and d.get('color') not in (None, (0.0, 0.0, 0.0)):
                x = p1.x
                if 71 < x < 554 and p1.y < 215:
                    marks.setdefault('a', []).append(cal['a']['a'] * x + cal['a']['b'])
                elif 71 < x < 281:
                    marks.setdefault('b', []).append(cal['b']['a'] * x + cal['b']['b'])
                elif 343 < x < 434:
                    marks.setdefault('c_left_fitcentre', []).append(cal['c_left']['a'] * x + cal['c_left']['b'])
                elif 464 < x < 554:
                    marks.setdefault('c_right_fitcentre', []).append(cal['c_right']['a'] * x + cal['c_right']['b'])

    out = dict(calibration=cal, y_cal_a=ya, y_cal_b=yb, markers_A=marks)
    w, f = spec['a']
    o = np.argsort(w)
    w, f = w[o], f[o]
    out['panel_a_npts'] = int(w.size)
    out['panel_a_range'] = [float(w.min()), float(w.max())]
    out['panel_a_median_step_A'] = float(np.median(np.diff(w)))
    res = {}
    res['OIII5007'] = fit1(w, f, 5002.5, 5012.0, 5007.3)
    res['OIII4959'] = fit1(w, f, 4955.5, 4963.5, 4959.4)
    # [N I] doublet: two Gaussians, common width, free separation
    m = (w > 5194.5) & (w < 5205.5)
    x, y = w[m], f[m]
    p0 = [35, 20, 5198.0, 0.3, np.median(y), 0, 2.36]
    p, cov = curve_fit(g2, x, y, p0=p0, maxfev=40000)
    e = np.sqrt(np.diag(cov))
    res['NI_doublet'] = dict(A1=p[0], A2=p[1], mu1=p[2], mu2=p[2] + p[6], sig=abs(p[3]), sep=p[6], e_mu1=e[2], e_sep=e[6],
                             ratio=p[0] / p[1], rms=float((y - g2(x, *p)).std()))
    # single-Gaussian on strong component only
    res['NI5198_single'] = fit1(w, f, 5196.6, 5199.3, 5198.0, 0.3)
    res['NI5200_single'] = fit1(w, f, 5199.4, 5201.6, 5200.3, 0.3)
    # other emission peaks in 5200-5300 to look for OH lines
    from scipy.signal import find_peaks
    mm = (w > 5205) & (w < 5300)
    pk, prop = find_peaks(f[mm], prominence=6)
    res['other_peaks_5205_5300'] = [(round(float(w[mm][i]), 2), round(float(prop['prominences'][j]), 1)) for j, i in enumerate(pk)]
    out['panel_a_fits'] = {k: ({a: float(b) for a, b in v.items()} if isinstance(v, dict) else v) for k, v in res.items()}

    # calibration-free velocities
    lamN = res['NI_doublet']['mu1']
    lamN2 = res['NI_doublet']['mu2']
    eN = res['NI_doublet']['e_mu1']
    vel = {}
    for line, lo in (('OIII5007', O3_AIR[1]), ('OIII4959', O3_AIR[0])):
        lamO = res[line]['mu']
        eO = res[line]['e_mu']
        for ref, l0n, lamn in (('NI5198', NI_AIR[0], lamN), ('NI5200', NI_AIR[1], lamN2)):
            vt = C * ((lamO / lamn) * (l0n / lo) - 1)
            ev = C * np.hypot(eO / lamO, eN / lamn)
            vel[f'{line}_vs_{ref}'] = dict(V_topo=float(vt), V_helio=float(vt + VCORR), err=float(ev))
        # apparent velocities on the plotted scale, under air and vacuum interpretations
        vel[f'{line}_plot_scale_air'] = float(C * (lamO / lo - 1))
        vel[f'{line}_plot_scale_vac'] = float(C * (lamO / (lo * n_air(lo * 1.00028)) - 1))
    for ref, l0a, l0v, lamn in (('NI5198', NI_AIR[0], NI_VAC[0], lamN), ('NI5200', NI_AIR[1], NI_VAC[1], lamN2)):
        vel[f'{ref}_plot_scale_air'] = float(C * (lamn / l0a - 1))
        vel[f'{ref}_plot_scale_vac'] = float(C * (lamn / l0v - 1))
    vel['NI_sep_measured_vs_true_A'] = [float(res['NI_doublet']['sep']), float(NI_AIR[1] - NI_AIR[0])]
    out['velocities'] = vel

    # panel b (P1-P4): [O III] 5007 position on the plotted scale
    pb = {}
    for k in ('b_P1', 'b_P2', 'b_P3', 'b_P4'):
        if k in spec:
            wb, fb = spec[k]
            o = np.argsort(wb)
            wb, fb = wb[o], fb[o]
            try:
                r5 = fit1(wb, fb, 5003.0, 5012.0, 5007.3)
                pb[k] = dict(mu=float(r5['mu']), e_mu=float(r5['e_mu']), sig=float(r5['sig']), A=float(r5['A']),
                             cont=float(r5['cont']), V_plot_air=float(C * (r5['mu'] / O3_AIR[1] - 1)))
            except Exception as ex:
                pb[k] = dict(err=str(ex))
    out['panel_b'] = pb
    OUTJ.write_text(json.dumps(out, indent=1, default=float))
    np.savez(OUTN, **{f'{k}_wave': v[0] for k, v in spec.items()}, **{f'{k}_flux': v[1] for k, v in spec.items()})
    print(json.dumps({k: out[k] for k in ('calibration', 'markers_A', 'panel_a_npts', 'panel_a_range', 'panel_a_median_step_A')}, indent=0, default=float))
    print(json.dumps(out['panel_a_fits'], indent=0, default=float))
    print(json.dumps(vel, indent=0))
    print(json.dumps(pb, indent=0))


if __name__ == '__main__':
    main()
