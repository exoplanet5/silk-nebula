"""w4_refit_geom: geometry for the refitted truncated shocks: smooth curved limb (wave-2 R_c = 58-66 pc) vs a FLAT
edge-on sheet/ribbon vs a RIPPLED sheet (tangencies of both orientations), for every acceptable southern model of
w4_refit_fit (likelihood A).  Profile/line-profile code from verify_shockmodel_analysis and the C IV Monte Carlo
from verify_shockmodel_rt are imported read-only.

Per model and geometry:
  closure : n0 from the observed ridge width (curved: projected FWHM; flat: FWHM of the emissivity profile j(s);
            rippled: width set by the ripple amplitude, so n0 is fixed by the brightness instead)
  brightness: curved -> I_pk predicted; flat -> LOS depth L = I_obs / j_max needed (G = I_obs / I_face);
            rippled -> n0 such that N_t tangencies (N_t f_overlap = 1, 2, 4) of radius R_r = 10, 20 pc give I_obs
  P/k = n0 mu m_H v^2 / k, t* = N*/(n0 v), t_cool
  profile asymmetry: model profile at the intrinsic FWHM 0.28 deg, blurred to the stacked 0.41 deg, split-normal fit
            sigma_W/sigma_E (W = post-shock side; observed 0.98 +- 0.05, w3_judge_skeptic)
  C IV scattering: tau1548 (face-on, scale-free), G_CIV (curved: projected; flat: L / Delta_CIV), LOS velocity span
            (curved: as wave 2; flat edge-on: 0, and 1 Doppler width for a few-degree tilt) -> S1548, S1551;
            FUV per R = full thin line list (w4_refit_grid UV lines) with C IV x S
  offsets (flat): emission centroids of FUV (after scattering) and [O II] relative to [O III], in deg (+ = east/front)
  line width: flat sheet edge-on (thermal only) and tilted by 5/10 deg; rippled R_r = 10/20 pc; curved 58 pc
Outputs: data/derived/w4_refit_geom.json
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_geom.py
"""
import sys, json, numpy as np
from scipy.optimize import curve_fit
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_analysis as SA
import verify_shockmodel_rt as RT
import verify_shockmodel_flow as F
import w4_refit_common as C

ROOT = '/Users/mickey/astronomy/silk-nebula'
PC, YR, K_B = 3.0857e18, 3.156e7, 1.380649e-16
O_MASS, C_MASS = 15.999 * 1.66054e-24, 12.011 * 1.66054e-24
MASS = F.Network('asplund09').mass_per_H
D_FID, D_SET = 282.0, (252.0, 282.0, 318.0)
W_FID, W_SET, W_STACK = 0.28, (0.15, 0.28, 0.38), 0.41
I_OBS = (1.62, 1.4, 2.5)          # R: DESI reference, range
HC = 6.62607015e-27 * 2.99792458e10


def track(T, tag):
    p = tag + '__'
    return {k[len(p):]: T[k] for k in T.files if k.startswith(p)}


def emiss(tr, i, key, lam, n0, n0ref=0.1):
    """photon emissivity (4pi, ph cm^-3 s^-1) on segment midpoints, depth s (cm) behind the front, at density n0."""
    s = tr['s'][:i + 1] * n0ref / n0
    Fc = tr[key][:i + 1] * n0 / n0ref / (HC / (lam * 1e-8))
    ds = np.diff(s)
    j = np.diff(Fc) / np.where(ds > 0, ds, np.inf)
    return 0.5 * (s[1:] + s[:-1]), j, ds


def fwhm_of(x, y):
    k = int(np.argmax(y)); h = y[k] / 2
    above = np.nonzero(y >= h)[0]
    return float(x[above.max()] - x[above.min()]) if above.size else np.nan


def flat_profile(tr, i, n0):
    sm, j, ds = emiss(tr, i, 'F_O3', 5006.84, n0)
    # resample on a uniform grid (top-hat-like profiles)
    xg = np.linspace(0, sm[-1] + ds[-1] / 2, 2000)
    jg = np.interp(xg, sm, j, left=j[0], right=0.0)
    return xg, jg


def split_normal_ratio(x, y):
    """fit A*exp(-(x-x0)^2/2s^2) with s = sW (x > x0, post-shock side) or sE (x < x0) + const; return sW/sE."""
    def f(x, A, x0, sE, sW, c):
        s = np.where(x < x0, sE, sW)
        return A * np.exp(-0.5 * ((x - x0) / s)**2) + c
    k = int(np.argmax(y))
    p0 = [y[k], x[k], 0.1, 0.1, 0.0]
    try:
        p, _ = curve_fit(f, x, y, p0=p0, maxfev=20000)
        return float(abs(p[3]) / abs(p[2]))
    except Exception:
        return np.nan


def blur_to(xdeg, y, fwhm_tot):
    """convolve a profile (x in deg, uniform) with a Gaussian so that the result has FWHM fwhm_tot."""
    f0 = fwhm_of(xdeg, y)
    sig = np.sqrt(max(fwhm_tot**2 - f0**2, 1e-6)) / 2.3548
    dx = xdeg[1] - xdeg[0]
    k = np.arange(-int(5 * sig / dx) - 1, int(5 * sig / dx) + 2) * dx
    g = np.exp(-0.5 * (k / sig)**2); g /= g.sum()
    pad = len(k)
    yp = np.concatenate([np.zeros(pad), y, np.zeros(pad)])
    xp = np.concatenate([xdeg[0] - dx * np.arange(pad, 0, -1), xdeg, xdeg[-1] + dx * np.arange(1, pad + 1)])
    return xp, np.convolve(yp, g, mode='same')


def asym_flat(tr, i, n0, d):
    xg, jg = flat_profile(tr, i, n0)
    xdeg = np.degrees(xg / (d * PC))          # 0 = front (east); increasing toward the post-shock (west) side
    xp, yp = blur_to(xdeg, jg / jg.max(), W_STACK)
    return split_normal_ratio(xp, yp)


def asym_curved(pr, d):
    hs, I = pr['hs'], pr['I']
    hu = np.linspace(0, hs.max() * 0.999, 4000)
    Iu = np.interp(hu, hs, I)
    xdeg = np.degrees(hu / (d * PC))
    m = xdeg < 1.5
    xp, yp = blur_to(xdeg[m], Iu[m] / Iu.max(), W_STACK)
    return split_normal_ratio(xp, yp)


def centroid(tr, i, key, n0, n0ref=0.1):
    s = tr['s'][:i + 1] * n0ref / n0
    dF = np.diff(tr[key][:i + 1])
    sm = 0.5 * (s[1:] + s[:-1])
    return float(np.sum(sm * dF) / np.sum(dF))


def suppression(tr, i, n0, G_C4, span):
    NC4 = float(tr['NCIV'][i])
    Tc4 = float(tr['wc4_T'][i])
    b = np.sqrt(2 * K_B * Tc4 / C_MASS) / 1e5
    t1548 = 1.497e-15 * 0.190 * 1548.2 * NC4 / b
    t1551 = 1.497e-15 * 0.0952 * 1550.8 * NC4 / b
    S1 = RT.edge_on_suppression(t1548, G_C4, span, nphot=5000)
    S2 = RT.edge_on_suppression(t1551, G_C4, span, nphot=5000)
    return dict(tau1548=float(t1548), b_kms=float(b), G_CIV=float(G_C4), span=float(span), S1548=float(S1),
                S1551=float(S2), S_eff=float(2 / 3 * S1 + 1 / 3 * S2))


def fuv_after(uv, S_eff):
    p = uv['FUV_parts']
    civ = p.get('CIV1549', 0.0)
    return float(uv['FUV_thin'] - civ + civ * S_eff), float(civ * S_eff)


def flat_line(tr, i, v, incl_deg):
    """[O III]5007 LOS profile of an edge-on flat sheet whose normal is tilted by incl from the sky plane."""
    dF = np.diff(tr['F_O3'][:i + 1])
    xm = 0.5 * (tr['x'][:i] + tr['x'][1:i + 1]); Tm = 0.5 * (tr['T'][:i] + tr['T'][1:i + 1])
    u = v * (1 - 1 / xm) * np.sin(np.radians(incl_deg))
    sig = np.sqrt(K_B * Tm / O_MASS) / 1e5
    vg = np.linspace(-120, 120, 2401)
    p = np.sum(dF[:, None] * np.exp(-0.5 * ((vg[None, :] - u[:, None]) / sig[:, None])**2) / sig[:, None], axis=0)
    p /= p.max()
    return fwhm_of(vg, p), float(np.sum(dF * u) / np.sum(dF)), float(np.sum(dF * sig) / np.sum(dF))


def model_geom(rec, tr, d=D_FID, W=W_FID, I_obs=I_OBS[0], full=True):
    v = rec['v']; i = rec['pick']['i']
    n0ref = 0.1
    o = {'tag': rec['tag'], 'v': v, 'xn': rec['xn'], 'b': rec['b'], 'logN': rec['pick']['logN'], 'T_star': rec['pick']['T'],
         'chi2': rec['pick']['chi2'], 'd': d, 'W': W}
    Wcm = np.radians(W) * d * PC
    t_star_n0 = float(tr['t'][i] * n0ref); t_end_n0 = float(rec['t_end_n0'])
    # ---------- flat sheet
    sm, j, ds = emiss(tr, i, 'F_O3', 5006.84, n0ref)
    xg, jg = flat_profile(tr, i, n0ref)
    f_ref = fwhm_of(xg, jg)
    n0 = n0ref * f_ref / Wcm
    xg, jg = flat_profile(tr, i, n0)
    I_face = float(tr['F_O3'][i] * n0 / n0ref / (HC / (5006.84e-8)) * 1e-6)
    jmax = float(jg.max())
    L = I_obs / (jmax * 1e-6)
    sstar = float(tr['s'][i] * n0ref / n0)
    fl = dict(n0=float(n0), P_k=float(n0 * MASS * (v * 1e5)**2 / K_B), I_face=I_face, G=float(I_obs / I_face),
              L_pc=float(L / PC), L_over_len=float(L / PC / (np.radians(16.94) * d)), sstar_pc=sstar / PC,
              t_star_yr=t_star_n0 / n0 / YR, t_cool_yr=t_end_n0 / n0 / YR)
    if full:
        fl['asym_WE'] = asym_flat(tr, i, n0, d)
        sC = float(tr['wc4_s'][i] * n0ref / n0)
        DeltaC = max(2 * sC, 1e-6)
        for lab, span in (('span0', 0.0), ('span1', 1.0)):
            sup = suppression(tr, i, n0, L / DeltaC, span)
            fuv, civ = fuv_after(rec['uv'], sup['S_eff'])
            fl['CIV_' + lab] = dict(sup, FUV=fuv, CIV_CU=civ)
        # offsets (deg, + = toward the front / east)
        c_o3 = centroid(tr, i, 'F_O3', n0)
        c_o2 = centroid(tr, i, 'F_O2', n0)
        c_c4 = centroid(tr, i, 'F_uv_CIV1549', n0)
        parts = rec['uv']['FUV_parts']
        S = fl['CIV_span0']['S_eff']
        # FUV-band centroid: C IV (x S) and the other lines (approx. each by its own centroid when available)
        cents = {'CIV1549': c_c4}
        for kk in ('OIII1663', 'SiIV1400', 'NIV1486', 'OIV1400', 'HeII1640', 'NIII1750', '2g'):
            try:
                cents[kk] = centroid(tr, i, 'F_uv_' + kk, n0)
            except Exception:
                cents[kk] = c_o3
        w = {k: parts.get(k, 0.0) * (S if k == 'CIV1549' else 1.0) for k in cents}
        c_fuv = sum(w[k] * cents[k] for k in cents) / max(sum(w.values()), 1e-30)
        w_thin = {k: parts.get(k, 0.0) for k in cents}
        c_fuv_thin = sum(w_thin[k] * cents[k] for k in cents) / max(sum(w_thin.values()), 1e-30)
        deg = lambda ds_: float(np.degrees(ds_ / (d * PC)))
        fl['offset_FUV_deg'] = deg(c_o3 - c_fuv)
        fl['offset_FUV_thin_deg'] = deg(c_o3 - c_fuv_thin)
        fl['offset_CIV_deg'] = deg(c_o3 - c_c4)
        fl['offset_OII_deg'] = deg(c_o3 - c_o2)
        for inc in (0.0, 5.0, 10.0):
            fw, vc, sg = flat_line(tr, i, v, inc)
            fl[f'line_i{int(inc)}'] = dict(fwhm=fw, v_shift=vc, sig_th=sg)
    o['flat'] = fl
    # ---------- curved limb (wave-2 closure W)
    Rc = d * np.radians(11.9) * PC
    pr = SA.profile(tr, i, n0ref, n0ref, Rc)
    n0c = n0ref * pr['fwhm'] / Wcm
    for _ in range(6):
        pr = SA.profile(tr, i, n0c, n0ref, Rc)
        n0c *= pr['fwhm'] / Wcm
    pr = SA.profile(tr, i, n0c, n0ref, Rc)
    cu = dict(n0=float(n0c), P_k=float(n0c * MASS * (v * 1e5)**2 / K_B), I_pk=pr['I_pk'], G=pr['G'],
              t_star_yr=t_star_n0 / n0c / YR)
    if full:
        cu['asym_WE'] = asym_curved(pr, d)
        lp = SA.line_profile(tr, i, n0c, n0ref, Rc, v, pr['h_pk'])
        cu['line_fwhm'] = lp['fwhm']
        prC = SA.profile(tr, i, n0c, n0ref, Rc, key='F_CIV', lam=1549.5)
        sC = float(tr['wc4_s'][i] * n0ref / n0c)
        Tc4 = float(tr['wc4_T'][i]); b = np.sqrt(2 * K_B * Tc4 / C_MASS) / 1e5
        span = v * float(tr['wc4_u'][i]) * (prC['G'] * 2 * sC) / Rc / b
        sup = suppression(tr, i, n0c, prC['G'], span)
        fuv, civ = fuv_after(rec['uv'], sup['S_eff'])
        cu['CIV'] = dict(sup, FUV=fuv, CIV_CU=civ)
        cu['offset_CIV_peak_deg'] = float(np.degrees((pr['h_pk'] - prC['h_pk']) / (d * PC)))
    o['curved'] = cu
    # ---------- rippled sheet: tangencies of radius R_r; brightness fixes n0 (width = ripple amplitude)
    if full:
        rp = {}
        for Rr in (10.0, 20.0):
            for Nt in (1.0, 2.0, 4.0):
                n0r = n0
                for _ in range(12):
                    prr = SA.profile(tr, i, n0r, n0ref, Rr * PC)
                    Ipred = Nt * prr['I_pk']
                    n0r *= (I_obs / Ipred)**(1 / 1.5)
                prr = SA.profile(tr, i, n0r, n0ref, Rr * PC)
                key = f'R{int(Rr)}_N{int(Nt)}'
                rp[key] = dict(n0=float(n0r), P_k=float(n0r * MASS * (v * 1e5)**2 / K_B), G1=prr['G'],
                               fwhm1_deg=float(np.degrees(prr['fwhm'] / (d * PC))), t_star_yr=t_star_n0 / n0r / YR)
            lpr = SA.line_profile(tr, i, rp[f'R{int(Rr)}_N2']['n0'], n0ref, Rr * PC, v,
                                  SA.profile(tr, i, rp[f'R{int(Rr)}_N2']['n0'], n0ref, Rr * PC)['h_pk'])
            rp[f'R{int(Rr)}_line_fwhm'] = lpr['fwhm']
            # C IV for a ripple tangency: projected G_CIV, gradient as for the curved limb, at N_t = 2
            n2 = rp[f'R{int(Rr)}_N2']['n0']
            prC = SA.profile(tr, i, n2, n0ref, Rr * PC, key='F_CIV', lam=1549.5)
            sC = float(tr['wc4_s'][i] * n0ref / n2)
            Tc4 = float(tr['wc4_T'][i]); b = np.sqrt(2 * K_B * Tc4 / C_MASS) / 1e5
            span = v * float(tr['wc4_u'][i]) * (prC['G'] * 2 * sC) / (Rr * PC) / b
            sup = suppression(tr, i, n2, prC['G'], span)
            rp[f'R{int(Rr)}_CIV'] = dict(sup, FUV=fuv_after(rec['uv'], sup['S_eff'])[0])
        o['ripple'] = rp
    return o


def rngs(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    return [float(v.min()), float(np.percentile(v, 16)), float(np.median(v)), float(np.percentile(v, 84)), float(v.max())] if v.size else None


if __name__ == '__main__':
    fit = json.load(open(ROOT + '/data/derived/w4_refit_fit.json'))
    T = np.load(ROOT + '/data/derived/w4_refit_tracks.npz')
    byt = {r['tag']: r for r in fit['runs']}
    VA = fit['variants']['A']
    mn, mr = VA['chi2_min_nat'], VA['chi2_min_R88']
    out = {'notes': __doc__, 'models': []}
    for tag in VA['acceptable']:
        r = dict(byt[tag])
        r['pick'] = r['A_nat'] if r['A_nat']['chi2'] - mn <= r['A_R88']['chi2'] - mr else r['A_R88']
        tr = track(T, tag)
        o = model_geom(r, tr)
        o['uv'] = {k: r['uv'][k] for k in ('FUV_thin', 'FUV_nonres', 'NUV')}
        o['uv_parts'] = r['uv']['FUV_parts']
        o['nuv_parts'] = r['uv']['NUV_parts']
        o['ratios'] = r['pick']['ratios']
        # distance / width variants (flat and curved closures, no RT)
        o['variants'] = {}
        for d in D_SET:
            for W in W_SET:
                g = model_geom(r, tr, d=d, W=W, full=False)
                o['variants'][f'd{int(d)}_W{W}'] = dict(flat=g['flat'], curved=g['curved'])
        out['models'].append(o)
        f, c = o['flat'], o['curved']
        print(f"{tag:28s} chi2={o['chi2']:.1f} | flat n0={f['n0']:.3f} P={f['P_k']:.2e} L={f['L_pc']:.1f}pc G={f['G']:.0f} "
              f"asym={f['asym_WE']:.2f} FUV0={f['CIV_span0']['FUV']:.0f} FUV1={f['CIV_span1']['FUV']:.0f} "
              f"offFUV={f['offset_FUV_deg']:+.3f} offOII={f['offset_OII_deg']:+.3f} FWHM0={f['line_i0']['fwhm']:.0f} | "
              f"curved n0={c['n0']:.3f} I={c['I_pk']:.2f} asym={c['asym_WE']:.2f} FUV={c['CIV']['FUV']:.0f} "
              f"FWHM={c['line_fwhm']:.0f} | NUV={o['uv']['NUV']:.0f} R10N2 P={o['ripple']['R10_N2']['P_k']:.2e} "
              f"FWHM={o['ripple']['R10_line_fwhm']:.0f}", flush=True)
    M = out['models']
    g = lambda path: rngs([eval('m' + path) for m in M])
    out['summary'] = {
        'n_models': len(M),
        'flat': {k: g(f"['flat']['{k}']") for k in ('n0', 'P_k', 'L_pc', 'L_over_len', 'G', 'sstar_pc', 't_star_yr', 't_cool_yr',
                                                 'asym_WE', 'offset_FUV_deg', 'offset_FUV_thin_deg', 'offset_CIV_deg', 'offset_OII_deg')},
        'flat_FUV_span0': g("['flat']['CIV_span0']['FUV']"), 'flat_FUV_span1': g("['flat']['CIV_span1']['FUV']"),
        'flat_S_span0': g("['flat']['CIV_span0']['S_eff']"), 'flat_GCIV': g("['flat']['CIV_span0']['G_CIV']"),
        'flat_tau1548': g("['flat']['CIV_span0']['tau1548']"),
        'flat_line_fwhm_i0': g("['flat']['line_i0']['fwhm']"), 'flat_line_fwhm_i5': g("['flat']['line_i5']['fwhm']"),
        'flat_line_fwhm_i10': g("['flat']['line_i10']['fwhm']"), 'flat_vshift_i5': g("['flat']['line_i5']['v_shift']"),
        'flat_vshift_i10': g("['flat']['line_i10']['v_shift']"),
        'curved': {k: g(f"['curved']['{k}']") for k in ('n0', 'P_k', 'I_pk', 'G', 'asym_WE', 'line_fwhm', 't_star_yr', 'offset_CIV_peak_deg')},
        'curved_FUV': g("['curved']['CIV']['FUV']"), 'curved_S': g("['curved']['CIV']['S_eff']"),
        'FUV_thin': g("['uv']['FUV_thin']"), 'FUV_nonres': g("['uv']['FUV_nonres']"), 'NUV': g("['uv']['NUV']"),
        'SiIV': g("['uv_parts']['SiIV1400']"), 'CIV_thin': g("['uv_parts']['CIV1549']"),
        'ripple': {k: g(f"['ripple']['{k}']['P_k']") for k in ('R10_N1', 'R10_N2', 'R10_N4', 'R20_N1', 'R20_N2', 'R20_N4')},
        'ripple_n0': {k: g(f"['ripple']['{k}']['n0']") for k in ('R10_N2', 'R20_N2')},
        'ripple_line_fwhm': {k: g(f"['ripple']['{k}_line_fwhm']") for k in ('R10', 'R20')},
        'ripple_FUV': {k: g(f"['ripple']['{k}_CIV']['FUV']") for k in ('R10', 'R20')},
        'by_distance_flat_P': {f'd{int(d)}': g(f"['variants']['d{int(d)}_W0.28']['flat']['P_k']") for d in D_SET},
        'by_distance_flat_L': {f'd{int(d)}': g(f"['variants']['d{int(d)}_W0.28']['flat']['L_pc']") for d in D_SET},
        'by_distance_flat_tstar': {f'd{int(d)}': g(f"['variants']['d{int(d)}_W0.28']['flat']['t_star_yr']") for d in D_SET},
        'by_width_flat_P': {f'W{W}': g(f"['variants']['d282_W{W}']['flat']['P_k']") for W in W_SET},
        'by_width_flat_L': {f'W{W}': g(f"['variants']['d282_W{W}']['flat']['L_pc']") for W in W_SET},
        'by_distance_curved_I': {f'd{int(d)}': g(f"['variants']['d{int(d)}_W0.28']['curved']['I_pk']") for d in D_SET},
    }
    json.dump(out, open(ROOT + '/data/derived/w4_refit_geom.json', 'w'), indent=1,
              default=lambda x: float(x) if np.isscalar(x) else None)
    for k, v in out['summary'].items():
        print(k, v)
