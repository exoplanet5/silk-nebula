"""verify_shockmodel_analysis: physical consequences of the truncated-shock fits (tasks b-f of the shockmodel agent).

Reads data/derived/verify_shockmodel_fit.json + verify_shockmodel_tracks.npz (tracks at n0_ref = 0.1; the flow is
scale-free in n0 at fixed b = B0/sqrt(n0), so s, t scale as 1/n0 and fluxes as n0).
For every acceptable model (chi2 <= global min + 4 with either the native or the 'R88-calibrated Hb' scoring):
  (b) edge-on brightness: projected profile of the model emissivity j(s) on a curved sheet of radius R_c
      (thin-shell chord z^2 = 2 R_c (h - s)); peak intensity I_pk and FWHM.  Two closures:
        W: n0 fixed by matching the observed filament FWHM (0.28 deg at d) -> predicted I_pk and the implied
           ram pressure rho0 v^2;
        P: n0 fixed by pressure balance rho0 v_s^2 = P_drive (Monogem) -> predicted I_pk and FWHM.
  (c) time since the shock entered the gas (t*) and full cooling time t(T=1e4 K) vs the Monogem age (1.1e5 yr).
  (d) C IV 1548/1551, O III] 1661/66, He II 1640 in the GALEX FUV band per rayleigh of [O III]5007, optically thin and
      with the Monte-Carlo resonance-scattering suppression of C IV for the edge-on geometry; compared with the
      measured 92 +- 10 CU per R (south, verify_xrayfuv).
  (e) [O III]5007 line profile at the projected peak: radial post-shock flow u = v_s (1 - 1/x) projected on the line
      of sight along the chord + thermal O broadening.
  (f) proper motion of the shock front.
  Also: [O II]3729/3726 predicted at the closure density; Ha surface brightness.
Outputs: data/derived/verify_shockmodel_analysis.json
"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_atomic as A
import verify_shockmodel_flow as F
import verify_shockmodel_rt as RT

ROOT = '/Users/mickey/astronomy/silk-nebula'
PC = 3.0857e18
YR = 3.156e7
K_B = 1.380649e-16
M_H = 1.6735575e-24
RAY = 1e6                       # photons cm^-2 s^-1 (4pi sr) per rayleigh -> I_R = F_ph * 1e-6 (face-on flux F_ph)
D_PC = {'lo': 250.0, 'fid': 280.0, 'hi': 320.0}
FWHM_DEG = 0.28
I_OBS = {'fid': 1.6, 'lo': 1.2, 'hi': 2.3}   # [O III]5007 ridge peak, R (DESI south median 1.6; E strand 1.4-2.3)
P_DRIVE = {'lo': 3e4, 'fid': 1e5, 'hi': 4e5}  # K cm^-3 (Monogem, see note)
AGE_MONOGEM = 1.1e5
FUV_MEAS = (92.0, 10.0)        # CU per R(5007), south, verify_xrayfuv
EBV = {'fg280': 0.006, 'sfd': 0.03}
O_MASS = 15.999 * 1.66054e-24
C_MASS = 12.011 * 1.66054e-24


def load():
    d = json.load(open(ROOT + '/data/derived/verify_shockmodel_fit.json'))
    T = np.load(ROOT + '/data/derived/verify_shockmodel_tracks.npz')
    return d, T


def track(T, tag):
    return {k.split('__', 1)[1]: T[k] for k in T.files if k.startswith(tag + '__')}


def galex_weff():
    a = np.loadtxt(ROOT + '/data/raw/verify_xrayfuv_lit/svo_GALEX.FUV.dat')
    lam, Ae = a[:, 0], a[:, 1]
    tot = np.trapezoid(Ae, lam)
    return lambda l: tot / np.interp(l, lam, Ae)


def ccm_k(lam_A, Rv=3.1):
    x = 1e4 / lam_A
    if x >= 3.3:
        Fa = -0.04473 * (x - 5.9)**2 - 0.009779 * (x - 5.9)**3 if x > 5.9 else 0.0
        Fb = 0.2130 * (x - 5.9)**2 + 0.1207 * (x - 5.9)**3 if x > 5.9 else 0.0
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67)**2 + 0.341) + Fa
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62)**2 + 0.263) + Fb
    else:
        y = x - 1.82
        a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    return (a + b / Rv) * Rv


def profile(tr, istar, n0, n0ref, Rc_cm, key='F_O3', lam=5006.84):
    """Projected intensity (R) vs depth h below the front for a sheet of curvature radius Rc, model truncated at istar."""
    s = tr['s'][:istar + 1] * n0ref / n0
    Fc = tr[key][:istar + 1] * n0 / n0ref / (A.H_PL * A.C_L / (lam * 1e-8))  # photon flux, cumulative
    ds = np.diff(s)
    jmid = np.diff(Fc) / np.where(ds > 0, ds, np.inf)       # photon emissivity (4pi) in each segment
    smid = 0.5 * (s[1:] + s[:-1])
    sstar = s[-1]
    I_face = Fc[-1] * 1e-6
    hs = np.concatenate([np.linspace(0, sstar, 120)[1:], sstar * np.logspace(0.0, 2.0, 120)[1:]])
    I = np.zeros(hs.size)
    for k, h in enumerate(hs):
        umin, umax = np.sqrt(max(h - sstar, 0.0)), np.sqrt(h)
        u = np.linspace(umin, umax, 400)
        sv = h - u**2
        jj = np.interp(sv, smid, jmid, left=jmid[0], right=0.0)
        jj[sv > sstar] = 0.0
        I[k] = 2 * np.sqrt(2 * Rc_cm) * np.trapezoid(jj, u) * 1e-6
    kp = int(np.argmax(I))
    half = I[kp] / 2
    inner = hs[kp:][np.argmax(I[kp:] < half)] if np.any(I[kp:] < half) else hs[-1]
    outer = hs[:kp + 1][np.argmax(I[:kp + 1] >= half)]
    return dict(I_pk=float(I[kp]), h_pk=float(hs[kp]), fwhm=float(inner - outer), I_face=float(I_face),
                G=float(I[kp] / I_face), sstar=float(sstar), hs=hs, I=I)


def line_profile(tr, istar, n0, n0ref, Rc_cm, v_s, h):
    """[O III]5007 LOS velocity profile at depth h (km/s grid), flow u = v_s(1-1/x), thermal O broadening."""
    s = tr['s'][:istar + 1] * n0ref / n0
    Fc = tr['F_O3'][:istar + 1]
    ds = np.diff(s)
    jmid = np.diff(Fc) / np.where(ds > 0, ds, np.inf)
    smid = 0.5 * (s[1:] + s[:-1])
    xm = 0.5 * (tr['x'][:istar] + tr['x'][1:istar + 1])
    Tm = 0.5 * (tr['T'][:istar] + tr['T'][1:istar + 1])
    sstar = s[-1]
    u = np.linspace(np.sqrt(max(h - sstar, 0.0)), np.sqrt(h), 600)
    sv = h - u**2
    jj = np.interp(sv, smid, jmid, left=jmid[0], right=0.0)
    jj[sv > sstar] = 0.0
    z = np.sqrt(2 * Rc_cm) * u
    uf = v_s * (1 - 1 / np.interp(sv, smid, xm))
    Tz = np.interp(sv, smid, Tm)
    vgrid = np.linspace(-150, 150, 1201)
    prof = np.zeros_like(vgrid)
    sig = np.sqrt(K_B * Tz / O_MASS) / 1e5
    for sign in (1, -1):
        vl = sign * uf * z / Rc_cm
        prof += np.sum(jj[:, None] * np.exp(-0.5 * ((vgrid[None, :] - vl[:, None]) / sig[:, None])**2) / sig[:, None], axis=0)
    prof /= prof.max()
    above = vgrid[prof >= 0.5]
    vlos_max = float(np.max(uf * z / Rc_cm))
    return dict(fwhm=float(above[-1] - above[0]), vlos_max=vlos_max, sig_th=float(np.average(sig, weights=jj + 1e-300)),
                v=vgrid, p=prof)


def o2_doublet(ne, T):
    at = F.atoms()['O2']
    return float(at.emiss(T, ne, 2, 1) / at.emiss(T, ne, 3, 1))


def analyse(rec, tr, n0ref, Weff, Rc_pc=None, d_pc=280.0, P=1e5, verbose=False, do_rt=True):
    v = rec['v']
    mass = F.Network('asplund09').mass_per_H if analyse.mass is None else analyse.mass
    analyse.mass = mass
    istar = rec['best']['i']
    Rc_pc = Rc_pc if Rc_pc is not None else d_pc * np.radians(11.9)
    Rc = Rc_pc * PC
    out = {'tag': rec['tag'], 'v': v, 'xn': rec['xn'], 'b': rec['b'], 'chi2': rec['best']['chi2'],
           'logN': float(np.log10(rec['best']['N'])), 'T_star': rec['best']['T'], 'Rc_pc': Rc_pc, 'd_pc': d_pc}
    fwhm_obs = np.radians(FWHM_DEG) * d_pc * PC
    # --- closure W: match the observed width
    pr_ref = profile(tr, istar, n0ref, n0ref, Rc)
    # FWHM scales exactly as 1/n0 only for fixed Rc/s; iterate
    n0 = n0ref * pr_ref['fwhm'] / fwhm_obs
    for _ in range(6):
        pr = profile(tr, istar, n0, n0ref, Rc)
        n0 *= pr['fwhm'] / fwhm_obs
    pr = profile(tr, istar, n0, n0ref, Rc)
    Pk = n0 * mass * (v * 1e5)**2 / K_B
    t_star = rec['best']['t_n0'] / n0 / YR
    t_cool = rec['t_end_n0'] / n0 / YR
    out['W'] = dict(n0=n0, I_pk=pr['I_pk'], I_face=pr['I_face'], G=pr['G'], fwhm_pc=pr['fwhm'] / PC,
                    sstar_pc=pr['sstar'] / PC, P_ram_k=Pk, t_star_yr=t_star, t_cool_yr=t_cool)
    lp = line_profile(tr, istar, n0, n0ref, Rc, v, pr['h_pk'])
    out['W']['line_fwhm'] = lp['fwhm']; out['W']['vlos_max'] = lp['vlos_max']; out['W']['sig_th'] = lp['sig_th']
    ne_o2 = tr['wo2_ne'][istar] * n0 / n0ref
    out['W']['ne_OII'] = float(ne_o2)
    out['W']['OII_doublet'] = o2_doublet(max(ne_o2, 1e-3), float(tr['wo2_T'][istar]))
    # --- closure P: pressure balance
    n0P = P * K_B / (mass * (v * 1e5)**2)
    prP = profile(tr, istar, n0P, n0ref, Rc)
    lpP = line_profile(tr, istar, n0P, n0ref, Rc, v, prP['h_pk'])
    out['P'] = dict(P_k=P, n0=n0P, I_pk=prP['I_pk'], I_face=prP['I_face'], G=prP['G'], fwhm_pc=prP['fwhm'] / PC,
                    fwhm_deg=float(np.degrees(prP['fwhm'] / (d_pc * PC))), sstar_pc=prP['sstar'] / PC,
                    t_star_yr=rec['best']['t_n0'] / n0P / YR, t_cool_yr=rec['t_end_n0'] / n0P / YR,
                    line_fwhm=lpP['fwhm'], vlos_max=lpP['vlos_max'],
                    ne_OII=float(tr['wo2_ne'][istar] * n0P / n0ref))
    out['P']['OII_doublet'] = o2_doublet(max(out['P']['ne_OII'], 1e-3), float(tr['wo2_T'][istar]))
    # --- (d) FUV
    hc = A.H_PL * A.C_L
    ph = lambda key, lam: tr[key][istar] / (hc / (lam * 1e-8))
    p5007 = ph('F_O3', 5006.84)
    # split C IV doublet by the model's own components (F_CIV is the sum); 1548:1551 ~ 2:1 when thin
    pC4 = ph('F_CIV', 1549.5)
    p1548, p1551 = pC4 * 2 / 3, pC4 / 3
    p1663 = ph('F_O3s', 1664.0)
    p1640 = ph('F_HeII1640', 1640.4)
    cu = lambda lam: 79577.47 / Weff(lam)
    thin = {'CIV': (p1548 + p1551) / p5007 * cu(1549.5), 'OIII]': p1663 / p5007 * cu(1664.0),
            'HeII1640': p1640 / p5007 * cu(1640.4)}
    out['FUV_thin_CU_per_R'] = thin
    # resonance scattering: C IV column (face-on) and the C IV-weighted geometry
    NC4 = float(tr['NCIV'][istar])
    Tc4 = float(tr['wc4_T'][istar])
    b_kms = np.sqrt(2 * K_B * Tc4 / C_MASS) / 1e5
    prC = profile(tr, istar, n0, n0ref, Rc, key='F_CIV', lam=1549.5)
    sC = float(tr['wc4_s'][istar] * n0ref / n0)          # C IV-weighted depth behind the front
    DeltaC = max(2 * sC, 1e-6)
    GC = prC['G']
    uC4 = v * float(tr['wc4_u'][istar])
    span = uC4 * (GC * DeltaC) / Rc / b_kms                 # LOS velocity span across the chord, Doppler units
    tau1548 = 1.497e-15 * 0.190 * 1548.2 * NC4 / b_kms
    tau1551 = 1.497e-15 * 0.0952 * 1550.8 * NC4 / b_kms
    S1548 = RT.edge_on_suppression(tau1548, GC, span, nphot=12000) if do_rt else np.nan
    S1551 = RT.edge_on_suppression(tau1551, GC, span, nphot=12000) if do_rt else np.nan
    # relative brightening of C IV vs [O III] at the [O III] peak position (different layer thickness)
    rel_G = GC / pr['G']
    civ_rt = (p1548 * S1548 + p1551 * S1551) / p5007 * cu(1549.5)
    out['FUV_RT'] = dict(NCIV=NC4, b_kms=float(b_kms), T_CIV=Tc4, tau1548_face=float(tau1548), G_CIV=GC,
                         G_OIII=pr['G'], span_dopp=float(span), S1548=S1548, S1551=S1551,
                         CIV_CU_per_R=float(civ_rt), CIV_peak_geom_factor=float(rel_G),
                         total_CU_per_R=float(civ_rt + thin['OIII]'] + thin['HeII1640']),
                         total_thin_CU_per_R=float(sum(thin.values())))
    for lab, e in EBV.items():
        att = 10**(-0.4 * e * (ccm_k(1549.5) - ccm_k(5006.84)))
        out['FUV_RT'][f'total_CU_per_R_ebv_{lab}'] = float(out['FUV_RT']['total_CU_per_R'] * att)
    # (f) proper motion
    out['pm_mas_yr'] = float(v / (4.74 * d_pc) * 1e3)
    out['eta5007'] = rec['best']['eta5007']
    out['ratios'] = rec['best']['ratios']
    # Ha brightness at the W closure
    out['W']['I_Ha_R'] = float(pr['I_pk'] / rec['best']['ratios']['O3_Ha'] * (6562.8 / 5006.84))   # photon units
    return out


analyse.mass = None

if __name__ == '__main__':
    d, T = load()
    n0ref = d['n0_ref']
    Weff = galex_weff()
    runs = [r for r in d['runs'] if not r['complete']]
    c2min = min(r['best']['chi2'] for r in runs)
    c2minR = min(r['best_hbR88']['chi2'] for r in runs)
    acc = [r for r in runs if r['best']['chi2'] <= c2min + 4 or r['best_hbR88']['chi2'] <= c2minR + 4]
    print(f'global chi2 min native {c2min:.2f}, R88-Hb {c2minR:.2f}; acceptable runs {len(acc)} of {len(runs)}')
    res = {'c2min': c2min, 'c2minR88': c2minR, 'acceptable': [], 'inputs': dict(D_PC=D_PC, FWHM_DEG=FWHM_DEG, I_OBS=I_OBS,
           P_DRIVE=P_DRIVE, AGE=AGE_MONOGEM, FUV_MEAS=FUV_MEAS, EBV=EBV,
           Weff_1549=float(Weff(1549.5)), Weff_1664=float(Weff(1664.0)), CU_per_R_1549=float(79577.47 / Weff(1549.5)),
           k1549=ccm_k(1549.5), k5007=ccm_k(5006.84))}
    for r in sorted(acc, key=lambda r: r['best']['chi2']):
        tr = track(T, r['tag'])
        o = analyse(r, tr, n0ref, Weff)
        # variants: pressure range, local ripple curvature, distance range
        o['P_lo'] = analyse(r, tr, n0ref, Weff, P=P_DRIVE['lo'], do_rt=False)['P']
        o['P_hi'] = analyse(r, tr, n0ref, Weff, P=P_DRIVE['hi'], do_rt=False)['P']
        o['Rc15'] = {k: analyse(r, tr, n0ref, Weff, Rc_pc=15.0, do_rt=False)[k] for k in ('W', 'P')}
        o['d250'] = analyse(r, tr, n0ref, Weff, d_pc=250.0, do_rt=False)['W']
        o['d320'] = analyse(r, tr, n0ref, Weff, d_pc=320.0, do_rt=False)['W']
        res['acceptable'].append(o)
        W, P = o['W'], o['P']
        print(f"{o['tag']:34s} chi2={o['chi2']:.1f} logN*={o['logN']:.2f} T*={o['T_star']:.3g} eta={o['eta5007']:.3f} | "
              f"W: n0={W['n0']:.3f} I={W['I_pk']:.2f}R G={W['G']:.0f} s*={W['sstar_pc']:.3f}pc P/k={W['P_ram_k']:.2e} "
              f"t*={W['t_star_yr']:.2e} tc={W['t_cool_yr']:.2e} FWHMv={W['line_fwhm']:.0f} neOII={W['ne_OII']:.2f} dbl={W['OII_doublet']:.3f} | "
              f"P1e5: n0={P['n0']:.3f} I={P['I_pk']:.2f}R fwhm={P['fwhm_deg']:.3f}deg t*={P['t_star_yr']:.2e} FWHMv={P['line_fwhm']:.0f} | "
              f"FUV thin {o['FUV_RT']['total_thin_CU_per_R']:.0f} RT {o['FUV_RT']['total_CU_per_R']:.0f} CU/R "
              f"(tau1548={o['FUV_RT']['tau1548_face']:.2f} G_C4={o['FUV_RT']['G_CIV']:.0f} S={o['FUV_RT']['S1548']:.2f}) "
              f"pm={o['pm_mas_yr']:.0f}", flush=True)
    json.dump(res, open(ROOT + '/data/derived/verify_shockmodel_analysis.json', 'w'), indent=1,
              default=lambda x: float(x) if np.isscalar(x) else None)
