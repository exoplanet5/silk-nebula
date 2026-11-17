"""verify_paper_energetics.py -- tasks (c) energetics, (d) [O III] self-consistency, (e) mag -> Rayleigh (tag: paper).

(c) Dopita (1979) n_[SII] = 45 n0 (Vs/100)^2; Cioffi, McKee & Bertschinger (1988) R_PDS, t_PDS;
    Boumis+2022 / Cioffi E = 6.8e43 n0^1.16 Vs^1.35 Rs^3.16 erg. Re-derived at 160 pc, 280 pc, 1.6 kpc and with n0
    tied to the revised brightness (complete-shock H-beta yield bound; [O III] yield fraction).
(d) own 5-level O++ atom: Aggarwal & Keenan 1999 collision strengths (AK99, to 2e5 K; cross-check Storey, Sochi &
    Badnell 2014 to 2.5e4 K), A-values Froese Fischer+2004 / Storey & Zeippen 2000 (PyNeb o_iii_atom_FFT04-SZ00),
    energies NIST (PyNeb levels file). H-alpha: case-B alpha_eff = 1.17e-13 T4^-0.9 (Haffner+2003 EM relation).
(e) AB / Vega, effective width, [O III] 4959 leakage, interference-filter angle shift (f/0.95 converging beam,
    collimated off-axis field, f/2.8 reflector), Gaia integrated starlight, sky-fraction of the claimed signal.
Run: ~/.venvs/astro313/bin/python -I scripts/verify_paper_energetics.py
"""
import sys, json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/derived/verify_paper_raw/pyneb_atomic'
sys.path.insert(0, str(ROOT / 'scripts'))
PC = 3.0857e18
DEG = np.pi / 180

# ------------------------------------------------------------------ (c) energetics
def vs_dopita(ne_sii, n0):
    return 100.0 * np.sqrt(ne_sii / (45.0 * n0))


def E_cioffi(n0, vs, rs):          # erg; vs km/s, rs pc  (xi_m = 1)
    return 6.8e43 * n0 ** 1.16 * vs ** 1.35 * rs ** 3.16


def R_pds(E, n0):
    return 14.0 * (E / 1e51) ** (2 / 7) * n0 ** (-3 / 7)


def t_age(E, n0, rs, plus=1 / 3):
    return 9.98e3 * (E / 1e51) ** (3 / 14) * n0 ** (-4 / 7) * ((rs / R_pds(E, n0)) ** (10 / 3) + plus)


def t_pds(E, n0):
    return 13.3e3 * (E / 1e51) ** (3 / 14) * n0 ** (-4 / 7)


def v_cioffi(E, n0, t):            # PDS-phase shell speed (Cioffi+88 eq. 3.29-like): v = (2/5) R_PDS/t_PDS (4t/3t_PDS - 1/3)^-7/10
    rp, tp = R_pds(E, n0) * PC, t_pds(E, n0) * 3.156e7
    return 0.4 * rp / tp * (4 * t * 3.156e7 / (3 * tp) - 1 / 3) ** (-0.7) / 1e5


def energetics():
    out = {}
    # paper reproduction
    n0, ne, rs = 1.0, 64.0, 30.0
    vs = vs_dopita(ne, n0)
    E = E_cioffi(n0, vs, rs)
    out['paper_repro'] = dict(Vs=round(vs, 1), E=float('%.3g' % E), R_PDS=round(R_pds(E, n0), 2),
                              age_plus_third=float('%.3g' % t_age(E, n0, rs)), age_plus_one_as_printed=float('%.3g' % t_age(E, n0, rs, 1.0)),
                              age_plus_third_RPDS18=float('%.3g' % (9.98e3 * (E / 1e51) ** (3 / 14) * ((rs / 18.0) ** (10 / 3) + 1 / 3))),
                              E_needed_for_RPDS18=float('%.3g' % (1e51 * (18.0 / 14.0) ** 3.5)),
                              v_cioffi_at_age=round(v_cioffi(E, n0, t_age(E, n0, rs)), 1))
    # error propagation like the paper: n0 = 1 +- 0.3, Vs = 119 +- 5, Rs fixed
    rng = np.random.default_rng(3)
    n0s = np.clip(rng.normal(1.0, 0.3, 20000), 0.05, None); vss = rng.normal(119.0, 5.0, 20000)
    Es = E_cioffi(n0s, vss, 30.0)
    out['paper_error_mc'] = dict(E_16_50_84=[float('%.3g' % q) for q in np.percentile(Es, [16, 50, 84])],
                                 age_16_50_84=[float('%.3g' % q) for q in np.percentile(t_age(Es, n0s, 30.0), [16, 50, 84])])
    theta_paper = 30.0 / 160.0 / DEG          # 10.74 deg
    theta_mono = 12.10                         # deg, spine distance from Knies+24 two-SNR Monogem centre (context teammate)
    grid = {}
    for d in (160.0, 280.0, 1600.0):
        for tname, th in (('paper_10.74deg', theta_paper), ('Monogem_12.1deg', theta_mono)):
            rs = d * th * DEG
            for ne_ in (64.0, 150.0):
                for n0 in (1.0, 0.3, 0.1, 0.03, 0.01):
                    vs = vs_dopita(ne_, n0)
                    E = E_cioffi(n0, vs, rs)
                    grid[f'd{int(d)}_{tname}_ne{int(ne_)}_n0{n0}'] = dict(Rs=round(rs, 1), Vs=round(vs, 0), E=float('%.3g' % E),
                                                                           R_over_RPDS=round(rs / R_pds(E, n0), 2), age=float('%.3g' % t_age(E, n0, rs)))
            for vs in (100.0, 150.0):
                for n0 in (1.0, 0.1, 0.03, 0.01, 0.003, 0.001):
                    E = E_cioffi(n0, vs, rs)
                    grid[f'd{int(d)}_{tname}_Vs{int(vs)}_n0{n0}'] = dict(Rs=round(rs, 1), Vs=vs, E=float('%.3g' % E),
                                                                          R_over_RPDS=round(rs / R_pds(E, n0), 2), age=float('%.3g' % t_age(E, n0, rs)))
    out['grid'] = grid
    # n0 tied to brightness (radiative shock), d = 280 pc, Rs = 59.1 pc
    d, rs = 280.0, 280.0 * theta_mono * DEG
    w = d * 0.28 * DEG
    Lb_geom = 2 * np.sqrt(2 * rs / w)
    yieldHb_min = 0.117          # case-B Hbeta photons per recombination (1e4 K); >= 1 recombination per swept-up H in a complete shock
    bright = {}
    for vs in (100.0, 150.0):
        v100 = vs / 100.0
        for Lb in (3.0, 10.0, round(Lb_geom, 1)):
            for IHb_name, IHb in (('SDSS_south_2sigma_0.057R', 0.057), ('DESI_implied_0.12R', 0.12)):
                n0max = IHb / (Lb * 1e-6 * 1e7 * v100 * yieldHb_min)      # >=1 recombination per H
                n0max2 = n0max / 1.6                                          # ~1.6 recombinations per H (DS96-like)
                E = E_cioffi(n0max, vs, rs)
                bright[f'Hb_complete_shock_Vs{int(vs)}_Lb{Lb}_{IHb_name}'] = dict(n0_max=float('%.3g' % n0max), n0_max_1p6rec=float('%.3g' % n0max2),
                                                                                    E_max=float('%.3g' % E), R_over_RPDS=round(rs / R_pds(E, n0max), 2),
                                                                                    age=float('%.3g' % t_age(E, n0max, rs)))
            for I5007 in (2.0, 4.0, 11.5):
                for eps in (0.02, 0.1):
                    # face-on 5007 photons per swept-up H = eps * (1/2) mu_H m_H Vs^2 / h nu(5007)
                    ph_per_H = eps * 0.5 * 1.4 * 1.6726e-24 * (vs * 1e5) ** 2 / 3.967e-12
                    for Lb in (3.0, 10.0, round(Lb_geom, 1)):
                        n0 = I5007 / (Lb * 1e-6 * 1e5 * vs * ph_per_H)
                        E = E_cioffi(n0, vs, rs)
                        bright[f'OIII_{I5007}R_eps{eps}_Vs{int(vs)}_Lb{Lb}'] = dict(n0=float('%.3g' % n0), E=float('%.3g' % E),
                                                                                    R_over_RPDS=round(rs / R_pds(E, n0), 2), age=float('%.3g' % t_age(E, n0, rs)))
    out['brightness_tied'] = bright
    out['brightness_tied_notes'] = dict(Rs_pc=round(rs, 1), strand_FWHM_pc=round(w, 3), Lb_geometric=round(Lb_geom, 1),
                                        face_on_Hb_per_n0V100_R=round(1e-6 * 1e7 * yieldHb_min, 3),
                                        face_on_5007_per_eps_n0_V100cubed_R=round(1e-6 * 1e7 * 0.5 * 1.4 * 1.6726e-24 * 1e14 / 3.967e-12, 1))
    return out


# ------------------------------------------------------------------ (d) O++ emissivity
E3 = np.array([0.0, 113.178, 306.174, 20273.27, 43185.74])
G3 = np.array([1, 3, 5, 5, 1])
HCK = 1.438777


def read_coll(fn, nlev=5):
    logT, om = None, {}
    for line in open(fn):
        if line.startswith('***') or not line.strip():
            continue
        p = line.split()
        i, j = int(p[0]), int(p[1])
        vals = np.array([float(x) for x in p[2:]])
        if i == 0 and j == 0:
            logT = vals
        elif i <= nlev and j <= nlev:
            om[(i - 1, j - 1)] = vals
    return logT, om


def read_A(fn, nlev=5):
    rows = []
    for line in open(fn):
        p = line.split()
        if not p or line.startswith('***') or p[0] == 'Aij' or p[0].startswith('1/s'):
            continue
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            pass
    return np.array([r[:nlev] for r in rows[:nlev]])


O3C = {'AK99': read_coll(RAW / 'o_iii_coll_AK99.dat.txt'), 'SSB14': read_coll(RAW / 'o_iii_coll_SSB14.dat.txt')}
O3A = read_A(RAW / 'o_iii_atom_FFT04-SZ00.dat.txt')


def o3_pop(ne, T, coll='AK99'):
    logT, om = O3C[coll]
    lt = np.clip(np.log10(T), logT[0], logT[-1])
    M = np.zeros((5, 5))
    for (i, j), vals in om.items():
        Om = np.interp(lt, logT, vals)
        qd = 8.629e-6 / (G3[j] * np.sqrt(T)) * Om
        qu = qd * G3[j] / G3[i] * np.exp(-(E3[j] - E3[i]) * HCK / T)
        M[j, i] += ne * qu; M[i, i] -= ne * qu
        M[i, j] += ne * qd; M[j, j] -= ne * qd
    for uu in range(5):
        for l in range(uu):
            M[l, uu] += O3A[uu, l]; M[uu, uu] -= O3A[uu, l]
    M[0, :] = 1.0
    rhs = np.zeros(5); rhs[0] = 1
    return np.linalg.solve(M, rhs)


def eps5007(T, ne=1.0, coll='AK99'):
    """photons s^-1 per (n_e n_O++) in 5007 (cm^3 s^-1)"""
    p = o3_pop(ne, T, coll)
    return p[3] * O3A[3, 2] / ne


def r4363(T, ne=1.0):
    p = o3_pop(ne, T)
    e5007 = p[3] * O3A[3, 2] * (E3[3] - E3[2])
    e4363 = p[4] * O3A[4, 3] * (E3[4] - E3[3])
    return e4363 / e5007


OH_SUN = 4.90e-4      # Asplund+2009, 12+log O/H = 8.69
NE_NH = 1.15          # n_e / n_H for H+ and partly doubly ionised He (range 1.1-1.2)


def oiii_selfconsistency():
    out = {'O_H': OH_SUN, 'ne_over_nH': NE_NH}
    Ts = [1e4, 1.5e4, 2e4, 3e4, 5e4]
    em = {}
    for T in Ts:
        e = eps5007(T)
        coef = 1e-6 * e * OH_SUN * PC            # R per (pc cm^-6 of n_e n_H) for x = 1
        aHa = 0.361 * (T / 1e4) ** -0.9           # R per (pc cm^-6 of n_e n_p)
        em[f'{T:.0f}'] = dict(eps5007_cm3s=float('%.4g' % e), eps_SSB14=(float('%.4g' % eps5007(T, coll='SSB14')) if T <= 2.5e4 else None),
                              I5007_per_EM_x1=round(coef, 3), OIII_over_Ha_photon_x1=round(coef / aHa, 2), OIII_over_Ha_photon_x0p3=round(0.3 * coef / aHa, 2),
                              r4363_5007=float('%.4g' % r4363(T)))
    out['emissivity'] = em
    # required EM and n_e at d = 280 pc
    d = 280.0
    w = d * 0.28 * DEG
    rs = d * 12.10 * DEG
    paths = {'filament_L=w': w, 'shell_tangent_L=2sqrt(2Rw)': 2 * np.sqrt(2 * rs * w - w * w), 'flat_sheet_L=50pc': 50.0}
    out['geometry'] = dict(distance_pc=d, strand_FWHM_pc=round(w, 3), shell_radius_pc=round(rs, 1), paths_pc={k: round(v, 2) for k, v in paths.items()})
    tab = {}
    for I in (2.0, 3.0, 4.0):
        for T in Ts:
            coef = em[f'{T:.0f}']['I5007_per_EM_x1']
            for x in (1.0, 0.3):
                emh = I / (coef * x)                 # pc cm^-6 (n_e n_H)
                row = dict(EM_neNH=round(emh, 3), EM_ne2=round(emh * NE_NH, 3))
                for pk, L in paths.items():
                    row[f'ne_{pk}'] = round(np.sqrt(emh * NE_NH / L), 3)
                row['Ha_R_from_same_gas'] = round(emh / NE_NH * 0.361 * (T / 1e4) ** -0.9 * NE_NH / NE_NH, 3)
                row['P_over_k_filament'] = float('%.3g' % (row['ne_filament_L=w'] * (1 + 1 / NE_NH + 0.1 / NE_NH) * T))
                row['P_over_k_shell'] = float('%.3g' % (row['ne_shell_tangent_L=2sqrt(2Rw)'] * (1 + 1 / NE_NH + 0.1 / NE_NH) * T))
                tab[f'I{I}_T{T:.0f}_x{x}'] = row
    out['table_280pc'] = tab
    # scaling with distance for one fiducial case
    out['ne_scaling'] = 'n_e proportional to d^-1/2 (w and R_s scale with d); I_5007 proportional to EM'
    # T needed for the SDSS-south [O III]/Halpha = 4.7 +- 1.0 (photon units assumed) for x = 1 and 0.3
    Tg = np.logspace(3.9, 4.9, 400)
    rat = np.array([1e-6 * eps5007(T) * OH_SUN * PC / (0.361 * (T / 1e4) ** -0.9) for T in Tg])
    out['T_for_OIII_Ha'] = {f'{r}_x{x}': round(float(np.interp(r / x, rat, Tg)), 0) for r in (3.7, 4.7, 5.7) for x in (1.0, 0.3)}
    r43 = np.array([r4363(T) for T in Tg])
    out['T_for_4363_5007'] = {str(r): round(float(np.interp(r, r43, Tg)), 0) for r in (0.022, 0.056, 0.077)}
    return out


# ------------------------------------------------------------------ (e) brightness conversion
LAM = 5006.84
C_A = 2.99792458e18
R_5007 = 1e6 / (4 * np.pi) * 6.62607e-27 * C_A / LAM / (180 / np.pi * 3600) ** 2     # erg s-1 cm-2 arcsec-2 per R


def mag_to_R(m, W=60.0, ab_offset=0.0):
    fnu = 10 ** (-0.4 * (m + ab_offset + 48.6))
    return fnu * C_A / LAM ** 2 * W / R_5007


def R_to_mag(I, W=60.0):
    flam = I * R_5007 / W
    return -2.5 * np.log10(flam * LAM ** 2 / C_A) - 48.6


def filt(lam, cwl, fwhm=60.0, N=3):
    return 1.0 / (1.0 + ((lam - cwl) / (fwhm / 2)) ** (2 * N))


def eff_curve(lam, thetas, weights, nstar, cwl0=5007.0, fwhm=60.0):
    T = np.zeros_like(lam)
    for th, wt in zip(thetas, weights):
        c = cwl0 * np.sqrt(1 - (np.sin(th) / nstar) ** 2)
        T += wt * filt(lam, c, fwhm)
    return T / np.sum(weights)


def conversion_cases():
    lam = np.arange(4700.0, 5300.0, 0.05)
    cases = {}
    # beam geometries
    r = np.linspace(0, 1, 401)[1:]
    geoms = {}
    for fno in (0.95, 2.8):
        tmax = np.arctan(1 / (2 * fno))
        geoms[f'converging_f/{fno}_on_axis'] = (np.arctan(r * np.tan(tmax)), r)
    for fa in (0, 5, 10, 13, 15):
        geoms[f'collimated_field_{fa}deg'] = (np.array([fa * DEG]), np.array([1.0]))
    for nstar in (1.8, 2.05):
        for cwl0 in (5007.0, 5015.0, 5030.0):
            for gname, (th, wt) in geoms.items():
                T = eff_curve(lam, th, wt, nstar, cwl0)
                integ = np.trapezoid(T, lam)
                t5007 = np.interp(5006.84, lam, T); t4959 = np.interp(4958.91, lam, T)
                S = t5007 + t4959 / 2.98
                true_over_paper = integ / (60.0 * S) if S > 0 else np.inf
                cases[f'n{nstar}_cwl{cwl0:.0f}_{gname}'] = dict(T5007=round(float(t5007), 3), T4959=round(float(t4959), 3),
                                                                 int_T=round(float(integ), 1), peak_lambda=round(float(lam[np.argmax(T)]), 1),
                                                                 true_over_paper=round(float(true_over_paper), 3))
    return cases


def isl():
    out = {}
    rad = 0.25
    area = np.pi * (rad * 3600) ** 2
    for k in ('P1', 'P2', 'P3', 'P4'):
        fn = ROOT / f'data/derived/verify_paper_raw/gaia_isl/isl_{k}_r{rad:.2f}.csv'
        if not fn.exists() or fn.stat().st_size < 30:
            out[k] = 'missing'
            continue
        a = np.genfromtxt(fn, delimiter=',', names=True)
        row = {}
        for glim in (14.0, 15.0, 16.0, 17.0, 18.0, 19.0):
            sel = a['gbin'] >= glim
            f = a['fsum'][sel].sum()
            row[f'SB_G_fainter_than_{glim:.0f}'] = round(float(-2.5 * np.log10(f / area)), 2)
        row['SB_all'] = round(float(-2.5 * np.log10(a['fsum'].sum() / area)), 2)
        row['n_stars'] = int(a['n'].sum())
        row['n_per_deg2_G17_21'] = round(float(a['n'][(a['gbin'] >= 17) & (a['gbin'] < 21)].sum() / (np.pi * rad ** 2)), 0)
        out[k] = row
    return out


def brightness():
    out = {}
    mags = {'P1': 23.60, 'P2': 23.87, 'P3': 23.85}
    out['AB_W60'] = {k: round(mag_to_R(m), 2) for k, m in mags.items()}
    out['W50_W70'] = {k: [round(mag_to_R(m, 50), 2), round(mag_to_R(m, 70), 2)] for k, m in mags.items()}
    out['Vega_offset_pm0.1'] = {k: [round(mag_to_R(m, ab_offset=0.1), 2), round(mag_to_R(m, ab_offset=-0.1), 2)] for k, m in mags.items()}
    out['with_4959_fully_transmitted'] = {k: round(mag_to_R(m) / (1 + 1 / 2.98), 2) for k, m in mags.items()}
    out['R_to_mag_AB_W60'] = {str(I): round(float(R_to_mag(I)), 2) for I in (1, 2, 3, 4, 10, 13)}
    out['R_per_cgs'] = R_5007
    for sky in (21.5, 22.0):
        out[f'fraction_of_sky_{sky}'] = {str(I): round(float(10 ** (-0.4 * (R_to_mag(I) - sky))), 4) for I in (2, 3, 4, 10, 13)}
        out[f'excess_10to13_minus_2to4_as_sky_fraction_{sky}'] = [round(float((10 ** (-0.4 * (R_to_mag(10) - sky)) - 10 ** (-0.4 * (R_to_mag(4) - sky)))), 4),
                                                                   round(float((10 ** (-0.4 * (R_to_mag(13) - sky)) - 10 ** (-0.4 * (R_to_mag(2) - sky)))), 4)]
    out['filter_cases'] = conversion_cases()
    out['gaia_isl'] = isl()
    # zero point sanity: e-/s for m_AB = 0 through a 50 mm f/0.95 lens (D = 52.6 mm) in a 60 A band
    D = 5.0 / 0.95
    area = np.pi * (D / 2) ** 2
    nph = 3631e-23 * C_A / LAM ** 2 * 60.0 / (6.62607e-27 * C_A / LAM) * area
    out['zeropoint_sanity'] = dict(aperture_cm2=round(area, 1), photons_per_s_mAB0=float('%.3g' % nph),
                                   ZP_eps_1e_per_s_QE0p7_T0p8=round(2.5 * np.log10(nph * 0.7 * 0.8), 2))
    return out


def make_fig(res_d, res_e):
    from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, FG, SUN, WARN, OK, ACCENT, DANGER
    import matplotlib.pyplot as plt
    fig, axs = new_sky_figure(13.5, 5.2, '[O III] self-consistency at 280 pc and the SWONS mag -> R conversion',
                              'own 5-level O++ (AK99 Ω) · solar O/H · strand FWHM 0.28° · filter: 6 nm, n* 2.05, 3-cavity shape', ncols=3)
    ax = axs[0]
    Tg = np.logspace(4, np.log10(5e4), 60)
    for I, ls in ((2.0, ':'), (4.0, '--'), (11.5, '-')):
        for x, col in ((1.0, OIII_C), (0.3, WARN)):
            emh = I / np.array([1e-6 * eps5007(T) * OH_SUN * PC * x for T in Tg])
            ne = np.sqrt(emh * NE_NH / res_d['geometry']['paths_pc']['shell_tangent_L=2sqrt(2Rw)'])
            ax.plot(Tg / 1e4, ne, color=col, ls=ls, lw=1.3, label=f'{I:g} R, x(O++)={x:g}')
    ax.axhspan(0.1, 3, color=OK, alpha=0.08, lw=0)
    ax.set_xscale('log'); ax.set_yscale('log')
    from matplotlib.ticker import NullFormatter
    ax.set_xticks([1, 2, 3, 5]); ax.set_xticklabels(['1', '2', '3', '5']); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0.07, 4)
    ax.set_xlabel('T_e (10$^4$ K)'); ax.set_ylabel('n_e (cm$^{-3}$), shell-tangent path 25 pc')
    ax.legend(fontsize=7, loc='upper right', ncol=2)
    label(ax, 1.05, 0.115, 'n_e 0.1–3 band', color=OK, size=8)
    ax.grid(alpha=0.3, which='both')
    ax = axs[1]
    lam = np.arange(4880.0, 5100.0, 0.1)
    r = np.linspace(0, 1, 401)[1:]
    for (lab, th, wt, col) in [('normal incidence', np.array([0.0]), np.array([1.0]), FG),
                               ('f/0.95 cone, on axis', np.arctan(r * np.tan(np.arctan(1 / 1.9))), r, DANGER),
                               ('f/2.8 cone', np.arctan(r * np.tan(np.arctan(1 / 5.6))), r, ACCENT),
                               ('collimated, 13° field', np.array([13 * DEG]), np.array([1.0]), WARN)]:
        ax.plot(lam, eff_curve(lam, th, wt, 2.05), color=col, lw=1.4, label=lab)
    ax.axvline(5006.84, color=OIII_C, lw=0.8); ax.axvline(4958.91, color=OIII_C, lw=0.8, ls=':')
    label(ax, 5009, 0.03, '5007', color=OIII_C, size=8); label(ax, 4961, 0.03, '4959', color=OIII_C, size=8)
    ax.set_ylim(0, 1.3); ax.set_xlabel('wavelength (Å)'); ax.set_ylabel('relative transmission')
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    ax.grid(alpha=0.3)
    ax = axs[2]
    Is = np.array([1, 2, 3, 4, 6, 10, 13])
    ax.plot(Is, R_to_mag(Is), color=FG, lw=1.2)
    ax.axhspan(23.6, 23.87, color=DANGER, alpha=0.2, lw=0)
    label(ax, 1.05, 23.75, 'SWONS Table 3', color=DANGER, size=8)
    ax.axvspan(2, 4, color=OK, alpha=0.12, lw=0)
    label(ax, 2.05, 25.95, 'DESI/SDSS/NSNS 2–4 R', color=OK, size=8)
    ax.set_xscale('log'); ax.invert_yaxis()
    ax.set_xticks([1, 2, 3, 4, 6, 10, 13]); ax.set_xticklabels(['1', '2', '3', '4', '6', '10', '13'])
    ax.set_xlabel('I(5007) (R)'); ax.set_ylabel('AB mag arcsec$^{-2}$ in 60 Å')
    ax2 = ax.twinx()
    fr = np.array([2, 3, 5, 10, 15])
    ticks = 21.5 - 2.5 * np.log10(fr / 100.0)
    ax2.set_yticks(ticks)
    ax2.set_yticklabels([f'{f:d}%' for f in fr])
    ax2.set_ylim(ax.get_ylim())
    ax2.set_ylabel('fraction of a 21.5 mag/arcsec$^2$ sky', color=FGDIM)
    ax.grid(alpha=0.3)
    finish(fig, str(ROOT / 'figures/verify_paper_oiii_conversion.png'))


def main():
    res = {'energetics': energetics(), 'oiii_selfconsistency': oiii_selfconsistency(), 'brightness_conversion': brightness()}
    (ROOT / 'data/derived/verify_paper_energetics.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res['energetics']['paper_repro'], indent=1))
    print(json.dumps(res['energetics']['paper_error_mc']))
    print(json.dumps(res['oiii_selfconsistency']['emissivity'], indent=1))
    print(json.dumps(res['oiii_selfconsistency']['T_for_OIII_Ha']), json.dumps(res['oiii_selfconsistency']['T_for_4363_5007']))
    print(json.dumps(res['oiii_selfconsistency']['geometry']))
    b = res['brightness_conversion']
    for k in ('AB_W60', 'W50_W70', 'Vega_offset_pm0.1', 'with_4959_fully_transmitted', 'R_to_mag_AB_W60', 'fraction_of_sky_21.5', 'fraction_of_sky_22.0',
              'excess_10to13_minus_2to4_as_sky_fraction_21.5', 'excess_10to13_minus_2to4_as_sky_fraction_22.0', 'gaia_isl', 'zeropoint_sanity'):
        print(k, json.dumps(b[k]))
    make_fig(res['oiii_selfconsistency'], b)


if __name__ == '__main__':
    main()
