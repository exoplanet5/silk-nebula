"""paper_review: shared physics helpers for scripts/paper_checks.py and scripts/paper_review_figures.py
(critical review of Wu et al. 2026, arXiv:2610.13329). Pure functions only; no file I/O.

Rotation curves are re-implemented from the kd package (T. V. Wenger, https://github.com/tvwenger/kd,
master branch read as plain text on 2026-11-16; kd code itself is never imported or executed).
"""
import numpy as np

C_KMS = 299792.458
U_STD, V_STDc, W_STD = 10.27, 15.32, 7.74   # IAU standard solar motion components used by kd (km/s)


def unit_vec(l, b):
    l, b = np.deg2rad(l), np.deg2rad(b)
    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])


# ---------------------------------------------------------------- rotation curves / kinematic distance
def theta_urc14(Rg, a1=241.0, a2=0.90, a3=1.46, R0=8.34):
    """Persic+1996 URC as implemented in kd/reid14_rotcurve.py (Reid+2014 parameters)."""
    x = Rg / (a2 * R0)
    beta = 0.72 + 0.44 * np.log10((a3 / 1.5) ** 5)
    vd2 = beta * 1.97 * x ** 1.22 / (x ** 2 + 0.78 ** 2) ** 1.43
    vh2 = (1 - beta) * (1 + a3 ** 2) * x ** 2 / (x ** 2 + a3 ** 2)
    return a1 * np.sqrt(np.clip(vd2 + vh2, 0, None))


def theta_urc19(Rg, a2=0.977, a3=1.623, R0=8.166):
    """URC as implemented in kd/reid19_rotcurve.py (Reid+2019 FORTRAN reproduction)."""
    rho = Rg / (a2 * R0)
    lam = (a3 / 1.5) ** 5
    ll = np.log10(lam)
    t1 = 200 * lam ** 0.41
    t2 = np.sqrt(0.8 + 0.49 * ll + 0.75 * np.exp(-0.4 * lam) / (0.47 + 2.25 * lam ** 0.4))
    t3 = (0.72 + 0.44 * ll) * 1.97 * rho ** 1.22 / (rho ** 2 + 0.61) ** 1.43
    t4 = 1.6 * np.exp(-0.4 * lam) * rho ** 2 / (rho ** 2 + 2.25 * lam ** 0.4)
    return t1 / t2 * np.sqrt(np.clip(t3 + t4, 0, None))


def vlsr_general(l, b, d, R0, theta_fn, sun_uvw, theta0=None):
    """IAU-standard V_LSR of gas on a circular orbit at heliocentric distance d (kpc); full 3-D
    projection (includes cos b). theta_fn(R): rotation curve; sun_uvw: solar peculiar motion
    relative to the circular speed at R0 assumed by that model."""
    d = np.asarray(d, float)
    cl, sl = np.cos(np.deg2rad(l)), np.sin(np.deg2rad(l))
    cb, sb = np.cos(np.deg2rad(b)), np.sin(np.deg2rad(b))
    X = d * cb * cl - R0
    Y = d * cb * sl
    Rg = np.maximum(np.hypot(X, Y), 1e-6)
    th = theta_fn(Rg)
    th0 = theta_fn(np.array(R0)) if theta0 is None else theta0
    vx, vy = th * Y / Rg, -th * X / Rg
    U, V, W = sun_uvw
    vhel = (vx - U) * cb * cl + (vy - th0 - V) * cb * sl + (-W) * sb
    return vhel + (U_STD * cb * cl + V_STDc * cb * sl + W_STD * sb)


def vlsr_kd14(l, b, d, a1=241.0, a2=0.90, a3=1.46, R0=8.34):
    """Exact emulation of kd reid14 calc_vlsr: V = R0 sin l (Theta/R - Theta0/R0) (no cos b)."""
    dcb = np.asarray(d, float) * np.cos(np.deg2rad(b))
    Rg = np.sqrt(R0 ** 2 + dcb ** 2 - 2 * R0 * dcb * np.cos(np.deg2rad(l)))
    Rg = np.maximum(Rg, 1e-6)
    return R0 * np.sin(np.deg2rad(l)) * (theta_urc14(Rg, a1, a2, a3, R0) / Rg
                                          - theta_urc14(np.array(R0), a1, a2, a3, R0) / R0)


MODELS = {
    'flat_R0_8.15_T236_stdsun': lambda l, b, d: vlsr_general(l, b, d, 8.15, lambda r: 236.0 + 0 * r, (U_STD, V_STDc, W_STD), 236.0),
    'flat_R0_8.34_T240_stdsun': lambda l, b, d: vlsr_general(l, b, d, 8.34, lambda r: 240.0 + 0 * r, (U_STD, V_STDc, W_STD), 240.0),
    'reid14_A5_linear': lambda l, b, d: vlsr_general(l, b, d, 8.34, lambda r: 240.0 - 0.2 * (r - 8.34), (10.7, 15.6, 8.9), 240.0),
    'reid19_A5_linear': lambda l, b, d: vlsr_general(l, b, d, 8.15, lambda r: 236.0 - 0.1 * (r - 8.15), (10.6, 10.7, 7.6), 236.0),
    'kd_reid14_URC': lambda l, b, d: vlsr_kd14(l, b, d),
    'kd_reid19_URC': lambda l, b, d: vlsr_general(l, b, d, 8.166, lambda r: theta_urc19(r), (10.449, 12.092, 7.729)),
}


# ---------------------------------------------------------------- SNR energetics
def vs_dopita(ne, n0):
    """Dopita (1979): n_e[S II] ~ 45 n0 (V_s/100 km/s)^2  ->  V_s (km/s)."""
    return 100.0 * np.sqrt(ne / (45.0 * n0))


def e_boumis(n0, vs, rs):
    """Cioffi+1988 PDS energy in the form of Boumis+2022 (erg); n0 cm^-3, vs km/s, rs pc."""
    return 6.8e43 * n0 ** 1.16 * vs ** 1.35 * rs ** 3.16


def r_pds(e51, n0):
    return 14.0 * e51 ** (2 / 7) * n0 ** (-3 / 7)


def age_cioffi(e51, n0, rs, rpds, plus=1 / 3):
    """t = (3/4) t_PDS [(R/R_PDS)^(10/3) + 1/3]; the paper prints '+1'."""
    return 9.98e3 * e51 ** (3 / 14) * n0 ** (-4 / 7) * ((rs / rpds) ** (10 / 3) + plus)


# ---------------------------------------------------------------- [S II] density
SA, SB_, SC = 0.4315, 2107.0, 627.1   # Sanders et al. 2016, ApJ 816, 23 ([S II], T = 1e4 K)


def ne_proxauf(Rr):
    """Proxauf+2014 fit as printed in the paper (eq. 5), T = 1e4 K."""
    return 10 ** (0.0543 * np.tan(-3.0553 * Rr + 2.8506) + 6.98 - 10.6905 * Rr + 9.9186 * Rr ** 2 - 3.5442 * Rr ** 3)


def ne_sanders(Rr):
    return (SC * Rr - SA * SB_) / (SA - Rr)


def R_sanders(ne):
    return SA * (SB_ + ne) / (SC + ne)


# ---------------------------------------------------------------- emission measures
OMEGA_1D, G_3P, T_EXC, BR5007, OH = 2.29, 9.0, 28870.0, 0.75, 4.9e-4


def q_oiii(T):
    """[O III] 3P -> 1D2 collisional excitation rate (cm^3 s^-1), Omega=2.29 (Osterbrock & Ferland 2006)."""
    return 8.629e-6 / (np.sqrt(T) * G_3P) * OMEGA_1D * np.exp(-T_EXC / T)


def em_from_oiii(I_R, T, x_opp=1.0, oh=OH):
    """Minimum EM (pc cm^-6) = int n_e n_p dl needed for I_R Rayleigh of [O III] 5007 (low-density limit)."""
    col = 1e6 * I_R / (q_oiii(T) * BR5007)        # int n_e n(O++) dl  [cm^-5]
    return col / (x_opp * oh) / 3.0857e18


def em_from_ha(I_R, T):
    """EM = 2.75 T4^0.9 I_Halpha (Haffner+2003)."""
    return 2.75 * (T / 1e4) ** 0.9 * I_R


def oiii_ha_photon_ratio(T, x_opp=1.0, oh=OH):
    """[O III]5007/H-alpha in Rayleigh (photon) units for gas with O++/O = x at temperature T."""
    return 3.0857e18 * 1e-6 * q_oiii(T) * BR5007 * x_opp * oh * 2.75 * (T / 1e4) ** 0.9


# ---------------------------------------------------------------- surface brightness
def mag_to_R(m, lam=5006.84, width=60.0):
    """AB mag/arcsec^2 of pure line emission spread over a filter of FWHM `width` (A) -> Rayleigh.
    Returns (I_R, erg/s/cm2/arcsec2 per Rayleigh)."""
    fnu = 10 ** (-0.4 * (m + 48.60))                    # erg/s/cm2/Hz/arcsec2 (AB)
    flam = fnu * C_KMS * 1e13 / lam ** 2                # erg/s/cm2/A/arcsec2
    F = flam * width                                     # erg/s/cm2/arcsec2
    e_ph = 6.62607015e-27 * 2.99792458e10 / (lam * 1e-8)
    one_R = 1e6 / (4 * np.pi) / 4.25451703e10 * e_ph     # erg/s/cm2/arcsec2 per Rayleigh
    return F / one_R, one_R
