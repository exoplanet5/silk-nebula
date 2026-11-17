"""verify_altmodel_atomic: small self-contained atomic-physics toolkit (no pyneb / ChiantiPy).

Atomic data (all parsed as text, never executed):
  * CHIANTI v7.1.4 (Dere+1997 A&AS 125,149; Landi+2013 ApJ 763,86): .elvlc energies/statistical weights,
    .wgfa A-values, .splups Burgess-Tully-scaled Maxwellian collision strengths (electron impact),
    .rrparams (Badnell / Verner-Ferland RR fits), .drparams (Badnell DR fits), chianti.ioneq (CIE check),
    continuum/verner_short.txt (Verner & Yakovlev 1995 / Verner+1996 outer-shell photoionisation fits).
  * Dere (2007, A&A 466, 771) Table 29 (VizieR J/A+A/466/771): spline fits to the scaled direct
    (+ excitation-autoionisation) collisional ionisation rate coefficients.
  * H I recombination lines: case-B fits of Draine (2011, "Physics of the ISM", eqs. 14.6-14.8),
    alpha_B = 2.54e-13 T4^(-0.8163-0.0208 ln T4); 4 pi j(Hbeta) = 1.24e-25 T4^(-0.942-0.031 ln T4) n_e n_p.
    Halpha/Hbeta = 2.86 T4^-0.07 (approximation to Hummer & Storey 1987 case B, +-3% for 0.5<T4<2).
  * He II 4686 recombination: j(4686)/j(Hbeta) = 12 n(He++)/n(H+) at 1e4 K (ratio of the Osterbrock &
    Ferland 2006 Table 4.5 / 4.2 emissivities, ~ +-20%); used only for upper-limit comparisons.
Collisional excitation by protons is ignored (affects only fine-structure IR lines).
Units: energies erg, T in K, rate coefficients cm^3 s^-1.
"""
import os
import functools
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import exp1

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CH = os.path.join(ROOT, "data/raw/verify_altmodel/chianti/dbase")
DERE = os.path.join(ROOT, "data/raw/verify_altmodel/dere2007/table29.tsv")

K_B = 1.380649e-16
H_PL = 6.62607015e-27
C_L = 2.99792458e10
EV = 1.602176634e-12
RYD_K = 157887.51  # 1 Rydberg / k_B in K
CM1_K = 1.4387769  # hc/k in K cm
SYMB = {1: "h", 2: "he", 6: "c", 7: "n", 8: "o", 10: "ne", 12: "mg", 14: "si", 16: "s", 18: "ar"}
ZOF = {v: k for k, v in SYMB.items()}

# Solar photospheric abundances (Asplund et al. 2009, ARA&A 47, 481), by number relative to H.
ABUND_A09 = {1: 1.0, 2: 0.0851, 6: 2.69e-4, 7: 6.76e-5, 8: 4.90e-4, 10: 8.51e-5, 12: 3.98e-5, 14: 3.24e-5,
             16: 1.32e-5, 18: 2.51e-6}


def _iname(Z, ion):
    return f"{SYMB[Z]}_{ion}"


def _path(Z, ion, ext):
    n = _iname(Z, ion)
    return os.path.join(CH, SYMB[Z], n, f"{n}.{ext}")


# ----------------------------------------------------------------------------------------------
# CHIANTI readers
# ----------------------------------------------------------------------------------------------
def read_elvlc(Z, ion):
    E, g, lab = [], [], []
    with open(_path(Z, ion, "elvlc")) as f:
        for line in f:
            p = line.split()
            if len(p) < 6 or p[0] == "-1":
                break
            try:
                int(p[0])
            except ValueError:
                break
            # columns: idx conf config 2S+1 L Lsym J mult Eobs(cm-1) Eobs(Ryd) Eth(cm-1) Eth(Ryd)
            vals = p[-4:]
            mult = int(p[-5])
            eobs, eth = float(vals[0]), float(vals[2])
            e = eobs if (eobs > 0 or int(p[0]) == 1) else eth
            E.append(e)
            g.append(mult)
            lab.append(" ".join(p[2:-5]))
    return np.array(E), np.array(g, float), lab


def read_wgfa(Z, ion, nlev):
    A = np.zeros((nlev, nlev))  # A[u,l]
    wl = np.zeros((nlev, nlev))
    with open(_path(Z, ion, "wgfa")) as f:
        for line in f:
            p = line.split()
            if not p or p[0] == "-1":
                break
            try:
                l1, l2 = int(p[0]), int(p[1])
                w, a = float(p[2]), float(p[4])
            except (ValueError, IndexError):
                break
            lo, up = min(l1, l2) - 1, max(l1, l2) - 1
            if up >= nlev:
                continue
            A[up, lo] += a
            wl[up, lo] = abs(w)
    return A, wl


def read_splups(Z, ion):
    out = []
    with open(_path(Z, ion, "splups")) as f:
        for line in f:
            if line.strip().startswith("-1"):
                break
            if len(line) < 40:
                continue
            ints = [int(line[3 * k:3 * k + 3]) for k in range(5)]
            rest = line[15:].rstrip("\n")
            nums = [float(rest[10 * k:10 * k + 10]) for k in range(len(rest) // 10) if rest[10 * k:10 * k + 10].strip()]
            gf, de, cc = nums[0], nums[1], nums[2]
            spl = np.array(nums[3:])
            out.append((ints[2] - 1, ints[3] - 1, ints[4], gf, de, cc, spl))
    return out


def _bt_descale(T, ttype, de_ryd, cc, spl):
    kte = (T / RYD_K) / de_ryd
    n = len(spl)
    xs = np.linspace(0, 1, n)
    cs = CubicSpline(xs, spl, bc_type="natural")
    if ttype in (1, 4):
        x = 1.0 - np.log(cc) / np.log(kte + cc)
    else:
        x = kte / (kte + cc)
    y = cs(np.clip(x, 0, 1))
    if ttype == 1:
        return y * np.log(kte + np.e)
    if ttype == 2:
        return y
    if ttype == 3:
        return y / (kte + 1.0)
    if ttype == 4:
        return y * np.log(kte + cc)
    if ttype == 5:
        return y / kte
    if ttype == 6:
        return 10.0 ** y
    raise ValueError(ttype)


class Ion:
    """Multi-level ion with CHIANTI electron-impact collision strengths."""

    def __init__(self, Z, ion, nmax=None):
        self.Z, self.ion = Z, ion
        E, g, lab = read_elvlc(Z, ion)
        if nmax:
            E, g, lab = E[:nmax], g[:nmax], lab[:nmax]
        self.n = len(E)
        self.E = E  # cm^-1
        self.g = g
        self.lab = lab
        self.A, self.wl = read_wgfa(Z, ion, self.n)
        self.trans = [t for t in read_splups(Z, ion) if t[1] < self.n and t[0] < self.n and t[2] <= 6]

    def rate_matrix(self, T, ne):
        n = self.n
        R = np.array(self.A, copy=True)  # R[u,l] radiative; we build total rate matrix M[i,j] = rate i->j
        M = np.zeros((n, n))
        M += R  # downward radiative (u->l)
        for (lo, up, tt, gf, de, cc, spl) in self.trans:
            if tt > 5:  # type 6 = proton rates in some files; skip
                continue
            ups = _bt_descale(T, tt, de, cc, spl)
            dE_K = (self.E[up] - self.E[lo]) * CM1_K
            if dE_K <= 0:
                dE_K = de * RYD_K
            qdn = 8.629e-6 * ups / (self.g[up] * np.sqrt(T))
            qup = qdn * self.g[up] / self.g[lo] * np.exp(-dE_K / T)
            M[lo, up] += ne * qup
            M[up, lo] += ne * qdn
        return M

    def populations(self, T, ne):
        M = self.rate_matrix(T, ne)
        n = self.n
        M[1:, 0] += 1e-20  # negligible leak to ground: keeps levels with no CHIANTI rates (e.g. Mg II) non-singular
        L = M.T - np.diag(M.sum(axis=1))  # dN_j/dt = sum_i N_i M[i,j] - N_j sum_k M[j,k]
        L[0, :] = 1.0
        b = np.zeros(n)
        b[0] = 1.0
        return np.linalg.solve(L, b)

    def line_emis(self, T, ne, up, lo):
        """energy emitted per ion per unit n_e: eps such that j = eps * n_e * n_ion (erg cm^3 s^-1)."""
        p = self.populations(T, ne)
        nu = (self.E[up] - self.E[lo]) * C_L
        return p[up] * self.A[up, lo] * H_PL * nu / ne

    def all_line_power(self, T, ne):
        """total radiated power per ion per n_e (erg cm^3 s^-1), all transitions."""
        p = self.populations(T, ne)
        dE = (self.E[:, None] - self.E[None, :]) * C_L * H_PL
        P = (p[:, None] * self.A * np.clip(dE, 0, None)).sum()
        return P / ne


@functools.lru_cache(maxsize=None)
def get_ion(Z, ion, nmax=None):
    return Ion(Z, ion, nmax)


def find_level(ionobj, lab_sub, J=None):
    for i, l in enumerate(ionobj.lab):
        if lab_sub in l:
            return i
    raise KeyError(lab_sub)


# Named forbidden / UV lines: (Z, ion, upper idx, lower idx) using CHIANTI v7 level numbering (0-based)
LINES = {
    "OIII5007": (8, 3, 3, 2), "OIII4959": (8, 3, 3, 1), "OIII4363": (8, 3, 4, 3), "OIII2321": (8, 3, 4, 1),
    "OIII1666": (8, 3, 5, 2), "OIII1661": (8, 3, 5, 1),
    "OII3729": (8, 2, 1, 0), "OII3726": (8, 2, 2, 0),  # 2D5/2 -> 4S (3729), 2D3/2 -> 4S (3726)
    "NII6583": (7, 2, 3, 2), "NII6548": (7, 2, 3, 1), "NII5755": (7, 2, 4, 3),
    "SII6716": (16, 2, 2, 0), "SII6731": (16, 2, 1, 0),
    "NeIII3869": (10, 3, 3, 0), "NeIII3968": (10, 3, 3, 1),
    "SIII9069": (16, 3, 3, 1), "SIII9531": (16, 3, 3, 2),
    "ArIII7136": (18, 3, 3, 0),
    "CIV1548": (6, 4, 2, 0), "CIV1551": (6, 4, 1, 0),
    "OVI1032": (8, 6, 2, 0), "OVI1038": (8, 6, 1, 0),
    "NV1239": (7, 5, 2, 0), "NV1243": (7, 5, 1, 0),
    "CIII1909": (6, 3, 2, 0), "CIII977": None,
    "HeII1640": None,
}


def check_line_levels():
    """print level labels / wavelengths for each named line (sanity check of indices)."""
    rows = []
    for k, v in LINES.items():
        if v is None:
            continue
        Z, ion, up, lo = v
        I = get_ion(Z, ion)
        rows.append((k, I.lab[up], I.lab[lo], I.wl[up, lo], I.A[up, lo]))
    return rows


def emis(name, T, ne):
    Z, ion, up, lo = LINES[name]
    return get_ion(Z, ion).line_emis(T, ne, up, lo)


def oii_7325_emis(T, ne):
    """[O II] 7320+7330 (2P -> 2D) per O+ per n_e."""
    I = get_ion(8, 2)
    tot = 0.0
    p = I.populations(T, ne)
    for up in (3, 4):
        for lo in (1, 2):
            if I.A[up, lo] > 0:
                tot += p[up] * I.A[up, lo] * H_PL * C_L * (I.E[up] - I.E[lo])
    return tot / ne


def oii_3727_emis(T, ne):
    return emis("OII3726", T, ne) + emis("OII3729", T, ne)


# ----------------------------------------------------------------------------------------------
# Hydrogen / helium recombination
# ----------------------------------------------------------------------------------------------
def alpha_B(T):
    T4 = np.asarray(T) / 1e4
    return 2.54e-13 * T4 ** (-0.8163 - 0.0208 * np.log(T4))


def alpha_A_H(T):
    """case-A RR of H (Hui & Gnedin 1997 fit)."""
    lam = 2 * 157807.0 / np.asarray(T)
    return 1.269e-13 * lam ** 1.503 / (1 + (lam / 0.522) ** 0.47) ** 1.923


def j_hbeta(T):
    """4 pi j(Hbeta)/(n_e n_p), erg cm^3 s^-1 (case B; Draine 2011 eq 14.8)."""
    T4 = np.asarray(T) / 1e4
    return 1.24e-25 * T4 ** (-0.942 - 0.031 * np.log(T4))


def j_halpha(T):
    T4 = np.asarray(T) / 1e4
    return j_hbeta(T) * 2.86 * T4 ** (-0.07)


def j_heii4686(T):
    """4 pi j(4686)/(n_e n_He++) ~ 12 j(Hbeta) (see module doc; T-dependence of the ratio ignored)."""
    return 12.0 * j_hbeta(T)


def h_collisional_balmer(T, ne):
    """Collisionally excited Hbeta, Halpha, Lyalpha per H0 per n_e (erg cm^3 s^-1) from the CHIANTI H I
    model (levels up to n=5, cascades included in the level solution; optically thin = case A-like)."""
    I = get_ion(1, 1)
    p = I.populations(T, ne)
    out = {}
    for name, nu_, nl_ in (("Hb", 4, 2), ("Ha", 3, 2), ("Lya", 2, 1)):
        tot = 0.0
        for u in range(I.n):
            if not I.lab[u].startswith(f"{nu_}"):
                continue
            for l in range(I.n):
                if I.lab[l].startswith(f"{nl_}") and I.A[u, l] > 0:
                    tot += p[u] * I.A[u, l] * H_PL * C_L * (I.E[u] - I.E[l])
        out[name] = tot / ne
    return out


# ----------------------------------------------------------------------------------------------
# Ionisation / recombination rates
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=1)
def _dere_table():
    rows = {}
    with open(DERE) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#") and l.strip()]
    for l in lines:
        p = l.split("\t")
        if len(p) < 30:
            continue
        try:
            Z, ion, nspl = int(p[1]), int(p[2]), int(p[3])
        except ValueError:
            continue
        I, tmin = float(p[4]), float(p[5])
        xs = [float(v) for v in p[6:26] if v.strip()][:nspl]
        rho = [float(v) for v in p[26:46] if v.strip()][:nspl]
        rows[(Z, ion)] = (I, tmin, np.array(xs), np.array(rho))
    return rows


def ci_rate(Z, ion, T):
    """Collisional ionisation rate coefficient of ion stage `ion` (1 = neutral), Dere (2007)."""
    tab = _dere_table()
    if (Z, ion) not in tab:
        return np.zeros_like(np.asarray(T, float))
    I, tmin, xs, rho = tab[(Z, ion)]
    T = np.asarray(T, float)
    t = K_B * T / (I * EV)
    x = 1.0 - np.log(2.0) / np.log(t + 2.0)
    cs = CubicSpline(xs, rho, bc_type="natural")
    r = cs(np.clip(x, xs[0], xs[-1]))
    R = 1e-6 * r * exp1(1.0 / t) / (np.sqrt(t) * I ** 1.5)
    return np.where(T >= 0.3 * tmin, R, R)  # (extrapolation below Tmin handled by E1 Arrhenius factor)


def _read_params(Z, ion, ext):
    p = _path(Z, ion, ext)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        L = [l for l in f]
    return L


@functools.lru_cache(maxsize=None)
def _rr_par(Z, ion):
    L = _read_params(Z, ion, "rrparams")
    if L is None:
        return None
    typ = int(L[0].split()[0])
    vals = [float(v) for v in L[1].split()]
    return typ, vals


@functools.lru_cache(maxsize=None)
def _dr_par(Z, ion):
    L = _read_params(Z, ion, "drparams")
    if L is None:
        return None
    typ = int(L[0].split()[0])
    E = [float(v) for v in L[1].split()[2:]]
    c = [float(v) for v in L[2].split()[2:]]
    return typ, np.array(E), np.array(c)


def rr_rate(Z, ion, T):
    """Radiative recombination coefficient of ion stage `ion` (recombining ion; ion>=2)."""
    T = np.asarray(T, float)
    par = _rr_par(Z, ion)
    if par is None:
        return np.zeros_like(T)
    typ, v = par
    if typ in (1, 2):
        A, B, T0, T1 = v[3], v[4], v[5], v[6]
        if typ == 2:
            Cc, T2 = v[7], v[8]
            B = B + Cc * np.exp(-T2 / T)
        s0, s1 = np.sqrt(T / T0), np.sqrt(T / T1)
        return A / (s0 * (1 + s0) ** (1 - B) * (1 + s1) ** (1 + B))
    if typ == 3:
        A, eta = v[2], v[3]
        return A * (T / 1e4) ** (-eta)
    raise ValueError(typ)


def dr_rate(Z, ion, T):
    T = np.asarray(T, float)
    par = _dr_par(Z, ion)
    if par is None:
        return np.zeros_like(T)
    typ, E, c = par
    if typ == 1:
        return T ** -1.5 * np.sum(c[:, None] * np.exp(-E[:, None] / T[None, :]) if T.ndim else c * np.exp(-E / T), axis=0)
    if typ == 2:  # Mazzotta-type: c E^-? (not used for our ions)
        return T ** -1.5 * np.sum(c[:, None] * np.exp(-E[:, None] / T[None, :]), axis=0)
    raise ValueError(typ)


def rec_rate(Z, ion, T):
    """total (RR+DR) recombination coefficient of stage `ion` -> ion-1."""
    T = np.atleast_1d(np.asarray(T, float))
    if Z == 1 and ion == 2:
        return alpha_A_H(T)
    return rr_rate(Z, ion, T) + dr_rate(Z, ion, T)


def cie_fractions(Z, T):
    T = np.atleast_1d(np.asarray(T, float))
    nst = Z + 1
    f = np.zeros((len(T), nst))
    for k, t in enumerate(T):
        logx = np.zeros(nst)
        for i in range(1, nst):
            C = ci_rate(Z, i, t)
            a = rec_rate(Z, i + 1, t)[0]
            logx[i] = logx[i - 1] + np.log(max(C, 1e-300)) - np.log(max(a, 1e-300))
        logx -= logx.max()
        x = np.exp(logx)
        f[k] = x / x.sum()
    return f


@functools.lru_cache(maxsize=1)
def chianti_ioneq():
    p = os.path.join(CH, "ioneq/chianti.ioneq")
    with open(p) as f:
        L = f.readlines()
    nT, nZ = [int(v) for v in L[0].split()[:2]]
    logT = np.array([float(v) for v in L[1].split()[:nT]])
    out = {}
    for l in L[2:]:
        if l.strip().startswith("-1"):
            break
        Z, ion = int(l[0:3]), int(l[3:6])
        rest = l[6:].rstrip("\n")
        out[(Z, ion)] = np.array([float(rest[10 * k:10 * k + 10]) for k in range(nT)])
    return logT, out


# ----------------------------------------------------------------------------------------------
# Photoionisation cross sections (Verner & Yakovlev 1995 outer-shell form, CHIANTI verner_short.txt)
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=1)
def _verner():
    d = {}
    with open(os.path.join(CH, "continuum/verner_short.txt")) as f:
        for l in f:
            p = l.split()
            if len(p) < 10:
                continue
            Z, N, n, ll = [int(v) for v in p[:4]]
            eth, e0, s0, ya, P, yw = [float(v) for v in p[4:10]]
            d[(Z, Z - N + 1)] = (eth, e0, s0, ya, P, yw, ll)
    return d


def phot_xs(Z, ion, E_eV):
    """photoionisation cross-section (cm^2) of stage `ion` (1 = neutral) at photon energy E (eV)."""
    E = np.asarray(E_eV, float)
    v = _verner().get((Z, ion))
    if v is None:
        return np.zeros_like(E)
    eth, e0, s0, ya, P, yw, ll = v
    y = E / e0
    Q = 5.5 + ll - 0.5 * P
    F = ((y - 1) ** 2 + yw ** 2) * y ** (-Q) * (1 + np.sqrt(y / ya)) ** (-P)
    return np.where(E >= eth, s0 * F * 1e-18, 0.0)


def ion_potential(Z, ion):
    return _verner()[(Z, ion)][0]
