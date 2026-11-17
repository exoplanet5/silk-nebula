"""w3_oii_common.py -- shared helpers for the wave-3 [O II] doublet forensic test (tag oii).

Own code. Line profiles, doublet fitting (linear amplitudes for fixed v/sigma, profile likelihood over v, sigma),
bootstrap-covariance chi^2, and the PyNeb-format 5-level atoms (O+, S+) used for the density/temperature curves.
Wavelengths are vacuum. Doublet ratios are ENERGY ratios (fluxes in erg) unless stated; for 3729/3726 and
6716/6731 photon and energy ratios differ by < 0.1 % and 0.2 %, respectively.
"""
import os
import numpy as np
from scipy.special import erf

ROOT = '/Users/mickey/astronomy/silk-nebula'
C = 299792.458


def air2vac(lam):
    """Edlen 1966 / IAU (Morton 1991) standard air -> vacuum (A), iterated."""
    lam = np.asarray(lam, float)
    vac = lam.copy()
    for _ in range(4):
        s2 = (1e4 / vac) ** 2
        n = 1 + 6.4328e-5 + 2.94981e-2 / (146 - s2) + 2.5540e-4 / (41 - s2)
        vac = lam * n
    return vac


# NIST / Osterbrock & Ferland air wavelengths -> vacuum
L_AIR = {'OII3726': 3726.032, 'OII3729': 3728.815, 'OIII4959': 4958.911, 'OIII5007': 5006.843,
         'SII6716': 6716.440, 'SII6731': 6730.810, 'OII7319': 7318.92, 'OII7320': 7319.99, 'OII7330': 7329.66,
         'OII7331': 7330.73, 'Ha': 6562.80, 'NII6583': 6583.45}
L = {k: float(air2vac(v)) for k, v in L_AIR.items()}


def gauss(x, mu, s):
    """Unit-integral Gaussian sampled at x (point sampling)."""
    return np.exp(-0.5 * ((x - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)


def gauss_binint(x, mu, s, px):
    """Unit-integral Gaussian averaged over pixels of width px centred on x."""
    a = (x - 0.5 * px - mu) / (np.sqrt(2) * s)
    b = (x + 0.5 * px - mu) / (np.sqrt(2) * s)
    return 0.5 * (erf(b) - erf(a)) / px


def profile(x, mu, s, mode='point', px=0.8):
    if mode == 'point':
        return gauss(x, mu, s)
    if mode == 'binint':
        return gauss_binint(x, mu, s, px)
    raise ValueError(mode)


def design(lam, lines, v, sig, npoly=2, mode='point', px=0.8, lam0=None):
    """Columns: one unit-integral profile per line (list of rest vac wavelengths), then a polynomial baseline.
    sig: scalar or list (per line)."""
    sig = np.broadcast_to(np.atleast_1d(np.asarray(sig, float)), (len(lines),))
    cols = [profile(lam, l0 * (1 + v / C), s, mode, px) for l0, s in zip(lines, sig)]
    lam0 = lam.mean() if lam0 is None else lam0
    for p in range(npoly):
        cols.append(((lam - lam0) / 30.0) ** p)
    return np.vstack(cols)


def lin_solve(X, y, Winv=None, w=None):
    """Generalized least squares: Winv = inverse covariance (n x n) or w = diagonal weights."""
    if Winv is None:
        Winv = np.diag(w if w is not None else np.ones(y.size))
    A = X @ Winv @ X.T
    b = X @ Winv @ y
    p = np.linalg.solve(A, b)
    r = y - p @ X
    return p, float(r @ Winv @ r), np.linalg.inv(A)


def fit_doublet(lam, y, lines, v, sig, Winv=None, w=None, npoly=2, mode='point', px=0.8):
    X = design(lam, lines, v, sig, npoly, mode, px)
    p, chi2, cov = lin_solve(X, y, Winv, w)
    return p, chi2, cov


def fit_doublet_ratio_fixed(lam, y, lines, v, sig, ratio, Winv=None, w=None, npoly=2, mode='point', px=0.8):
    """Doublet with amplitude ratio fixed (line2/line1 = ratio): one amplitude + baseline."""
    sig = np.broadcast_to(np.atleast_1d(np.asarray(sig, float)), (len(lines),))
    g = profile(lam, lines[0] * (1 + v / C), sig[0], mode, px) + ratio * profile(lam, lines[1] * (1 + v / C), sig[1],
                                                                                    mode, px)
    X = np.vstack([g] + [((lam - lam.mean()) / 30.0) ** p for p in range(npoly)])
    p, chi2, cov = lin_solve(X, y, Winv, w)
    return p, chi2


def shrink_cov(bs, lam_sel=None, alpha=0.2):
    """Bootstrap covariance with shrinkage towards its diagonal (bs: [nb, npix])."""
    S = np.cov(bs, rowvar=False)
    D = np.diag(np.diag(S))
    return (1 - alpha) * S + alpha * D


def scan_best(lam, y, lines, Winv, vgrid, sgrid, npoly=2, mode='point', px=0.8, ratio=None, sig_ratio=None):
    """Grid search over (v, sigma); free amplitudes (or fixed ratio). sig_ratio: sigma(line2)/sigma(line1).
    Returns best (v, s), chi2 grid, amplitude array at best."""
    chi = np.full((len(vgrid), len(sgrid)), np.inf)
    best = (np.inf, None, None, None)
    for i, v in enumerate(vgrid):
        for j, s in enumerate(sgrid):
            sg = [s, s * (sig_ratio if sig_ratio else 1.0)]
            if ratio is None:
                p, c2, _ = fit_doublet(lam, y, lines, v, sg, Winv=Winv, npoly=npoly, mode=mode, px=px)
            else:
                p, c2 = fit_doublet_ratio_fixed(lam, y, lines, v, sg, ratio, Winv=Winv, npoly=npoly, mode=mode,
                                                px=px)
            chi[i, j] = c2
            if c2 < best[0]:
                best = (c2, v, s, p)
    return best, chi


# ------------------------------------------------------------------ PyNeb-format n-level atoms
PDIR = os.path.join(ROOT, 'data/raw/verify_shockmodel_atomic/pyneb')
H_PL, C_L = 6.62607015e-27, 2.99792458e10


def _eval_j(s):
    if '/' in s:
        a, b = s.split('/')
        return float(a) / float(b)
    return float(s)


def _read_levels(fn, n):
    E, J = [], []
    for line in open(os.path.join(PDIR, fn)):
        s = [x.strip() for x in line.split('|')]
        if len(s) >= 4 and s[2] != '' and s[3] != '':
            try:
                e = float(s[3])
                j = _eval_j(s[2])
            except ValueError:
                continue
            E.append(e)
            J.append(j)
        if len(E) == n:
            break
    return np.array(E), 2 * np.array(J) + 1


def _read_A(fn):
    rows = []
    for line in open(os.path.join(PDIR, fn)):
        if line.startswith('***') or line.startswith('Aij') or '1/s' in line:
            continue
        s = line.split()
        if s:
            rows.append([float(x) for x in s])
    return np.array(rows)


def _read_coll(fn):
    """Returns log10(T/K) grid and {(i, j): Upsilon array}; honours the file's T_UNIT header."""
    lt, om, unit = None, {}, 'log(K)'
    for line in open(os.path.join(PDIR, fn)):
        if line.startswith('***'):
            if 'T_UNIT' in line:
                unit = line.split('T_UNIT')[1].strip()
            continue
        s = line.split()
        if len(s) < 3:
            continue
        i, j = int(s[0]), int(s[1])
        v = np.array([float(x) for x in s[2:]])
        if i == 0 and j == 0:
            lt = v
        else:
            om[(i, j)] = v
    if unit.startswith('K/10000'):
        lt = np.log10(lt * 1e4)
    elif unit.strip() == 'K':
        lt = np.log10(lt)
    return lt, om


class Atom:
    """5-level atom; collisional excitation/de-excitation by electrons + radiative decay (statistical equilibrium)."""

    def __init__(self, levfile, afile, collfile, nlev=5):
        A = _read_A(afile)
        self.n = min(nlev, A.shape[0])
        self.A = A[:self.n, :self.n]
        self.E, self.g = _read_levels(levfile, self.n)
        self.lt, self.om = _read_coll(collfile)

    def pops(self, T, ne):
        n = self.n
        lT = np.clip(np.log10(T), self.lt[0], self.lt[-1])
        R = np.zeros((n, n))
        kT = T / 1.4387769
        for (i, j), v in self.om.items():
            if i > n or j > n:
                continue
            ups = np.interp(lT, self.lt, v)
            a, b = i - 1, j - 1
            dE = self.E[b] - self.E[a]
            qd = 8.629e-6 * ups / (self.g[b] * np.sqrt(T))
            qu = qd * self.g[b] / self.g[a] * np.exp(-dE / kT)
            R[a, b] += ne * qu
            R[b, a] += ne * qd
        for j in range(n):
            for i in range(j):
                R[j, i] += self.A[j, i]
        M = R.T.copy()
        M[np.diag_indices(n)] -= R.sum(axis=1)
        M[0, :] = 1.0
        rhs = np.zeros(n)
        rhs[0] = 1.0
        return np.linalg.solve(M, rhs)

    def emis(self, T, ne, j, i):
        """Energy emissivity per ion per n_e (erg cm^3 s^-1), upper j -> lower i (1-based)."""
        p = self.pops(T, ne)
        return p[j - 1] * self.A[j - 1, i - 1] * H_PL * C_L * (self.E[j - 1] - self.E[i - 1]) / ne


def oii_atoms():
    """O+ with two collision-strength sets x two A-value sets (PyNeb data files on disk)."""
    out = {}
    for cn, cf in (('Kal09', 'o_ii_coll_Kal09.dat'), ('T07', 'o_ii_coll_T07.dat')):
        for an, af in (('FFT04', 'o_ii_atom_FFT04.dat'), ('Z82-WFD96', 'o_ii_atom_Z82-WFD96.dat')):
            out[f'{cn}+{an}'] = Atom('o_ii_levels.dat', af, cf, 5)
    return out


def sii_atom():
    return Atom('s_ii_levels.dat', 's_ii_atom_IFF05.dat', 's_ii_coll_TZ10.dat', 5)


def oii_ratio(at, T, ne):
    """[O II] 3729/3726 energy ratio. Levels (1-based): 1 4S3/2, 2 2D5/2, 3 2D3/2, 4 2P3/2, 5 2P1/2."""
    return at.emis(T, ne, 2, 1) / at.emis(T, ne, 3, 1)


def oii_7325_3727(at, T, ne):
    """([O II] 7319+7320+7330+7331) / (3726+3729), energy."""
    num = sum(at.emis(T, ne, j, i) for j in (4, 5) for i in (2, 3))
    den = at.emis(T, ne, 2, 1) + at.emis(T, ne, 3, 1)
    return num / den


def sii_ratio(at, T, ne):
    """[S II] 6716/6731 energy ratio. Levels: 1 4S3/2, 2 2D3/2, 3 2D5/2 (S+ ordering: 2D3/2 below 2D5/2)."""
    # identify levels by wavelength: 6716 from 2D5/2, 6731 from 2D3/2
    lam = {j: 1e8 / (at.E[j - 1] - at.E[0]) for j in (2, 3)}
    j6716 = min(lam, key=lambda j: abs(lam[j] - 6718.3))
    j6731 = min(lam, key=lambda j: abs(lam[j] - 6732.7))
    return at.emis(T, ne, j6716, 1) / at.emis(T, ne, j6731, 1)


def invert(fun, target, ne_grid=np.logspace(-2, 5, 701)):
    """n_e at which fun(ne) = target (monotonic decreasing in ne); NaN outside the range."""
    vals = np.array([fun(n) for n in ne_grid])
    if target >= vals[0]:
        return 0.0
    if target <= vals[-1]:
        return np.inf
    return float(np.interp(-target, -vals, ne_grid))


def grid_fit(Ys, Wi, builders):
    """Vectorised GLS over many data vectors. Ys: [n, npix] (row 0 = data, rest = bootstrap replicates);
    builders: iterable of (params, X[k, npix]). For every row returns the best chi2, its params and coefficients."""
    q = np.einsum('ip,pq,iq->i', Ys, Wi, Ys)
    YW = Ys @ Wi
    best_c = np.full(Ys.shape[0], np.inf)
    best_par = [None] * Ys.shape[0]
    best_p = [None] * Ys.shape[0]
    for par, X in builders:
        A = X @ Wi @ X.T
        Ai = np.linalg.inv(A)
        B = YW @ X.T
        P = B @ Ai
        c2 = q - np.sum(B * P, 1)
        better = c2 < best_c
        if better.any():
            best_c[better] = c2[better]
            for i in np.where(better)[0]:
                best_par[i] = par
                best_p[i] = P[i]
    return best_c, best_par, best_p


def chi2_curve(Ys, Wi, X_list):
    """chi2 of every row for every design in X_list -> [len(X_list), n]."""
    q = np.einsum('ip,pq,iq->i', Ys, Wi, Ys)
    YW = Ys @ Wi
    out = []
    for X in X_list:
        A = X @ Wi @ X.T
        B = YW @ X.T
        out.append(q - np.sum(B * (B @ np.linalg.inv(A)), 1))
    return np.array(out)
