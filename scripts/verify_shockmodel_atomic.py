"""verify_shockmodel_atomic: atomic data and rate coefficients for the shockmodel 1-D post-shock code.

All data are parsed (never executed) from text files downloaded on 2026-11-16 into
data/raw/verify_shockmodel_atomic/ (URLs + byte counts in the MANIFEST text files / note):

  * Collisional ionization: Voronov (1997, ADNDT 65, 1) fits, from D. Verner's cfit.f DATA statements.
  * Radiative recombination: Badnell (2006) fits (amdpp.phys.strath.ac.uk/tamoc/RR/clist_K, ground level M=1);
    fallback Verner rrfit.f (Verner & Ferland 1996; Pequignot+91; Shull & van Steenberg 82; ...).
  * Dielectronic recombination: Badnell et al. (2003+) fits (tamoc/DR/clist_K, M=1). Sequences up to Al-like only
    (N<=13): no DR for S+ (P-like) and S2+ (Si-like); Fe is not tracked (CIE cooling only).
  * H charge transfer: Kingdon & Ferland (1996) fits + Stancil+99 (O+/O), Barragan+06 (O++, N++),
    transcribed from the Cloudy c17 source listing (atmdat_char_tran.cpp, doxygen HTML) -> parsed numbers.
  * Cooling: Gnat & Ferland (2012, ApJS 199, 20) ion-by-ion cooling efficiencies (erg cm^3 s^-1, T=1e4-1e8 K).
  * Forbidden/intercombination/resonance lines: PyNeb atomic data files (A-values, collision strengths, levels).
  * H I / He II recombination lines: Pequignot, Petitjean & Boisson (1991) fits (PyNeb *_rec_P91.func).

Conventions: element Z, ion charge q = 0..Z.  Rates:
  Cion(Z,q,T)  q -> q+1 (electron impact), cm^3/s
  Arec(Z,q,T)  q -> q-1 (RR + DR), cm^3/s   (q >= 1)
  CTrec(Z,q,T) X^q+ + H0 -> X^(q-1)+ + H+     CTion(Z,q,T) X^q+ + H+ -> X^(q+1)+ + H0
"""
import os, re
import numpy as np

ROOT = '/Users/mickey/astronomy/silk-nebula'
ADIR = os.path.join(ROOT, 'data/raw/verify_shockmodel_atomic')
PDIR = os.path.join(ADIR, 'pyneb')
K_B = 1.380649e-16
EV_K = 11604.518
H_PL = 6.62607e-27
C_L = 2.99792458e10

ELEMENTS = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Mg': 12, 'Si': 14, 'S': 16, 'Fe': 26}
GF_FILE = {1: 'A', 2: 'B', 6: 'F', 7: 'G', 8: 'H', 10: 'J', 12: 'L', 14: 'N', 16: 'P', 26: 'Z'}

# ---------------------------------------------------------------- abundances (by number, rel. to H)
ABUND = {
    # Asplund et al. 2009 (ARA&A 47, 481) photospheric
    'asplund09': {'He': 10**(10.93-12), 'C': 10**(8.43-12), 'N': 10**(7.83-12), 'O': 10**(8.69-12),
                  'Ne': 10**(7.93-12), 'Mg': 10**(7.60-12), 'Si': 10**(7.51-12), 'S': 10**(7.12-12),
                  'Fe': 10**(7.50-12)},
    # Allen et al. 2008 Table 1 'Solar' (log by number; used only to validate against their Table 6)
    'allen08': {'He': 10**-1.01, 'C': 10**-3.44, 'N': 10**-3.95, 'O': 10**-3.07, 'Ne': 10**-3.91,
                'Mg': 10**-4.42, 'Si': 10**-4.45, 'S': 10**-4.79, 'Fe': 10**-4.63},
    # Asplund et al. 2005 (used by Gnat & Sternberg 2007) -- for the GS07 validation
    'asplund05': {'He': 10**(10.93-12), 'C': 10**(8.39-12), 'N': 10**(7.78-12), 'O': 10**(8.66-12),
                  'Ne': 10**(7.84-12), 'Mg': 10**(7.53-12), 'Si': 10**(7.51-12), 'S': 10**(7.14-12),
                  'Fe': 10**(7.45-12)},
}
# 'ism_dep': Asplund09 with refractory depletion typical of warm diffuse gas (Mg, Si, Fe gas fraction 0.2)
ABUND['ism_dep'] = dict(ABUND['asplund09'])
for _e in ('Mg', 'Si', 'Fe'):
    ABUND['ism_dep'][_e] *= 0.2


# ---------------------------------------------------------------- Voronov collisional ionization
def _load_voronov():
    cf = {}
    pat = re.compile(r"DATA\(CF\(I,\s*(\d+),\s*(\d+)\),I=1,5\)/\s*([^/]+)/")
    with open(os.path.join(ADIR, 'cfit.f')) as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                iz, ine = int(m.group(1)), int(m.group(2))
                vals = [float(v) for v in m.group(3).replace('E', 'e').split(',')]
                cf[(iz, ine)] = vals
    return cf


_CF = _load_voronov()


def Cion(Z, q, T):
    """Voronov 1997 direct collisional ionization rate of X^q+ (q=0..Z-1), cm^3 s^-1."""
    if q >= Z:
        return np.zeros_like(np.asarray(T, float))
    dE, P, A, X, Kx = _CF[(Z, Z - q)]
    U = dE / (np.asarray(T, float) / EV_K)
    with np.errstate(over='ignore', under='ignore'):
        c = A * (1.0 + P * np.sqrt(U)) / (X + U) * U**Kx * np.exp(-U)
    return np.where(U > 80.0, 0.0, c)


def ion_potential_eV(Z, q):
    return _CF[(Z, Z - q)][0]


# ---------------------------------------------------------------- Badnell RR / DR
def _load_badnell_rr():
    rr = {}
    with open(os.path.join(ADIR, 'badnell_RR_clist_K')) as fh:
        for line in fh:
            s = line.split()
            if len(s) < 8 or not s[0].isdigit():
                continue
            Z, N, M, W = map(int, s[:4])
            if M != 1:
                continue
            v = [float(x) for x in s[4:]]
            A, B, T0, T1 = v[:4]
            C, T2 = (v[4], v[5]) if len(v) >= 6 else (0.0, 0.0)
            rr[(Z, N)] = (A, B, T0, T1, C, T2)
    return rr


def _load_badnell_dr():
    c, e = {}, {}
    block = None
    with open(os.path.join(ADIR, 'badnell_DR_clist_K')) as fh:
        for line in fh:
            if 'C1' in line and 'Z  N' in line:
                block = 'c'; continue
            if 'E1' in line and 'Z  N' in line:
                block = 'e'; continue
            s = line.split()
            if block is None or len(s) < 5 or not s[0].isdigit():
                continue
            Z, N, M, W = map(int, s[:4])
            if M != 1:
                continue
            vals = [float(x) for x in s[4:]]
            (c if block == 'c' else e)[(Z, N)] = np.array(vals)
    dr = {}
    for k in c:
        if k in e:
            n = min(len(c[k]), len(e[k]))
            dr[k] = (c[k][:n], e[k][:n])
    return dr


def _load_verner_rr():
    """Parse Verner rrfit.f data statements. Keys use in = number of electrons of the *product* ion."""
    rrec, rnew, fe = {}, {}, {}
    txt = open(os.path.join(ADIR, 'rrfit.f')).read()
    # join continuation lines (column-6 '&' or '*')
    txt = re.sub(r"\n\s{5}[&*]", " ", txt)
    for m in re.finditer(r"data\s*\(rrec\(i,\s*(\d+),\s*(\d+)\),i=1,2\)\s*/([^/]+)/", txt, re.I):
        rrec[(int(m.group(1)), int(m.group(2)))] = [float(x.replace('d', 'e')) for x in m.group(3).split(',')]
    for m in re.finditer(r"data\s*\(rnew\(i,\s*(\d+),\s*(\d+)\),i=1,4\)\s*/([^/]+)/", txt, re.I):
        rnew[(int(m.group(1)), int(m.group(2)))] = [float(x.replace('d', 'e')) for x in m.group(3).split(',')]
    return rrec, rnew


_RR = _load_badnell_rr()
_DR = _load_badnell_dr()
_VRREC, _VRNEW = _load_verner_rr()


def RR(Z, q, T):
    """Radiative recombination X^q+ -> X^(q-1)+ (q>=1)."""
    T = np.asarray(T, float)
    N = Z - q  # electrons of recombining ion
    if (Z, N) in _RR:
        A, B, T0, T1, C, T2 = _RR[(Z, N)]
        Bp = B + C * np.exp(-T2 / T) if C != 0 else B
        st0, st1 = np.sqrt(T / T0), np.sqrt(T / T1)
        return A / (st0 * (1 + st0)**(1 - Bp) * (1 + st1)**(1 + Bp))
    ine = N + 1  # Verner: electrons of product ion
    if (Z, ine) in _VRNEW and (ine <= 3 or ine == 11 or 6 <= Z <= 8 or Z == 10):
        a, b, T0, T1 = _VRNEW[(Z, ine)]
        tt = np.sqrt(T / T0)
        return a / (tt * (tt + 1)**(1 - b) * (1 + np.sqrt(T / T1))**(1 + b))
    if (Z, ine) in _VRREC:
        a, b = _VRREC[(Z, ine)]
        return a / (T * 1e-4)**b
    raise KeyError(f'no RR data for Z={Z} q={q}')


def DR(Z, q, T):
    T = np.asarray(T, float)
    N = Z - q
    if (Z, N) not in _DR:
        return np.zeros_like(T)
    c, E = _DR[(Z, N)]
    with np.errstate(under='ignore'):
        return T**-1.5 * np.sum(c[:, None] * np.exp(-E[:, None] / np.atleast_1d(T)[None, :]), axis=0).reshape(T.shape)


def Arec(Z, q, T):
    if q == 0:
        return np.zeros_like(np.asarray(T, float))
    return RR(Z, q, T) + DR(Z, q, T)


def has_DR(Z, q):
    return (Z, Z - q) in _DR


# ---------------------------------------------------------------- charge transfer with H
def _load_kingdon():
    raw = open(os.path.join(ADIR, 'cloudy_ct_src.txt')).read()
    raw = re.sub(r"(?m)^\s*\d+ ", "", raw)          # strip doxygen listing line numbers
    raw = re.sub(r"/\*.*?\*/", "", raw, flags=re.S)
    ion, rec = {}, {}
    pat = re.compile(r"_itmp\d+\[\]\s*=\s*\{([^}]*)\};.*?(CTIonData|CTRecombData)\[(\d+)\]\[(\d+)\]", re.S)
    for m in pat.finditer(raw):
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(1))]
        # guard: the regex may run across a later block; only accept if the array block is short
        if len(nums) not in (7, 8):
            continue
        nel, istage = int(m.group(3)), int(m.group(4))
        (ion if m.group(2) == 'CTIonData' else rec)[(nel + 1, istage)] = nums
    return ion, rec


_KI, _KR = _load_kingdon()


def CTrec(Z, q, T):
    """X^q+ + H0 -> X^(q-1)+ + H+ (q>=1); cm^3 s^-1."""
    T = np.asarray(T, float)
    if q == 0 or Z == 1:
        return np.zeros_like(T)
    if Z == 8 and q == 1:  # Stancil+99 fit as in Cloudy c17
        a, b, c, d, e, f = 2.3344302e-10, 2.3651505e-10, -1.3146803e-10, 2.9979994e-11, -2.8577012e-12, 1.1963502e-13
        l = np.log(np.clip(T, 10., 1e4))   # fit range of Stancil+99 tables; clamp above 1e4 K
        return ((((f*l + e)*l + d)*l + c)*l + b)*l + a
    if Z == 8 and q == 2:  # Barragan+06
        return np.where(T <= 1500., 0.5337e-9*(T/100.)**-0.076, 0.4344e-9 + 0.6340e-9*np.log10(np.clip(T, 1500., 5e4)/1500.)**2.06)
    if Z == 7 and q == 2:
        return np.where(T <= 1500., 0.8692e-9*(T/1500.)**0.17,
                        np.where(T <= 2e4, 0.9703e-9*(T/1e4)**0.058,
                                 1.0101e-9 + 1.4589e-9*np.log10(np.clip(T, 2e4, 5e4)/2e4)**2.06))
    key = (Z, q - 1)
    if q - 1 > 3:  # Dalgarno estimate used by Cloudy for highly charged ions
        return np.full_like(T, 1.92e-9 * q)
    if key not in _KR or _KR[key][0] <= 0:
        return np.zeros_like(T)
    a, b, c, d, tmin, tmax = _KR[key][:6]
    t = np.clip(T, tmin, tmax) * 1e-4
    return a * 1e-9 * t**b * (1 + c * np.exp(d * t))


def CTion(Z, q, T):
    """X^q+ + H+ -> X^(q+1)+ + H0; cm^3 s^-1 (only q=0 has data)."""
    T = np.asarray(T, float)
    if Z == 1 or q > 0:
        return np.zeros_like(T)
    if Z == 8:
        a, b, c, d, e, f, g = -7.6767404e-14, -3.7282001e-13, -1.488594e-12, -3.6606214e-12, 2.0699463e-12, -2.6139493e-13, 1.1580844e-14
        l = np.log(np.clip(T, 200., 1e4))
        return np.maximum((((((g*l + f)*l + e)*l + d)*l + c)*l + b)*l + a, 0.0)
    key = (Z, 0)
    if key not in _KI or _KI[key][0] <= 0:
        return np.zeros_like(T)
    a, b, c, d, tmin, tmax, dE = _KI[key][:7]
    t = np.clip(T, tmin, tmax) * 1e-4
    return a * 1e-9 * t**b * (1 + c * np.exp(d * t)) * np.exp(-dE * 1e4 / T)


# ---------------------------------------------------------------- Gnat & Ferland 2012 cooling
def _load_gf12(Z):
    path = os.path.join(ADIR, f'gnat_ferland2012_datafile1{GF_FILE[Z]}.txt')
    rows = []
    for line in open(path):
        s = line.split()
        if len(s) == Z + 3:
            try:
                rows.append([float(x) for x in s])
            except ValueError:
                pass
    a = np.array(rows)
    return np.log10(a[:, 0]), np.log10(np.maximum(a[:, 1:Z + 2], 1e-40)), np.log10(np.maximum(a[:, Z + 2], 1e-40))


_GF = {Z: _load_gf12(Z) for Z in GF_FILE}


def cool_ion(Z, T):
    """Array (Z+1,) of Gnat&Ferland ion cooling efficiencies at scalar T (clamped to 1e4..1e8)."""
    lt, lc, _ = _GF[Z]
    x = np.clip(np.log10(T), lt[0], lt[-1])
    return 10**np.array([np.interp(x, lt, lc[:, i]) for i in range(Z + 1)])


def cool_cie(Z, T):
    lt, _, lcie = _GF[Z]
    return 10**np.interp(np.clip(np.log10(T), lt[0], lt[-1]), lt, lcie)


# ---------------------------------------------------------------- n-level atoms (PyNeb files)
def _read_levels(fn, n):
    E = []
    J = []
    for line in open(os.path.join(PDIR, fn)):
        s = [x.strip() for x in line.split('|')]
        if len(s) >= 4 and s[2] != '' and s[3] != '':
            try:
                e = float(s[3]); j = eval_j(s[2])
            except ValueError:
                continue
            E.append(e); J.append(j)
        if len(E) == n:
            break
    return np.array(E), 2 * np.array(J) + 1


def eval_j(s):
    if '/' in s:
        a, b = s.split('/')
        return float(a) / float(b)
    return float(s)


def _read_A(fn):
    rows = []
    for line in open(os.path.join(PDIR, fn)):
        if line.startswith('***') or line.startswith('Aij') or '1/s' in line:
            continue
        s = line.split()
        if s:
            rows.append([float(x) for x in s])
    return np.array(rows)  # A[j, i] = A(j -> i), j upper (row), i lower (col)


def _read_coll(fn):
    lt, om = None, {}
    for line in open(os.path.join(PDIR, fn)):
        if line.startswith('***'):
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
    return lt, om


class Atom:
    """Simple n-level atom: populations from collisional excitation / de-excitation + radiative decay."""

    def __init__(self, name, levfile, afile, collfile, nlev):
        self.name = name
        A = _read_A(afile)
        self.n = min(nlev, A.shape[0])
        self.A = A[:self.n, :self.n]
        self.E, self.g = _read_levels(levfile, self.n)
        self.lt, self.om = _read_coll(collfile)
        self.lt_range = (self.lt[0], self.lt[-1])

    def pops(self, T, ne):
        n = self.n
        lT = np.clip(np.log10(T), *self.lt_range)
        R = np.zeros((n, n))
        kT_cm = T / 1.438777  # kT in cm^-1
        for (i, j), v in self.om.items():
            if i > n or j > n:
                continue
            ups = np.interp(lT, self.lt, v)
            if ups <= 0:
                continue
            a, b = i - 1, j - 1
            dE = self.E[b] - self.E[a]
            qd = 8.629e-6 * ups / (self.g[b] * np.sqrt(T))
            qu = qd * self.g[b] / self.g[a] * np.exp(-dE / kT_cm)
            R[a, b] += ne * qu   # a -> b
            R[b, a] += ne * qd
        for j in range(n):
            for i in range(j):
                R[j, i] += self.A[j, i]
        M = R.T.copy()
        M[np.diag_indices(n)] -= R.sum(axis=1)
        M[0, :] = 1.0
        rhs = np.zeros(n); rhs[0] = 1.0
        return np.linalg.solve(M, rhs)

    def emiss(self, T, ne, j, i):
        """Emissivity per ion per electron, erg cm^3 s^-1 (4pi), for transition upper j -> lower i (1-based)."""
        p = self.pops(T, ne)
        nu_cm = self.E[j - 1] - self.E[i - 1]
        return p[j - 1] * self.A[j - 1, i - 1] * H_PL * C_L * nu_cm / ne

    def wl(self, j, i):
        return 1e8 / (self.E[j - 1] - self.E[i - 1])


def load_atoms():
    return {
        'O3': Atom('O III', 'o_iii_levels.dat', 'o_iii_atom_FFT04-SZ00.dat', 'o_iii_coll_AK99.dat', 6),
        'O2': Atom('O II', 'o_ii_levels.dat', 'o_ii_atom_FFT04.dat', 'o_ii_coll_Kal09.dat', 5),
        'N2': Atom('N II', 'n_ii_levels.dat', 'n_ii_atom_FFT04.dat', 'n_ii_coll_T11.dat', 5),
        'S2': Atom('S II', 's_ii_levels.dat', 's_ii_atom_IFF05.dat', 's_ii_coll_TZ10.dat', 5),
        'Ne3': Atom('Ne III', 'ne_iii_levels.dat', 'ne_iii_atom_FFT04.dat', 'ne_iii_coll_McL11.dat', 5),
        'C4': Atom('C IV', 'c_iv_levels.dat', 'c_iv_atom_WFD96.dat', 'c_iv_coll_AK04.dat', 3),
    }


# ---------------------------------------------------------------- recombination lines (PPB91 fits)
def _ppb91(fn):
    out = {}
    for line in open(os.path.join(PDIR, fn)):
        s = line.split()
        if len(s) >= 8 and s[2] in ('A', 'B'):
            try:
                out[(float(s[0]), s[2])] = tuple(float(x) for x in s[3:7])
            except ValueError:
                pass
    return out


_PH = _ppb91('h_i_rec_P91.func')
_PHE = _ppb91('he_ii_rec_P91.func')


def alpha_eff(line, T, case='B'):
    """Effective recombination coefficient (cm^3/s): 1e-13 z a t^b/(1+c t^d), t=T/(1e4 z^2) (PPB91)."""
    if line in (4861.0, 6563.0):
        a, b, c, d = _PH[(line, case)]; z = 1.0
    else:
        a, b, c, d = _PHE[(line, 'B')]; z = 2.0
    t = np.asarray(T, float) / (1e4 * z * z)
    return 1e-13 * z * a * t**b / (1 + c * t**d)


def rec_emiss(line, T, case='B'):
    """erg cm^3 s^-1 per recombining ion per electron."""
    return alpha_eff(line, T, case) * H_PL * C_L / (line * 1e-8)


if __name__ == '__main__':
    # quick sanity print
    T = np.array([1e4, 3e4, 1e5, 3e5])
    print('Voronov H', Cion(1, 0, T))
    print('O2+ rec', Arec(8, 2, T), 'O3+ rec', Arec(8, 3, T))
    print('CT O+ + H', CTrec(8, 1, T), ' O + H+', CTion(8, 0, T), ' O++ + H', CTrec(8, 2, T))
    print('Kingdon keys rec', sorted(_KR.keys())[:40])
    print('Hb alpha_B(1e4)', alpha_eff(4861.0, 1e4), ' HeII4686', alpha_eff(4686.0, 1e4))
    at = load_atoms()
    for k, a in at.items():
        print(k, a.n, a.E, a.g)
    o3 = at['O3']
    for TT in (1e4, 3e4, 1e5):
        r = o3.emiss(TT, 1.0, 5, 4) / o3.emiss(TT, 1.0, 4, 3)
        print('4363/5007 at', TT, r, 'wl', o3.wl(4, 3), o3.wl(5, 4), o3.wl(6, 3))
    o2 = at['O2']
    for ne in (1, 100, 1000):
        print('[OII]3729/3726 ne', ne, o2.emiss(1e4, ne, 2, 1) / o2.emiss(1e4, ne, 3, 1), o2.wl(2, 1), o2.wl(3, 1))
    print('missing DR (q>=1) for tracked:', [(Z, q) for Z in (1, 2, 6, 7, 8, 10, 12, 14, 16) for q in range(1, Z + 1) if not has_DR(Z, q) and Z > 1][:30])


# ---------------------------------------------------------------- H I collisional excitation of Balmer lines
CDIR = os.path.join(ADIR, 'chianti')


def _load_h1():
    """CHIANTI h_1 (Anderson et al. 2000/2002 Upsilons; hydrogenic A-values), levels 1-25 (n <= 5)."""
    lev = {}
    for line in open(os.path.join(CDIR, 'h_1.elvlc')):
        s = line.split()
        if len(s) >= 7 and s[0].isdigit():
            i = int(s[0]); conf = s[1]; J = float(s[4]); E = float(s[5])
            lev[i] = (conf, 2 * J + 1, E)
        if line.strip().startswith('-1'):
            break
    A = np.zeros((26, 26))
    for line in open(os.path.join(CDIR, 'h_1.wgfa')):
        s = line.split()
        if line.strip().startswith('-1'):
            break
        if len(s) >= 5:
            try:
                i, j = int(s[0]), int(s[1]); wl = float(s[2]); a = float(s[4])
            except ValueError:
                continue
            if wl > 0 and j <= 25:
                A[j, i] += a        # A[upper, lower]
    ups = {}
    L = open(os.path.join(CDIR, 'h_1.scups')).read().split('\n')
    k = 0
    while k < len(L):
        s = L[k].split()
        if not s or s[0] == '-1':
            break
        i, j = int(s[0]), int(s[1]); dE = float(s[2]); hi = float(s[4]) if s[4] != '-1' else None
        nsp, typ, C = int(s[5]), int(s[6]), float(s[7])
        xs = np.array([float(v) for v in L[k + 1].split()])
        ys = np.array([float(v) for v in L[k + 2].split()])
        ups[(i, j)] = (dE, nsp, typ, C, xs, ys)
        k += 3
    return lev, A, ups


def _bt_upsilon(dE, typ, C, xs, ys, T):
    kt = T / 157887.5 / dE   # kT/dE with kT in Ryd
    if typ in (1, 4):
        x = 1 - np.log(C) / np.log(kt + C)
    else:
        x = kt / (kt + C)
    y = np.interp(x, xs, ys)
    if typ == 1:
        return y * np.log(kt + np.e)
    if typ == 2:
        return y
    if typ == 3:
        return y / (kt + 1)
    if typ == 4:
        return y * np.log(kt + C)
    return y


_H1 = None


def h_balmer_coll(T):
    """Collisional excitation (from H0 1s, low density, case A cascades through n<=5) photon-production
    coefficients for Halpha and Hbeta, cm^3 s^-1 per H0 per electron; and for Ly-alpha (for reference)."""
    global _H1
    if _H1 is None:
        _H1 = _load_h1()
    lev, A, ups = _H1
    n = 25
    q = np.zeros(n + 1)
    for j in range(2, n + 1):
        if (1, j) not in ups:
            continue
        dE, nsp, typ, C, xs, ys = ups[(1, j)]
        U = _bt_upsilon(dE, typ, C, xs, ys, T)
        q[j] = 8.629e-6 * U / (2.0 * np.sqrt(T)) * np.exp(-dE * 157887.5 / T)
    # level n-quantum number
    nq = {i: int(re.match(r'(\d+)', lev[i][0]).group(1)) for i in lev if i <= n}
    Atot = A[:n + 1, :n + 1].sum(axis=1)
    # cascade: population flux P_j = q_j + sum_{k>j} P_k A_kj/Atot_k  (photons through each level)
    P = q.copy()
    for k in range(n, 1, -1):
        if Atot[k] <= 0:
            continue
        for j in range(1, k):
            P[j] += P[k] * A[k, j] / Atot[k] if j > 1 else 0.0
    ha = sum(P[k] * A[k, j] / Atot[k] for k in range(2, n + 1) for j in range(1, k)
             if Atot[k] > 0 and nq[k] == 3 and nq[j] == 2)
    hb = sum(P[k] * A[k, j] / Atot[k] for k in range(2, n + 1) for j in range(1, k)
             if Atot[k] > 0 and nq[k] == 4 and nq[j] == 2)
    lya = sum(P[k] * A[k, 1] / Atot[k] for k in range(2, n + 1) if Atot[k] > 0 and nq[k] == 2)
    return ha, hb, lya
