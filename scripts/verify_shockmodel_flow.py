"""verify_shockmodel_flow: simplified 1-D steady-flow non-equilibrium radiative-shock model (shockmodel agent).

Physics (approximations, see notes/verify_shockmodel.md):
  * plane-parallel steady flow behind a perpendicular MHD shock (transverse B frozen in, B ∝ rho);
    exact mass / momentum / energy flux conservation; Rankine-Hugoniot jump from the same relations.
  * T_e = T_ion immediately behind the shock (no Coulomb-equilibration delay).
  * non-equilibrium ionization of H, He, C, N, O, Ne, Mg, Si, S (all stages): Voronov collisional ionization,
    Badnell RR + DR (Verner fallback), H charge transfer (Kingdon & Ferland + updates).  Fe in CIE (cooling only).
  * radiative cooling = n_e * sum_ions n_ion * Lambda_ion(T) (Gnat & Ferland 2012 ion-by-ion efficiencies).
  * NO radiative transfer: no photoionization of the post-shock gas by its own EUV / no precursor computed
    self-consistently (the pre-shock ionization state is an input).  Optically thin lines.
  * optional 'complete' tail: when T reaches T_floor (default 1e4 K) T is held fixed (mimicking photo-heating)
    and the gas recombines isothermally until x(H+) < 1e-2.

Integrator: operator-split; explicit energy update + backward-Euler (implicit) ionization per element,
adaptive step limited to dT/T < 0.4 % and |df| < 0.03.
"""
import sys, numpy as np
from scipy.optimize import brentq
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_atomic as A

M_H = 1.6735575e-24
K_B = 1.380649e-16
AMU = {'H': 1.008, 'He': 4.0026, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'Ne': 20.18, 'Mg': 24.305,
       'Si': 28.085, 'S': 32.06, 'Fe': 55.845}
TRACK = ('H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S')
LT_GRID = np.arange(3.0, 7.5001, 0.002)


class Network:
    def __init__(self, abund='asplund09', track=TRACK, cie=('Fe',), h_case='B'):
        ab = A.ABUND[abund] if isinstance(abund, str) else abund
        self.ab = {'H': 1.0, **ab}
        self.track = list(track)
        self.cie = list(cie)
        self.Z = {e: A.ELEMENTS[e] for e in self.track + self.cie}
        self.slices, i = {}, 0
        for e in self.track:
            self.slices[e] = slice(i, i + self.Z[e] + 1); i += self.Z[e] + 1
        self.nstate = i
        T = 10**LT_GRID
        self.tab = {}
        for e in self.track:
            Z = self.Z[e]
            C = np.array([A.Cion(Z, q, T) for q in range(Z + 1)])
            R = np.array([A.Arec(Z, q, T) for q in range(Z + 1)])
            if e == 'H' and h_case == 'B':
                R[1] = _alphaB_H(T)
            CTr = np.array([A.CTrec(Z, q, T) for q in range(Z + 1)])
            CTi = np.array([A.CTion(Z, q, T) for q in range(Z + 1)])
            L = np.array([A.cool_ion(Z, t) for t in T]).T
            self.tab[e] = (C, R, CTr, CTi, L)
        self.tab_cie = {e: np.array([A.cool_cie(self.Z[e], t) for t in T]) for e in self.cie}
        self.mass_per_H = M_H * sum(self.ab[e] * AMU[e] / 1.008 for e in self.ab)  # g per H nucleus
        self.n_nuc_per_H = sum(self.ab[e] for e in self.ab)

    # ---- interpolation helpers
    def _w(self, T):
        x = (np.log10(T) - LT_GRID[0]) / (LT_GRID[1] - LT_GRID[0])
        i = int(np.clip(np.floor(x), 0, len(LT_GRID) - 2))
        w = np.clip(x - i, 0.0, 1.0)
        return i, w

    def coeffs(self, e, T):
        i, w = self._w(T)
        return [arr[:, i] * (1 - w) + arr[:, i + 1] * w for arr in self.tab[e]]

    def ne_per_H(self, f):
        s = 0.0
        for e in self.track:
            fq = f[self.slices[e]]
            s += self.ab[e] * np.dot(np.arange(len(fq)), fq)
        for e in self.cie:   # assume Fe ~ singly/doubly ionized: count 2 electrons (negligible)
            s += self.ab[e] * 2.0
        return s

    def ntot_per_H(self, f):
        return self.n_nuc_per_H + self.ne_per_H(f)

    def cooling(self, f, T, nH):
        """erg cm^-3 s^-1."""
        ne = nH * self.ne_per_H(f)
        L = 0.0
        i, w = self._w(T)
        for e in self.track:
            Lq = self.tab[e][4][:, i] * (1 - w) + self.tab[e][4][:, i + 1] * w
            L += self.ab[e] * np.dot(Lq, f[self.slices[e]])
        for e in self.cie:
            L += self.ab[e] * (self.tab_cie[e][i] * (1 - w) + self.tab_cie[e][i + 1] * w)
        return ne * nH * L

    def implicit_ion_step(self, f, T, nH, dt):
        """Backward-Euler update of all tracked ionization fractions at fixed T, n_e, n(H0), n(H+)."""
        fn = f.copy()
        ne = nH * self.ne_per_H(f)
        fH = f[self.slices['H']]
        nH0, nHp = nH * fH[0], nH * fH[1]
        for e in self.track:
            C, R, CTr, CTi, _ = self.coeffs(e, T)
            ion = ne * C + nHp * CTi
            rec = ne * R + nH0 * CTr
            if e == 'H':
                ion = ne * C
                rec = ne * R
            n = len(ion)
            M = np.zeros((n, n))
            idx = np.arange(n)
            M[idx, idx] = -(ion + rec)
            M[idx[1:], idx[:-1]] += ion[:-1]       # q-1 -> q
            M[idx[:-1], idx[1:]] += rec[1:]        # q+1 -> q
            Amat = np.eye(n) - dt * M
            sl = self.slices[e]
            x = np.linalg.solve(Amat, f[sl])
            x = np.clip(x, 0.0, None)
            fn[sl] = x / x.sum()
        return fn

    def cie_state(self, T, nH=1.0, ct=True):
        """Collisional ionization equilibrium fractions at T (iterating n_e and H0 for charge transfer)."""
        f = np.zeros(self.nstate)
        for e in self.track:
            sl = self.slices[e]; tmp = np.zeros(sl.stop - sl.start); tmp[-1] = 1; f[sl] = tmp
        for it in range(30):
            ne = nH * self.ne_per_H(f)
            fH = f[self.slices['H']]
            nH0, nHp = nH * fH[0], nH * fH[1]
            fnew = f.copy()
            for e in self.track:
                C, R, CTr, CTi, _ = self.coeffs(e, T)
                ion = ne * C + (nHp * CTi if ct and e != 'H' else 0)
                rec = ne * R + (nH0 * CTr if ct and e != 'H' else 0)
                n = len(ion)
                # equilibrium: f_{q+1}/f_q = ion_q / rec_{q+1}
                lr = np.zeros(n)
                for q in range(n - 1):
                    lr[q + 1] = lr[q] + np.log(max(ion[q], 1e-300)) - np.log(max(rec[q + 1], 1e-300))
                lr -= lr.max()
                x = np.exp(lr); fnew[self.slices[e]] = x / x.sum()
            if np.max(np.abs(fnew - f)) < 1e-9:
                f = fnew; break
            f = fnew
        return f

    def state_from(self, spec, T0=1e4):
        """Build an initial fraction vector from dict {element: {q: frac}}; unspecified elements -> CIE at T0."""
        f = self.cie_state(T0)
        for e, d in spec.items():
            if e not in self.slices:
                continue
            sl = self.slices[e]
            x = np.zeros(sl.stop - sl.start)
            for q, v in d.items():
                x[q] = v
            f[sl] = x / x.sum()
        return f


def _alphaB_H(T):
    # case-B total recombination coefficient of H (Storey & Hummer 1995 total, via the PyNeb SH95 caseB table
    # is density-dependent; use the standard fit alpha_B = 2.59e-13 (T/1e4)^-0.7 (Osterbrock & Ferland 2006))
    return 2.59e-13 * (np.asarray(T) / 1e4)**-0.7


# ------------------------------------------------------------------ preshock ionization states
PRESHOCK = {
    # neutral atomic gas (HI wall-like): C, Mg, Si, S, Fe singly ionized by the ISRF
    'neutral': ({'H': {0: 0.999, 1: 0.001}, 'He': {0: 1.0}, 'C': {1: 1.0}, 'N': {0: 1.0}, 'O': {0: 0.999, 1: 0.001},
                 'Ne': {0: 1.0}, 'Mg': {1: 1.0}, 'Si': {1: 1.0}, 'S': {1: 1.0}}, 100.0),
    # warm ionized medium-like photoionized gas (O mostly O+)
    'wim': ({'H': {0: 0.05, 1: 0.95}, 'He': {0: 0.6, 1: 0.4}, 'C': {1: 0.7, 2: 0.3}, 'N': {0: 0.05, 1: 0.85, 2: 0.1},
             'O': {0: 0.05, 1: 0.85, 2: 0.1}, 'Ne': {0: 0.2, 1: 0.8}, 'Mg': {1: 0.3, 2: 0.7}, 'Si': {1: 0.3, 2: 0.7},
             'S': {1: 0.7, 2: 0.3}}, 8000.0),
    # strongly preionized (precursor of a >~120 km/s shock, Shull & McKee 1979-like): H+, He+, O+/O++
    'precursor': ({'H': {1: 1.0}, 'He': {1: 0.9, 2: 0.1}, 'C': {2: 0.6, 3: 0.4}, 'N': {1: 0.3, 2: 0.6, 3: 0.1},
                   'O': {1: 0.5, 2: 0.5}, 'Ne': {1: 0.5, 2: 0.5}, 'Mg': {2: 1.0}, 'Si': {2: 0.5, 3: 0.5},
                   'S': {2: 0.6, 3: 0.4}}, 1.2e4),
    # fully ionized but low ionization: H+ with metals in CIE at 1e4 K (O mostly O+?)
}


class Shock:
    def __init__(self, net, vs_kms, n0, B0_uG=1.0, pre='wim', T0=None, Tfloor=1e4, complete=False,
                 Nmax=1e21, tmax=None, xHp_stop=1e-2):
        self.net, self.vs, self.n0, self.B0 = net, vs_kms * 1e5, n0, B0_uG * 1e-6
        spec, Tpre = PRESHOCK[pre] if isinstance(pre, str) else pre
        self.T0 = Tpre if T0 is None else T0
        self.f0 = net.state_from(spec, T0=max(self.T0, 1e4))
        self.Tfloor, self.complete, self.Nmax, self.xHp_stop = Tfloor, complete, Nmax, xHp_stop
        self.tmax = tmax
        rho0 = n0 * net.mass_per_H
        self.rho0 = rho0
        P0 = n0 * net.ntot_per_H(self.f0) * K_B * self.T0
        self.b = self.B0**2 / (8 * np.pi)
        self.J = rho0 * self.vs
        self.Pi = rho0 * self.vs**2 + P0 + self.b
        self.E0 = 0.5 * self.vs**2 + 2.5 * P0 / rho0 + 2 * self.b / rho0
        xs = np.linspace(1.0001, 4.5, 4000)
        Ex = self.Efun(xs)
        self.x_sonic = xs[np.argmax(Ex)]
        # max compression where P -> 0 (magnetic support)
        self.x_max = self._xmax()
        self.x1 = brentq(lambda x: self.Efun(x) - self.E0, self.x_sonic, min(4.0 + 1e-6, self.x_max - 1e-9) + 0.0)

    def _xmax(self):
        if self.b <= 0:
            return 1e9
        g = lambda x: self.Pi - self.J * self.vs / x - self.b * x * x
        hi = 2.0
        while g(hi) > 0:
            hi *= 2
        return brentq(g, 1.0001, hi)

    def P_of_x(self, x):
        return self.Pi - self.J * self.vs / x - self.b * x * x

    def Efun(self, x):
        v = self.vs / x
        return 0.5 * v * v + 2.5 * self.P_of_x(x) / (self.rho0 * x) + 2 * self.b * x / self.rho0

    def x_of_E(self, E, xprev):
        lo = max(self.x_sonic, 1.0001)
        hi = self.x_max * (1 - 1e-12) if self.b > 0 else 1e8
        if self.Efun(hi) > E:   # cannot cool further at this E (magnetically supported limit)
            return hi
        return brentq(lambda x: self.Efun(x) - E, max(lo, xprev * 0.999), hi, xtol=1e-10 * xprev)

    def run(self, dT_rel=0.004, df_max=0.03, verbose=False):
        net = self.net
        f = self.f0.copy()
        E = self.E0
        x = self.x1
        nH = self.n0 * x
        P = self.P_of_x(x)
        T = P / (nH * net.ntot_per_H(f) * K_B)
        self.Ts = T
        t, s = 0.0, 0.0
        rec = []
        L = net.cooling(f, T, nH)
        dt = 1e-4 * (1.5 * P) / L
        isothermal = False
        nstep = 0
        while True:
            nstep += 1
            if not isothermal:
                L = net.cooling(f, T, nH)
                for _ in range(40):
                    Enew = E - L / (self.rho0 * x) * dt
                    xn = self.x_of_E(Enew, x)
                    nHn = self.n0 * xn
                    fn = net.implicit_ion_step(f, T, 0.5 * (nH + nHn), dt)
                    Tn = self.P_of_x(xn) / (nHn * net.ntot_per_H(fn) * K_B)
                    if abs(Tn - T) / T > dT_rel or np.max(np.abs(fn - f)) > df_max:
                        dt *= 0.5; continue
                    break
            else:
                Tn = self.Tfloor
                for _ in range(40):
                    fn = net.implicit_ion_step(f, Tn, nH, dt)
                    if np.max(np.abs(fn - f)) > df_max:
                        dt *= 0.5; continue
                    break
                nt = net.ntot_per_H(fn)
                xn = brentq(lambda xx: self.P_of_x(xx) - self.n0 * xx * nt * K_B * Tn, max(x * 0.5, 1.0001),
                            self.x_max * (1 - 1e-12) if self.b > 0 else 1e9)
                nHn = self.n0 * xn
                Enew = self.Efun(xn)
            v = self.vs / x
            vn = self.vs / xn
            ds = 0.5 * (v + vn) * dt
            rec.append((t, s, T, nH, x, f, dt, ds))
            t += dt; s += ds
            f, T, x, nH, E = fn, Tn, xn, nHn, Enew
            dfmax_used = None
            dt *= 1.4
            N = self.n0 * self.vs * t
            if self.tmax is not None and t >= self.tmax:
                break
            if N >= self.Nmax:
                break
            if not isothermal and T <= self.Tfloor:
                if not self.complete:
                    break
                isothermal = True
            if isothermal and f[net.slices['H']][1] < self.xHp_stop:
                break
            if nstep > 200000:
                raise RuntimeError('too many steps')
        rec.append((t, s, T, nH, x, f, 0.0, 0.0))
        self.nstep = nstep
        self.t = np.array([r[0] for r in rec])
        self.s = np.array([r[1] for r in rec])
        self.T = np.array([r[2] for r in rec])
        self.nH = np.array([r[3] for r in rec])
        self.x = np.array([r[4] for r in rec])
        self.f = np.array([r[5] for r in rec])
        self.ds = np.array([r[7] for r in rec])
        self.dt = np.array([r[6] for r in rec])
        self.N = self.n0 * self.vs * self.t   # swept column of H nuclei behind the shock (cm^-2)
        self.ne = self.nH * np.array([net.ne_per_H(ff) for ff in self.f])
        return self

    def frac(self, el, q):
        return self.f[:, self.net.slices[el].start + q]


# ------------------------------------------------------------------ emission post-processing
LINES = {
    # name: (kind, ion-element, charge, atom-key, upper, lower)  wavelengths from the level energies
    'OIII5007': ('cf', 'O', 2, 'O3', 4, 3), 'OIII4959': ('cf', 'O', 2, 'O3', 4, 2),
    'OIII4363': ('cf', 'O', 2, 'O3', 5, 4), 'OIII1666': ('cf', 'O', 2, 'O3', 6, 3), 'OIII1661': ('cf', 'O', 2, 'O3', 6, 2),
    'OII3729': ('cf', 'O', 1, 'O2', 2, 1), 'OII3726': ('cf', 'O', 1, 'O2', 3, 1),
    'OII7320': ('cf', 'O', 1, 'O2', 4, 2), 'OII7330': ('cf', 'O', 1, 'O2', 4, 3),
    'OII7319b': ('cf', 'O', 1, 'O2', 5, 2), 'OII7330b': ('cf', 'O', 1, 'O2', 5, 3),
    'NII6583': ('cf', 'N', 1, 'N2', 4, 3), 'NII6548': ('cf', 'N', 1, 'N2', 4, 2), 'NII5755': ('cf', 'N', 1, 'N2', 5, 4),
    'SII6716': ('cf', 'S', 1, 'S2', 3, 1), 'SII6731': ('cf', 'S', 1, 'S2', 2, 1),
    'NeIII3869': ('cf', 'Ne', 2, 'Ne3', 4, 1),
    'CIV1548': ('cf', 'C', 3, 'C4', 3, 1), 'CIV1551': ('cf', 'C', 3, 'C4', 2, 1),
    'Hb': ('rec', 'H', 1, 4861.0, None, None), 'Ha': ('rec', 'H', 1, 6563.0, None, None),
    'HeII4686': ('rec', 'He', 2, 4686.0, None, None),
}
LAM = {'OIII5007': 5006.84, 'OIII4959': 4958.91, 'OIII4363': 4363.21, 'OIII1666': 1666.15, 'OIII1661': 1660.81,
       'OII3729': 3728.82, 'OII3726': 3726.03, 'NII6583': 6583.45, 'NII6548': 6548.05, 'NII5755': 5754.6,
       'SII6716': 6716.44, 'SII6731': 6730.82, 'NeIII3869': 3868.76, 'CIV1548': 1548.19, 'CIV1551': 1550.77,
       'Hb': 4861.33, 'Ha': 6562.80, 'HeII4686': 4685.68, 'OII7320': 7320.0, 'OII7330': 7330.0,
       'OII7319b': 7319.0, 'OII7330b': 7330.0}

_ATOMS = None


def atoms():
    global _ATOMS
    if _ATOMS is None:
        _ATOMS = A.load_atoms()
    return _ATOMS


def emissivities(sh, case='B', thin=1):
    """Per-step volume emissivities (erg cm^-3 s^-1, 4pi) for all LINES; returns dict of arrays."""
    at = atoms()
    net = sh.net
    out = {}
    n = len(sh.T)
    for name, (kind, el, q, key, up, lo) in LINES.items():
        j = np.zeros(n)
        nion = sh.nH * net.ab[el] * sh.frac(el, q)
        if kind == 'cf':
            atom = at[key]
            for k in range(0, n, thin):
                if nion[k] <= 0 or sh.ne[k] <= 0:
                    continue
                j[k] = sh.ne[k] * nion[k] * atom.emiss(sh.T[k], sh.ne[k], up, lo)
            if thin > 1:
                good = np.arange(0, n, thin)
                j = np.interp(np.arange(n), good, j[good])
        else:
            cs = case if el == 'H' else 'B'
            j = sh.ne * nion * A.rec_emiss(key, np.maximum(sh.T, 1e3), cs)
        out[name] = j
    # He II 1640: recombination, energy ratio I(1640)/I(4686) ~ 6.8 (case B, ~1e4 K; adopted, see note)
    out['HeII1640'] = 6.8 * out['HeII4686']
    # collisional excitation of H I Balmer lines from H0 (CHIANTI/Anderson+2000 Upsilons, case-A cascades)
    qa, qb = _hcoll_tables()
    lT = np.log10(np.maximum(sh.T, 1e4))
    nH0 = sh.nH * sh.frac('H', 0)
    ha_c = sh.ne * nH0 * 10**np.interp(lT, _HC_LT, qa) * A.H_PL * A.C_L / (LAM['Ha'] * 1e-8)
    hb_c = sh.ne * nH0 * 10**np.interp(lT, _HC_LT, qb) * A.H_PL * A.C_L / (LAM['Hb'] * 1e-8)
    out['Ha_rec'], out['Hb_rec'] = out['Ha'].copy(), out['Hb'].copy()
    out['Ha_coll'], out['Hb_coll'] = ha_c, hb_c
    out['Ha'] = out['Ha'] + ha_c
    out['Hb'] = out['Hb'] + hb_c
    return out


_HC_LT = np.arange(4.0, 7.01, 0.02)
_HC = None


def _hcoll_tables():
    global _HC
    if _HC is None:
        v = np.array([A.h_balmer_coll(10**t)[:2] for t in _HC_LT])
        _HC = (np.log10(np.maximum(v[:, 0], 1e-40)), np.log10(np.maximum(v[:, 1], 1e-40)))
    return _HC


def cumulative(sh, j):
    """Cumulative face-on flux (erg cm^-2 s^-1) vs swept column sh.N, trapezoid in ds."""
    F = {}
    for k, arr in j.items():
        seg = 0.5 * (arr[:-1] + arr[1:]) * np.diff(sh.s)
        F[k] = np.concatenate([[0.0], np.cumsum(seg)])
    return F


def photons(Ferg, name):
    return Ferg / (A.H_PL * A.C_L / (LAM[name] * 1e-8))


if __name__ == '__main__':
    import time
    net = Network('asplund09')
    t0 = time.time()
    sh = Shock(net, 120, 0.1, B0_uG=1.0, pre='wim', complete=True).run()
    print('steps', sh.nstep, 'time', time.time() - t0, 'Ts', sh.Ts, 'x1', sh.x1, 'final T', sh.T[-1], 'N', sh.N[-1])
    j = emissivities(sh, thin=3)
    F = cumulative(sh, j)
    for k in ('OIII5007', 'OII3729', 'OII3726', 'Hb', 'Ha', 'NII6583', 'HeII4686', 'CIV1548', 'OIII4363'):
        print(k, F[k][-1] / F['Hb'][-1])
