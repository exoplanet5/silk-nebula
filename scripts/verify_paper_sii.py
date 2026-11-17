"""verify_paper_sii.py -- independent 5-level S+ density diagnostic (tag: paper).

Atomic data (read as plain-text tables, transcribed from the original papers by the PyNeb team,
https://github.com/Morisset/PyNeb_devel/tree/master/pyneb/atomic_data, fetched 2026-11-16, see
data/derived/verify_paper_raw/pyneb_atomic/fetch_log.txt):
  energies   NIST ASD (s_ii_levels.dat): 4S3/2 0, 2D3/2 14852.94, 2D5/2 14884.73, 2P1/2 24524.83, 2P3/2 24571.54 cm-1
  Omega      Tayal & Zatsarinny 2010, ApJS 188, 32 (s_ii_coll_TZ10.dat); Ramsbottom, Bell & Stafford 1996, ADNDT 63, 57 (RBS96)
  A-values   Tayal & Zatsarinny 2010 (s_ii_atom_TZ10-IFF05.dat), Mendoza & Zeippen 1982, MNRAS 198, 127 (VVF96-MZ82a),
             Rynkun, Gaigalas & Jonsson 2019, A&A 623, A155 (RGJ19)
The level-population solver is my own code. Compared with the paper's eq. 5 (Proxauf, Oettl & Kimeswenger 2014).
Run: ~/.venvs/astro313/bin/python -I scripts/verify_paper_sii.py
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/derived/verify_paper_raw/pyneb_atomic'
sys.path.insert(0, str(ROOT / 'scripts'))

E = np.array([0.0, 14852.94, 14884.73, 24524.83, 24571.54])  # cm-1
G = np.array([4, 4, 6, 2, 4])
HCK = 1.438777  # cm K


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
        if not p or line.startswith('***') or p[0] in ('Aij',) or p[0].startswith('1/s'):
            continue
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            continue
    A = np.array([r[:nlev] for r in rows[:nlev]])
    return A  # A[upper, lower]


COLL = {'TZ10': read_coll(RAW / 's_ii_coll_TZ10.dat.txt'), 'RBS96': read_coll(RAW / 's_ii_coll_RBS96.dat.txt')}
AVAL = {'TZ10': read_A(RAW / 's_ii_atom_TZ10-IFF05.dat.txt'), 'MZ82': read_A(RAW / 's_ii_atom_VVF96-MZ82a.dat.txt'),
        'RGJ19': read_A(RAW / 's_ii_atom_RGJ19.dat.txt')}


def populations(ne, T, coll='TZ10', aval='TZ10'):
    logT, om = COLL[coll]
    A = AVAL[aval]
    n = 5
    M = np.zeros((n, n))
    for (i, j), vals in om.items():
        Om = np.interp(np.log10(T), logT, vals)
        dE = (E[j] - E[i]) * HCK / T
        qd = 8.629e-6 / (G[j] * np.sqrt(T)) * Om        # j -> i (de-excitation)
        qu = qd * G[j] / G[i] * np.exp(-dE)             # i -> j
        M[j, i] += ne * qu; M[i, i] -= ne * qu
        M[i, j] += ne * qd; M[j, j] -= ne * qd
    for u in range(n):
        for l in range(u):
            M[l, u] += A[u, l]; M[u, u] -= A[u, l]
    M[0, :] = 1.0
    rhs = np.zeros(n); rhs[0] = 1.0
    return np.linalg.solve(M, rhs)


def ratio(ne, T, **kw):
    """I(6716)/I(6731) in energy units: level 3 (2D5/2) -> 1 over level 2 (2D3/2) -> 1."""
    p = populations(ne, T, **kw)
    A = AVAL[kw.get('aval', 'TZ10')]
    return (p[2] * A[2, 0] * E[2]) / (p[1] * A[1, 0] * E[1])


def ne_from_R(R, T, **kw):
    f = lambda lg: ratio(10 ** lg, T, **kw) - R
    lo, hi = -2.0, 6.0
    if f(lo) < 0:   # above low-density limit
        return np.nan if R > ratio(1e-2, T, **kw) + 1e-9 else 1e-2
    if f(hi) > 0:
        return np.inf
    return 10 ** brentq(f, lo, hi, xtol=1e-6)


def proxauf(R):
    return 10 ** (0.0543 * np.tan(-3.0553 * R + 2.8506) + 6.98 - 10.6905 * R + 9.9186 * R ** 2 - 3.5442 * R ** 3)


def main():
    Rs = [1.03, 1.10, 1.14, 1.19]
    region = {1.10: 'P1', 1.14: 'P2', 1.19: 'P3', 1.03: 'P4'}
    combos = [('TZ10', 'TZ10'), ('TZ10', 'MZ82'), ('TZ10', 'RGJ19'), ('RBS96', 'MZ82')]
    res = {'atomic_data': {'energies_cm-1': E.tolist(), 'g': G.tolist(),
                           'A_6716_6731': {k: [float(v[2, 0]), float(v[1, 0])] for k, v in AVAL.items()},
                           'Omega_1e4K_12_13_23': {k: [float(np.interp(4.0, c[0], c[1][(0, 1)])), float(np.interp(4.0, c[0], c[1][(0, 2)])),
                                                       float(np.interp(4.0, c[0], c[1][(1, 2)]))] for k, c in COLL.items()}},
           'limits': {}, 'ne': {}, 'R_for_ne': {}}
    for T in (8000.0, 1.0e4):
        for c, a in combos:
            key = f'{c}+{a}@{int(T)}K'
            res['limits'][key] = dict(R_low=round(ratio(1e-3, T, coll=c, aval=a), 4), R_at_1=round(ratio(1.0, T, coll=c, aval=a), 4),
                                      R_high=round(ratio(1e7, T, coll=c, aval=a), 4))
            res['ne'][key] = {f'{R:.2f}({region[R]})': round(ne_from_R(R, T, coll=c, aval=a), 1) for R in Rs}
            res['R_for_ne'][key] = {f'{n}': round(ratio(n, T, coll=c, aval=a), 4) for n in (1, 3, 10, 64, 100, 150, 300)}
    res['proxauf_eq5'] = {f'{R:.2f}({region[R]})': round(float(proxauf(R)), 1) for R in Rs}
    res['proxauf_eq5_T8000_sqrtT'] = {k: round(v * np.sqrt(0.8), 1) for k, v in res['proxauf_eq5'].items()}
    rr = np.linspace(1.0, 1.45, 4501)
    pp = proxauf(rr)
    res['proxauf_R_for_ne'] = {str(n): round(float(rr[np.argmin(np.abs(pp - n))]), 3) for n in (64, 100, 150)}
    # proxauf vs my 1e4 K TZ10+TZ10 curve
    grid = np.logspace(0, 4, 81)
    mine = np.array([ratio(n, 1e4) for n in grid])
    res['proxauf_vs_mine_1e4'] = {f'{n:.0f}': [round(float(r), 4), round(float(rr[np.argmin(np.abs(pp - n))]), 4)] for n, r in zip(grid[::10], mine[::10])}
    # Monte Carlo on Table 2 intensities (I6716 +- 0.02, I6731 +- 0.02)
    rng = np.random.default_rng(7)
    tab = {'P1': (0.14, 0.13, 1.10), 'P2': (0.17, 0.15, 1.14), 'P3': (0.11, 0.10, 1.19), 'P4': (0.10, 0.10, 1.03)}
    Rl = ratio(1e-3, 1e4)
    nev = np.vectorize(lambda R: ne_from_R(R, 1e4) if R < Rl else 1.0)
    mc = {}
    for k, (a, b, Rt) in tab.items():
        R = rng.normal(a, 0.02, 4000) / rng.normal(b, 0.02, 4000)
        R = R[(R > 0.3) & (R < 2.5)]
        q = np.percentile(R, [16, 50, 84])
        nq = nev(np.clip(q[::-1], 0.46, 1.6))
        mc[k] = dict(R_from_intensities=round(a / b, 3), R_table=Rt, R_16_50_84=[round(float(x), 3) for x in q],
                     ne_at_R_84_50_16=[round(float(x), 1) for x in nq], frac_above_LDL=round(float(np.mean(R >= Rl)), 3))
    res['table2_mc_1e4K_TZ10'] = mc
    # Can per-fibre noise turn region ratios of 1.03-1.19 into a per-fibre median n_e of 64 (paper eq. 5)?
    # Fibres with R outside eq. 5's monotonic range (R >= 1.435, i.e. above the low-density limit) are dropped, as a
    # density cannot be computed for them; the median of the remaining eq.-5 densities is returned.
    sim = {}
    for sigR in (0.03, 0.15, 0.30):
        row = {}
        for Rtrue in (1.03, 1.10, 1.19, 1.30, 1.34, 1.38, 1.42):
            Rf = rng.normal(Rtrue, sigR, 20000)
            ok = (Rf > 0.45) & (Rf < 1.435)
            ne_f = proxauf(Rf[ok])
            row[f'{Rtrue:.2f}'] = dict(median_ne=round(float(np.median(ne_f)), 1), frac_dropped=round(float(1 - ok.mean()), 3),
                                       mean_R_kept=round(float(Rf[ok].mean()), 3))
        sim[f'sigmaR_{sigR}'] = row
    res['per_fibre_median_sim_eq5'] = sim
    out = ROOT / 'data/derived/verify_paper_sii.json'
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    make_fig(res)


def make_fig(res):
    from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, SII_C, NII_C, FGDIM, FG, SUN, WARN, OK, ACCENT
    import matplotlib.pyplot as plt
    fig, ax = new_sky_figure(9, 5.6, '[S II] 6716/6731 density diagnostic',
                             'own 5-level S+ atom (TZ10 Ω; TZ10 / MZ82 A) vs paper eq. 5 (Proxauf+14) · Table 2 ratios')
    ne = np.logspace(-0.5, 4.5, 300)
    for (c, a, T, col, ls) in [('TZ10', 'TZ10', 1e4, ACCENT, '-'), ('TZ10', 'TZ10', 8000, ACCENT, '--'),
                               ('TZ10', 'MZ82', 1e4, OK, '-'), ('RBS96', 'MZ82', 1e4, WARN, ':')]:
        ax.plot(ne, [ratio(n, T, coll=c, aval=a) for n in ne], color=col, ls=ls, lw=1.4,
                label=f'Ω {c} + A {a}, {T/1e3:.0f}k K')
    rr = np.linspace(1.0, 1.44, 400)
    ax.plot(proxauf(rr), rr, color=SUN, lw=1.6, label='paper eq. 5 (Proxauf+14)')
    ax.axvspan(64, 150, color=HA_C, alpha=0.18, lw=0)
    label(ax, 68, 0.62, 'paper\n64–150', color=HA_C, size=8)
    for R, k in [(1.10, 'P1'), (1.14, 'P2'), (1.19, 'P3'), (1.03, 'P4')]:
        ax.axhline(R, color=SII_C, lw=0.7, alpha=0.6)
        label(ax, 1.3, R + 0.008, f'{k} R={R:.2f}', color=SII_C, size=7.5)
    ax.set_xscale('log'); ax.set_xlim(0.5, 3e4); ax.set_ylim(0.4, 1.55)
    ax.set_xlabel('n_e (cm$^{-3}$)'); ax.set_ylabel('I(6716)/I(6731)')
    ax.legend(loc='lower left', fontsize=7.5, ncol=1)
    ax.grid(alpha=0.3)
    finish(fig, str(ROOT / 'figures/verify_paper_sii.png'))


if __name__ == '__main__':
    main()
