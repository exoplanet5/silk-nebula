"""verify_paper_kd.py -- independent kinematic-distance check of Wu et al. 2026 Table 4 (tag: paper).

Own implementation (not a re-run of wave-1 paper_review_lib):
  * Persic+1996 universal rotation curve with the Reid+2014 A5 parameters as tabulated in the
    kd package source (a1=241, a2=0.90, a3=1.46, R0=8.34; read, not executed, from
    https://github.com/tvwenger/kd  kd/reid14_rotcurve.py) and a flat curve.
  * Two projections: 'kd' (V = R0 sin l (Theta/R - Theta0/R0), no cos b, as in kd) and
    '3d' (x cos b, the physically correct projection for a circular orbit at height z).
  * IAU standard solar motion (U,V,W) = (10.27, 15.32, 7.74) km/s (kd_utils.py defaults).
  * Traditional distance: 1-pc grid 0.001-30 kpc, outer-Galaxy branch, tolerance 0.1 km/s.
  * Monte Carlo distance mimicking kd.pdf_kd: resample a1,a2,a3,R0 and V (Gaussian), grid floor
    0.01 kpc, unsolved samples dropped, boundary-corrected (linear-combination) Gaussian KDE with
    a lower bound at 0, mode + 68.3 % minimum-width interval on a 1000-point grid.
Outputs data/derived/verify_paper_kd.json and figures/verify_paper_kd.png.
Run: ~/.venvs/astro313/bin/python -I scripts/verify_paper_kd.py
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
rng = np.random.default_rng(20261117)

# ---------------- paper inputs (Wu+2026 Tables 1 and 4; arXiv v1 == AJ text) ----------------
REG = {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72), 'P4': (115.20, 16.18)}
# line: (V_helio, err, V_LSR, MC mode, MC -err, MC +err, traditional)   [distances as printed]
T4 = {
 'P1': {'OIII': (26.8, 4.3, 15.5, 1.63, 0.84, 1.01, 1.79), 'SII': (11.8, 2.5, -1.0, 0.06, 0.05, 1.17, 0.16),
        'NII': (11.8, 1.6, -0.9, 0.11, 0.10, 1.12, 0.18), 'Ha': (16.3, 1.2, 4.7, 0.47, 0.45, 1.18, 0.61)},
 'P2': {'OIII': (26.8, 4.1, 14.5, 1.58, 0.83, 1.23, 1.85), 'SII': (12.0, 1.6, 0.5, 0.21, 0.20, 1.06, 0.17),
        'NII': (12.0, 1.6, -0.5, 0.08, 0.07, 1.35, 0.17), 'Ha': (17.0, 1.2, 4.0, 0.31, 0.31, 1.44, 0.54)},
 'P3': {'OIII': (23.8, 3.3, 11.2, 0.98, 0.69, 1.08, 1.20), 'SII': (11.7, 1.9, 1.7, 0.21, 0.21, 1.12, 0.23),
        'NII': (11.7, 1.9, 0.4, 0.09, 0.08, 1.02, 0.14), 'Ha': (18.4, 1.4, 6.1, 0.34, 0.33, 1.20, 0.58)},
 'P4': {'SII': (12.7, 2.1, 1.1, 0.22, 0.21, 1.61, 0.22), 'NII': (13.3, 2.3, 1.0, 0.28, 0.27, 0.93, 0.49),
        'Ha': (19.7, 1.4, 4.5, 0.54, 0.52, 1.18, 0.80)},
}

# ---------------- rotation curves (own code) ----------------
R14 = dict(a1=241.0, a2=0.90, a3=1.46, R0=8.34)
R14_ERR = dict(a1=8.0, a2=0.06, a3=0.16, R0=0.16)
USTD, VSTD, WSTD = 10.27, 15.32, 7.74


def theta_urc(R, a1, a2, a3, R0):
    """Persic+96 URC (eqs 8-11 of Persic, Salucci & Stel 1996) in the Reid+14 parametrisation."""
    x = R / (a2 * R0)
    beta = 0.72 + 0.44 * np.log10((a3 / 1.5) ** 5)
    vd2 = beta * 1.97 * x ** 1.22 / (x ** 2 + 0.78 ** 2) ** 1.43
    vh2 = (1 - beta) * (1 + a3 ** 2) * x ** 2 / (x ** 2 + a3 ** 2)
    tot = vd2 + vh2
    return a1 * np.sqrt(np.where(tot > 0, tot, np.nan))


def rgal(l, b, d, R0):
    dc = d * np.cos(np.radians(b))
    return np.sqrt(R0 ** 2 + dc ** 2 - 2 * R0 * dc * np.cos(np.radians(l)))


def vlsr_model(l, b, d, curve='urc', proj='kd', a1=241.0, a2=0.90, a3=1.46, R0=8.34, theta0=None):
    R = np.maximum(rgal(l, b, d, R0), 1e-6)
    if curve == 'urc':
        th, th0 = theta_urc(R, a1, a2, a3, R0), theta_urc(R0, a1, a2, a3, R0)
    else:  # flat
        th0 = theta0 if theta0 is not None else theta_urc(R0, a1, a2, a3, R0)
        th = th0
    v = R0 * np.sin(np.radians(l)) * (th / R - th0 / R0)
    if proj == '3d':
        v = v * np.cos(np.radians(b))
    return v


def helio_to_lsr(l, b, vh):
    lr, br = np.radians(l), np.radians(b)
    return vh + (USTD * np.cos(lr) + VSTD * np.sin(lr)) * np.cos(br) + WSTD * np.sin(br)


GRID = np.arange(0.001, 30.0 + 0.001, 0.001)


def trad_kd(l, b, v, tol=0.1, **kw):
    vg = vlsr_model(l, b, GRID, **kw)
    i = np.nanargmin(np.abs(vg - v))
    return GRID[i] if abs(vg[i] - v) < tol else np.nan


def lc_kde(samples, x, lower=0.0):
    """Gaussian KDE with linear-combination boundary correction at `lower` (Jones 1993),
    bandwidth h = sigma n^-1/5 (Scott covariance, as pyqt_fit's default)."""
    s = samples[np.isfinite(samples)]
    n = s.size
    h = np.std(s, ddof=1) * n ** (-0.2)
    uu = (x - lower) / h
    a0, a1 = norm.cdf(uu), -norm.pdf(uu)
    a2 = a0 - uu * norm.pdf(uu)
    out = np.empty_like(x)
    for j0 in range(0, x.size, 100):
        sl = slice(j0, j0 + 100)
        z = (x[sl, None] - s[None, :]) / h
        k = norm.pdf(z)
        num = ((a2[sl, None] - a1[sl, None] * z) * k).sum(1)
        out[sl] = num / (n * h * (a0[sl] * a2[sl] - a1[sl] ** 2))
    return np.clip(out, 0, None)


def hpd(samples, alpha=0.683, nb=1000):
    s = samples[np.isfinite(samples)]
    x = np.linspace(s.min(), s.max(), nb)
    p = lc_kde(s, x)
    mode = x[np.argmax(p)]
    order = np.argsort(p)[::-1]
    c = np.cumsum(p[order] / p.sum())
    k = np.searchsorted(c, alpha) + 1
    sel = x[order[:k]]
    return mode, sel.min(), sel.max()


def mc_kd(l, b, v, verr, n=10000, dmin=0.01, dmax=12.0, tol=0.1, proj='kd'):
    """kd.pdf_kd emulation; grid truncated at 12 kpc (V(d) is monotonic here and all solutions < 6 kpc)."""
    g = np.arange(dmin, dmax + 0.001, 0.001)
    a1 = rng.normal(R14['a1'], R14_ERR['a1'], n); a2 = rng.normal(R14['a2'], R14_ERR['a2'], n)
    a3 = rng.normal(R14['a3'], R14_ERR['a3'], n); R0 = rng.normal(R14['R0'], R14_ERR['R0'], n)
    vs = rng.normal(v, verr, n)
    out = np.full(n, np.nan)
    for j0 in range(0, n, 400):
        sl = slice(j0, j0 + 400)
        vg = vlsr_model(l, b, g[None, :], proj=proj, a1=a1[sl, None], a2=a2[sl, None], a3=a3[sl, None], R0=R0[sl, None])
        dv = np.abs(vg - vs[sl, None])
        i = np.nanargmin(dv, axis=1)
        ok = dv[np.arange(i.size), i] < tol
        out[sl] = np.where(ok, g[i], np.nan)
    return out


def main():
    res = {'inputs': 'Wu+2026 Table 1 positions; Table 4 velocities and distances as printed (AJ HTML == arXiv v1)',
           'rotation_curves': {'R14_URC': R14, 'R14_URC_err': R14_ERR, 'solar_motion_std': [USTD, VSTD, WSTD]}}
    lb = {}
    for k, (ra, dec) in REG.items():
        c = SkyCoord(ra * u.deg, dec * u.deg).galactic
        lb[k] = (c.l.deg, c.b.deg)
    res['lb'] = lb

    # 1) rotation-curve velocity at key distances and local gradient
    keyd = [0.16, 0.28, 0.30, 1.0, 1.6]
    vtab = {}
    for k in ['P1', 'P2', 'P3']:
        l, b = lb[k]
        row = {}
        for name, kw in {'URC_kd': dict(curve='urc', proj='kd'), 'URC_3d': dict(curve='urc', proj='3d'),
                         'flat_kd': dict(curve='flat', proj='kd'), 'flat_3d': dict(curve='flat', proj='3d'),
                         'flat_R0_8.15_T236_3d': dict(curve='flat', proj='3d', R0=8.15, theta0=236.0)}.items():
            row[name] = {f'{d:.2f}kpc': round(float(vlsr_model(l, b, d, **kw)), 2) for d in keyd}
            row[name]['grad_0_2kpc'] = round(float(vlsr_model(l, b, 2.0, **kw) / 2.0), 2)
        vtab[k] = row
    res['vlsr_vs_distance'] = vtab

    # 2) traditional distances vs Table 4, both unit hypotheses
    trad = {}
    chi = {'kpc': [], '100pc': []}
    chiHa = {'kpc_Vtable': [], '100pc_Vtable': [], 'kpc_Vhelio': [], '100pc_Vhelio': []}
    for k, lines in T4.items():
        l, b = lb[k]
        trad[k] = {}
        for ln, (vh, ve, vl, mcm, mcn, mcp, tr) in lines.items():
            vl_from_vh = helio_to_lsr(l, b, vh)
            d_tab_v = {name: trad_kd(l, b, vl, **kw) for name, kw in
                       {'URC_kd': {}, 'URC_3d': dict(proj='3d'), 'flat_3d': dict(curve='flat', proj='3d')}.items()}
            d_from_vh = trad_kd(l, b, vl_from_vh)
            v_need_kpc = float(vlsr_model(l, b, tr))           # V that kd would need for d = tr kpc
            v_need_100 = float(vlsr_model(l, b, tr / 10.0))    # ... for d = tr x 100 pc
            trad[k][ln] = dict(V_helio=vh, V_LSR_table=vl, V_LSR_from_helio_std=round(float(vl_from_vh), 2),
                               table_trad=tr, my_trad_from_table_V={kk: (None if np.isnan(vv) else round(float(vv), 3)) for kk, vv in d_tab_v.items()},
                               my_trad_from_Vhelio_URC_kd=None if np.isnan(d_from_vh) else round(float(d_from_vh), 3),
                               V_needed_if_kpc=round(v_need_kpc, 2), V_needed_if_100pc=round(v_need_100, 2))
            if ln == 'OIII':
                chi['kpc'].append((vl - v_need_kpc) / ve)
                chi['100pc'].append((vl - v_need_100) / ve)
            if ln == 'Ha':
                chiHa['kpc_Vtable'].append((vl - v_need_kpc) / ve)
                chiHa['100pc_Vtable'].append((vl - v_need_100) / ve)
                chiHa['kpc_Vhelio'].append((vl_from_vh - v_need_kpc) / ve)
                chiHa['100pc_Vhelio'].append((vl_from_vh - v_need_100) / ve)
    res['traditional'] = trad
    res['OIII_unit_test'] = {h: dict(z=[round(float(z), 2) for z in zz], chi2=round(float(np.sum(np.square(zz))), 2),
                                     p_value=float(1 - __import__('scipy').stats.chi2.cdf(np.sum(np.square(zz)), 3)))
                             for h, zz in chi.items()}

    res['Ha_unit_test'] = {h: dict(z=[round(float(z), 2) for z in zz], chi2=round(float(np.sum(np.square(zz))), 2),
                                   p_value=float(1 - __import__('scipy').stats.chi2.cdf(np.sum(np.square(zz)), 4)))
                           for h, zz in chiHa.items()}

    # 3) Monte Carlo (kd.pdf_kd emulation) for all Table 4 entries
    mc = {}
    samples_store = {}
    for k, lines in T4.items():
        l, b = lb[k]
        mc[k] = {}
        for ln, (vh, ve, vl, mcm, mcn, mcp, tr) in lines.items():
            s = mc_kd(l, b, vl, ve)
            frac_bad = float(np.mean(~np.isfinite(s)))
            mode, lo, hi = hpd(s)
            mc[k][ln] = dict(table=(mcm, mcm - mcn, mcm + mcp), table_lower_bound=round(mcm - mcn, 3),
                             mine=(round(float(mode), 3), round(float(lo), 3), round(float(hi), 3)),
                             frac_no_solution=round(frac_bad, 3), median_valid=round(float(np.nanmedian(s)), 3))
            samples_store[f'{k}_{ln}'] = s
    res['monte_carlo'] = mc
    lbs = [v['table_lower_bound'] for kk in mc.values() for v in kk.values() if v['table'][0] < 0.6]
    res['red_line_MC_lower_bounds'] = dict(values=lbs, note='kd.pdf_kd uses dist_min = 0.01 kpc; HPD lower bound = grid floor when the mode is near 0')

    # 4) per-fibre-median hypothesis for the Halpha traditional distances (digitized Fig 5, geometry tag)
    from matplotlib.path import Path as MPath
    tr = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
    dig = np.genfromtxt(ROOT / 'data/derived/geometry_lamost_halpha_vlsr_fig5_digitized.csv', delimiter=',', skip_header=1)
    fib = {}
    for k in REG:
        poly = np.array(tr['regions'][k]['drawn_in_fig1b']['outline'])
        inside = MPath(poly).contains_points(dig[:, :2])
        v = dig[inside, 2]; cen = dig[inside, 4]
        l, b = lb[k]
        vg = vlsr_model(l, b, GRID)
        d = np.interp(v, vg, GRID, left=np.nan)  # V(d) monotonic; v below V(1 pc) -> no solution
        d[cen > 0] = 9.0   # >=15 km/s: far beyond median anyway
        d[cen < 0] = np.nan  # <=0 km/s: no outer-Galaxy solution in kd
        fib[k] = dict(n=int(inside.sum()), n_censored_low=int((cen < 0).sum()), n_censored_high=int((cen > 0).sum()),
                      median_V=round(float(np.median(v)), 2), mean_V=round(float(np.mean(v)), 2),
                      median_d_valid=round(float(np.nanmedian(d)), 3), frac_valid=round(float(np.mean(np.isfinite(d))), 3),
                      table_Ha_trad=T4[k]['Ha'][6], table_Ha_VLSR=T4[k]['Ha'][2],
                      d_at_median_V=round(float(np.interp(np.median(v), vg, GRID, left=np.nan)), 3))
    res['per_fibre_Ha_test'] = fib
    res['per_fibre_note'] = ('Digitized Fig 5 colours are clipped at 0 and 15 km/s; fibres at <=0 km/s have no outer-Galaxy solution '
                             'in kd and are dropped, so a per-fibre median distance is biased high relative to d(median V).')

    out = ROOT / 'data/derived/verify_paper_kd.json'
    out.write_text(json.dumps(res, indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item()))
    np.savez_compressed(ROOT / 'data/derived/verify_paper_kd_mc_samples.npz', **samples_store)
    print(json.dumps({k: res[k] for k in ['lb', 'OIII_unit_test', 'Ha_unit_test', 'red_line_MC_lower_bounds', 'per_fibre_Ha_test']}, indent=1))
    for k in trad:
        for ln in trad[k]:
            t = trad[k][ln]; m = mc[k][ln]
            print(k, ln, 'Vtab', t['V_LSR_table'], 'Vh->lsr', t['V_LSR_from_helio_std'], 'trad tab', t['table_trad'],
                  'mine', t['my_trad_from_table_V'], 'fromVh', t['my_trad_from_Vhelio_URC_kd'],
                  '| Vneed kpc', t['V_needed_if_kpc'], '100pc', t['V_needed_if_100pc'],
                  '| MC tab', m['table'], 'mine', m['mine'], 'bad', m['frac_no_solution'])
    print(json.dumps(vtab, indent=0)[:1500])
    make_fig(lb, res, samples_store)


def make_fig(lb, res, samples):
    from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, SII_C, NII_C, FGDIM, FG, SUN, WARN, OK
    import matplotlib.pyplot as plt
    fig, axs = new_sky_figure(13, 5.2, 'Wu+2026 Table 4 · kinematic-distance units',
                              'own Reid+14 URC (kd form) · P1 l=203.8 b=19.5 · points: Table 4 distances read as kpc (filled) or x100 pc (open)',
                              ncols=2)
    ax = axs[0]
    l, b = lb['P1']
    d = np.linspace(0.001, 3.0, 600)
    for name, kw, ls in [('URC (kd form)', {}, '-'), ('URC x cos b', dict(proj='3d'), '--'), ('flat 240', dict(curve='flat', proj='3d'), ':')]:
        ax.plot(d, vlsr_model(l, b, d, **kw), color=FG, ls=ls, lw=1.2, label=name)
    ax.fill_between(d, vlsr_model(l, b, d) - 7, vlsr_model(l, b, d) + 7, color=FGDIM, alpha=0.15, lw=0, label='±7 km/s cloud motion')
    cols = {'OIII': OIII_C, 'Ha': HA_C, 'SII': SII_C, 'NII': NII_C}
    done = set()
    for k, lines in T4.items():
        for ln, (vh, ve, vl, mcm, mcn, mcp, tr) in lines.items():
            c = cols[ln]
            lab = ln if ln not in done else None
            done.add(ln)
            ax.errorbar(tr, vl, yerr=ve, fmt='o', color=c, ms=5, capsize=0, label=lab)
            ax.errorbar(tr / 10, vl, yerr=ve, fmt='o', mfc='none', mec=c, color=c, ms=5, capsize=0, alpha=0.8)
    ax.set_xlim(0, 2.6); ax.set_ylim(-6, 22)
    ax.set_xlabel('distance (kpc)'); ax.set_ylabel('V_LSR (km/s)')
    ax.axhline(0, color='#2a323c', lw=0.8)
    label(ax, 0.27, 17.0, 'open = x100 pc reading:\n[O III] would need V≈1.8', color=WARN, size=8)
    label(ax, 1.3, 3.5, 'filled = kpc reading:\n[O III] on the curve (chi2 0.05)', color=OK, size=8)
    ax.legend(loc='lower right', ncol=2, fontsize=7.5)
    ax.grid(alpha=0.3)
    # panel b: MC PDFs P1 OIII and lower-bound signature
    ax = axs[1]
    for key, c in [('P1_OIII', OIII_C), ('P3_OIII', '#9be7ff'), ('P1_SII', SII_C), ('P4_Ha', HA_C)]:
        s = samples[key]; s = s[np.isfinite(s)]
        x = np.linspace(0, 4, 400)
        ax.plot(x, lc_kde(s, x), color=c, lw=1.4, label=key.replace('_', ' '))
        k, ln = key.split('_')
        mcm = T4[k][ln][3]
        ax.axvline(mcm, color=c, ls=':', lw=1)
    ax.axvline(0.01, color=SUN, lw=1, ls='-')
    label(ax, 0.05, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 1, 'kd grid floor 0.01 kpc', color=SUN, size=8)
    label(ax, 2.2, ax.get_ylim()[1] * 0.6, 'solid: own kd.pdf_kd emulation\ndotted: Table 4 MC modes as kpc', color=FGDIM, size=8)
    ax.set_xlim(0, 4); ax.set_xlabel('Monte Carlo kinematic distance (kpc)'); ax.set_ylabel('PDF')
    ax.legend(loc='upper right', fontsize=7.5, ncol=2)
    ax.grid(alpha=0.3)
    finish(fig, str(ROOT / 'figures/verify_paper_kd.png'))


if __name__ == '__main__':
    main()
