"""Fit and plot the LAMOST LRS / MRS residual stacks (topic spectra_sdss_lamost).

Inputs : data/derived/spectra_sdss_lamost_{lrs,mrs}_resid_stacks.json, ..._resid_stars.csv
Outputs: data/derived/spectra_sdss_lamost_lamost_stackfit.json
         figures/spectra_sdss_lamost_lamost_stacks.png
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, NII_C, SII_C, FGDIM, FG, WARN

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
R = {'lrs': 1800.0, 'mrs': 7500.0}
LINESETS = {'OIII': ['OIII5007'], 'Hb': ['Hb'], 'Ha': ['NII6548', 'Ha', 'NII6583'], 'SII': ['SII6716', 'SII6731']}


def model_factory(lines, sig_A, free_v=True):
    def f(x, v, *amps):
        y = np.zeros_like(x)
        for ln, a in zip(lines, amps):
            c = VAC[ln] * (1 + v / C_KMS)
            y += a * np.exp(-0.5 * ((x - c) / sig_A) ** 2) / (np.sqrt(2 * np.pi) * sig_A)
        return y
    return f


def fit_stack(g, m, e, lines, res, v0=20.0):
    sig_A = np.mean(g) / res / 2.3548
    sig_A = np.hypot(sig_A, np.mean(g) * 19.0 / C_KMS)
    ok = np.isfinite(m) & np.isfinite(e) & (e > 0)
    # local baseline: subtract median outside lines
    far = ok.copy()
    for ln in lines:
        far &= np.abs(g - VAC[ln]) > 5 * sig_A
    base = np.nanmedian(m[far]) if far.sum() > 5 else 0.0
    y = m - base
    f = model_factory(lines, sig_A)
    best = None
    for vs in (-40.0, 0.0, 20.0, 40.0):
        try:
            p, c = curve_fit(f, g[ok], y[ok], p0=[vs] + [1.0] * len(lines), sigma=e[ok], absolute_sigma=True, maxfev=20000)
            chi = np.sum(((y[ok] - f(g[ok], *p)) / e[ok]) ** 2)
            if best is None or chi < best[2]:
                best = (p, c, chi)
        except RuntimeError:
            pass
    p, c, chi = best
    dof = ok.sum() - len(p)
    sc = np.sqrt(max(chi / dof, 1))
    out = dict(v_helio=float(p[0]), v_err=float(np.sqrt(c[0, 0]) * sc), chi2_red=float(chi / dof), baseline=float(base))
    for ln, a, i in zip(lines, p[1:], range(1, len(p))):
        out['F_' + ln] = float(a)
        out['F_' + ln + '_err'] = float(np.sqrt(c[i, i]) * sc)
    # fixed-velocity amplitude (v = +20 km/s) for comparison
    p2, c2 = curve_fit(lambda x, *a: f(x, v0, *a), g[ok], y[ok], p0=[1.0] * len(lines), sigma=e[ok], absolute_sigma=True)
    chi2 = np.sum(((y[ok] - f(g[ok], v0, *p2)) / e[ok]) ** 2)
    sc2 = np.sqrt(max(chi2 / (ok.sum() - len(p2)), 1))
    for ln, a, i in zip(lines, p2, range(len(p2))):
        out['Ffix_' + ln] = float(a)
        out['Ffix_' + ln + '_err'] = float(np.sqrt(c2[i, i]) * sc2)
    free = lines == ['OIII5007']
    return out, (g, y, e, f(g, *p) if free else f(g, v0, *p2))


def wmean(v, e):
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    v, e = v[ok], e[ok]
    if len(v) < 5:
        return np.nan, np.nan, 0
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med))
    k = np.abs(v - med) < 5 * mad
    w = 1 / e[k] ** 2
    mu = np.sum(w * v[k]) / np.sum(w)
    chi = np.sum(w * (v[k] - mu) ** 2) / max(k.sum() - 1, 1)
    return float(mu), float(np.sqrt(1 / np.sum(w)) * np.sqrt(max(chi, 1))), int(k.sum())


def main():
    out = {}
    curves = {}
    for mode in ('lrs', 'mrs'):
        sp = ROOT / f'data/derived/spectra_sdss_lamost_{mode}_resid_stacks.json'
        if not sp.exists():
            continue
        st = json.loads(sp.read_text())
        T = Table.read(ROOT / f'data/derived/spectra_sdss_lamost_{mode}_resid_stars.csv')
        out[mode] = {'per_star_weighted': {}, 'stack_fits': {}}
        for col in [c for c in T.colnames if c.startswith('I_') and not c.endswith('_err')]:
            ln = col[2:]
            on = T['sample'] == 'on'
            a = wmean(np.array(T[col][on]), np.array(T[col + '_err'][on]))
            b = wmean(np.array(T[col][~on]), np.array(T[col + '_err'][~on]))
            out[mode]['per_star_weighted'][ln] = dict(on=a[0], on_err=a[1], n_on=a[2], off=b[0], off_err=b[1], n_off=b[2],
                                                      diff=a[0] - b[0], diff_err=float(np.hypot(a[1], b[1])))
        for win, lines in LINESETS.items():
            if f'{win}_on' not in st:
                continue
            g = np.array(st[f'{win}_on']['grid'])
            mon, eon = np.array(st[f'{win}_on']['mean']), np.array(st[f'{win}_on']['err'])
            mof, eof = np.array(st[f'{win}_off']['mean']), np.array(st[f'{win}_off']['err'])
            # restrict to the neighbourhood of the lines
            lo = min(VAC[l] for l in lines) - (12 if mode == 'lrs' else 4)
            hi = max(VAC[l] for l in lines) + (12 if mode == 'lrs' else 4)
            m = (g > lo) & (g < hi)
            for tag, mm, ee in [('on', mon, eon), ('off', mof, eof), ('on-off', mon - mof, np.hypot(eon, eof))]:
                try:
                    r, cur = fit_stack(g[m], mm[m], ee[m], lines, R[mode])
                except Exception as ex:          # noqa
                    r, cur = {'error': str(ex)}, None
                out[mode]['stack_fits'][f'{win}_{tag}'] = r
                curves[(mode, win, tag)] = cur
    # MRS [OIII]: narrow-window Gaussian (free centre and width) + linear baseline, per plate and combined.
    # OFF-strand stars carry a NEGATIVE residual (super-sky contains strand [OIII]); its centre gives the
    # heliocentric velocity of the [OIII] in the sky fibres.
    sp = ROOT / 'data/derived/spectra_sdss_lamost_mrs_resid_stacks.json'
    if sp.exists():
        st = json.loads(sp.read_text())
        out.setdefault('mrs', {})['oiii_dip_fits'] = {}
        lc0 = VAC['OIII5007']

        def gmod(x, a, c, s, b0, b1):
            return a * np.exp(-0.5 * ((x - c) / s) ** 2) / (np.sqrt(2 * np.pi) * s) + b0 + b1 * (x - lc0)
        for key in [k for k in st if k.startswith('OIII_')]:
            g = np.array(st[key]['grid']); m = np.array(st[key]['mean']); e = np.array(st[key]['err'])
            w = (g > 5006.9) & (g < 5010.3) & np.isfinite(m) & (e > 0)
            best = None
            for a0 in (-1.0, 1.0):
                try:
                    p, c = curve_fit(gmod, g[w], m[w], p0=[a0, lc0 + 0.3, 0.4, 0, 0], sigma=e[w], absolute_sigma=True,
                                     bounds=([-50, 5007.6, 0.2, -10, -10], [50, 5009.4, 1.2, 10, 10]), maxfev=20000)
                    chi = float(np.sum(((m[w] - gmod(g[w], *p)) / e[w]) ** 2))
                    if best is None or chi < best[2]:
                        best = (p, c, chi)
                except (RuntimeError, ValueError):
                    pass
            if best is None:
                continue
            p, c, chi = best
            dof = int(w.sum()) - 5
            sc = np.sqrt(max(chi / dof, 1))
            out['mrs']['oiii_dip_fits'][key] = dict(
                n_stars=st[key]['n'], F=float(p[0]), F_err=float(np.sqrt(c[0, 0]) * sc),
                centre_A=float(p[1]), v_helio=float((p[1] / lc0 - 1) * C_KMS), v_err=float(np.sqrt(c[1, 1]) * sc / lc0 * C_KMS),
                sigma_A=float(p[2]), sigma_kms=float(p[2] / lc0 * C_KMS), chi2_red=float(chi / dof),
                snr=float(abs(p[0]) / (np.sqrt(c[0, 0]) * sc)))
            curves[('mrsdip', key)] = (g[w], m[w], e[w], gmod(g[w], *p))
    (ROOT / 'data/derived/spectra_sdss_lamost_lamost_stackfit.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:6000])

    # figure: rows = LRS [OIII], LRS Halpha/[NII], MRS [OIII], MRS Halpha/[NII]; on (colour), off (grey), on-off
    panels = [('lrs', 'OIII', OIII_C, 'LRS  [O III] 5007'), ('lrs', 'Hb', HB_C, 'LRS  Hβ'),
              ('lrs', 'Ha', NII_C, 'LRS  Hα + [N II]'), ('lrs', 'SII', SII_C, 'LRS  [S II]'),
              ('mrs', 'OIII', OIII_C, 'MRS  [O III] 5007'), ('mrs', 'Ha', NII_C, 'MRS  Hα + [N II]')]
    fig, axes = new_sky_figure(12, 9.5, 'LAMOST DR12 · nebular residuals in stellar spectra on the Silk strand',
                               'stack of star residuals (heliocentric, vacuum), R/Å; ON |d|<0.2–0.25°, OFF 0.6–2.5°',
                               nrows=3, ncols=2)
    for ax, (mode, win, col, ttl) in zip(axes.ravel(), panels):
        for tag, c, lw, z in [('off', FGDIM, 0.9, 1), ('on', col, 1.2, 2)]:
            cur = curves.get((mode, win, tag))
            if cur is None:
                continue
            g, y, e, mod = cur
            ax.fill_between(g, y - e, y + e, color=c, alpha=0.18, lw=0, zorder=z)
            ax.plot(g, y, color=c, lw=lw, label=tag.upper(), zorder=z + 1)
        cur = curves.get((mode, win, 'on-off'))
        if cur is not None and mode == 'mrs':
            dip = curves.get(('mrsdip', f'{win}_off'))
            if dip is not None:
                ax.plot(dip[0], dip[3], color=WARN, lw=1.1, ls='--', label='fit OFF dip', zorder=4)
                r = out['mrs']['oiii_dip_fits'][f'{win}_off']
                label(ax, 0.02, 0.06, f"OFF dip: v_helio = {r['v_helio']:+.1f}±{r['v_err']:.1f} km/s, "
                      f"σ = {r['sigma_kms']:.0f} km/s", color=WARN, size=7.5, transform=ax.transAxes)
            else:
                label(ax, 0.02, 0.06, 'sky-subtraction mismatch dominates; not fitted', color=WARN, size=7,
                      transform=ax.transAxes)
            cur = None
        if cur is not None:
            g, y, e, mod = cur
            ax.plot(g, mod, color=WARN, lw=1.0, ls='--', label='fit ON−OFF', zorder=4)
            r = out[mode]['stack_fits'][f'{win}_on-off']
            key = [k for k in r if k.startswith('F_') and not k.endswith('_err')]
            if mode == 'mrs':
                txt = 'derivative-shaped: sky-subtraction mismatch, no velocity'
            elif win == 'OIII':
                txt = f"ON−OFF {r['F_OIII5007']:.2f}±{r['F_OIII5007_err']:.2f} R, v = {r['v_helio']:+.0f}±{r['v_err']:.0f} km/s"
            else:
                ks = [k for k in r if k.startswith('Ffix_') and not k.endswith('_err')]
                txt = 'ON−OFF ' + ', '.join(f"{k[5:]} {r[k]:.2f}±{r[k + '_err']:.2f}" for k in ks) + ' R'
            label(ax, 0.02, 0.06, txt, color=WARN, size=7.5, transform=ax.transAxes)
        for ln in LINESETS[win]:
            ax.axvline(VAC[ln] * (1 + 20 / C_KMS), color=FGDIM, lw=0.5, ls=':')
        ax.axhline(0, color=FGDIM, lw=0.5)
        label(ax, 0.02, 0.88, ttl, color=col, size=9, transform=ax.transAxes)
        ax.set_ylabel('R / Å')
        ax.grid(alpha=0.2)
        ax.legend(fontsize=6.5, loc='upper right')
    for ax in axes[-1]:
        ax.set_xlabel('vacuum wavelength, heliocentric (Å)')
    finish(fig, ROOT / 'figures/spectra_sdss_lamost_lamost_stacks.png')


if __name__ == '__main__':
    main()
