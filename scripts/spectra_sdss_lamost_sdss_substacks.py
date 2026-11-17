"""Regional ON-minus-OFF stacks of SDSS/BOSS sky-fibre spectra with a wavelength zero-point check
(topic tag spectra_sdss_lamost).

What it adds to spectra_sdss_lamost_sdss_stack.py
  (1) Wavelength zero point.  In every plate-mjd-spectrograph group the mean TOTAL sky (flux+sky) of the
      sky fibres is fitted around the airglow lines [O I] 5577 (blue camera) and [O I] 6300, 6364 (red
      camera).  These are at rest in the topocentric frame, so in the heliocentric spec-lite frame they must
      sit at lambda_vac * (1 + v_corr/c) = lambda_vac * (1 - HELIO_RV/c) (SDSS HELIO_RV = -v_corr).  The offset dv_blue / dv_red (km/s) is the residual wavelength
      error of the pipeline at that epoch; the line velocities of each stack are corrected by the
      ON-fibre-weighted mean offset of the groups that enter it.
  (2) Stacks per region: all, south (Dec < 20, LAMOST MRS-N field / paper regions), north (Dec >= 20),
      and E strand vs W strand (Dec 15-20).  For each: [O III] 4959/5007, Hbeta, Halpha, [N II], [S II]
      intensities (diffuse-corrected Rayleigh), free-velocity fits for [O III] and Halpha+[N II], ratios and
      2-sigma limits.
  (3) Velocity null test: the same stacking with ON/OFF swapped for fibres 1.0-1.15 deg from the strand
      (should give no line).

Inputs : data/derived/spectra_sdss_lamost_sdss_skyfibres.npz, ..._sdss_groups.csv, data/derived/silk_trace.json
Outputs: data/derived/spectra_sdss_lamost_sdss_substacks.json, ..._sdss_zeropoints.csv,
         figures/spectra_sdss_lamost_sdss_substacks.png
usage  : python -I spectra_sdss_lamost_sdss_substacks.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS, air_to_vac, rayleigh_per_fibre_flux, FIBRE_DIAM
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, NII_C, SII_C, FGDIM, WARN, FG
import spectra_sdss_lamost_trace as TR

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
ON, OFF = 0.15, (0.45, 2.0)
AIRGLOW = {'OI5577': (float(air_to_vac(5577.339)), 'blue'), 'OI6300': (float(air_to_vac(6300.304)), 'red'),
           'OI6364': (float(air_to_vac(6363.776)), 'red')}
SETS = [('OIII', ['OIII4959', 'OIII5007'], 4940, 5030, 'blue'),
        ('Hb', ['Hb'], 4835, 4890, 'blue'),
        ('HaNII', ['NII6548', 'Ha', 'NII6583'], 6530, 6600, 'red'),
        ('SII', ['SII6716', 'SII6731'], 6700, 6750, 'red')]


def load():
    with np.load(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz') as z:
        d = {k: z[k] for k in z.files}
    G = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_groups.csv')
    return d, G


def gfit(x, y, e, c0, sig0):
    f = lambda xx, a, c, s, b0, b1: a * np.exp(-0.5 * ((xx - c) / s) ** 2) + b0 + b1 * (xx - c0)
    p, cv = curve_fit(f, x, y, p0=[np.nanmax(y) - np.nanmedian(y), c0, sig0, np.nanmedian(y), 0], sigma=e,
                      absolute_sigma=True, maxfev=20000)
    chi = np.sum(((y - f(x, *p)) / e) ** 2) / max(len(x) - 5, 1)
    return p, np.sqrt(np.diag(cv)) * np.sqrt(max(chi, 1))


def zeropoints(d, keys_arr, groups, lam):
    """Airglow-line centroid offsets (km/s) per group."""
    tot = d['flux'] + d['sky']          # not ivar-masked: SDSS-I flags bright sky lines (and_mask) at 5577/6300
    rows = []
    for k in groups:
        m = keys_arr == k
        t = np.where(d['wdisp'][m] > 0, tot[m], np.nan)
        mean = np.nanmedian(t, 0)
        err = 1.2533 * np.nanstd(t, 0) / np.sqrt(np.maximum(np.sum(np.isfinite(t), 0), 1))
        helio = float(np.nanmedian(d['helio_rv'][m]))
        r = dict(group=k, helio_rv=helio, n=int(m.sum()))
        for nm, (l0, cam) in AIRGLOW.items():
            lc = l0 * (1 - helio / C_KMS)    # topocentric rest; SDSS HELIO_RV = -v_corr (astropy), see sdss_fit TOPO note
            w = (lam > lc - 9) & (lam < lc + 9) & np.isfinite(mean) & (err > 0)
            if w.sum() < 12:
                r[nm] = np.nan; r[nm + '_err'] = np.nan
                continue
            try:
                p, pe = gfit(lam[w], mean[w], np.maximum(err[w], 1e-3), lc, 1.6)
                r[nm] = float((p[1] / lc - 1) * C_KMS)
                r[nm + '_err'] = float(pe[1] / lc * C_KMS)
            except (RuntimeError, ValueError):
                r[nm] = np.nan; r[nm + '_err'] = np.nan
        rows.append(r)
    return rows


def stack(d, lam, sel_groups, keys_arr, on_mask, off_mask, fpsf):
    diffs, wts, non, used = [], [], 0, []
    for k in sel_groups:
        m = keys_arr == k
        on, off = m & on_mask, m & off_mask
        if on.sum() < 1 or off.sum() < 4:
            continue
        diam = FIBRE_DIAM[str(d['instrument'][np.where(m)[0][0]])]
        conv = 1e-17 / rayleigh_per_fibre_flux(lam, diam) * fpsf[k]

        def wmean(idx):
            w = np.where(d['ivar'][idx] > 0, d['ivar'][idx], 0)
            return np.sum(d['flux'][idx] * w, 0) / np.maximum(np.sum(w, 0), 1e-30), np.sum(w, 0)
        fon, won = wmean(on)
        foff, woff = wmean(off)
        good = (won > 0) & (woff > 0)
        var = (1 / np.maximum(won, 1e-30) + 1 / np.maximum(woff, 1e-30)) * conv ** 2
        diffs.append(np.where(good, (fon - foff) * conv, np.nan))
        wts.append(np.where(good, 1 / var, 0))
        non += int(on.sum())
        used.append((k, int(on.sum())))
    if not diffs:
        return None
    D, W = np.array(diffs), np.array(wts)
    med = np.nanmedian(D, 0)
    sd = 1.4826 * np.nanmedian(np.abs(D - med), 0)
    W = np.where(np.isfinite(D) & (np.abs(D - med) < 5 * np.maximum(sd, 1e-6)), W, 0)
    st = np.nansum(np.nan_to_num(D) * W, 0) / np.maximum(W.sum(0), 1e-30)
    er = 1 / np.sqrt(np.maximum(W.sum(0), 1e-30))
    lf = ((lam > 4890) & (lam < 4940)) | ((lam > 5030) & (lam < 5065))
    rchi = np.nanmean(((st[lf] - np.nanmedian(st[lf])) / er[lf]) ** 2)
    er = er * np.sqrt(max(rchi, 1))
    return dict(stack=st, err=er, n_on=non, n_groups=len(used), used=used, err_inflation=float(np.sqrt(max(rchi, 1))))


def fit_set(lam, st, er, lines, lo, hi, sig_pix, free_v):
    m = (lam > lo) & (lam < hi) & (er > 0) & np.isfinite(st)
    x, y, e = lam[m], st[m], er[m]
    sig = np.hypot(sig_pix * x.mean(), x.mean() * 19.0 / C_KMS)

    def f(xx, v, b, *a):
        out = b + 0 * xx
        for ln, aa in zip(lines, a):
            c = VAC[ln] * (1 + v / C_KMS)
            out = out + aa * np.exp(-0.5 * ((xx - c) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
        return out
    res = {}
    p, c = curve_fit(lambda xx, b, *a: f(xx, 20.0, b, *a), x, y, p0=[0] + [0.3] * len(lines), sigma=e,
                     absolute_sigma=True)
    chi = np.sum(((y - f(x, 20.0, *p)) / e) ** 2) / max(len(x) - len(p), 1)
    for i, ln in enumerate(lines):
        res[ln] = dict(I=float(p[i + 1]), err=float(np.sqrt(c[i + 1, i + 1]) * np.sqrt(max(chi, 1))))
    model = f(x, 20.0, *p)
    if free_v:
        try:
            pv, cv = curve_fit(f, x, y, p0=[20.0] + list(p), sigma=e, absolute_sigma=True, maxfev=20000)
            chv = np.sum(((y - f(x, *pv)) / e) ** 2) / max(len(x) - len(pv), 1)
            res['v_raw'] = float(pv[0])
            res['v_raw_err'] = float(np.sqrt(cv[0, 0]) * np.sqrt(max(chv, 1)))
            model = f(x, *pv)
        except RuntimeError:
            pass
    return res, (x, y, e, model)


def ratio(a, ea, b, eb):
    r = a / b if b > 0 else np.nan
    er = abs(r) * np.hypot(ea / a, eb / b) if (a != 0 and b > 0) else np.nan
    lim2 = a / max(b + 2 * eb, 1e-6)        # lower limit on a/b using b's 2-sigma upper bound
    return dict(value=float(r), err=float(er), lower_2sigma=float(lim2))


def main():
    d, G = load()
    lam = 10 ** d['loglam']
    keys_arr = np.array([f'{p}-{m}-{s}' for p, m, s in zip(d['plate'], d['mjd'], d['spectrograph'])])
    gkey = lambda r: f"{int(r['plate'])}-{int(r['mjd'])}-{int(r['sp'])}"
    bad = {gkey(r) for r in G if (r['sb5040_ABmag'] < 20.4) or (r['I_OIII_cal_err'] > 0.6)}
    fpsf = {gkey(r): float(r['f_psf']) for r in G}
    good_groups = [gkey(r) for r in G if gkey(r) not in bad]
    print(len(good_groups), 'good groups of', len(G))

    zp = zeropoints(d, keys_arr, good_groups, lam)
    Z = Table(rows=[[r['group'], r['helio_rv'], r['n']] + [r[k] for k in ('OI5577', 'OI5577_err', 'OI6300', 'OI6300_err',
                                                                         'OI6364', 'OI6364_err')] for r in zp],
              names=['group', 'helio_rv', 'n_sky', 'dv_OI5577', 'dv_OI5577_err', 'dv_OI6300', 'dv_OI6300_err',
                     'dv_OI6364', 'dv_OI6364_err'])
    Z.write(ROOT / 'data/derived/spectra_sdss_lamost_sdss_zeropoints.csv', overwrite=True)
    zb = {r['group']: r['OI5577'] for r in zp}
    zr = {r['group']: np.nanmean([r['OI6300'], r['OI6364']]) for r in zp}
    print('zero points: blue median %.2f (16-84%%: %.2f..%.2f) km/s; red median %.2f (%.2f..%.2f)' % (
        np.nanmedian(list(zb.values())), *np.nanpercentile(list(zb.values()), [16, 84]),
        np.nanmedian(list(zr.values())), *np.nanpercentile(list(zr.values()), [16, 84])))

    geo = TR.load()
    dist, inside, useW = TR.nearest_strand_distance(d['ra'], d['dec'], geo)
    dE, inE = TR.signed_distance(d['ra'], d['dec'], geo['E'])
    dW, inW = TR.signed_distance(d['ra'], d['dec'], geo['W'])
    adist = np.abs(dist)
    dec = d['dec']
    on_all = inside & (adist < ON)
    off_all = inside & (adist > OFF[0]) & (adist < OFF[1])
    regions = {
        'all': (on_all, off_all),
        'south_dec_lt_20': (on_all & (dec < 20), off_all & (dec < 20)),
        'north_dec_ge_20': (on_all & (dec >= 20), off_all & (dec >= 20)),
        'E_strand_dec15_20': (inE & (np.abs(dE) < ON) & (np.abs(dW) > 0.3) & (dec >= 15) & (dec < 20),
                              off_all & (dec >= 15) & (dec < 20)),
        'W_strand_dec15_20': (inW & (np.abs(dW) < ON) & (np.abs(dE) > 0.3) & (dec >= 15) & (dec < 20),
                              off_all & (dec >= 15) & (dec < 20)),
        'null_1.0-1.15deg': (inside & (adist > 1.0) & (adist < 1.15), inside & (adist > 1.3) & (adist < 2.5)),
    }
    wd = np.nanmedian(np.where(d['wdisp'] > 0, d['wdisp'], np.nan), 0)
    out = {'on_deg': ON, 'off_deg': OFF, 'trace_version': geo['version'], 'airglow_vac_A': {k: v[0] for k, v in AIRGLOW.items()},
           'zeropoint_summary': dict(blue_median=float(np.nanmedian(list(zb.values()))),
                                     blue_p16_p84=[float(x) for x in np.nanpercentile(list(zb.values()), [16, 84])],
                                     red_median=float(np.nanmedian(list(zr.values()))),
                                     red_p16_p84=[float(x) for x in np.nanpercentile(list(zr.values()), [16, 84])]),
           'units': 'Rayleigh (diffuse-corrected: x f_psf); velocities heliocentric km/s', 'regions': {}}
    curves = {}
    for reg, (onm, offm) in regions.items():
        S = stack(d, lam, good_groups, keys_arr, onm, offm, fpsf)
        if S is None:
            continue
        w_on = {k: n for k, n in S['used']}
        tot_on = sum(w_on.values())
        zpb = float(np.nansum([zb.get(k, np.nan) * n for k, n in w_on.items() if np.isfinite(zb.get(k, np.nan))]) /
                    max(sum(n for k, n in w_on.items() if np.isfinite(zb.get(k, np.nan))), 1))
        zpr = float(np.nansum([zr.get(k, np.nan) * n for k, n in w_on.items() if np.isfinite(zr.get(k, np.nan))]) /
                    max(sum(n for k, n in w_on.items() if np.isfinite(zr.get(k, np.nan))), 1))
        # scatter of zero points among the contributing groups -> systematic error of the mean
        zbl = np.array([zb.get(k, np.nan) for k in w_on]); zrl = np.array([zr.get(k, np.nan) for k in w_on])
        r = dict(n_on=S['n_on'], n_groups=S['n_groups'], err_inflation=S['err_inflation'], zp_blue=zpb, zp_red=zpr,
                 zp_blue_scatter=float(np.nanstd(zbl)), zp_red_scatter=float(np.nanstd(zrl)),
                 dec_range_on=[float(np.min(dec[onm])), float(np.max(dec[onm]))] if onm.sum() else None, lines={})
        for nm, lines, lo, hi, cam in SETS:
            sig_pix = np.nanmedian(wd[(lam > lo) & (lam < hi)]) * 1e-4 * np.log(10)
            res, cur = fit_set(lam, S['stack'], S['err'], lines, lo, hi, sig_pix, free_v=nm in ('OIII', 'HaNII'))
            if 'v_raw' in res:
                zp_ = zpb if cam == 'blue' else zpr
                res['v_helio'] = res['v_raw'] - zp_
                res['v_helio_err_stat'] = res['v_raw_err']
                zsc = r['zp_blue_scatter'] if cam == 'blue' else r['zp_red_scatter']
                res['v_helio_err_zp'] = float(zsc / np.sqrt(max(S['n_groups'], 1)))
            r['lines'][nm] = res
            curves[(reg, nm)] = cur
        L = {}
        for nm in r['lines']:
            for k, v in r['lines'][nm].items():
                if isinstance(v, dict):
                    L[k] = v
        r['ratios'] = {
            'OIII5007/Hb': ratio(L['OIII5007']['I'], L['OIII5007']['err'], L['Hb']['I'], L['Hb']['err']),
            'OIII5007/Ha': ratio(L['OIII5007']['I'], L['OIII5007']['err'], L['Ha']['I'], L['Ha']['err']),
            'OIII5007/4959': ratio(L['OIII5007']['I'], L['OIII5007']['err'], L['OIII4959']['I'], L['OIII4959']['err']),
            'NII6583/Ha': ratio(L['NII6583']['I'], L['NII6583']['err'], L['Ha']['I'], L['Ha']['err']),
            'SII/Ha': ratio(L['SII6716']['I'] + L['SII6731']['I'], np.hypot(L['SII6716']['err'], L['SII6731']['err']),
                            L['Ha']['I'], L['Ha']['err']),
            'Ha/Hb': ratio(L['Ha']['I'], L['Ha']['err'], L['Hb']['I'], L['Hb']['err']),
        }
        r['Hb_upper_2sigma'] = float(max(L['Hb']['I'], 0) + 2 * L['Hb']['err'])
        out['regions'][reg] = r
        o = r['lines']['OIII']; h = r['lines']['HaNII']
        print(f"{reg:20s} nOn {S['n_on']:4d} ng {S['n_groups']:3d} | OIII {L['OIII5007']['I']:.2f}+-{L['OIII5007']['err']:.2f} "
              f"Hb {L['Hb']['I']:.3f}+-{L['Hb']['err']:.3f} Ha {L['Ha']['I']:.2f}+-{L['Ha']['err']:.2f} "
              f"NII {L['NII6583']['I']:.2f}+-{L['NII6583']['err']:.2f} SII {L['SII6716']['I'] + L['SII6731']['I']:.2f} | "
              f"v[OIII] {o.get('v_helio', np.nan):+.1f}+-{o.get('v_raw_err', np.nan):.1f} (zp {zpb:+.1f}) "
              f"v[NII/Ha] {h.get('v_helio', np.nan):+.1f}+-{h.get('v_raw_err', np.nan):.1f} (zp {zpr:+.1f})")
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_substacks.json').write_text(json.dumps(out, indent=1, default=float))

    # ------------------------------------------------------------------ figure: south vs north vs null
    fig, axes = new_sky_figure(12.5, 8.6, 'SDSS sky fibres · Silk strand spectrum, south vs north',
                               f'ON (<{ON}° from strand) minus same-plate OFF ({OFF[0]}–{OFF[1]}°), diffuse-corrected; '
                               'heliocentric vacuum λ', nrows=2, ncols=2)
    styles = [('south_dec_lt_20', 'Dec < 20°', None, 1.4), ('north_dec_ge_20', 'Dec ≥ 20°', '#b0bec5', 1.0),
              ('null_1.0-1.15deg', 'null (1.0–1.15°)', '#5f6b76', 0.8)]
    ttl = {'OIII': ('[O III] 4959, 5007', OIII_C), 'Hb': ('Hβ', HB_C), 'HaNII': ('Hα + [N II]', NII_C),
           'SII': ('[S II] 6716, 6731', SII_C)}
    for ax, (nm, lines, lo, hi, cam) in zip(axes.ravel(), SETS):
        name, col = ttl[nm]
        for reg, lab, c, lw in styles:
            cur = curves.get((reg, nm))
            if cur is None:
                continue
            x, y, e, mod = cur
            cc = col if c is None else c
            if reg.startswith('south'):
                ax.fill_between(x, y - e, y + e, color=cc, alpha=0.2, lw=0)
            ax.plot(x, y, color=cc, lw=lw, drawstyle='steps-mid', label=lab)
        for ln in lines:
            ax.axvline(VAC[ln] * (1 + 20 / C_KMS), color=FGDIM, lw=0.5, ls=':')
        ax.axhline(0, color=FGDIM, lw=0.5)
        label(ax, 0.02, 0.9, name, color=col, size=9.5, transform=ax.transAxes)
        rs = out['regions']['south_dec_lt_20']['lines'][nm]
        rn = out['regions']['north_dec_ge_20']['lines'][nm]
        key = lines[-1]
        txt = f"{key.replace('OIII', '[O III] ').replace('NII', '[N II] ').replace('SII', '[S II] ')}: " \
              f"S {rs[key]['I']:.2f}±{rs[key]['err']:.2f}, N {rn[key]['I']:.2f}±{rn[key]['err']:.2f} R"
        if nm == 'HaNII':
            txt = f"Hα: S {rs['Ha']['I']:.2f}±{rs['Ha']['err']:.2f}, N {rn['Ha']['I']:.2f}±{rn['Ha']['err']:.2f} R"
        label(ax, 0.02, 0.05, txt, color=WARN, size=7.5, transform=ax.transAxes)
        ax.set_ylabel('R / Å')
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7, loc='upper right')
    for ax in axes[-1]:
        ax.set_xlabel('vacuum wavelength, heliocentric (Å)')
    finish(fig, ROOT / 'figures/spectra_sdss_lamost_sdss_substacks.png')


if __name__ == '__main__':
    main()
