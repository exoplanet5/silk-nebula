"""Differential stack of SDSS/BOSS sky-fibre spectra: fibres ON an [OIII] strand minus fibres OFF it in the
same plate-spectrograph (topic spectra_sdss_lamost).  Airglow, zodiacal light and the pipeline sky model
cancel; what remains is the strand's own spectrum (diffuse-corrected, Rayleigh per Angstrom).

Inputs : data/derived/spectra_sdss_lamost_sdss_skyfibres.npz, ..._sdss_fibres.csv, ..._sdss_groups.csv
Outputs: data/derived/spectra_sdss_lamost_sdss_diffstack.json, figures/spectra_sdss_lamost_sdss_diffstack.png
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS, rayleigh_per_fibre_flux, FIBRE_DIAM
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, NII_C, SII_C, FGDIM, WARN
import spectra_sdss_lamost_trace as TR

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
ON, OFF = 0.15, (0.45, 2.0)
PANELS = [('[O III] 4959/5007', 4940, 5030, ['OIII4959', 'OIII5007'], OIII_C),
          ('Hβ', 4835, 4890, ['Hb'], HB_C),
          ('Hα + [N II]', 6530, 6600, ['NII6548', 'Ha', 'NII6583'], NII_C),
          ('[S II] 6716/6731', 6700, 6750, ['SII6716', 'SII6731'], SII_C)]


def main():
    with np.load(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz') as z:
        d = {k: z[k] for k in z.files}          # load once (NpzFile re-reads arrays on every access)
    lam = 10 ** d['loglam']
    G = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_groups.csv')
    bad = {(int(r['plate']), int(r['mjd']), int(r['sp'])) for r in G
           if (r['sb5040_ABmag'] < 20.4) or (r['I_OIII_cal_err'] > 0.6)}
    fpsf = {(int(r['plate']), int(r['mjd']), int(r['sp'])): float(r['f_psf']) for r in G}
    geo = TR.load()
    dist, inside, _ = TR.nearest_strand_distance(d['ra'], d['dec'], geo)
    keys = list(zip(d['plate'].tolist(), d['mjd'].tolist(), d['spectrograph'].tolist()))
    diffs, wts, ngrp, non = [], [], 0, 0
    for k in sorted(set(keys)):
        if k in bad or k not in fpsf:
            continue
        m = np.array([kk == k for kk in keys]) & inside
        on = m & (np.abs(dist) < ON)
        off = m & (np.abs(dist) > OFF[0]) & (np.abs(dist) < OFF[1])
        if on.sum() < 1 or off.sum() < 4:
            continue
        diam = FIBRE_DIAM[str(d['instrument'][np.where(m)[0][0]])]
        conv = 1e-17 / rayleigh_per_fibre_flux(lam, diam) * fpsf[k]          # flux per fibre -> R/A, diffuse-corrected
        def wmean(idx):
            f = d['flux'][idx]
            iv = d['ivar'][idx]
            w = np.where(iv > 0, iv, 0)
            return np.sum(f * w, 0) / np.maximum(np.sum(w, 0), 1e-30), np.sum(w, 0)
        fon, won = wmean(on)
        foff, woff = wmean(off)
        diff = (fon - foff) * conv
        var = (1 / np.maximum(won, 1e-30) + 1 / np.maximum(woff, 1e-30)) * conv ** 2
        good = (won > 0) & (woff > 0)
        diffs.append(np.where(good, diff, np.nan))
        wts.append(np.where(good, 1 / var, 0))
        ngrp += 1
        non += int(on.sum())
    D = np.array(diffs)
    W = np.array(wts)
    # robust: 5-sigma clip per pixel against the median
    med = np.nanmedian(D, 0)
    sd = 1.4826 * np.nanmedian(np.abs(D - med), 0)
    W = np.where(np.isfinite(D) & (np.abs(D - med) < 5 * np.maximum(sd, 1e-6)), W, 0)
    stack = np.nansum(np.nan_to_num(D) * W, 0) / np.maximum(W.sum(0), 1e-30)
    err = 1 / np.sqrt(np.maximum(W.sum(0), 1e-30))
    # inflate errors by the scatter-based reduced chi2 in line-free regions
    lf = (lam > 4890) & (lam < 4940) | (lam > 5030) & (lam < 5065)
    rchi = np.nanmean(((stack[lf] - np.nanmedian(stack[lf])) / err[lf]) ** 2)
    err *= np.sqrt(max(rchi, 1))
    res = {'n_groups': ngrp, 'n_on_fibres': non, 'on_deg': ON, 'off_deg': OFF, 'err_inflation': float(np.sqrt(max(rchi, 1))),
           'lines': {}}
    wd = np.nanmedian(np.where(d['wdisp'] > 0, d['wdisp'], np.nan), 0)
    fig, axes = new_sky_figure(12, 7.6, 'SDSS sky fibres · ON-strand minus OFF-strand spectrum',
                               f'{non} fibres within {ON}° of an [O III] strand, {ngrp} plate-spectrographs; '
                               'same-plate OFF fibres 0.45–2° away; diffuse-corrected', nrows=2, ncols=2)
    for ax, (ttl, lo, hi, lines, col) in zip(axes.ravel(), PANELS):
        m = (lam > lo) & (lam < hi) & (err > 0) & np.isfinite(stack)
        x, y, e = lam[m], stack[m], err[m]
        sig = np.nanmedian(wd[m]) * 1e-4 * np.log(10) * x.mean()
        sig = np.hypot(sig, x.mean() * 19 / C_KMS)

        def f(xx, b, *a):
            out = b + 0 * xx
            for ln, aa in zip(lines, a):
                c = VAC[ln] * (1 + 20 / C_KMS)
                out = out + aa * np.exp(-0.5 * ((xx - c) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
            return out
        p, c = curve_fit(f, x, y, p0=[0] + [0.5] * len(lines), sigma=e, absolute_sigma=True)
        chi = np.sum(((y - f(x, *p)) / e) ** 2) / max(len(x) - len(p), 1)
        sc = np.sqrt(max(chi, 1))
        for i, ln in enumerate(lines):
            res['lines'][ln] = dict(I=float(p[i + 1]), err=float(np.sqrt(c[i + 1, i + 1]) * sc))
        # free common velocity (heliocentric) for line sets with a significant detection

        def fv(xx, v, b, *a):
            out = b + 0 * xx
            for ln, aa in zip(lines, a):
                cc = VAC[ln] * (1 + v / C_KMS)
                out = out + aa * np.exp(-0.5 * ((xx - cc) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
            return out
        try:
            pv, cv = curve_fit(fv, x, y, p0=[20.0] + list(p), sigma=e, absolute_sigma=True, maxfev=20000)
            chv = np.sum(((y - fv(x, *pv)) / e) ** 2) / max(len(x) - len(pv), 1)
            res['lines'][lines[-1]]['v_helio_free'] = float(pv[0])
            res['lines'][lines[-1]]['v_helio_free_err'] = float(np.sqrt(cv[0, 0]) * np.sqrt(max(chv, 1)))
        except RuntimeError:
            pass
        ax.fill_between(x, y - e, y + e, color=col, alpha=0.2, lw=0)
        ax.plot(x, y, color=col, lw=1.2, drawstyle='steps-mid')
        ax.plot(x, f(x, *p), color=WARN, lw=1.0, ls='--')
        for ln in lines:
            ax.axvline(VAC[ln] * (1 + 20 / C_KMS), color=FGDIM, lw=0.5, ls=':')
        ax.axhline(0, color=FGDIM, lw=0.5)
        txt = ', '.join(f"{ln.replace('OIII', '').replace('NII', 'N').replace('SII', 'S')}: {res['lines'][ln]['I']:.2f}±{res['lines'][ln]['err']:.2f}"
                        for ln in lines)
        label(ax, 0.02, 0.88, ttl, color=col, size=9.5, transform=ax.transAxes)
        label(ax, 0.02, 0.05, txt + ' R', color=WARN, size=7.5, transform=ax.transAxes)
        ax.set_ylabel('R / Å')
        ax.grid(alpha=0.2)
    for ax in axes[-1]:
        ax.set_xlabel('vacuum wavelength, heliocentric (Å)')
    L = res['lines']
    res['ratios'] = {
        'OIII5007/Hb': L['OIII5007']['I'] / L['Hb']['I'] if L['Hb']['I'] > 0 else None,
        'OIII5007/Hb_lower_2sigma': L['OIII5007']['I'] / max(L['Hb']['I'] + 2 * L['Hb']['err'], 1e-3),
        'OIII5007/Ha': L['OIII5007']['I'] / L['Ha']['I'] if L['Ha']['I'] > 0 else None,
        'OIII5007/OIII4959': L['OIII5007']['I'] / L['OIII4959']['I'] if L['OIII4959']['I'] > 0 else None,
        'NII6583/Ha': L['NII6583']['I'] / L['Ha']['I'] if L['Ha']['I'] > 0 else None,
        'SII/Ha': (L['SII6716']['I'] + L['SII6731']['I']) / L['Ha']['I'] if L['Ha']['I'] > 0 else None,
        'SII6716/6731': L['SII6716']['I'] / L['SII6731']['I'] if L['SII6731']['I'] > 0 else None,
    }
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_diffstack.json').write_text(json.dumps(res, indent=1))
    np.savez(ROOT / 'data/derived/spectra_sdss_lamost_sdss_diffstack.npz', lam=lam, stack=stack, err=err)
    print(json.dumps(res, indent=1))
    finish(fig, ROOT / 'figures/spectra_sdss_lamost_sdss_diffstack.png')


if __name__ == '__main__':
    main()
