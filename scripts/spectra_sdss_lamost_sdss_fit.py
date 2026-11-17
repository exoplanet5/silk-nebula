"""Absolute nebular-line intensities from SDSS DR17 sky-fibre spectra (topic tag spectra_sdss_lamost).

Method
------
(1) For every plate-mjd-spectrograph group, the mean total sky spectrum of its SKY fibres
    (flux + pipeline sky = what each sky fibre saw) is fitted, line by line, with
        sky(l) = S_conv(l; v) * (a0 + a1 x + a2 x^2) [+ broad lamp bump] [+ OH lines] + sum_k F_k G_k(l)
    where S_conv is the Kurucz solar spectrum (zodiacal light / scattered sunlight / integrated starlight
    all have near-solar absorption spectra) broadened to the SDSS LSF (wdisp), G_k unit-area Gaussians of
    the instrumental width at the vacuum line wavelengths for V_helio = +20 km/s, and F_k the line fluxes
    (1e-17 erg/s/cm^2 per fibre).  The fit is linear in all amplitudes; the solar-template velocity is
    chosen on a grid (-60..+60 km/s).  [OIII]4959 is tied to 5007/2.98.
(2) Per sky fibre, the pipeline residual (flux column) is fitted with a local linear baseline plus the
    same Gaussians -> Delta F (fibre minus spectrograph mean).  Absolute per-fibre intensity
    = group-mean F + Delta F - <Delta F>_group.
(3) Fluxes per fibre -> Rayleigh using the fibre solid angle (3" SDSS-I/II, 2" BOSS).  The SDSS flux
    calibration is tuned to point sources (PSF photometry), so diffuse light is over-estimated by
    1/f_psf, f_psf = point-source fraction inside the fibre.  We report the as-calibrated value I_cal
    and I = f_psf * I_cal with f_psf from a Gaussian PSF of the plate's median seeing (systematic ~20%).

Output: data/derived/spectra_sdss_lamost_sdss_groups.csv, ..._sdss_fibres.csv, ..._sdss_calib.json

TOPO note (fixed 2026-11-16 18:10): the spec-lite header keyword HELIO_RV has the opposite sign of the
astropy heliocentric correction v_corr (checked on 12 plates, |difference| < 1.2 km/s), and the [O I] 5577
centroid follows lambda0 * (1 - HELIO_RV/c) (slope -2.03 vs the earlier assumption).  Topocentric features
(airglow, OH, scattered sunlight) are therefore placed at lambda0 * (1 - HELIO_RV/c).  The nebular lines,
fibre residuals and on/off differences never depended on this; only the group-level airglow/OH modelling did.
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.optimize import least_squares

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS, sun_at, gauss, linfit, flux_to_rayleigh, FIBRE_DIAM

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
V_NEB = 20.0                      # km/s heliocentric, adopted line centre (Galactic; [OIII] +25, red +12..+20)
PIX_KMS = 1e-4 * np.log(10) * C_KMS   # 69.03 km/s per SDSS pixel
SIG_INTR = 19.0                   # km/s intrinsic sigma ([OIII] FWHM ~45 km/s, paper Table 2)
VGRID = np.arange(-60, 61, 10.0)

WIN = {  # fitting windows (vacuum A); 'mask' = lamp/telluric pieces to exclude
    'OIII': dict(lo=4935.0, hi=5045.0, lines=['OIII'], bump=True),
    'Hb': dict(lo=4836.0, hi=4895.0, lines=['Hb'], bump=False),
    'RED': dict(lo=6520.0, hi=6760.0, lines=['Ha', 'NII', 'SII6716', 'SII6731'], bump=False),
}
NULL = {'OIII': 5026.0, 'Hb': 4880.0}      # fake-line positions for null tests (no known features)


def line_profile(name, lam, sig_pix_kms, v=V_NEB):
    """Unit-flux profile (per A) of a (possibly tied) nebular line set."""
    sig = np.sqrt(sig_pix_kms ** 2 + SIG_INTR ** 2)
    def g(l0):
        c = l0 * (1 + v / C_KMS)
        return gauss(lam, c, c * sig / C_KMS)
    if name == 'OIII':
        return g(VAC['OIII5007']) + g(VAC['OIII4959']) / 2.98
    if name == 'NII':
        return g(VAC['NII6583']) + g(VAC['NII6548']) / 2.94
    if name == 'Hb':
        return g(VAC['Hb'])
    if name == 'Ha':
        return g(VAC['Ha'])
    if name == 'SII6716':
        return g(VAC['SII6716'])
    if name == 'SII6731':
        return g(VAC['SII6731'])
    if name.startswith('null'):
        return g(float(name[4:]))
    raise KeyError(name)


# ------------------------------------------------------------------------------------ groups
def group_means(d, lam):
    tot = d['flux'] + d['sky']
    keys = sorted(set(zip(d['plate'].tolist(), d['mjd'].tolist(), d['spectrograph'].tolist())))
    out = []
    for k in keys:
        m = (d['plate'] == k[0]) & (d['mjd'] == k[1]) & (d['spectrograph'] == k[2])
        if m.sum() < 6:
            continue
        t = np.where(d['ivar'][m] > 0, tot[m], np.nan)
        # reject fibres whose residual continuum is deviant (object light in a sky fibre)
        rc = np.nanmedian(np.where(d['ivar'][m] > 0, d['flux'][m], np.nan)[:, (lam > 5030) & (lam < 5060)], axis=1)
        mad = 1.4826 * np.nanmedian(np.abs(rc - np.nanmedian(rc)))
        good = np.abs(rc - np.nanmedian(rc)) < max(4 * mad, 0.3)
        t = t[good]
        med = np.nanmedian(t, 0)
        sd = 1.4826 * np.nanmedian(np.abs(t - med), 0)
        clip = np.abs(t - med) < 4 * np.maximum(sd, 1e-3)
        tc = np.where(clip, t, np.nan)
        mean = np.nanmean(tc, 0)
        nn = np.sum(np.isfinite(tc), 0)
        err = np.nanstd(tc, 0) / np.sqrt(np.maximum(nn, 1))
        wd = np.nanmedian(np.where(d['wdisp'][m] > 0, d['wdisp'][m], np.nan), 0)
        sel = np.where(m)[0]
        out.append(dict(plate=k[0], mjd=k[1], sp=k[2], n=int(good.sum()), mean=mean, err=err, wdisp=wd,
                        ra=float(np.mean(d['ra'][m])), dec=float(np.mean(d['dec'][m])),
                        instrument=str(d['instrument'][sel[0]]), seeing=float(d['seeing50'][sel[0]]),
                        helio=float(np.nanmedian(d['helio_rv'][m])), airmass=float(d['airmass'][sel[0]]),
                        idx=sel[good]))
    return out


# ------------------------------------------------------------------------------- calibrations
def calibrate(groups, lam):
    """Global shape parameters: lamp bump near 4980 A and OH line centres in the red window."""
    # stack in the topocentric frame, normalised
    stk = []
    for g in groups:
        lt = lam / (1 - g['helio'] / C_KMS)       # helio-frame -> topocentric (SDSS HELIO_RV = -v_corr; see TOPO note)
        y = np.interp(lam, lt, g['mean'])
        stk.append(y / np.nanmedian(y[(lam > 6610) & (lam < 6690)]))
    stk = np.nanmedian(np.array(stk), 0)
    wd = np.nanmedian([g['wdisp'] for g in groups], 0)
    sig_kms = np.nanmedian(wd[(lam > 6500)]) * PIX_KMS
    # --- red window: find OH peaks
    m = (lam > 6500) & (lam < 6780)
    x, y = lam[m], stk[m]
    S = sun_at(x, sig_kms)
    cont = np.polyval(np.polyfit(x, y / S, 2, w=None), x) * S
    r = y - cont
    from scipy.signal import find_peaks
    pk, _ = find_peaks(r, height=0.04 * np.nanmedian(y), distance=2)
    cands = x[pk]
    neb = [VAC['Ha'], VAC['NII6583'], VAC['NII6548'], VAC['SII6716'], VAC['SII6731']]
    oh0 = [c for c in cands if min(abs(c - n) for n in neb) > 2.5]
    sigA = np.median(wd[m]) * 1e-4 * np.log(10) * 6600

    def model(p):
        n = len(oh0)
        cen = p[:n]; amp = p[n:2 * n]; c = p[2 * n:2 * n + 3]; neba = p[2 * n + 3:]
        mod = S * np.polyval(c, x - 6640)
        for ci, ai in zip(cen, amp):
            mod += ai * gauss(x, ci, sigA)
        for ni, ai in zip(neb, neba):
            mod += ai * gauss(x, ni, sigA)
        return mod
    p0 = np.concatenate([oh0, np.full(len(oh0), 0.5), [0, 0, np.nanmedian(y / S)], np.full(len(neb), 0.2)])
    res = least_squares(lambda p: model(p) - y, p0)
    oh_cen = res.x[:len(oh0)]
    oh_amp = res.x[len(oh0):2 * len(oh0)]
    keep = oh_amp > 0.05
    # --- lamp bump around 4980 A
    m2 = (lam > 4935) & (lam < 5045)
    x2, y2 = lam[m2], stk[m2] / np.nanmedian(stk[m2])
    S2 = sun_at(x2, np.nanmedian(wd[m2]) * PIX_KMS)
    s5 = np.median(wd[m2]) * 1e-4 * np.log(10) * 5000

    def model2(p):
        return (S2 * np.polyval(p[:3], x2 - 4990) + p[3] * gauss(x2, p[4], p[5])
                + p[6] * (gauss(x2, VAC['OIII5007'], s5) + gauss(x2, VAC['OIII4959'], s5) / 2.98))
    r2 = least_squares(lambda p: model2(p) - y2, [0, 0, 1, 10, 4982, 8, 0.1],
                       bounds=([-np.inf] * 3 + [0, 4970, 3, -5], [np.inf] * 3 + [100, 4995, 20, 5]))
    cal = dict(oh_centres_topo=[float(c) for c in oh_cen[keep]], oh_rel_amp=[float(a) for a in oh_amp[keep]],
               bump_centre=float(r2.x[4]), bump_sigma=float(r2.x[5]), bump_amp_rel=float(r2.x[3]),
               red_nebular_rel_amp=[float(a) for a in res.x[2 * len(oh0) + 3:]])
    return cal, stk


# --------------------------------------------------------------------------------- group fits
def fit_group(g, lam, cal, win, extra_null=None):
    w = WIN[win]
    m = (lam > w['lo']) & (lam < w['hi']) & np.isfinite(g['mean']) & (g['err'] > 0)
    x, y, e = lam[m], g['mean'][m], g['err'][m]
    # add a 1% floor to the per-pixel error (template mismatch)
    e = np.sqrt(e ** 2 + (0.01 * np.abs(y)) ** 2)
    iv = 1 / e ** 2
    sig_kms = float(np.nanmedian(g['wdisp'][m]) * PIX_KMS)
    lines = list(w['lines']) + ([extra_null] if extra_null else [])
    best = None
    for v in VGRID:
        S = sun_at(x, sig_kms, v_kms=v - g['helio'])   # zodiacal/solar features ~ at rest in the topocentric frame
        xx = (x - x.mean()) / (x.max() - x.min())
        cols = [S, S * xx, S * xx ** 2]
        names = ['c0', 'c1', 'c2']
        if w['bump']:
            cols.append(gauss(x, cal['bump_centre'], cal['bump_sigma']) * (1 + 0 * x)); names.append('bump')
            cols.append(gauss(x, cal['bump_centre'], cal['bump_sigma']) * xx); names.append('bump1')
        if win == 'RED':
            sA = sig_kms / C_KMS
            for i, c in enumerate(cal['oh_centres_topo']):
                ch = c * (1 - g['helio'] / C_KMS)
                cols.append(gauss(x, ch, ch * sA)); names.append(f'oh{i}')
        for ln in lines:
            cols.append(line_profile(ln, x, sig_kms)); names.append(ln)
        coef, cov, chi2, model = linfit(y, iv, cols)
        if best is None or chi2 < best[2]:
            best = (coef, cov, chi2, model, v, names, x, y, e)
    coef, cov, chi2, model, v, names, x, y, e = best
    out = {'v_sun': v, 'chi2': chi2, 'dof': len(x) - len(coef)}
    rchi = max(chi2 / out['dof'], 1.0)
    for ln in lines:
        i = names.index(ln)
        out[ln] = coef[i]
        out[ln + '_err'] = np.sqrt(cov[i, i] * rchi)
    if win == 'RED':
        oh = [coef[names.index(f'oh{i}')] for i in range(len(cal['oh_centres_topo']))]
        out['oh_sum'] = float(np.sum(oh))
        j = int(np.argmin(np.abs(np.array(cal['oh_centres_topo']) - 6555.6)))
        out['oh6555'] = float(oh[j])
    out['_model'] = (x, y, e, model)
    return out


# ---------------------------------------------------------------------------- per-fibre fits
FWIN = {'OIII': [(4944.0, 4976.0), (4994.0, 5024.0)], 'Hb': [(4848.0, 4878.0)],
        'Ha': [(6545.0, 6605.0)], 'SII': [(6705.0, 6745.0)]}


def fit_fibre(lam, flux, ivar, wdisp, helio, cal):
    res = {}
    for key, wins in FWIN.items():
        m = np.zeros_like(lam, bool)
        for a, b in wins:
            m |= (lam > a) & (lam < b)
        m &= ivar > 0
        if m.sum() < 15:
            continue
        x, y, iv = lam[m], flux[m], ivar[m]
        sig_kms = float(np.nanmedian(wdisp[m]) * PIX_KMS)
        xx = (x - x.mean()) / 30.0
        cols = [np.ones_like(x), xx]
        if key == 'OIII':
            # separate baselines per sub-window
            sub = x > 4990
            cols = [sub * 1.0, (~sub) * 1.0, xx]
            names = ['b0', 'b1', 'b2', 'OIII']
            cols.append(line_profile('OIII', x, sig_kms))
        elif key == 'Hb':
            names = ['b0', 'b1', 'Hb']
            cols.append(line_profile('Hb', x, sig_kms))
        elif key == 'Ha':
            names = ['b0', 'b1', 'Ha', 'NII']
            cols += [line_profile('Ha', x, sig_kms), line_profile('NII', x, sig_kms)]
        else:
            names = ['b0', 'b1', 'SII6716', 'SII6731']
            cols += [line_profile('SII6716', x, sig_kms), line_profile('SII6731', x, sig_kms)]
        coef, cov, chi2, model = linfit(y, iv, cols)
        rchi = max(chi2 / max(len(x) - len(coef), 1), 1.0)
        for nm in names[len(names) - (2 if key in ('Ha', 'SII') else 1):]:
            i = names.index(nm)
            res['d' + nm] = coef[i]
            res['d' + nm + '_err'] = np.sqrt(cov[i, i] * rchi)
    return res


def psf_fraction(seeing, diam):
    s = seeing / 2.3548
    return 1 - np.exp(-(diam / 2) ** 2 / (2 * s ** 2))


def main():
    with np.load(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz') as z:
        d = {k: z[k] for k in z.files}          # load once (NpzFile re-reads arrays on every access)
    lam = 10 ** d['loglam']
    groups = group_means(d, lam)
    print(len(groups), 'groups')
    cal, stk = calibrate(groups, lam)
    print(json.dumps({k: (v if not isinstance(v, list) else [round(a, 2) for a in v]) for k, v in cal.items()}))
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_calib.json').write_text(json.dumps(cal, indent=1))
    np.save(ROOT / 'data/derived/spectra_sdss_lamost_sdss_globalstack.npy', np.vstack([lam, stk]))

    rows = []
    models = {}
    for g in groups:
        r = dict(plate=g['plate'], mjd=g['mjd'], sp=g['sp'], n_sky=g['n'], ra=g['ra'], dec=g['dec'],
                 instrument=g['instrument'], seeing=g['seeing'], airmass=g['airmass'], helio_rv=g['helio'])
        diam = FIBRE_DIAM[g['instrument']]
        r['f_psf'] = psf_fraction(g['seeing'] if np.isfinite(g['seeing']) and g['seeing'] > 0.5 else 1.5, diam)
        for win in ['OIII', 'Hb', 'RED']:
            fr = fit_group(g, lam, cal, win)
            models[(g['plate'], g['mjd'], g['sp'], win)] = fr.pop('_model')
            if win in NULL:
                fn = fit_group(g, lam, cal, win, extra_null=f'null{NULL[win]}')
                fr['null'] = fn[f'null{NULL[win]}']
            for k, v in fr.items():
                r[f'{win}_{k}' if k in ('v_sun', 'chi2', 'dof', 'null', 'oh_sum', 'oh6555') else k] = v
            if win == 'RED':
                r['I_OH6555_cal'] = float(flux_to_rayleigh(r['RED_oh6555'], 6555.6, diam))
        # continuum level near 5007 for reference (1e-17 erg/s/cm2/A per fibre)
        mm = (lam > 5030) & (lam < 5045)
        r['cont5040'] = float(np.nanmedian(g['mean'][mm]))
        # surface brightness of the continuum (AB mag/arcsec^2 at 5040 A), as calibrated
        fl_arc = r['cont5040'] * 1e-17 / (np.pi * (diam / 2) ** 2)
        r['sb5040_ABmag'] = float(-2.5 * np.log10(max(fl_arc, 1e-30) * 5040 ** 2 / 2.99792458e18 / 3.631e-20))
        # airglow [O I] 5577 (absolute-calibration sanity check; typical 100-400 R)
        m5 = (lam > 5562) & (lam < 5596) & np.isfinite(g['mean'])
        x5 = lam[m5]
        s5 = float(np.nanmedian(g['wdisp'][m5]) * PIX_KMS)
        c5 = 5578.89 * (1 - g['helio'] / C_KMS)
        co, cv, _, _ = linfit(g['mean'][m5], 1 / np.maximum(g['err'][m5], 1e-3) ** 2,
                              [np.ones_like(x5), x5 - 5578, gauss(x5, c5, c5 * s5 / C_KMS)])
        r['I_OI5577_cal'] = float(flux_to_rayleigh(co[2], 5578.89, diam))
        # airglow [N I] 5198/5200 doublet (vac 5199.35/5201.71, blended at R~1800; profile 0.6/0.4)
        mn = (lam > 5183) & (lam < 5218) & np.isfinite(g['mean']) & (g['err'] > 0)
        if mn.sum() > 20:
            xn = lam[mn]
            sn = float(np.nanmedian(g['wdisp'][mn]) * PIX_KMS)
            Sn = sun_at(xn, sn, v_kms=-g['helio'])
            xx = (xn - 5200) / 35
            hv = 1 - g['helio'] / C_KMS
            prof = 0.6 * gauss(xn, 5199.35 * hv, 5199.35 * sn / C_KMS) + 0.4 * gauss(xn, 5201.71 * hv, 5201.71 * sn / C_KMS)
            co, cv, _, _ = linfit(g['mean'][mn], 1 / g['err'][mn] ** 2, [Sn, Sn * xx, Sn * xx ** 2, prof])
            r['I_NI5199_cal'] = float(flux_to_rayleigh(co[3], 5200.3, diam))
        else:
            r['I_NI5199_cal'] = np.nan
        for ln, l0 in [('OIII', VAC['OIII5007']), ('Hb', VAC['Hb']), ('Ha', VAC['Ha']), ('NII', VAC['NII6583']),
                       ('SII6716', VAC['SII6716']), ('SII6731', VAC['SII6731'])]:
            r[f'I_{ln}_cal'] = float(flux_to_rayleigh(r[ln], l0, diam))
            r[f'I_{ln}_cal_err'] = float(flux_to_rayleigh(r[ln + '_err'], l0, diam))
        for win, l0 in [('OIII', VAC['OIII5007']), ('Hb', VAC['Hb'])]:
            r[f'I_{win}_null_cal'] = float(flux_to_rayleigh(r[f'{win}_null'], l0, diam))
        rows.append(r)
    T = Table(rows)
    T.write(ROOT / 'data/derived/spectra_sdss_lamost_sdss_groups.csv', overwrite=True)

    # per-fibre
    gl = {(r['plate'], r['mjd'], r['sp']): r for r in rows}
    frow = []
    for g in groups:
        r0 = gl[(g['plate'], g['mjd'], g['sp'])]
        tmp = []
        for i in g['idx']:
            fr = fit_fibre(lam, d['flux'][i], d['ivar'][i], d['wdisp'][i], g['helio'], cal)
            fr.update(plate=g['plate'], mjd=g['mjd'], sp=g['sp'], fiber=int(d['fiber'][i]),
                      ra=float(d['ra'][i]), dec=float(d['dec'][i]), instrument=g['instrument'])
            tmp.append(fr)
        for key, gl_key, l0 in [('OIII', 'OIII', VAC['OIII5007']), ('Hb', 'Hb', VAC['Hb']), ('Ha', 'Ha', VAC['Ha']),
                                ('NII', 'NII', VAC['NII6583']), ('SII6716', 'SII6716', VAC['SII6716']),
                                ('SII6731', 'SII6731', VAC['SII6731'])]:
            vals = np.array([t.get('d' + key, np.nan) for t in tmp])
            errs = np.array([t.get('d' + key + '_err', np.nan) for t in tmp])
            ok = np.isfinite(vals) & (errs > 0)
            mu = np.sum(vals[ok] / errs[ok] ** 2) / np.sum(1 / errs[ok] ** 2) if ok.sum() > 2 else 0.0
            diam = FIBRE_DIAM[g['instrument']]
            for t, v_, e_ in zip(tmp, vals, errs):
                t[f'I_{key}_cal'] = float(flux_to_rayleigh(r0[gl_key] + v_ - mu, l0, diam))
                t[f'I_{key}_cal_err'] = float(flux_to_rayleigh(e_, l0, diam))
                t[f'dI_{key}_cal'] = float(flux_to_rayleigh(v_ - mu, l0, diam))
        for t in tmp:
            t['f_psf'] = r0['f_psf']
        frow += tmp
    keys = ['plate', 'mjd', 'sp', 'fiber', 'ra', 'dec', 'instrument', 'f_psf'] + sorted(set().union(*[t.keys() for t in frow]) - {'plate', 'mjd', 'sp', 'fiber', 'ra', 'dec', 'instrument', 'f_psf'})
    F = Table(rows=[{k: t.get(k, np.nan) for k in keys} for t in frow], names=keys)
    F.write(ROOT / 'data/derived/spectra_sdss_lamost_sdss_fibres.csv', overwrite=True)
    print('wrote', len(T), 'groups,', len(F), 'fibres')


if __name__ == '__main__':
    main()
