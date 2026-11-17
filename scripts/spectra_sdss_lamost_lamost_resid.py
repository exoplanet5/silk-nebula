"""Nebular residual lines in LAMOST DR12 v2.0 stellar spectra along the Silk strand (topic spectra_sdss_lamost).

LAMOST subtracts a super-sky built from the sky fibres of each spectrograph (LRS: ~19-20 per spectrograph;
MRS-N plate: 134-205 'sky' fibres per spectrograph = the nebular fibres).  A star on the strand therefore
carries a residual  I(strand at star) - <I over that spectrograph's sky fibres>.

Method (per window): every spectrum is continuum-normalised and shifted to its stellar rest frame (catalogue RV);
for each target a data-driven template = median of the K nearest OFF-strand stars (chi^2 in the window with the
nebular-line pixels masked) is shifted back to the target frame; target = template*(a+b x) + nebular Gaussians.
Residual spectra (target/template - 1, heliocentric frame) are stacked for ON (|d|<0.2-0.25 deg) and OFF stars.
Equivalent widths are converted to Rayleigh with the star's photometric flux density (PS1 g/r for LRS,
Gaia G/BP/RP for MRS) and the fibre aperture fraction for a Gaussian PSF of the plate seeing (3.3" fibre).

usage: python -I spectra_sdss_lamost_lamost_resid.py lrs|mrs|mrsmn
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, C_KMS, gauss, linfit, rayleigh_per_fibre_flux

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
V0 = 20.0
K = 25
MASK_TEMPLATES = False

CFG = {
    'lrs': dict(sample=ROOT / 'data/raw/lamost_dr12/lrs_strand_sample.csv', dir=ROOT / 'data/raw/lamost_dr12/lrs/lrs_strand_sample',
                prefix='lrs', dlog=1e-4, R=1800.0,
                windows={'OIII': (4925.0, 5045.0, ['OIII4959', 'OIII5007']), 'Hb': (4820.0, 4905.0, ['Hb']),
                         'Ha': (6510.0, 6620.0, ['NII6548', 'Ha', 'NII6583']), 'SII': (6680.0, 6770.0, ['SII6716', 'SII6731'])}),
    'mrs': dict(sample=ROOT / 'data/raw/lamost_dr12/mrs_strand_sample.csv', dir=ROOT / 'data/raw/lamost_dr12/mrs/mrs_strand_sample',
                prefix='mrs', dlog=1e-5, R=7500.0,
                windows={'OIII': (4990.0, 5030.0, ['OIII5007']), 'Ha': (6545.0, 6600.0, ['NII6548', 'Ha', 'NII6583']),
                         'SII': (6705.0, 6745.0, ['SII6716', 'SII6731'])}),
    # the paper's MRS-N plate: its 83 public star spectra as targets ('on'), OFF stars of other plates as templates
    'mrsmn': dict(sample=ROOT / 'data/derived/spectra_sdss_lamost_mnplate_sample.csv',
                  dir=[ROOT / 'data/raw/lamost_dr12/mrs/mrs_strand_sample', ROOT / 'data/raw/lamost_dr12/mrs/MN074209N160100N01'],
                  prefix='mrsmn', dlog=1e-5, R=7500.0,
                  windows={'OIII': (4990.0, 5030.0, ['OIII5007'])}),
}


def read_spec(f, mode):
    h = fits.open(f, memmap=False)
    hd = h[0].header
    out = {'seeing': float(hd.get('SEEING', np.nan) or np.nan), 'helio_rv': float(hd.get('HELIO_RV', np.nan) or np.nan),
           'nskies': hd.get('NSKIES'), 'spid': hd.get('SPID'), 'planid': hd.get('PLANID')}
    if mode == 'lrs':
        d = h[1].data
        out['w'] = [np.array(d['WAVELENGTH'][0], float)]
        out['f'] = [np.array(d['FLUX'][0], float)]
        out['iv'] = [np.array(d['IVAR'][0], float) * (np.array(d['ANDMASK'][0]) == 0)]
    else:
        out['w'], out['f'], out['iv'] = [], [], []
        for ext in ('COADD_B', 'COADD_R'):
            d = h[ext].data
            if len(d) == 0 or len(d['WAVELENGTH'][0]) == 0:
                continue
            out['w'].append(np.array(d['WAVELENGTH'][0], float))
            out['f'].append(np.array(d['FLUX'][0], float))
            out['iv'].append(np.array(d['IVAR'][0], float) * (np.array(d['ANDMASK'][0]) == 0))
    return out


def window_arrays(spec, lo, hi, grid):
    for w, f, iv in zip(spec['w'], spec['f'], spec['iv']):
        if w.min() < lo and w.max() > hi:
            m = (w > lo - 20) & (w < hi + 20)
            fi = np.interp(grid, w[m], f[m])
            ii = np.interp(grid, w[m], iv[m])
            ii[ii < 0] = 0
            return fi, ii
    return None, None


def line_mask(grid, lines, rv, sig_A):
    """mask (True = keep) nebular-line pixels; grid is in the star rest frame."""
    keep = np.ones_like(grid, bool)
    for ln in lines:
        lc = VAC[ln] * (1 + V0 / C_KMS) / (1 + rv / C_KMS)
        keep &= np.abs(grid - lc) > 3.5 * sig_A + 0.5
    return keep


def phot_flambda(row, lam, mode):
    """Total stellar flux density (erg/s/cm2/A) at lam from catalogue photometry (approximate)."""
    if mode == 'lrs':
        g, r = float(row['mag_ps_g']), float(row['mag_ps_r'])
        if not (10 < g < 25 and 10 < r < 25):
            gg = float(row['gaia_g_mean_mag'])
            if not 5 < gg < 25:
                return np.nan
            g = r = gg
        lg, lr = 4810.0, 6170.0
        fg = 3631e-23 * 10 ** (-0.4 * g) * 2.99792458e18 / lg ** 2
        fr = 3631e-23 * 10 ** (-0.4 * r) * 2.99792458e18 / lr ** 2
        # power-law interpolation in lambda
        a = np.log(fr / fg) / np.log(lr / lg)
        return fg * (lam / lg) ** a
    G = float(row['gaia_g_mean_mag'])
    if not 5 < G < 25:
        return np.nan
    # Gaia G Vega zero point: f_lambda ~ 2.5e-9 erg/s/cm2/A at 6400 A (pivot 6218 A); F-K stars: f(5008)/f(6400) ~ 1.0-1.2
    return 2.5e-9 * 10 ** (-0.4 * G) * 1.1


def aperture_fraction(seeing, diam=3.3):
    s = (seeing if np.isfinite(seeing) and seeing > 0.5 else 3.5) / 2.3548
    return 1 - np.exp(-(diam / 2) ** 2 / (2 * s ** 2))


def main():
    mode = sys.argv[1]
    cfg = CFG[mode]
    S = Table.read(cfg['sample'], format='ascii.csv')
    dirs = cfg['dir'] if isinstance(cfg['dir'], list) else [cfg['dir']]
    files = {int(Path(f).name.split('-')[1].split('.')[0]): f for dd in dirs for f in glob.glob(str(dd / '*.fits.gz'))}
    S = S[[int(o) in files for o in S['obsid']]]
    if mode == 'lrs':
        rv = np.array(S['z'], float) * C_KMS
    else:
        rv = np.array(S['rv_br1'], float)
        alt = np.array(S['rv_lasp1'], float)
        rv = np.where(rv > -999, rv, np.where(alt > -999, alt, np.nan))
    ok = np.isfinite(rv) & (np.abs(rv) < 400)
    S, rv = S[ok], rv[ok]
    print(mode, len(S), 'spectra with RV;', 'on', int(np.sum(S['sample'] == 'on')), 'off', int(np.sum(S['sample'] == 'off')))
    specs = [read_spec(files[int(o)], mode) for o in S['obsid']]
    results = {k: [] for k in cfg['windows']}
    stacks = {}
    meta = []
    for i, (row, sp) in enumerate(zip(S, specs)):
        meta.append(dict(obsid=int(row['obsid']), sample=str(row['sample']), dist=float(row['dist']), ra=float(row['ra']),
                         dec=float(row['dec']), rv=float(rv[i]), seeing=sp['seeing'], planid=str(sp['planid']),
                         spid=sp['spid'], nskies=sp['nskies']))
    for win, (lo, hi, lines) in cfg['windows'].items():
        grid = np.exp(np.arange(np.log(lo), np.log(hi), cfg['dlog'] * np.log(10)))   # rest/helio common grid
        sig_A = np.mean(grid) / cfg['R'] / 2.3548
        N = len(S)
        REST = np.full((N, len(grid)), np.nan)
        HEL = np.full((N, len(grid)), np.nan)
        HIV = np.zeros((N, len(grid)))
        for i, sp in enumerate(specs):
            f, iv = window_arrays(sp, lo, hi, grid)
            if f is None or np.sum(iv > 0) < 0.7 * len(grid):
                continue
            med = np.nanmedian(f[iv > 0])
            if not np.isfinite(med) or med <= 0:
                continue
            HEL[i] = f / med
            HIV[i] = iv * med ** 2
            REST[i] = np.interp(grid, grid / (1 + rv[i] / C_KMS), f / med)       # to stellar rest frame
        off = np.where((np.array(S['sample']) == 'off') & np.isfinite(REST[:, 0]))[0]
        # optional: blank each template star's own nebular-line pixels.  Off by default: at R~1800 the blanked
        # window (+-6 A) covers the stellar Balmer cores and interpolating across them biases Hb/Ha badly.
        RESTM = REST.copy()
        if MASK_TEMPLATES:
            for j in range(N):
                RESTM[j, ~line_mask(grid, lines, rv[j], sig_A)] = np.nan
        resid_stack = []
        for i in range(N):
            if not np.isfinite(HEL[i, 0]):
                continue
            keep = line_mask(grid, lines, rv[i], sig_A)
            pool = off[off != i]
            if mode == 'mrs':      # same plate templates first (same instrument setting), else all
                same = pool[np.array([meta[j]['planid'] == meta[i]['planid'] for j in pool])]
                pool = same if len(same) > 3 * K else pool
            d2 = np.nanmean((REST[pool][:, keep] - REST[i][keep]) ** 2, axis=1)
            nb = pool[np.argsort(d2)[:K]]
            tmpl_rest = np.nanmedian(RESTM[nb], axis=0)
            bad = ~np.isfinite(tmpl_rest)
            if bad.any():
                tmpl_rest[bad] = np.interp(grid[bad], grid[~bad], tmpl_rest[~bad])
            tmpl = np.interp(grid, grid * (1 + rv[i] / C_KMS), tmpl_rest)        # back to target (helio) frame
            y, ivy = HEL[i], HIV[i]
            xx = (grid - grid.mean()) / (hi - lo)
            cols = [tmpl, tmpl * xx]
            for ln in lines:
                c = VAC[ln] * (1 + V0 / C_KMS)
                cols.append(gauss(grid, c, c / cfg['R'] / 2.3548 * np.sqrt(1 + (cfg['R'] * 19 / C_KMS) ** 2)))
            good = ivy > 0
            coef, cov, chi2, model = linfit(y[good], ivy[good], [c[good] for c in cols])
            rchi = max(chi2 / max(good.sum() - len(coef), 1), 1)
            cont = coef[0] * tmpl + coef[1] * tmpl * xx
            res = dict(i=i)
            for j, ln in enumerate(lines):
                res[ln] = coef[2 + j]                          # in units of (continuum-normalised flux) x A = EW
                res[ln + '_err'] = np.sqrt(cov[2 + j, 2 + j] * rchi)
            res['rchi2'] = chi2 / max(good.sum() - len(coef), 1)
            results[win].append(res)
            r = np.where(good, y / np.where(cont > 0, cont, np.nan) - 1, np.nan)
            resid_stack.append((i, r, ivy * cont ** 2))
        # stacks in the heliocentric frame (EW units, i.e. fraction of continuum)
        plans = sorted({m['planid'] for m in meta}) if mode == 'mrs' else []
        for samp, plan in [(s_, None) for s_ in ('on', 'off')] + [(s_, p_) for p_ in plans for s_ in ('on', 'off')]:
            idx = [k for k, (i, r, w) in enumerate(resid_stack) if meta[i]['sample'] == samp
                   and (plan is None or meta[i]['planid'] == plan)]
            if len(idx) < 10:
                continue
            R_ = np.array([resid_stack[k][1] for k in idx])
            W_ = np.array([resid_stack[k][2] for k in idx])
            # convert each residual to Rayleigh-per-A units using the star's photometric scaling
            scale = []
            for k in idx:
                i = resid_stack[k][0]
                fl = phot_flambda(S[i], np.mean(grid), mode)
                fap = aperture_fraction(meta[i]['seeing'])
                scale.append(fl * fap / rayleigh_per_fibre_flux(np.mean(grid), 3.3))
            scale = np.array(scale)
            okk = np.isfinite(scale)
            Rr = R_[okk] * scale[okk][:, None]
            Wr = W_[okk] / scale[okk][:, None] ** 2
            # robust: clip 4-sigma per pixel
            med = np.nanmedian(Rr, 0)
            mad = 1.4826 * np.nanmedian(np.abs(Rr - med), 0)
            Wr = np.where(np.abs(Rr - med) < 5 * mad, Wr, 0)
            Wr = np.where(np.isfinite(Rr), Wr, 0)
            num = np.nansum(np.nan_to_num(Rr) * Wr, 0)
            den = np.sum(Wr, 0)
            mean = num / np.where(den > 0, den, np.nan)
            # error from bootstrap over stars
            rng = np.random.default_rng(1)
            bs = []
            for _ in range(100):
                b = rng.integers(0, len(Rr), len(Rr))
                bs.append(np.nansum(np.nan_to_num(Rr[b]) * Wr[b], 0) / np.where(Wr[b].sum(0) > 0, Wr[b].sum(0), np.nan))
            err = np.nanstd(bs, 0)
            stacks[f'{win}_{samp}' + (f'_{plan}' if plan else '')] = dict(grid=grid.tolist(), mean=mean.tolist(), err=err.tolist(), n=int(okk.sum()))
        print(win, 'done', len(results[win]))
    # per-star table
    rows = []
    for k, m in enumerate(meta):
        r = dict(m)
        fl_ok = True
        for win, (lo, hi, lines) in cfg['windows'].items():
            rr = [x for x in results[win] if x['i'] == k]
            if not rr:
                continue
            rr = rr[0]
            for ln in lines:
                lam = VAC[ln]
                fl = phot_flambda(S[k], lam, mode)
                fap = aperture_fraction(m['seeing'])
                conv = fl * fap / rayleigh_per_fibre_flux(lam, 3.3)
                r['EW_' + ln] = rr[ln]
                r['EW_' + ln + '_err'] = rr[ln + '_err']
                r['I_' + ln] = rr[ln] * conv
                r['I_' + ln + '_err'] = rr[ln + '_err'] * conv
        rows.append(r)
    keys = sorted(set().union(*[r.keys() for r in rows]))
    T = Table(rows=[{k: r.get(k, np.nan) for k in keys} for r in rows], names=keys)
    T.write(ROOT / f'data/derived/spectra_sdss_lamost_{cfg["prefix"]}_resid_stars.csv', overwrite=True)
    (ROOT / f'data/derived/spectra_sdss_lamost_{cfg["prefix"]}_resid_stacks.json').write_text(json.dumps(stacks))
    print('wrote', len(T))


if __name__ == '__main__':
    main()
