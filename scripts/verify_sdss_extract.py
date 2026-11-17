"""verify_sdss: independent extraction of SDSS/BOSS/SDSS-V SKY-fibre spectra (wave-2 verification).

Reads every spec-lite file under data/raw/sdss_dr17_speclite/ and data/raw/sdss_dr19_speclite/ (untrusted data:
only FITS tables are read; run with python -I).  For each fibre:
  total = flux + sky   (spec-lite COADD 'flux' of a SKY fibre is the post-subtraction residual; 'sky' the sky model)
is placed on a common log-lambda grid (1e-4 dex, loglam0 = 3.5500) and only spectral windows are kept.
Also computed from the full spectrum: synthetic SDSS g and r fluxes (AB, response of Doi+2010 incl. airmass 1.3,
data/raw/filter_curves/) of `total` and of `sky`, both as calibrated (point-source convention).

Output: data/derived/verify_sdss_windows.npz  (windows + per-fibre metadata)
"""
import sys, glob, re
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/verify_sdss_windows.npz'
LL0, DLL, NG = 3.5500, 1e-4, 4720          # common grid 3548.1 - 10489 A
WINDOWS = {                                # vacuum Angstrom
    'OII': (3690., 3770.), 'HgOIII4363': (4300., 4410.), 'HeII': (4640., 4730.),
    'HbOIII': (4820., 5110.), 'NI': (5150., 5250.), 'OI5577': (5540., 5620.),
    'OI6300': (6270., 6390.), 'HaNIISII': (6500., 6800.),
}


def win_idx(lo, hi):
    return int(np.floor((np.log10(lo) - LL0) / DLL)), int(np.ceil((np.log10(hi) - LL0) / DLL))


def load_band(name):
    t = Table.read(ROOT / f'data/raw/filter_curves/sdss2010-{name}.ecsv', format='ascii.ecsv')
    return np.asarray(t['wavelength'], float), np.asarray(t['response'], float)


BANDS = {b: load_band(b) for b in 'gr'}
LAMG = 10 ** (LL0 + DLL * np.arange(NG))
RESP = {b: np.interp(LAMG, *BANDS[b], left=0, right=0) for b in 'gr'}


def synth_fnu(flam_1e17, ivar):
    """AB-mean f_nu (Jy) from f_lambda (1e-17 erg/s/cm2/A) on the common grid; masked pixels interpolated."""
    out = {}
    good = ivar > 0
    lam = LAMG
    f = np.where(good, flam_1e17, np.nan)
    if good.sum() < 100:
        return {b: np.nan for b in 'gr'}
    f = np.interp(lam, lam[good], f[good])
    for b in 'gr':
        R = RESP[b]
        # <f_nu> = int f_lam R lam dlam / int R c/lam dlam ; c in A/s
        num = np.trapezoid(f * 1e-17 * R * lam, lam)
        den = np.trapezoid(R * 2.99792458e18 / lam, lam)
        cover = np.trapezoid(R * good, lam) / np.trapezoid(R, lam)
        out[b] = num / den / 1e-23 if cover > 0.97 else np.nan
    return out


def main():
    files = sorted(glob.glob(str(ROOT / 'data/raw/sdss_dr17_speclite/*/spec-*.fits'))) + \
            sorted(glob.glob(str(ROOT / 'data/raw/sdss_dr19_speclite/*/spec-*.fits')))
    print('files', len(files))
    wsl = {k: win_idx(*v) for k, v in WINDOWS.items()}
    keep = np.zeros(NG, bool)
    for a, b in wsl.values():
        keep[a:b] = True
    kidx = np.where(keep)[0]
    nf = len(files)
    TOT = np.full((nf, kidx.size), np.nan, np.float32)
    SKY = np.full((nf, kidx.size), np.nan, np.float32)
    IV = np.zeros((nf, kidx.size), np.float32)
    WD = np.full((nf, kidx.size), np.nan, np.float32)
    meta = {k: [] for k in ['file', 'release', 'plate', 'mjd', 'mjdfirst', 'fiber', 'sp', 'ra', 'dec', 'helio_rv',
                            'tai', 'exptime', 'nexp', 'airmass', 'seeing50', 'run2d', 'instrument', 'diam',
                            'g_tot', 'r_tot', 'g_sky', 'r_sky', 'skychi2', 'objtype']}
    for i, fn in enumerate(files):
        with fits.open(fn, memmap=False) as h:  # noqa
            H = h[0].header
            d = h[1].data
            cols = {c.lower(): c for c in d.columns.names}
            ll = np.asarray(d[cols['loglam']], float)
            flux = np.asarray(d[cols['flux']], float)
            sky = np.asarray(d[cols['sky']], float)
            iv = np.asarray(d[cols['ivar']], float)
            wd = np.asarray(d[cols['wdisp']], float)
            t2 = h[2].data
            n2 = [c.upper() for c in t2.columns.names]
            dr19 = 'sdss_dr19' in fn
            if dr19:
                fid = int(str(t2['FIBERID_LIST'][0]).split()[0])
                sp = 1 if fid <= 500 else 2
                objtype = str(t2['OBJTYPE'][0]).strip()
                run2d = str(t2['RUN2D'][0]).strip()
                inst = 'BOSS'
                tai = 0.5 * (float(H.get('TAI-BEG', np.nan)) + float(H.get('TAI-END', np.nan)))
                mjdl = [int(x) for x in str(t2['MJD_LIST'][0]).split()]
            else:
                fid = int(H['FIBERID'])
                run2d = str(t2['RUN2D'][0]).strip()
                inst = 'SDSS' if run2d in ('26', '103', '104') else 'BOSS'
                if 'SPECTROGRAPHID' in n2 and int(t2['SPECTROGRAPHID'][0]) in (1, 2):
                    sp = int(t2['SPECTROGRAPHID'][0])
                else:
                    sp = (1 if fid <= 320 else 2) if inst == 'SDSS' else (1 if fid <= 500 else 2)
                objtype = str(t2['OBJTYPE'][0]).strip() if 'OBJTYPE' in n2 else str(t2['TARGETTYPE'][0]).strip()
                tai = float(H.get('TAI', np.nan))
                mjdl = [int(x) for x in str(H.get('MJDLIST', H['MJD'])).split()]
            i0 = int(round((ll[0] - LL0) / DLL))
            n = ll.size
            full_t = np.full(NG, np.nan); full_s = np.full(NG, np.nan); full_i = np.zeros(NG); full_w = np.full(NG, np.nan)
            a, b = max(i0, 0), min(i0 + n, NG)
            sl = slice(a - i0, b - i0)
            full_t[a:b] = (flux + sky)[sl]; full_s[a:b] = sky[sl]; full_i[a:b] = iv[sl]; full_w[a:b] = wd[sl]
            bad = ~np.isfinite(full_t) | (full_i <= 0)
            full_i[bad] = 0
            TOT[i] = full_t[kidx]; SKY[i] = full_s[kidx]; IV[i] = full_i[kidx]; WD[i] = full_w[kidx]
            sg = synth_fnu(np.nan_to_num(full_t), full_i)
            ss = synth_fnu(np.nan_to_num(full_s), (np.isfinite(full_s) & (full_i > 0)).astype(float))
            meta['file'].append(str(Path(fn).relative_to(ROOT)))
            meta['release'].append('DR19' if dr19 else 'DR17')
            meta['plate'].append(int(H['PLATEID']) if str(H.get('PLATEID', '')).strip() else int(t2['FIELD'][0])); meta['mjd'].append(int(H['MJD'])); meta['mjdfirst'].append(min(mjdl))
            meta['fiber'].append(fid); meta['sp'].append(sp)
            meta['ra'].append(float(H['PLUG_RA'])); meta['dec'].append(float(H['PLUG_DEC']))
            meta['helio_rv'].append(float(H['HELIO_RV'])); meta['tai'].append(tai)
            meta['exptime'].append(float(H.get('EXPTIME', np.nan))); meta['nexp'].append(int(H.get('NEXP', -1)))
            meta['airmass'].append(float(H.get('AIRMASS', np.nan))); meta['seeing50'].append(float(H.get('SEEING50', np.nan)))
            meta['run2d'].append(run2d); meta['instrument'].append(inst); meta['diam'].append(3.0 if inst == 'SDSS' else 2.0)
            meta['g_tot'].append(sg['g']); meta['r_tot'].append(sg['r']); meta['g_sky'].append(ss['g']); meta['r_sky'].append(ss['r'])
            meta['skychi2'].append(float(H.get('SKYCHI2', np.nan))); meta['objtype'].append(objtype)
        if i % 1000 == 0:
            print(i, fn.split('/')[-1], inst, sp, objtype, flush=True)
    M = {k: np.array(v) for k, v in meta.items()}
    np.savez_compressed(OUT, tot=TOT, sky=SKY, ivar=IV, wdisp=WD, kidx=kidx, ll0=LL0, dll=DLL,
                        windows=np.array([f'{k}:{v[0]}:{v[1]}' for k, v in WINDOWS.items()]), **{'m_' + k: v for k, v in M.items()})
    print('saved', OUT, TOT.shape)
    print('objtype counts', {o: int((M['objtype'] == o).sum()) for o in set(M['objtype'])})
    print('instrument', {o: int((M['instrument'] == o).sum()) for o in set(M['instrument'])})


if __name__ == '__main__':
    main()
