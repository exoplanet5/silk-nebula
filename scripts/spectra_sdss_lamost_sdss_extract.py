"""Extract the nebular-line windows of every downloaded SDSS DR17 SKY-fibre spec-lite file into one
compact archive (topic tag spectra_sdss_lamost).

Input : data/raw/sdss_dr17_speclite/<plate>/spec-<plate>-<mjd>-<fiber>.fits  (see *_get_sdss_sky.py)
        data/raw/sdss_dr17_skyserver/platex_box.csv (seeing, airmass, exptime per plate-mjd)
Output: data/derived/spectra_sdss_lamost_sdss_skyfibres.npz
        arrays on the common SDSS log-lambda grid (1e-4 dex pixels; vacuum, heliocentric):
        flux (residual after pipeline sky subtraction), sky (pipeline sky model), ivar, wdisp,
        plus per-fibre metadata.
"""
import glob
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
WINDOWS = [(4800.0, 5070.0), (5175.0, 5225.0), (5560.0, 5600.0), (6280.0, 6320.0), (6490.0, 6790.0)]
IDX = np.concatenate([np.arange(int(np.ceil(np.log10(a) * 1e4)), int(np.floor(np.log10(b) * 1e4)) + 1)
                      for a, b in WINDOWS])
LOGLAM = IDX * 1e-4


def main():
    files = sorted(glob.glob(str(ROOT / 'data/raw/sdss_dr17_speclite/*/spec-*.fits')))
    files += sorted(glob.glob(str(ROOT / 'data/raw/sdss_dr19_speclite/*/spec-*.fits')))   # SDSS-V DR19 (field, catalogid)
    print(len(files), 'files')
    px = Table.read(ROOT / 'data/raw/sdss_dr17_skyserver/platex_box.csv', format='ascii.csv', comment='#')
    see = {(int(r['plate']), int(r['mjd'])): (float(r['seeing50']), float(r['airmass']), float(r['exptime']))
           for r in px}
    n = len(files)
    npx = len(IDX)
    flux = np.zeros((n, npx), np.float32)
    sky = np.zeros((n, npx), np.float32)
    ivar = np.zeros((n, npx), np.float32)
    wdisp = np.zeros((n, npx), np.float32)
    meta = {k: [] for k in ['plate', 'mjd', 'fiber', 'ra', 'dec', 'objtype', 'instrument', 'spectrograph',
                            'seeing50', 'airmass', 'exptime', 'nexp', 'tai', 'skychi2', 'run2d', 'helio_rv']}
    keep = np.zeros(n, bool)
    for i, f in enumerate(files):
        try:
            h = fits.open(f, memmap=False)
            d = h[1].data
            hd = h[0].header
            s2 = h[2].data
        except Exception as e:                       # corrupt download
            print('skip', f, e)
            continue
        dr19 = 'LOGLAM' in d.names and 'loglam' not in d.names
        if dr19:
            d = {k.lower(): d[k] for k in d.names}
        li = np.round(d['loglam'] * 1e4).astype(int)
        pos = np.searchsorted(li, IDX)
        ok = (pos < len(li)) & (li[np.clip(pos, 0, len(li) - 1)] == IDX)
        pp = pos[ok]
        flux[i, ok] = d['flux'][pp]
        sky[i, ok] = d['sky'][pp]
        ivar[i, ok] = d['ivar'][pp] * (d['and_mask'][pp] == 0)
        wdisp[i, ok] = d['wdisp'][pp]
        if dr19:   # SDSS-V: plate -> field, fibre -> catalogid, single APO BOSS spectrograph (500 fibres)
            plate, mjd, fib = int(s2['FIELD'][0]), int(hd['MJD']), int(s2['CATALOGID'][0])
            boss = True
        else:
            plate, mjd, fib = int(hd['PLATEID']), int(hd['MJD']), int(s2['FIBERID'][0])
            boss = len(d['loglam']) > 4200
        meta['plate'].append(plate); meta['mjd'].append(mjd); meta['fiber'].append(fib)
        meta['ra'].append(float(hd['PLUG_RA'])); meta['dec'].append(float(hd['PLUG_DEC']))
        okey = 'OBJTYPE' if 'OBJTYPE' in s2.names else 'TARGETTYPE'
        meta['objtype'].append(str(s2[okey][0]).strip())
        meta['instrument'].append('BOSS' if boss else 'SDSS')
        meta['spectrograph'].append(1 if (dr19 or fib <= (500 if boss else 320)) else 2)
        sa = see.get((plate, mjd), (np.nan, np.nan, np.nan))
        if dr19:
            sa = (float(hd.get('SEEING50', np.nan)), float(hd.get('AIRMASS', np.nan)), float(hd.get('EXPTIME', np.nan)))
        meta['seeing50'].append(sa[0]); meta['airmass'].append(sa[1]); meta['exptime'].append(sa[2])
        meta['nexp'].append(int(hd.get('NEXP', 0))); meta['tai'].append(float(hd.get('TAI', np.nan) or np.nan))
        meta['skychi2'].append(float(hd.get('SKYCHI2', np.nan)))
        meta['run2d'].append('v6_1_3' if dr19 else str(hd.get('RUN2D', '')))
        meta['helio_rv'].append(float(hd.get('HELIO_RV', np.nan)))
        keep[i] = True
    sel = keep
    out = dict(loglam=LOGLAM, flux=flux[sel], sky=sky[sel], ivar=ivar[sel], wdisp=wdisp[sel])
    for k, v in meta.items():
        out[k] = np.array(v)
    np.savez_compressed(ROOT / 'data/derived/spectra_sdss_lamost_sdss_skyfibres.npz', **out)
    print('saved', sel.sum(), 'fibres;', len(set(zip(out['plate'], out['mjd']))), 'plate-mjds')


if __name__ == '__main__':
    sys.exit(main())
