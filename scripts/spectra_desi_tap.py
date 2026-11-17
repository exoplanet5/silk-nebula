"""spectra_desi_tap.py -- Data Lab TAP queries for the DESI DR1 Silk Nebula residual-emission study.

Downloads (CSV, to data/raw/datalab_tap/):
  zpix_box.csv         desi_dr1.zpix rows in RA 111-122, Dec +6..+28 (healpix coadds = SPARCL spectra)
  fiberassign_box.csv  desi_dr1.fiberassign rows (all fa_types incl. sky) for tiles touching the box
  exposure_box.csv     desi_dr1.exposure rows for those tiles
  tiles_box.csv        desi_dr1.tile rows around the box

Service: https://datalab.noirlab.edu/tap/sync (TAP, ADQL). Accessed 2026-11-16.
Usage: ~/.venvs/astro313/bin/python -I scripts/spectra_desi_tap.py [zpix|fiberassign|exposure|all]
"""
import sys
import time
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/datalab_tap'
OUT.mkdir(parents=True, exist_ok=True)
TAP = 'https://datalab.noirlab.edu/tap/sync'
BOX = dict(ra=(111.0, 122.0), dec=(6.0, 28.0))


def tap(query, outfile, timeout=900, tries=3):
    for k in range(tries):
        try:
            r = requests.post(TAP, data=dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=query),
                              timeout=timeout, stream=True)
            r.raise_for_status()
            n = 0
            with open(outfile, 'wb') as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
                    n += len(chunk)
            head = open(outfile, 'rb').read(300)
            if b'QUERY_STATUS' in head and b'ERROR' in head:
                raise RuntimeError(head.decode(errors='replace'))
            print(f'{outfile.name}: {n/1e6:.1f} MB')
            return outfile
        except Exception as e:  # retry politely
            print('TAP error', k, e)
            time.sleep(10 * (k + 1))
    raise RuntimeError('TAP failed: ' + query[:200])


ZCOLS = ('targetid, survey, program, spectype, subtype, z, zerr, zwarn, deltachi2, chi2, '
         'mean_fiber_ra, mean_fiber_dec, glon, glat, healpix, coadd_fiberstatus, coadd_numexp, '
         'coadd_numnight, coadd_numtile, coadd_exptime, mean_mjd, min_mjd, max_mjd, '
         'mean_psf_to_fiber_specflux, mean_delta_x, mean_delta_y, tsnr2_lrg, tsnr2_bgs, '
         'zcat_primary, desi_target, bgs_target, mws_target, scnd_target')


def q_zpix():
    q = (f'SELECT {ZCOLS} FROM desi_dr1.zpix WHERE mean_fiber_ra BETWEEN {BOX["ra"][0]} AND {BOX["ra"][1]} '
         f'AND mean_fiber_dec BETWEEN {BOX["dec"][0]} AND {BOX["dec"][1]}')
    tap(q, OUT / 'zpix_box.csv')


def tile_ids():
    import csv
    p = OUT / 'tiles_box.csv'
    if not p.exists():
        tap('SELECT tileid, survey, program, faprgrm, tilera, tiledec, lastnight, nexp, exptime, efftime_spec '
            'FROM desi_dr1.tile WHERE tilera BETWEEN 106 AND 127 AND tiledec BETWEEN 2 AND 32 ORDER BY tilera', p)
    return [int(r['tileid']) for r in csv.DictReader(open(p))]


def q_fiberassign():
    ids = tile_ids()
    q = ('SELECT targetid, tileid, petal_loc, device_loc, location, fiber, fa_type, fa_target, fiberstatus, '
         'target_ra, target_dec, fiberassign_x, fiberassign_y, priority FROM desi_dr1.fiberassign '
         f'WHERE tileid IN ({",".join(map(str, ids))})')
    tap(q, OUT / 'fiberassign_box.csv', timeout=1800)


def q_exposure():
    ids = tile_ids()
    q = ('SELECT expid, tileid, night, mjd, exptime, efftime_spec, airmass, seeing_etc, transparency_gfa, '
         'sky_mag_g_spec, sky_mag_r_spec, tilera, tiledec FROM desi_dr1.exposure '
         f'WHERE tileid IN ({",".join(map(str, ids))})')
    tap(q, OUT / 'exposure_box.csv')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('zpix', 'all'):
        q_zpix()
    if what in ('fiberassign', 'all'):
        q_fiberassign()
    if what in ('exposure', 'all'):
        q_exposure()
