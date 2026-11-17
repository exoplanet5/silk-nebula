"""w4_family_tap.py -- wave-4 `family`: Data Lab TAP queries (DESI DR1) for Monogem-periphery [O III] filaments
outside the Silk box: zpix coadds, fiberassign rows and exposure/GFA rows of the tiles that touch the regions.

Service: https://datalab.noirlab.edu/tap/sync (ADQL). Outputs in data/raw/w4_family_tap/ (+ .query.txt with date).
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_tap.py [zpix|fiberassign|exposure|all] REGION
REGION: south  (RA 99-116, Dec -2..+11.5: Monogem S limb + G206.6+6.1)
        nw     (RA 92-103, Dec 12..27)
        east   (RA 117-127, Dec 2..30: fainter filaments east of the Silk)
"""
import csv
import sys
import time
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/w4_family_tap'
OUT.mkdir(parents=True, exist_ok=True)
TAP = 'https://datalab.noirlab.edu/tap/sync'
BOXES = {'south': ((99.0, 116.0), (-2.0, 11.5)), 'nw': ((92.0, 103.0), (12.0, 27.0)),
         'east': ((117.0, 127.0), (2.0, 30.0))}
ZCOLS = ('targetid, survey, program, spectype, subtype, z, zerr, zwarn, deltachi2, '
         'mean_fiber_ra, mean_fiber_dec, glon, glat, healpix, coadd_fiberstatus, coadd_numexp, '
         'coadd_numtile, coadd_exptime, mean_mjd, min_mjd, max_mjd, mean_psf_to_fiber_specflux, zcat_primary')


def tap(query, outfile, timeout=1800, tries=3):
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
            Path(str(outfile) + '.query.txt').write_text(f'{TAP}\n{time.strftime("%Y-%m-%dT%H:%M:%S")}\n{query}\n')
            print(f'{outfile.name}: {n/1e6:.2f} MB', flush=True)
            return outfile
        except Exception as e:
            print('TAP error', k, e, flush=True)
            time.sleep(10 * (k + 1))
    raise RuntimeError('TAP failed: ' + query[:200])


def tiles_touching(region, pad=1.65):
    (ra0, ra1), (d0, d1) = BOXES[region]
    import numpy as np
    out = []
    for r in csv.DictReader(open(OUT / 'tiles_monogem.csv')):
        ra, dec = float(r['tilera']), float(r['tiledec'])
        c = np.cos(np.radians(dec))
        if ra0 - pad / c <= ra <= ra1 + pad / c and d0 - pad <= dec <= d1 + pad:
            out.append(int(r['tileid']))
    return out


def main():
    what, region = sys.argv[1], sys.argv[2]
    (ra0, ra1), (d0, d1) = BOXES[region]
    if what in ('zpix', 'all'):
        tap(f'SELECT {ZCOLS} FROM desi_dr1.zpix WHERE mean_fiber_ra BETWEEN {ra0} AND {ra1} '
            f'AND mean_fiber_dec BETWEEN {d0} AND {d1}', OUT / f'zpix_{region}.csv')
    ids = tiles_touching(region)
    print(region, 'tiles', len(ids))
    if what in ('fiberassign', 'all') and ids:
        tap('SELECT targetid, tileid, petal_loc, location, fiber, fa_type, fiberstatus, target_ra, target_dec '
            f'FROM desi_dr1.fiberassign WHERE tileid IN ({",".join(map(str, ids))})',
            OUT / f'fiberassign_{region}.csv')
    if what in ('exposure', 'all') and ids:
        tap('SELECT expid, tileid, night, mjd, exptime, efftime_spec, airmass, seeing_etc, seeing_gfa, '
            'fiber_fracflux_gfa, transparency_gfa, survey, program FROM desi_dr1.exposure '
            f'WHERE tileid IN ({",".join(map(str, ids))})', OUT / f'exposure_{region}.csv')


if __name__ == '__main__':
    main()
