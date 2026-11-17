"""Query SDSS DR17 SkyServer (SQL REST) for every spectroscopic fibre in the Silk Nebula region.

Topic tag: spectra_sdss_lamost.  Access date 2026-11-16.
Endpoint: https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch (table SpecObjAll).
The server rejects long multi-column result sets (spurious 'Arithmetic overflow'), so the columns
are fetched in groups and joined on (plate, mjd, fiberid).

Output: data/raw/sdss_dr17_skyserver/specobjall_box.csv  (+ the raw group CSVs)
"""
import io
import sys
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table, join, vstack

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/sdss_dr17_skyserver'
OUT.mkdir(parents=True, exist_ok=True)
URL = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'

RA = (108.5, 124.5)          # search box padded by one plate radius (1.49 deg) beyond RA 111-122, Dec 6-28
DEC_EDGES = [4.5, 9.5, 14.5, 19.5, 24.5, 29.5]
KEY = 'p.plate, p.mjd, p.fiberid'
GROUPS = {
    'a': 'p.ra, p.dec, p.targettype, p.sourcetype, p.class, p.subclass',
    'b': 'p.z, p.zwarning, p.snmedian, p.spectrographid',
    'c': 'p.run2d, p.survey, p.programname, p.instrument',
    'd': 'p.xfocal, p.yfocal, p.spectroSkyFlux_g, p.spectroSkyFlux_r',
}


def sql(q, tries=4):
    for i in range(tries):
        try:
            r = requests.get(URL, params={'cmd': q, 'format': 'csv'}, timeout=300)
            if r.status_code == 200 and r.text.startswith('#Table1'):
                return Table.read(r.text.split('\n', 1)[1], format='ascii.csv', fast_reader=False)
            print('  bad response', r.status_code, r.text[:200])
        except requests.RequestException as e:
            print('  error', e)
        time.sleep(5 * (i + 1))
    raise RuntimeError('query failed: ' + q[:120])


def main():
    strips = []
    for d0, d1 in zip(DEC_EDGES[:-1], DEC_EDGES[1:]):
        hi = d1 - 1e-7 if d1 < DEC_EDGES[-1] else d1
        where = f'from specobjall p where p.ra between {RA[0]} and {RA[1]} and p.dec between {d0} and {hi:.7f}'
        tabs = []
        for g, cols in GROUPS.items():
            t = sql(f'select {KEY}, {cols} {where}')
            t.write(OUT / f'specobjall_dec{d0:04.1f}_{g}.csv', overwrite=True)
            tabs.append(t)
        m = tabs[0]
        for t in tabs[1:]:
            m = join(m, t, keys=['plate', 'mjd', 'fiberid'])
        print(f'Dec {d0}-{d1}: {len(m)} fibres')
        strips.append(m)
    allt = vstack(strips)
    allt.write(OUT / 'specobjall_box.csv', overwrite=True)
    print('total', len(allt), 'sky fibres', int(np.sum(allt['targettype'] == 'SKY')))


if __name__ == '__main__':
    sys.exit(main())
