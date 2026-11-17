"""Download SDSS-V DR19 BOSS spec-lite files of SKY fibres near the Silk Nebula strand.
Topic tag spectra_sdss_lamost; access date 2026-11-16.

Catalogue: SkyServer DR19 table spAll (objtype='SKY'), saved as
           data/raw/sdss_dr19_skyserver/dr19_spall_sky_box.csv
Files    : https://dr19.sdss.org/sas/dr19/spectro/boss/redux/v6_1_3/spectra/lite/<field6>/<mjd>/spec-<field6>-<mjd>-<catalogid>.fits
           (data.sdss5.org redirects to data.sdss.org, where DR19 paths 404; dr19.sdss.org serves them)
"""
import concurrent.futures as cf
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/sdss_dr19_speclite'
BASE = 'https://dr19.sdss.org/sas/dr19/spectro/boss/redux/v6_1_3/spectra/lite'


def fetch(r):
    field, mjd, cid = int(r['field']), int(r['mjd']), int(r['catalogid'])
    d = OUT / f'{field:06d}'
    d.mkdir(parents=True, exist_ok=True)
    name = f'spec-{field:06d}-{mjd}-{cid}.fits'
    f = d / name
    if f.exists() and f.stat().st_size > 50000:
        return 'cached'
    for i in range(4):
        try:
            x = requests.get(f'{BASE}/{field:06d}/{mjd}/{name}', timeout=90)
            if x.status_code == 200 and x.content[:6] == b'SIMPLE':
                f.write_bytes(x.content)
                return 'ok'
            if x.status_code == 404:
                return '404'
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    return 'fail'


def main():
    t = Table.read(ROOT / 'data/raw/sdss_dr19_skyserver/dr19_spall_sky_box.csv', format='ascii.csv', comment='#')
    m = (t['fiber_ra'] > 112) & (t['fiber_ra'] < 121) & (t['fiber_dec'] > 12) & (t['fiber_dec'] < 24)
    rows = list(t[m])
    print(len(rows), 'files')
    st = {}
    with cf.ThreadPoolExecutor(4) as ex:
        for i, s in enumerate(ex.map(fetch, rows)):
            st[s] = st.get(s, 0) + 1
            if i % 200 == 0:
                print(i, st, flush=True)
    print('done', st)


if __name__ == '__main__':
    main()
