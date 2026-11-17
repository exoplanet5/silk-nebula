"""Download SDSS DR17 spec-lite files of SKY fibres (and a few objects) in the Silk Nebula region.

Topic tag: spectra_sdss_lamost.  Access date 2026-11-16.
SAS DR17 paths (data model: https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/):
  SDSS-I/II (run2d 26, 103, 104):  https://data.sdss.org/sas/dr17/sdss/spectro/redux/<run2d>/spectra/lite/<pppp>/spec-<pppp>-<mjd>-<ffff>.fits
  BOSS/eBOSS (run2d v5_13_2):     https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/<plate>/spec-<plate>-<mjd>-<ffff>.fits
Each spec-lite COADD HDU has columns flux, loglam (vacuum, heliocentric), ivar, and_mask, or_mask,
wdisp, sky, model; flux and sky in 1e-17 erg/s/cm^2/A per fibre.  For a SKY fibre,
flux = sky-subtracted residual and flux + sky = the fibre's total sky spectrum.

usage: python -I spectra_sdss_lamost_get_sdss_sky.py [--objects N]
"""
import argparse
import concurrent.futures as cf
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
CAT = ROOT / 'data/raw/sdss_dr17_skyserver/specobjall_box.csv'
OUT = ROOT / 'data/raw/sdss_dr17_speclite'
BOX = dict(ra=(110.0, 123.0), dec=(5.0, 29.0))


def url_for(plate, mjd, fiber, run2d):
    run2d = str(run2d)
    if run2d.startswith('v5'):
        return (f'https://data.sdss.org/sas/dr17/eboss/spectro/redux/{run2d}/spectra/lite/'
                f'{plate:04d}/spec-{plate:04d}-{mjd}-{fiber:04d}.fits')
    return (f'https://data.sdss.org/sas/dr17/sdss/spectro/redux/{run2d}/spectra/lite/'
            f'{plate:04d}/spec-{plate:04d}-{mjd}-{fiber:04d}.fits')


def fetch(row):
    plate, mjd, fiber, run2d = int(row['plate']), int(row['mjd']), int(row['fiberid']), row['run2d']
    d = OUT / f'{plate:04d}'
    d.mkdir(parents=True, exist_ok=True)
    f = d / f'spec-{plate:04d}-{mjd}-{fiber:04d}.fits'
    if f.exists() and f.stat().st_size > 50000:
        return 'cached', f
    u = url_for(plate, mjd, fiber, run2d)
    for i in range(4):
        try:
            r = requests.get(u, timeout=90)
            if r.status_code == 200 and r.content[:6] == b'SIMPLE':
                f.write_bytes(r.content)
                return 'ok', f
            if r.status_code == 404:
                return '404 ' + u, None
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    return 'fail ' + u, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--objects', type=int, default=0,
                    help='also fetch N bright STAR fibres per plate-mjd (object residual method)')
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    t = Table.read(CAT)
    inbox = ((t['ra'] >= BOX['ra'][0]) & (t['ra'] <= BOX['ra'][1]) &
             (t['dec'] >= BOX['dec'][0]) & (t['dec'] <= BOX['dec'][1]))
    sel = t[inbox & (t['targettype'] == 'SKY')]
    rows = list(sel)
    if a.objects:
        st = t[inbox & (t['class'] == 'STAR') & (t['targettype'] != 'SKY')]
        for key in set(zip(st['plate'], st['mjd'])):
            s = st[(st['plate'] == key[0]) & (st['mjd'] == key[1])]
            s = s[np.argsort(-np.nan_to_num(s['snmedian']))][:a.objects]
            rows += list(s)
    # fetch plates nearest the filament centre first
    c0 = np.array([116.68, 16.42])
    cen = {}
    for r in rows:
        cen.setdefault((int(r['plate']), int(r['mjd'])), []).append((float(r['ra']), float(r['dec'])))
    cen = {k: np.mean(v, axis=0) for k, v in cen.items()}

    def pdist(r):
        ra, dec = cen[(int(r['plate']), int(r['mjd']))]
        return ((ra - c0[0]) * np.cos(np.radians(16.4))) ** 2 + (dec - c0[1]) ** 2
    rows.sort(key=pdist)
    print(len(rows), 'files to fetch')
    stats = {}
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for i, (st_, f) in enumerate(ex.map(fetch, rows)):
            k = st_.split()[0]
            stats[k] = stats.get(k, 0) + 1
            if k not in ('ok', 'cached'):
                print(st_)
            if i % 500 == 0:
                print(i, stats, flush=True)
    print('done', stats)


if __name__ == '__main__':
    main()
