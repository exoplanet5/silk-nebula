"""Download LAMOST DR12 v2.0 spectra (public international release, end of Sept 2026) for the
Silk Nebula field.  Topic tag spectra_sdss_lamost; access date 2026-11-16.

Endpoints (no login needed for public data):
  SQL      : https://www.lamost.org/dr12/v2.0/sql/q   (form fields sql=..., output.fmt=csv)
  LRS FITS : https://www.lamost.org/dr12/v2.0/spectrum/fits/<obsid>
  MRS FITS : https://www.lamost.org/dr12/v2.0/medspectrum/fits/<obsid>

usage:
  python -I spectra_sdss_lamost_get_lamost.py mrs-plan PLANID [--max N]
  python -I spectra_sdss_lamost_get_lamost.py mrs-list CSV            (column obsid)
  python -I spectra_sdss_lamost_get_lamost.py lrs-list CSV            (column obsid)
"""
import concurrent.futures as cf
import sys
import time
from pathlib import Path

import requests
from astropy.table import Table, unique

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/lamost_dr12'
BASE = 'https://www.lamost.org/dr12/v2.0'


def sql(q, out):
    for i in range(4):
        try:
            r = requests.post(f'{BASE}/sql/q', data={'sql': q, 'output.fmt': 'csv'}, timeout=600)
            if r.status_code == 200 and not r.text.lstrip().startswith('<'):
                Path(out).write_text(r.text.lstrip('﻿'))
                return Table.read(out, format='ascii.csv')
            print('bad', r.status_code, r.text[:200])
        except requests.RequestException as e:
            print('err', e)
        time.sleep(5 * (i + 1))
    raise RuntimeError(q)


def fetch(args):
    kind, obsid, outdir = args
    outdir.mkdir(parents=True, exist_ok=True)
    f = outdir / f'{kind}-{obsid}.fits.gz'
    if f.exists() and f.stat().st_size > 10000:
        return 'cached'
    url = f'{BASE}/{"medspectrum" if kind == "mrs" else "spectrum"}/fits/{obsid}'
    for i in range(4):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 200 and r.content[:2] == b'\x1f\x8b':
                f.write_bytes(r.content)
                return 'ok'
            if r.status_code == 404:
                return '404'
        except requests.RequestException:
            pass
        time.sleep(3 * (i + 1))
    return 'fail'


def run(kind, obsids, outdir, workers=4):
    stats = {}
    with cf.ThreadPoolExecutor(workers) as ex:
        for i, s in enumerate(ex.map(fetch, [(kind, o, outdir) for o in obsids])):
            stats[s] = stats.get(s, 0) + 1
            if i % 200 == 0:
                print(i, stats, flush=True)
    print('done', stats)


def main():
    mode = sys.argv[1]
    if mode == 'mrs-plan':
        plan = sys.argv[2]
        t = sql(f"select obsid, mobsid, planid, obsdate, lmjm, spid, fiberid, ra, dec, snr, fibertype, band, "
                f"subproject, gaia_g_mean_mag, gaia_bp_mean_mag, gaia_rp_mean_mag, rv_br1, rv_br1_err, rv_lasp1 "
                f"from med_catalogue where planid = '{plan}'", RAW / f'med_catalogue_{plan}.csv')
        u = unique(t, keys='obsid')
        run('mrs', [int(o) for o in u['obsid']], RAW / 'mrs' / plan)
    elif mode in ('mrs-list', 'lrs-list'):
        t = Table.read(sys.argv[2], format='ascii.csv')
        u = unique(t, keys='obsid')
        kind = mode.split('-')[0]
        outdir = RAW / kind / Path(sys.argv[2]).stem
        run(kind, [int(o) for o in u['obsid']], outdir, workers=int(sys.argv[3]) if len(sys.argv) > 3 else 4)


if __name__ == '__main__':
    main()
