"""Download NSNS DR0.2 FITS HiPS tiles (linear intensity, Rayleighs) covering the Silk Nebula region.

NSNS = Northern Sky Narrowband Survey by S. Ziegenbalg, https://www.simg.de/nebulae3/dr0_2 (CC-BY-NC-SA 4.0).
HiPS IDs ivo://simg.de/P/NSNS/DR0_2/{oiii,halpha,sii,vc,sn-oiii,...}; tiles 512x512, max order 6.

Usage: python imaging_optical_fetch_nsns.py <order> <hips> [<hips> ...]
"""
import sys, json, time, datetime as dt
from pathlib import Path
import numpy as np
import healpy as hp
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/nsns_hips'
BASE = 'https://simg.de/nebulae3/dr0_2/'
BOX = dict(ra=(106.0, 127.0), dec=(2.0, 31.0))


def tiles_for_box(order):
    nside = 2 ** order
    ra0, ra1 = BOX['ra']; de0, de1 = BOX['dec']
    corners = [(ra0, de0), (ra1, de0), (ra1, de1), (ra0, de1)]
    vec = np.array([hp.ang2vec(r, d, lonlat=True) for r, d in corners])
    return sorted(hp.query_polygon(nside, vec, inclusive=True, nest=True).tolist())


def main(order, names):
    log = []
    for name in names:
        d = OUT / name / f'Norder{order}'
        d.mkdir(parents=True, exist_ok=True)
        props = OUT / name / 'properties'
        if not props.exists():
            r = requests.get(BASE + name + '/properties', timeout=60)
            props.write_bytes(r.content)
        for n in tiles_for_box(order):
            D = (n // 10000) * 10000
            url = f'{BASE}{name}/Norder{order}/Dir{D}/Npix{n}.fits'
            out = d / f'Npix{n}.fits'
            if out.exists() and out.stat().st_size > 1000:
                continue
            for k in range(3):
                try:
                    r = requests.get(url, timeout=120)
                    if r.status_code == 200 and r.content[:6] == b'SIMPLE':
                        out.write_bytes(r.content); status = 'ok'; break
                    status = f'http {r.status_code}'
                    if r.status_code == 404:
                        break
                except Exception as e:
                    status = f'error {e}'
                time.sleep(3 * (k + 1))
            log.append(dict(url=url, status=status, accessed=dt.datetime.now(dt.timezone.utc).isoformat()))
            print(status, url)
            time.sleep(0.3)
    (OUT / f'fetch_log_order{order}_{"_".join(names)}.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2:])
