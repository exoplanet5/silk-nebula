"""Pan-STARRS1 DR2 stack cutouts (g, r) across the Silk Nebula.

Service: https://ps1images.stsci.edu/cgi-bin/ps1filenames.py + fitscut.cgi (MAST). Stack images are
'stk.<b>.unconv' skycell stacks (0.25"/pix). We request 25' (6000 px) cutouts block-averaged by fitscut to
600 px (2.5"/pix) with format=fits.
Usage: python imaging_optical_fetch_ps1.py test | strips
"""
import sys, json, time, datetime as dt
from pathlib import Path
import numpy as np
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/ps1_dr2'
FN = 'https://ps1images.stsci.edu/cgi-bin/ps1filenames.py'
CUT = 'https://ps1images.stsci.edu/cgi-bin/fitscut.cgi'


def get(url, params, timeout=300, tries=3):
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            print('http', r.status_code, r.url)
        except Exception as e:
            print('err', e)
        time.sleep(5 * (k + 1))
    return None


def cutout(ra, dec, band, name, size=6000, out_size=600):
    d = OUT / 'cutouts'; d.mkdir(parents=True, exist_ok=True)
    out = d / f'{name}_{band}.fits'
    if out.exists() and out.stat().st_size > 10000:
        return out
    r = get(FN, dict(ra=ra, dec=dec, filters=band, type='stack'))
    if r is None:
        return None
    lines = [l.split() for l in r.text.strip().splitlines()[1:]]
    if not lines:
        print('no file', ra, dec); return None
    fname = lines[0][7]
    p = dict(ra=ra, dec=dec, size=size, format='fits', red=fname, output_size=out_size)
    r2 = get(CUT, p, timeout=600)
    if r2 is None or r2.content[:6] != b'SIMPLE':
        print('bad cutout', ra, dec, band, None if r2 is None else r2.content[:200]); return None
    out.write_bytes(r2.content)
    (d / f'{name}_{band}.json').write_text(json.dumps(dict(
        url=r2.url, filename=fname, accessed=dt.datetime.now(dt.timezone.utc).isoformat()), indent=1))
    return out


def strips():
    """Cross-filament strips: 12' cutouts (2880 px -> 288 px, 2.5"/pix) every 8' from ridge-1.7 to ridge+1.7 deg
    (great-circle offsets) at the Dec of P3, P1, P2 and the northern arc. Cutouts are clipped by PS1 to the
    skycell containing their centre; overlapping steps fill the gaps when mosaicked."""
    rid = json.loads((ROOT / 'data/derived/imaging_optical_nsns_ridge.json').read_text())
    STRIPS = {'P3': 14.72, 'P1': 16.60, 'P2': 17.74, 'N22': 22.3}
    for nm, dec in STRIPS.items():
        r0 = float(np.polyval(rid['poly_ra_of_dec'], dec))
        offs = np.arange(-1.7, 1.7001, 8 / 60)
        for k, off in enumerate(offs):
            ra = r0 + off / np.cos(np.radians(dec))
            for b in 'gr':
                o = cutout(round(ra, 5), dec, b, f'strip_{nm}_{k:02d}', size=2880, out_size=288)
                print(nm, k, b, o, flush=True)
                time.sleep(0.3)


if __name__ == '__main__':
    if sys.argv[1] == 'test':
        print(cutout(116.53, 14.72, 'g', 'test_P3'))
        print(cutout(112.2625, 13.2464, 'g', 'medusa_A21'))
        print(cutout(112.2625, 13.2464, 'r', 'medusa_A21'))
    else:
        strips()
