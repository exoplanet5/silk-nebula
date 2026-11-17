"""Fetch DSS1/DSS2 blue and red tiles (SkyView, 6"/pix, Clip resampling) along the Silk Nebula ridge.

Tiles are 2x2 deg TAN images centred at ridge RA -/+ 1 deg (two columns) for Dec = 11, 13, ..., 27.
Bands: DSS2 Blue (POSS-II IIIaJ+GG385, includes [OIII] 5007), DSS2 Red (POSS-II IIIaF+RG610, includes Halpha),
       DSS1 Blue (POSS-I 103aO), DSS1 Red (POSS-I 103aE+red plexiglass, includes Halpha).
Output: data/raw/skyview/<band>/dss_<band>_d<dec>_c<col>.fits + .json (query URL, access time, plates used)
Usage:  python imaging_optical_fetch_dss.py <band> [<band> ...]      e.g. 'DSS2 Blue'
"""
import sys, json, time, datetime as dt
from pathlib import Path
import numpy as np
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/skyview'
BASE = 'https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl'


def ridge_ra(dec):
    r = json.loads((ROOT / 'data/derived/imaging_optical_nsns_ridge.json').read_text())
    return float(np.polyval(r['poly_ra_of_dec'], np.clip(dec, *r['dec_range'])))


def tiles():
    out = []
    for dec in range(11, 28, 2):
        ra0 = ridge_ra(dec)
        for col, off in (('E', +1.0), ('W', -1.0)):
            out.append((dec, col, ra0 + off / np.cos(np.radians(dec))))
    return out


def fetch_band(band):
    d = OUT / band.replace(' ', '_'); d.mkdir(parents=True, exist_ok=True)
    for dec, col, ra in tiles():
        name = f"dss_{band.replace(' ', '_')}_d{dec:02d}_{col}"
        out = d / f'{name}.fits'
        if out.exists() and out.stat().st_size > 100000:
            continue
        q = dict(Survey=band, Position=f'{ra:.4f},{dec:.4f}', Size='2,2', Pixels='1200,1200', Projection='Tan',
                 Coordinates='J2000', Sampler='Clip', Return='FITS')
        for k in range(3):
            try:
                r = requests.get(BASE, params=q, timeout=900)
                if r.status_code == 200 and r.content[:6] == b'SIMPLE' and len(r.content) > 100000:
                    out.write_bytes(r.content)
                    txt = r.content[:200000].decode('ascii', 'replace')
                    plates = sorted(set(x.split('/')[-1].strip() for x in
                                        [txt[i:i + 80] for i in range(0, len(txt), 80)] if 'Used image' in x))
                    (d / f'{name}.json').write_text(json.dumps(dict(
                        survey=band, url=r.url, accessed=dt.datetime.now(dt.timezone.utc).isoformat(),
                        center=[ra, dec], plates=plates), indent=1))
                    print('ok', name, plates, flush=True)
                    break
                print('bad', name, r.status_code, len(r.content), flush=True)
            except Exception as e:
                print('err', name, e, flush=True)
            time.sleep(10 * (k + 1))
        time.sleep(1)


if __name__ == '__main__':
    for b in sys.argv[1:]:
        fetch_band(b)
