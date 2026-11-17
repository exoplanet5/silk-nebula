"""Fetch wide-field narrowband / photographic images of the Silk Nebula region from NASA SkyView.

Usage:  python imaging_optical_fetch_skyview.py <set>
  set = halpha  -> H-Alpha Comp (Finkbeiner 2003), VTSS Ha/Bg, SHASSA H/C/CC/Sm  (wide maps)
        dss     -> DSS1/DSS2 Blue+Red strips (see DSS_TILES)

All output is written to data/raw/skyview/<survey>/ as FITS, with a JSON log of the exact query URL
and access time. SkyView: https://skyview.gsfc.nasa.gov (runquery.pl batch interface).
"""
import sys, json, time, datetime as dt
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/skyview'
BASE = 'https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl'

# Wide maps centred on the filament (TAN projection, J2000)
WIDE = dict(Position='116.5,17.0', Size='17,24', Projection='Tan', Coordinates='J2000')
HALPHA_SETS = [
    ('H-Alpha Comp', dict(Pixels='510,720', **WIDE)),            # 2'/pix, native 6' FWHM
    ('VTSS Ha',      dict(Pixels='510,720', **WIDE)),
    ('VTSS Bg',      dict(Pixels='510,720', **WIDE)),
    # SHASSA covers Dec < +15: south end of the filament; 48"/pix ~ native 47.6"
    ('SHASSA H',  dict(Position='116.5,10.5', Size='13,11', Pixels='975,825', Projection='Tan', Coordinates='J2000')),
    ('SHASSA C',  dict(Position='116.5,10.5', Size='13,11', Pixels='975,825', Projection='Tan', Coordinates='J2000')),
    ('SHASSA CC', dict(Position='116.5,10.5', Size='13,11', Pixels='975,825', Projection='Tan', Coordinates='J2000')),
    ('SHASSA Sm', dict(Position='116.5,10.5', Size='13,11', Pixels='975,825', Projection='Tan', Coordinates='J2000')),
]


def fetch(survey, params, outdir, name, sampler='Clip', tries=3):
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f'{name}.fits'
    if out.exists() and out.stat().st_size > 10000:
        print('exists', out); return out
    q = dict(Survey=survey, Return='FITS', Sampler=sampler, **params)
    for k in range(tries):
        try:
            r = requests.get(BASE, params=q, timeout=600)
            ok = r.status_code == 200 and r.content[:6] == b'SIMPLE'
            if ok:
                out.write_bytes(r.content)
                log = dict(survey=survey, url=r.url, accessed=dt.datetime.utcnow().isoformat() + 'Z',
                           bytes=len(r.content), sampler=sampler)
                (outdir / f'{name}.json').write_text(json.dumps(log, indent=1))
                print('ok', survey, out, len(r.content)); return out
            print('bad response', survey, r.status_code, r.content[:200])
        except Exception as e:
            print('error', survey, e)
        time.sleep(5 * (k + 1))
    return None


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'halpha'
    if which == 'halpha':
        for s, p in HALPHA_SETS:
            fetch(s, p, OUT / s.replace(' ', '_'), s.replace(' ', '_') + '_wide')
