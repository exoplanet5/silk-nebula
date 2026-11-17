"""verify_xrayfuv: fetch wide-field HI4PI N(HI) and SFD E(B-V) around the Monogem Ring from SkyView.

Galactic CAR, centre (l,b)=(203,+13), 50 x 50 deg, 0.2 deg pixels. Logged in MANIFEST.jsonl.
"""
import hashlib, json, time, urllib.parse
from pathlib import Path
import requests

OUT = Path('/Users/mickey/astronomy/silk-nebula/data/raw/verify_xrayfuv_skyview')
SV = 'https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl'
SURVEYS = {'RASS Background 1': 'rass_r1', 'RASS Background 2': 'rass_r2', 'HI4PI': 'hi4pi_nhi', 'SFD Dust Map': 'sfd_ebv', 'RASS Background 4': 'rass_r4',
           'RASS Background 5': 'rass_r5', 'RASS Background 6': 'rass_r6', 'RASS Background 7': 'rass_r7'}
for sv, tag in SURVEYS.items():
    p = {'Position': '203,13', 'Survey': sv, 'Coordinates': 'Galactic', 'Projection': 'Car',
         'Pixels': '250,250', 'Size': '50,50', 'Sampler': 'Clip', 'ClipIntensive': '', 'Return': 'FITS'}
    dest = OUT / f'wide_{tag}_l203b13_car.fits'
    if dest.exists():
        continue
    r = requests.get(SV, params=p, timeout=900)
    b = r.content
    ok = b[:6] == b'SIMPLE'
    dest = dest if ok else dest.with_suffix('.INVALID.txt')
    dest.write_bytes(b)
    with open(OUT / 'MANIFEST.jsonl', 'a') as fo:
        fo.write(json.dumps({'access_date': time.strftime('%Y-%m-%d'), 'url': SV + '?' + urllib.parse.urlencode(p),
                             'file': str(dest.relative_to(OUT.parents[2])), 'bytes': len(b),
                             'sha256': hashlib.sha256(b).hexdigest(), 'ok': ok}) + '\n')
    print(sv, ok, len(b))
