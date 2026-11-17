"""verify_xrayfuv: fetch the eRASS1 Main catalogue (Merloni+2024, VizieR J/A+A/682/A34/erass1-m)
over the Monogem Ring region from the CfA VizieR mirror, in RA strips, for point-source masking
and an empirical exposure map (MLExp1). Logged in MANIFEST.jsonl.
"""
import hashlib, json, time, urllib.parse
from pathlib import Path
import requests

OUT = Path('/Users/mickey/astronomy/silk-nebula/data/raw/verify_xrayfuv_erass1cat')
BASE = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
COLS = 'RA_ICRS,DE_ICRS,GLON,GLAT,ELON,ELAT,MJD,EXT,EXTLike,DetLike0,MLcts1,MLRate1,MLFlux1,MLBkg1,MLExp1'
for ra0 in range(78, 136, 4):
    dest = OUT / f'erass1m_ra{ra0:03d}_{ra0+4:03d}.tsv'
    if dest.exists() and dest.stat().st_size > 1000:
        continue
    p = {'-source': 'J/A+A/682/A34/erass1-m', '-out': COLS, '-out.max': 'unlimited',
         'RA_ICRS': f'{ra0}..{ra0+4}', 'DE_ICRS': '-16..44', '-oc.form': 'd'}
    for attempt in range(3):
        try:
            r = requests.get(BASE, params=p, timeout=900)
            b = r.content
            break
        except Exception as e:
            print('retry', ra0, e); time.sleep(10)
    dest.write_bytes(b)
    n = sum(1 for ln in b.decode('utf-8', 'replace').splitlines() if ln and ln[0].isdigit())
    with open(OUT / 'MANIFEST.jsonl', 'a') as fo:
        fo.write(json.dumps({'access_date': time.strftime('%Y-%m-%d'), 'url': BASE + '?' + urllib.parse.urlencode(p),
                             'file': str(dest.relative_to(OUT.parents[2])), 'bytes': len(b), 'rows': n,
                             'sha256': hashlib.sha256(b).hexdigest()}) + '\n')
    print(ra0, len(b), n, flush=True)
