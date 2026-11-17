"""AstroBin public search API: plate-solved images whose field could contain the Silk Nebula, published before the
SWONS detection (Nov 2023). Endpoint https://www.astrobin.com/api/v2/images/image-search/ (no key needed), filters
coord_ra_min/max, coord_dec_min/max, field_radius_min, date_published_max, text. Saves raw JSON pages to
data/raw/astrobin_search/ and a filtered candidate table to data/derived/imaging_optical_astrobin_candidates.csv.
Usage: python imaging_optical_astrobin_search.py <tag> <text or ''> <field_radius_min> [date_published_max]
"""
import sys, json, time, datetime as dt
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/astrobin_search'
API = 'https://www.astrobin.com/api/v2/images/image-search/'
HDR = {'User-Agent': 'Mozilla/5.0 (research script; low rate)'}


def main(tag, text, rmin, dmax='2023-11-01'):
    OUT.mkdir(parents=True, exist_ok=True)
    p = dict(coord_ra_min=100, coord_ra_max=133, coord_dec_min=0, coord_dec_max=36, field_radius_min=rmin,
             date_published_max=dmax, page_size=100)
    if text:
        p['text'] = text
    page = 1; allr = []; log = []
    while True:
        p['page'] = page
        r = requests.get(API, params=p, headers=HDR, timeout=120)
        log.append(dict(url=r.url, status=r.status_code, accessed=dt.datetime.now(dt.timezone.utc).isoformat()))
        if r.status_code != 200:
            print('http', r.status_code, r.text[:200]); break
        d = r.json()
        (OUT / f'{tag}_p{page:02d}.json').write_text(r.text)
        allr += d['results']
        print(page, d['count'], len(allr), flush=True)
        if not d.get('next') or page >= 40:
            break
        page += 1; time.sleep(1.5)
    (OUT / f'{tag}_log.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), *(sys.argv[4:5]))
