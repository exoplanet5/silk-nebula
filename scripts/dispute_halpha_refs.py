"""dispute_halpha_refs.py -- verify the references cited in Section 4 of the Comment (tag: halpha).

Queries https://api.crossref.org/works/<DOI>, the arXiv export API and the DataCite API, and writes the returned
bibliographic fields to data/derived/dispute_halpha_refs.json. Responses are parsed as data only.
The Kurucz (2005) solar flux atlas has no DOI; its identity is taken from the local file header
data/raw/verify_velocity_solar/kurucz2005_header.txt (downloaded 2026-11-16 by verify_velocity).
Run: ~/.venvs/astro313/bin/python -I scripts/dispute_halpha_refs.py
"""
import json
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/dispute_halpha_refs.json'
UA = {'User-Agent': 'generic-research-client/1.0'}

CROSSREF = {
    'Wu2026': '10.3847/1538-3881/aeb11c',
    'Zhang2021': '10.1088/1674-4527/21/11/280',
    'Wu2022': '10.1088/1674-4527/ac7387',
    'Ren2021': '10.1088/1674-4527/21/3/51',
    'Haffner2003': '10.1086/378850',
    'ZhangLY2025': '10.3847/1538-3881/ae0f94',
    'Wenger2018': '10.3847/1538-4357/aaaec8',
    'Reid2014': '10.1088/0004-637X/783/2/130',
    'Astropy2022': '10.3847/1538-4357/ac7c74',
}
ARXIV = {'Wu2026': '2610.13329', 'Zhang2021': '2108.08021', 'Wu2022': '2110.04547', 'Ren2021': '2008.06236',
         'ZhangLY2025': '2510.08950'}
DATACITE = {'NIST_ASD': '10.18434/T4W30F'}


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=30) as r:
        return r.read().decode('utf-8', 'replace')


def main():
    out = {'_retrieved': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    for k, doi in CROSSREF.items():
        try:
            m = json.loads(get('https://api.crossref.org/works/' + doi))['message']
            out[k] = dict(source='crossref', doi=m.get('DOI'),
                          authors=[(a.get('family', '') + ', ' + a.get('given', '')).strip(', ')
                                   for a in m.get('author', [])[:6]] or [a.get('name') for a in m.get('author', [])[:3]],
                          title=re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0])).strip(),
                          journal=(m.get('container-title') or [''])[0], volume=m.get('volume'),
                          page=m.get('page') or m.get('article-number'),
                          year=(m.get('issued', {}).get('date-parts') or [[None]])[0][0])
        except Exception as e:                                  # noqa: BLE001
            out[k] = dict(source='crossref', doi=doi, error=str(e))
        time.sleep(0.5)
    try:
        xml = get('https://export.arxiv.org/api/query?id_list=' + ','.join(ARXIV.values()) + '&max_results=10')
        for entry in xml.split('<entry>')[1:]:
            aid = re.search(r'<id>http://arxiv.org/abs/([^<v]+)', entry).group(1)
            title = re.sub(r'\s+', ' ', re.search(r'<title>(.*?)</title>', entry, re.S).group(1)).strip()
            doi = re.search(r'<arxiv:doi>(.*?)</arxiv:doi>', entry)
            key = [k for k, v in ARXIV.items() if v == aid]
            if key:
                out[key[0]]['arxiv'] = dict(id=aid, title=title, doi=doi.group(1) if doi else None)
    except Exception as e:                                      # noqa: BLE001
        out['arxiv_error'] = str(e)
    for k, doi in DATACITE.items():
        try:
            a = json.loads(get('https://api.datacite.org/dois/' + doi))['data']['attributes']
            out[k] = dict(source='datacite', doi=a['doi'], title=a['titles'][0]['title'],
                          creators=[c.get('name') for c in a['creators']], year=a.get('publicationYear'),
                          note='Version cited (5.12, 2024) taken from the "how to reference" line of the NIST ASD '
                               'H I output pages retrieved 2026-11-16')
        except Exception as e:                                  # noqa: BLE001
            out[k] = dict(source='datacite', doi=doi, error=str(e))
    hdr = (ROOT / 'data/raw/verify_velocity_solar/kurucz2005_header.txt').read_text().splitlines()[:4]
    out['Kurucz2005_atlas'] = dict(source='local file header', header=hdr,
                                   url='http://kurucz.harvard.edu/sun/fluxatlas2005/solarfluxintwl.asc')
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in out.items():
        print(k, json.dumps(v, ensure_ascii=False)[:300])


if __name__ == '__main__':
    main()
