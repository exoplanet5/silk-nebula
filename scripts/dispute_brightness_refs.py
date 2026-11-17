"""dispute_brightness_refs.py -- verify the references cited in dispute section 5 against Crossref and arXiv.

Queries https://api.crossref.org/works/<DOI> and the arXiv export API; writes the returned bibliographic fields
(authors, year, title, journal, volume, page) to data/derived/dispute_brightness_refs.json. Responses are parsed
as data only. Run: ~/.venvs/astro313/bin/python -I scripts/dispute_brightness_refs.py
"""
import json, re, time, urllib.request, urllib.parse
from pathlib import Path

OUT = Path('/Users/mickey/astronomy/silk-nebula/data/derived/dispute_brightness_refs.json')
UA = {'User-Agent': 'generic-research-client/1.0'}

DOIS = {
    'Wu2026': '10.3847/1538-3881/aeb11c',
    'Fesen2024': '10.3847/1538-4365/ad410a',
    'Ziegenbalg2025': '10.3847/2515-5172/adfec7',
    'Abdurrouf2022': '10.3847/1538-4365/ac4414',
    'GaiaDR3_2023': '10.1051/0004-6361/202243940',
    'BlantonRoweis2007': '10.1086/510127',
    'OkeGunn1983': '10.1086/160817',
    'Jordi2010': '10.1051/0004-6361/201015441',
    'Martin2005': '10.1086/426387',
    'Morrissey2007': '10.1086/520512',
    'Raymond1979': '10.1086/190562',
    'CoxRaymond1985': '10.1086/163649',
    'Raymond1988': '10.1086/165945',
    'SFD1998': '10.1086/305772',
    'Hanuschik2003': '10.1051/0004-6361:20030885',
    'Almeida2023_DR18': '10.3847/1538-4365/acda98',
}
SEARCH = {'DopitaSutherland1996': 'Dopita Sutherland 1996 Spectral signatures of fast shocks Astrophysical Journal Supplement 102 161'}
ARXIV = ['2610.13329', '2608.26666', '2503.14745', '2507.07093']


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=40) as f:
        return f.read().decode('utf-8', 'replace')


def summ(m):
    return dict(doi=m.get('DOI'), authors=[a.get('family', a.get('name', '?')) for a in m.get('author', [])[:4]],
                n_authors=len(m.get('author', [])), year=(m.get('issued', {}).get('date-parts') or [[None]])[0][0],
                title=re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0])).strip(),
                journal=(m.get('container-title') or [''])[0], volume=m.get('volume'),
                page=m.get('page') or m.get('article-number'))


out = {'crossref': {}, 'arxiv': {}, 'accessed': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
for k, doi in DOIS.items():
    try:
        out['crossref'][k] = summ(json.loads(get('https://api.crossref.org/works/' + urllib.parse.quote(doi)))['message'])
    except Exception as e:
        out['crossref'][k] = {'doi': doi, 'error': str(e)}
    time.sleep(0.5)
for k, q in SEARCH.items():
    try:
        u = 'https://api.crossref.org/works?' + urllib.parse.urlencode({'query.bibliographic': q, 'rows': 1})
        out['crossref'][k] = summ(json.loads(get(u))['message']['items'][0])
    except Exception as e:
        out['crossref'][k] = {'query': q, 'error': str(e)}
try:
    t = get('https://export.arxiv.org/api/query?id_list=' + ','.join(ARXIV))
    for e in t.split('<entry>')[1:]:
        i = re.search(r'<id>(.*?)</id>', e).group(1)
        out['arxiv'][i] = dict(title=' '.join(re.search(r'<title>(.*?)</title>', e, re.S).group(1).split()),
                               authors=re.findall(r'<name>(.*?)</name>', e)[:4],
                               published=re.search(r'<published>(.*?)</published>', e).group(1))
except Exception as e:
    out['arxiv']['error'] = str(e)
json.dump(out, open(OUT, 'w'), indent=1)
for k, v in out['crossref'].items():
    print(k, '|', v.get('authors'), v.get('year'), '|', (v.get('title') or '')[:70], '|', v.get('journal'), v.get('volume'), v.get('page'), v.get('doi'), v.get('error', ''))
for k, v in out['arxiv'].items():
    print(k, v)
