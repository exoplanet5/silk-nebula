"""dispute_references_verify.py -- verify every work cited in the six Comment sections (tag: references).

For each reference the sections cite (notes/dispute/{distance,kinematic,oiii,halpha,brightness,other}.md), this
script queries:
  * https://api.crossref.org/works/<DOI>               (journal articles with a DOI)
  * https://export.arxiv.org/api/query?id_list=<id>    (arXiv identifiers; also returns journal_ref / DOI)
  * https://api.datacite.org/dois/<DOI>                (NIST ASD, Zenodo software DOIs)
  * https://api.crossref.org/works?query.bibliographic= (works without a DOI, to confirm none is registered,
                                                         and arXiv-only works, to look for a journal version)
and compares the returned first-author surname, year, volume and first page / article ID with what the sections
print (EXPECTED below). Responses are parsed as JSON / XML text only; nothing is executed.

Output: data/derived/dispute_references_records.json
Run:    ~/.venvs/astro313/bin/python -I scripts/dispute_references_verify.py
"""
import json
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/dispute_references_records.json'
UA = {'User-Agent': 'generic-research-client/1.0 (reference check)'}

# key: (doi, arxiv, first-author family as printed, year, volume, first page/article as printed)
EXPECTED = {
    'Abdurrouf2022': ('10.3847/1538-4365/ac4414', None, "Abdurro'uf", 2022, '259', '35'),
    'AllendePrieto2001': ('10.1086/322874', 'astro-ph/0106360', 'Allende Prieto', 2001, '556', 'L63'),
    'Almeida2023': ('10.3847/1538-4365/acda98', None, 'Almeida', 2023, '267', '44'),
    'Anderson2012': ('10.1088/0004-637X/754/1/62', '1205.4228', 'Anderson', 2012, '754', '62'),
    'Asplund2004': ('10.1051/0004-6361:20034328', 'astro-ph/0312290', 'Asplund', 2004, '417', '751'),
    'Astropy2022': ('10.3847/1538-4357/ac7c74', None, 'Astropy Collaboration', 2022, '935', '167'),
    'Blanton2007': ('10.1086/510127', 'astro-ph/0606170', 'Blanton', 2007, '133', '734'),
    'Boumis2022': ('10.1093/mnras/stac412', None, 'Boumis', 2022, '512', '1658'),
    'Bovy2017': ('10.1093/mnrasl/slx027', '1610.07610', 'Bovy', 2017, '468', 'L63'),
    'Cioffi1988': ('10.1086/166834', None, 'Cioffi', 1988, '334', '252'),
    'Condon1998a': ('10.1017/S0074180900128177', None, 'Condon', 1998, '179', '19'),
    'Condon1998b': ('10.1086/300337', None, 'Condon', 1998, '115', '1693'),
    'Cox1985': ('10.1086/163649', None, 'Cox', 1985, '298', '651'),
    'Cui2012': ('10.1088/1674-4527/12/9/003', None, 'Cui', 2012, '12', '1197'),
    'DESI2025': (None, '2503.14745', 'DESI Collaboration', 2025, None, None),
    'Dopita1979': ('10.1086/190594', None, 'Dopita', 1979, '40', '455'),
    'Dopita1996': ('10.1086/192255', None, 'Dopita', 1996, '102', '161'),
    'Dopita2016': ('10.1007/s10509-016-2657-8', '1601.01337', 'Dopita', 2016, '361', '61'),
    'Fesen2024': ('10.3847/1538-4365/ad410a', None, 'Fesen', 2024, '272', '36'),
    'Fesen2026': (None, '2608.26666', 'Fesen', 2026, None, None),
    'Gaensler2008': ('10.1071/AS08004', '0808.2550', 'Gaensler', 2008, '25', '184'),
    'Gaia2023': ('10.1051/0004-6361/202243940', None, 'Gaia Collaboration', 2023, '674', 'A1'),
    'Guy2023': ('10.3847/1538-3881/acb212', None, 'Guy', 2023, '165', '144'),
    'Haffner1999': ('10.1086/307734', 'astro-ph/9904143', 'Haffner', 1999, '523', '223'),
    'Haffner2003': ('10.1086/378850', 'astro-ph/0309117', 'Haffner', 2003, '149', '405'),
    'Haffner2009': ('10.1103/RevModPhys.81.969', '0901.0941', 'Haffner', 2009, '81', '969'),
    'Hanuschik2003': ('10.1051/0004-6361:20030885', None, 'Hanuschik', 2003, '407', '1157'),
    'Haslam1982': (None, None, 'Haslam', 1982, '47', '1'),
    'Intema2017': ('10.1051/0004-6361/201628536', None, 'Intema', 2017, '598', 'A78'),
    'Jordi2010': ('10.1051/0004-6361/201015441', None, 'Jordi', 2010, '523', 'A48'),
    'Katz2023': ('10.1051/0004-6361/202244220', None, 'Katz', 2023, '674', 'A5'),
    'Knies2024': ('10.1051/0004-6361/202348834', '2401.17289', 'Knies', 2024, '688', 'A90'),
    'Kopsacheili2020': ('10.1093/mnras/stz2594', None, 'Kopsacheili', 2020, '491', '889'),
    'Kurucz2005': (None, None, 'Kurucz', 2005, '8', '189'),
    'Luo2015': ('10.1088/1674-4527/15/8/002', None, 'Luo', 2015, '15', '1095'),
    'Madsen2006': ('10.1086/508441', 'astro-ph/0609558', 'Madsen', 2006, '652', '401'),
    'Marino2013': ('10.1051/0004-6361/201321956', '1307.5316', 'Marino', 2013, '559', 'A114'),
    'Martin2005': ('10.1086/426387', None, 'Martin', 2005, '619', 'L1'),
    'Mendoza1982': ('10.1093/mnras/198.1.127', None, 'Mendoza', 1982, '198', '127'),
    'Morrissey2007': ('10.1086/520512', None, 'Morrissey', 2007, '173', '682'),
    'Morton2000': ('10.1086/317349', None, 'Morton', 2000, '130', '403'),
    'Oke1983': ('10.1086/160817', None, 'Oke', 1983, '266', '713'),
    'Persic1996': ('10.1093/mnras/278.1.27', None, 'Persic', 1996, '281', '27'),
    'Pettini2004': ('10.1111/j.1365-2966.2004.07591.x', 'astro-ph/0401128', 'Pettini', 2004, '348', 'L59'),
    'Proxauf2014': ('10.1051/0004-6361/201322581', None, 'Proxauf', 2014, '561', 'A10'),
    'Ramsbottom1996': ('10.1006/adnd.1996.0009', None, 'Ramsbottom', 1996, '63', '57'),
    'Raymond1979': ('10.1086/190562', None, 'Raymond', 1979, '39', '1'),
    'Raymond1988': ('10.1086/165945', None, 'Raymond', 1988, '324', '869'),
    'Redfield2008': ('10.1086/524002', '0709.4480', 'Redfield', 2008, '673', '283'),
    'Reid2014': ('10.1088/0004-637X/783/2/130', '1401.5377', 'Reid', 2014, '783', '130'),
    'Reid2019': ('10.3847/1538-4357/ab4a11', '1910.03357', 'Reid', 2019, '885', '131'),
    'Ren2021': ('10.1088/1674-4527/21/3/51', '2008.06236', 'Ren', 2021, '21', '051'),
    'Reynolds1977': ('10.1086/154908', None, 'Reynolds', 1977, '211', '115'),
    'Rynkun2019': ('10.1051/0004-6361/201834931', None, 'Rynkun', 2019, '623', 'A155'),
    'Sabin2013': ('10.1093/mnras/stt160', None, 'Sabin', 2013, '431', '279'),
    'Schlegel1998': ('10.1086/305772', None, 'Schlegel', 1998, '500', '525'),
    'Schoenrich2010': ('10.1111/j.1365-2966.2010.16253.x', None, 'Schönrich', 2010, '403', '1829'),
    'SDSS2025': (None, '2507.07093', 'SDSS Collaboration', 2025, None, None),
    'Stark1989': ('10.1086/167334', None, 'Stark', 1989, '339', '763'),
    'Storey2000': ('10.1046/j.1365-8711.2000.03184.x', None, 'Storey', 2000, '312', '813'),
    'Tayal2010': ('10.1088/0067-0049/188/1/32', None, 'Tayal', 2010, '188', '32'),
    'Wenger2018': ('10.3847/1538-4357/aaaec8', '1802.04203', 'Wenger', 2018, '856', '52'),
    'Wu2022': ('10.1088/1674-4527/ac7387', '2110.04547', 'Wu', 2022, '22', '075015'),
    'Wu2026': ('10.3847/1538-3881/aeb11c', '2610.13329', 'Wu', 2026, '172', '298'),
    'Zhang2021': ('10.1088/1674-4527/21/11/280', '2108.08021', 'Zhang', 2021, '21', '280'),
    'ZhangLY2025': ('10.3847/1538-3881/ae0f94', '2510.08950', 'Zhang', 2025, '170', '309'),
    'Ziegenbalg2025': ('10.3847/2515-5172/adfec7', None, 'Ziegenbalg', 2025, '9', '227'),
}
DATACITE = {'NIST_ASD': '10.18434/T4W30F', 'kd_zenodo': '10.5281/zenodo.1166001'}
# works with no DOI or only an arXiv ID: search Crossref for a registered record / journal version
BIBSEARCH = {
    'Haslam1982': 'Haslam Salter Stoffel Wilson 1982 408 MHz all-sky continuum survey atlas of contour maps',
    'Kurucz2005': 'Kurucz 2005 New atlases for solar flux irradiance central intensity limb intensity',
    'DESI2025': 'Data Release 1 of the Dark Energy Spectroscopic Instrument',
    'SDSS2025': 'The Nineteenth Data Release of the Sloan Digital Sky Survey',
    'Fesen2026': 'Optical Discovery of New and Candidate Galactic Supernova Remnants Plus Optical Imaging of the '
                 'Monogem Supernova Remnant',
}


def fold(s):
    s = unicodedata.normalize('NFKD', s or '')
    return re.sub(r'[^a-z]', '', ''.join(c for c in s if not unicodedata.combining(c)).lower())


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40) as r:
                return r.read().decode('utf-8', 'replace')
        except Exception:                                       # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(2 + 3 * i)


def clean(t):
    return re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', t or '')).strip()


def year_of(m, field):
    dp = (m.get(field) or {}).get('date-parts') or [[None]]
    return dp[0][0]


def crossref_record(m):
    auth = m.get('author', [])
    names = []
    for a in auth:
        if 'family' in a:
            names.append((a.get('family', '') + ', ' + a.get('given', '')).strip(', '))
        else:
            names.append(a.get('name', ''))
    return dict(doi=m.get('DOI'), type=m.get('type'), title=clean((m.get('title') or [''])[0]),
                subtitle=clean((m.get('subtitle') or [''])[0]),
                journal=(m.get('container-title') or [''])[0],
                short_journal=(m.get('short-container-title') or [''])[0],
                volume=m.get('volume'), issue=m.get('issue'), page=m.get('page'),
                article_number=m.get('article-number'), n_authors=len(auth), authors=names[:6],
                year_issued=year_of(m, 'issued'), year_print=year_of(m, 'published-print'),
                year_online=year_of(m, 'published-online'),
                updates=[u.get('type') for u in m.get('update-to', [])],
                relation=sorted((m.get('relation') or {}).keys()))


def arxiv_records(ids):
    out = {}
    for chunk in [ids[i:i + 10] for i in range(0, len(ids), 10)]:
        xml = get('https://export.arxiv.org/api/query?id_list=' + ','.join(chunk) + '&max_results=20')
        for entry in xml.split('<entry>')[1:]:
            idm = re.search(r'<id>https?://arxiv.org/abs/(.+?)</id>', entry)
            if not idm:
                continue
            aid = re.sub(r'v\d+$', '', idm.group(1))
            out[aid] = dict(
                id_version=idm.group(1),
                title=clean(re.search(r'<title>(.*?)</title>', entry, re.S).group(1)),
                authors=[clean(a) for a in re.findall(r'<name>(.*?)</name>', entry)][:6],
                n_authors=len(re.findall(r'<name>', entry)),
                published=re.search(r'<published>(.*?)</published>', entry).group(1),
                updated=re.search(r'<updated>(.*?)</updated>', entry).group(1),
                journal_ref=clean((re.search(r'<arxiv:journal_ref[^>]*>(.*?)</arxiv:journal_ref>', entry, re.S)
                                   or [None, ''])[1]) or None,
                doi=clean((re.search(r'<arxiv:doi[^>]*>(.*?)</arxiv:doi>', entry, re.S) or [None, ''])[1]) or None,
                comment=clean((re.search(r'<arxiv:comment[^>]*>(.*?)</arxiv:comment>', entry, re.S)
                               or [None, ''])[1]) or None)
        time.sleep(3.5)
    return out


def compare(key, exp, rec):
    doi, _, fam, yr, vol, pg = exp
    issues = []
    first = (rec.get('authors') or [''])[0]
    if fold(fam) not in fold(first) and fold(first.split(',')[0]) not in fold(fam):
        issues.append(f'first author: printed "{fam}", record "{first}"')
    yrs = {rec.get('year_print'), rec.get('year_issued')} - {None}
    if yr not in yrs:
        issues.append(f'year: printed {yr}, record print/issued {sorted(yrs)} (online {rec.get("year_online")})')
    if vol and str(rec.get('volume')) != vol:
        issues.append(f'volume: printed {vol}, record {rec.get("volume")}')
    first_page = (rec.get('page') or '').split('-')[0] or (rec.get('article_number') or '')
    if pg and first_page.lstrip('0') != pg.lstrip('0') and (rec.get('article_number') or '').lstrip('0') != pg.lstrip('0'):
        issues.append(f'page/article: printed {pg}, record page {rec.get("page")} article {rec.get("article_number")}')
    return issues


def main():
    out = dict(_retrieved=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), crossref={}, arxiv={}, datacite={},
               bibsearch={}, comparison={})
    for key, exp in EXPECTED.items():
        doi = exp[0]
        if not doi:
            continue
        try:
            m = json.loads(get('https://api.crossref.org/works/' + urllib.parse.quote(doi, safe='/:;()')))['message']
            rec = crossref_record(m)
            out['crossref'][key] = rec
            out['comparison'][key] = compare(key, exp, rec)
        except Exception as e:                                  # noqa: BLE001
            out['crossref'][key] = dict(doi=doi, error=str(e))
            out['comparison'][key] = [f'crossref lookup failed: {e}']
        time.sleep(0.4)
    ids = [e[1] for e in EXPECTED.values() if e[1]]
    try:
        ax = arxiv_records(ids)
    except Exception as e:                                      # noqa: BLE001
        ax = {'_error': str(e)}
    for key, exp in EXPECTED.items():
        if exp[1]:
            out['arxiv'][key] = ax.get(exp[1], dict(error='not returned'))
    for key, doi in DATACITE.items():
        try:
            a = json.loads(get('https://api.datacite.org/dois/' + doi))['data']['attributes']
            out['datacite'][key] = dict(doi=a.get('doi'), title=(a.get('titles') or [{}])[0].get('title'),
                                        creators=[c.get('name') for c in a.get('creators', [])][:6],
                                        publisher=a.get('publisher'), year=a.get('publicationYear'),
                                        version=a.get('version'), url=a.get('url'),
                                        related=[(r.get('relationType'), r.get('relatedIdentifier'))
                                                 for r in a.get('relatedIdentifiers', [])][:8])
        except Exception as e:                                  # noqa: BLE001
            out['datacite'][key] = dict(doi=doi, error=str(e))
    for key, q in BIBSEARCH.items():
        try:
            items = json.loads(get('https://api.crossref.org/works?rows=3&query.bibliographic='
                                   + urllib.parse.quote(q)))['message']['items']
            out['bibsearch'][key] = [dict(doi=i.get('DOI'), title=clean((i.get('title') or [''])[0]),
                                          journal=(i.get('container-title') or [''])[0],
                                          first_author=((i.get('author') or [{}])[0].get('family')
                                                        or (i.get('author') or [{}])[0].get('name')),
                                          year=year_of(i, 'issued'), volume=i.get('volume'),
                                          page=i.get('page') or i.get('article-number')) for i in items]
        except Exception as e:                                  # noqa: BLE001
            out['bibsearch'][key] = dict(error=str(e))
        time.sleep(0.4)
    # journal versions of the arXiv-only works: title search restricted to the AAS (10.3847) prefix
    out['journal_version_search'] = {}
    for key, q in (('DESI2025', 'Data Release 1 of the Dark Energy Spectroscopic Instrument'),
                   ('SDSS2025', 'The Nineteenth Data Release of the Sloan Digital Sky Survey'),
                   ('Fesen2026', 'Optical Discovery of New and Candidate Galactic Supernova Remnants')):
        try:
            items = json.loads(get('https://api.crossref.org/works?rows=5&filter=prefix:10.3847,'
                                   'from-pub-date:2025-01-01&query.title=' + urllib.parse.quote(q)))['message']['items']
            out['journal_version_search'][key] = [
                dict(doi=i.get('DOI'), title=clean((i.get('title') or [''])[0]),
                     journal=(i.get('container-title') or [''])[0], year=year_of(i, 'issued'),
                     volume=i.get('volume'), page=i.get('page') or i.get('article-number')) for i in items]
        except Exception as e:                                  # noqa: BLE001
            out['journal_version_search'][key] = dict(error=str(e))
        time.sleep(0.4)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print('written', OUT)
    for k, v in out['comparison'].items():
        if v:
            print('MISMATCH', k, v)
    for k, v in out['arxiv'].items():
        print('ARXIV', k, v.get('id_version'), '|', v.get('journal_ref'), '|', v.get('doi'), '|',
              (v.get('title') or v.get('error', ''))[:90])


if __name__ == '__main__':
    main()
