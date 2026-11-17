"""dispute_references_manual.py -- record the checks of the references that have no Crossref DOI (tag: references).

Fetches, and stores as data:
  * the NIST ASD N I level page (version string, recommended citation, the 2D* J=5/2, 3/2 level energies);
  * HTTP status of the NRAO VLA OSS 2025B resolution page and the QHYCCD web note;
  * Crossref records of the A&A commentary that reprints the Haslam et al. (1982) citation, and of the
    data-credit papers that were considered but NOT added (Gaia mission, eROSITA DR1, astropy 2013/2018);
  * the kd commit used (project copy data/derived/verify_paper_raw/kd_commits.json).
Web-search confirmations (ADS records of Haslam 1982 and Kurucz 2005) are recorded with their URLs.
Responses are parsed as text/JSON only.

Output: data/derived/dispute_references_manual.json
Run:    ~/.venvs/astro313/bin/python -I scripts/dispute_references_manual.py
"""
import json
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/dispute_references_manual.json'
UA = {'User-Agent': 'generic-research-client/1.0 (reference check)'}
NIST_URL = ('https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=N+I&submit=Retrieve+Data&units=0'
            '&format=1&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&j_out=on&temp=')


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40) as r:
        return r.status, r.read().decode('utf-8', 'replace')


def crossref(doi):
    _, s = get('https://api.crossref.org/works/' + doi)
    m = json.loads(s)['message']
    a = m.get('author', [])
    return dict(doi=m['DOI'], first_author=(a[0].get('family') or a[0].get('name')) if a else None,
                second_author=a[1].get('family') if len(a) > 1 else None, n_authors=len(a),
                year=m['issued']['date-parts'][0][0], journal=(m.get('container-title') or [''])[0],
                volume=m.get('volume'), page=m.get('page') or m.get('article-number'),
                title=re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0])).strip())


def main():
    out = dict(_retrieved=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    st, s = get(NIST_URL)
    i = s.find('NIST ASD Team (')
    cite = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', s[i - 80:i + 420]).replace('&nbsp;', ' ')).strip()
    lv = re.findall(r'\|\s*(5/2|3/2)\s*\|\s*(1922\d\.\d+|1923\d\.\d+)', s)
    out['NIST_ASD'] = dict(url=NIST_URL, http=st, recommended_citation=cite,
                           version=(re.search(r'ver\. ([0-9.]+)', cite) or [None, None])[1],
                           N_I_2D_levels_cm1=lv,
                           section_values_cm1={'5/2': 19224.464, '3/2': 19233.177})
    for key, url in (('NRAO_OSS_2025B', 'https://science.nrao.edu/facilities/vla/docs/manuals/oss2025b/'
                                       'performance/resolution'),
                     ('QHYCCD_2024', 'https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-'
                                     'with-qhyccd-cameras/')):
        try:
            st, body = get(url)
            t = re.search(r'<title>(.*?)</title>', body, re.S)
            out[key] = dict(url=url, http=st, page_title=re.sub(r'\s+', ' ', t.group(1)).strip() if t else None)
        except Exception as e:                                  # noqa: BLE001
            out[key] = dict(url=url, error=str(e))
    qh = (ROOT / 'sources/web/discovery_qhyccd_pr.md').read_text().splitlines()[:4]
    out['QHYCCD_2024']['archived_copy'] = dict(path='sources/web/discovery_qhyccd_pr.md', header=qh)
    out['Haslam1982'] = dict(
        bibcode='1982A&AS...47....1H', doi=None,
        ads_url='https://ui.adsabs.harvard.edu/abs/1982A&AS...47....1H/abstract',
        confirmed_by=['ADS record (web search result title: "A 408-MHZ All-Sky Continuum Survey. II. The Atlas of '
                      'Contour Maps", C. G. T. Haslam, 1982)',
                      'A&A commentary whose title reprints the citation "Haslam C. G. T., Salter C. J., Stoffel H., '
                      'and Wilson W. E., 1982, A&AS, 47, 1"'],
        commentary=crossref('10.1051/0004-6361/200912161'))
    out['Kurucz2005'] = dict(
        bibcode='2005MSAIS...8..189K', doi=None,
        ads_url='https://ui.adsabs.harvard.edu/abs/2005MSAIS...8..189K/abstract',
        confirmed_by=['ADS record (web search result: "New atlases for solar flux, irradiance, central intensity, '
                      'and ...", R. L. Kurucz, 2005)', 'kurucz.harvard.edu/sun.html lists FLUXATLAS2005'])
    commits = json.loads((ROOT / 'data/derived/verify_paper_raw/kd_commits.json').read_text())
    out['kd'] = dict(commit=commits[0]['sha'], date=commits[0]['commit']['committer']['date'],
                     readme_request='If you use this work, please reference https://zenodo.org/record/1166001',
                     zenodo_doi='10.5281/zenodo.1166001 (DataCite: "tvwenger/kd v1.0", 2018)')
    out['considered_not_added'] = {}
    for key, doi in (('Gaia2016_mission', '10.1051/0004-6361/201629272'),
                     ('Merloni2024_eROSITA_DR1', '10.1051/0004-6361/202347165'),
                     ('Astropy2013', '10.1051/0004-6361/201322068'),
                     ('Astropy2018', '10.3847/1538-3881/aabc4f')):
        try:
            out['considered_not_added'][key] = crossref(doi)
        except Exception as e:                                  # noqa: BLE001
            out['considered_not_added'][key] = dict(doi=doi, error=str(e))
        time.sleep(0.4)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps({k: (v if k != 'considered_not_added' else list(v)) for k, v in out.items()},
                     ensure_ascii=False)[:2500])


if __name__ == '__main__':
    main()
