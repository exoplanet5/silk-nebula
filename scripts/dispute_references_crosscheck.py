"""dispute_references_crosscheck.py -- match in-text citations in the six Comment sections to the verified list.

Reads notes/dispute/{distance,kinematic,oiii,halpha,brightness,other}.md. For each section it extracts
  (a) author-year citations from the body (text before '**References used in this section.**'), and
  (b) the first author + year of each entry in that section's reference list,
then matches both against MASTER (the verified entries of notes/dispute/references.md; DOIs/arXiv IDs checked by
scripts/dispute_references_verify.py -> data/derived/dispute_references_records.json).
Reports: citations with no verified entry, reference-list entries never cited in the body, and informal
shorthand (e.g. 'Reid+14').

Output: data/derived/dispute_references_crosscheck.json
Run:    ~/.venvs/astro313/bin/python -I scripts/dispute_references_crosscheck.py
"""
import json
import re
import unicodedata
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
SECTIONS = ['distance', 'kinematic', 'oiii', 'halpha', 'brightness', 'other']
OUT = ROOT / 'data/derived/dispute_references_crosscheck.json'

# key -> (first-author surname as cited, year)
MASTER = {
    'Abdurrouf2022': ("Abdurro'uf", 2022), 'AllendePrieto2001': ('Allende Prieto', 2001),
    'Almeida2023': ('Almeida', 2023), 'Anderson2012': ('Anderson', 2012), 'Asplund2004': ('Asplund', 2004),
    'Astropy2022': ('Astropy Collaboration', 2022), 'Blanton2007': ('Blanton', 2007), 'Boumis2022': ('Boumis', 2022),
    'Bovy2017': ('Bovy', 2017), 'Cioffi1988': ('Cioffi', 1988), 'Condon1998': ('Condon', 1998),
    'Cox1985': ('Cox', 1985), 'Cui2012': ('Cui', 2012), 'DESI2026': ('DESI Collaboration', 2026),
    'Dopita1979': ('Dopita', 1979), 'Dopita1996': ('Dopita', 1996), 'Dopita2016': ('Dopita', 2016),
    'Fesen2024': ('Fesen', 2024), 'Fesen2026': ('Fesen', 2026), 'Gaensler2008': ('Gaensler', 2008),
    'Gaia2023': ('Gaia Collaboration', 2023), 'Guy2023': ('Guy', 2023), 'Haffner1999': ('Haffner', 1999),
    'Haffner2003': ('Haffner', 2003), 'Haffner2009': ('Haffner', 2009), 'Hanuschik2003': ('Hanuschik', 2003),
    'Haslam1982': ('Haslam', 1982), 'Intema2017': ('Intema', 2017), 'Jordi2010': ('Jordi', 2010),
    'Katz2023': ('Katz', 2023), 'Knies2024': ('Knies', 2024), 'Kopsacheili2020': ('Kopsacheili', 2020),
    'Kramida2024': ('Kramida', 2024), 'Kurucz2005': ('Kurucz', 2005), 'Luo2015': ('Luo', 2015),
    'Madsen2006': ('Madsen', 2006), 'Marino2013': ('Marino', 2013), 'Martin2005': ('Martin', 2005),
    'Mendoza1982': ('Mendoza', 1982), 'Morrissey2007': ('Morrissey', 2007), 'Morton2000': ('Morton', 2000),
    'NRAO2025': ('NRAO', 2025), 'Oke1983': ('Oke', 1983), 'Persic1996': ('Persic', 1996),
    'Pettini2004': ('Pettini', 2004), 'Proxauf2014': ('Proxauf', 2014), 'QHYCCD2024': ('QHYCCD', 2024),
    'Ramsbottom1996': ('Ramsbottom', 1996), 'Raymond1979': ('Raymond', 1979), 'Raymond1988': ('Raymond', 1988),
    'Redfield2008': ('Redfield', 2008), 'Reid2014': ('Reid', 2014), 'Reid2019': ('Reid', 2019),
    'Ren2021': ('Ren', 2021), 'Reynolds1977': ('Reynolds', 1977), 'Rynkun2019': ('Rynkun', 2019),
    'Sabin2013': ('Sabin', 2013), 'Schlegel1998': ('Schlegel', 1998), 'Schoenrich2010': ('Schönrich', 2010),
    'SDSS2026': ('SDSS Collaboration', 2026), 'Stark1989': ('Stark', 1989), 'Storey2000': ('Storey', 2000),
    'Tayal2010': ('Tayal', 2010), 'Wenger2018': ('Wenger', 2018), 'Wu2022': ('Wu', 2022), 'Wu2026': ('Wu', 2026),
    'Zhang2021': ('Zhang', 2021), 'ZhangLY2025': ('Zhang', 2025), 'Ziegenbalg2025': ('Ziegenbalg', 2025),
}
# names that appear only inside verification notes, not as cited works
NOTE_ONLY = {('sato', 2008), ('wielebinski', 2009)}

NAME = r"[A-Z][A-Za-z'’À-ſ\-]+(?: (?:Prieto|Collaboration))?"
CIT = re.compile(rf"(?<![\w.])({NAME})(?: et al\.|,? (?:&|and) {NAME})?,? \(?((?:19|20)\d\d)\b")
SHORT = re.compile(r"\b([A-Z][a-z]+)\+(\d\d)\b")
REFLINE = re.compile(r"^- (.+?) \(((?:19|20)\d\d)\)", re.M)
STOP = {'Table', 'Figure', 'Fig', 'Sect', 'Section', 'Eq', 'Dec', 'DR', 'Commit', 'Crossref', 'January', 'February',
        'November', 'Note', 'Sun', 'Moon', 'December', 'September', 'October'}


def fold(s):
    s = unicodedata.normalize('NFKD', s)
    return re.sub(r"[^a-z]", '', ''.join(c for c in s if not unicodedata.combining(c)).lower())


LOOKUP = {}
for k, (n, y) in MASTER.items():
    LOOKUP.setdefault((fold(n), y), []).append(k)
    LOOKUP.setdefault((fold(n.split()[0]), y), []).append(k)


def match(name, year):
    return sorted(set(LOOKUP.get((fold(name), year), []) + LOOKUP.get((fold(name.split()[0]), year), [])))


def main():
    out = {}
    for sec in SECTIONS:
        txt = (ROOT / f'notes/dispute/{sec}.md').read_text()
        body, _, refs = txt.partition('**References used in this section.**')
        body_hits, unmatched = {}, []
        for m in CIT.finditer(body):
            name, year = m.group(1), int(m.group(2))
            if name in STOP or name.isupper() and name not in ('NRAO', 'QHYCCD', 'DESI', 'SDSS'):
                continue
            keys = match(name, year)
            if keys:
                for k in keys:
                    body_hits[k] = body_hits.get(k, 0) + 1
            elif (fold(name), year) not in NOTE_ONLY:
                unmatched.append(f'{name} {year}')
        ref_keys, ref_unmatched = [], []
        for m in REFLINE.finditer(refs):
            first = re.split(r",| \(|&", m.group(1))[0].strip()
            keys = match(first, int(m.group(2)))
            (ref_keys.extend(keys) if keys else ref_unmatched.append(f'{first} {m.group(2)}'))
        # an author with several works in one year range is ambiguous only if both are in MASTER for that year
        cited_only_in_list = sorted(set(ref_keys) - set(body_hits))
        cited_not_in_list = sorted(set(body_hits) - set(ref_keys))
        out[sec] = dict(body_citations=dict(sorted(body_hits.items())), body_unmatched=sorted(set(unmatched)),
                        reflist_keys=sorted(set(ref_keys)), reflist_unmatched=ref_unmatched,
                        in_reflist_not_cited_in_body=cited_only_in_list,
                        cited_in_body_not_in_reflist=cited_not_in_list,
                        informal_shorthand=sorted({f'{a}+{b}' for a, b in SHORT.findall(body)}))
    allcited = sorted({k for v in out.values() for k in list(v['body_citations']) + v['reflist_keys']})
    out['_summary'] = dict(n_master=len(MASTER), n_used=len(allcited),
                           master_not_used=sorted(set(MASTER) - set(allcited)))
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for sec in SECTIONS:
        v = out[sec]
        print(f"== {sec}: body {len(v['body_citations'])} works; list {len(v['reflist_keys'])}")
        for f in ('body_unmatched', 'reflist_unmatched', 'in_reflist_not_cited_in_body',
                  'cited_in_body_not_in_reflist', 'informal_shorthand'):
            if v[f]:
                print('  ', f, v[f])
    print(out['_summary'])


if __name__ == '__main__':
    main()
