#!/usr/bin/env python
"""spectra_ism: fetch interstellar-absorption, RM and pulsar catalogues from VizieR.

Mirror: https://vizier.cfa.harvard.edu/viz-bin/asu-tsv (vizier.cds.unistra.fr asu-tsv returned the CDS
"Page unavailable" 404 page on 2026-11-16). Tables are small; full tables are fetched (-out.max=unlimited)
or cones where the table is large. Output: data/raw/spectra_ism_vizier/<tag>.tsv
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_fetch_vizier.py [tag ...]
"""
import sys
import time

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, fetch  # noqa: E402

VIZ = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
OUT = RAW / 'spectra_ism_vizier'

# (tag, source, extra parameters); cones are centred on the filament mid-point (116.0, +20.0)
CONE = {'-c': '116.0 +20.0', '-c.rd': '20'}
JOBS = [
    ('welsh2010_naicaii', 'J/A+A/510/A54', {}),
    ('redfield_linsky2008_lism', 'J/ApJ/673/283', {}),
    ('redfield_linsky2002_caii', 'J/ApJS/139/439', {}),
    ('gudennavar2012_ism_columns', 'J/ApJS/199/8', {}),
    ('taylor2009_nvss_rm_cone', 'J/ApJ/702/1230', CONE),
    ('megier2009_caii', 'J/A+A/507/833', {}),
    ('hunter2006_caii_nai', 'J/MNRAS/367/1478', {}),
    ('sfeir1999_nai', 'J/A+A/346/785', {}),
    ('zasowski2015_apogee_dib_cone', 'J/ApJ/798/35', CONE),
    ('psrcat_bpsr_cone', 'B/psr', {'-c': '116.0 +20.0', '-c.rd': '35'}),
]


def run(tags=None):
    OUT.mkdir(parents=True, exist_ok=True)
    for tag, src, extra in JOBS:
        if tags and tag not in tags:
            continue
        params = {'-source': src, '-out.all': '', '-out.max': 'unlimited', '-oc.form': 'd'}
        params.update(extra)
        dest = OUT / f'{tag}.tsv'
        for attempt in range(4):
            ok, n, st = fetch(VIZ, dest, params=params, timeout=180, retries=1)
            txt = dest.read_text(errors='replace') if ok else ''
            bad = ('Error=' in txt and 'QUERY_STATUS=ERROR' in txt) or ('#RESOURCE=\n#Name:\n' in txt and '#Table' not in txt)
            if ok and not bad:
                nrow = sum(1 for ln in txt.splitlines() if ln and not ln.startswith('#'))
                print(f'{tag:32s} {src:18s} OK {n} B, ~{nrow} lines')
                break
            print(f'{tag:32s} {src:18s} attempt {attempt} failed (status {st}, bad={bad})')
            time.sleep(15 + 10 * attempt)
        time.sleep(3)


if __name__ == '__main__':
    run(sys.argv[1:] or None)
