"""verify_desi_tap.py -- (wave-2 DESI verification) GFA fibre-acceptance columns of the DESI DR1 exposures used
by the wave-1 sample, for the point-source -> surface-brightness conversion check.

Service: NOIRLab Data Lab TAP https://datalab.noirlab.edu/tap/sync (ADQL), table desi_dr1.exposure.
Output : data/raw/verify_desi_tap/exposure_gfa_box.csv  (+ query and date in exposure_gfa_box.query.txt)
Usage  : ~/.venvs/astro313/bin/python -I scripts/verify_desi_tap.py
"""
import csv
import time
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/verify_desi_tap'
TAP = 'https://datalab.noirlab.edu/tap/sync'

ids = sorted({int(r['tileid']) for r in csv.DictReader(open(ROOT / 'data/raw/datalab_tap/tiles_box.csv'))})
q = ('SELECT expid, tileid, night, mjd, exptime, efftime_spec, airmass, seeing_etc, seeing_gfa, airmass_gfa, '
     'transparency_gfa, fiber_fracflux_gfa, fiber_fracflux_elg_gfa, fiber_fracflux_bgs_gfa, fiberfac_gfa, '
     'fiberfac_elg_gfa, sky_mag_ab_gfa, sky_mag_g_spec, survey, program FROM desi_dr1.exposure '
     f'WHERE tileid IN ({",".join(map(str, ids))})')
r = requests.post(TAP, data=dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=q), timeout=600)
r.raise_for_status()
(OUT / 'exposure_gfa_box.csv').write_bytes(r.content)
(OUT / 'exposure_gfa_box.query.txt').write_text(f'{TAP}\n{time.strftime("%Y-%m-%dT%H:%M:%S")}\n{q}\n')
print(r.content[:600].decode(errors='replace'))
print(len(r.content), 'bytes')
