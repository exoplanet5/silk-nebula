"""SDSS DR18 imaging for the Silk Nebula region (SkyServer SQL + SAS frame files).

Step 1 (fields): list photometric fields in RA 109-124, Dec 6-30 -> data/raw/sdss_dr18/fields_box.csv
Step 2 (frames): download g and r 'frame' files (calibrated, sky-subtracted; sky model in HDU2) for runs/camcols
                 crossing the filament -> data/raw/sdss_dr18/frames/
SkyServer: https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch
SAS:       https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/<run>/<camcol>/frame-<b>-<run6>-<camcol>-<field4>.fits.bz2
Usage: python imaging_optical_fetch_sdss.py fields | frames
"""
import sys, csv, io, json, time, datetime as dt
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/sdss_dr18'
SQL = 'https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch'


def fields():
    OUT.mkdir(parents=True, exist_ok=True)
    q = ('SELECT run, camcol, field, rerun, ra, dec, raMin, raMax, decMin, decMax, score, quality, '
         'sky_g, sky_r, mjd_g FROM Field WHERE ra BETWEEN 109 AND 124 AND dec BETWEEN 6 AND 30')
    r = requests.get(SQL, params=dict(cmd=q, format='csv'), timeout=300)
    (OUT / 'fields_box.csv').write_text(r.text)
    (OUT / 'fields_box.json').write_text(json.dumps(dict(url=r.url, accessed=dt.datetime.now(dt.timezone.utc).isoformat()), indent=1))
    lines = [l for l in r.text.splitlines() if l and l[0].isdigit()]
    print('fields:', len(lines))
    print(r.text[:600])


if __name__ == '__main__':
    {'fields': fields}[sys.argv[1]]()
