#!/usr/bin/env python
"""spectra_ism: download selected archival 1D high-resolution spectra of early-type stars near the
Silk Nebula (selected from data/derived/spectra_ism_archive_candidates.csv; highest-S/N product per star
covering Ca II K and/or Na I D).

ESO phase-3 products:  https://dataportal.eso.org/dataPortal/file/<dp_id>
SOPHIE s1d:            http://atlas.obs-hp.fr/sophie/sophie.cgi?c=i&a=mime:application/fits&o=sophie:[s1d,<seq>]
ELODIE s1d:            http://atlas.obs-hp.fr/elodie/fE.cgi?n=e500&c=i&z=s1d&a=mime:application/fits&o=elodie:<night>/<num>
Files go to data/raw/spectra_ism_spectra/; manifest data/raw/spectra_ism_spectra/manifest.csv
"""
import csv
import sys
import time

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, fetch  # noqa: E402

OUT = RAW / 'spectra_ism_spectra'
OUT.mkdir(parents=True, exist_ok=True)

ESO = [  # star, dp_id, note
    ('HD 61717', 'ADP.2024-03-09T06:44:09.481', 'UVES red 4726-6835'),
    ('HD 61717', 'ADP.2024-03-09T06:44:09.485', 'UVES blue 3282-4563'),
    ('HD 63432', 'ADP.2024-03-09T06:34:44.956', 'UVES red'),
    ('HD 63432', 'ADP.2024-03-09T06:34:44.960', 'UVES blue'),
    ('HD 63774', 'ADP.2024-03-09T07:02:16.839', 'UVES red'),
    ('HD 63774', 'ADP.2024-03-09T07:02:16.843', 'UVES blue'),
    ('HD 68332', 'ADP.2024-03-09T07:02:16.832', 'UVES red'),
    ('HD 68332', 'ADP.2024-03-09T07:02:16.836', 'UVES blue'),
    ('HD 53588', 'ADP.2024-03-09T06:47:54.823', 'UVES red'),
    ('HD 53588', 'ADP.2024-03-09T06:47:54.827', 'UVES blue'),
    ('HD 52124', 'ADP.2024-03-09T06:44:09.495', 'UVES red'),
    ('HD 52124', 'ADP.2024-03-09T06:44:09.499', 'UVES blue'),
    ('11 CMi', 'ADP.2020-09-14T18:52:16.793', 'UVES red 4727-6835'),
    ('11 CMi', 'ADP.2020-09-14T18:50:06.638', 'UVES blue 3105-3983'),
    ('BR CMi', 'ADP.2014-09-23T11:02:27.347', 'HARPS'),
    ('8 Cnc', 'ADP.2014-09-16T11:04:03.200', 'HARPS'),
    ('HD 67959', 'ADP.2014-09-26T16:53:39.010', 'HARPS'),
    ('HD 52479', 'ADP.2014-09-26T16:53:53.353', 'HARPS'),
    ('BD+17 1832', 'ADP.2021-01-22T01:14:15.071', 'HARPS'),
]
SOPHIE = [('HD 61038', 1894720), ('HD 64705', 1894780), ('HD 66552', 1894750), ('HD 55869', 1742770),
          ('HD 57930', 2347140), ('HD 52710', 1510420), ('64 Gem', 767240), ('HD 53588', 2072860),
          ('HD 67959', 1705680), ('HD 65900', 1705670), ('8 Cnc', 310410)]
ELODIE = [('HD 71557', '20000129/0020'), ('HD 59605', '20020328/0019'), ('HD 58144', '20020330/0021'),
          ('HD 54217', '20011128/0029'), ('HD 71262', '20001215/0021'), ('HD 71310', '20001216/0020'),
          ('HD 59473', '20011221/0015'), ('HD 56359', '20011219/0027'), ('HD 64934', '20041205/0032'),
          ('64 Gem', '20050130/0036')]

man = []
for star, dp, note in ESO:
    dest = OUT / f"eso_{star.replace(' ', '').replace('+', 'p')}_{dp}.fits"
    url = f'https://dataportal.eso.org/dataPortal/file/{dp}'
    ok, n, st = (True, dest.stat().st_size, 'cached') if dest.exists() else fetch(url, dest, timeout=300, min_bytes=10000)
    man.append(dict(star=star, archive='ESO', id=dp, note=note, url=url, file=dest.name, ok=ok, bytes=n))
    print(star, dp, ok, n)
    time.sleep(1)
for star, seq in SOPHIE:
    dest = OUT / f"sophie_{star.replace(' ', '')}_{seq}_s1d.fits"
    url = f'http://atlas.obs-hp.fr/sophie/sophie.cgi?c=i&a=mime:application/fits&o=sophie:[s1d,{seq}]'
    ok, n, st = (True, dest.stat().st_size, 'cached') if dest.exists() else fetch(url, dest, timeout=300, min_bytes=10000)
    man.append(dict(star=star, archive='SOPHIE', id=str(seq), note='SOPHIE HR s1d', url=url, file=dest.name, ok=ok, bytes=n))
    print(star, seq, ok, n)
    time.sleep(1)
for star, num in ELODIE:
    dest = OUT / f"elodie_{star.replace(' ', '')}_{num.replace('/', '_')}_s1d.fits"
    url = f'http://atlas.obs-hp.fr/elodie/fE.cgi?n=e500&c=i&z=s1d&a=mime:application/fits&o=elodie:{num}'
    ok, n, st = (True, dest.stat().st_size, 'cached') if dest.exists() else fetch(url, dest, timeout=300, min_bytes=10000)
    man.append(dict(star=star, archive='ELODIE', id=num, note='ELODIE s1d', url=url, file=dest.name, ok=ok, bytes=n))
    print(star, num, ok, n)
    time.sleep(1)
with open(OUT / 'manifest.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(man[0].keys()))
    w.writeheader()
    w.writerows(man)
print('total MB', sum(m['bytes'] for m in man) / 1e6)
