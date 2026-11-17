"""Filter AstroBin search pages (data/raw/astrobin_search/<tag>_p*.json) to images whose plate-solved field covers
part of the Silk Nebula [OIII] ridge (NSNS ridge polynomial, Dec 12-27). Coverage = ridge points inside the circle of
0.75 x fieldRadius around the field centre AND inside the RA/Dec min-max box.
Output: data/derived/imaging_optical_astrobin_candidates_<tag>.csv
Usage: python imaging_optical_astrobin_filter.py <tag>"""
import sys, glob, json, csv
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from imaging_optical_common import RAW, DER
from imaging_optical_nsns import ridge_spine

tag = sys.argv[1]
sra, sdec = ridge_spine(step=0.1)
m = (sdec >= 12) & (sdec <= 27); sra, sdec = sra[m], sdec[m]
rows = []
seen = set()
for f in sorted(glob.glob(str(RAW / f'astrobin_search/{tag}_p*.json'))):
    for r in json.loads(open(f).read())['results']:
        if r['hash'] in seen or r.get('coordRaCenter') is None or r.get('fieldRadius') is None:
            continue
        seen.add(r['hash'])
        ra0, de0, fr = float(r['coordRaCenter']), float(r['coordDecCenter']), float(r['fieldRadius'])
        d = np.degrees(np.arccos(np.clip(np.sin(np.radians(de0)) * np.sin(np.radians(sdec)) +
                                         np.cos(np.radians(de0)) * np.cos(np.radians(sdec)) * np.cos(np.radians(sra - ra0)), -1, 1)))
        inbox = (sra >= float(r['coordRaMin'] or -1)) & (sra <= float(r['coordRaMax'] or 999)) & \
                (sdec >= float(r['coordDecMin'] or -99)) & (sdec <= float(r['coordDecMax'] or 99))
        cov = (d <= 0.75 * fr) & inbox
        if cov.sum() == 0:
            continue
        txt = (r.get('text') or '').lower()
        rows.append(dict(hash=r['hash'], url=f"https://www.astrobin.com/{r['hash']}/", title=r['title'], user=r.get('userDisplayName'),
                         published=r['published'], ra=ra0, dec=de0, field_radius=fr, pixel_scale=r.get('pixelScale'),
                         integration_h=round(float(r.get('integration') or 0) / 3600, 2),
                         ridge_dec_covered=f'{sdec[cov].min():.1f}-{sdec[cov].max():.1f}', ridge_len_deg=round(0.1 * cov.sum(), 1),
                         oiii_in_text=('oiii' in txt or 'o iii' in txt or '[o iii]' in txt), hd=r.get('hdThumbnail')))
rows.sort(key=lambda x: -x['ridge_len_deg'])
with open(DER / f'imaging_optical_astrobin_candidates_{tag}.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ['hash'])
    w.writeheader(); w.writerows(rows)
print(len(seen), 'solved images;', len(rows), 'cover the ridge')
for x in rows[:40]:
    print(x['hash'], x['published'][:10], f"{x['ra']:.1f} {x['dec']:+.1f} r={x['field_radius']:.1f}", x['ridge_dec_covered'], x['ridge_len_deg'], x['integration_h'], '|', x['title'][:60], '|', x['user'])
