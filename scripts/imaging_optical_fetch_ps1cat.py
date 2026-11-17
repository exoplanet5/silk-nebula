"""PS1 DR2 mean-object photometry (MAST catalogs API) of faint stars near each DSS tile centre.

Used to calibrate DSS photographic 'density volume' against stellar flux (g for blue plates, r for red plates).
API: https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv
Output: data/raw/ps1_dr2/catalog/ps1cat_<tile>.csv
"""
import sys, json, time, datetime as dt
from pathlib import Path
import requests
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_optical_fetch_dss import tiles

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/ps1_dr2/catalog'
API = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log = []
    for dec, col, ra in tiles():
        name = f'd{dec:02d}_{col}'
        out = OUT / f'ps1cat_{name}.csv'
        if out.exists() and out.stat().st_size > 1000:
            continue
        p = {'ra': f'{ra:.4f}', 'dec': f'{dec:.4f}', 'radius': 0.3, 'nDetections.gte': 6,
             'gMeanPSFMag.gte': 17.0, 'gMeanPSFMag.lte': 21.0,
             'columns': '[objID,raMean,decMean,gMeanPSFMag,rMeanPSFMag,gMeanKronMag,rMeanKronMag,nDetections]',
             'pagesize': 20000}
        for k in range(3):
            try:
                r = requests.get(API, params=p, timeout=300)
                if r.status_code == 200 and 'raMean' in r.text[:500]:
                    out.write_text(r.text); break
                print('bad', name, r.status_code, r.text[:200])
            except Exception as e:
                print('err', name, e)
            time.sleep(5 * (k + 1))
        log.append(dict(tile=name, url=r.url, accessed=dt.datetime.now(dt.timezone.utc).isoformat(),
                        nrows=r.text.count('\n') - 1))
        print(name, log[-1]['nrows'], flush=True)
    (OUT / 'query_log.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    main()
