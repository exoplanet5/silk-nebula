"""Gaia DR3 stars (G < 17.5) along the Silk Nebula ridge for star masking of DSS/PS1/SDSS images.

Source: NOIRLab Astro Data Lab TAP, table gaia_dr3.gaia_source (https://datalab.noirlab.edu/tap).
Region: Dec 9.5-28.5 in 1-deg Dec slices, RA = ridge(dec) -/+ 2.6 deg/cos(dec).
Output: data/raw/gaia_dr3_datalab/gaia_g175_dec<d>.csv + query log JSON
"""
import json, time, datetime as dt
from pathlib import Path
import numpy as np
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/gaia_dr3_datalab'
TAP = 'https://datalab.noirlab.edu/tap/sync'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rid = json.loads((ROOT / 'data/derived/imaging_optical_nsns_ridge.json').read_text())
    log = []
    for d0 in range(9, 29):
        out = OUT / f'gaia_g175_dec{d0:02d}.csv'
        if out.exists() and out.stat().st_size > 1000:
            continue
        dmid = d0 + 0.5
        r0 = float(np.polyval(rid['poly_ra_of_dec'], np.clip(dmid, *rid['dec_range'])))
        hw = 2.6 / np.cos(np.radians(dmid + 0.5))
        q = (f'SELECT ra, dec, phot_g_mean_mag FROM gaia_dr3.gaia_source WHERE dec >= {d0} AND dec < {d0 + 1} '
             f'AND ra BETWEEN {r0 - hw:.3f} AND {r0 + hw:.3f} AND phot_g_mean_mag < 17.5')
        for k in range(3):
            try:
                r = requests.get(TAP, params=dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=q), timeout=600)
                if r.status_code == 200 and r.text.startswith('ra'):
                    out.write_text(r.text); break
                print('bad', d0, r.status_code, r.text[:300])
            except Exception as e:
                print('err', d0, e)
            time.sleep(10 * (k + 1))
        n = r.text.count('\n') - 1
        log.append(dict(query=q, url=TAP, accessed=dt.datetime.now(dt.timezone.utc).isoformat(), nrows=n))
        print(d0, n, flush=True)
    (OUT / 'query_log.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    main()
