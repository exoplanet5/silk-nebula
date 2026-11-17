"""Query MAST CAOM (TAP) for GALEX tile background maps (-fd/-nd-skybg) covering the Silk field.

MAST JSON API (mast.stsci.edu/api/v0/invoke) is blocked by the proxy; TAP (mast.stsci.edu/vo-tap) works.
Sibling products (-int, -rrhr) sit in the same directory with the same stem.
Output: data/raw/galex_mast/galex_skybg.ecsv
"""
import collections
import numpy as np
import pyvo
from astropy.table import vstack, unique
tap = pyvo.dal.TAPService('https://mast.stsci.edu/vo-tap/api/v0.1/caom')
parts = []
for d0 in range(5, 28, 3):
    q = f"""SELECT obs_id, target_name, s_ra, s_dec, t_exptime, em_min, access_url, access_estsize
    FROM ivoa.obscore WHERE obs_collection='GALEX' AND s_ra BETWEEN 109.0 AND 124.0
    AND s_dec >= {d0} AND s_dec < {d0 + 3} AND access_url LIKE '%d-skybg.fits.gz'"""
    r = tap.search(q, maxrec=20000).to_table()
    print(d0, len(r))
    parts.append(r)
t = unique(vstack(parts), keys='access_url')
t['band'] = ['FUV' if '-fd-' in u else 'NUV' for u in t['access_url']]
t['survey'] = [str(x).split('_')[0] for x in t['target_name']]
t['tile'] = [u.split('/')[-1].replace('-fd-skybg.fits.gz', '').replace('-nd-skybg.fits.gz', '') for u in t['access_url']]
t['is_coadd'] = ['/01-main/' in u for u in t['access_url']]
t.write('/Users/mickey/astronomy/silk-nebula/data/raw/galex_mast/galex_skybg.ecsv', overwrite=True, format='ascii.ecsv')
print(len(t), collections.Counter(zip(t['band'], t['survey'])).most_common(12))
print(collections.Counter(t['is_coadd']))
print(t['access_url'][0])
print(np.percentile(t['t_exptime'][t['survey'] == 'AIS'], [10, 50, 90]))
