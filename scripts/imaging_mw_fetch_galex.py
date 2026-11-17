"""Download GALEX GR6/7 coadd background maps (-fd/-nd-skybg) for tiles centred in RA 112.5-120.5, Dec 10.5-22.5.

The pipeline skybg map is the source-free (sigma-clipped) diffuse background estimate, in counts/s/pixel
(1.5" pixels), i.e. already star-masked at the pipeline level. URLs from MAST CAOM TAP (galex_skybg.ecsv).
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from astropy.table import Table
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_fetch import fetch, RAW
t = Table.read(RAW / 'galex_mast' / 'galex_skybg.ecsv')
t = t[t['is_coadd'] & (t['s_ra'] > 112.5) & (t['s_ra'] < 120.5) & (t['s_dec'] > 10.5) & (t['s_dec'] < 22.5)]
print(len(t), 'files')

def get(row):
    url = row['access_url'].replace('http://', 'https://')
    return fetch(url, RAW / 'galex_mast' / 'skybg' / url.split('/')[-1], 'galex_mast',
                 'GALEX GR6/GR7 pipeline coadd skybg (MAST galex.stsci.edu)', timeout=300,
                 note=f"{row['survey']} {row['band']} exptime {row['t_exptime']:.0f}s")

with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(get, t))
print('done', sum(r is not None for r in res), 'of', len(res))
