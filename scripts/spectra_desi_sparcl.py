"""spectra_desi_sparcl.py -- minimal REST client for NOIRLab SPARCL (no sparclclient, JSON only, never pickle).

Endpoints were taken from the sparclclient source (github.com/astro-datalab/sparclclient, client.py:
find() -> POST {api}/find/?limit=&offset=&sort= with body {"outfields":[...], "search":[[field, v...]]};
retrieve() -> POST {api}/spectras/?include=a,b&format=json&dataset_list=DESI-DR1 with body [sparcl_id,...])
and the server's OpenAPI schema https://astrosparcl.datalab.noirlab.edu/api/docs/?format=openapi
(SPARCL API version 14.2, basePath /api; the client's older /sparc/ base path returns 404).
Accessed 2026-11-16.

Functions
  find_box(ra0, ra1, dec0, dec1, out_csv)  -> CSV of sparcl_id, targetid, survey, program, ra, dec, spectype, redshift
  retrieve(ids, include)                  -> (list of record dicts, compressed bytes, seconds)
"""
import csv
import gzip
import json
import time
import requests

API = 'https://astrosparcl.datalab.noirlab.edu/api'
DR = 'DESI-DR1'
FIND_FIELDS = ['sparcl_id', 'targetid', 'survey', 'program', 'ra', 'dec', 'spectype', 'redshift', 'specprimary']


def find_box(ra0, ra1, dec0, dec1, out_csv, page=20000, timeout=600):
    rows, offset = [], 0
    while True:
        body = {'outfields': FIND_FIELDS,
                'search': [['data_release', DR], ['ra', str(ra0), str(ra1)], ['dec', str(dec0), str(dec1)]]}
        for k in range(5):
            try:
                r = requests.post(f'{API}/find/?limit={page}&offset={offset}&sort=sparcl_id', json=body,
                                  timeout=timeout)
                r.raise_for_status()
                recs = r.json()[1:]
                break
            except Exception as e:
                print('find retry', k, e, flush=True)
                time.sleep(15 * (k + 1))
        else:
            raise RuntimeError('find failed')
        rows += recs
        print(f'find offset {offset}: +{len(recs)} (total {len(rows)})', flush=True)
        if len(recs) < page:
            break
        offset += page
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIND_FIELDS, extrasaction='ignore')
        w.writeheader()
        for rec in rows:
            w.writerow(rec)
    return len(rows)


def retrieve(ids, include=('sparcl_id', 'flux', 'ivar', 'model'), timeout=(20, 420), session=None):
    """POST ids to /spectras/ with format=json. Returns (records, n_compressed_bytes, seconds).

    The response is gzip-compressed JSON; it is parsed with json only (no pickle)."""
    s = session or requests
    url = f'{API}/spectras/?include={",".join(include)}&format=json&dataset_list={DR}'
    t0 = time.time()
    r = s.post(url, json=list(ids), headers={'Accept-Encoding': 'gzip', 'Content-Type': 'application/json'},
               timeout=timeout, stream=True)
    r.raise_for_status()
    raw = r.raw.read(decode_content=False)
    nbytes = len(raw)
    if r.headers.get('Content-Encoding', '').lower() == 'gzip':
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    meta = data[0]
    if not meta.get('status', {}).get('success', False):
        raise RuntimeError(f'SPARCL status: {meta.get("status")}')
    return data[1:], nbytes, time.time() - t0


if __name__ == '__main__':
    import sys
    out = sys.argv[1]
    ra0, ra1, dec0, dec1 = map(float, sys.argv[2:6])
    print(find_box(ra0, ra1, dec0, dec1, out))
