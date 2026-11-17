"""Parallel downloader for the SDSS frame files listed in <dir>/urls.txt (verify_precovery).
Uses urllib (system proxy settings from the environment); writes <dir>/download_log.jsonl with url, status,
bytes and UTC time per file. Files already present with nonzero size are skipped.
Usage: python -I verify_precovery_download.py <dir> [nthreads]
"""
import sys, json, time, urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

d = Path(sys.argv[1]); nt = int(sys.argv[2]) if len(sys.argv) > 2 else 8
urls = [u.strip() for u in (d / 'urls.txt').read_text().splitlines() if u.strip()]
log = open(d / 'download_log.jsonl', 'a')


def get(u):
    out = d / u.rsplit('/', 1)[1]
    if out.exists() and out.stat().st_size > 0:
        return dict(url=u, status='exists', bytes=out.stat().st_size)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(u, timeout=300) as r:
                data = r.read()
            tmp = out.with_suffix('.part')
            tmp.write_bytes(data); tmp.rename(out)
            return dict(url=u, status=200, bytes=len(data), utc=datetime.now(timezone.utc).isoformat())
        except Exception as e:
            err = str(e); time.sleep(3)
    return dict(url=u, status='error', error=err, utc=datetime.now(timezone.utc).isoformat())


t0 = time.time()
with ThreadPoolExecutor(nt) as ex:
    for i, r in enumerate(ex.map(get, urls)):
        log.write(json.dumps(r) + '\n'); log.flush()
        if i % 25 == 0:
            print(i, len(urls), r.get('status'), f'{time.time() - t0:.0f}s', flush=True)
print('done', time.time() - t0)
