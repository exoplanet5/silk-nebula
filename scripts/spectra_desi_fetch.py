"""spectra_desi_fetch.py -- stream DESI DR1 spectra from SPARCL and keep only the +-40 A line windows.

Usage:  python -I scripts/spectra_desi_fetch.py <selection.csv> <tag> [--workers 3] [--batch 25] [--sigma]
  selection.csv: needs a sparcl_id column (rows from data/raw/sparcl_desi_dr1/find_box.csv).
  Output: data/derived/spectra_desi_windows/<tag>_<batchno>.npz with arrays
          sparcl_id (str), flux, ivar, model (float32, [n, len(WIN_IDX)]), and wave_sigma if --sigma.
  Log:    data/derived/spectra_desi_windows/<tag>_log.jsonl (bytes transferred, seconds per batch).
Resumable: ids already present in existing <tag>_*.npz files are skipped.
Downloaded JSON is parsed with json only; nothing is unpickled or executed.
"""
import argparse
import csv
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_sparcl import retrieve  # noqa: E402
from spectra_desi_lines import WIN_IDX, NPIX  # noqa: E402

OUTDIR = Path('/Users/mickey/astronomy/silk-nebula/data/derived/spectra_desi_windows')
OUTDIR.mkdir(parents=True, exist_ok=True)
lock = threading.Lock()


def done_ids(tag):
    s = set()
    for p in OUTDIR.glob(f'{tag}_*.npz'):
        with np.load(p, allow_pickle=False) as d:
            s.update(d['sparcl_id'].tolist())
    return s


def work(bno, ids, tag, include, tries=6):
    sess = requests.Session()
    for k in range(tries):
        try:
            recs, nbytes, dt = retrieve(ids, include=include, session=sess)
            break
        except Exception as e:
            msg = (type(e).__name__ + ": " + str(e))[-400:]
            with lock:
                print(f'batch {bno} try {k} failed: {msg}', flush=True)
            time.sleep(min(30 * (k + 1), 150))
    else:
        return bno, 0, 0, 0.0
    out = {'sparcl_id': [], 'flux': [], 'ivar': [], 'model': [], 'wave_sigma': []}
    for r in recs:
        f = np.asarray(r['flux'], np.float32)
        if f.size != NPIX:
            continue
        out['sparcl_id'].append(r['sparcl_id'])
        out['flux'].append(f[WIN_IDX])
        out['ivar'].append(np.asarray(r['ivar'], np.float32)[WIN_IDX])
        out['model'].append(np.asarray(r['model'], np.float32)[WIN_IDX])
        if 'wave_sigma' in r:
            out['wave_sigma'].append(np.asarray(r['wave_sigma'], np.float32)[WIN_IDX])
    arrs = {k: np.array(v) for k, v in out.items() if len(v)}
    arrs['sparcl_id'] = np.array(out['sparcl_id'], dtype='U36')
    np.savez_compressed(OUTDIR / f'{tag}_{bno:05d}.npz', **arrs)
    with lock:
        with open(OUTDIR / f'{tag}_log.jsonl', 'a') as fh:
            fh.write(json.dumps(dict(batch=bno, n_req=len(ids), n_got=len(out['sparcl_id']), bytes=nbytes,
                                     sec=round(dt, 1), t=time.strftime('%Y-%m-%dT%H:%M:%S'))) + '\n')
    return bno, len(out['sparcl_id']), nbytes, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('selection')
    ap.add_argument('tag')
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--batch', type=int, default=25)
    ap.add_argument('--sigma', action='store_true')
    ap.add_argument('--max-gb', type=float, default=12.0)
    a = ap.parse_args()
    ids = [r['sparcl_id'] for r in csv.DictReader(open(a.selection))]
    have = done_ids(a.tag)
    todo = [i for i in ids if i not in have]
    existing = [int(p.stem.split('_')[-1]) for p in OUTDIR.glob(f'{a.tag}_*.npz')]
    b0 = max(existing) + 1 if existing else 0
    batches = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
    include = ['sparcl_id', 'flux', 'ivar', 'model'] + (['wave_sigma'] if a.sigma else [])
    print(f'{len(ids)} selected, {len(have)} done, {len(todo)} to fetch in {len(batches)} batches', flush=True)
    tot_b, tot_n, t0 = 0, 0, time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(work, b0 + j, b, a.tag, include) for j, b in enumerate(batches)]
        for f in as_completed(futs):
            bno, n, nb, dt = f.result()
            tot_b += nb
            tot_n += n
            el = time.time() - t0
            print(f'batch {bno}: {n} spectra {nb/1e6:.1f} MB {dt:.0f}s | total {tot_n} spectra '
                  f'{tot_b/1e9:.2f} GB {el/60:.1f} min ({tot_n/max(el,1):.2f}/s)', flush=True)
            if tot_b > a.max_gb * 1e9:
                print('transfer cap reached; cancelling remaining batches', flush=True)
                for g in futs:
                    g.cancel()
                break


if __name__ == '__main__':
    main()
