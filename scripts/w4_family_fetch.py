"""w4_family_fetch.py -- wave-4 `family`: fetch the DESI DR1 coadds in data/derived/w4_family_sel_queue.csv from
SPARCL, reusing the wave-1 resumable fetcher (scripts/spectra_desi_fetch.py: same +-40 A windows, JSON-only parsing,
nothing unpickled) with its output redirected to data/derived/w4_family_windows/ (tag 'fam').
Log: data/derived/w4_family_windows/fam_log.jsonl (bytes, seconds per batch).
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_fetch.py [--workers 3] [--minutes 40] [--queue CSV]
"""
import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/w4_family_windows'
SEL = ROOT / 'data/derived/w4_family_sel_queue.csv'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--minutes', type=float, default=40.0)
    ap.add_argument('--queue', default=str(SEL))
    a = ap.parse_args()
    import spectra_desi_fetch as f
    OUT.mkdir(parents=True, exist_ok=True)
    f.OUTDIR = OUT
    tag = 'fam'
    ids = [r['sparcl_id'] for r in csv.DictReader(open(a.queue))]
    have = f.done_ids(tag)
    todo = [i for i in ids if i not in have]
    existing = [int(p.stem.split('_')[-1]) for p in OUT.glob(f'{tag}_*.npz')]
    b0 = max(existing) + 1 if existing else 0
    batches = [todo[i:i + 25] for i in range(0, len(todo), 25)]
    print(f'{len(ids)} selected, {len(have)} done, {len(batches)} batches', flush=True)
    t0, tot = time.time(), 0
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(f.work, b0 + j, b, tag, ['sparcl_id', 'flux', 'ivar', 'model']) for j, b in
                enumerate(batches)]
        for fu in as_completed(futs):
            bno, n, nb, dt = fu.result()
            tot += n
            el = (time.time() - t0) / 60
            print(f'batch {bno}: {n} spectra {nb / 1e6:.1f} MB {dt:.0f}s | total {tot} in {el:.1f} min', flush=True)
            if el > a.minutes:
                print('time limit reached; cancelling remaining batches', flush=True)
                for g in futs:
                    g.cancel()
                break


if __name__ == '__main__':
    main()
