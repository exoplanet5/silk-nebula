"""verify_desi_fetch.py -- (wave-2 DESI verification) extend DESI DR1 coverage near P3 / Dec 12-15.5 and add
control fibres > 5 deg from every strand, reusing the wave-1 resumable SPARCL fetcher logic.

The wave-1 fetcher (scripts/spectra_desi_fetch.py) writes into data/derived/spectra_desi_windows/, which the
wave-1 scripts glob; to leave wave-1 products untouched, this wrapper imports that module and redirects its
OUTDIR to data/derived/verify_desi_windows/ before calling its own `work()` (same +-40 A windows, same JSON-only
parsing, nothing unpickled).

Usage:
  python -I scripts/verify_desi_fetch.py select          -> data/derived/verify_desi_sel_extend.csv
  python -I scripts/verify_desi_fetch.py fetch [--workers 3] [--minutes 60]
Selection (GALAXY/QSO, ZWARN=0, not already fetched in wave 1):
  near : Dec 12-15.5 and d_any < 1.0 deg, plus all within 0.7 deg of P3 (116.45,+14.72), sorted by d_any;
  ctrl : d_any > 5 deg anywhere in the SPARCL box (RA 111-122, Dec 6-28), random order (seed 7), up to 1200;
  queue interleaves near/ctrl batches 1:1 so a stopped run still has both.
"""
import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
OUT = D / 'verify_desi_windows'
SEL = D / 'verify_desi_sel_extend.csv'
P3 = (116.45, 14.72)


def select():
    cat = list(csv.DictReader(open(D / 'spectra_desi_catalog.csv')))
    have = set(r['sparcl_id'] for r in csv.DictReader(open(D / 'spectra_desi_fiber_table.csv')))
    near, ctrl = [], []
    for r in cat:
        if r['spectype'] not in ('GALAXY', 'QSO') or int(r['zwarn']) != 0 or r['sparcl_id'] in have:
            continue
        ra, dec, d = float(r['ra']), float(r['dec']), float(r['d_any'])
        dp3 = np.hypot((ra - P3[0]) * np.cos(np.radians(P3[1])), dec - P3[1])
        if (12.0 <= dec < 15.5 and d < 1.0) or dp3 < 0.7:
            r['set'] = 'near'
            near.append(r)
        elif d > 5.0:
            r['set'] = 'ctrl'
            ctrl.append(r)
    near.sort(key=lambda r: float(r['d_any']))
    rng = np.random.default_rng(7)
    ctrl = [ctrl[i] for i in rng.permutation(len(ctrl))[:1200]]
    q, B = [], 25
    for k in range(0, max(len(near), len(ctrl)), B):
        q += near[k:k + B] + ctrl[k:k + B]
    with open(SEL, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(q[0].keys()))
        w.writeheader()
        w.writerows(q)
    print(f'near {len(near)}  ctrl {len(ctrl)}  queue {len(q)} -> {SEL}')


def fetch(workers, minutes):
    import spectra_desi_fetch as f
    OUT.mkdir(parents=True, exist_ok=True)
    f.OUTDIR = OUT                       # redirect output away from the wave-1 directory
    tag = 'ext'
    ids = [r['sparcl_id'] for r in csv.DictReader(open(SEL))]
    have = f.done_ids(tag)
    todo = [i for i in ids if i not in have]
    existing = [int(p.stem.split('_')[-1]) for p in OUT.glob(f'{tag}_*.npz')]
    b0 = max(existing) + 1 if existing else 0
    batches = [todo[i:i + 25] for i in range(0, len(todo), 25)]
    print(f'{len(ids)} selected, {len(have)} done, {len(batches)} batches to fetch', flush=True)
    t0, tot = time.time(), 0
    with ThreadPoolExecutor(workers) as ex:
        futs = [ex.submit(f.work, b0 + j, b, tag, ['sparcl_id', 'flux', 'ivar', 'model']) for j, b in
                enumerate(batches)]
        for fu in as_completed(futs):
            bno, n, nb, dt = fu.result()
            tot += n
            el = (time.time() - t0) / 60
            print(f'batch {bno}: {n} spectra {nb / 1e6:.1f} MB {dt:.0f}s | total {tot} in {el:.1f} min', flush=True)
            if el > minutes:
                print('time limit reached; cancelling remaining batches', flush=True)
                for g in futs:
                    g.cancel()
                break


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['select', 'fetch'])
    ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--minutes', type=float, default=60.0)
    a = ap.parse_args()
    select() if a.mode == 'select' else fetch(a.workers, a.minutes)
