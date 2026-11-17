"""w4_family_select.py -- wave-4 `family`: choose DESI DR1 coadds to fetch from SPARCL for the S1 south-limb
filament, the G206.6+6.1 rim (control) and the E1 filament, by angular distance to the NSNS ridge traces
(data/derived/w4_family_traces.json, ridge points with S/N >= 4, linearly densified to 0.005 deg).

Inputs: data/raw/w4_family_sparcl/find_south.csv (SPARCL find, RA 99-116, Dec -2..11.5),
        data/raw/w4_family_tap/zpix_south.csv, data/derived/spectra_desi_catalog.csv (wave-1 box, for E1),
        windows already on disk (data/derived/spectra_desi_windows, verify_desi_windows).
Selection: S1 d < 0.30 deg; 200 random far controls in the south box; E1 the 700 nearest not on disk.
Output: data/derived/w4_family_sel_queue.csv (queue, near-first, interleaved with controls),
        data/derived/w4_family_sel_summary.json
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_select.py [--ext]
  --ext (second pass, after the first fetch): S1 0.30-0.45 deg and the next 800 E1 coadds within 0.12 deg
        -> data/derived/w4_family_sel_queue_ext.csv
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'


def uv(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def ridge_tree(name, sn_min=4.0, closed=False, maxgap=0.3):
    tr = json.load(open(D / 'w4_family_traces.json'))['filaments'][name]['ridge']
    p = [(r['ra'], r['dec']) for r in tr if r['sn'] >= sn_min]
    p = np.array(p)
    pts = []
    n = len(p)
    for k in range(n if closed else n - 1):
        a, b = p[k], p[(k + 1) % n]
        dd = np.hypot((b[0] - a[0]) * np.cos(np.radians(a[1])), b[1] - a[1])
        if dd > maxgap:            # do not bridge gaps in the traced ridge
            pts.append(a[None]); continue
        m = max(int(dd / 0.005), 1)
        f = np.arange(m)[:, None] / m
        pts.append(a[None] * (1 - f) + b[None] * f)
    pts.append(p[-1:])
    P = np.concatenate(pts)
    return cKDTree(uv(P[:, 0], P[:, 1])), P


def dist(tree, ra, dec):
    d, _ = tree.query(uv(ra, dec))
    return np.degrees(2 * np.arcsin(np.clip(d / 2, 0, 1)))


def have_ids():
    s = set()
    for sub in ('spectra_desi_windows', 'verify_desi_windows', 'w4_family_windows'):
        for f in (D / sub).glob('*_[0-9][0-9][0-9][0-9][0-9].npz'):
            with np.load(f, allow_pickle=False) as d:
                s.update(d['sparcl_id'].tolist())
    return s


def main():
    ext = '--ext' in sys.argv
    have = have_ids()
    rows = list(csv.DictReader(open(ROOT / 'data/raw/w4_family_sparcl/find_south.csv')))
    zp = {}
    for r in csv.DictReader(open(ROOT / 'data/raw/w4_family_tap/zpix_south.csv')):
        zp[(r['targetid'], r['survey'], r['program'])] = r
    ra = np.array([float(r['ra']) for r in rows]); dec = np.array([float(r['dec']) for r in rows])
    tS, _ = ridge_tree('S1'); tG, _ = ridge_tree('G206', closed=True, maxgap=0.35)
    dS, dG = dist(tS, ra, dec), dist(tG, ra, dec)
    ok = np.array([(r['targetid'], r['survey'], r['program']) in zp and
                   int(zp[(r['targetid'], r['survey'], r['program'])]['zwarn']) == 0 and
                   r['specprimary'] in ('True', 'true', '1') for r in rows])
    for r, a, b in zip(rows, dS, dG):
        r['d_S1'], r['d_G206'] = f'{a:.4f}', f'{b:.4f}'
    nearS = [r for r, o, a in zip(rows, ok, dS) if o and a < (0.45 if ext else 0.30) and r['sparcl_id'] not in have]
    nearG = [r for r, o, a, b in zip(rows, ok, dG, dS) if o and a < 0.30 and b >= 0.40 and r['sparcl_id'] not in have]
    far = [r for r, o, a, b in zip(rows, ok, dS, dG) if o and a > 0.8 and b > 0.8 and r['sparcl_id'] not in have]
    nearS.sort(key=lambda r: float(r['d_S1'])); nearG.sort(key=lambda r: float(r['d_G206']))
    for r in nearS: r['set'] = 'S1'
    for r in nearG: r['set'] = 'G206'
    rng = np.random.default_rng(11)
    far = [far[i] for i in rng.permutation(len(far))[:200]]
    if ext:
        far = []                   # controls already fetched in the first pass
    for r in far: r['set'] = 'ctrl_south'
    # E1 from the wave-1 catalogue (main survey, box RA 111-122)
    tE, _ = ridge_tree('E1')
    cat = list(csv.DictReader(open(D / 'spectra_desi_catalog.csv')))
    cra = np.array([float(r['ra']) for r in cat]); cdec = np.array([float(r['dec']) for r in cat])
    m = (cra > 117.5) & (cra < 121.5) & (cdec > 12.5) & (cdec < 17.6)
    dE = np.full(len(cat), 99.0); dE[m] = dist(tE, cra[m], cdec[m])
    nearE_all = [(c, d) for c, d in zip(cat, dE) if d < 0.40 and c['spectype'] in ('GALAXY', 'QSO', 'STAR')
                 and int(c['zwarn']) == 0]
    nearE_have = [c for c, d in nearE_all if c['sparcl_id'] in have]
    nearE = [dict(sparcl_id=c['sparcl_id'], targetid=c['targetid'], survey=c['survey'], program=c['program'],
                  ra=c['ra'], dec=c['dec'], spectype=c['spectype'], redshift=c['redshift'], specprimary='True',
                  d_S1='', d_G206='', d_E1=f'{d:.4f}', set='E1')
             for c, d in sorted(nearE_all, key=lambda x: x[1]) if c['sparcl_id'] not in have][:700]
    if ext:
        nearE = nearE[:0] + [r for r in nearE][:0] or [dict(sparcl_id=c['sparcl_id'], targetid=c['targetid'], survey=c['survey'], program=c['program'], ra=c['ra'], dec=c['dec'], spectype=c['spectype'], redshift=c['redshift'], specprimary='True', d_S1='', d_G206='', d_E1=f'{d:.4f}', set='E1') for c, d in sorted(nearE_all, key=lambda x: x[1]) if c['sparcl_id'] not in have and d < 0.12][:800]
    print('S1 <0.40:', len(nearS), ' G206 <0.30:', len(nearG), ' far ctrl:', len(far))
    print('E1 <0.40 in catalogue:', len(nearE_all), ' already on disk:', len(nearE_have), ' to fetch:', len(nearE))
    for r in nearS + nearG + far:
        r['d_E1'] = ''
    # queue: interleave in batches of 25: S1, E1, G206, ctrl
    B = 25
    lists = [nearS, nearE, nearG, far]
    q, k = [], 0
    while any(k < len(L) for L in lists):
        for L in lists:
            q += L[k:k + B]
        k += B
    fields = ['sparcl_id', 'targetid', 'survey', 'program', 'ra', 'dec', 'spectype', 'redshift', 'specprimary',
              'set', 'd_S1', 'd_G206', 'd_E1']
    with open(D / ('w4_family_sel_queue_ext.csv' if ext else 'w4_family_sel_queue.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        w.writeheader(); w.writerows(q)
    summ = dict(n_queue=len(q), S1=len(nearS), G206=len(nearG), ctrl_south=len(far), E1_fetch=len(nearE),
                E1_on_disk=len(nearE_have), S1_hist_d=np.histogram([float(r['d_S1']) for r in nearS], bins=[0, .05, .1, .2, .3])[0].tolist(),
                E1_max_d_fetched=max(float(r['d_E1']) for r in nearE) if nearE else None,
                note='S1 d<0.30 deg; 200 far controls (>0.8 deg from S1/G206); E1: the 700 nearest coadds not on disk; G206.6+6.1 has no DESI DR1 coadds (min distance 3.1 deg)')
    json.dump(summ, open(D / ('w4_family_sel_summary_ext.json' if ext else 'w4_family_sel_summary.json'), 'w'), indent=1)
    print(summ)


if __name__ == '__main__':
    main()
