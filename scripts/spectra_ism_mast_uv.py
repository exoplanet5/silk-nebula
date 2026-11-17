#!/usr/bin/env python
"""spectra_ism: inventory of archival UV spectra (FUSE, IUE, HST STIS/COS/GHRS) around the Silk Nebula.

MAST CAOM ObsCore TAP: https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync (accessed 2026-11-16; the old vao.stsci.edu CAOMTAP URL returned 404).
Box RA 104-128, Dec +4..+34 (same as the optical archive search). Output:
data/raw/spectra_ism_mast/mast_uv_obscore_box.csv, data/derived/spectra_ism_mast_uv_targets.csv
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_mast_uv.py
"""
import csv
import io
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, fetch, filament_coords  # noqa: E402
from spectra_ism_absorption_catalogs_lib import spine_sep  # noqa: E402

OUT = RAW / 'spectra_ism_mast'
OUT.mkdir(parents=True, exist_ok=True)
# A first unfiltered query (mast_uv_obscore_box.csv) hit the 100000-row cap, dominated by solar-system targets
# (Jupiter/Io/comets crossing Gemini). Two narrower queries are used instead; moving targets are dropped by name.
BASE = """SELECT obs_collection, instrument_name, target_name, s_ra, s_dec, em_min, em_max, em_res_power, t_exptime,
t_min, obs_id, dataproduct_type, access_url FROM ivoa.obscore
WHERE s_ra BETWEEN 104 AND 128 AND s_dec BETWEEN 4 AND 34 AND dataproduct_type = 'spectrum' AND em_min < 3.2e-7 AND """
QUERIES = {'mast_uv_fuse_iue_hut.csv': BASE + "obs_collection IN ('FUSE', 'IUE', 'HUT')",
           'mast_uv_hst.csv': BASE + "obs_collection = 'HST' AND (instrument_name LIKE 'COS%' OR instrument_name LIKE 'STIS%' OR instrument_name LIKE 'GHRS%' OR instrument_name LIKE 'FOS%')"}
SOLSYS = ('JUPITER', 'IO-', 'IO_', 'CALLISTO', 'GANYMEDE', 'EUROPA', 'HARTLEY', 'BORRELLY', 'COMET', 'SATURN', 'MARS',
          'VENUS', 'URANUS', 'NEPTUNE', 'CCDFLAT', 'NONE', 'ANY', 'TITAN', 'MOON', 'TEMPEL', 'ENCKE', 'ISON', 'P/')
import time  # noqa: E402
rows = []
for fn, q in QUERIES.items():
    dest = OUT / fn
    if '--cached' not in sys.argv or not dest.exists():
        ok, n, st = fetch('https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync', dest,
                          params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': q}, timeout=300)
        print(fn, ok, n, st)
        time.sleep(2)
    if dest.exists():
        rr = list(csv.DictReader(io.StringIO(dest.read_text(errors='replace'))))
        print(fn, 'rows', len(rr))
        rows += [r for r in rr if not any(k in r['target_name'].upper() for k in SOLSYS)]
print('UV spectra rows:', len(rows))
grp = defaultdict(list)
for r in rows:
    try:
        ra, de = float(r['s_ra']), float(r['s_dec'])
    except ValueError:
        continue
    grp[(r['target_name'].strip().upper(), round(ra, 2), round(de, 2))].append(r)
out = []
for (tn, ra, de), rr in grp.items():
    em_min = min(float(x['em_min'] or 'nan') for x in rr) * 1e10 if rr[0]['em_min'] else np.nan
    em_max = max(float(x['em_max'] or 'nan') for x in rr) * 1e10 if rr[0]['em_max'] else np.nan
    insts = sorted(set(f"{x['obs_collection']}/{x['instrument_name']}" for x in rr))
    cls = sorted(set((x.get('target_classification') or '')[:40] for x in rr))
    covers_ovi = any(x['obs_collection'] == 'FUSE' for x in rr)
    covers_civ = any((float(x['em_min'] or 1) * 1e10 < 1548) and (float(x['em_max'] or 0) * 1e10 > 1551) for x in rr)
    res = max(float(x['em_res_power'] or 0) for x in rr)
    out.append(dict(target=tn, ra=ra, dec=de, n_obs=len(rr), instruments='|'.join(insts), classification='|'.join(cls), lam_min_A=round(em_min, 1),
                    lam_max_A=round(em_max, 1), max_res_power=res, covers_OVI_1032=covers_ovi, covers_CIV_1548=covers_civ,
                    total_exptime_s=round(sum(float(x['t_exptime'] or 0) for x in rr), 1)))
ra = np.array([o['ra'] for o in out])
de = np.array([o['dec'] for o in out])
sep = spine_sep(ra, de)
dp, dsp, bey = filament_coords(ra, de)
for o, s_, d_, ds_ in zip(out, sep, dp, dsp):
    o['sep_spine_deg'] = round(float(s_), 3)
    o['dperp_deg'] = round(float(d_), 3)
    o['dec_spine'] = round(float(ds_), 2)
out.sort(key=lambda o: o['sep_spine_deg'])
with open(DER / 'spectra_ism_mast_uv_targets.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
for o in out[:60]:
    print(f"{o['target'][:22]:22s} {o['ra']:8.2f} {o['dec']:+7.2f} sep={o['sep_spine_deg']:5.2f} n={o['n_obs']:3d} "
          f"{o['instruments'][:50]:50s} {o['lam_min_A']:7.0f}-{o['lam_max_A']:7.0f} R={o['max_res_power']:.0f} OVI={o['covers_OVI_1032']} CIV={o['covers_CIV_1548']}")
