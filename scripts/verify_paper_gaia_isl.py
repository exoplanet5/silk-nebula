"""verify_paper_gaia_isl.py -- integrated starlight (Gaia DR3 G) per magnitude bin near P1-P4 (tag: paper).

Cone query (G magnitudes only) of Gaia DR3 (I/355/gaiadr3) at the VizieR mirror vizier.cfa.harvard.edu; the NOIRLab
Data Lab aggregated TAP query returned HTTP 502 on 2026-11-16 (log kept as query_log_datalab_failed_502.txt).
Output: data/derived/verify_paper_raw/gaia_isl/isl_<P>_r<rad>.csv (gbin, n, fsum) and query_log.txt
Run: ~/.venvs/astro313/bin/python -I scripts/verify_paper_gaia_isl.py
"""
import datetime, sys
from pathlib import Path
import requests

OUT = Path('/Users/mickey/astronomy/silk-nebula/data/derived/verify_paper_raw/gaia_isl')
OUT.mkdir(parents=True, exist_ok=True)
URL = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
REG = {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72), 'P4': (115.20, 16.18)}
RAD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
import numpy as np

for k, (ra, dec) in REG.items():
    params = {'-source': 'I/355/gaiadr3', '-c': f'{ra}+{dec}', '-c.rd': f'{RAD}', '-out': 'Gmag', '-out.max': 'unlimited'}
    t0 = datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        r = requests.get(URL, params=params, timeout=300)
        txt, status = r.text, r.status_code
    except Exception as e:
        txt, status = '', f'error {e!r}'
    g = []
    for line in txt.splitlines():
        try:
            g.append(float(line.strip()))
        except ValueError:
            pass
    g = np.array(g)
    bins = np.arange(3.0, 22.0, 0.5)
    rows = ['gbin,n,fsum']
    for b0 in bins:
        sel = (g >= b0) & (g < b0 + 0.5)
        rows.append(f'{b0:.1f},{int(sel.sum())},{np.sum(10 ** (-0.4 * g[sel])):.6e}')
    fn = OUT / f'isl_{k}_r{RAD:.2f}.csv'
    fn.write_text('\n'.join(rows) + '\n')
    with open(OUT / 'query_log.txt', 'a') as f:
        f.write(f'{t0} {URL} params={params} status={status} n_stars={g.size} -> {fn.name}\n')
    print(k, status, g.size, np.round(np.percentile(g, [1, 50, 99]), 2) if g.size else '')
