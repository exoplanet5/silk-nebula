"""spectra_ism: parser for multi-table VizieR asu-tsv files (data only; no code from the files is run).

read_asu(path) -> dict {table_name: list of row dicts}. Blank fields become ''.
"""
from pathlib import Path


def read_asu(path):
    lines = Path(path).read_text(errors='replace').splitlines()
    tables, cur, hdr, state = {}, None, None, 0
    for ln in lines:
        if ln.startswith('#Name:'):
            cur = ln.split(':', 1)[1].strip()
            hdr, state = None, 1
            continue
        if ln.startswith('#') or cur is None:
            continue
        if not ln.strip():
            hdr, state = None, 1
            continue
        if state == 1:                       # header row
            hdr = ln.split('\t')
            state = 2
            continue
        if state == 2:                       # units row
            state = 3
            continue
        if state == 3 and set(ln.replace('\t', '').strip()) <= {'-'}:
            state = 4
            continue
        if state >= 3 and hdr is not None:
            vals = ln.split('\t')
            tables.setdefault(cur, []).append({h.strip(): (vals[i].strip() if i < len(vals) else '') for i, h in enumerate(hdr)})
            state = 4
    return tables


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')
