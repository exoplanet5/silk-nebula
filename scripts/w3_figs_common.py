"""Shared helpers for the wave-3 composite figures (tag `figs`).

Only reads existing project files; writes nothing by itself. Figures use scripts/silk_style.py
(dark palette, Menlo, haloed labels, credit via finish()).
"""
import sys
import json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
DER = ROOT / 'data/derived'
RAW = ROOT / 'data/raw'
FIG = ROOT / 'figures'
sys.path.insert(0, str(ROOT / 'scripts'))
from silk_style import *  # noqa: E402,F401,F403
import silk_style as ss  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

C_KMS = 299792.458


def load_trace():
    """silk_trace.json v1.0 -> dict with E/W/W2 (N,2) arrays, Ha patch polygons, regions, LAMOST field."""
    tr = json.loads((DER / 'silk_trace.json').read_text())

    def pts(p):
        return np.array([[q['ra'], q['dec']] if isinstance(q, dict) else [q[0], q[1]] for q in p], float)
    out = {'tr': tr,
           'E': pts(tr['strands']['E']['points']),
           'W': pts(tr['strands']['W']['points']),
           'W2': pts(tr['strands']['W2']['points']),
           'patch': np.array(tr['halpha_patch'], float),
           'patch30': np.array(tr['halpha_patch_info']['outer_polygon_30pct'], float),
           'regions': {k: (v['ra'], v['dec']) for k, v in tr['regions'].items()},
           'lamost': (tr['lamost_field']['ra'], tr['lamost_field']['dec'], tr['lamost_field']['radius_deg'])}
    return out


def small_circle(ra0, dec0, r_deg, n=361):
    """Points of a small circle of angular radius r about (ra0, dec0)."""
    ra0, dec0, r = np.radians(ra0), np.radians(dec0), np.radians(r_deg)
    pa = np.linspace(0, 2 * np.pi, n)
    dec = np.arcsin(np.sin(dec0) * np.cos(r) + np.cos(dec0) * np.sin(r) * np.cos(pa))
    ra = ra0 + np.arctan2(np.sin(pa) * np.sin(r) * np.cos(dec0), np.cos(r) - np.sin(dec0) * np.sin(dec))
    return np.degrees(ra), np.degrees(dec)


def ang_sep(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def json_dump(obj, path):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    Path(path).write_text(json.dumps(obj, indent=1, default=conv))


def panel_tag(ax, text, x=0.012, y=0.985, size=9.5, color='#cfd6dd'):
    """Short bold panel label at the upper left inside an axes."""
    return ax.text(x, y, text, transform=ax.transAxes, ha='left', va='top', fontsize=size,
                   fontweight='bold', color=color, path_effects=[ss.HALO])
