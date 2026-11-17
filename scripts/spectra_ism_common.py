"""spectra_ism: shared helpers (paths, velocity frames, filament-relative coordinates, HTTP).

All network fetches go through fetch() which records URL, status, bytes and access date in
data/raw/spectra_ism_fetch_log.jsonl (one JSON object per line). Downloaded bytes are written
to data/raw/<source>/ and only ever parsed as data (FITS/CSV/VOTable/text), never executed.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import ROOT, REGIONS, CENTER, load_trace, gal  # noqa: E402,F401

RAW = ROOT / 'data/raw'
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'
LOG = RAW / 'spectra_ism_fetch_log.jsonl'
ACCESS_DATE = '2026-11-16'
UA = {'User-Agent': 'silk-nebula-research/1.0 (python-requests)'}

# Standard solar motion used by Wu+2026 (20 km/s toward RA 18h03m50.29s, Dec +30d00m16.8s, J2000)
APEX_STD = SkyCoord('18h03m50.29s', '+30d00m16.8s')
V_SUN_STD = 20.0


def vlsr_minus_vhelio(ra, dec):
    """V_LSR - V_helio (km/s) for the standard solar motion (the convention of WHAM and Wu+2026)."""
    c = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg)
    return V_SUN_STD * np.cos(c.separation(APEX_STD).rad)


def fetch(url, dest, params=None, timeout=120, retries=3, method='GET', data=None, min_bytes=1):
    """Download url -> dest (Path). Returns (ok, nbytes, status). Logs every attempt."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    status, n, err = None, 0, ''
    for k in range(retries):
        try:
            if method == 'POST':
                r = requests.post(url, data=data or params, timeout=timeout, headers=UA)
            else:
                r = requests.get(url, params=params, timeout=timeout, headers=UA)
            status = r.status_code
            if r.ok and len(r.content) >= min_bytes:
                dest.write_bytes(r.content)
                n = len(r.content)
                break
            err = f'HTTP {status}, {len(r.content)} bytes'
        except Exception as ex:  # network errors are recorded, not raised
            err = repr(ex)[:300]
        time.sleep(2 + 3 * k)
    rec = {'url': url, 'params': params if method == 'GET' else None, 'post': data if method == 'POST' else None,
           'dest': str(dest.relative_to(ROOT)) if str(dest).startswith(str(ROOT)) else str(dest),
           'status': status, 'bytes': n, 'error': err if n == 0 else '', 'access_date': ACCESS_DATE,
           'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    with open(LOG, 'a') as fh:
        fh.write(json.dumps(rec) + '\n')
    return n > 0, n, status


# ------------------------------------------------------------------------------------------
# Filament-relative coordinates from the geometry teammate's trace (E strand = spine)
# ------------------------------------------------------------------------------------------
def spine_arrays():
    tr = load_trace()
    if tr is None:  # provisional straight-line fallback through P3-P1-P2
        dec = np.linspace(12, 28.25, 66)
        ra = np.interp(dec, [14.72, 16.60, 17.74], [116.45, 116.65, 116.21])
        return ra, dec, None
    sp = tr['spine']
    return np.array([p['ra'] for p in sp]), np.array([p['dec'] for p in sp]), tr


def filament_coords(ra, dec):
    """Signed perpendicular offset (deg, + = east/larger RA) from the spine and the Dec of the
    nearest spine point; points beyond the spine Dec range get offsets measured from the end
    point (flag beyond=True)."""
    sra, sdec, _ = spine_arrays()
    ra, dec = np.atleast_1d(ra).astype(float), np.atleast_1d(dec).astype(float)
    # densify the spine
    t = np.linspace(0, 1, 2000)
    fd = np.interp(t, np.linspace(0, 1, len(sdec)), sdec)
    fr = np.interp(t, np.linspace(0, 1, len(sra)), sra)
    cs = SkyCoord(fr * u.deg, fd * u.deg)
    cp = SkyCoord(ra * u.deg, dec * u.deg)
    idx, sep, _ = cp.match_to_catalog_sky(cs)
    sgn = np.sign((ra - fr[idx]) * np.cos(np.deg2rad(dec)))
    sgn[sgn == 0] = 1
    beyond = (idx == 0) | (idx == len(t) - 1)
    return sgn * sep.deg, fd[idx], beyond


def in_polygon(ra, dec, poly):
    """Even-odd point-in-polygon in (RA, Dec) plane (fine for the small, convex-ish envelope)."""
    from matplotlib.path import Path as MPath
    p = MPath(np.asarray(poly))
    return p.contains_points(np.c_[np.atleast_1d(ra), np.atleast_1d(dec)])


def draw_trace(ax, color='#4fc3f7', lw=1.0, alpha=0.9, envelope=True, regions=True, lb=False, label=True):
    """Overlay spine, envelope and P1-P4 on an RA/Dec (or l,b) axis."""
    from silk_style import label as hlabel
    sra, sdec, tr = spine_arrays()
    if lb:
        x, y = gal(sra, sdec)
    else:
        x, y = sra, sdec
    ax.plot(x, y, '-', color=color, lw=lw, alpha=alpha, label='[O III] spine' if label else None)
    if envelope and tr is not None:
        pol = np.asarray(tr['envelope_polygon'])
        px, py = (gal(pol[:, 0], pol[:, 1]) if lb else (pol[:, 0], pol[:, 1]))
        ax.plot(np.r_[px, px[0]], np.r_[py, py[0]], ':', color=color, lw=0.7, alpha=0.7)
    if regions:
        for k, (r, d) in REGIONS.items():
            x0, y0 = (gal(r, d) if lb else (r, d))
            x0, y0 = float(np.atleast_1d(x0)[0]), float(np.atleast_1d(y0)[0])
            ax.plot(x0, y0, 'o', ms=4, mfc='none', mec='#ffffff' if k != 'P4' else '#ffb84f', mew=1.0)
            hlabel(ax, x0, y0, ' ' + k, size=7.5, va='bottom')
