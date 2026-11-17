"""verify_xrayfuv shared helpers (independent of wave-1 code): sky geometry, trace placement,
perpendicular offsets, catalogue reading.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'

# Monogem centres (RA, Dec deg)
CENTRES = {
    'PSR': (104.9506, 14.2393),       # PSR B0656+14 (ATNF); used by wave-1 imaging_mw
    'K1': (107.08, 13.50),            # Knies+2024 single-SNR centre, R = 11.97 deg
    'K2': (104.15, 15.30),            # Knies+2024 two-SNR Monogem centre, R = 11.00 deg
}


def radec_to_vec(ra, dec):
    ra = np.radians(ra); dec = np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def vec_to_radec(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    ra = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360
    dec = np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1)))
    return ra, dec


def sep_pa(ra0, dec0, ra, dec):
    """Angular separation (deg) and position angle (deg E of N) of (ra,dec) from (ra0,dec0)."""
    r0, d0, r, d = map(np.radians, (ra0, dec0, ra, dec))
    cs = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
    sep = np.degrees(np.arccos(np.clip(cs, -1, 1)))
    pa = np.degrees(np.arctan2(np.sin(r - r0) * np.cos(d),
                               np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0))) % 360
    return sep, pa


def offset_point(ra0, dec0, sep, pa):
    """Point at separation sep (deg) and PA (deg) from (ra0, dec0)."""
    r0, d0, s, p = map(np.radians, (ra0, dec0, sep, pa))
    d = np.arcsin(np.sin(d0) * np.cos(s) + np.cos(d0) * np.sin(s) * np.cos(p))
    r = r0 + np.arctan2(np.sin(p) * np.sin(s) * np.cos(d0), np.cos(s) - np.sin(d0) * np.sin(d))
    return np.degrees(r) % 360, np.degrees(d)


def load_spine():
    t = json.loads((DER / 'silk_trace.json').read_text())
    sp = t['spine']
    ra = np.array([p['ra'] for p in sp]); dec = np.array([p['dec'] for p in sp])
    return ra, dec, t


def densify(ra, dec, step=0.01):
    """Densify a polyline along great circles to ~step deg spacing; returns ra, dec, s (deg)."""
    v = radec_to_vec(ra, dec)
    out = [v[0]]
    for a, b in zip(v[:-1], v[1:]):
        ang = np.arccos(np.clip(np.dot(a, b), -1, 1))
        n = max(int(np.ceil(np.degrees(ang) / step)), 1)
        for k in range(1, n + 1):
            f = k / n
            w = (np.sin((1 - f) * ang) * a + np.sin(f * ang) * b) / np.sin(ang) if ang > 0 else b
            out.append(w)
    out = np.array(out)
    seg = np.degrees(np.arccos(np.clip(np.sum(out[1:] * out[:-1], axis=1), -1, 1)))
    s = np.concatenate([[0], np.cumsum(seg)])
    r, d = vec_to_radec(out)
    return r, d, s


def rotate_about(ra, dec, ra_c, dec_c, theta_deg):
    """Rotate points about the axis through (ra_c, dec_c) by theta (deg, right-hand rule)."""
    k = radec_to_vec(ra_c, dec_c)
    v = radec_to_vec(ra, dec)
    th = np.radians(theta_deg)
    vr = v * np.cos(th) + np.cross(k, v) * np.sin(th) + np.outer(v @ k, k) * (1 - np.cos(th))
    return vec_to_radec(vr)


def translate(ra, dec, dra_deg, ddec_deg, ra_ref, dec_ref):
    """Rigidly move a curve: rotation taking (ra_ref,dec_ref) to the point offset by
    (dra*cos dec, ddec) in the tangent plane at the reference point."""
    ra2, dec2 = offset_point(ra_ref, dec_ref, np.hypot(dra_deg, ddec_deg),
                             np.degrees(np.arctan2(dra_deg, ddec_deg)))
    a = radec_to_vec(ra_ref, dec_ref); b = radec_to_vec(ra2, dec2)
    axis = np.cross(a, b); n = np.linalg.norm(axis)
    if n < 1e-12:
        return ra, dec
    axis /= n; th = np.arcsin(min(n, 1.0))
    v = radec_to_vec(ra, dec)
    vr = v * np.cos(th) + np.cross(axis, v) * np.sin(th) + np.outer(v @ axis, axis) * (1 - np.cos(th))
    return vec_to_radec(vr)


class CurveCoords:
    """Signed perpendicular offset x (deg; x>0 on the right-hand side when walking from the
    first to the last point, i.e. EAST for the S->N Silk spine) and along-curve s for pixels."""
    def __init__(self, ra, dec, step=0.01):
        from scipy.spatial import cKDTree
        r, d, s = densify(ra, dec, step)
        self.v = radec_to_vec(r, d); self.s = s; self.ra = r; self.dec = d
        tang = np.gradient(self.v, axis=0)
        self.tang = tang / np.linalg.norm(tang, axis=1, keepdims=True)
        self.tree = cKDTree(self.v)

    def coords(self, pv, maxdist_deg=4.0):
        dist, idx = self.tree.query(pv, distance_upper_bound=2 * np.sin(np.radians(maxdist_deg) / 2))
        ok = np.isfinite(dist)
        idx_c = np.where(ok, idx, 0)
        t = self.v[idx_c]; tg = self.tang[idx_c]
        # sign: walking N along the spine, the east side is on the LEFT when viewed on the sky
        # (RA increases to the left). n = t x tangent points to the left of travel (east for S->N).
        nrm = np.cross(t, tg)
        sgn = -np.sign(np.sum(pv * nrm, axis=1))  # verified: P4 (west) -> negative
        ang = 2 * np.degrees(np.arcsin(np.clip(dist / 2, 0, 1)))
        x = np.where(ok, sgn * ang, np.nan)
        s = np.where(ok, self.s[idx_c], np.nan)
        end = (idx_c == 0) | (idx_c == len(self.s) - 1)
        x[end] = np.nan
        return x, s


def read_vizier_tsv(path, cols_float=None):
    lines = [ln for ln in Path(path).read_text(errors='replace').splitlines() if ln and not ln.startswith('#')]
    if len(lines) < 4:
        return {}
    hdr = lines[0].split('\t')
    rows = [ln.split('\t') for ln in lines[3:]]
    rows = [r for r in rows if len(r) == len(hdr)]
    out = {}
    for j, h in enumerate(hdr):
        col = [r[j].strip() for r in rows]
        try:
            out[h] = np.array([float(c) if c else np.nan for c in col])
        except ValueError:
            out[h] = np.array(col)
    return out
