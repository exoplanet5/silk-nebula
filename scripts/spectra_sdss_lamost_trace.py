"""Filament geometry helpers for the spectra_sdss_lamost analysis (reads the geometry teammate's
data/derived/silk_trace.json v1.0; falls back to silk_geometry's provisional values)."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import silk_geometry as sg

ROOT = Path('/Users/mickey/astronomy/silk-nebula')


def load():
    tr = sg.load_trace()
    out = {'version': None}
    if tr is None:
        sp = np.array([[116.45, 14.72], [116.65, 16.60], [116.21, 17.74]])
        out.update(E=sp, W=None, envelope=None, halpha=None)
        return out
    out['version'] = tr.get('version')
    out['E'] = np.array([[p['ra'], p['dec']] for p in tr['spine']])
    w = tr.get('strands', {}).get('W')
    out['W'] = np.array([[p['ra'], p['dec']] for p in w['points']]) if w else None
    out['envelope'] = np.array(tr['envelope_polygon']) if tr.get('envelope_polygon') else None
    out['halpha'] = np.array(tr['halpha_patch']) if tr.get('halpha_patch') else None
    out['hwhm'] = tr.get('half_width_deg')
    out['regions'] = {k: (v['ra'], v['dec']) for k, v in tr['regions'].items()}
    return out


def _tangent(ra, dec, ra0, dec0):
    """Gnomonic-ish small-field projection (deg) about (ra0, dec0)."""
    x = (np.asarray(ra) - ra0) * np.cos(np.radians(dec0))
    y = np.asarray(dec) - dec0
    return x, y


def signed_distance(ra, dec, poly):
    """Distance (deg) of points to a polyline; sign + = east (larger RA) of the line at the same Dec.
    Also returns 'inside_span' = the point's Dec lies within the Dec range of the polyline."""
    ra = np.atleast_1d(ra).astype(float)
    dec = np.atleast_1d(dec).astype(float)
    best = np.full(ra.shape, np.inf)
    for (a0, d0), (a1, d1) in zip(poly[:-1], poly[1:]):
        dm = 0.5 * (d0 + d1)
        px, py = _tangent(ra, dec, a0, dm)
        bx, by = _tangent(a1, d1, a0, dm)
        L2 = bx * bx + by * by
        t = np.clip((px * bx + py * by) / L2, 0, 1)
        best = np.minimum(best, np.hypot(px - t * bx, py - t * by))
    o = np.argsort(poly[:, 1])
    ra_line = np.interp(dec, poly[o, 1], poly[o, 0])
    sign = np.where(ra >= ra_line, 1.0, -1.0)
    inside = (dec >= poly[:, 1].min()) & (dec <= poly[:, 1].max())
    return best * sign, inside


def nearest_strand_distance(ra, dec, geo):
    dE, inE = signed_distance(ra, dec, geo['E'])
    if geo.get('W') is None:
        return dE, inE, np.zeros(len(dE), int)
    dW, inW = signed_distance(ra, dec, geo['W'])
    useW = (np.abs(dW) < np.abs(dE)) & inW
    d = np.where(useW, dW, dE)
    return d, inE | inW, useW.astype(int)


def point_in_poly(ra, dec, poly):
    from matplotlib.path import Path as MPath
    return MPath(poly).contains_points(np.column_stack([ra, dec]))
