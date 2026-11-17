"""verify_hiwham: shared helpers (trace geometry, rigid trace transforms for nulls, LSR conversion, GALFA mosaic).

Independent re-implementation for the wave-2 HI/WHAM verification (tag hiwham). Reads only
data/derived/silk_trace.json (geometry) and raw survey files; does not import any wave-1 analysis module.

Conventions
-----------
* Signed perpendicular offset d (deg): d > 0 on the CONVEX / outer side of the trace. For the real trace this is
  EAST (larger RA), as in wave 1. For null copies (rotated / translated / mirrored) d > 0 is the side that maps
  from the original east side under the rigid motion (mirror copies get their sign flipped so that d > 0 stays
  the convex side of the copy).
* Along-spine coordinate s (deg) from the south end of the E spine.
* V_LSR: standard solar motion 20 km/s toward RA 18h, Dec +30 (B1900), the convention of both WHAM and GALFA.
"""
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord, FK4
from astropy.time import Time
import astropy.units as u
from scipy.spatial import cKDTree

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw'
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'
SCR = Path('/private/tmp/claude-502/-Users-mickey-astronomy-silk-nebula/c7d7dd5c-f92f-4bd3-b591-acbd37ba6cef/scratchpad/hiwham')
SCR.mkdir(parents=True, exist_ok=True)

TRACE = json.loads((DER / 'silk_trace.json').read_text())


# ----------------------------------------------------------------------------------------------- velocities
_APEX = SkyCoord(ra=270 * u.deg, dec=30 * u.deg, frame=FK4(equinox=Time('B1900'))).icrs


def vlsr_minus_vhelio(ra, dec):
    """V_LSR - V_helio (km/s) for the standard solar motion (20 km/s to 18h,+30 B1900)."""
    c = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg)
    return 20.0 * np.cos(c.separation(_APEX).rad)


# ----------------------------------------------------------------------------------------------- trace geometry
def strand(name='E'):
    pts = TRACE['strands'][name]['points'] if name != 'spine' else TRACE['spine']
    ra = np.array([p['ra'] for p in pts], float)
    dec = np.array([p['dec'] for p in pts], float)
    o = np.argsort(dec)
    return ra[o], dec[o]


def tp_fwd(ra, dec, ra0, dec0):
    """Gnomonic projection (deg); xi increases to the EAST."""
    r, d = np.radians(ra), np.radians(dec)
    r0, d0 = np.radians(ra0), np.radians(dec0)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
    xi = np.cos(d) * np.sin(r - r0) / cosc
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc
    return np.degrees(xi), np.degrees(eta)


def tp_inv(xi, eta, ra0, dec0):
    x, y = np.radians(xi), np.radians(eta)
    r0, d0 = np.radians(ra0), np.radians(dec0)
    rho = np.hypot(x, y)
    c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        dec = np.arcsin(np.cos(c) * np.sin(d0) + np.where(rho > 0, y * np.sin(c) * np.cos(d0) / rho, 0))
        ra = r0 + np.arctan2(x * np.sin(c), rho * np.cos(d0) * np.cos(c) - y * np.sin(d0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(dec)


def densify(ra, dec, step=0.01):
    """Densify a polyline in its own tangent plane (deg step along arc)."""
    ra0, dec0 = np.mean(ra), np.mean(dec)
    x, y = tp_fwd(ra, dec, ra0, dec0)
    sl = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    sd = np.arange(0, sl[-1] + 1e-9, step)
    xd, yd = np.interp(sd, sl, x), np.interp(sd, sl, y)
    rd, dd = tp_inv(xd, yd, ra0, dec0)
    return rd, dd


class Trace:
    """A (possibly transformed) copy of a strand polyline; computes (d, s) for arbitrary sky positions.

    The copy is defined by a rigid motion of the ORIGINAL polyline in the tangent plane of the original
    centroid: rotate by `rot` deg (counter-clockwise, i.e. N through E... in (xi,eta) math sense), optional
    E-W mirror, then re-deproject at a new centre (ra_c, dec_c). Identity: rot=0, mirror=False, centre=None.
    """

    def __init__(self, ra, dec, rot=0.0, mirror=False, centre=None, step=0.01):
        rd, dd = densify(ra, dec, step)
        self.ra0, self.dec0 = float(np.mean(ra)), float(np.mean(dec))
        x, y = tp_fwd(rd, dd, self.ra0, self.dec0)
        if mirror:
            x = -x
        th = np.radians(rot)
        xr = np.cos(th) * x - np.sin(th) * y
        yr = np.sin(th) * x + np.cos(th) * y
        self.mirror = mirror
        self.rot = rot
        self.rac, self.decc = (self.ra0, self.dec0) if centre is None else centre
        self.ra, self.dec = tp_inv(xr, yr, self.rac, self.decc)
        # work in the tangent plane of the copy's own centre
        self.x, self.y = tp_fwd(self.ra, self.dec, self.rac, self.decc)
        self.s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(self.x), np.diff(self.y)))])
        self.L = float(self.s[-1])
        tx = np.gradient(self.x); ty = np.gradient(self.y)
        n = np.hypot(tx, ty)
        self.tx, self.ty = tx / n, ty / n
        self.tree = cKDTree(np.c_[self.x, self.y])

    def coords(self, ra, dec, dmax=4.0):
        """Return d (deg, + = convex side), s (deg), inside (bool: projection falls within the polyline)."""
        px, py = tp_fwd(ra, dec, self.rac, self.decc)
        sh = np.shape(px)
        px, py = np.ravel(px), np.ravel(py)
        dist, k = self.tree.query(np.c_[px, py], distance_upper_bound=dmax + 0.05)
        ok = np.isfinite(dist)
        k = np.where(ok, k, 0)
        vx, vy = px - self.x[k], py - self.y[k]
        along = vx * self.tx[k] + vy * self.ty[k]        # residual along-tangent (small)
        cross = self.tx[k] * vy - self.ty[k] * vx        # >0: pixel to the LEFT of S->N direction = WEST
        d = -cross
        if self.mirror:
            d = -d
        s = self.s[k] + along
        inside = ok & (k > 0) & (k < len(self.x) - 1)
        d = np.where(ok, d, np.nan)
        s = np.where(ok, s, np.nan)
        return d.reshape(sh), s.reshape(sh), inside.reshape(sh)


def spine_dec_of_s(s):
    ra, dec = strand('E')
    tr = Trace(ra, dec)
    return np.interp(s, tr.s, tr.dec)


# ----------------------------------------------------------------------------------------------- GALFA mosaic
GALFA_BLOCK = 2            # arcmin block averaging of the 1' GALFA lattice
GALFA_RA = (108.0, 125.0)
GALFA_DEC = (4.0, 30.0)


def galfa_cube(rebuild=False):
    """Return (cube[v, y, x] K, v_kms, ra_1d, dec_1d) on a 2' plate-carree lattice (cached in scratchpad).

    Built from the nine raw DR2 'wide' cutouts by index arithmetic (header CRVAL/CRPIX/CDELT, linear RA/Dec
    axes on GALFA's global 1' lattice), then block-averaged 2x2 (nanmean). Missing tiles stay NaN.
    """
    cache = SCR / 'galfa_mosaic_2arcmin.npz'
    if cache.exists() and not rebuild:
        z = np.load(cache)
        return z['cube'], z['v'], z['ra'], z['dec']
    from astropy.io import fits
    import glob
    files = sorted(glob.glob(str(RAW / 'galfa_hi_dr2' / 'galfa_dr2_wide_*.fits')))
    nx = int(round((GALFA_RA[1] - GALFA_RA[0]) * 60))
    ny = int(round((GALFA_DEC[1] - GALFA_DEC[0]) * 60))
    # 1' lattice centres: RA_i = 125 - (i+0.5)/60 ; Dec_j = 4 + (j+0.5)/60
    cube = None
    for f in files:
        with fits.open(f) as h:
            hd = h[0].header
            d = h[0].data.astype(np.float32)
        if cube is None:
            nv = d.shape[0]
            v = (np.arange(nv) + 1 - hd['CRPIX3']) * hd['CDELT3'] / 1e3 + hd['CRVAL3'] / 1e3
            cube = np.full((nv, ny, nx), np.nan, np.float32)
        ra_pix0 = hd['CRVAL1'] + (1 - hd['CRPIX1']) * hd['CDELT1']
        de_pix0 = hd['CRVAL2'] + (1 - hd['CRPIX2']) * hd['CDELT2']
        i0 = int(round((GALFA_RA[1] - ra_pix0) * 60 - 0.5))
        j0 = int(round((de_pix0 - GALFA_DEC[0]) * 60 - 0.5))
        # sanity: lattice alignment
        assert abs((GALFA_RA[1] - ra_pix0) * 60 - 0.5 - i0) < 0.02, f
        assert abs((de_pix0 - GALFA_DEC[0]) * 60 - 0.5 - j0) < 0.02, f
        assert hd['CDELT1'] < 0 and hd['CDELT2'] > 0
        tny, tnx = d.shape[1:]
        xs = slice(max(i0, 0), min(i0 + tnx, nx)); ys = slice(max(j0, 0), min(j0 + tny, ny))
        sub = d[:, ys.start - j0:ys.stop - j0, xs.start - i0:xs.stop - i0]
        cur = cube[:, ys, xs]
        cube[:, ys, xs] = np.where(np.isfinite(sub), sub, cur)
    b = GALFA_BLOCK
    nv = cube.shape[0]
    cb = np.nanmean(cube.reshape(nv, ny // b, b, nx // b, b), axis=(2, 4)).astype(np.float32)
    ra = GALFA_RA[1] - (np.arange(nx // b) + 0.5) * b / 60.0
    dec = GALFA_DEC[0] + (np.arange(ny // b) + 0.5) * b / 60.0
    np.savez(cache, cube=cb, v=v, ra=ra, dec=dec)
    return cb, v, ra, dec


def nhi_slice(cube, v, lo, hi):
    """Optically thin N(HI) (cm^-2) in [lo, hi] km/s with fractional channel weights."""
    dv = abs(v[1] - v[0])
    e0, e1 = v - dv / 2, v + dv / 2
    w = np.clip(np.minimum(e1, hi) - np.maximum(e0, lo), 0, None)   # km/s of each channel inside window
    k = w > 0
    out = 1.823e18 * np.tensordot(w[k], np.nan_to_num(cube[k]), axes=(0, 0))
    bad = ~np.isfinite(cube[k]).all(axis=0)
    out[bad] = np.nan
    return out


def binned_profile(d, val, edges, sel=None, stat='mean', min_n=5):
    """Binned profile of val vs d (nan-aware)."""
    m = np.isfinite(d) & np.isfinite(val)
    if sel is not None:
        m &= sel
    idx = np.digitize(d[m], edges) - 1
    nb = len(edges) - 1
    ok = (idx >= 0) & (idx < nb)
    idx, x = idx[ok], val[m][ok]
    n = np.bincount(idx, minlength=nb)
    if stat == 'mean':
        s = np.bincount(idx, weights=x, minlength=nb)
        p = np.where(n >= min_n, s / np.maximum(n, 1), np.nan)
    else:
        p = np.full(nb, np.nan)
        order = np.argsort(idx, kind='stable')
        xs, ids = x[order], idx[order]
        bounds = np.searchsorted(ids, np.arange(nb + 1))
        for i in range(nb):
            if n[i] >= min_n:
                p[i] = np.median(xs[bounds[i]:bounds[i + 1]])
    return p, n
