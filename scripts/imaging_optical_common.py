"""Helpers for the imaging_optical analyses (Silk Nebula precovery search).

- target WCS grids (TAN) for wide maps
- HiPS (NSNS DR0.2) tile mosaicking via the tiles' own HPX WCS
- HEALPix (Finkbeiner 2003, LAMBDA) sampling onto a WCS grid
- filament trace geometry: along / across coordinates relative to the spine
"""
import sys, glob, json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
sys.path.insert(0, str(ROOT / 'scripts'))
import silk_geometry as sg  # noqa: E402

RAW = ROOT / 'data/raw'
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'


def tan_wcs(ra0, dec0, width_deg, height_deg, pix_arcmin):
    """North-up, east-left TAN WCS; returns (wcs, (ny, nx))."""
    nx = int(round(width_deg * 60 / pix_arcmin)); ny = int(round(height_deg * 60 / pix_arcmin))
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [(nx + 1) / 2, (ny + 1) / 2]
    w.wcs.cdelt = [-pix_arcmin / 60, pix_arcmin / 60]
    return w, (ny, nx)


def wide_grid(pix_arcmin=2.0):
    return tan_wcs(116.5, 17.0, 17.0, 24.0, pix_arcmin)


def pixel_radec(w, shape):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    ra, dec = w.wcs_pix2world(xx, yy, 0)
    return ra, dec


# ---------------------------------------------------------------- HiPS
def hips_mosaic(name, order, w, shape):
    """Reproject all downloaded NSNS HiPS tiles of <name>/Norder<order> onto (w, shape).
    Tiles are 512x512 with an HPX WCS in their header; BSCALE is applied by astropy."""
    from reproject import reproject_interp
    acc = np.zeros(shape); cnt = np.zeros(shape)
    for f in sorted(glob.glob(str(RAW / f'nsns_hips/{name}/Norder{order}/Npix*.fits'))):
        with fits.open(f) as h:
            data = h[0].data.astype(float)
            tw = WCS(h[0].header)
        arr, fp = reproject_interp((data, tw), w, shape_out=shape, order='bilinear')
        good = (fp > 0) & np.isfinite(arr)
        acc[good] += arr[good]; cnt[good] += 1
    out = np.full(shape, np.nan)
    m = cnt > 0
    out[m] = acc[m] / cnt[m]
    return out


# ---------------------------------------------------------------- HEALPix
def healpix_sample(map_, w, shape, coord='G', nest=False):
    import healpy as hp
    ra, dec = pixel_radec(w, shape)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    if coord == 'G':
        lon, lat = c.galactic.l.deg, c.galactic.b.deg
    else:
        lon, lat = ra, dec
    nside = hp.npix2nside(map_.size)
    return hp.get_interp_val(map_, lon, lat, lonlat=True, nest=nest)


# ---------------------------------------------------------------- trace
def spine():
    """Return (ra, dec, sigma, source_str, preliminary_flag) of the current trace spine."""
    t = sg.load_trace()
    if t is None:
        raise RuntimeError('no trace file')
    sp = t['spine']
    ra = np.array([p['ra'] for p in sp]); dec = np.array([p['dec'] for p in sp])
    sig = np.array([p.get('sigma_deg', np.nan) for p in sp])
    return ra, dec, sig, t.get('version', '?'), bool(t.get('preliminary', True))


def trace_info():
    t = sg.load_trace()
    return t


def along_across(ra, dec, sra=None, sdec=None, nfine=4000):
    """Signed perpendicular distance (deg, + = east of spine) and along-spine distance (deg from the
    southernmost spine point) for each (ra, dec). Uses a densely interpolated spine and local tangent plane."""
    if sra is None:
        sra, sdec, *_ = spine()
    o = np.argsort(sdec); sra, sdec = sra[o], sdec[o]
    # dense spine by interpolation in dec
    dd = np.linspace(sdec.min(), sdec.max(), nfine)
    rr = np.interp(dd, sdec, sra)
    # local projection centred at spine middle
    ra0 = np.median(rr); dec0 = np.median(dd)
    def proj(a, d):
        a = np.radians(a); d = np.radians(d); a0 = np.radians(ra0); d0 = np.radians(dec0)
        cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
        x = np.cos(d) * np.sin(a - a0) / cosc
        y = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
        return np.degrees(x), np.degrees(y)   # x increases to the EAST
    sx, sy = proj(rr, dd)
    seg = np.hypot(np.diff(sx), np.diff(sy)); s_cum = np.concatenate([[0], np.cumsum(seg)])
    px, py = proj(np.ravel(ra), np.ravel(dec))
    from scipy.spatial import cKDTree
    tree = cKDTree(np.c_[sx, sy])
    dist, idx = tree.query(np.c_[px, py])
    idx = np.clip(idx, 1, len(sx) - 2)
    tx = sx[idx + 1] - sx[idx - 1]; ty = sy[idx + 1] - sy[idx - 1]
    tn = np.hypot(tx, ty); tx /= tn; ty /= tn
    dx = px - sx[idx]; dy = py - sy[idx]
    across = dx * ty - dy * tx          # + when point is to the east (x larger) for a northward tangent
    along = s_cum[idx] + dx * tx + dy * ty
    # flag points beyond the spine ends
    beyond = (along < -0.05) | (along > s_cum[-1] + 0.05)
    shp = np.shape(ra)
    return across.reshape(shp), along.reshape(shp), beyond.reshape(shp)


def overlay_trace(ax, w=None, color='#4fc3f7', lw=1.0, ls='--', alpha=0.9, halfwidth=True, transform=None,
                  label='Silk trace'):
    """Plot the spine (and +-half-width) in world coords on a WCSAxes (or with a given transform)."""
    t = sg.load_trace()
    sra, sdec, sig, ver, prelim = spine()
    o = np.argsort(sdec)
    tr = transform if transform is not None else ax.get_transform('world')
    ax.plot(sra[o], sdec[o], ls=ls, color=color, lw=lw, alpha=alpha, transform=tr, label=label)
    if halfwidth and t is not None and 'half_width_deg' in t:
        hw = t['half_width_deg']
        cosd = np.cos(np.radians(sdec[o]))
        for s in (-1, 1):
            ax.plot(sra[o] + s * hw / cosd, sdec[o], ls=':', color=color, lw=0.7, alpha=0.6, transform=tr)


def mark_regions(ax, transform=None, color='#ffffff', size=8, with_p4=True):
    from silk_style import label as slabel
    tr = transform if transform is not None else ax.get_transform('world')
    for k, (ra, dec) in sg.REGIONS.items():
        if k == 'P4' and not with_p4:
            continue
        c = '#ffb84f' if k == 'P4' else color
        ax.plot([ra], [dec], marker='o', ms=4, mfc='none', mec=c, mew=1.0, transform=tr)
        slabel(ax, ra - 0.25, dec, k, color=c, size=size, transform=tr, va='center')


def save_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=float))


# ---------------------------------------------------------------- v1.0 geometry (teammate file, read-only)
def trace_v1():
    """Return dict with spine (E strand) ra/dec, W strand ra/dec (may be empty), Halpha patch polygon, version."""
    t = sg.load_trace()
    sp = t['spine']
    out = dict(version=t.get('version'), preliminary=t.get('preliminary', True),
               ra=np.array([p['ra'] for p in sp]), dec=np.array([p['dec'] for p in sp]))
    W = t.get('strands', {}).get('W', {}).get('points', [])
    out['w_ra'] = np.array([p['ra'] for p in W]); out['w_dec'] = np.array([p['dec'] for p in W])
    out['patch'] = np.array(t.get('halpha_patch', []), float)
    return out


def overlay_v1(ax, tr, spine=True, wstrand=True, patch=True, color='#4fc3f7', lw=0.8):
    from silk_style import HA_C
    g = trace_v1()
    if spine:
        o = np.argsort(g['dec'])
        ax.plot(g['ra'][o], g['dec'][o], '-', color=color, lw=lw, alpha=0.85, transform=tr,
                label=f"[OIII] spine v{g['version']}")
    if wstrand and len(g['w_ra']):
        o = np.argsort(g['w_dec'])
        ax.plot(g['w_ra'][o], g['w_dec'][o], ':', color=color, lw=lw, alpha=0.85, transform=tr, label='W strand')
    if patch and len(g['patch']):
        p = np.vstack([g['patch'], g['patch'][:1]])
        ax.plot(p[:, 0], p[:, 1], '-', color=HA_C, lw=0.8, transform=tr, label='Hα patch')


def in_patch(ra, dec):
    from matplotlib.path import Path as MPath
    g = trace_v1()
    return MPath(g['patch']).contains_points(np.c_[np.ravel(ra), np.ravel(dec)]).reshape(np.shape(ra))
