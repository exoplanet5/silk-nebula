"""Common grids, reprojection and trace-relative profile tools for the imaging_mw study.

Grids
  EQ : plate-carree RA/Dec, RA 108-125, Dec +4..+30, 0.05 deg (3') pixels (x = RA decreasing to the right)
  GAL: plate-carree l/b,   l 185-225, b 0..+35,    0.10 deg pixels (x = l decreasing to the right)

Trace: read via silk_geometry.load_trace() at call time (so the latest teammate trace is used);
falls back to the P3-P1-P2 provisional spine if no trace file is present.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import REGIONS, CENTER, ROOT, load_trace, MONOGEM_PULSAR, GEMINGA, gal  # noqa: E402

warnings.simplefilter('ignore', FITSFixedWarning)
RAW = ROOT / 'data/raw'
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'
DER.mkdir(exist_ok=True, parents=True)


def car_wcs(frame, lon0, lon1, lat0, lat1, pix):
    """Plate-carree WCS with CRVAL2 = 0 so that pixel rows are lines of constant latitude."""
    nx = int(round((lon1 - lon0) / pix)); ny = int(round((lat1 - lat0) / pix))
    w = WCS(naxis=2)
    if frame == 'icrs':
        w.wcs.ctype = ['RA---CAR', 'DEC--CAR']; w.wcs.radesys = 'ICRS'
    else:
        w.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']
    lonc = 0.5 * (lon0 + lon1)
    w.wcs.crval = [lonc, 0.0]
    w.wcs.cdelt = [-pix, pix]
    w.wcs.crpix = [nx / 2 + 0.5, 1 - (lat0 + pix / 2) / pix]
    return w, (ny, nx)


EQ = dict(frame='icrs', lon=(108.0, 125.0), lat=(4.0, 30.0), pix=0.05)
GALG = dict(frame='galactic', lon=(185.0, 225.0), lat=(0.0, 35.0), pix=0.10)
EQ_WCS, EQ_SHAPE = car_wcs('icrs', *EQ['lon'], *EQ['lat'], EQ['pix'])
GAL_WCS, GAL_SHAPE = car_wcs('galactic', *GALG['lon'], *GALG['lat'], GALG['pix'])
EQ_EXTENT = [EQ['lon'][1], EQ['lon'][0], EQ['lat'][0], EQ['lat'][1]]     # imshow extent (RA decreasing)
GAL_EXTENT = [GALG['lon'][1], GALG['lon'][0], GALG['lat'][0], GALG['lat'][1]]


def grid_coords(which='eq'):
    w, shp = (EQ_WCS, EQ_SHAPE) if which == 'eq' else (GAL_WCS, GAL_SHAPE)
    yy, xx = np.mgrid[0:shp[0], 0:shp[1]]
    lon, lat = w.wcs_pix2world(xx, yy, 0)
    return lon, lat


def reproject_to(path_or_hdu, which='eq', hdu=0, method='interp', nan_value=None, order='bilinear'):
    """Reproject a FITS image onto the EQ or GAL grid. Returns float32 array."""
    from reproject import reproject_interp, reproject_adaptive
    if isinstance(path_or_hdu, (str, Path)):
        hl = fits.open(path_or_hdu)
        h = hl[hdu]
        data, header = h.data.astype(float), h.header
    else:
        data, header = path_or_hdu
    if nan_value is not None:
        data = np.where(data == nan_value, np.nan, data)
    w, shp = (EQ_WCS, EQ_SHAPE) if which == 'eq' else (GAL_WCS, GAL_SHAPE)
    win = WCS(header).celestial
    if method == 'adaptive':
        out, _ = reproject_adaptive((data, win), w, shape_out=shp, conserve_flux=False, kernel='gaussian',
                                    boundary_mode='ignore')
    else:
        out, _ = reproject_interp((data, win), w, shape_out=shp, order=order)
    return out.astype(np.float32)


def save_grid(arr, name, which='eq', bunit='', note=''):
    w = EQ_WCS if which == 'eq' else GAL_WCS
    hdr = w.to_header()
    hdr['BUNIT'] = bunit
    hdr['COMMENT'] = note[:68]
    out = DER / 'imaging_mw_grids'
    out.mkdir(exist_ok=True)
    p = out / f'imaging_mw_{which}_{name}.fits'
    fits.PrimaryHDU(arr.astype(np.float32), header=hdr).writeto(p, overwrite=True)
    return p


def load_grid(name, which='eq'):
    p = DER / 'imaging_mw_grids' / f'imaging_mw_{which}_{name}.fits'
    return fits.getdata(p).astype(float)


# ------------------------------------------------------------------------------------------------
# Trace handling

def get_spine():
    """Return (ra, dec, meta) arrays for the current filament spine, densified to 0.05 deg."""
    tr = load_trace()
    meta = {'source': 'silk_trace.json', 'version': None, 'preliminary': None}
    if tr is not None and tr.get('spine'):
        ra = np.array([p['ra'] for p in tr['spine']]); dec = np.array([p['dec'] for p in tr['spine']])
        meta['version'] = tr.get('version'); meta['preliminary'] = tr.get('preliminary')
        meta['half_width_deg'] = tr.get('half_width_deg')
        meta['raw'] = tr
    else:
        dk = np.array([REGIONS['P3'][1], REGIONS['P1'][1], REGIONS['P2'][1]])
        rk = np.array([REGIONS['P3'][0], REGIONS['P1'][0], REGIONS['P2'][0]])
        dec = np.linspace(11.5, 21.5, 21); ra = np.polyval(np.polyfit(dk, rk, 2), dec)
        meta['source'] = 'fallback quadratic P3-P1-P2'
    order = np.argsort(dec)
    ra, dec = ra[order], dec[order]
    # densify along arc length
    c = SkyCoord(ra * u.deg, dec * u.deg)
    seg = c[:-1].separation(c[1:]).deg
    s = np.concatenate([[0], np.cumsum(seg)])
    sd = np.arange(0, s[-1] + 1e-9, 0.05)
    return np.interp(sd, s, ra), np.interp(sd, s, dec), meta


def trace_coords(lon, lat, which='eq', spine=None):
    """Signed perpendicular offset (deg; + = east/increasing RA) and along-spine coordinate (deg from south end)
    for every grid pixel, using a local tangent-plane approximation around the spine."""
    ra_s, dec_s, _ = spine if spine is not None else get_spine()
    if which == 'eq':
        ra, dec = lon, lat
    else:
        c = SkyCoord(l=lon.ravel() * u.deg, b=lat.ravel() * u.deg, frame='galactic').icrs
        ra, dec = c.ra.deg.reshape(lon.shape), c.dec.deg.reshape(lon.shape)
    # project everything to a tangent plane at the trace centre
    ra0, dec0 = np.mean(ra_s), np.mean(dec_s)
    def tp(r, d):
        r, d = np.radians(r), np.radians(d); r0, d0 = np.radians(ra0), np.radians(dec0)
        cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
        x = np.cos(d) * np.sin(r - r0) / cosc
        y = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc
        return np.degrees(x), np.degrees(y)
    xs, ys = tp(ra_s, dec_s)
    xp, yp = tp(ra, dec)
    seglen = np.hypot(np.diff(xs), np.diff(ys))
    s_cum = np.concatenate([[0], np.cumsum(seglen)])
    best_d = np.full(xp.shape, np.inf); best_off = np.zeros(xp.shape); best_s = np.zeros(xp.shape)
    for i in range(len(xs) - 1):
        ax, ay, bx, by = xs[i], ys[i], xs[i + 1], ys[i + 1]
        dx, dy = bx - ax, by - ay; L2 = dx * dx + dy * dy
        t = np.clip(((xp - ax) * dx + (yp - ay) * dy) / L2, 0, 1)
        px, py = ax + t * dx, ay + t * dy
        d = np.hypot(xp - px, yp - py)
        # sign: + if pixel is to the EAST. In this tangent plane x increases with RA (east), so for a
        # south->north segment cross > 0 means the pixel is at smaller x = WEST  -> sign -1.
        # (checked: P4 at RA 115.20, west of the spine at RA ~116.6, must get a negative offset)
        cross = dx * (yp - ay) - dy * (xp - ax)
        sgn = np.where(cross > 0, -1.0, 1.0)
        upd = d < best_d
        best_d[upd] = d[upd]; best_off[upd] = (sgn * d)[upd]; best_s[upd] = (s_cum[i] + t * np.sqrt(L2))[upd]
    # beyond the ends the "offset" is a distance from the end point; flag those
    beyond = (best_s <= 1e-6) | (best_s >= s_cum[-1] - 1e-6)
    return best_off, best_s, beyond, s_cum[-1]


def spine_s_of(ra, dec, spine=None):
    """Along-spine coordinate of a given sky position (for P1-P3)."""
    off, s, _, _ = trace_coords(np.array([[ra]]), np.array([[dec]]), 'eq', spine)
    return float(s[0, 0]), float(off[0, 0])


def binned_profile(img, off, s, s_lo, s_hi, edges, mask=None, stat='median', beyond=None):
    """Median (or mean) of img in bins of perpendicular offset for pixels with s_lo <= s <= s_hi.
    Returns centres, value, error (1.2533*sigma/sqrt(N) for the median, robust sigma via MAD), N."""
    sel = np.isfinite(img) & (s >= s_lo) & (s <= s_hi)
    if mask is not None:
        sel &= ~mask
    if beyond is not None:
        sel &= ~beyond
    cen = 0.5 * (edges[1:] + edges[:-1])
    val = np.full(len(cen), np.nan); err = np.full(len(cen), np.nan); n = np.zeros(len(cen), int)
    for i in range(len(cen)):
        m = sel & (off >= edges[i]) & (off < edges[i + 1])
        x = img[m]
        n[i] = x.size
        if x.size < 3:
            continue
        if stat == 'median':
            val[i] = np.median(x)
            sig = 1.4826 * np.median(np.abs(x - val[i]))
            err[i] = 1.2533 * sig / np.sqrt(x.size)
        else:
            val[i] = np.mean(x); err[i] = np.std(x) / np.sqrt(x.size)
    return cen, val, err, n


def disk_stat(img, lon, lat, lon0, lat0, radius, which='eq'):
    """Median and robust sigma of img within `radius` deg of (lon0, lat0) on the grid."""
    c0 = SkyCoord(lon0 * u.deg, lat0 * u.deg, frame='icrs' if which == 'eq' else 'galactic')
    c = SkyCoord(lon * u.deg, lat * u.deg, frame='icrs' if which == 'eq' else 'galactic')
    m = (c.separation(c0).deg < radius) & np.isfinite(img)
    x = img[m]
    if x.size == 0:
        return np.nan, np.nan, 0
    med = np.median(x)
    return med, 1.4826 * np.median(np.abs(x - med)), x.size


# ------------------------------------------------------------------------------------------------
# Plot helpers (use the shared silk_style)

def overlay_eq(ax, spine=None, regions=True, color='#ffffff', lw=1.0, halfwidth=True, label_regions=True,
               lamost=True, ls='-'):
    from silk_style import label, FGDIM
    ra_s, dec_s, meta = spine if spine is not None else get_spine()
    ax.plot(ra_s, dec_s, color=color, lw=lw, ls=ls, alpha=0.9, zorder=5)
    if halfwidth and meta.get('half_width_deg'):
        hw = meta['half_width_deg'] / np.cos(np.radians(dec_s))
        ax.plot(ra_s + hw, dec_s, color=color, lw=0.5, ls=':', alpha=0.6, zorder=5)
        ax.plot(ra_s - hw, dec_s, color=color, lw=0.5, ls=':', alpha=0.6, zorder=5)
    tr = meta.get('raw') or {}
    for sk in ['W', 'W2']:
        pts = (tr.get('strands') or {}).get(sk, {}).get('points') or []
        if pts:
            ax.plot([p['ra'] for p in pts], [p['dec'] for p in pts], color=color, lw=lw * 0.7, ls='--', alpha=0.8, zorder=5)
    if tr.get('halpha_patch'):
        hp_ = np.array(tr['halpha_patch'] + [tr['halpha_patch'][0]])
        ax.plot(hp_[:, 0], hp_[:, 1], color='#ff6f6f', lw=0.8, zorder=5)
    if regions:
        for k, (r, d) in REGIONS.items():
            mk = 'x' if k == 'P4' else 'o'
            ax.plot(r, d, mk, ms=4, mfc='none', mec='#ffd54f', mew=1.0, zorder=6)
            if label_regions:
                label(ax, r - 0.25, d, k, color='#ffd54f', size=7, va='center', ha='left', zorder=7)
    if lamost:
        t = np.linspace(0, 2 * np.pi, 200)
        ra0, dec0, rr = 116.3, 16.3, 2.5
        tr = meta.get('raw') or {}
        if tr.get('lamost_field'):
            ra0, dec0, rr = tr['lamost_field']['ra'], tr['lamost_field']['dec'], tr['lamost_field']['radius_deg']
        ax.plot(ra0 + rr * np.cos(t) / np.cos(np.radians(dec0 + rr * np.sin(t))), dec0 + rr * np.sin(t),
                color=FGDIM, lw=0.6, ls='--', zorder=4)


def spine_gal(spine=None):
    ra_s, dec_s, _ = spine if spine is not None else get_spine()
    l, b = gal(ra_s, dec_s)
    return l, b


def eq_axes(ax, title=None):
    ax.set_xlim(EQ_EXTENT[0], EQ_EXTENT[1]); ax.set_ylim(EQ_EXTENT[2], EQ_EXTENT[3])
    ax.set_aspect(1 / np.cos(np.radians(17.0)))
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    if title:
        ax.set_title(title, loc='left', fontsize=9)


def gal_axes(ax, title=None):
    ax.set_xlim(GAL_EXTENT[0], GAL_EXTENT[1]); ax.set_ylim(GAL_EXTENT[2], GAL_EXTENT[3])
    ax.set_aspect(1 / np.cos(np.radians(17.5)))
    ax.set_xlabel('l (deg)'); ax.set_ylabel('b (deg)')
    if title:
        ax.set_title(title, loc='left', fontsize=9)


def stretch(img, lo=1, hi=99.5, sel=None):
    x = img[np.isfinite(img)] if sel is None else img[sel & np.isfinite(img)]
    return np.nanpercentile(x, lo), np.nanpercentile(x, hi)
