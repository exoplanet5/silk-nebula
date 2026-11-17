"""verify_precovery step 1: per-field templates, geometry, masks and observing conditions (cached).

For every SDSS DR18 field in the box (fields_box_full.csv) compute, from my own code:
  * footprint-averaged (mean and median of 9x7 samples over the 13.5' x 9.8' frame, oriented along the scan)
    values of the ridge templates: gauss028 (E+W, FWHM 0.28 deg, unit), gauss028_EWW2, gaussE028,
    gaussfw0.60, gaussfw1.00 (broad tests), desi (DESI forward-model amplitudes, R), fibre (SDSS sky-fibre
    strand amplitudes, R), each also shifted in RA by s (deg, at fixed Dec) for null / injection tests;
  * signed distance of the field centre to the E spine and to the nearest strand (real and shifted);
  * bright-star mask (Tycho-2 VT<7 from data/raw/vizier_tycho2 + 5 bright stars missing from that query,
    coordinates J2000 from the Yale BSC values typed here), planet proximity (astropy builtin ephemeris);
  * Moon altitude / separation / illumination at mjd_g (APO), ecliptic latitude, elongation from the Sun,
    IRIS 100 um (imaging_mw grid) and Finkbeiner Halpha at field centre.
Output: data/derived/verify_precovery_fields.npz
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body, get_sun, GeocentricTrueEcliptic
import astropy.units as u
from verify_precovery_lib import (load_fields, load_strands, make_model, scan_frames, footprint_points, DER, RAW)

SHIFTS = [round(s, 2) for s in np.r_[np.arange(-6.0, -0.75, 0.2), np.arange(0.8, 6.01, 0.2)]]
TEMPLATES = ['gauss028', 'gauss028_EWW2', 'gaussE028', 'gaussfw0.60', 'gaussfw1.00', 'desi', 'fibre']
SHIFT_TEMPLATES = ['gauss028', 'desi']

BRIGHT_EXTRA = [  # name, RA, Dec (J2000, deg), V  -- stars absent from (or without positions in) the Tycho-2 query
    ('Pollux', 116.3290, 28.0262, 1.14), ('Procyon', 114.8255, 5.2250, 0.34), ('Castor', 113.6496, 31.8883, 1.58),
    ('Gomeisa', 111.7877, 8.2893, 2.89), ('delta Gem', 110.0307, 21.9823, 3.53)]


def tycho_bright(vmax=7.0):
    rows = [l.rstrip('\n').split('\t') for l in open(RAW / 'vizier_tycho2/tyc2_VT11_box.tsv') if l.strip() and not l.startswith('#')]
    out = []
    for r in rows[3:]:
        try:
            ra, de, vt = float(r[0]), float(r[1]), float(r[3])
        except ValueError:
            continue
        if vt < vmax:
            out.append((ra, de, vt))
    out += [(r, d, v) for _, r, d, v in BRIGHT_EXTRA]
    return np.array(out)


def star_mask(a, stars):
    from scipy.spatial import cKDTree
    from verify_precovery_lib import uvec
    t = cKDTree(uvec(a['ra'], a['dec']))
    bad = np.zeros(len(a), bool); near = np.full(len(a), 99.0)
    for ra, de, v in stars:
        r = 1.0 if v < 1.5 else 0.5 if v < 3 else 0.25 if v < 5 else 0.12
        idx = t.query_ball_point(uvec(ra, de), np.radians(r))
        bad[idx] = True
    return bad


def sample_grid(path, ra, dec):
    with fits.open(path) as h:
        w = WCS(h[0].header); d = h[0].data.astype(float)
    x, y = w.world_to_pixel_values(ra, dec)
    ix = np.round(x).astype(int); iy = np.round(y).astype(int)
    ok = (ix >= 0) & (iy >= 0) & (ix < d.shape[1]) & (iy < d.shape[0])
    out = np.full(len(ra), np.nan); out[ok] = d[iy[ok], ix[ok]]
    return out


def main():
    a = load_fields()
    N = len(a)
    strands, trace = load_strands()
    frames = scan_frames(a)
    ra_fp, dec_fp = footprint_points(a, frames=frames)
    out = dict(run=a['run'].astype(int), camcol=a['camcol'].astype(int), field=a['field'].astype(int),
               ra=a['ra'], dec=a['dec'], mjd=a['mjd_g'], quality=a['quality'].astype(int))
    # geometry: distances to E spine and to nearest strand (E/W)
    dE, pdE, oE, _ = strands['E'].geometry(a['ra'], a['dec'])
    dW, pdW, oW, _ = strands['W'].geometry(a['ra'], a['dec'])
    out.update(dE=dE, decE=pdE, overE=oE, dW=dW, overW=oW)
    # templates
    for k in TEMPLATES:
        m = make_model(k, strands)
        v = m(ra_fp, dec_fp)
        out[f'T_{k}_mean'] = v.mean(1); out[f'T_{k}_med'] = np.median(v, 1)
        out[f'T_{k}_centre'] = m(a['ra'], a['dec'])
        print('template', k, 'max', np.nanmax(out[f'T_{k}_mean']))
    for k in SHIFT_TEMPLATES:
        arr = np.zeros((len(SHIFTS), N))
        for i, s in enumerate(SHIFTS):
            m = make_model(k, strands, ra_shift=s)
            arr[i] = m(ra_fp, dec_fp).mean(1)
        out[f'TS_{k}'] = arr
    dEs = np.zeros((len(SHIFTS), N)); oEs = np.zeros((len(SHIFTS), N))
    for i, s in enumerate(SHIFTS):
        d, _, o, _ = strands['E'].geometry(a['ra'] - s, a['dec'])
        dEs[i] = d; oEs[i] = o
    out.update(shifts=np.array(SHIFTS), dE_shift=dEs, overE_shift=oEs)
    # cross-ridge profile basis: fraction of footprint samples in bins of signed distance to the E spine
    dfp, _, ofp, _ = strands['E'].geometry(ra_fp, dec_fp)
    edges = np.arange(-2.5, 2.5001, 0.2)
    prof = np.zeros((len(edges) - 1, N))
    for j in range(len(edges) - 1):
        prof[j] = ((dfp >= edges[j]) & (dfp < edges[j + 1])).mean(1)
    out.update(prof_edges=edges, prof_basis=prof)
    # stars
    stars = tycho_bright()
    out['star_bad'] = star_mask(a, stars)
    print('bright-star masked fields', out['star_bad'].sum(), 'of', N, 'stars', len(stars))
    # observing conditions per field at mjd_g (TAI ~ UTC to 1 min; irrelevant here)
    apo = EarthLocation(lat=32.7803 * u.deg, lon=-105.8203 * u.deg, height=2788 * u.m)
    t = Time(a['mjd_g'], format='mjd', scale='tai')
    c = SkyCoord(a['ra'] * u.deg, a['dec'] * u.deg)
    moon = get_body('moon', t, apo)
    sun = get_body('sun', t, apo)
    aa = AltAz(obstime=t, location=apo)
    out['moon_alt'] = moon.transform_to(aa).alt.deg
    out['moon_sep'] = c.separation(moon, origin_mismatch='ignore').deg if hasattr(c, 'separation') else np.nan
    elong = sun.separation(moon, origin_mismatch='ignore').rad
    out['moon_illum'] = (1 - np.cos(elong)) / 2
    out['sun_alt'] = sun.transform_to(aa).alt.deg
    ecl = c.transform_to(GeocentricTrueEcliptic(equinox=t))
    secl = sun.transform_to(GeocentricTrueEcliptic(equinox=t))
    out['ecl_beta'] = ecl.lat.deg
    out['ecl_dlam'] = ((ecl.lon.deg - secl.lon.deg + 180) % 360) - 180    # 0 = Sun direction; +-180 antisolar
    out['alt'] = c.transform_to(aa).alt.deg
    # planets: minimum separation from field centre
    psep = np.full(N, 99.0)
    for p in ('venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune'):
        b = get_body(p, t, apo)
        psep = np.minimum(psep, c.separation(b, origin_mismatch='ignore').deg)
    out['planet_sep'] = psep
    # maps
    grid = DER / 'imaging_mw_grids'
    out['iris100'] = sample_grid(grid / 'imaging_mw_eq_iris100.fits', a['ra'], a['dec'])
    out['halpha_f03'] = sample_grid(grid / 'imaging_mw_eq_halpha_f03.fits', a['ra'], a['dec'])
    # sky estimators and other columns
    for col in a.dtype.names:
        if col.startswith('sky') or col.startswith('airmass') or col.startswith('sigPix') or col.startswith('nMgy') \
                or col in ('psfWidth_g', 'nStars', 'nTotal', 'nBrightObj_g', 'score'):
            out['f_' + col] = a[col]
    np.savez_compressed(DER / 'verify_precovery_fields.npz', **out)
    print('planet within 1 deg of a field:', int((psep < 1).sum()), ' moon up (alt>0) fields:', int((out['moon_alt'] > 0).sum()))


if __name__ == '__main__':
    main()
