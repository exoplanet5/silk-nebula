"""spectra_ism: loaders for archival 1D high-resolution spectra (ESO phase-3 UVES/HARPS, SOPHIE s1d, ELODIE s1d).

load(path) -> dict(wave_air_A, flux, err, frame, berv, R, date, inst, object, ra, dec)
  frame 'topo' : wavelengths topocentric (UVES phase-3 SPECSYS=TOPOCENT); v_helio = v_topo + berv
  frame 'bary' : wavelengths already barycentric (HARPS phase-3 SPECSYS=BARYCENT; SOPHIE s1d)
  ELODIE s1d are topocentric. Both OHP frames were verified empirically with telluric O2 gamma-band lines
  (scripts/spectra_ism_frame_test.py -> data/derived/spectra_ism_frame_test.csv).
  If an ESO product lacks a BERV keyword, it is computed with astropy for Paranal/La Silla at DATE-OBS (+ EXPTIME/2).
Only FITS data are read (astropy.io.fits); nothing in the files is executed.
"""
import numpy as np
from astropy.io import fits

C = 299792.458
LINES = {'CaK': 3933.663, 'CaH': 3968.469, 'NaD2': 5889.951, 'NaD1': 5895.924}   # air wavelengths (A)
FVAL = {'CaK': 0.6267, 'CaH': 0.3116, 'NaD2': 0.6408, 'NaD1': 0.3201}          # Morton 2003 oscillator strengths
LAM_VAC = {'CaK': 3934.777, 'CaH': 3969.591, 'NaD2': 5891.583, 'NaD1': 5897.558}


def load(path):
    h = fits.open(path)
    h0 = h[0].header
    name = str(path).split('/')[-1]
    if name.startswith('eso_'):
        d = h[1].data
        w = np.asarray(d['WAVE'][0], float)
        cols = d.columns.names
        fcol = 'FLUX_REDUCED' if 'FLUX_REDUCED' in cols else 'FLUX'
        ecol = 'ERR_REDUCED' if 'ERR_REDUCED' in cols else 'ERR'
        f = np.asarray(d[fcol][0], float)
        e = np.asarray(d[ecol][0], float)
        if 'FLUX' in cols and fcol != 'FLUX':
            f2 = np.asarray(d['FLUX'][0], float)
            if np.isfinite(f2).sum() > 0.9 * len(f2) and np.nanmedian(f2) > 0:
                pass  # keep FLUX_REDUCED (not flux-calibrated, no response ripple from calibration)
        if np.nanmedian(w) < 1000:
            w = w * 10.0
        spec = h0.get('SPECSYS', 'TOPOCENT')
        frame = 'bary' if spec.startswith('BARY') else 'topo'
        berv = h0.get('HIERARCH ESO QC VRAD BARYCOR', h0.get('HIERARCH ESO DRS BERV', np.nan))
        R = h0.get('SPEC_RES', np.nan)
        inst = h0.get('INSTRUME', '?')
    else:
        f = np.asarray(h[0].data, float)
        n = f.size
        crval, cdelt = h0['CRVAL1'], h0['CDELT1']
        crpix = h0.get('CRPIX1', 1.0)
        w = crval + (np.arange(n) + 1 - crpix) * cdelt
        if name.startswith('elodie_'):
            e = np.asarray(h['NOISE'].data, float) if 'NOISE' in h else np.full(n, np.nan)
            berv = h0.get('BERV', np.nan)
            R, inst = 42000.0, 'ELODIE'
        else:
            e = np.full(n, np.nan)
            berv = h0.get('HIERARCH OHP DRS BERV', np.nan)
            mode = h0.get('HIERARCH OHP INS FIBER MODE', '')
            R = 75000.0 if 'HR' in str(mode) or not mode else 39000.0
            inst = 'SOPHIE'
        frame = 'topo' if name.startswith('elodie_') else 'bary'
    if not np.isfinite(berv):
        berv = compute_berv(h0)
    good = np.isfinite(w) & np.isfinite(f)
    return dict(wave=w[good], flux=f[good], err=e[good] if e is not None else None, frame=frame, berv=float(berv),
                R=float(R), date=h0.get('DATE-OBS', ''), inst=inst, object=h0.get('OBJECT', ''),
                ra=h0.get('RA', np.nan), dec=h0.get('DEC', np.nan), header=h0)


def compute_berv(h0):
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy.time import Time
    import astropy.units as u
    lat, lon = h0.get('HIERARCH ESO TEL GEOLAT', -24.6253), h0.get('HIERARCH ESO TEL GEOLON', -70.4022)
    hgt = h0.get('HIERARCH ESO TEL GEOELEV', 2635.0)
    loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=hgt * u.m)
    t = Time(h0['DATE-OBS']) + 0.5 * float(h0.get('EXPTIME', 0.0)) * u.s
    c = SkyCoord(float(h0['RA']) * u.deg, float(h0['DEC']) * u.deg)
    return float(c.radial_velocity_correction('barycentric', obstime=t, location=loc).to(u.km / u.s).value)


def window(sp, line, vmin=-400, vmax=400):
    """Return heliocentric velocity, flux, err in a window around an interstellar line."""
    lam0 = LINES[line]
    v = (sp['wave'] / lam0 - 1) * C
    if sp['frame'] == 'topo':
        v = v + sp['berv']
    m = (v > vmin) & (v < vmax)
    e = sp['err'][m] if sp['err'] is not None else None
    return v[m], sp['flux'][m], e
