"""verify_distance: shared readers for the UVES / HARPS / SOPHIE / ELODIE products on disk.

Independent of the wave-1 spectra_ism_* code (no shared functions).  All inputs are
untrusted FITS files read with astropy (no pickles).  Run scripts with `python -I`.
"""
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
SPEC = ROOT / 'data/raw/spectra_ism_spectra'
C = 299792.458

# air rest wavelengths (Angstrom); f-values (Morton 2003)
LINES = {
    'CaK': (3933.6614, 0.6267), 'CaH': (3968.4673, 0.3116),
    'NaD2': (5889.9510, 0.6408), 'NaD1': (5895.9242, 0.3201),
    'TiII3384': (3383.7588, 0.358), 'CaI4227': (4226.7280, 1.77),
    'CHp4232': (4232.548, 0.00545), 'CH4300': (4300.3132, 0.00506),
    'CN3874': (3874.608, 0.0342), 'FeI3860': (3859.9114, 0.0217),
    'NaI3302': (3302.369, 0.00899),
}


def read_uves(path):
    """Return dict(wave, flux, err, header) for an ESO phase-3 1D spectrum."""
    with fits.open(path, memmap=False) as h:
        H = h[0].header.copy()
        d = h[1].data
        names = d.columns.names
        w = np.array(d['WAVE'][0], float)
        fk = 'FLUX' if 'FLUX' in names else names[1]
        ek = 'ERR' if 'ERR' in names else names[2]
        f = np.array(d[fk][0], float)
        e = np.array(d[ek][0], float)
        fr = np.array(d['FLUX_REDUCED'][0], float) if 'FLUX_REDUCED' in names else f
        er = np.array(d['ERR_REDUCED'][0], float) if 'ERR_REDUCED' in names else e
    return dict(wave=w, flux=f, err=e, flux_red=fr, err_red=er, header=H,
                specsys=H.get('SPECSYS'), berv=H.get('HIERARCH ESO QC VRAD BARYCOR'),
                mjd=H.get('MJD-OBS'), R=H.get('SPEC_RES'))


def window(sp, lam0, vlo, vhi, frame_shift=0.0, use_reduced=True):
    """Cut a window around lam0 in velocity (km/s) relative to the spectrum's frame.

    frame_shift is added to the native-frame velocity (e.g. BERV for topocentric -> bary)."""
    w = sp['wave']
    v = (w / lam0 - 1) * C + frame_shift
    m = (v > vlo) & (v < vhi) & np.isfinite(sp['flux'])
    f = sp['flux_red'][m] if use_reduced else sp['flux'][m]
    e = sp['err_red'][m] if use_reduced else sp['err'][m]
    return v[m], f, e
