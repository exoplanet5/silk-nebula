"""Shared helpers for the spectra_sdss_lamost analysis: line lists, unit conversion, solar template,
linear line+continuum fitting.

All wavelengths VACUUM (SDSS, LAMOST and the Kurucz/CALSPEC solar model are vacuum).
Vacuum rest wavelengths from the NIST ASD (Ritz) values converted with the IAU/Morton (2000) formula.
"""
import numpy as np
from pathlib import Path
from astropy.io import fits

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
C_KMS = 299792.458

# air -> vacuum (Morton 2000 / IAU), sufficient to 0.01 A in the optical
def air_to_vac(lam_air):
    lam_air = np.asarray(lam_air, float)
    s2 = (1e4 / lam_air) ** 2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    return lam_air * n


# air wavelengths (A) of the nebular lines (standard values, e.g. Osterbrock & Ferland 2006)
AIR = {'Hb': 4861.33, 'OIII4959': 4958.91, 'OIII5007': 5006.84, 'NI5200': 5199.0,
       'NII6548': 6548.05, 'Ha': 6562.80, 'NII6583': 6583.45, 'SII6716': 6716.44, 'SII6731': 6730.82,
       'OI6300': 6300.30, 'HeI5876': 5875.62}
VAC = {k: float(air_to_vac(v)) for k, v in AIR.items()}

# fibre diameters (arcsec)
FIBRE_DIAM = {'SDSS': 3.0, 'BOSS': 2.0, 'LAMOST': 3.3}
ARCSEC2_SR = (np.pi / 180 / 3600) ** 2


def rayleigh_per_fibre_flux(lam_A, diam_arcsec):
    """Flux per fibre (erg/s/cm^2) produced by 1 Rayleigh of a line at lam_A filling the fibre."""
    omega = np.pi * (diam_arcsec / 2) ** 2 * ARCSEC2_SR
    e_ph = 6.62607015e-27 * 2.99792458e18 / lam_A
    return 1e6 / (4 * np.pi) * omega * e_ph


def flux_to_rayleigh(f_1e17, lam_A, diam_arcsec):
    return np.asarray(f_1e17) * 1e-17 / rayleigh_per_fibre_flux(lam_A, diam_arcsec)


# ---------------------------------------------------------------- solar template
_SUN = None


def load_sun():
    """Kurucz solar model (CALSPEC sun_mod_001.fits, R=5000, vacuum) -> (lam, flux)."""
    global _SUN
    if _SUN is None:
        d = fits.open(ROOT / 'data/raw/calspec_sun/sun_mod_001.fits')[1].data
        m = (d['WAVELENGTH'] > 3500) & (d['WAVELENGTH'] < 10500)
        _SUN = (np.array(d['WAVELENGTH'][m], float), np.array(d['FLUX'][m], float))
    return _SUN


def sun_at(lam, sigma_kms, v_kms=0.0, r_model=5000.0):
    """Solar spectrum broadened to a Gaussian LSF of sigma_kms (total) and evaluated at lam (A)."""
    ls, fs = load_sun()
    sig_model = C_KMS / r_model / 2.3548
    sig = np.sqrt(max(sigma_kms ** 2 - sig_model ** 2, 1.0))
    dv = 5.0                                     # km/s grid
    lo, hi = np.log(lam.min() * 0.995), np.log(lam.max() * 1.005)
    grid = np.exp(np.arange(lo, hi, dv / C_KMS))
    f = np.interp(grid, ls, fs)
    k = np.arange(-int(5 * sig / dv) - 1, int(5 * sig / dv) + 2) * dv
    ker = np.exp(-0.5 * (k / sig) ** 2)
    ker /= ker.sum()
    fc = np.convolve(f, ker, mode='same')
    return np.interp(lam / (1 + v_kms / C_KMS), grid, fc)


def gauss(lam, cen, sig):
    return np.exp(-0.5 * ((lam - cen) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)


def linfit(y, ivar, cols):
    """Weighted linear least squares; returns (coef, cov, chi2, model)."""
    A = np.vstack(cols).T
    w = np.sqrt(np.clip(ivar, 0, None))
    Aw = A * w[:, None]
    yw = y * w
    coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.pinv(Aw.T @ Aw)
    model = A @ coef
    chi2 = float(np.sum(((y - model) ** 2) * ivar))
    return coef, cov, chi2, model
