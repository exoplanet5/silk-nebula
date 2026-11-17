"""spectra_ism: angular separation from the [O III] spine (shared helper)."""
import sys
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import spine_arrays  # noqa: E402


def spine_sep(ra, dec):
    sra, sdec, _ = spine_arrays()
    t = np.linspace(0, 1, 1500)
    cs = SkyCoord(np.interp(t, np.linspace(0, 1, len(sra)), sra) * u.deg,
                  np.interp(t, np.linspace(0, 1, len(sdec)), sdec) * u.deg)
    _, sep, _ = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg).match_to_catalog_sky(cs)
    return sep.deg
