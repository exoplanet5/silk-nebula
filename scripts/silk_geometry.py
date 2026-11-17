"""Key published coordinates for the Silk Nebula (provisional geometry shared by all scripts).

Sources: Wu et al. 2026, arXiv:2610.13329 (AJ, doi:10.3847/1538-3881/aeb11c), Sections I-III;
QHYCCD/SWONS press note 2024-07-12; Ziegenbalg NSNS page simg.de/nebulae3/gem-01.html.
A refined filament trace, if derived, lives in data/derived/silk_trace.json (see load_trace()).
"""
import json
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')

CENTER = (116.68, 16.42)            # RA, Dec deg (J2000), paper Sect. I; (l, b) = (203.96, +19.44)
REGIONS = {                          # paper Table 1
    'P1': (116.65, 16.60),           # bright [OIII]
    'P2': (116.21, 17.74),           # bright [OIII]
    'P3': (116.45, 14.72),           # bright [OIII]
    'P4': (115.20, 16.18),           # nebula-free background reference
}
LAMOST_FIELD = {'center_star': 'HIP 37797', 'date': '2024-02-01', 'diameter_deg': 5.0,
                'exposures': '3 x 900 s', 'n_spectra': 2684, 'n_oiii': 279}
NSNS_ARC_CENTER = (115.45, 19.933)  # 07:41:48 +19:56, ~17 deg [OIII] arc (Ziegenbalg NSNS)
EXTENT = {'paper_oiii_deg': (1.0, 10.0), 'press_oiii_deg': (1.5, 15.0), 'paper_halpha_north_deg': (0.5, 1.0)}
SEARCH_BOX = {'ra': (111.0, 122.0), 'dec': (6.0, 28.0)}   # generous provisional box

# Published line velocities (km/s) and widths, paper Tables 2 & 4
V_HELIO = {'OIII': (23.8, 26.8), 'Ha': (16.3, 19.7), 'NII': (11.7, 13.3), 'SII': (11.7, 12.7)}
V_LSR = {'OIII': (11.2, 15.5), 'Ha': (4.0, 6.1), 'NII': (-0.9, 1.0), 'SII': (-1.0, 1.7)}
FWHM_INTRINSIC = {'OIII': (43.1, 48.9), 'Ha': (31.0, 32.4), 'NII': (26.1, 30.0), 'SII': (29.0, 32.9)}

# Nearby objects frequently invoked
MONOGEM_PULSAR = {'name': 'PSR B0656+14', 'radec': (104.9506, 14.2393), 'dist_pc': 288}
GEMINGA = {'name': 'Geminga (PSR J0633+1746)', 'radec': (98.4756, 17.7703), 'dist_pc': 250}


def gal(ra, dec):
    c = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg).galactic
    return c.l.deg, c.b.deg


def load_trace():
    """Return the refined trace dict if data/derived/silk_trace.json exists, else None."""
    p = ROOT / 'data/derived/silk_trace.json'
    return json.loads(p.read_text()) if p.exists() else None


if __name__ == '__main__':
    for k, (ra, dec) in {'CENTER': CENTER, **REGIONS}.items():
        l, b = gal(ra, dec)
        print(f'{k:7s} RA {ra:8.3f} Dec {dec:+7.3f}  l {l[0]:8.3f} b {b[0]:+7.3f}')
