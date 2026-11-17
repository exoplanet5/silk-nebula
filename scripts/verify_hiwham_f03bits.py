"""verify_hiwham: which surveys feed the Finkbeiner (2003) H-alpha composite along the Silk spine?

Bit definitions are taken from the FITS header of Halpha_mask_fwhm06_1024.fits (bit 0 WHAM, 1 VTSS, 2 SHASSA,
3 STAR, 4 SATUR, 5 BRT_OBJ, 6 BIG_OBJ, 7 HI_VEL). Samples the bitmask at every spine point and in a +-3 deg
band (0.05 deg steps across the spine), and reports the Dec ranges covered by each survey.
Output: data/derived/verify_hiwham_f03bits.json
"""
import json
import sys

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import RAW, DER, strand, Trace, tp_inv

mk = fits.getdata(RAW / 'lambda_finkbeiner_halpha/Halpha_mask_fwhm06_1024.fits', 1)['MASK'].astype(int).ravel()
NAMES = ['WHAM', 'VTSS', 'SHASSA', 'STAR', 'SATUR', 'BRT_OBJ', 'BIG_OBJ', 'HI_VEL']
ra, dec = strand('E')
c = SkyCoord(ra * u.deg, dec * u.deg).galactic
bits = mk[hp.ang2pix(1024, c.l.deg, c.b.deg, lonlat=True, nest=True)]
out = {'n_spine_points': int(len(ra)), 'spine': {}}
for k, nm in enumerate(NAMES):
    on = ((bits >> k) & 1) == 1
    out['spine'][nm] = {'fraction': float(on.mean()),
                        'dec_min': float(dec[on].min()) if on.any() else None,
                        'dec_max': float(dec[on].max()) if on.any() else None}
# band +-3 deg around the spine, densely resampled along the spine
T = Trace(ra, dec)
band = {}
for nm_k, nm in enumerate(NAMES[:3]):
    band[nm] = []
for dd in np.arange(-3.0, 3.01, 0.25):
    # point at signed offset dd along the local normal (d > 0 = east/convex side, as in Trace.coords)
    r2, d2 = tp_inv(T.x + dd * T.ty, T.y - dd * T.tx, T.rac, T.decc)
    chk, _, _ = T.coords(r2, d2, dmax=4.0)
    assert np.nanmedian(np.abs(chk - dd)) < 0.05, (dd, np.nanmedian(chk))
    cg = SkyCoord(r2 * u.deg, d2 * u.deg).galactic
    b2 = mk[hp.ang2pix(1024, cg.l.deg, cg.b.deg, lonlat=True, nest=True)]
    for k, nm in enumerate(NAMES[:3]):
        on = ((b2 >> k) & 1) == 1
        band[nm].append({'d': float(dd), 'fraction': float(on.mean()),
                         'dec_max': float(d2[on].max()) if on.any() else None})
out['band'] = band
(DER / 'verify_hiwham_f03bits.json').write_text(json.dumps(out, indent=1))
print(json.dumps(out['spine'], indent=1))
for nm in NAMES[:3]:
    print(nm, [(b['d'], round(b['fraction'], 2), None if b['dec_max'] is None else round(b['dec_max'], 2)) for b in band[nm]])
