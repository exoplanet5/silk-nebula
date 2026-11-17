"""Reproject the NVSS / TGSS core cutouts (1' TAN, 6x10 deg) onto the EQ grid (3', flux-averaging adaptive)."""
import sys
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import RAW, reproject_to, save_grid
for tag, name, note in [('NVSS', 'nvss', 'NVSS 1.4 GHz (interferometer; scales > ~15 arcmin resolved out)'),
                        ('TGSS_ADR1', 'tgss', 'TGSS ADR1 150 MHz (interferometer)')]:
    h = fits.open(RAW / 'skyview_imaging_mw' / f'core_{tag}.fits')
    d = h[0].data.astype(float); d[d <= -32767] = np.nan
    print(tag, 'native rms (MAD) Jy/beam', 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d))), 'median', np.nanmedian(d))
    a = reproject_to((d, h[0].header), 'eq', method='adaptive')
    save_grid(a, name, 'eq', 'Jy/beam', note)
