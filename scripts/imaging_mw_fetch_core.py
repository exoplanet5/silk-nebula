"""SkyView cutouts of interferometric radio surveys (NVSS 1.4 GHz, TGSS ADR1 150 MHz) over the Silk core
(TAN centred RA 116.5, Dec 16.5, 6 x 10 deg, 1' pixels, flux-averaging Clip sampler), and an AKARI retry."""
import sys, time
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_fetch import fetch, RAW, SV
for name, tag in [('NVSS', 'NVSS'), ('TGSS ADR1', 'TGSS_ADR1'), ('AKARI WIDE-S', 'AKARI_WIDE-S')]:
    params = {'Position': '116.5,16.5', 'Survey': name, 'Coordinates': 'J2000', 'Projection': 'Tan',
              'Pixels': '360,600', 'Size': '6,10', 'Sampler': 'Clip', 'ClipIntensive': '', 'Return': 'FITS'}
    fetch(SV, RAW / 'skyview_imaging_mw' / f'core_{tag}.fits', 'skyview_imaging_mw', 'SkyView current; survey ' + name,
          params=params, timeout=1800, note='TAN 6x10 deg @ 1 arcmin, Clip intensive')
    time.sleep(3)
