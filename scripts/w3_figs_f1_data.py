"""F1 inputs: a wide NSNS DR0.2 [O III] mosaic (CAR, 1 arcmin) and a smoothed eRASS1 0.2-0.6 keV grid.

Inputs (all on disk, no downloads):
  data/raw/nsns_hips/oiii/Norder3/Npix*.fits   NSNS DR0.2 [O III] HiPS tiles, Rayleigh (fetched by imaging_optical)
  data/derived/verify_xrayfuv_erass_pixels.npz  eRASS1 DR1 0.2-0.4 + 0.4-0.6 keV rate pixels (verify_xrayfuv), with
                                                the verify_xrayfuv point-source mask (catmask) and validity flag (good)
Output: data/derived/w3_figs_f1_maps.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_f1_data.py
"""
import sys
import numpy as np
from astropy.wcs import WCS
from scipy import ndimage
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import DER
from imaging_optical_common import hips_mosaic

RA0, RA1, DE0, DE1, PIX = 107.5, 124.5, 7.5, 30.5, 1.0 / 60   # deg


def car_wcs():
    nx = int(round((RA1 - RA0) / PIX)); ny = int(round((DE1 - DE0) / PIX))
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    w.wcs.crval = [0.5 * (RA0 + RA1), 0.0]
    w.wcs.cdelt = [-PIX, PIX]
    w.wcs.crpix = [(nx + 1) / 2, 1 - DE0 / PIX + 0.5]
    return w, (ny, nx)


def main():
    w, shp = car_wcs()
    yy, xx = np.mgrid[0:shp[0], 0:shp[1]]
    ra, de = w.wcs_pix2world(xx[0], np.zeros(shp[1]), 0)[0], w.wcs_pix2world(np.zeros(shp[0]), yy[:, 0], 0)[1]
    print('grid RA %.3f..%.3f Dec %.3f..%.3f' % (ra.max(), ra.min(), de.min(), de.max()))
    o3 = hips_mosaic('oiii', 3, w, shp)
    print('NSNS [O III] finite fraction %.3f, median %.3f R' % (np.isfinite(o3).mean(), np.nanmedian(o3)))
    # eRASS1 0.2-0.6 keV on a 0.1 deg grid, point sources masked, normalised-convolution smoothing
    e = np.load(DER / 'verify_xrayfuv_erass_pixels.npz')
    ok = e['good'] & ~e['catmask']
    val = (e['b0204'] + e['b0406']).astype(float)
    ok &= np.isfinite(val)
    step = 0.1
    re = np.arange(RA0 - 3, RA1 + 3 + 1e-6, step); dd = np.arange(DE0 - 3, DE1 + 3 + 1e-6, step)
    S, _, _ = np.histogram2d(e['dec'][ok], e['ra'][ok], bins=[dd, re], weights=val[ok])
    N, _, _ = np.histogram2d(e['dec'][ok], e['ra'][ok], bins=[dd, re])
    sig = 0.35 / step
    xs = ndimage.gaussian_filter(S, sig) / np.maximum(ndimage.gaussian_filter(N, sig), 1e-6)
    xs[ndimage.gaussian_filter(N, sig) < 0.2] = np.nan
    np.savez_compressed(DER / 'w3_figs_f1_maps.npz', oiii=o3.astype(np.float32), ra=ra, dec=de,
                        x_ra=0.5 * (re[1:] + re[:-1]), x_dec=0.5 * (dd[1:] + dd[:-1]), xray=xs.astype(np.float32),
                        note=np.array('NSNS DR0.2 [O III] (R) on CAR 1 arcmin; eRASS1 0.2-0.6 keV cts/s/deg2 '
                                      '(verify_xrayfuv pixels, catmask applied), 0.1 deg bins smoothed sigma 0.35 deg'))
    print('saved', DER / 'w3_figs_f1_maps.npz', 'xray range', np.nanpercentile(xs, [5, 50, 95]))


if __name__ == '__main__':
    main()
