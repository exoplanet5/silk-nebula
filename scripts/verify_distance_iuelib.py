"""verify_distance: minimal reader for IUE NEWSIPS high-dispersion MXHI files (untrusted FITS, read-only).

NEWSIPS MXHI: one table row per echelle order; wavelength = WAVELENGTH + DELTAW*i (vacuum, heliocentric
correction applied by NEWSIPS, see HISTORY 'NET CORRECTION VECTOR TO HELIOCENTRIC VELOCITY');
ABS_CAL = ripple-corrected absolutely calibrated flux; NOISE in FN units (scaled to ABS_CAL here);
QUALITY < 0 = flagged pixel.
"""
from pathlib import Path
import numpy as np
from astropy.io import fits

C = 299792.458
# vacuum wavelengths and f-values (Morton 2003)
UVLINES = {
    'CIV1548': (1548.204, 0.1899), 'CIV1550': (1550.781, 0.09475),
    'SiIV1393': (1393.7602, 0.513), 'SiIV1402': (1402.7729, 0.254),
    'NV1238': (1238.821, 0.156), 'NV1242': (1242.804, 0.0777),
    'SiII1526': (1526.7070, 0.133), 'SiII1304': (1304.3702, 0.0928), 'SiII1260': (1260.4221, 1.18),
    'CII1334': (1334.5323, 0.1278), 'CIIs1335': (1335.7077, 0.115),
    'SII1250': (1250.578, 0.00543), 'SII1253': (1253.805, 0.0109), 'SII1259': (1259.518, 0.0166),
    'AlII1670': (1670.7886, 1.74), 'FeII1608': (1608.4511, 0.0577), 'OI1302': (1302.1685, 0.048),
    'AlIII1854': (1854.7184, 0.559), 'AlIII1862': (1862.7910, 0.278), 'SiIII1206': (1206.500, 1.63),
    'NiII1370': (1370.132, 0.0588), 'CI1656': (1656.928, 0.149),
}


def read_mxhi(path):
    with fits.open(path, memmap=False) as h:
        H = h[0].header
        d = h[1].data
        orders = []
        for row in d:
            n = int(row['NPOINTS'])
            w = row['WAVELENGTH'] + row['DELTAW'] * np.arange(n)
            sp = int(row['STARTPIX']) - 1          # STARTPIX is 1-based (checked against NET != 0)
            flux = np.array(row['ABS_CAL'][sp:sp + n], float)
            net = np.array(row['NET'][sp:sp + n], float)
            noise = np.array(row['NOISE'][sp:sp + n], float)
            q = np.array(row['QUALITY'][sp:sp + n], int)
            with np.errstate(divide='ignore', invalid='ignore'):
                scale = np.where(net != 0, flux / net, np.nan)
            err = np.abs(noise * scale)
            bad = (flux == 0) | ~np.isfinite(flux)
            flux[bad] = np.nan
            err[bad] = np.nan
            orders.append(dict(order=int(row['ORDER']), wave=w, flux=flux, err=err, q=q))
        ap = str(H.get('APERTURE', '')).strip().upper()
        pre = 'S' if ap.startswith('SMALL') else 'L'
        meta = dict(target=str(H.get(pre + 'TARGET') or H.get('LTARGET') or H.get('STARGET') or '').strip(),
                    jd=H.get(pre + 'JD-MID'), aperture=ap, radvelo=H.get(pre + 'RADVELO'),
                    camera=H.get('CAMERA'), image=H.get('IMAGE'))
    return orders, meta


def extract(orders, lam0, vwin=700.0, edge_frac=0.06):
    """Velocity-window spectra of a line from every order that contains it (away from order edges)."""
    out = []
    for o in orders:
        w = o['wave']
        if len(w) < 50:
            continue
        lo = w[0] + edge_frac * (w[-1] - w[0])
        hi = w[-1] - edge_frac * (w[-1] - w[0])
        if not (lo < lam0 < hi):
            continue
        v = (w / lam0 - 1) * C
        m = (np.abs(v) < vwin) & np.isfinite(o['flux'])
        out.append(dict(order=o['order'], v=v[m], f=o['flux'][m], e=o['err'][m], q=o['q'][m],
                        centrality=abs(lam0 - 0.5 * (w[0] + w[-1])) / (w[-1] - w[0])))
    return out
