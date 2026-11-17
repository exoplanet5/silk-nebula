"""Fetch Edenhofer et al. (2024) 3D-dust line-of-sight profiles toward the Silk Nebula via HTTP Range requests.

The full product mean_and_std_healpix.fits (Zenodo record 10658339, v1.0.x, 3.25 GB) is NOT downloaded.
Layout (read from its FITS headers, 2026-11-16): HDU1 'MEAN' float32 image [516 radii x 786432 NEST pixels,
nside 256], HDU2 'STD.' same shape, HDU3 radial pixel centres (pc), HDU4 radial pixel boundaries (pc).
Units: dust extinction density in "E of Zhang, Green & Rix (2023)" per pc.
For each line of sight we read the 4 contiguous nside-256 sub-pixels of the enclosing nside-128 NEST pixel
(~0.46 deg) for all 516 radii (16 bytes per radius) with multipart byte-range requests.

Output (raw fetched values, float32): data/raw/context_dust/edenhofer24_los_profiles.npz
Usage: ~/.venvs/astro313/bin/python -I scripts/context_dust_fetch.py
"""
import sys, re, json, time
from pathlib import Path
import numpy as np
import requests
import healpy as hp
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import REGIONS, CENTER, ROOT, gal

URL = 'https://zenodo.org/records/10658339/files/mean_and_std_healpix.fits?download=1'
NPIX, NRAD = 786432, 516
OFF_MEAN = 5760                       # 2 header blocks
OFF_STD = 1623202560 + 2880           # HDU2 header at 1623202560 (one block)
OFF_RC = 3246402240 + 2880            # HDU3 data (516 x float32 big-endian)
OFF_RB = 3246402240 + 2880 + 2880 + 2880  # HDU4 data (517 x float32)
OUT = ROOT / 'data/raw/context_dust/edenhofer24_los_profiles.npz'


def multirange(ranges, tries=4):
    hdr = {'Range': 'bytes=' + ','.join(f'{a}-{b}' for a, b in ranges)}
    for k in range(tries):
        try:
            r = requests.get(URL, headers=hdr, timeout=300)
            if r.status_code != 206:
                raise RuntimeError(f'status {r.status_code}')
            ct = r.headers.get('content-type', '')
            if len(ranges) == 1 and 'multipart' not in ct:
                return [r.content]
            body = r.content
            # multipart/byteranges parsing (boundary from header, or detect)
            m = re.search(rb'--([^\r\n]+)\r?\n', body)
            bnd = m.group(1)
            parts = body.split(b'--' + bnd)
            out = {}
            for p in parts:
                mm = re.search(rb'Content-Range: bytes (\d+)-(\d+)/\d+\r?\n\r?\n', p)
                if not mm:
                    continue
                a, b = int(mm.group(1)), int(mm.group(2))
                start = mm.end()
                out[a] = p[start:start + (b - a + 1)]
            return [out[a] for a, b in ranges]
        except Exception as e:
            print('retry', k, e)
            time.sleep(5 * (k + 1))
    raise RuntimeError('multirange failed')


def los_list():
    pts = {k: v for k, v in REGIONS.items()}
    pts['CENTER'] = CENTER
    # samples along the refined E-strand spine (data/derived/silk_trace.json v1.0, geometry teammate)
    from silk_geometry import load_trace
    sp = load_trace()['spine']
    sra = np.array([q['ra'] for q in sp]); sde = np.array([q['dec'] for q in sp])
    for dec in np.arange(12.0, 28.01, 1.0):
        pts[f'arc_dec{dec:04.1f}'] = (round(float(np.interp(dec, sde, sra)), 3), round(dec, 2))
    pts['east2_P1'] = (118.75, 16.60)
    pts['west2_P1'] = (114.55, 16.60)
    pts['east4_P1'] = (120.85, 16.60)
    pts['west4_P1'] = (112.45, 16.60)
    pts['Monogem_ctr'] = (107.08, 13.50)
    pts['WD_J0736+1444'] = (114.0329, 14.7475)
    return pts


if __name__ == '__main__':
    rc = np.frombuffer(multirange([(OFF_RC, OFF_RC + NRAD * 4 - 1)])[0], dtype='>f4').astype(float)
    rb = np.frombuffer(multirange([(OFF_RB, OFF_RB + (NRAD + 1) * 4 - 1)])[0], dtype='>f4').astype(float)
    print('radii', rc[:3], rc[-3:], 'bounds', rb[:2], rb[-2:])
    pts = los_list()
    names, means, stds, pix = [], [], [], []
    for name, (ra, dec) in pts.items():
        l, b = gal(ra, dec)
        p128 = hp.ang2pix(128, l[0], b[0], lonlat=True, nest=True)
        sub0 = p128 * 4                                  # 4 nested children at nside 256 are contiguous
        prof = {}
        for tag, off in (('mean', OFF_MEAN), ('std', OFF_STD)):
            vals = []
            rows = list(range(NRAD))
            for i in range(0, NRAD, 172):
                chunk = rows[i:i + 172]
                rngs = [(off + (r * NPIX + sub0) * 4, off + (r * NPIX + sub0) * 4 + 15) for r in chunk]
                vals += [np.frombuffer(x, dtype='>f4') for x in multirange(rngs)]
            prof[tag] = np.array(vals)  # (516, 4)
        names.append(name); means.append(prof['mean']); stds.append(prof['std']); pix.append(sub0)
        print(f'{name:16s} l={l[0]:7.2f} b={b[0]:+6.2f} nest128={p128}  E(<1.25kpc)~'
              f'{np.sum(prof["mean"].mean(1) * np.diff(rb)):.3f}')
    np.savez(OUT, names=np.array(names), ra=np.array([pts[n][0] for n in names]), dec=np.array([pts[n][1] for n in names]),
             mean=np.array(means, dtype='f4'), std=np.array(stds, dtype='f4'), nest256_first=np.array(pix),
             r_centers=rc, r_bounds=rb,
             meta=json.dumps({'url': URL, 'zenodo_record': 10658339, 'file': 'mean_and_std_healpix.fits',
                              'accessed': '2026-11-16', 'units': 'E (Zhang, Green & Rix 2023) per pc',
                              'pixels': 'nside256 NEST, 4 children of the nside128 pixel containing each LOS'}))
    print('saved', OUT)
