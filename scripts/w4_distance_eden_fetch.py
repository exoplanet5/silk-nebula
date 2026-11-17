"""w4_distance: read an Edenhofer et al. (2024) 3D-dust sub-cube around the Silk ridge and PGCC G199.09+19.95 by
HTTP Range requests (the 3.25 GB file is NOT downloaded).

File: Zenodo record 10658339, mean_and_std_healpix.fits (layout from scripts/context_dust_fetch.py, which read the
FITS headers on 2026-11-16): HDU1 MEAN float32 [516 radii x 786432 NEST pixels, nside 256]; HDU2 STD same;
HDU3 radial centres; HDU4 radial bounds. Units: E (Zhang, Green & Rix 2023) per pc.

Region: all nside-16 NEST parents whose centre or corners fall in Galactic l 189-211, b 13.5-27.5 (covers RA ~112-121,
Dec ~10-30). Each nside-16 parent holds 256 contiguous nside-256 children -> 1 KiB per radius.
Radii: indices with 90 <= r <= 750 pc. MEAN for all parents; STD only for the parents containing the PGCC and the
south-ridge point P1.
Output: data/raw/w4_distance_eden/eden_subcube.npz (+ fetch_log.tsv)
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_eden_fetch.py
"""
import re, time, json, datetime
from pathlib import Path
import numpy as np
import requests
import healpy as hp

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUTD = ROOT / 'data/raw/w4_distance_eden'
OUTD.mkdir(parents=True, exist_ok=True)
LOG = OUTD / 'fetch_log.tsv'
URL = 'https://zenodo.org/records/10658339/files/mean_and_std_healpix.fits?download=1'
NPIX, NRAD = 786432, 516
OFF_MEAN = 5760
OFF_STD = 1623202560 + 2880
OFF_RC = 3246402240 + 2880
OFF_RB = 3246402240 + 2880 + 2880 + 2880
NBYTES = 0


def log(msg):
    new = not LOG.exists()
    with open(LOG, 'a') as f:
        if new:
            f.write('time_utc\turl\tnote\n')
        f.write(f'{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{URL}\t{msg}\n')


def multirange(ranges, tries=5):
    global NBYTES
    hdr = {'Range': 'bytes=' + ','.join(f'{a}-{b}' for a, b in ranges)}
    for k in range(tries):
        try:
            r = requests.get(URL, headers=hdr, timeout=300)
            if r.status_code != 206:
                raise RuntimeError(f'status {r.status_code}')
            NBYTES += len(r.content)
            ct = r.headers.get('content-type', '')
            if len(ranges) == 1 and 'multipart' not in ct:
                return [r.content]
            body = r.content
            m = re.search(rb'--([^\r\n]+)\r?\n', body)
            bnd = m.group(1)
            out = {}
            for p in body.split(b'--' + bnd):
                mm = re.search(rb'Content-Range: bytes (\d+)-(\d+)/\d+\r?\n\r?\n', p)
                if not mm:
                    continue
                a, b = int(mm.group(1)), int(mm.group(2))
                out[a] = p[mm.end():mm.end() + (b - a + 1)]
            return [out[a] for a, b in ranges]
        except Exception as e:
            print('retry', k, e)
            time.sleep(5 * (k + 1))
    raise RuntimeError('multirange failed')


def parents16():
    ps = set()
    for l in np.arange(189, 211.01, 0.5):
        for b in np.arange(13.5, 27.51, 0.5):
            ps.add(int(hp.ang2pix(16, l, b, lonlat=True, nest=True)))
    return sorted(ps)


if __name__ == '__main__':
    t0 = time.time()
    rc = np.frombuffer(multirange([(OFF_RC, OFF_RC + NRAD * 4 - 1)])[0], dtype='>f4').astype(float)
    rb = np.frombuffer(multirange([(OFF_RB, OFF_RB + (NRAD + 1) * 4 - 1)])[0], dtype='>f4').astype(float)
    ridx = np.where((rc >= 90) & (rc <= 750))[0]
    par = parents16()
    pg = int(hp.ang2pix(16, 199.094, 19.947, lonlat=True, nest=True))
    p1 = int(hp.ang2pix(16, 207.07, 18.33, lonlat=True, nest=True))   # approx (l,b) of P1 region
    std_par = sorted({pg, p1})
    print(len(par), 'nside16 parents;', len(ridx), 'radii', rc[ridx[0]], rc[ridx[-1]])
    mean = np.zeros((len(par), len(ridx), 256), dtype='f4')
    std = np.zeros((len(std_par), len(ridx), 256), dtype='f4')
    for tag, off, plist, arr in (('mean', OFF_MEAN, par, mean), ('std', OFF_STD, std_par, std)):
        for ip, p in enumerate(plist):
            sub0 = p * 256
            for i in range(0, len(ridx), 170):
                chunk = ridx[i:i + 170]
                rngs = [(off + (int(r) * NPIX + sub0) * 4, off + (int(r) * NPIX + sub0) * 4 + 1023) for r in chunk]
                vals = multirange(rngs)
                for j, v in enumerate(vals):
                    arr[ip, i + j] = np.frombuffer(v, dtype='>f4')
            print(tag, ip + 1, '/', len(plist), f'{NBYTES / 1e6:.1f} MB', f'{time.time() - t0:.0f}s')
    np.savez_compressed(OUTD / 'eden_subcube.npz', mean=mean, std=std, parents16=np.array(par),
                        std_parents16=np.array(std_par), ridx=ridx, r_centers=rc, r_bounds=rb,
                        meta=json.dumps({'url': URL, 'zenodo_record': 10658339, 'accessed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                         'units': 'E (Zhang, Green & Rix 2023) per pc',
                                         'layout': 'nside256 NEST children p*256 ... p*256+255 of nside16 parent p',
                                         'bytes_downloaded': NBYTES}))
    log(f'range reads: {len(par)} nside16 parents mean, {len(std_par)} std, radii {rc[ridx[0]]:.1f}-{rc[ridx[-1]]:.1f} pc, {NBYTES} bytes')
    print('done', NBYTES / 1e6, 'MB', time.time() - t0, 's')
