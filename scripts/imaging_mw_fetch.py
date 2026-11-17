"""Download multiwavelength imaging for the Silk Nebula imaging_mw study.

Usage:  ~/.venvs/astro313/bin/python -I scripts/imaging_mw_fetch.py <group> [<group> ...]
Groups: skyview_eq, skyview_gal, erosita, uvbkgd, galfa, lab, wssa (see functions below).

Every file is logged to data/raw/<archive>/MANIFEST.jsonl with URL, release, bytes, sha256 and
access date. Downloaded content is treated as data only (FITS / text), never executed.
"""
import hashlib
import json
import sys
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw'
ACCESS_DATE = '2026-11-16'
UA = {'User-Agent': 'silk-nebula-research-script/1.0 (python-requests)'}


def log(archive, rec):
    p = RAW / archive / 'MANIFEST.jsonl'
    p.parent.mkdir(parents=True, exist_ok=True)
    rec = {'access_date': ACCESS_DATE, 'utc': datetime.now(timezone.utc).isoformat(timespec='seconds'), **rec}
    with p.open('a') as f:
        f.write(json.dumps(rec) + '\n')


def fetch(url, dest, archive, release, params=None, timeout=600, retries=3, min_bytes=2000, note=''):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > min_bytes:
        print('exists', dest.name)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    full = url + ('?' + urllib.parse.urlencode(params) if params else '')
    for k in range(retries):
        try:
            with requests.get(url, params=params, stream=True, timeout=timeout, headers=UA) as r:
                r.raise_for_status()
                h = hashlib.sha256()
                tmp = dest.with_suffix(dest.suffix + '.part')
                n = 0
                with tmp.open('wb') as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk); h.update(chunk); n += len(chunk)
                if n < min_bytes:
                    raise IOError(f'only {n} bytes')
                tmp.rename(dest)
                log(archive, {'url': full, 'file': str(dest.relative_to(ROOT)), 'bytes': n,
                              'sha256': h.hexdigest(), 'release': release, 'note': note})
                print(f'ok {dest.name} {n/1e6:.1f} MB')
                return dest
        except Exception as e:  # network / HTTP errors: retry politely
            print(f'retry {k+1}/{retries} {dest.name}: {e}')
            time.sleep(10 * (k + 1))
    log(archive, {'url': full, 'file': str(dest.relative_to(ROOT)), 'status': 'FAILED', 'release': release})
    return None


# ----------------------------------------------------------------------------------------------
# SkyView (https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl). TAN cutouts, later reprojected
# onto the common grids by imaging_mw_common.py.
SV = 'https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl'
# survey short name, pixel size (arcmin), sampler
EQ_SURVEYS = [
    ('RASS Background 1', 12, 'NN'), ('RASS Background 2', 12, 'NN'), ('RASS Background 3', 12, 'NN'),
    ('RASS Background 4', 12, 'NN'), ('RASS Background 5', 12, 'NN'), ('RASS Background 6', 12, 'NN'),
    ('RASS Background 7', 12, 'NN'),
    ('RASS-Int Soft', 4, 'Clip'), ('RASS-Int Hard', 4, 'Clip'), ('RASS-Int Broad', 4, 'Clip'),
    ('RASS-Cnt Soft', 4, 'Clip'), ('RASS-Cnt Broad', 4, 'Clip'),
    ('IRIS 100', 2, 'Clip'), ('IRIS 60', 2, 'Clip'), ('IRIS 25', 2, 'Clip'), ('IRIS 12', 2, 'Clip'),
    ('SFD100m', 2, 'Clip'), ('SFD Dust Map', 2, 'Clip'),
    ('Planck 857 I', 2, 'Clip'), ('Planck 545 I', 2, 'Clip'), ('Planck 353 PI', 10, 'Clip'),
    ('Planck 353 PA', 10, 'NN'),
    ('AKARI WIDE-S', 2, 'Clip'), ('AKARI N160', 2, 'Clip'),
    ('HI4PI', 6, 'Clip'), ('EBHIS', 6, 'Clip'),
    ('408MHz', 10, 'Clip'), ('1420MHz (Bonn)', 10, 'Clip'), ('4850MHz', 2, 'Clip'),
    ('Planck 030 I', 10, 'Clip'), ('Planck 030 PI', 10, 'Clip'), ('WMAP K', 10, 'Clip'),
    ('DRAO 22 MHz', 20, 'Clip'), ('0035MHz', 20, 'Clip'),
    ('H-Alpha Comp', 3, 'Clip'), ('CO', 7.5, 'NN'),
    ('Fermi 3', 10, 'Clip'), ('Fermi 4', 10, 'Clip'), ('Fermi 5', 10, 'Clip'),
]
GAL_SURVEYS = [
    ('RASS Background 1', 12, 'NN'), ('RASS Background 2', 12, 'NN'), ('RASS Background 4', 12, 'NN'),
    ('RASS Background 5', 12, 'NN'), ('RASS Background 6', 12, 'NN'), ('RASS Background 7', 12, 'NN'),
    ('IRIS 100', 6, 'Clip'), ('Planck 857 I', 6, 'Clip'), ('HI4PI', 8, 'Clip'), ('H-Alpha Comp', 6, 'Clip'),
    ('408MHz', 15, 'Clip'), ('1420MHz (Bonn)', 15, 'Clip'), ('Planck 030 PI', 15, 'Clip'),
    ('Planck 030 I', 15, 'Clip'), ('WMAP K', 15, 'Clip'), ('Fermi 3', 15, 'Clip'), ('Fermi 4', 15, 'Clip'),
    ('Fermi 5', 15, 'Clip'), ('Planck 353 PA', 15, 'NN'), ('Planck 353 PI', 15, 'Clip'),
]


def skyview(field):
    if field == 'eq':
        pos, coords, size, surveys = '116.5,17.0', 'J2000', (21.0, 29.0), EQ_SURVEYS
    else:
        pos, coords, size, surveys = '205.0,17.5', 'Galactic', (44.0, 39.0), GAL_SURVEYS
    for name, pix_arcmin, sampler in surveys:
        nx = int(round(size[0] * 60 / pix_arcmin)); ny = int(round(size[1] * 60 / pix_arcmin))
        params = {'Position': pos, 'Survey': name, 'Coordinates': coords, 'Projection': 'Tan',
                  'Pixels': f'{nx},{ny}', 'Size': f'{size[0]},{size[1]}', 'Sampler': sampler,
                  'Return': 'FITS'}
        if sampler == 'Clip':
            params['ClipIntensive'] = ''  # treat pixel values as intensive (average, not sum)
        tag = name.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')
        fetch(SV, RAW / 'skyview_imaging_mw' / f'{field}_{tag}.fits', 'skyview_imaging_mw',
              'SkyView current (2026-11 access); survey ' + name, params=params, timeout=900,
              note=f'TAN {size} deg @ {pix_arcmin} arcmin, sampler {sampler}')
        time.sleep(2)


def erosita():
    base = 'https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/HalfSkyMaps_dr1/maps_zheng/'
    for f in ['e01_TM8_0.200_0.400keV_3arcm_c947_AIT_sum_Ratesmap.fits',
              'e01_TM8_0.400_0.600keV_3arcm_c947_AIT_sum_Ratesmap.fits',
              'e01_TM8_0.600_1.000keV_3arcm_c947_AIT_sum_Ratesmap.fits',
              'e01_TM8_1.000_2.300keV_3arcm_c947_AIT_sum_Ratesmap.fits',
              'e01_TM8_2.300_5.000keV_3arcm_c947_AIT_sum_Ratesmap.fits',
              'rat_0.534_0.616_keV_tm8.fits', 'rat_0.614_0.694_keV_tm8.fits']:
        fetch(base + f, RAW / 'erosita_dr1' / f, 'erosita_dr1',
              'eROSITA-DE DR1 eRASS1 half-sky maps (Zheng et al. 2024, A&A 681, A77)', timeout=1200)


def uvbkgd():
    base = 'https://archive.stsci.edu/missions/hlsp/uv-bkgd/aitoff/'
    for f in ['hlsp_uv-bkgd_galex_map_allsky_fuv_v1_skymap.fits', 'hlsp_uv-bkgd_galex_map_allsky_nuv_v1_skymap.fits']:
        fetch(base + f, RAW / 'mast_uvbkgd' / f, 'mast_uvbkgd',
              'MAST HLSP uv-bkgd v1 (Murthy 2014, ApJS 213, 32), GALEX GR6+GR7', timeout=1200)


def galfa():
    api = 'https://galfa-dr2-cutout-dptnqjtpkq-uc.a.run.app/cutout'
    p = {'product': 'wide', 'ra': 116.5, 'dec': 17.0, 'width': 17.0, 'height': 26.0, 'v_lo': -30, 'v_hi': 40}
    fetch(api, RAW / 'galfa_hi_dr2' / 'galfa_dr2_wide_ra116.5_dec17_17x26_vm30p40.fits', 'galfa_hi_dr2',
          'GALFA-HI DR2 (Peek et al. 2018, ApJS 234, 2) wide cubes via www.galfa-hi.space cutout service',
          params=p, timeout=3600)


def lab():
    fetch('https://lambda.gsfc.nasa.gov/data/foregrounds/HI/lab_healpix.tar.bz2',
          RAW / 'lambda_imaging_mw' / 'lab_healpix.tar.bz2', 'lambda_imaging_mw',
          'LAMBDA LAB HI survey HEALPix (Kalberla et al. 2005; Land & Slosar 2007 repackaging)', timeout=1800)


if __name__ == '__main__':
    for g in sys.argv[1:]:
        {'skyview_eq': lambda: skyview('eq'), 'skyview_gal': lambda: skyview('gal'), 'erosita': erosita,
         'uvbkgd': uvbkgd, 'galfa': galfa, 'lab': lab}[g]()
