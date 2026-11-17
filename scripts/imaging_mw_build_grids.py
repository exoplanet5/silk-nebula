"""Reproject all downloaded imaging onto the common EQ (RA/Dec, 3') and GAL (l/b, 6') grids.

Output: data/derived/imaging_mw_grids/imaging_mw_{eq,gal}_<name>.fits and a catalogue
data/derived/imaging_mw_grids/grid_catalog.json (name -> source file, unit, conversion notes).
Run: ~/.venvs/astro313/bin/python -I scripts/imaging_mw_build_grids.py [skyview erosita uvbkgd lab]
"""
import json
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_common import RAW, DER, reproject_to, save_grid, grid_coords, EQ_WCS, GAL_WCS, EQ_SHAPE, GAL_SHAPE

CAT_P = DER / 'imaging_mw_grids' / 'grid_catalog.json'
(DER / 'imaging_mw_grids').mkdir(exist_ok=True, parents=True)
CAT = json.loads(CAT_P.read_text()) if CAT_P.exists() else {}


def add(name, which, arr, src, unit, note=''):
    p = save_grid(arr, name, which, unit, note)
    CAT[f'{which}:{name}'] = {'file': str(p.relative_to(DER.parent.parent)), 'source': str(src), 'unit': unit,
                              'note': note}
    print('saved', which, name, f'{np.nanmedian(arr):.4g}', unit)


SKY = {  # skyview tag -> (grid name, unit, note)
    'RASS_Background_1': ('rass_r1', '1e-6 cts/s/arcmin2', 'Snowden+97 R1 (0.11-0.284 keV)'),
    'RASS_Background_2': ('rass_r2', '1e-6 cts/s/arcmin2', 'Snowden+97 R2 (0.14-0.284 keV)'),
    'RASS_Background_3': ('rass_r3', '1e-6 cts/s/arcmin2', 'Snowden+97 R3 (0.2-0.83 keV)'),
    'RASS_Background_4': ('rass_r4', '1e-6 cts/s/arcmin2', 'Snowden+97 R4 (0.44-1.01 keV)'),
    'RASS_Background_5': ('rass_r5', '1e-6 cts/s/arcmin2', 'Snowden+97 R5 (0.56-1.21 keV)'),
    'RASS_Background_6': ('rass_r6', '1e-6 cts/s/arcmin2', 'Snowden+97 R6 (0.73-1.56 keV)'),
    'RASS_Background_7': ('rass_r7', '1e-6 cts/s/arcmin2', 'Snowden+97 R7 (1.05-2.04 keV)'),
    'RASS-Int_Soft': ('rass_int_soft', 'cts/s/pix(SkyView)', 'RASS3 soft 0.1-0.4 keV intensity, point sources included'),
    'RASS-Int_Hard': ('rass_int_hard', 'cts/s/pix(SkyView)', 'RASS3 hard 0.5-2 keV intensity'),
    'RASS-Int_Broad': ('rass_int_broad', 'cts/s/pix(SkyView)', 'RASS3 broad 0.1-2.4 keV intensity'),
    'IRIS_100': ('iris100', 'MJy/sr', 'IRIS 100um (Miville-Deschenes & Lagache 2005)'),
    'IRIS_60': ('iris60', 'MJy/sr', 'IRIS 60um'),
    'IRIS_25': ('iris25', 'MJy/sr', 'IRIS 25um'),
    'IRIS_12': ('iris12', 'MJy/sr', 'IRIS 12um'),
    'SFD100m': ('sfd100', 'MJy/sr', 'SFD98 100um'),
    'SFD_Dust_Map': ('sfd_ebv', 'mag', 'SFD98 E(B-V)'),
    'Planck_857_I': ('planck857', 'MJy/sr', 'Planck PR 857 GHz'),
    'Planck_545_I': ('planck545', 'MJy/sr', 'Planck PR 545 GHz'),
    'Planck_353_PI': ('planck353_pi', 'K_CMB', 'Planck 353 GHz polarized intensity'),
    'Planck_353_PA': ('planck353_pa', 'rad', 'Planck 353 GHz polarization angle (SkyView)'),
    'HI4PI': ('hi4pi_nhi', 'cm-2', 'HI4PI NHI |v|<600 km/s'),
    'EBHIS': ('ebhis_nhi', 'cm-2', 'EBHIS NHI'),
    '408MHz': ('haslam408', 'K', 'Haslam+82 408 MHz'),
    '1420MHz_Bonn': ('stockert1420', 'mK', 'Reich 1986 Stockert 1420 MHz'),
    '4850MHz': ('gb6_4850', 'Jy/beam', 'GB6 4850 MHz (large scales filtered)'),
    'Planck_030_I': ('planck030', 'K_CMB', 'Planck 30 GHz I'),
    'Planck_030_PI': ('planck030_pi', 'K_CMB', 'Planck 30 GHz polarized intensity'),
    'WMAP_K': ('wmapk', 'mK', 'WMAP 9yr K band (23 GHz)'),
    'DRAO_22_MHz': ('drao22', 'K', 'DRAO 22 MHz'),
    '0035MHz': ('gtee35', 'K?', 'Gauribidanur 34.5 MHz'),
    'H-Alpha_Comp': ('halpha_f03', 'R', 'Finkbeiner 2003 Halpha composite (WHAM/VTSS/SHASSA)'),
    'Fermi_3': ('fermi3', 'ph/cm2/s/sr', 'Fermi-LAT band 3 (SkyView)'),
    'Fermi_4': ('fermi4', 'ph/cm2/s/sr', 'Fermi-LAT band 4 (SkyView)'),
    'Fermi_5': ('fermi5', 'ph/cm2/s/sr', 'Fermi-LAT band 5 (SkyView)'),
    'NVSS': ('nvss', 'Jy/beam', 'NVSS 1.4 GHz (interferometer; >~15 arcmin scales resolved out)'),
    'TGSS_ADR1': ('tgss', 'Jy/beam', 'TGSS ADR1 150 MHz (interferometer)'),
    'AKARI_WIDE-S': ('akari90', 'MJy/sr', 'AKARI FIS WIDE-S 90um'),
}


def do_skyview():
    for which in ['eq', 'gal']:
        for tag, (name, unit, note) in SKY.items():
            f = RAW / 'skyview_imaging_mw' / f'{which}_{tag}.fits'
            if not f.exists():
                continue
            try:
                hd = fits.open(f)
            except Exception as e:
                print('skip', f.name, e); continue
            d = hd[0].data.astype(float)
            d[d <= -32767] = np.nan
            arr = reproject_to((d, hd[0].header), which)
            add(name, which, arr, f.relative_to(RAW), unit, note)


def do_erosita():
    base = RAW / 'erosita_dr1'
    bands = {'0.200_0.400': 'ero_0204', '0.400_0.600': 'ero_0406', '0.600_1.000': 'ero_0610',
             '1.000_2.300': 'ero_1023', '2.300_5.000': 'ero_2350'}
    for k, name in bands.items():
        f = base / f'e01_TM8_{k}keV_3arcm_c947_AIT_sum_Ratesmap.fits'
        hd = fits.open(f, memmap=True)
        d = np.array(hd[0].data, dtype=float)
        d[d == 0] = np.nan            # zero = outside the eROSITA-DE hemisphere / no exposure
        for which in ['eq', 'gal']:
            arr = reproject_to((d, hd[0].header), which) * 400.0   # cts/s per 9 arcmin2 -> cts/s/deg2
            add(name, which, arr, f.relative_to(RAW), 'cts/s/deg2 (TM8)',
                f'eRASS1 {k} keV rate map (Zheng+24), point sources included')
    for k, name in {'0.534_0.616': 'ero_ovii', '0.614_0.694': 'ero_oviii'}.items():
        f = base / f'rat_{k}_keV_tm8.fits'
        hd = fits.open(f)
        d = hd[0].data.astype(float)  # counts/exposure: zeros are genuine Poisson outcomes, keep them
        d[~np.isfinite(d)] = np.nan
        for which in ['eq', 'gal']:
            arr = reproject_to((d, hd[0].header), which)
            add(name, which, arr, f.relative_to(RAW), 'cts/s/deg2 (TM8)', f'eRASS1 narrow band {k} keV (Zheng+24)')


def do_uvbkgd():
    for band in ['fuv', 'nuv']:
        f = RAW / 'mast_uvbkgd' / f'hlsp_uv-bkgd_galex_map_allsky_{band}_v1_skymap.fits'
        hd = fits.open(f)
        for i, suffix in [(0, ''), (1, '_hdu1')]:
            d = hd[i].data.astype(float)
            d[(d <= 0) | ~np.isfinite(d)] = np.nan
            hdr = hd[0].header if i == 1 else hd[i].header   # HDU1 has no WCS of its own
            for which in ['eq', 'gal']:
                arr = reproject_to((d, hdr), which, order='nearest-neighbor')
                add(f'murthy_{band}{suffix}', which, arr, f.relative_to(RAW), 'CU (ph/cm2/s/sr/A)' if i == 0 else 'HDU1',
                    'Murthy 2014 GALEX diffuse Aitoff map, 0.1 deg pixels')


def do_lab():
    import healpy as hp
    from reproject import reproject_from_healpix
    base = RAW / 'lambda_imaging_mw' / 'lab_extract' / 'LAB_Healpix'
    for lo in [-20, -10, 0, 10, 20, 30]:
        f = base / f'LAB_cut{lo:+d}.fits'.replace('+0', '+0')
        m = hp.read_map(str(f), field=0, verbose=False) if 'verbose' in hp.read_map.__code__.co_varnames else hp.read_map(str(f), field=0)
        m = np.where(m == hp.UNSEEN, np.nan, m)
        for which, w, shp in [('eq', EQ_WCS, EQ_SHAPE), ('gal', GAL_WCS, GAL_SHAPE)]:
            arr, _ = reproject_from_healpix((m, 'galactic'), w, shape_out=shp, nested=False)
            nhi = arr * 10.0 * 1.823e18      # mean Tb (K) over a 10 km/s slice -> NHI (optically thin)
            add(f'lab_v{lo:+d}', which, nhi.astype(np.float32), f.relative_to(RAW), 'cm-2',
                f'LAB NHI V_LSR {lo}..{lo + 10} km/s (from mean Tb x 10 km/s)')


if __name__ == '__main__':
    jobs = sys.argv[1:] or ['skyview', 'erosita', 'uvbkgd', 'lab']
    for j in jobs:
        {'skyview': do_skyview, 'erosita': do_erosita, 'uvbkgd': do_uvbkgd, 'lab': do_lab}[j]()
        CAT_P.write_text(json.dumps(CAT, indent=1))
