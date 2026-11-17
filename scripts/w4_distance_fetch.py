"""w4_distance: small, logged catalogue fetches for the distance close-out (wave 4, tag `distance`).

Writes ONLY into data/raw/w4_distance_cat/ (new directory) and appends to its fetch_log.tsv.
Downloaded content is treated as data; run with `python -I`.

Fetches
  pgcc      VizieR J/A+A/594/A28/pgcc cone 2 deg around PGCC G199.09+19.95 (SIMBAD position 115.28 +20.833)
  zhang     (tried; mirror returned no rows) VizieR J/MNRAS/524/1855/xpparams (Zhang, Green & Rix 2023 XP extinctions) cone 1.0 deg around the PGCC
            and around three control / ridge positions
  starhorse VizieR I/354/starhorse2021 (Anders+2022) cones 1.0 deg, AV/distance percentiles
  gaiapgcc  Data Lab gaia_dr3.gaia_source cone 1.0 deg around the PGCC (parallax, GSP-Phot A0/E(BP-RP)) for the join
  hot       Data Lab gaia_dr3.gaia_source: blue (bp_rp < 0.10) stars with parallax > 1.1 mas in a box enclosing both
            strands (RA 112.6-118.4, Dec 11.3-29.0); distances to the strands are computed later (w4_distance_targets.py)
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_fetch.py [pgcc zhang gaiapgcc hot]
"""
import sys, time, datetime, urllib.parse
from pathlib import Path
import requests
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/w4_distance_cat'
OUT.mkdir(parents=True, exist_ok=True)
LOG = OUT / 'fetch_log.tsv'
VIZ = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
TAP = 'https://datalab.noirlab.edu/tap/sync'
PGCC = (115.28, 20.833)
# control / context cones for the Zhang+23 extinction-distance test (RA, Dec)
CONES = {
    'pgcc': PGCC,
    'ridgeS_P1': (116.65, 16.60),       # south spine (P1)
    'east4_Dec21': (120.00, 21.00),     # ~4 deg east of the ridge at Dec 21 (dust 'rising east')
    'north_Dec27': (114.70, 27.00),     # north arc
    'west3_Dec21': (112.30, 21.00),     # ~3 deg west of the patch (control)
}


def log(url, fname, status, nbytes, note=''):
    new = not LOG.exists()
    with open(LOG, 'a') as f:
        if new:
            f.write('time_utc\turl\tfile\tstatus\tbytes\tnote\n')
        f.write(f'{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{url}\t{fname}\t{status}\t{nbytes}\t{note}\n')


def get(url, fname, params=None, tries=3, timeout=300):
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            full = r.url
            (OUT / fname).write_bytes(r.content)
            log(full, fname, r.status_code, len(r.content))
            print(fname, r.status_code, len(r.content))
            if r.status_code == 200:
                return True
        except Exception as e:  # network hiccup
            log(url, fname, 'ERR', 0, repr(e)[:200])
            print('retry', k, e)
        time.sleep(4 * (k + 1))
    return False


def viz(source, ra, dec, rad, fname, extra=None, cols='-out.all'):
    p = {'-source': source, '-c': f'{ra} {dec}', '-c.r': str(rad), '-c.u': 'deg', '-out.max': 'unlimited',
         '-oc.form': 'd', '-out.add': '_r'}
    if cols == '-out.all':
        p['-out.all'] = ''
    else:
        p['-out'] = cols
    if extra:
        p.update(extra)
    return get(VIZ, fname, params=p)


def tap(query, fname):
    p = {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query}
    return get(TAP, fname, params=p, timeout=600)


GAIA_COLS = ('source_id, ra, dec, l, b, parallax, parallax_error, parallax_over_error, pmra, pmdec, ruwe, '
             'phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, radial_velocity, '
             'teff_gspphot, logg_gspphot, mh_gspphot, azero_gspphot, ag_gspphot, ebpminrp_gspphot, '
             'distance_gspphot, distance_gspphot_lower, distance_gspphot_upper, classprob_dsc_combmod_star, '
             'has_xp_continuous')

if __name__ == '__main__':
    what = sys.argv[1:] or ['pgcc', 'starhorse', 'gaiapgcc', 'hot']
    if 'pgcc' in what:
        viz('J/A+A/594/A28/pgcc', *PGCC, 2.0, 'pgcc_cone2deg.tsv')
    if 'zhang' in what:   # NOTE 2026-11-17: the CfA mirror returned the header but NO rows for this table
        for name, (ra, dec) in CONES.items():
            viz('J/MNRAS/524/1855/xpparams', ra, dec, 1.0, f'zhang23_{name}_r1.0.tsv')
    if 'starhorse' in what:  # Anders+2022 StarHorse (Gaia EDR3 + PS1/2MASS/WISE): AV and distance percentiles
        cols = ('Source,RA_ICRS,DE_ICRS,teff50,logg50,dist16,dist50,dist84,AV16,AV50,AV84,fidelity,FlagOut,'
                'FlagPhoto,AG50,GMAG0,BP-RP0')
        for name, (ra, dec) in CONES.items():
            viz('I/354/starhorse2021', ra, dec, 1.0, f'starhorse_{name}_r1.0.tsv')  # -out.all (a column list returned no rows)
    if 'gaiapgcc' in what:
        for name, (ra, dec) in CONES.items():
            q = (f"SELECT {GAIA_COLS} FROM gaia_dr3.gaia_source WHERE "
                 f"ra BETWEEN {ra - 1.0 / np.cos(np.radians(dec)):.4f} AND {ra + 1.0 / np.cos(np.radians(dec)):.4f} AND dec BETWEEN {dec - 1.0} AND {dec + 1.0} AND phot_g_mean_mag < 17.8 AND parallax > 0.3 AND parallax < 100")
            tap(q, f'gaia_{name}_r1.0.csv')
    if 'hot' in what:
        q = (f"SELECT {GAIA_COLS} FROM gaia_dr3.gaia_source WHERE ra BETWEEN 112.6 AND 118.4 AND "
             f"dec BETWEEN 11.3 AND 29.0 AND parallax > 1.1 AND bp_rp < 0.10")
        tap(q, 'gaia_hot_box.csv')
