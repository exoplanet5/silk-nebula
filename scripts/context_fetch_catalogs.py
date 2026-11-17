"""Fetch catalog extracts for the Silk Nebula astrophysical-context study (topic tag: context).

All results are saved as raw files under data/raw/context_*/ with a provenance log
(data/raw/context_catalogs/context_fetch_log.json: URL, catalogue, access date).

Sources (accessed 2026-11-16):
  * VizieR CfA mirror  https://vizier.cfa.harvard.edu/viz-bin/asu-tsv   (vizier.cds.unistra.fr was down / "Page unavailable")
      VII/297  Green (2025) Galactic SNR catalogue
      J/MNRAS/508/3877/maincat  Gentile Fusillo+ (2021) Gaia EDR3 white dwarfs
      J/A+A/662/A40/hotsd, knownhsd  Culpan+ (2022) hot subdwarfs
      J/A+A/635/A193/sdcatdr2  Geier (2020) known hot subdwarfs (with spectroscopic Teff)
      V/84/main, V/84/dist  Acker+ (1992) Strasbourg-ESO PNe
  * SIMBAD TAP CfA mirror  https://simbad.cfa.harvard.edu/simbad/sim-tap  (nebulae, hot stars)
Usage:  ~/.venvs/astro313/bin/python -I scripts/context_fetch_catalogs.py
"""
import json, time, sys
from pathlib import Path
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/raw/context_catalogs'
OUT.mkdir(parents=True, exist_ok=True)
VIZ = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
SIMBAD_TAP = 'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
CEN = (116.68, 16.42)
LOG = OUT / 'context_fetch_log.json'
log = json.loads(LOG.read_text()) if LOG.exists() else []


def get(url, params, fname, tries=3, method='get'):
    for k in range(tries):
        try:
            r = (requests.post if method == 'post' else requests.get)(url, params=params if method == 'get' else None,
                                                                     data=params if method == 'post' else None, timeout=300)
            r.raise_for_status()
            (OUT / fname).write_bytes(r.content)
            log.append({'file': str(OUT / fname), 'url': r.url if method == 'get' else url, 'params': params,
                        'accessed': '2026-11-16', 'bytes': len(r.content)})
            print(f'{fname}: {len(r.content)} bytes')
            return
        except Exception as e:
            print('retry', fname, e)
            time.sleep(5 * (k + 1))
    print('FAILED', fname)


def viz(source, fname, extra=None, cone=None, out='**'):
    p = {'-source': source, '-out.max': 'unlimited', '-out.add': '_r', '-oc.form': 'd', '-out': out}
    if out == '**':
        p['-out.all'] = ''
        p.pop('-out')
    if cone:
        p.update({'-c': f'{CEN[0]} {CEN[1]}', '-c.r': str(cone), '-c.u': 'deg'})
    if extra:
        p.update(extra)
    get(VIZ, p, fname)


def tap(query, fname):
    get(SIMBAD_TAP, {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query}, fname, method='post')


if __name__ == '__main__':
    which = sys.argv[1:] or ['green', 'wd', 'sd', 'pn', 'simbad']
    if 'green' in which:
        viz('VII/297/snrs', 'green2025_VII297_all.tsv')
    if 'wd' in which:
        # hot WD candidates within 20 deg; Pwd>0.75, parallax > 1 mas (d < ~1 kpc) ; TeffH > 20 kK
        viz('J/MNRAS/508/3877/maincat', 'gf21_wd_cone20_hot.tsv', cone=20,
            extra={'Pwd': '>0.75', 'Plx': '>1', 'TeffH': '>20000'},
            out='WDJname,GaiaEDR3,RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,e_pmRA,pmDE,e_pmDE,Pwd,Gmag,BP-RP,GMAG,meanAV,'
                'TeffH,e_TeffH,loggH,MassH,TeffHe,e_TeffHe,loggHe,MassHe,rgeo,b_rgeo,B_rgeo,SDSS12,umag,gmag,GLON,GLAT')
    if 'sd' in which:
        viz('J/A+A/662/A40/hotsd', 'culpan22_hotsd_cone20.tsv', cone=20)
        viz('J/A+A/662/A40/knownhsd', 'culpan22_knownhsd_cone20.tsv', cone=20)
        viz('J/A+A/635/A193/sdcatdr2', 'geier20_sdcat_cone20.tsv', cone=20)
    if 'pn' in which:
        viz('V/84/main', 'acker92_pn_cone25.tsv', cone=25)
        viz('V/84/dist', 'acker92_pn_dist_all.tsv')
    if 'simbad' in which:
        box = "ra BETWEEN 100 AND 135 AND dec BETWEEN -2 AND 36"
        tap("SELECT main_id, ra, dec, otype, sp_type, plx_value, plx_err, pmra, pmdec, rvz_radvel, galdim_majaxis, "
            "galdim_minaxis, nbref FROM basic WHERE " + box +
            " AND otype IN ('HII','SNR','SR?','PN','PN?','RNe','ISM','Cld','GNe','DNe','EmO','MoC','sh','bub','SFR',"
            "'HH','glb','cor','CGb','HVC','HI','LSB','out','Neb')", 'simbad_nebulae_box.csv')
        cone = ("CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %.3f, %.3f, 25)) = 1" % CEN)
        tap("SELECT main_id, ra, dec, otype, sp_type, plx_value, plx_err, pmra, pmdec, rvz_radvel, flux.flux AS vmag "
            "FROM basic LEFT JOIN flux ON flux.oidref = basic.oid AND flux.filter = 'V' WHERE " + cone +
            " AND (sp_type LIKE 'O%' OR sp_type LIKE 'B0%' OR sp_type LIKE 'B1%' OR sp_type LIKE 'B2%' "
            " OR sp_type LIKE 'sdO%' OR sp_type LIKE 'DO%' OR sp_type LIKE 'DAO%' OR sp_type LIKE 'PG1159%' "
            " OR sp_type LIKE 'WR%' OR sp_type LIKE 'W%' OR sp_type LIKE 'sdB%' OR sp_type LIKE 'DA%' OR otype IN ('WD*','WD?','HS*','HS?','PN'))",
            'simbad_hotstars_cone25.csv')
    LOG.write_text(json.dumps(log, indent=1))
