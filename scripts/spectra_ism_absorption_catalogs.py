#!/usr/bin/env python
"""spectra_ism: literature Na I D / Ca II K interstellar absorption toward stars around the Silk Nebula.

Inputs (VizieR CfA mirror https://vizier.cfa.harvard.edu/viz-bin/asu-tsv, fetched 2026-11-16 by
spectra_ism_fetch_vizier.py and an ad-hoc call logged in data/raw/spectra_ism_fetch_log.jsonl):
  welsh2010_naicaii.tsv   J/A+A/510/A54   Welsh et al. 2010 (Na I + Ca II compilation, 1857 stars)
  lallement2003_nai.tsv   J/A+A/411/447   Lallement et al. 2003 (Na I D2, 1005 sightlines)
  sfeir1999_nai.tsv       J/A+A/346/785   Sfeir et al. 1999 (Na I, OHP)
  megier2009_caii.tsv     J/A+A/507/833   Megier et al. 2009 (Ca II, OB stars)
  hunter2006_caii_nai.tsv J/MNRAS/367/1478 Hunter et al. 2006 (UV Na I / Ca II components)
  gudennavar2012_ism_columns.tsv J/ApJS/199/8 (UV/optical column compilation; C IV etc.)
Distances: SIMBAD parallaxes (CfA SIMBAD TAP, basic.plx_value; mostly Gaia DR3) fetched here;
catalogue (mostly Hipparcos) distances kept for reference.
Output: data/derived/spectra_ism_abs_catalog_stars.csv, spectra_ism_abs_catalog_summary.json
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_absorption_catalogs.py
"""
import csv
import io
import json
import re
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, fetch, filament_coords, spine_arrays  # noqa: E402
from spectra_ism_vizier_parse import read_asu, fnum  # noqa: E402

V = RAW / 'spectra_ism_vizier'
RMAX = 12.0          # deg from the spine
DMAX = 1500.0        # pc


def spine_sep(ra, dec):
    sra, sdec, _ = spine_arrays()
    t = np.linspace(0, 1, 1500)
    cs = SkyCoord(np.interp(t, np.linspace(0, 1, len(sra)), sra) * u.deg,
                  np.interp(t, np.linspace(0, 1, len(sdec)), sdec) * u.deg)
    _, sep, _ = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg).match_to_catalog_sky(cs)
    return sep.deg


def hdkey(name):
    s = re.sub(r'\s+', ' ', str(name).strip())
    m = re.match(r'^(HD|HIP|BD)\s*([+-]?\d[\d ]*)', s)
    return f'{m.group(1)} {m.group(2).replace(" ", "")}' if m else s


stars = {}


def add(key, ra, dec, **kw):
    r = stars.setdefault(key, {'name': key, 'ra': ra, 'dec': dec, 'src': []})
    for k, v in kw.items():
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            continue
        if k not in r or r[k] in ('', None):
            r[k] = v
    return r


# Welsh et al. 2010
W = read_asu(V / 'welsh2010_naicaii.tsv')['J/A+A/510/A54/table1']
for r in W:
    ra, de = fnum(r['_RA']), fnum(r['_DE'])
    if not np.isfinite(ra):
        continue
    nm = r['Name'].strip()
    key = ('HD ' + nm) if nm.isdigit() else hdkey(nm)
    s = add(key, ra, de, d_cat=fnum(r['Dist']), l_d_cat=r['l_Dist'],
            EW_D2_mA=fnum(r['EWD2']),  # VizieR unit 0.1 pm = 1 mA
            l_EW_D2=r['l_EWD2'], logN_NaI=fnum(r['N(NaI)']), l_logN_NaI=r['l_N(NaI)'],
            EW_K_mA=fnum(r['EWCaII']), l_EW_K=r['l_EWCaII'], logN_CaII=fnum(r['N(CaII)']), l_logN_CaII=r['l_N(CaII)'])
    s['src'].append('W10')

# Lallement et al. 2003 (both tables)
L = read_asu(V / 'lallement2003_nai.tsv')
for tname, rows in L.items():
    for r in rows:
        ra = fnum(r.get('_RA.icrs', r.get('_RA', '')))
        de = fnum(r.get('_DE.icrs', r.get('_DE', '')))
        if not np.isfinite(ra):
            continue
        key = hdkey(r.get('Star', ''))
        ew = fnum(r.get('W(D2)', ''))       # VizieR labels the unit 0.1 nm, but values are mA (e.g. 36.5 for HD 141513)
        s = add(key, ra, de, d_cat=fnum(r.get('Dist', '')), EW_D2_mA=ew, l_EW_D2=r.get('l_W(D2)', ''))
        s['src'].append('L03')

# Sfeir et al. 1999
S = read_asu(V / 'sfeir1999_nai.tsv')
for tname, rows in S.items():
    for r in rows:
        ra, de = fnum(r.get('_RA.icrs', '')), fnum(r.get('_DE.icrs', ''))
        if not np.isfinite(ra):
            continue
        hd = r.get('HD', '').strip()
        key = f'HD {hd}' if hd.isdigit() else hdkey(r.get('Name', hd))
        s = add(key, ra, de, d_cat=fnum(r.get('Dist', '')), EW_D2_mA=fnum(r.get('EWD2', '')), l_EW_D2=r.get('l_EWD2', ''))
        s['src'].append('S99')

# Megier et al. 2009 (Ca II, distant OB stars)
for r in read_asu(V / 'megier2009_caii.tsv')['J/A+A/507/833/table1']:
    ra, de = fnum(r['_RA.icrs']), fnum(r['_DE.icrs'])
    plx = fnum(r['plx'])
    s = add(hdkey(r['Name']), ra, de, d_cat=1000 / plx if plx > 0 else fnum(r['DCaII']), EW_K_mA=fnum(r['EWK']))
    s['src'].append('M09')

# Hunter et al. 2006 component table (LSR velocities)
H = read_asu(V / 'hunter2006_caii_nai.tsv')
hcomp = {}
for r in H.get('J/MNRAS/367/1478/table3', []):
    hcomp.setdefault('HD ' + r['HD'].strip(), []).append(r)
for r in H.get('J/MNRAS/367/1478/stars', []):
    ra, de = fnum(r['_RA']), fnum(r['_DE'])
    s = add('HD ' + r['HD'].strip(), ra, de, d_cat=fnum(r['Dist']), logN_CaII_3393=fnum(r['logNCaIIp']))
    s['src'].append('H06')

# Gudennavar et al. 2012: UV columns (C IV) and N(HI)
G = read_asu(V / 'gudennavar2012_ism_columns.tsv')
guv = []
for tname in ('J/ApJS/199/8/table1', 'J/ApJS/199/8/table2'):
    for r in G.get(tname, []):
        ra, de = fnum(r['_RA']), fnum(r['_DE'])
        if not np.isfinite(ra):
            continue
        guv.append(dict(name=r['Name'].strip(), ra=ra, dec=de, sp=r['SpType'].strip(), dist_kpc=fnum(r['Dist']),
                        N_CIV=fnum(r['N(CIV)']), N_HI=fnum(r['N(HI)']), N_H2=fnum(r['N(H2)']), N_FeII=fnum(r['N(FeII)']),
                        N_OI=fnum(r['N(OI)']), EBV=fnum(r['E(B-V)']), table=tname.split('/')[-1]))

# --------------------------------------------------------------------------------------------
keys = list(stars)
ra = np.array([stars[k]['ra'] for k in keys])
dec = np.array([stars[k]['dec'] for k in keys])
sep = spine_sep(ra, dec)
sel = [k for k, s_ in zip(keys, sep) if s_ < RMAX]
print('stars with literature Na I / Ca II within', RMAX, 'deg of the spine:', len(sel))

# SIMBAD parallaxes for the selection (CfA mirror)
ids = sorted(set(sel))
q = "SELECT i.id, b.main_id, b.plx_value, b.plx_err, b.plx_bibcode, b.sp_type, b.rvz_radvel FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE i.id IN (" + \
    ','.join("'" + x.replace("'", "") + "'" for x in ids) + ')'
dest = RAW / 'spectra_ism_vizier' / 'simbad_plx_abs_catalog_stars.csv'
if not dest.exists():
    fetch('https://simbad.cfa.harvard.edu/simbad/sim-tap/sync', dest,
          params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': q}, timeout=300)
simb = {}
if dest.exists():
    for r in csv.DictReader(io.StringIO(dest.read_text(errors='replace'))):
        simb[re.sub(r'\s+', ' ', r['id'].strip())] = r
print('SIMBAD matches', len(simb))

out = []
for k in sel:
    s = stars[k]
    sb = simb.get(k, {})
    plx, eplx = fnum(sb.get('plx_value')), fnum(sb.get('plx_err'))
    d_plx = 1000 / plx if plx > 0 and eplx / plx < 0.2 else np.nan
    d_best = d_plx if np.isfinite(d_plx) else fnum(s.get('d_cat'))
    dp, dsp, bey = filament_coords(s['ra'], s['dec'])
    out.append(dict(name=k, main_id=sb.get('main_id', ''), ra=round(s['ra'], 5), dec=round(s['dec'], 5),
                    sep_spine_deg=round(float(spine_sep([s['ra']], [s['dec']])[0]), 3), dperp_deg=round(float(dp[0]), 3),
                    dec_spine=round(float(dsp[0]), 2), beyond_spine_end=bool(bey[0]), sp=sb.get('sp_type', ''),
                    plx_mas=plx, eplx_mas=eplx, plx_bib=sb.get('plx_bibcode', ''), d_plx_pc=round(d_plx, 1) if np.isfinite(d_plx) else '',
                    d_cat_pc=s.get('d_cat', ''), d_best_pc=round(d_best, 1) if np.isfinite(d_best) else '',
                    EW_D2_mA=s.get('EW_D2_mA', ''), l_EW_D2=s.get('l_EW_D2', ''), logN_NaI=s.get('logN_NaI', ''),
                    l_logN_NaI=s.get('l_logN_NaI', ''), EW_K_mA=s.get('EW_K_mA', ''), l_EW_K=s.get('l_EW_K', ''),
                    logN_CaII=s.get('logN_CaII', ''), l_logN_CaII=s.get('l_logN_CaII', ''), src='+'.join(sorted(set(s['src'])))))
out.sort(key=lambda r: r['sep_spine_deg'])
with open(DER / 'spectra_ism_abs_catalog_stars.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)

# UV (C IV etc.) compilation entries near the filament
gra = np.array([g['ra'] for g in guv])
gde = np.array([g['dec'] for g in guv])
gsep = spine_sep(gra, gde)
guv_near = [dict(g, sep_spine_deg=round(float(s_), 2)) for g, s_ in zip(guv, gsep) if s_ < 15]
guv_near.sort(key=lambda g: g['sep_spine_deg'])

summ = {'n_catalog_stars_within_12deg': len(out),
        'n_with_d_lt_1500': int(sum(1 for r in out if r['d_best_pc'] != '' and r['d_best_pc'] < DMAX)),
        'gudennavar_within_15deg': guv_near,
        'hunter2006_components_within_15deg': {k: hcomp[k] for k in hcomp
                                               if k in stars and spine_sep([stars[k]['ra']], [stars[k]['dec']])[0] < 15}}
json.dump(summ, open(DER / 'spectra_ism_abs_catalog_summary.json', 'w'), indent=1, default=float)

print(f"{'name':12s} {'sep':>5s} {'dperp':>6s} {'decsp':>6s} {'d_plx':>7s} {'d_cat':>6s} {'EWD2':>7s} {'logNa':>6s} {'EWK':>6s} {'logCa':>6s} src")
for r in out:
    if r['d_best_pc'] == '' or r['d_best_pc'] > DMAX:
        continue
    print(f"{r['name']:12s} {r['sep_spine_deg']:5.2f} {r['dperp_deg']:+6.2f} {r['dec_spine']:6.2f} {str(r['d_plx_pc']):>7s} "
          f"{str(r['d_cat_pc']):>6s} {r['l_EW_D2']}{str(r['EW_D2_mA']):>6s} {r['l_logN_NaI']}{str(r['logN_NaI']):>5s} "
          f"{r['l_EW_K']}{str(r['EW_K_mA']):>5s} {r['l_logN_CaII']}{str(r['logN_CaII']):>5s} {r['src']}")
print('Gudennavar UV entries within 15 deg:', len(guv_near))
for g in guv_near[:30]:
    print(g)
