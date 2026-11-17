"""Parse the ATNF Pulsar Catalogue text database (psrcat.db, v2.9.0) and extract pulsars near the Silk Nebula.

Input : data/raw/context_psr/psrcat_tar/psrcat.db  (downloaded 2026-11-16 from
        https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz ; only the text DB is read,
        the C code in the package is never compiled or executed)
Output: data/derived/context_pulsars_30deg.csv  (all pulsars within 30 deg of the Silk centre)
Usage : ~/.venvs/astro313/bin/python -I scripts/context_psrcat.py
"""
import sys, re
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import CENTER, ROOT

DB = ROOT / 'data/raw/context_psr/psrcat_tar/psrcat.db'
KEYS = ['PSRJ', 'PSRB', 'RAJ', 'DECJ', 'PMRA', 'PMDEC', 'PX', 'DM', 'P0', 'P1', 'F0', 'F1', 'DIST_A', 'DIST_DM',
        'DIST_DM1', 'TYPE', 'ASSOC', 'POSEPOCH', 'NAME', 'ELONG', 'ELAT', 'PMELONG', 'PMELAT']


def _err(valstr, errstr):
    """psrcat uncertainty is in units of the last quoted digit of the value (unless it contains '.' or 'e')."""
    try:
        e = float(errstr)
    except Exception:
        return np.nan
    if any(c in errstr for c in '.eE'):
        return e
    v = valstr.lower()
    mant = v.split('e')[0]
    exp = int(v.split('e')[1]) if 'e' in v else 0
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    return e * 10 ** (exp - dec)


def parse(path=DB):
    recs, cur = [], {}
    for line in path.read_text(errors='replace').splitlines():
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            if cur:
                recs.append(cur)
            cur = {}
            continue
        parts = line.split()
        if not parts:
            continue
        k = parts[0]
        if k in KEYS and k not in cur:
            if k in ('TYPE', 'ASSOC'):
                cur[k] = ' '.join(parts[1:])
                continue
            v = parts[1] if len(parts) > 1 else ''
            cur[k] = v
            if len(parts) > 2 and k in ('PMRA', 'PMDEC', 'PX', 'DM', 'DIST_A', 'PMELONG', 'PMELAT'):
                cur['e_' + k] = _err(v, parts[2])
    if cur:
        recs.append(cur)
    return recs


def fnum(x):
    try:
        return float(str(x).replace('D', 'e'))
    except Exception:
        return np.nan


def table(recs):
    rows = []
    for r in recs:
        pm_from_ecl = None
        try:
            if 'RAJ' in r and 'DECJ' in r:
                c = SkyCoord(r['RAJ'], r['DECJ'], unit=(u.hourangle, u.deg))
            elif 'ELONG' in r and 'ELAT' in r:
                ce = SkyCoord(lon=fnum(r['ELONG']) * u.deg, lat=fnum(r['ELAT']) * u.deg,
                              pm_lon_coslat=(fnum(r.get('PMELONG', 0)) if 'PMELONG' in r else 0) * u.mas / u.yr,
                              pm_lat=(fnum(r.get('PMELAT', 0)) if 'PMELAT' in r else 0) * u.mas / u.yr,
                              frame='barycentrictrueecliptic')
                c = SkyCoord(ce.icrs.ra, ce.icrs.dec)
                if 'PMELONG' in r and 'PMELAT' in r:
                    pm_from_ecl = (ce.icrs.pm_ra_cosdec.value, ce.icrs.pm_dec.value)
            else:
                continue
        except Exception:
            continue
        P0 = fnum(r.get('P0')) if 'P0' in r else 1 / fnum(r.get('F0', np.nan))
        if 'P1' in r:
            P1 = fnum(r['P1'])
        else:
            F0, F1 = fnum(r.get('F0', np.nan)), fnum(r.get('F1', np.nan))
            P1 = -F1 / F0 ** 2
        age = P0 / (2 * P1) / 3.15576e7 if P1 > 0 else np.nan
        edot = 4 * np.pi ** 2 * 1e45 * P1 / P0 ** 3 if P1 > 0 else np.nan
        px, epx = fnum(r.get('PX')), r.get('e_PX', np.nan)
        dA, dDM = fnum(r.get('DIST_A')), fnum(r.get('DIST_DM'))
        # catalogue logic for best distance: DIST_A, else significant PX, else YMW16 DM distance
        if np.isfinite(dA):
            d, dsrc = dA, 'DIST_A'
        elif np.isfinite(px) and np.isfinite(epx) and px / epx > 3:
            d, dsrc = 1 / px, 'PX'
        else:
            d, dsrc = dDM, 'DM(YMW16)'
        rows.append(dict(psrj=r.get('PSRJ'), psrb=r.get('PSRB', ''), ra=c.ra.deg, dec=c.dec.deg,
                         l=c.galactic.l.deg, b=c.galactic.b.deg,
                         sep_silk=c.separation(SkyCoord(*CENTER, unit='deg')).deg,
                         P0_s=P0, P1=P1, age_yr=age, edot=edot, dm=fnum(r.get('DM')), e_dm=r.get('e_DM', np.nan),
                         pmra=fnum(r.get('PMRA')) if pm_from_ecl is None else pm_from_ecl[0],
                         e_pmra=r.get('e_PMRA', r.get('e_PMELONG', np.nan)),
                         pmdec=fnum(r.get('PMDEC')) if pm_from_ecl is None else pm_from_ecl[1],
                         e_pmdec=r.get('e_PMDEC', r.get('e_PMELAT', np.nan)),
                         px=px, e_px=epx, dist_a=dA, dist_dm_ymw16=dDM, dist_dm_ne2001=fnum(r.get('DIST_DM1')),
                         dist_best_kpc=d, dist_src=dsrc, type=r.get('TYPE', ''), assoc=r.get('ASSOC', '')))
    return rows


if __name__ == '__main__':
    import csv
    rows = table(parse())
    print('pulsars parsed:', len(rows))
    sel = sorted([r for r in rows if r['sep_silk'] < 30], key=lambda r: r['sep_silk'])
    out = ROOT / 'data/derived/context_pulsars_30deg.csv'
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(sel[0].keys()))
        w.writeheader()
        for r in sel:
            w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    print('within 30 deg:', len(sel), '->', out)
    for r in sel:
        if (np.isfinite(r['dist_best_kpc']) and r['dist_best_kpc'] < 1.0) or 'J0633+1746' in r['psrj']:
            print(f"{r['psrj']:12s} {r['psrb']:9s} sep={r['sep_silk']:5.1f} l={r['l']:6.1f} b={r['b']:+5.1f} "
                  f"d={r['dist_best_kpc']:.3f}({r['dist_src']}) PX={r['px']}+-{r['e_px']} DM={r['dm']} "
                  f"pm=({r['pmra']},{r['pmdec']})+-({r['e_pmra']},{r['e_pmdec']}) age={r['age_yr']:.3g} "
                  f"Edot={r['edot']:.2g} type={r['type']}")
