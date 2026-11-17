#!/usr/bin/env python
"""spectra_ism: inventory of archival high-resolution stellar spectra toward early-type stars near the
Silk Nebula, for interstellar Ca II K/H and Na I D absorption.

Queries (all logged to data/raw/spectra_ism_fetch_log.jsonl):
  * SIMBAD TAP (CfA mirror)  https://simbad.cfa.harvard.edu/simbad/sim-tap/sync : O/B/A stars, V < 10
  * ESO ObsCore TAP          https://archive.eso.org/tap_obs/sync : phase-3 1D spectra (UVES, FEROS, HARPS,
                             ESPRESSO, UVES_SQUAD)
  * ELODIE archive (OHP)     http://atlas.obs-hp.fr/elodie/fE.cgi  region listing, a=csv
  * SOPHIE archive (OHP)     http://atlas.obs-hp.fr/sophie/sophie.cgi region listing, a=csv
Box: RA 104-128 deg, Dec +4..+34 deg (covers stars within ~10 deg of the [O III] spine).
Output: data/raw/spectra_ism_archives/*, data/derived/spectra_ism_archive_candidates.csv
"""
import csv
import io
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, fetch, filament_coords, spine_arrays  # noqa: E402

OUT = RAW / 'spectra_ism_archives'
OUT.mkdir(parents=True, exist_ok=True)
BOX = (104.0, 128.0, 4.0, 34.0)


def tap(url, query, dest):
    ok, n, st = fetch(url, dest, params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query},
                      timeout=300, retries=3)
    if not ok:
        print('TAP failed', url, st)
        return []
    return list(csv.DictReader(io.StringIO(dest.read_text(errors='replace'))))


def run_queries():
    if '--cached' in sys.argv:
        rd = lambda p: list(csv.DictReader(io.StringIO((OUT / p).read_text(errors='replace'))))  # noqa: E731
        return rd('simbad_OBA_V10_box.csv'), rd('eso_obscore_box.csv')
    # SIMBAD early-type stars
    q = f"""SELECT b.main_id, b.ra, b.dec, b.sp_type, b.plx_value, b.plx_err, b.rvz_radvel, b.otype, f.V AS vmag
FROM basic AS b JOIN allfluxes AS f ON b.oid = f.oidref
WHERE b.ra BETWEEN {BOX[0]} AND {BOX[1]} AND b.dec BETWEEN {BOX[2]} AND {BOX[3]}
AND (b.sp_type LIKE 'O%' OR b.sp_type LIKE 'B%' OR b.sp_type LIKE 'A%') AND f.V < 10"""
    sim = tap('https://simbad.cfa.harvard.edu/simbad/sim-tap/sync', q, OUT / 'simbad_OBA_V10_box.csv')
    print('SIMBAD O/B/A V<10:', len(sim))
    # ESO phase 3 spectra
    q = f"""SELECT obs_collection, instrument_name, target_name, s_ra, s_dec, em_min, em_max, em_res_power, snr,
t_exptime, t_min, dp_id, access_url FROM ivoa.ObsCore
WHERE s_ra BETWEEN {BOX[0]} AND {BOX[1]} AND s_dec BETWEEN {BOX[2]} AND {BOX[3]} AND dataproduct_type='spectrum'
AND obs_collection IN ('UVES','FEROS','HARPS','ESPRESSO','UVES_SQUAD')"""
    eso = tap('https://archive.eso.org/tap_obs/sync', q, OUT / 'eso_obscore_box.csv')
    print('ESO phase-3 spectra:', len(eso))
    # ELODIE / SOPHIE region listings (RA in hours)
    fetch('http://atlas.obs-hp.fr/elodie/fE.cgi', OUT / 'elodie_region.txt',
          params={'ob': 'objname,dataset,imanum', 'r': '[06 56 +04],[08 32 +34]', 'a': 'csv'})
    fetch('http://atlas.obs-hp.fr/sophie/sophie.cgi', OUT / 'sophie_region.txt',
          params={'n': 'sophies', 'ob': 'ra,seq', 'r': '[06 56 +04],[08 32 +34]', 'a': 'csv'})
    return sim, eso


def parse_ohp(path):
    """Parse the OHP 'a=csv' listing (tab separated, '#' comments, '$'-prefixed header)."""
    rows = []
    if not path.exists():
        return rows
    hdr = None
    for ln in path.read_text(errors='replace').splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        parts = ln.split('\t')
        if ln.startswith('$'):
            hdr = [p.strip().replace('$link[', '').replace(']', '').replace("'", '').lstrip('$') for p in ln.split()]
            continue
        if hdr is None:
            continue
        rows.append(dict(zip(hdr, [p.strip() for p in parts])))
    return rows


def j2000_to_deg(s):
    s = s.strip().lstrip('J')
    try:
        c = SkyCoord(f'{s[0:2]}h{s[2:4]}m{s[4:8]}s', f'{s[8:11]}d{s[11:13]}m{s[13:]}s')
        return c.ra.deg, c.dec.deg
    except Exception:
        return np.nan, np.nan


if __name__ == '__main__':
    sim, eso = run_queries()
    elo = parse_ohp(OUT / 'elodie_region.txt')
    sop = parse_ohp(OUT / 'sophie_region.txt')
    print('ELODIE rows', len(elo), 'SOPHIE rows', len(sop))
    if sop:
        print('SOPHIE columns:', list(sop[0].keys())[:20])

    sra = np.array([float(r['ra']) for r in sim])
    sde = np.array([float(r['dec']) for r in sim])
    csim = SkyCoord(sra * u.deg, sde * u.deg)
    dperp, dsp, beyond = filament_coords(sra, sde)
    spra, spdec, _ = spine_arrays()
    cspine = SkyCoord(np.interp(np.linspace(0, 1, 800), np.linspace(0, 1, len(spra)), spra) * u.deg,
                      np.interp(np.linspace(0, 1, 800), np.linspace(0, 1, len(spdec)), spdec) * u.deg)
    _, dmin, _ = csim.match_to_catalog_sky(cspine)

    def match(rows, rakey, deckey, conv=None, rad=10.0):
        if not rows:
            return [[] for _ in sim]
        if conv:
            pos = np.array([conv(r) for r in rows], float)
        else:
            pos = np.array([[float(r[rakey] or 'nan'), float(r[deckey] or 'nan')] for r in rows])
        good = np.isfinite(pos).all(1)
        c = SkyCoord(pos[good, 0] * u.deg, pos[good, 1] * u.deg)
        idx = np.where(good)[0]
        out = [[] for _ in sim]
        i1, i2, sep, _ = c.search_around_sky(csim, rad * u.arcsec)
        for a, b in zip(i1, i2):
            out[a].append(rows[idx[b]])
        return out

    m_eso = match(eso, 's_ra', 's_dec')
    m_elo = match(elo, None, None, conv=lambda r: j2000_to_deg(r.get('j2000', '')), rad=30.0)

    m_sop = match(sop, None, None, conv=lambda r: j2000_to_deg(r.get('j2000', '')), rad=30.0) if sop else [[] for _ in sim]

    rows = []
    for i, r in enumerate(sim):
        plx = float(r['plx_value']) if r['plx_value'] else np.nan
        ne, nl, ns = len(m_eso[i]), len(m_elo[i]), len(m_sop[i])
        best_eso = sorted(m_eso[i], key=lambda x: -float(x['snr'] or 0))[:1]
        best_elo = sorted(m_elo[i], key=lambda x: -float(x.get('sn', 0) or 0))[:1]
        rows.append(dict(name=r['main_id'], ra=round(sra[i], 5), dec=round(sde[i], 5), sp=r['sp_type'], V=r['vmag'],
                         plx=plx, dist_pc=round(1000 / plx, 1) if plx > 0 else '', dperp_deg=round(float(dperp[i]), 3),
                         dmin_spine_deg=round(float(dmin.deg[i]), 3), dec_spine=round(float(dsp[i]), 2),
                         n_eso=ne, eso_best=(f"{best_eso[0]['obs_collection']}|{best_eso[0]['dp_id']}|snr={best_eso[0]['snr']}"
                                             if best_eso else ''),
                         n_elodie=nl, elodie_best=(f"{best_elo[0].get('dataset')}/{best_elo[0].get('imanum')}|sn={best_elo[0].get('sn')}"
                                                   if best_elo else ''),
                         n_sophie=ns,
                         sophie_best=(lambda b: f"{b[0].get('seq')}|{b[0].get('mode')}|sn26={b[0].get('sn26')}" if b else '')(
                             sorted(m_sop[i], key=lambda x: -float(x.get('sn26', 0) or 0))[:1])))
    rows.sort(key=lambda x: x['dmin_spine_deg'])
    with open(DER / 'spectra_ism_archive_candidates.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    sel = [r for r in rows if (r['n_eso'] or r['n_elodie'] or r['n_sophie']) and r['dist_pc'] != '' and r['dist_pc'] < 1500]
    print('O/B/A stars with archival high-res spectra and d<1.5 kpc:', len(sel))
    for r in sel[:80]:
        print(f"{r['name']:22s} {r['sp']:10s} V={r['V']:>5s} d={r['dist_pc']:7} dmin={r['dmin_spine_deg']:5.2f} "
              f"dperp={r['dperp_deg']:+6.2f} ESO={r['n_eso']:3d} {r['eso_best'][:45]:45s} ELO={r['n_elodie']} {r['elodie_best']} SOP={r['n_sophie']} {r['sophie_best']}")
