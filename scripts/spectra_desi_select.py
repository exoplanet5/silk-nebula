"""spectra_desi_select.py -- join SPARCL ids with desi_dr1.zpix metadata and build download selections.

Inputs : data/raw/sparcl_desi_dr1/find_box.csv (SPARCL find), data/raw/datalab_tap/zpix_box.csv (TAP),
         data/derived/silk_trace.json (geometry teammate) via silk_geometry.load_trace().
Outputs: data/derived/spectra_desi_catalog.csv   one row per SPARCL record with zpix columns and
                                                  dRA_spine (deg, + = east of the E spine at that Dec),
                                                  d_spine (great-circle deg to spine polyline), d_any (to E/W/W2)
         data/derived/spectra_desi_sel_<name>.csv selections (pilot, strip, control)
Usage  : python -I scripts/spectra_desi_select.py catalog|pilot|strip <frac_in> <frac_out>|control <frac>
"""
import csv
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import ROOT, REGIONS, load_trace  # noqa: E402

CAT = ROOT / 'data/derived/spectra_desi_catalog.csv'
RNG = np.random.default_rng(20261116)


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def dist_to_polyline(ra, dec, pra, pdec, nsub=20):
    """Great-circle distance (deg) from points to a densely resampled polyline."""
    t = np.linspace(0, 1, nsub, endpoint=False)
    sra = np.concatenate([pra[i] + t * (pra[i + 1] - pra[i]) for i in range(len(pra) - 1)] + [pra[-1:]])
    sdec = np.concatenate([pdec[i] + t * (pdec[i + 1] - pdec[i]) for i in range(len(pdec) - 1)] + [pdec[-1:]])
    P = unitvec(sra, sdec)
    out = np.empty(len(ra))
    X = unitvec(ra, dec)
    for i in range(0, len(ra), 20000):
        d = X[i:i + 20000] @ P.T
        out[i:i + 20000] = np.degrees(np.arccos(np.clip(d.max(1), -1, 1)))
    return out


def trace_arrays():
    t = load_trace()
    sp = t['spine']
    E = (np.array([p['ra'] for p in sp]), np.array([p['dec'] for p in sp]))
    strands = {'E': E}
    for k in ('W', 'W2'):
        if k in t.get('strands', {}):
            pts = t['strands'][k]['points']
            strands[k] = (np.array([p['ra'] for p in pts]), np.array([p['dec'] for p in pts]))
    return t, strands


def build_catalog():
    zp = {}
    with open(ROOT / 'data/raw/datalab_tap/zpix_box.csv') as fh:
        for r in csv.DictReader(fh):
            zp[(r['targetid'], r['survey'], r['program'])] = r
    rows = []
    with open(ROOT / 'data/raw/sparcl_desi_dr1/find_box.csv') as fh:
        for r in csv.DictReader(fh):
            z = zp.get((r['targetid'], r['survey'], r['program']))
            if z is None:
                continue
            rows.append({**{k: r[k] for k in ('sparcl_id', 'targetid', 'survey', 'program', 'spectype',
                                              'redshift')},
                         'ra': z['mean_fiber_ra'], 'dec': z['mean_fiber_dec'],
                         **{k: z[k] for k in ('zwarn', 'deltachi2', 'subtype', 'coadd_numexp', 'coadd_numtile',
                                              'coadd_exptime', 'mean_mjd', 'min_mjd', 'max_mjd',
                                              'mean_psf_to_fiber_specflux', 'healpix', 'coadd_fiberstatus',
                                              'zcat_primary', 'glon', 'glat')}})
    ra = np.array([float(r['ra']) for r in rows])
    dec = np.array([float(r['dec']) for r in rows])
    t, S = trace_arrays()
    dE = dist_to_polyline(ra, dec, *S['E'])
    dany = np.minimum.reduce([dist_to_polyline(ra, dec, *S[k]) for k in S])
    # signed RA offset from the E spine at the same Dec (interpolated; NaN outside spine Dec range)
    sra, sdec = S['E']
    o = np.argsort(sdec)
    ra_sp = np.interp(dec, sdec[o], sra[o], left=np.nan, right=np.nan)
    dra = (ra - ra_sp) * np.cos(np.radians(dec))
    for r, a, b, c in zip(rows, dra, dE, dany):
        r['dRA_spine'] = f'{a:.4f}'
        r['d_spine'] = f'{b:.4f}'
        r['d_any'] = f'{c:.4f}'
    with open(CAT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'catalog: {len(rows)} rows -> {CAT}')


def load_catalog():
    with open(CAT) as fh:
        return list(csv.DictReader(fh))


def write_sel(name, rows):
    p = ROOT / f'data/derived/spectra_desi_sel_{name}.csv'
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'{name}: {len(rows)} -> {p}')


def extragal(r):
    return r['spectype'] in ('GALAXY', 'QSO') and r['zwarn'] == '0'


def sep(ra, dec, ra0, dec0):
    return np.degrees(np.arccos(np.clip(unitvec(ra, dec) @ unitvec(ra0, dec0), -1, 1)))


def pilot():
    """E-W strip through P1 and P4 (Dec 16.10-16.70, RA 114.6-118.6) + discs r=0.2 deg on P2, P3."""
    rows = load_catalog()
    out = []
    for r in rows:
        if not extragal(r):
            continue
        ra, dec = float(r['ra']), float(r['dec'])
        instrip = (16.10 <= dec <= 16.70) and (114.6 <= ra <= 118.6)
        indisc = any(sep(ra, dec, *REGIONS[k]) < 0.2 for k in ('P2', 'P3'))
        if instrip or indisc:
            out.append(r)
    write_sel('pilot', out)


def strip(frac_in, frac_out, exclude=('pilot',)):
    """Within 3 deg of any strand: keep frac_in of extragalactic spectra with d_any<1 deg, frac_out for 1-3 deg."""
    rows = load_catalog()
    skip = set()
    for e in exclude:
        p = ROOT / f'data/derived/spectra_desi_sel_{e}.csv'
        if p.exists():
            skip.update(r['sparcl_id'] for r in csv.DictReader(open(p)))
    out = []
    for r in rows:
        if not extragal(r) or r['sparcl_id'] in skip:
            continue
        d = float(r['d_any'])
        if d < 1.0 and RNG.random() < frac_in:
            out.append(r)
        elif 1.0 <= d < 3.0 and RNG.random() < frac_out:
            out.append(r)
    write_sel('strip', out)


def control(frac):
    """Controls > 5 deg from every strand (within the box)."""
    rows = load_catalog()
    out = [r for r in rows if extragal(r) and float(r['d_any']) > 5.0 and RNG.random() < frac]
    write_sel('control', out)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'catalog':
        build_catalog()
    elif cmd == 'pilot':
        pilot()
    elif cmd == 'strip':
        strip(float(sys.argv[2]), float(sys.argv[3]))
    elif cmd == 'control':
        control(float(sys.argv[2]))
