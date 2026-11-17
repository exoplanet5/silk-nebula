"""spectra_desi_rvcheck.py -- DESI DR1 velocity zero-point check: MWS VAC VRAD vs Gaia DR3 radial_velocity.

Purpose: the DESI [O III] residual velocity is compared with LAMOST MRS. This checks that DESI DR1 spectra
(barycentric vacuum wavelengths) carry no large velocity zero-point error in this field.

Inputs (Data Lab TAP, https://datalab.noirlab.edu/tap/sync, accessed 2026-11-16; see spectra_desi_tap.py):
  data/raw/datalab_tap/mws_box.csv              desi_dr1.mws: targetid, target_ra, target_dec, vrad, vrad_err,
                                                rvs_warn, sn_b, teff, logg, feh  (RA 111-122, Dec +6..+28)
  data/raw/datalab_tap/gaia_dr3_rv_box_withrv.csv  gaia_dr3.gaia_source rows with radial_velocity in the same box
Output: data/derived/spectra_desi_rvcheck.json
Usage : python -I scripts/spectra_desi_rvcheck.py
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/datalab_tap'


def xyz(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def load(p, cols):
    rows = list(csv.DictReader(open(p)))
    out = {}
    for c in cols:
        v = [r[c] for r in rows]
        out[c] = np.array([float(x) if x not in ('', 'nan', 'None') else np.nan for x in v])
    return out


def robust(x):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return float(med), float(mad), float(1.2533 * mad / np.sqrt(len(x)))


def main():
    m = load(RAW / 'mws_box.csv', ['target_ra', 'target_dec', 'vrad', 'vrad_err', 'rvs_warn', 'sn_b', 'teff'])
    g = load(RAW / 'gaia_dr3_rv_box_withrv.csv', ['ra', 'dec', 'radial_velocity', 'radial_velocity_error',
                                                  'phot_g_mean_mag', 'rv_template_teff'])
    tree = cKDTree(xyz(g['ra'], g['dec']))
    d, j = tree.query(xyz(m['target_ra'], m['target_dec']), distance_upper_bound=2 * np.sin(np.radians(1 / 3600) / 2))
    ok = np.isfinite(d) & (j < len(g['ra']))
    jj = j[ok]
    dv = m['vrad'][ok] - g['radial_velocity'][jj]
    edv = np.hypot(m['vrad_err'][ok], g['radial_velocity_error'][jj])
    G = g['phot_g_mean_mag'][jj]
    good = (m['rvs_warn'][ok] == 0) & np.isfinite(dv) & (m['vrad_err'][ok] < 2) & \
        (g['radial_velocity_error'][jj] < 3) & (np.abs(dv) < 50) & (m['sn_b'][ok] > 10)
    out = {'n_mws': int(len(m['vrad'])), 'n_gaia_rv': int(len(g['ra'])), 'n_matched_1arcsec': int(ok.sum()),
           'n_good': int(good.sum())}
    med, sig, emed = robust(dv[good])
    out['dv_DESI_minus_Gaia_all'] = dict(median=med, mad_sigma=sig, err_median=emed)
    bins = {}
    for lo, hi in ((8, 11), (11, 12), (12, 13), (13, 14), (14, 15), (15, 16)):
        s = good & (G >= lo) & (G < hi)
        if s.sum() > 20:
            a, b, c = robust(dv[s])
            bins[f'G{lo}-{hi}'] = dict(n=int(s.sum()), median=a, mad_sigma=b, err_median=c)
    out['dv_by_G'] = bins
    bright = good & (G < 12)
    if bright.sum() > 20:
        a, b, c = robust(dv[bright])
        out['dv_G_lt_12'] = dict(n=int(bright.sum()), median=a, mad_sigma=b, err_median=c,
                                 note='Gaia DR3 RV zero point is best for G_RVS < 11 (Katz et al. 2023)')
    out['note'] = ('DESI MWS VAC VRAD (RVSpecFit, barycentric) minus Gaia DR3 radial_velocity for stars matched '
                   'within 1 arcsec; rvs_warn==0, vrad_err<2, Gaia err<3, |dv|<50 km/s, sn_b>10')
    (ROOT / 'data/derived/spectra_desi_rvcheck.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()


def three_way():
    """DESI zpix STAR cz (by survey/program) vs LAMOST LRS cz (teammate download, VizieR/LAMOST DR12 catalogue in
    data/raw/lamost_dr12/lrs_catalogue_strandbox.csv) vs LAMOST MRS rv_br1 (med_catalogue_box.csv) vs Gaia DR3."""
    C = 299792.458
    zp = [r for r in csv.DictReader(open(RAW / 'zpix_box.csv')) if r['spectype'] == 'STAR' and r['zwarn'] == '0']
    D = dict(ra=np.array([float(r['mean_fiber_ra']) for r in zp]), dec=np.array([float(r['mean_fiber_dec']) for r in zp]),
             v=np.array([float(r['z']) * C for r in zp]), prog=np.array([r['survey'] + '/' + r['program'] for r in zp]))
    L = ROOT / 'data/raw/lamost_dr12'
    lr = [r for r in csv.DictReader(open(L / 'lrs_catalogue_strandbox.csv', encoding='utf-8-sig'))
          if r['class'] == 'STAR' and float(r['snrg']) > 15]
    LR = dict(ra=np.array([float(r['ra']) for r in lr]), dec=np.array([float(r['dec']) for r in lr]),
              v=np.array([float(r['z']) * C for r in lr]))
    mr = [r for r in csv.DictReader(open(L / 'med_catalogue_box.csv', encoding='utf-8-sig'))
          if r['band'] == 'B' and float(r['rv_br1']) > -9000 and 0 < float(r['rv_br1_err']) < 2 and float(r['snr']) > 15]
    MR = dict(ra=np.array([float(r['ra']) for r in mr]), dec=np.array([float(r['dec']) for r in mr]),
              v=np.array([float(r['rv_br1']) for r in mr]))
    g = load(RAW / 'gaia_dr3_rv_box_withrv.csv', ['ra', 'dec', 'radial_velocity', 'radial_velocity_error'])
    G = dict(ra=g['ra'], dec=g['dec'], v=g['radial_velocity'])
    G_ok = g['radial_velocity_error'] < 3
    G = {k: v[G_ok] for k, v in G.items()}
    r1 = 2 * np.sin(np.radians(1.5 / 3600) / 2)

    def match(A, B, cut=60):
        t = cKDTree(xyz(B['ra'], B['dec']))
        d, j = t.query(xyz(A['ra'], A['dec']), distance_upper_bound=r1)
        ok = np.isfinite(d) & (j < len(B['ra']))
        dv = A['v'][ok] - B['v'][j[ok]]
        s = np.abs(dv) < cut
        return ok, dv, s

    out = {}
    for p in np.unique(D['prog']):
        sel = D['prog'] == p
        A = {k: v[sel] for k, v in D.items()}
        for nm, B in (('LRS', LR), ('MRS', MR), ('Gaia', G)):
            ok, dv, s = match(A, B)
            if s.sum() >= 5:
                a, b, c = robust(dv[s])
                out[f'DESI_{p}_minus_{nm}'] = dict(n=int(s.sum()), median=a, mad_sigma=b, err_median=c)
    for nm, A in (('LRS', LR), ('MRS', MR)):
        ok, dv, s = match(A, G)
        if s.sum() >= 5:
            a, b, c = robust(dv[s])
            out[f'{nm}_minus_Gaia'] = dict(n=int(s.sum()), median=a, mad_sigma=b, err_median=c)
    out['note'] = ('cz differences (km/s) of stars matched within 1.5 arcsec; |dv|<60; LAMOST LRS snrg>15 (z*c), '
                   'MRS blue-band rv_br1 with err<2 and snr>15; Gaia err<3')
    return out


if __name__ == '__main__':
    o = three_way()
    p = ROOT / 'data/derived/spectra_desi_rvcheck.json'
    j = json.loads(p.read_text())
    j['three_way'] = o
    p.write_text(json.dumps(j, indent=1))
    print(json.dumps(o, indent=1))
