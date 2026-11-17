"""verify_velocity: LAMOST MRS blue-arm vs red-arm velocity zero point on the paper's plate.

Test: do the LAMOST 2D/1D products that the MRS-N pipeline starts from (arc-lamp wavelength solution,
vacuum, heliocentric correction) carry a blue-minus-red velocity offset on plate MN074209N160100N01
(2024-02-01), and is the heliocentric correction right (absolute RVs vs Gaia DR3)?

Inputs (all downloaded by this agent, logged in data/raw/verify_velocity_lamost/download_log.txt, or read-only
wave-1 files):
  data/raw/verify_velocity_lamost/med_catalogue_rv_MN074209N160100N01.csv   (LAMOST DR12 v2.0 med_catalogue)
  data/raw/verify_velocity_lamost/med_catalogue_rv_2024-01-26_2024-02-07_coaddB_snr15.csv
  data/raw/verify_velocity_lamost/med_catalogue_rv_box_coaddB.csv
  data/raw/datalab_tap/gaia_dr3_rv_box_withrv.csv  (wave-1, Gaia DR3 RVs in the box)
  data/raw/lamost_dr12/mrs/MN074209N160100N01/*.fits.gz (wave-1; headers only)

rv_b0 / rv_r0: cross-correlation RVs (483 Kurucz templates) from the B (4950-5350) and R (6300-6800) arms
separately (LAMOST MRS data-production description).  rv_b1/rv_r1: after per-spectrograph zero-point.

Run: ~/.venvs/astro313/bin/python -I scripts/verify_velocity_lamost_arms.py
Output: data/derived/verify_velocity_lamost_arms.json
"""
import glob
import json
import re
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/verify_velocity_lamost'
OUT = ROOT / 'data/derived/verify_velocity_lamost_arms.json'
PLAN = 'MN074209N160100N01'


def rstat(d, e=None):
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return dict(n=0)
    m = np.median(d)
    mad = 1.4826 * np.median(np.abs(d - m))
    # 3-sigma clipped mean
    k = np.abs(d - m) < 3 * max(mad, 0.3)
    dm = d[k]
    rng = np.random.default_rng(1)
    boot = [np.median(rng.choice(d, len(d))) for _ in range(2000)]
    return dict(n=int(len(d)), median=float(m), mad_sigma=float(mad), clipped_mean=float(dm.mean()),
                clipped_sem=float(dm.std(ddof=1) / np.sqrt(len(dm))) if len(dm) > 1 else None,
                median_boot_sigma=float(np.std(boot)))


def load(f):
    t = Table.read(f, format='ascii.csv')
    t.rename_column(t.colnames[0], t.colnames[0].lstrip('﻿'))
    return t


def good(t):
    m = (t['rv_b0'] > -9000) & (t['rv_r0'] > -9000)
    for c in ('rv_b_flag', 'rv_r_flag'):
        if c in t.colnames:
            m &= (t[c] == 0)
    return t[m]


def main():
    res = {}
    # ---------------------------------------------------------------- MN plate
    mn = load(RAW / f'med_catalogue_rv_{PLAN}.csv')
    mnc = good(mn[(mn['coadd'] == 1) & (mn['band'] == 'B')])
    d = mnc['rv_b0'] - mnc['rv_r0']
    res['MN_plate_coadd_b0_minus_r0_all'] = rstat(d)
    for cut in (5, 10, 15):
        s = mnc[mnc['snr'] > cut]
        res[f'MN_plate_coadd_b0_minus_r0_snr>{cut}'] = rstat(s['rv_b0'] - s['rv_r0'])
    res['MN_plate_coadd_b1_minus_b0'] = rstat(mnc['rv_b1'] - mnc['rv_b0'])
    res['MN_plate_coadd_r1_minus_r0'] = rstat(mnc['rv_r1'] - mnc['rv_r0'])
    # per spectrograph
    persp = {}
    for sp in np.unique(mnc['spid']):
        s = mnc[mnc['spid'] == sp]
        persp[int(sp)] = rstat(s['rv_b0'] - s['rv_r0'])
    res['MN_plate_coadd_b0_minus_r0_per_spid'] = persp
    # single exposures: B rows with coadd=0 carry per-exposure B RV; match R rows by obsid+lmjm
    se = good(mn[(mn['coadd'] == 0) & (mn['band'] == 'B')])
    perexp = {}
    for lm in np.unique(se['lmjm']):
        s = se[se['lmjm'] == lm]
        perexp[str(lm)] = rstat(s['rv_b0'] - s['rv_r0'])
    res['MN_plate_single_exposure_b0_minus_r0'] = perexp

    # ---------------------------------------------------------------- Gaia DR3 absolute check
    g = Table.read(ROOT / 'data/raw/datalab_tap/gaia_dr3_rv_box_withrv.csv', format='ascii.csv')
    gc = SkyCoord(g['ra'] * u.deg, g['dec'] * u.deg)

    def gaia_cmp(t, label):
        c = SkyCoord(t['ra'] * u.deg, t['dec'] * u.deg)
        idx, sep, _ = c.match_to_catalog_sky(gc)
        m = (sep < 2 * u.arcsec) & (g['radial_velocity_error'][idx] < 5)
        out = {}
        for col in ('rv_b0', 'rv_r0', 'rv_b1', 'rv_r1', 'rv_br1'):
            if col in t.colnames:
                out[col + '_minus_gaia'] = rstat(np.asarray(t[col][m]) - np.asarray(g['radial_velocity'][idx][m]))
        res[label] = out
        return m

    gaia_cmp(mnc, 'MN_plate_vs_gaia')

    # ---------------------------------------------------------------- comparison plates
    box = good(load(RAW / 'med_catalogue_rv_box_coaddB.csv'))
    box = box[box['snr'] > 10]
    perplate = []
    for p in np.unique(box['planid']):
        s = box[box['planid'] == p]
        if len(s) < 20:
            continue
        st = rstat(s['rv_b0'] - s['rv_r0'])
        perplate.append(dict(planid=str(p), obsdate=str(s['obsdate'][0]), **st))
    res['box_plates_b0_minus_r0'] = perplate
    med = np.array([q['median'] for q in perplate])
    res['box_plates_summary'] = dict(nplates=len(med), median_of_medians=float(np.median(med)),
                                     min=float(med.min()), max=float(med.max()),
                                     frac_abs_gt_3=float(np.mean(np.abs(med) > 3)))
    gaia_box = {}
    for p in np.unique(box['planid']):
        s = box[box['planid'] == p]
        c = SkyCoord(s['ra'] * u.deg, s['dec'] * u.deg)
        idx, sep, _ = c.match_to_catalog_sky(gc)
        m = (sep < 2 * u.arcsec) & (g['radial_velocity_error'][idx] < 5)
        if m.sum() < 15:
            continue
        gaia_box[str(p)] = dict(
            obsdate=str(s['obsdate'][0]),
            b0=rstat(np.asarray(s['rv_b0'][m]) - np.asarray(g['radial_velocity'][idx][m]))['median'],
            r0=rstat(np.asarray(s['rv_r0'][m]) - np.asarray(g['radial_velocity'][idx][m]))['median'],
            n=int(m.sum()))
    res['box_plates_minus_gaia_median'] = gaia_box

    night = good(load(RAW / 'med_catalogue_rv_2024-01-26_2024-02-07_coaddB_snr15.csv'))
    pn = []
    for p in np.unique(night['planid']):
        s = night[night['planid'] == p]
        st = rstat(s['rv_b0'] - s['rv_r0'])
        pn.append(dict(planid=str(p), obsdate=str(s['obsdate'][0]), **st))
    res['nearby_nights_plates_b0_minus_r0'] = pn

    # ---------------------------------------------------------------- heliocentric correction check
    loc = EarthLocation(lon=117.58 * u.deg, lat=40.39 * u.deg, height=960 * u.m)
    files = sorted(glob.glob(str(ROOT / f'data/raw/lamost_dr12/mrs/{PLAN}/*.fits.gz')))
    hv = []
    for f in files[:83]:
        with fits.open(f) as h:
            H = h[0].header
            hv.append((H['RA'], H['DEC'], H['HELIO_RV'], H['VACUUM'], H['HELIO']))
            exps = [(e.name, e.header['DATE-OBS'], e.header.get('EXPTIME')) for e in h[3:]]
    hv = np.array([(a, b, c) for a, b, c, _, _ in hv], float)
    vc = []
    for name, dobs, ex in exps:
        hh = re.match(r'(\d{4}-\d\d-\d\d)T(\d+):(\d+):(\d+)', dobs).groups()
        t = Time(f'{hh[0]}T{int(hh[1]):02d}:{int(hh[2]):02d}:{int(hh[3]):02d}', scale='utc')
        t = t + (ex or 900) / 2 * u.s
        c = SkyCoord(116.2525 * u.deg, 15.895 * u.deg)
        vc.append(dict(ext=name, mid_utc=t.isot,
                       v_bary=float(c.radial_velocity_correction('barycentric', obstime=t, location=loc).to_value(u.km / u.s))))
    res['helio'] = dict(header_HELIO_RV_median=float(np.median(hv[:, 2])), header_HELIO_RV_range=[float(hv[:, 2].min()), float(hv[:, 2].max())],
                        header_VACUUM=True, header_HELIO=True, exposures=vc,
                        note='LAMOST HELIO_RV = -v_corr(astropy); topocentric-rest features sit at V_helio = v_corr ~ -8.3 km/s')
    OUT.write_text(json.dumps(res, indent=1))
    # short print
    for k in ('MN_plate_coadd_b0_minus_r0_all', 'MN_plate_coadd_b0_minus_r0_snr>10', 'MN_plate_coadd_b1_minus_b0',
              'MN_plate_coadd_r1_minus_r0', 'box_plates_summary'):
        print(k, res[k])
    print('per spid', {k: (v['n'], round(v.get('median', np.nan), 2), round(v.get('mad_sigma', np.nan), 2)) for k, v in persp.items()})
    print('per exp', {k: (v['n'], round(v['median'], 2)) for k, v in perexp.items()})
    print('vs gaia', json.dumps(res['MN_plate_vs_gaia'], indent=0))
    print('helio', res['helio'])
    print('box vs gaia', gaia_box)


if __name__ == '__main__':
    main()
