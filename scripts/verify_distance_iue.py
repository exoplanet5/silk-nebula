"""verify_distance (b): IUE SWP high-dispersion (NEWSIPS MXHI) absorption toward hot stars behind / around the Silk.

Stars (SIMBAD CfA mirror, data/raw/verify_distance_meta/simbad_uv_candidates.csv; Gaia DR3 parallaxes):
  HD 60848  O8:V:pe  plx 0.523 mas (~1.9 kpc)  2.26 deg W of the spine (concave/Monogem side), 13 SWP-HI
  NGC 2392  O6f CSPN plx 0.545 (~1.8 kpc)      4.0 deg W                                         7 SWP-HI
  NGC 2371  CSPN     plx 0.581 (~1.7 kpc)      2.9 deg W (N end)                                 1 SWP-HI
  HD 58050  B2Ve     plx 2.21 (452 pc)         5.3 deg W                                         4 SWP-HI
  HD 52559  B2IV-V   plx 1.23 (811 pc)         12.7 deg W (Monogem interior)                     2 SWP-HI
  HD 51354  B3ne     plx 1.79 (559 pc)         11.6 deg W (Monogem interior, near centre)        2 SWP-HI
  HD 53929  B9III    plx 4.45 (225 pc)         12.0 deg W                                        2 SWP-HI
  HD 66665  B1V      plx 0.563 (~1.8 kpc)      7.4 deg E (convex side, outside the arc) = control 3 SWP-HI
Files: data/raw/verify_distance_iue/*.mxhi.gz (fetch log in the same folder).

Method: per spectrum and line, take every echelle order containing the line (away from order ends), normalise by
a quadratic continuum fitted at |v| = 90..300 km/s (C IV/N V: 90..220), co-add orders and spectra with inverse-
variance weights after aligning each spectrum on the mean low-ion (S II 1250/1253/1259, Si II 1304/1526,
Al II 1670, Fe II 1608) centroid of that star.  Apparent-optical-depth (AOD) columns and EWs in velocity windows,
error from the propagated NEWSIPS noise inflated by the empirical continuum scatter.  IUE LSF FWHM ~ 25 km/s.
Outputs: data/derived/verify_distance_iue.json, data/derived/verify_distance_iue_profiles.npz
Run: ~/.venvs/astro313/bin/python -I scripts/verify_distance_iue.py
"""
import glob
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_distance_iuelib import read_mxhi, extract, UVLINES, C  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula'
ALIAS = {'HD60848': 'HD 60848', 'NGC2392': 'NGC 2392', 'HR 2633': 'HD 52559', 'HR 2817': 'HD 58050', 'NGC2371': 'NGC 2371'}
STARINFO = {
    'HD 60848': dict(sp='O8:V:pe', plx=0.5227, eplx=0.0548, rv=5.47, dperp=-2.26, role='behind, W (Monogem side)'),
    'NGC 2392': dict(sp='O6f CSPN', plx=0.5451, eplx=0.0467, rv=84.2, dperp=-4.04, role='behind, W'),
    'NGC 2371': dict(sp='CSPN', plx=0.5812, eplx=0.0459, rv=21.0, dperp=-2.88, role='behind, W, N end'),
    'HD 58050': dict(sp='B2Ve', plx=2.2123, eplx=0.0542, rv=36.9, dperp=-5.28, role='W'),
    'HD 52559': dict(sp='B2IV-V', plx=1.2329, eplx=0.0414, rv=106.7, dperp=-12.68, role='Monogem interior'),
    'HD 51354': dict(sp='B3ne', plx=1.7888, eplx=0.0539, rv=None, dperp=-11.61, role='Monogem interior'),
    'HD 53929': dict(sp='B9III', plx=4.4509, eplx=0.0495, rv=11.2, dperp=-11.97, role='Monogem interior, near'),
    'HD 66665': dict(sp='B1V', plx=0.5632, eplx=0.0635, rv=16.9, dperp=7.41, role='control, E (outside arc)'),
}
LOW = ['SII1250', 'SII1253', 'SII1259', 'SiII1304', 'SiII1526', 'AlII1670', 'FeII1608']
HIGH = ['SiIV1393', 'SiIV1402', 'CIV1548', 'CIV1550', 'NV1238', 'NV1242', 'AlIII1854', 'AlIII1862']
ALL = LOW + HIGH + ['CII1334', 'SiII1260', 'OI1302']
VG = np.arange(-400, 400.01, 2.0)
WINDOWS = {'core_m50_p60': (-50, 60), 'local_m10_p20': (-10, 20), 'silk_p5_p30': (5, 30),
           'neg_m60_m20': (-60, -20), 'hv_neg_m200_m70': (-200, -70), 'hv_pos_p70_p200': (70, 200)}


def norm_line(ex, lam0, line):
    """Return normalised flux & error on VG for one spectrum (inverse-variance co-add of its orders)."""
    num = np.zeros_like(VG)
    den = np.zeros_like(VG)
    inner = 220.0 if line.startswith(('CIV', 'NV')) else 300.0
    for e in ex:
        good = (e['q'] >= 0) & np.isfinite(e['f']) & np.isfinite(e['e']) & (e['e'] > 0)
        v, f, er = e['v'][good], e['f'][good], e['e'][good]
        if len(v) < 40:
            continue
        cm = (np.abs(v) > 90) & (np.abs(v) < inner)
        if cm.sum() < 15:
            continue
        p = np.polyfit(v[cm] / 100, f[cm], 2)
        cont = np.polyval(p, v / 100)
        fn, en = f / cont, er / np.abs(cont)
        scat = np.std(fn[cm] - 1)
        en = np.sqrt(en ** 2 + 0 * scat ** 2)
        fi = np.interp(VG, v, fn, left=np.nan, right=np.nan)
        ei = np.interp(VG, v, en, left=np.nan, right=np.nan)
        # empirical noise: max(pipeline, continuum scatter)
        ei = np.maximum(ei, scat)
        w = np.where(np.isfinite(fi) & np.isfinite(ei) & (ei > 0), 1 / ei ** 2, 0)
        num += np.nan_to_num(fi) * w
        den += w
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / den, 1 / np.sqrt(den)


def centroid(fn, en, lo=-60, hi=70):
    m = (VG > lo) & (VG < hi) & np.isfinite(fn)
    a = np.clip(1 - fn[m], 0, None)
    if a.sum() <= 0:
        return np.nan, np.nan
    vc = np.sum(VG[m] * a) / a.sum()
    ev = np.sqrt(np.sum(((VG[m] - vc) * en[m]) ** 2)) / a.sum()
    return vc, ev


def aod(fn, en, lam, fval, lo, hi):
    m = (VG >= lo) & (VG <= hi) & np.isfinite(fn)
    if m.sum() < 3:
        return dict(W=np.nan, eW=np.nan, logN=np.nan, logN_3sig=np.nan)
    x = np.clip(fn[m], 0.03, None)
    dv = VG[1] - VG[0]
    W = np.sum(1 - fn[m]) * dv / C * lam * 1000
    eW = np.sqrt(np.sum(en[m] ** 2)) * dv / C * lam * 1000 * np.sqrt(4.0)   # ~4 resampled px per IUE resel.
    tau = -np.log(x)
    N = 3.768e14 / (fval * lam) * np.sum(tau) * dv
    eN = 3.768e14 / (fval * lam) * np.sqrt(np.sum((en[m] / x) ** 2)) * dv * 2.0
    N3 = 1.13e17 * (max(W, 0) + 3 * eW) / (fval * lam ** 2)
    sat = bool(np.nanmin(fn[m]) < 0.15)
    return dict(W=round(float(W), 1), eW=round(float(eW), 1), logN=round(float(np.log10(N)), 2) if N > 0 else None,
                elogN=round(float(eN / N / np.log(10)), 2) if N > 0 else None,
                logN_3sig_upper=round(float(np.log10(N3)), 2), saturated=sat)


files = sorted(glob.glob(f'{ROOT}/data/raw/verify_distance_iue/*.mxhi.gz'))
by_star = {}
for f in files:
    orders, meta = read_mxhi(f)
    name = ALIAS.get(meta['target'], meta['target'])
    by_star.setdefault(name, []).append((f.split('/')[-1].split('.')[0], orders, meta))

OUT = {'stars': {}, 'windows_kms': WINDOWS, 'note': 'velocities heliocentric, NEWSIPS zero point; IUE resolution ~25 km/s'}
PROF = {}
for star, lst in by_star.items():
    # per-spectrum normalised profiles, low-ion centroid per spectrum
    per = []
    for sid, orders, meta in lst:
        prof = {}
        for ln in ALL:
            lam, fv = UVLINES[ln]
            ex = extract(orders, lam, vwin=450)
            if ex:
                prof[ln] = norm_line(ex, lam, ln)
        cents = []
        for ln in LOW:
            if ln in prof:
                vc, ev = centroid(*prof[ln], lo=-50, hi=60)
                if np.isfinite(vc) and np.isfinite(ev) and ev < 6:
                    cents.append((vc, ev))
        vlow = (np.average([c[0] for c in cents], weights=[1 / c[1] ** 2 for c in cents])
                if cents else np.nan)
        per.append(dict(id=sid, aperture=meta['aperture'], jd=meta['jd'], radvelo=meta['radvelo'], vlow=vlow,
                        nlow=len(cents), prof=prof))
    vlows = np.array([p['vlow'] for p in per])
    vref = np.nanmedian(vlows)
    # co-add with alignment to the median low-ion centroid
    coadd = {}
    for ln in ALL:
        num = np.zeros_like(VG)
        den = np.zeros_like(VG)
        for p in per:
            if ln not in p['prof']:
                continue
            fn, en = p['prof'][ln]
            sh = (p['vlow'] - vref) if np.isfinite(p['vlow']) else 0.0
            fi = np.interp(VG, VG - sh, fn, left=np.nan, right=np.nan)
            ei = np.interp(VG, VG - sh, en, left=np.nan, right=np.nan)
            w = np.where(np.isfinite(fi) & np.isfinite(ei) & (ei > 0), 1 / ei ** 2, 0)
            num += np.nan_to_num(fi) * w
            den += w
        if den.max() > 0:
            with np.errstate(invalid='ignore', divide='ignore'):
                coadd[ln] = (num / den, 1 / np.sqrt(den))
    res = {}
    for ln, (fn, en) in coadd.items():
        lam, fv = UVLINES[ln]
        r = {k: aod(fn, en, lam, fv, *w) for k, w in WINDOWS.items()}
        vc, ev = centroid(fn, en, -50, 60)
        r['centroid_m50_p60'] = [round(float(vc), 1), round(float(ev), 1)]
        r['median_err_cont'] = round(float(np.nanmedian(en[(np.abs(VG) > 90) & (np.abs(VG) < 200)])), 3)
        res[ln] = r
        PROF[f'{star}|{ln}|f'] = fn
        PROF[f'{star}|{ln}|e'] = en
    info = STARINFO.get(star, {})
    d = 1000 / info['plx'] if info.get('plx') else None
    OUT['stars'][star] = dict(info=info, dist_pc=round(d) if d else None, n_spectra=len(per),
                              spectra=[dict(id=p['id'], ap=p['aperture'], jd=p['jd'], vlow=round(float(p['vlow']), 1)
                                            if np.isfinite(p['vlow']) else None, nlow=p['nlow']) for p in per],
                              vlow_ref=round(float(vref), 1) if np.isfinite(vref) else None,
                              vlow_scatter=round(float(np.nanstd(vlows)), 1), lines=res)
    print(star, 'n', len(per), 'vlow', np.round(vlows, 1), 'ref', round(float(vref), 1))
    for ln in ('SII1253', 'SiII1304', 'AlIII1854', 'SiIV1393', 'SiIV1402', 'CIV1548', 'CIV1550', 'NV1238'):
        if ln in res:
            r = res[ln]
            print(f"   {ln:10s} core W={r['core_m50_p60']['W']}+-{r['core_m50_p60']['eW']} logN={r['core_m50_p60']['logN']}"
                  f"+-{r['core_m50_p60']['elogN']} (<{r['core_m50_p60']['logN_3sig_upper']}) cen={r['centroid_m50_p60']}"
                  f" hv-={r['hv_neg_m200_m70']['W']}+-{r['hv_neg_m200_m70']['eW']} hv+={r['hv_pos_p70_p200']['W']}+-{r['hv_pos_p70_p200']['eW']}")

# derived: high-minus-low centroid offsets and per-kpc columns
for star, s in OUT['stars'].items():
    L = s['lines']
    off = {}
    for hi in ('SiIV1393', 'SiIV1402', 'CIV1548', 'CIV1550', 'AlIII1854'):
        if hi in L and 'SII1253' in L:
            a, ea = L[hi]['centroid_m50_p60']
            b, eb = L['SII1253']['centroid_m50_p60']
            off[hi + '-SII1253'] = [round(a - b, 1), round(float(np.hypot(ea, eb)), 1)]
    s['centroid_offsets'] = off
    if s['dist_pc']:
        pk = {}
        for hi in ('SiIV1393', 'CIV1548', 'NV1238', 'AlIII1854'):
            if hi in L and L[hi]['core_m50_p60']['logN']:
                pk[hi] = round(float(L[hi]['core_m50_p60']['logN'] - np.log10(s['dist_pc'] / 1000)), 2)
        s['logN_per_kpc_core'] = pk

np.savez_compressed(f'{ROOT}/data/derived/verify_distance_iue_profiles.npz', VG=VG, **{k.replace(' ', '_'): v for k, v in PROF.items()})
open(f'{ROOT}/data/derived/verify_distance_iue.json', 'w').write(json.dumps(OUT, indent=1, default=float))
print('saved')

# ------------------------------------------------------------------ refined doublet fits for the clean O-star sightline(s)
from scipy.optimize import least_squares  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402
FINE = np.arange(-400, 400.01, 0.5)
LSF_FWHM = 25.0


def doublet_fit(star, l1, l2, cont_in=50.0, cont_out=170.0, nboot=300, seed=1):
    """Joint fit of a doublet: common v, b, logN; local linear continua (independent per line) re-fitted from
    the raw co-add over cont_in<|v|<cont_out.  Errors by bootstrap over resolution elements of the residuals."""
    rng = np.random.default_rng(seed)
    f1, e1 = PROF[f'{star}|{l1}|f'], PROF[f'{star}|{l1}|e']
    f2, e2 = PROF[f'{star}|{l2}|f'], PROF[f'{star}|{l2}|e']
    m = (np.abs(VG) < cont_out) & np.isfinite(f1) & np.isfinite(f2)
    v = VG[m]
    y1, y2, s1, s2 = f1[m], f2[m], e1[m], e2[m]
    (lam1, fv1), (lam2, fv2) = UVLINES[l1], UVLINES[l2]
    sig = LSF_FWHM / 2.3548 / 0.5

    def model(p, v=v):
        lN, b, v0, a1, b1, a2, b2 = p
        t1 = 1.4974e-15 * 10 ** lN * fv1 * lam1 / b * np.exp(-((FINE - v0) / b) ** 2)
        t2 = 1.4974e-15 * 10 ** lN * fv2 * lam2 / b * np.exp(-((FINE - v0) / b) ** 2)
        g1 = np.interp(v, FINE, gaussian_filter1d(np.exp(-t1), sig))
        g2 = np.interp(v, FINE, gaussian_filter1d(np.exp(-t2), sig))
        return (a1 + b1 * v / 100) * g1, (a2 + b2 * v / 100) * g2

    def res(p, y1=y1, y2=y2):
        m1, m2 = model(p)
        return np.concatenate([(m1 - y1) / s1, (m2 - y2) / s2])
    p0 = [13.0, 15.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    lo = [11.0, 3.0, -60, 0.5, -1, 0.5, -1]
    hi = [15.5, 60.0, 60, 1.5, 1, 1.5, 1]
    r = least_squares(res, p0, bounds=(lo, hi))
    m1, m2 = model(r.x)
    chi2r = 2 * r.cost / (2 * len(v) - 7)
    boots = []
    r1, r2 = y1 - m1, y2 - m2
    L = int(round(LSF_FWHM / 2.0))
    for _ in range(nboot):
        idx = rng.integers(0, len(v) // L, len(v) // L + 1)
        bb1 = np.concatenate([r1[i * L:(i + 1) * L] for i in idx])[:len(v)]
        bb2 = np.concatenate([r2[i * L:(i + 1) * L] for i in idx])[:len(v)]
        rb = least_squares(res, r.x, bounds=(lo, hi), kwargs=dict(y1=m1 + bb1, y2=m2 + bb2))
        boots.append(rb.x[:3])
    boots = np.array(boots)
    return dict(logN=round(float(r.x[0]), 2), e_logN=round(float(np.std(boots[:, 0])), 2),
                b_kms=round(float(r.x[1]), 1), e_b=round(float(np.std(boots[:, 1])), 1),
                v_helio=round(float(r.x[2]), 1), e_v=round(float(np.std(boots[:, 2])), 1),
                chi2_red=round(float(chi2r), 2), n_pix=int(len(v)))


REFINED = {}
for star in ('HD 60848', 'HD 58050'):
    REFINED[star] = {}
    for l1, l2 in (('SiIV1393', 'SiIV1402'), ('CIV1548', 'CIV1550'), ('NV1238', 'NV1242'), ('AlIII1854', 'AlIII1862'),
                   ('SII1253', 'SII1259')):
        if f'{star}|{l1}|f' in PROF and f'{star}|{l2}|f' in PROF:
            try:
                REFINED[star][l1[:-4]] = doublet_fit(star, l1, l2)
            except Exception as ex:  # noqa: BLE001
                REFINED[star][l1[:-4]] = str(ex)
    print(star, json.dumps(REFINED[star]))
OUT['refined_doublet_fits'] = REFINED
OUT['refined_note'] = ('joint doublet fits (common logN, b, v; Gaussian tau convolved with 25 km/s IUE LSF; local linear '
                       'continua over |v|<170 km/s); S II uses 1253+1259; N V likely upper limit (check logN error)')
open(f'{ROOT}/data/derived/verify_distance_iue.json', 'w').write(json.dumps(OUT, indent=1, default=float))
