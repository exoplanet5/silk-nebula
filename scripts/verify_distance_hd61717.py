"""verify_distance (a): independent vetting of the narrow Ca II component at V_helio ~ -17 km/s toward HD 61717.

Inputs (on disk, untrusted FITS): ESO phase-3 UVES spectra of HD 61717 (prog. 112.25VN.001, 2023-11-12),
the other UVES stars of the same programme, 11 CMi (UVES 2000) and three HARPS stars used as telluric references.
Gaia DR3 / APOGEE DR17 metadata fetched to data/raw/verify_distance_meta/ (see notes/verify_distance.md).

Steps
 1. Frame check: telluric lines of HD 61717 vs 11 CMi in native (topocentric) wavelengths (BERVs differ by 56 km/s).
 2. Instrumental resolution from unresolved telluric O2 lines in the red arm.
 3. Photospheric lines: one- vs two-component rotational fits (Mg II 4481, Fe II 5018, Fe II 5316, Si II 6347,
    Ti II 4468) -> SB2 test, component velocities, v sin i.
 4. Ca II K+H: simultaneous Gaussian-tau components convolved with the LSF on native pixels, empirical noise,
    block-bootstrap + LSF + continuum-order Monte Carlo; per-line AOD (saturation test); LSF scan.
 5. Na I D2/D1 at -17: telluric template = median of three HARPS spectra moved to the topocentric frame and
    degraded to R(HD 61717); joint D2+D1 fit with the Ca II -17 velocity (b free) + local components.
 6. Matched-filter limits at -17 (and +12) for Ti II 3384, Ca I 4227, CH+ 4232, CH 4300, CN 3874, Fe I 3860;
    shell-star narrow cores (Ti II 3759/3761, Fe II 4233, Fe II 5169) at -17 km/s.
 7. Uniqueness: matched filter at the same heliocentric and the same topocentric velocity in the programme stars.
Outputs: data/derived/verify_distance_hd61717.json, data/derived/verify_distance_hd61717_profiles.npz
Run: ~/.venvs/astro313/bin/python -I scripts/verify_distance_hd61717.py
"""
import json
import sys
import glob
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_distance_uves import read_uves, SPEC, ROOT, C, LINES  # noqa: E402
from astropy.coordinates import SkyCoord, EarthLocation  # noqa: E402
from astropy.time import Time  # noqa: E402
from astropy.io import fits  # noqa: E402
import astropy.units as u  # noqa: E402

OUT = {}
rng = np.random.default_rng(61717)
PARANAL = EarthLocation.from_geodetic(-70.4045 * u.deg, -24.6268 * u.deg, 2635 * u.m)

blue = read_uves(SPEC / 'eso_HD61717_ADP.2024-03-09T06:44:09.485.fits')
red = read_uves(SPEC / 'eso_HD61717_ADP.2024-03-09T06:44:09.481.fits')
BERV = float(blue['berv'])
ref_red = read_uves(SPEC / 'eso_11CMi_ADP.2020-09-14T18:52:16.793.fits')
H = ref_red['header']
berv_ref = SkyCoord(H['RA'] * u.deg, H['DEC'] * u.deg).radial_velocity_correction(
    'barycentric', obstime=Time(H['MJD-OBS'], format='mjd'), location=PARANAL).to('km/s').value
OUT['frames'] = dict(hd61717_specsys=blue['specsys'], hd61717_berv=BERV, hd61717_mjd=blue['mjd'],
                     hd61717_date=blue['header'].get('DATE-OBS'), hd61717_prog=blue['header'].get('PROG_ID'),
                     hd61717_airmass=red['header'].get('HIERARCH ESO TEL AIRM START'),
                     ref_11cmi_specsys=ref_red['specsys'], ref_11cmi_berv_astropy=round(berv_ref, 3),
                     ref_11cmi_airmass=H.get('HIERARCH ESO TEL AIRM START'))


def resample(w, f, grid):
    return np.interp(grid, w, f, left=np.nan, right=np.nan)


def contnorm(w, f, deg=3, clip=(1.5, 3.0), it=6):
    m = np.isfinite(f)
    x = (w - w.mean()) / (np.ptp(w) / 2)
    for _ in range(it):
        p = np.polyfit(x[m], f[m], deg)
        r = f - np.polyval(p, x)
        s = np.nanstd(r[m])
        m = m & (r > -clip[0] * s) & (r < clip[1] * s)
    return f / np.polyval(p, x)


# ---------------------------------------------------------------- 1. frame check (telluric cross-correlation)
segments = [(5885.0, 5888.6), (5890.3, 5895.0), (5896.5, 5900.0), (6276.0, 6300.0), (6300.5, 6320.0)]
lags = np.arange(-15, 15.01, 0.1)
cc_tot = np.zeros_like(lags)
for a, b in segments:
    g = np.arange(a, b, 0.005)
    d1 = 1 - contnorm(g, resample(red['wave'], red['flux_red'], g))
    f2 = contnorm(g, resample(ref_red['wave'], ref_red['flux_red'], g))
    for i, L in enumerate(lags):
        f2s = 1 - np.interp(g, g * (1 + L / C), f2)
        ok = np.isfinite(d1) & np.isfinite(f2s)
        cc_tot[i] += np.sum((d1[ok] - d1[ok].mean()) * (f2s[ok] - f2s[ok].mean()))
ipk = int(np.argmax(cc_tot))
pp = np.polyfit(lags[ipk - 5:ipk + 6], cc_tot[ipk - 5:ipk + 6], 2)
lag_pk = -pp[1] / (2 * pp[0])
OUT['frames']['telluric_xcorr_lag_hd61717_minus_11cmi_kms'] = round(float(lag_pk), 2)
OUT['frames']['expected_lag_if_products_were_barycentric_kms'] = round(BERV - berv_ref, 2)

# ---------------------------------------------------------------- 2. LSF from telluric lines (red arm)
g = np.arange(6270, 6330, 0.004)
fr = contnorm(g, resample(red['wave'], red['flux_red'], g))
fr2 = contnorm(g, resample(ref_red['wave'], ref_red['flux_red'], g))
pk, _ = find_peaks(1 - fr, height=0.04, distance=40)
fwhm_kms = []
for k in pk:
    lam = g[k]
    if np.nanmin(1 - fr2[(g > lam - 0.03) & (g < lam + 0.03)]) < 0.02:   # must be telluric (in 11 CMi too)
        continue
    m = (g > lam - 0.12) & (g < lam + 0.12)
    x, y = g[m], fr[m]

    def mod(p):
        return p[0] + p[4] * (x - lam) - p[1] * np.exp(-0.5 * ((x - p[2]) / p[3]) ** 2)
    r = least_squares(lambda p: mod(p) - y, [np.median(y), 0.05, lam, 0.03, 0.0],
                      bounds=([0, 0, lam - 0.05, 0.005, -1], [2, 1, lam + 0.05, 0.2, 1]))
    if 0.02 < r.x[1] < 0.6:
        fwhm_kms.append(2.3548 * r.x[3] / r.x[2] * C)
fwhm_kms = np.array(fwhm_kms)
fwhm_red_inst = np.sqrt(max(np.median(fwhm_kms) ** 2 - 1.2 ** 2, 1.0))   # remove ~1.2 km/s intrinsic O2 width
FWHM_BLUE = fwhm_red_inst * red['R'] / blue['R']
OUT['lsf'] = dict(n_telluric_lines=int(len(fwhm_kms)), fwhm_kms_median=round(float(np.median(fwhm_kms)), 2),
                  fwhm_kms_p16_84=[round(float(x), 2) for x in np.percentile(fwhm_kms, [16, 84])],
                  R_red_header=red['R'], R_blue_header=blue['R'], fwhm_blue_header_kms=round(C / blue['R'], 2),
                  fwhm_blue_adopted_kms=round(float(FWHM_BLUE), 2),
                  note='measured telluric FWHM minus ~1.2 km/s intrinsic, scaled to the blue arm by the header R ratio')

# ---------------------------------------------------------------- 3. photospheric lines: SB2 test


def rot_kernel(vsini, dv, eps=0.6):
    kv = np.arange(-vsini, vsini + dv / 2, dv)
    x = kv / vsini
    k = np.clip(2 * (1 - eps) * np.sqrt(np.clip(1 - x ** 2, 0, None)) + 0.5 * np.pi * eps * (1 - x ** 2), 0, None)
    return k / k.sum()


def rot_prof(v, rv, vsini, sub=None):
    """unit-area rotational profile (with 4 km/s intrinsic Gaussian) on grid v; sub=[(offset,weight)]"""
    dv = v[1] - v[0]
    base = np.zeros_like(v)
    for off, wgt in (sub or [(0.0, 1.0)]):
        base += wgt * np.exp(-0.5 * ((v - rv - off) / 4.0) ** 2)
    if vsini > 2 * dv:
        base = np.convolve(base, rot_kernel(vsini, dv), mode='same')
    return base / (base.sum() * dv)


PH_LINES = [('MgII4481', 4481.126, [(0.0, 10 ** 0.749), ((4481.150 - 4481.126) / 4481.126 * C, 10 ** -0.553),
                                     ((4481.325 - 4481.126) / 4481.126 * C, 10 ** 0.594)]),
            ('FeII5018', 5018.440, None), ('FeII5316', 5316.615, None), ('SiII6347', 6347.109, None),
            ('TiII4468', 4468.507, None)]
VG = np.arange(-450, 450.01, 1.0)
ph_data = []
for nm, lam0, sub in PH_LINES:
    sp = blue if lam0 < 4600 else red
    v = (sp['wave'] / lam0 - 1) * C + sp['berv']
    m = np.isfinite(sp['flux_red'])
    f = np.interp(VG, v[m], sp['flux_red'][m])
    edge = np.abs(VG) > 320
    p = np.polyfit(VG[edge], f[edge], 1)
    f = f / np.polyval(p, VG)
    ph_data.append((nm, sub, f))
FIT_V = np.abs(VG) < 330


def ph_model(p, ncomp):
    out = []
    for j, (nm, sub, f) in enumerate(ph_data):
        y = np.ones_like(VG)
        for c in range(ncomp):
            rv, vs = p[2 * c], p[2 * c + 1]
            eq = p[2 * ncomp + j * ncomp + c]
            y = y - eq * rot_prof(VG, rv, vs, sub)
        out.append(y)
    return out


def ph_fit(ncomp, starts):
    best = None
    for s in starts:
        nl = len(ph_data)
        p0 = list(s) + [3.0] * (nl * ncomp)
        lo = [-300, 15] * ncomp + [0] * (nl * ncomp)
        hi = [300, 400] * ncomp + [100] * (nl * ncomp)

        def res(p):
            mods = ph_model(p, ncomp)
            return np.concatenate([(mo - f)[FIT_V] for mo, (_, _, f) in zip(mods, ph_data)])
        r = least_squares(res, p0, bounds=(lo, hi))
        if best is None or r.cost < best.cost:
            best = r
    n = sum(FIT_V.sum() for _ in ph_data)
    return best, 2 * best.cost, n


r1, chi1, n1 = ph_fit(1, [(30, 150), (40, 230), (10, 300)])
r2, chi2_, n2 = ph_fit(2, [(-110, 80, 130, 80), (-90, 60, 140, 60), (-130, 100, 120, 100)])
rms1 = np.sqrt(chi1 / n1)
rms2 = np.sqrt(chi2_ / n2)
bic1 = n1 * np.log(chi1 / n1) + len(r1.x) * np.log(n1)
bic2 = n2 * np.log(chi2_ / n2) + len(r2.x) * np.log(n2)
ncl = len(ph_data)
eqw2 = np.array(r2.x[4:]).reshape(ncl, 2)
OUT['stellar'] = dict(
    one_component=dict(rv=round(float(r1.x[0]), 1), vsini=round(float(r1.x[1]), 1), rms=round(float(rms1), 5),
                       BIC=round(float(bic1), 1)),
    two_component=dict(rv1=round(float(r2.x[0]), 1), vsini1=round(float(r2.x[1]), 1), rv2=round(float(r2.x[2]), 1),
                       vsini2=round(float(r2.x[3]), 1), rms=round(float(rms2), 5), BIC=round(float(bic2), 1),
                       eqw_ratio_comp2_over_comp1_per_line={nm: round(float(e[1] / e[0]), 2) if e[0] > 0 else None
                                                            for (nm, _, _), e in zip(ph_data, eqw2)}),
    delta_BIC_1minus2=round(float(bic1 - bic2), 1),
    lines=[nm for nm, _, _ in ph_data],
    note='equal-velocity-for-all-lines rotational profiles (limb darkening 0.6); a large positive delta_BIC favours SB2')
OUT['stellar']['catalogue'] = dict(
    simbad=dict(sp='A2', V=7.57, rv=41.824, e_rv=0.106, vsini=93.3, vsini_ref='2020AJ....160..120J (APOGEE DR16)'),
    gaia_dr3=dict(source_id=866945131588948864, plx=3.6711, e_plx=0.0670, ruwe=2.143, rv=35.71, e_rv=1.32,
                  vbroad=160.4, teff_gspphot=10050.8, logg=3.9965, dist_gspphot=230.5, nss_flag=0,
                  file='data/raw/verify_distance_meta/gaia_dr3_hd61717.tsv'),
    apogee_dr17='data/raw/verify_distance_meta/apogee_dr17_hd61717_{star,visits}.csv')
# APOGEE visits
try:
    rows = [ln.split(',') for ln in (ROOT / 'data/raw/verify_distance_meta/apogee_dr17_hd61717_visits.csv').read_text().splitlines()[2:]]
    vh = np.array([float(r[3]) for r in rows])
    mj = np.array([float(r[1]) for r in rows])
    OUT['stellar']['apogee_visits'] = dict(n=int(len(vh)), mjd_range=[float(mj.min()), float(mj.max())],
                                           vhelio_min=float(vh.min()), vhelio_max=float(vh.max()),
                                           vhelio_mean=round(float(vh.mean()), 2), vhelio_std=round(float(vh.std(ddof=1)), 2),
                                           median_vrelerr=0.32, starflag='RV_SUSPECT + SUSPECT_BROAD_LINES (visit bits 17, 20)')
except Exception as ex:  # noqa: BLE001
    OUT['stellar']['apogee_visits'] = str(ex)

# ---------------------------------------------------------------- 4. Ca II K + H optical-depth fit on native pixels


def tau_gauss(v, N, b, v0, lam, f):
    return 1.4974e-15 * N * f * lam / b * np.exp(-((v - v0) / b) ** 2)


FINE = np.arange(-100, 100.001, 0.1)


def native_window(sp, lam0, vlo=-85, vhi=85):
    v = (sp['wave'] / lam0 - 1) * C + sp['berv']
    m = (v > vlo) & (v < vhi) & np.isfinite(sp['flux_red'])
    return v[m], sp['flux_red'][m]


vK, fK = native_window(blue, LINES['CaK'][0])
vH, fH = native_window(blue, LINES['CaH'][0])


def emp_noise(v, f):
    m = (np.abs(v) > 40) & (np.abs(v) < 85)
    x, y = v[m], f[m]
    sel = x > 0
    r = np.concatenate([y[sel] - np.polyval(np.polyfit(x[sel], y[sel], 2), x[sel]),
                        y[~sel] - np.polyval(np.polyfit(x[~sel], y[~sel], 2), x[~sel])])
    return np.std(r) / np.median(y)


sK = np.median(fK[np.abs(vK) > 50])
sH = np.median(fH[np.abs(vH) > 50])
fK, fH = fK / sK, fH / sH
nK, nH = emp_noise(vK, fK), emp_noise(vH, fH)
pipe_err_ratio = float(np.median(blue['err_red'] / blue['flux_red']))
NCOMP = 3


def ca_model(p, fwhm, vK=vK, vH=vH, ncont=3, only=None):
    tK = np.zeros_like(FINE)
    tH = np.zeros_like(FINE)
    for i in range(NCOMP):
        if only is not None and i not in only:
            continue
        lN, b, v0 = p[3 * i:3 * i + 3]
        tK += tau_gauss(FINE, 10 ** lN, b, v0, *LINES['CaK'])
        tH += tau_gauss(FINE, 10 ** lN, b, v0, *LINES['CaH'])
    sig = fwhm / 2.3548 / 0.1
    eK = np.interp(vK, FINE, gaussian_filter1d(np.exp(-tK), sig, mode='nearest'))
    eH = np.interp(vH, FINE, gaussian_filter1d(np.exp(-tH), sig, mode='nearest'))
    cK = np.polyval(p[3 * NCOMP:3 * NCOMP + ncont], vK / 100)
    cH = np.polyval(p[3 * NCOMP + ncont:3 * NCOMP + 2 * ncont], vH / 100)
    return cK * eK, cH * eH, cK, cH


def fit_ca(fwhm, fK=fK, fH=fH, ncont=3, x0=None):
    p0 = [11.5, 1.5, -17.0, 11.2, 4.0, 9.0, 11.3, 3.0, 14.5] + [0] * (ncont - 1) + [1.0] + [0] * (ncont - 1) + [1.0]
    lo = [9.5, 0.2, -22, 9.5, 0.5, 0, 9.5, 0.5, 11] + [-5] * (2 * ncont)
    hi = [13.5, 8, -12, 13, 12, 13, 13, 10, 22] + [5] * (2 * ncont)
    if x0 is not None and len(x0) == len(p0):
        p0 = list(np.clip(x0, np.array(lo) + 1e-6, np.array(hi) - 1e-6))

    def res(p):
        mK, mH, _, _ = ca_model(p, fwhm, ncont=ncont)
        return np.concatenate([(mK - fK) / nK, (mH - fH) / nH])
    return least_squares(res, p0, bounds=(lo, hi), x_scale='jac')


r0 = fit_ca(FWHM_BLUE)
p_best = r0.x
chi2 = 2 * r0.cost
dof = len(vK) + len(vH) - len(p_best)
mK0, mH0, cK0, cH0 = ca_model(p_best, FWHM_BLUE)
mc = []
for it in range(300):
    fw = FWHM_BLUE * rng.uniform(0.85, 1.15)
    nc = int(rng.integers(2, 5))
    out = []
    for f_, m_ in ((fK, mK0), (fH, mH0)):
        rr = f_ - m_
        L = 6
        nb = len(rr) // L + 1
        idx = rng.integers(0, len(rr) // L, nb)
        rb = np.concatenate([rr[i * L:(i + 1) * L] for i in idx])[:len(rr)]
        out.append(m_ + rb)
    r = fit_ca(fw, fK=out[0], fH=out[1], ncont=nc)
    mc.append(r.x[:9])
mc = np.array(mc)
ca_res = []
for i in range(NCOMP):
    d = {}
    for j, nm in enumerate(['logN', 'b', 'v']):
        d[nm] = round(float(p_best[3 * i + j]), 3)
        d['e_' + nm] = round(float(np.std(mc[:, 3 * i + j])), 3)
    ca_res.append(d)
b17, eb17 = p_best[1], np.std(mc[:, 1])
OUT['caII_fit'] = dict(components=ca_res, chi2=round(float(chi2), 1), dof=int(dof),
                       noise_per_native_px_K=round(float(nK), 5), noise_per_native_px_H=round(float(nH), 5),
                       pipeline_err_over_flux=round(pipe_err_ratio, 4),
                       fwhm_lsf_kms=round(float(FWHM_BLUE), 2), n_mc=len(mc),
                       T_max_K_thermal=round(float(40 * 1.6605e-24 * (b17 * 1e5) ** 2 / (2 * 1.380649e-16))),
                       T_max_K_thermal_b_plus_2sig=round(float(40 * 1.6605e-24 * ((b17 + 2 * eb17) * 1e5) ** 2 / (2 * 1.380649e-16))),
                       method='K+H simultaneous Gaussian-tau, LSF-convolved on 0.1 km/s grid, native pixels, empirical noise; '
                              'MC = 6-px block bootstrap of residuals x LSF +-15% x continuum order 2-4')
bscan = []
for Rb in (55000, 65030, 70000, 75000, 85000):
    r = fit_ca(C / Rb, x0=p_best)
    bscan.append(dict(R=Rb, logN=round(float(r.x[0]), 3), b=round(float(r.x[1]), 3), v=round(float(r.x[2]), 3),
                      chi2=round(float(2 * r.cost), 1)))
OUT['caII_fit']['lsf_scan'] = bscan


def aod(v, f, cont, lam, fval, lo, hi, noise):
    m = (v > lo) & (v < hi)
    x = np.clip(f[m] / cont[m], 0.02, None)
    dv = np.gradient(v)[m]
    tau = -np.log(x)
    N = 3.768e14 / (fval * lam) * np.sum(tau * dv)
    W = np.sum((1 - x) * dv) / C * lam * 1000
    eW = noise * np.sqrt(m.sum()) * np.median(dv) / C * lam * 1000 * 1.3   # 1.3: residual pixel correlation
    return N, W, eW


NKa, WK, eWK = aod(vK, fK, cK0, *LINES['CaK'], -24, -10, nK)
NHa, WH, eWH = aod(vH, fH, cH0, *LINES['CaH'], -24, -10, nH)
OUT['caII_aod_m24_m10'] = dict(W_K_mA=round(WK, 2), eW_K_mA=round(eWK, 2), W_H_mA=round(WH, 2), eW_H_mA=round(eWH, 2),
                               logNa_K=round(float(np.log10(NKa)), 3), logNa_H=round(float(np.log10(NHa)), 3),
                               W_K_over_W_H=round(WK / WH, 3))
# ambient components AOD, -5..+30
NKb, WKb, _ = aod(vK, fK, cK0, *LINES['CaK'], -6, 30, nK)
NHb, WHb, _ = aod(vH, fH, cH0, *LINES['CaH'], -6, 30, nH)
OUT['caII_aod_m6_p30'] = dict(W_K_mA=round(WKb, 2), W_H_mA=round(WHb, 2), logNa_K=round(float(np.log10(NKb)), 3),
                              logNa_H=round(float(np.log10(NHb)), 3))

# ---------------------------------------------------------------- 5. Na I D: HARPS telluric template


def harps_topo(fn):
    with fits.open(SPEC / fn, memmap=False) as h:
        Hh = h[0].header
        d = h[1].data
        w = np.array(d['WAVE'][0], float)
        f = np.array(d['FLUX'][0], float)
    bv = float(Hh['HIERARCH ESO DRS BERV'])
    # barycentric -> topocentric:  lambda_topo = lambda_bary / (1 + BERV/c)
    return w / (1 + bv / C), f, bv, float(Hh.get('HIERARCH ESO TEL AIRM START', np.nan))


TELREF = ['eso_8Cnc_ADP.2014-09-16T11:04:03.200.fits', 'eso_HD52479_ADP.2014-09-26T16:53:53.353.fits',
          'eso_HD67959_ADP.2014-09-26T16:53:39.010.fits']
GT = np.arange(5883.0, 5902.0, 0.004)
R_TARGET = red['R'] if red['R'] else 74450.0
sig_conv_kms = np.sqrt((C / R_TARGET) ** 2 - (C / 115000.0) ** 2) / 2.3548
tel_stack = []
tel_meta = []
for fn in TELREF:
    w, f, bv, am = harps_topo(fn)
    fi = contnorm(GT, resample(w, f, GT), deg=3)
    # mask each reference star's own Na I (helio -20..+45 -> topo = helio - BERV)
    for lam0 in (LINES['NaD2'][0], LINES['NaD1'][0]):
        vt = (GT / lam0 - 1) * C
        mm = (vt + bv > -20) & (vt + bv < 45)
        fi[mm] = np.nan
    sig_px = sig_conv_kms / C * 5892 / 0.004
    ok = np.isfinite(fi)
    fs = gaussian_filter1d(np.where(ok, fi, 1.0), sig_px)
    fs[~ok] = np.nan
    tau = -np.log(np.clip(fs, 0.05, None))
    tel_stack.append(tau / np.nanmax(np.where(np.isfinite(tau), tau, 0)[(GT > 5891.5) & (GT < 5891.8)]))  # norm on 5891.65 line
    tel_meta.append(dict(file=fn, berv=round(bv, 2), airmass=am))
tel_stack = np.array(tel_stack)
TAU_T = np.nanmedian(tel_stack, axis=0)
TAU_T = np.where(np.isfinite(TAU_T), TAU_T, 0.0)
f61 = contnorm(GT, resample(red['wave'], red['flux_red'], GT), deg=3)
v_hel = {k: (GT / LINES[k][0] - 1) * C + BERV for k in ('NaD2', 'NaD1')}
noise_na = np.nanstd(np.diff(f61[(GT > 5899) & (GT < 5901.5)])) / np.sqrt(2)
# native-pixel noise
wn = red['wave']
mn = (wn > 5899.0) & (wn < 5901.5)
xn, yn = wn[mn], red['flux_red'][mn]
noise_na_native = float(np.std(yn - np.polyval(np.polyfit(xn, yn, 2), xn)) / np.median(yn))
px_native = float(np.median(np.diff(wn[mn])))
# fit window: both lines -60..+50 helio; model = cont * exp(-alpha*TAU_T - tau_ISM)
fitm = ((v_hel['NaD2'] > -60) & (v_hel['NaD2'] < 50)) | ((v_hel['NaD1'] > -60) & (v_hel['NaD1'] < 50))
fitm &= np.isfinite(f61)


def na_model(p, with_m17=True):
    alpha, lN17, b17_, lNa, ba, va, lNb, bb, vb = p[:9]
    tau = alpha * TAU_T
    for k in ('NaD2', 'NaD1'):
        lam, fv = LINES[k]
        vv = v_hel[k]
        if with_m17:
            tau = tau + tau_gauss(vv, 10 ** lN17, b17_, p_best[2], lam, fv)
        tau = tau + tau_gauss(vv, 10 ** lNa, ba, va, lam, fv) + tau_gauss(vv, 10 ** lNb, bb, vb, lam, fv)
    sig_px = (C / R_TARGET) / 2.3548 / C * 5892 / 0.004
    return gaussian_filter1d(np.exp(-tau), sig_px) * (p[9] + p[10] * (GT - 5892))


p0 = [1.0, 10.2, 1.7, 10.3, 3.0, 5.0, 11.0, 2.5, 13.0, 1.0, 0.0]
lo = [0.0, 8.0, 1.69, 8.0, 0.5, -5, 8.0, 0.5, 9, 0.9, -0.1]
hi = [5.0, 12.5, 1.71, 12.5, 8.0, 10, 12.5, 8.0, 20, 1.1, 0.1]   # b(-17) fixed at 1.7 km/s (Ca II b x 1-1.3)
rna = least_squares(lambda p: (na_model(p) - f61)[fitm] / noise_na, p0, bounds=(lo, hi))
# Monte Carlo for N(Na I, -17)
mcna = []
base = na_model(rna.x)
resid = (f61 - base)
for it in range(200):
    L = 40
    idx = rng.integers(0, len(GT) // L, len(GT) // L + 1)
    rb = np.concatenate([resid[i * L:(i + 1) * L] for i in idx])[:len(GT)]
    # telluric template perturbation: use one reference at random 1/3 of the time
    r = least_squares(lambda p: (na_model(p) - (base + rb))[fitm] / noise_na, rna.x, bounds=(lo, hi))
    mcna.append(r.x[:3])
mcna = np.array(mcna)
# per-line EW at -17 after dividing the fitted telluric + local components
tel_only = gaussian_filter1d(np.exp(-rna.x[0] * TAU_T), (C / R_TARGET) / 2.3548 / C * 5892 / 0.004)
pn = rna.x.copy()
pn[1] = 5.0
loc_only = na_model(pn, with_m17=False) / (pn[9] + pn[10] * (GT - 5892))
resid_17 = f61 / (rna.x[9] + rna.x[10] * (GT - 5892)) / loc_only
ew17 = {}
for k in ('NaD2', 'NaD1'):
    lam, fv = LINES[k]
    m = (v_hel[k] > -22) & (v_hel[k] < -12)
    W = float(np.sum(1 - resid_17[m]) * 0.004 * 1000)
    # empirical: scatter of 10-km/s box sums of the residual spectrum at line-free velocities
    vk = v_hel[k]
    boxes = []
    for vc in np.arange(-75, -30, 10.0).tolist() + np.arange(30, 55, 10.0).tolist():
        mb = (vk > vc - 5) & (vk < vc + 5) & np.isfinite(resid_17)
        if mb.sum() > 20:
            boxes.append(np.sum(1 - resid_17[mb]) * 0.004 * 1000)
    eW = float(np.std(boxes, ddof=1)) if len(boxes) > 3 else np.nan
    Wtel = float(np.sum(1 - tel_only[m]) * 0.004 * 1000)
    ew17[k] = dict(W_mA=round(W, 2), eW_mA=round(eW, 2), W_telluric_removed_mA=round(Wtel, 2))
OUT['naI'] = dict(telluric_refs=tel_meta, telluric_scale_alpha=round(float(rna.x[0]), 3),
                  m17=dict(logN=round(float(rna.x[1]), 3), e_logN=round(float(np.std(mcna[:, 1])), 3),
                           b=round(float(rna.x[2]), 2), e_b=round(float(np.std(mcna[:, 2])), 2),
                           v_fixed=round(float(p_best[2]), 2)),
                  local=[dict(logN=round(float(rna.x[3]), 3), b=round(float(rna.x[4]), 2), v=round(float(rna.x[5]), 2)),
                         dict(logN=round(float(rna.x[6]), 3), b=round(float(rna.x[7]), 2), v=round(float(rna.x[8]), 2))],
                  per_line_ew_m17=ew17, noise_per_native_px_5899_5901=round(noise_na_native, 5),
                  logN_from_ew=dict({k: round(float(np.log10(1.13e17 * max(ew17[k]['W_mA'], 1e-3) / (LINES[k][1] * LINES[k][0] ** 2))), 3)
                                     for k in ('NaD2', 'NaD1')}),
                  chi2_red=round(float(2 * rna.cost / (fitm.sum() - len(rna.x))), 2),
                  note='joint D2+D1 fit; the -17 component velocity is fixed to the Ca II value; HARPS telluric template '
                       'median of 3 stars (each star\'s own Na I masked), degraded to R=74,450')
NCa = 10 ** p_best[0]
NNa = 10 ** rna.x[1]
eNa = np.std(10 ** mcna[:, 1])
OUT['ratio_m17'] = dict(logN_CaII=round(float(p_best[0]), 3), logN_NaI=round(float(rna.x[1]), 3),
                        CaII_over_NaI=round(float(NCa / NNa), 1),
                        CaII_over_NaI_68=[round(float(NCa / (NNa + eNa)), 1), round(float(NCa / max(NNa - eNa, 1e8)), 1)])
# local ratio for comparison
NCa_loc = 10 ** p_best[3] + 10 ** p_best[6]
NNa_loc = 10 ** rna.x[3] + 10 ** rna.x[6]
OUT['ratio_local'] = dict(logN_CaII=round(float(np.log10(NCa_loc)), 3), logN_NaI=round(float(np.log10(NNa_loc)), 3),
                          CaII_over_NaI=round(float(NCa_loc / NNa_loc), 2))

# ---------------------------------------------------------------- 6. matched-filter limits for other species


def matched(sp, lam0, vc, b, fwhm, vlo=-90, vhi=90, mask=((-25, 25),), highpass=False):
    """Return EW (mA) and 1-sigma error of a weak line with fixed shape (Gaussian tau, b, LSF) at vc.
    highpass=True divides by a 15-km/s running median first, so only narrow (FWHM < ~8 km/s) features survive."""
    v = (sp['wave'] / lam0 - 1) * C + sp['berv']
    m = (v > vlo) & (v < vhi) & np.isfinite(sp['flux_red'])
    v, f = v[m], sp['flux_red'][m]
    cm = np.ones_like(v, bool)
    for a, bb in mask:
        cm &= ~((v > a) & (v < bb))
    cm &= np.abs(v - vc) > 3 * fwhm
    if highpass:
        from scipy.ndimage import median_filter
        npx = int(round(15.0 / np.median(np.diff(v)))) | 1
        fn = f / median_filter(f, size=npx, mode='nearest')
        cm = np.abs(v - vc) > 3 * fwhm
        for a, bb in mask:
            cm &= ~((v > a) & (v < bb))
    else:
        p = np.polyfit(v[cm] / 100, f[cm], 3)
        fn = f / np.polyval(p, v / 100)
    noise = np.std(fn[cm] - 1)
    sig = np.sqrt(b ** 2 / 2 + (fwhm / 2.3548) ** 2)
    shape = np.exp(-0.5 * ((v - vc) / sig) ** 2)                # depth profile per unit peak depth
    dv = np.gradient(v)
    area_per_depth = np.sum(shape * dv)
    sel = np.abs(v - vc) < 4 * sig
    A = np.sum(((1 - fn) * shape)[sel]) / np.sum((shape ** 2)[sel])
    eA = noise / np.sqrt(np.sum((shape ** 2)[sel])) * 1.3
    W = A * area_per_depth / C * lam0 * 1000
    eW = eA * area_per_depth / C * lam0 * 1000
    return float(W), float(eW), float(noise)


others = {}
for key in ('TiII3384', 'CaI4227', 'CHp4232', 'CH4300', 'CN3874', 'FeI3860', 'CaK', 'CaH'):
    lam0, fv = LINES[key]
    W, eW, nz = matched(blue, lam0, p_best[2], max(b17, 1.0), FWHM_BLUE)
    Wl, eWl, _ = matched(blue, lam0, 12.0, 4.0, FWHM_BLUE, mask=((-25, -8),))
    Nup = 1.13e17 * (max(W, 0) + 3 * eW) / (fv * lam0 ** 2)
    Nl = 1.13e17 * Wl / (fv * lam0 ** 2)
    others[key] = dict(W_m17_mA=round(W, 2), eW_m17_mA=round(eW, 2), logN_m17_3sig_upper=round(float(np.log10(Nup)), 2),
                       W_local_p12_mA=round(Wl, 2), eW_local_mA=round(eWl, 2),
                       logN_local=round(float(np.log10(Nl)), 2) if Wl > 3 * eWl else None,
                       noise_per_px=round(nz, 4))
OUT['other_lines'] = others
tiup = others['TiII3384']['logN_m17_3sig_upper']
OUT['TiII_over_CaII'] = dict(m17_3sig_upper=round(float(10 ** (tiup - p_best[0])), 2),
                             local=round(float(10 ** (others['TiII3384']['logN_local'] - np.log10(NCa_loc))), 2)
                             if others['TiII3384']['logN_local'] else None)
shell = {}
for nm, lam0 in (('TiII3759', 3759.291), ('TiII3761', 3761.320), ('FeII4233', 4233.172), ('FeII5169', 5169.033)):
    sp = blue if lam0 < 4600 else red
    fw = FWHM_BLUE if sp is blue else fwhm_red_inst
    W, eW, _ = matched(sp, lam0, p_best[2], max(b17, 1.0), fw, mask=())
    shell[nm] = dict(W_m17_mA=round(W, 2), eW_mA=round(eW, 2))
OUT['shell_test_m17'] = shell

# ---------------------------------------------------------------- 7. uniqueness across programme stars
uniq = []
for fn in sorted(glob.glob(str(SPEC / 'eso_HD*_ADP.2024-*.fits'))):
    sp = read_uves(fn)
    if not (sp['wave'].min() < 3933 < sp['wave'].max()):
        continue
    star = sp['header']['OBJECT']
    vt = p_best[2] - BERV + sp['berv']
    row = dict(star=star, berv=round(float(sp['berv']), 2), mjd=round(float(sp['mjd']), 3))
    for ln in ('CaK', 'CaH'):
        lam0 = LINES[ln][0]
        Wh, eWh, _ = matched(sp, lam0, p_best[2], max(b17, 1.0), FWHM_BLUE, mask=((-5, 25),))
        Wt, eWt, _ = matched(sp, lam0, vt, max(b17, 1.0), FWHM_BLUE, mask=((-5, 25),))
        Wn, eWn, _ = matched(sp, lam0, p_best[2], max(b17, 1.0), FWHM_BLUE, mask=((-5, 25),), highpass=True)
        row[ln] = dict(W_same_vhelio=round(Wh, 2), eW=round(eWh, 2), W_same_vtopo=round(Wt, 2), eW_topo=round(eWt, 2),
                       vtopo_equiv_helio=round(float(vt), 2), W_narrow_highpass=round(Wn, 2), eW_narrow=round(eWn, 2))
    uniq.append(row)
OUT['uniqueness_programme_stars'] = uniq

np.savez_compressed(ROOT / 'data/derived/verify_distance_hd61717_profiles.npz', vK=vK, fK=fK, vH=vH, fH=fH,
                    mK=mK0, mH=mH0, cK=cK0, cH=cH0, GT=GT, f61=f61, tel=tel_only, na_model=na_model(rna.x),
                    vD2=v_hel['NaD2'], vD1=v_hel['NaD1'], resid17=resid_17, VG=VG,
                    ph_f=np.array([f for _, _, f in ph_data]), ph_m1=np.array(ph_model(r1.x, 1)),
                    ph_m2=np.array(ph_model(r2.x, 2)), ph_names=np.array([n for n, _, _ in ph_data]))
(ROOT / 'data/derived/verify_distance_hd61717.json').write_text(json.dumps(OUT, indent=1, default=float))
print(json.dumps({k: OUT[k] for k in OUT if k not in ('uniqueness_programme_stars',)}, indent=1, default=float))
print(json.dumps(OUT['uniqueness_programme_stars'], default=float))
