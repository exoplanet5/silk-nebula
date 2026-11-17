"""verify_hiwham: WHAM-SS DR1 H-alpha on vs off the Silk spine at the revised [O III] velocity, with nulls.

Independent re-implementation (raw WHAM DR1 table + silk_trace.json only).

1. Beams in RA 90-142, Dec -2..+42. Each spectrum is put on a common 2 km/s grid (interpolation of the native
   ~2.05 km/s channels) and characterised by
     * I_tot      : INTEN column (-80..+80 km/s, survey value)
     * I_win      : integral over the filament window  v_c +- 25 km/s, v_c = V_helio(+13.5) + (V_LSR - V_helio)
                    at the beam (V_LSR ~ +0.3 south .. +4 north)
     * m1, sig    : first moment and dispersion over -45..+45 km/s after a linear baseline from |v| 55-80
     * single-Gaussian + linear-baseline fit of the main component over -50..+50 (I_g, v_g, FWHM_g)
2. On/off: on = |d| < 0.5 deg, off = 1.5 < |d| < 3 deg (both sides), matched in along-spine coordinate
   (|ds| < 0.75 deg). Paired differences averaged over beams, for south (spine Dec <= 17.6), north, all.
3. Template regression: beam-averaged unit-peak Gaussian filament (FWHM 0.276 deg) on the E spine, D_b, computed
   on a 1-deg-diameter top-hat disc (437 sample points). Fit  y_b = quadratic(d, s) + I_peak * D_b  over beams
   with |d| < 3.5 deg. I_peak is the intrinsic (undiluted) filament H-alpha peak.
4. Nulls: the same statistics for rigid copies of the trace shifted in Galactic longitude by 6..45 deg (both
   signs, same latitude band, i.e. rotations about the Galactic pole) and in latitude by -4..+4 deg, with
   rotations of 0/+-15/180 deg and E-W mirrors; copies must avoid the real |d|<3.5 band.
5. Injection test: a synthetic filament (I_peak = 1 R, Gaussian in velocity at v_c, FWHM 40 km/s) is added to
   the real spectra with the beam dilution D_b and recovered through the same pipeline.
Outputs: data/derived/verify_hiwham_wham_beams.csv, verify_hiwham_wham_results.json, verify_hiwham_wham_spec.npz
"""
import json
import sys

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import (RAW, DER, strand, Trace, vlsr_minus_vhelio, tp_inv)

rng = np.random.default_rng(161116)
VG = np.arange(-90.0, 90.01, 2.0)
V_HELIO_OIII = 13.5          # E strand +13.0 (DESI), DESI ridge +14.5, SDSS +12.5: adopt +13.5 +- 1
HALF_WIN = 25.0
FWHM_IP = 12.0
SIGF = 0.276 / 2.3548

W = fits.open(RAW / 'wham_ss_dr1/wham-ss-DR1-v161116-170912.fits')[1].data
gl, gb = W['GAL-LON'].astype(float), W['GAL-LAT'].astype(float)
c = SkyCoord(gl * u.deg, gb * u.deg, frame='galactic').icrs
ra_all, dec_all = c.ra.deg, c.dec.deg
sel = np.where((gl > 125) & (gl < 285) & (gb > -4) & (gb < 50))[0]
ra, dec, l, b = ra_all[sel], dec_all[sel], gl[sel], gb[sel]
NB = len(sel)
spec = np.full((NB, len(VG)), np.nan)
var = np.full((NB, len(VG)), np.nan)
for j, i in enumerate(sel):
    vv = W['VELOCITY'][i].astype(float); ff = W['DATA'][i].astype(float); ee = W['VARIANCE'][i].astype(float)
    m = np.isfinite(vv) & np.isfinite(ff) & np.isfinite(ee)
    o = np.argsort(vv[m]); vv, ff, ee = vv[m][o], ff[m][o], ee[m][o]
    ok = (VG >= vv[0]) & (VG <= vv[-1])
    spec[j, ok] = np.interp(VG[ok], vv, ff)
    # native channel ~2.05 km/s ~ grid step: keep per-channel variance
    var[j, ok] = np.interp(VG[ok], vv, ee)
I_tot = W['INTEN'][sel].astype(float)
E_tot = W['ERROR'][sel].astype(float)
star = W['STAR'][sel].astype(bool)
corr = vlsr_minus_vhelio(ra, dec)
v_c = V_HELIO_OIII + corr
print('beams', NB, 'v_c range', v_c.min(), v_c.max())


def measure(sp, vr):
    """window integral, moments, single-gaussian main component for spectra array sp (n, nv)."""
    n = sp.shape[0]
    out = np.full((n, 10), np.nan)
    for j in range(n):
        f = sp[j]
        m = np.isfinite(f)
        if m.sum() < 50:
            continue
        w = m & (np.abs(VG - v_c[j]) <= HALF_WIN)
        out[j, 0] = np.sum(f[w]) * 2.0
        out[j, 1] = np.sqrt(np.nansum(vr[j][w])) * 2.0
        bl = m & (np.abs(VG) > 55) & (np.abs(VG) < 80)
        if bl.sum() > 6:
            pc = np.polyfit(VG[bl], f[bl], 1)
            g = f - np.polyval(pc, VG)
        else:
            g = f.copy()
        mm = m & (np.abs(VG) <= 45)
        m0 = np.sum(g[mm]) * 2.0
        if m0 > 0:
            m1 = np.sum(g[mm] * VG[mm]) * 2.0 / m0
            m2 = np.sqrt(max(np.sum(g[mm] * (VG[mm] - m1) ** 2) * 2.0 / m0, 0))
            out[j, 2], out[j, 3] = m1, m2
        fm = m & (np.abs(VG) <= 50)
        try:
            p, _ = curve_fit(lambda x, a, x0, s, c0, c1: a * np.exp(-0.5 * ((x - x0) / s) ** 2) + c0 + c1 * x,
                             VG[fm], f[fm], p0=[np.nanmax(f[fm]), 0.0, 13.0, 0.0, 0.0],
                             sigma=np.sqrt(np.clip(vr[j][fm], 1e-8, None)),
                             bounds=([0, -30, 4, -0.5, -0.01], [50, 30, 40, 0.5, 0.01]), maxfev=5000)
            out[j, 4] = p[0] * p[2] * np.sqrt(2 * np.pi)
            out[j, 5] = p[1]
            out[j, 6] = 2.3548 * p[2]
        except Exception:
            pass
        out[j, 7] = np.sum(f[m & (VG >= -80) & (VG <= 80)]) * 2.0
        out[j, 8] = np.sum(f[m & (VG >= -30) & (VG <= -2)]) * 2.0          # blue wing window (a posteriori)
        out[j, 9] = np.sum(f[m & (np.abs(VG - v_c[j]) <= 12)]) * 2.0       # narrow core window at v_c
    return out


MEAS = measure(spec, var)
COLS = ['I_win', 'eI_win', 'm1', 'sig', 'I_g', 'v_g', 'fwhm_g', 'I_sum', 'I_blue', 'I_core']
print('measured')

# ------------------------------------------------------------------ disc sampling for beam dilution
rr = np.sqrt(rng.uniform(0, 1, 437)) * 0.5
th = rng.uniform(0, 2 * np.pi, 437)
# use a deterministic sunflower pattern instead of random for stability
k = np.arange(437) + 0.5
rr = 0.5 * np.sqrt(k / 437); th = np.pi * (1 + 5 ** 0.5) * k
DX, DY = rr * np.cos(th), rr * np.sin(th)


def dilution(tr, bra, bdec):
    """Beam-averaged unit-peak Gaussian filament template for beams (bra, bdec) and trace copy tr."""
    D = np.zeros(len(bra))
    for j in range(len(bra)):
        pr, pd = tp_inv(DX, DY, bra[j], bdec[j])
        d, s, ins = tr.coords(pr, pd, dmax=1.5)
        g = np.where(ins & np.isfinite(d), np.exp(-0.5 * (np.nan_to_num(d, nan=9) / SIGF) ** 2), 0.0)
        D[j] = g.mean()
    return D


ra_E, dec_E = strand('E')
REAL = Trace(ra_E, dec_E)
s_split = float(np.interp(17.6, REAL.dec, REAL.s))


def analyse(tr, Y, need_D=True):
    """On-off paired differences and template regression for a trace copy. Y: (NB, ny) beam quantities."""
    d, s, ins = tr.coords(ra, dec, dmax=4.0)
    near = ins & np.isfinite(d) & (np.abs(d) < 3.5)
    idx = np.where(near)[0]
    res = {}
    on = idx[np.abs(d[idx]) < 0.5]
    off = idx[(np.abs(d[idx]) > 1.5) & (np.abs(d[idx]) < 3.0)]
    diffs, sv = [], []
    for j in on:
        mE = off[(d[off] > 0) & (np.abs(s[off] - s[j]) < 0.75)]
        mW = off[(d[off] < 0) & (np.abs(s[off] - s[j]) < 0.75)]
        parts = [np.nanmean(Y[mm], axis=0) for mm in (mE, mW) if len(mm)]
        if not parts:
            continue
        diffs.append(Y[j] - np.nanmean(parts, axis=0)); sv.append(s[j])
    diffs = np.array(diffs); sv = np.array(sv)
    for name, m in (('south', sv <= s_split), ('north', sv > s_split), ('all', np.ones(len(sv), bool))):
        if len(sv) and m.sum():
            res[f'onoff_{name}'] = np.nanmean(diffs[m], axis=0)
            res[f'n_{name}'] = int(m.sum())
        else:
            res[f'onoff_{name}'] = np.full(Y.shape[1], np.nan); res[f'n_{name}'] = 0
    if need_D:
        D = dilution(tr, ra[idx], dec[idx])
        res['D_on_median'] = float(np.median(D[np.abs(d[idx]) < 0.5])) if (np.abs(d[idx]) < 0.5).any() else np.nan
        for name, m in (('south', s[idx] <= s_split + 0.5), ('north', s[idx] > s_split - 0.5), ('all', np.ones(len(idx), bool))):
            dd, ss, DD = d[idx][m], s[idx][m], D[m]
            X = np.c_[np.ones(m.sum()), dd, ss, dd ** 2, ss ** 2, dd * ss, DD]
            coefs = []
            for q in range(Y.shape[1]):
                y = Y[idx][m, q]
                ok = np.isfinite(y)
                if ok.sum() < 12:
                    coefs.append(np.nan); continue
                cc, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
                coefs.append(cc[-1])
            res[f'reg_{name}'] = np.array(coefs)
            res[f'nreg_{name}'] = int(m.sum())
        res['D'] = D; res['idx'] = idx
    res['d'] = d; res['s'] = s
    return res


Y = MEAS[:, [0, 2, 3, 4, 5, 6, 7]].copy()
YN = ['I_win', 'm1', 'sig', 'I_g', 'v_g', 'fwhm_g', 'I_sum']
Y = np.c_[Y, I_tot, MEAS[:, 8], MEAS[:, 9]]; YN += ['I_tot', 'I_blue', 'I_core']
Y[star] = np.nan            # beams flagged as containing a bright star are excluded
print('star-flagged beams excluded:', int(star.sum()))
REALR = analyse(REAL, Y)
print('real on-off south', dict(zip(YN, np.round(REALR['onoff_south'], 3))), REALR['n_south'])
print('real on-off north', dict(zip(YN, np.round(REALR['onoff_north'], 3))), REALR['n_north'])
print('real reg south', dict(zip(YN, np.round(REALR['reg_south'], 3))))
print('real reg north', dict(zip(YN, np.round(REALR['reg_north'], 3))))
print('D on median', REALR['D_on_median'])

# ------------------------------------------------------------------ nulls
real_band = np.zeros(NB, bool)
real_band[REALR['idx']] = True
cE = SkyCoord(ra_E * u.deg, dec_E * u.deg).galactic
nulls = []
specs = []
for dl in list(np.arange(-66, -7.9, 3.0)) + list(np.arange(8, 66.1, 3.0)):
    for db in (-6.0, -3.0, 0.0, 3.0, 6.0, 9.0):
        specs.append((dl, db))
for dl, db in specs:
    cg = SkyCoord(l=(cE.l.deg + dl) * u.deg, b=(cE.b.deg + db) * u.deg, frame='galactic').icrs
    bc = float(np.mean(cE.b.deg) + db)
    if bc < 12 or bc > 30:
        continue
    for rot, mir in ((0, False), (0, True), (180, False), (180, True), (20, False), (-20, True)):
        tr = Trace(cg.ra.deg, cg.dec.deg, rot=rot, mirror=mir)
        tg = SkyCoord(tr.ra * u.deg, tr.dec * u.deg).galactic
        if tg.b.deg.min() < 2 or tg.b.deg.max() > 45 or tg.l.deg.min() < 130 or tg.l.deg.max() > 280:
            continue
        d, s, ins = tr.coords(ra, dec, dmax=4.0)
        band = ins & (np.abs(d) < 3.5)
        if band.sum() < 40 or (band & real_band).sum() > 0.05 * band.sum():
            continue
        nulls.append(((dl, db, rot, mir), tr))
nulls_all = nulls
# thin: one placement per (dl, db) cell for each orientation is kept; all are used but noted as correlated
print('null copies', len(nulls))
NR = []
for (meta, tr) in nulls:
    r = analyse(tr, Y)
    NR.append((meta, r))
print('nulls analysed')


def collect(key):
    return np.array([r[key] for _, r in NR])


summary = {'columns': YN, 'n_null': len(NR), 'v_helio_oiii_adopted': V_HELIO_OIII, 'half_window_kms': HALF_WIN,
           'D_on_median_real': REALR['D_on_median'], 's_split': s_split}
for key in ['onoff_south', 'onoff_north', 'onoff_all', 'reg_south', 'reg_north', 'reg_all']:
    nv = collect(key)
    summary[key] = {}
    for q, nm in enumerate(YN):
        x = nv[:, q]; x = x[np.isfinite(x)]
        summary[key][nm] = {'real': float(REALR[key][q]), 'null_mean': float(np.mean(x)), 'null_std': float(np.std(x)),
                            'null_p16': float(np.percentile(x, 16)), 'null_p84': float(np.percentile(x, 84)),
                            'z': float((REALR[key][q] - np.mean(x)) / np.std(x)),
                            'frac_null_ge_real': float(np.mean(x >= REALR[key][q])), 'n': int(len(x))}
for nm in ['n_south', 'n_north', 'n_all', 'nreg_south', 'nreg_north']:
    summary[nm] = REALR[nm]

# ------------------------------------------------------------------ injection test (I_peak = 1 R)
D_real = np.zeros(NB); D_real[REALR['idx']] = REALR['D']
fw_inj = 40.0
prof = np.exp(-0.5 * ((VG[None, :] - v_c[:, None]) / (fw_inj / 2.3548)) ** 2)
prof /= (prof.sum(axis=1, keepdims=True) * 2.0)
spec_inj = spec + D_real[:, None] * prof * 1.0
MI = measure(spec_inj, var)
Yi = np.c_[MI[:, [0, 2, 3, 4, 5, 6, 7]], I_tot + D_real * 1.0, MI[:, 8], MI[:, 9]]
Yi[star] = np.nan
RI = analyse(REAL, Yi)
summary['injection_1R'] = {'onoff_south': dict(zip(YN, (RI['onoff_south'] - REALR['onoff_south']).tolist())),
                           'onoff_north': dict(zip(YN, (RI['onoff_north'] - REALR['onoff_north']).tolist())),
                           'reg_south': dict(zip(YN, (RI['reg_south'] - REALR['reg_south']).tolist())),
                           'reg_north': dict(zip(YN, (RI['reg_north'] - REALR['reg_north']).tolist())),
                           'fwhm_kms': fw_inj}
# window capture fraction for different assumed line widths
cap = {}
for fw in (25.0, 32.0, 40.0, 48.0, 60.0):
    sgm = fw / 2.3548
    from scipy.special import erf
    cap[str(fw)] = float(erf(HALF_WIN / (sgm * np.sqrt(2))))
summary['window_capture_fraction'] = cap

# ------------------------------------------------------------------ per-beam table for the real trace
import csv
dR, sR = REALR['d'], REALR['s']
with open(DER / 'verify_hiwham_wham_beams.csv', 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['ra', 'dec', 'l', 'b', 'd_deg', 's_deg', 'D_dilution', 'v_c_lsr', 'I_tot', 'eI_tot', 'star'] + COLS)
    for j in REALR['idx']:
        wr.writerow([f'{ra[j]:.3f}', f'{dec[j]:.3f}', f'{l[j]:.3f}', f'{b[j]:.3f}', f'{dR[j]:.3f}', f'{sR[j]:.3f}',
                     f'{D_real[j]:.4f}', f'{v_c[j]:.2f}', f'{I_tot[j]:.4f}', f'{E_tot[j]:.4f}', int(star[j])] +
                    [f'{x:.4f}' for x in MEAS[j]])
# stacks for figure: on / off mean spectra (south / north), aligned in LSR
idx = REALR['idx']
dd, ss = dR[idx], sR[idx]
stk = {}
for name, m in (('south', ss <= s_split), ('north', ss > s_split)):
    on = idx[m & (np.abs(dd) < 0.5)]
    off = idx[m & (np.abs(dd) > 1.5) & (np.abs(dd) < 3.0)]
    stk[f'on_{name}'] = np.nanmean(spec[on], axis=0); stk[f'off_{name}'] = np.nanmean(spec[off], axis=0)
    stk[f'non_{name}'] = len(on); stk[f'noff_{name}'] = len(off)
np.savez_compressed(DER / 'verify_hiwham_wham_spec.npz', vgrid=VG, D_real=D_real, ra=ra, dec=dec, d=dR, s=sR,
                    I_tot=I_tot, meas=MEAS, v_c=v_c,
                    null_onoff_south=collect('onoff_south'), null_onoff_north=collect('onoff_north'),
                    null_reg_south=collect('reg_south'), null_reg_north=collect('reg_north'),
                    null_meta=np.array([list(m[:2]) + [m[2], int(m[3])] for m, _ in NR], float),
                    cols=np.array(YN), **{k: np.asarray(v) for k, v in stk.items()})
(DER / 'verify_hiwham_wham_results.json').write_text(json.dumps(summary, indent=1))
for key in ['onoff_south', 'onoff_north', 'reg_south', 'reg_north']:
    print(key)
    for nm in YN:
        x = summary[key][nm]
        print(f"   {nm:7s} real {x['real']:+8.3f}  null {x['null_mean']:+8.3f} +- {x['null_std']:.3f}  z {x['z']:+.2f}")
print('injection', json.dumps(summary['injection_1R'], indent=0)[:1500])
