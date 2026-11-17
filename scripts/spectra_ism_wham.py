#!/usr/bin/env python
"""spectra_ism: WHAM Sky Survey DR1 velocity-resolved H-alpha toward the Silk Nebula.

Input (read-only): data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits  (WHAM-SS DR1 v161116,
  2017-09-12 build; teammate download; byte-identical (91,200,960 B, same primary header md5) to the
  file served at https://uwmadison.box.com/shared/static/4kccrw2bgad7muss3z2po7rezklenxxz.fits, the
  remote URL hard-coded in whampy/skySurvey.py, github.com/Deech08/whampy; checked 2026-11-16).
  Teammate NSNS [O III] map (R) data/derived/imaging_optical_nsns_maps.npz, used only to weight
  WHAM beams by their beam-averaged [O III] surface brightness.

Outputs: data/derived/spectra_ism_wham_beams.csv, spectra_ism_wham_results.json,
         spectra_ism_wham_stacks.npz
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_wham.py
"""
import json
import sys

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import (ROOT, DER, REGIONS, vlsr_minus_vhelio, filament_coords,  # noqa: E402
                                spine_arrays)

RNG = np.random.default_rng(1116)
FWHM_IP = 12.0                     # WHAM instrumental FWHM, km/s (Haffner et al. 2003)
SIG_IP = FWHM_IP / 2.3548
VGRID = np.arange(-90.0, 80.01, 2.0)

# ------------------------------------------------------------------------------------------
W = fits.open(ROOT / 'data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits')[1].data
l_all, b_all = W['GAL-LON'].astype(float), W['GAL-LAT'].astype(float)
c_all = SkyCoord(l_all * u.deg, b_all * u.deg, frame='galactic').icrs
ra_all, dec_all = c_all.ra.deg, c_all.dec.deg
sel = np.where((ra_all > 100) & (ra_all < 133) & (dec_all > -4) & (dec_all < 38))[0]
print('WHAM beams in RA 100-133, Dec -4..38:', len(sel))

ra, dec, l, b = ra_all[sel], dec_all[sel], l_all[sel], b_all[sel]
N = len(sel)
spec = np.full((N, len(VGRID)), np.nan)
svar = np.full((N, len(VGRID)), np.nan)
for j, i in enumerate(sel):
    v = W['VELOCITY'][i].astype(float)
    f = W['DATA'][i].astype(float)
    e = W['VARIANCE'][i].astype(float)
    m = np.isfinite(v) & np.isfinite(f) & np.isfinite(e)
    o = np.argsort(v[m])
    vv, ff, ee = v[m][o], f[m][o], e[m][o]
    ok = (VGRID >= vv.min()) & (VGRID <= vv.max())
    spec[j, ok] = np.interp(VGRID[ok], vv, ff)
    svar[j, ok] = np.interp(VGRID[ok], vv, ee)
inten = W['INTEN'][sel].astype(float)
ierr = W['ERROR'][sel].astype(float)
star = W['STAR'][sel].astype(bool)

# Filament coordinates
dperp, dec_sp, beyond = filament_coords(ra, dec)

# Beam-averaged NSNS [O III] (teammate map, Rayleigh) over the 1-deg-diameter beam
nm = np.load(ROOT / 'data/derived/imaging_optical_nsns_maps.npz')
mra, mdec, moiii = nm['ra'], nm['dec'], nm['oiii']
mcoord = SkyCoord(mra.ravel() * u.deg, mdec.ravel() * u.deg)
o_beam = np.full(N, np.nan)
o_cov = np.zeros(N)
inside_map = (ra > mra.min() + 0.5) & (ra < mra.max() - 0.5) & (dec > mdec.min() + 0.5) & (dec < mdec.max() - 0.5)
for j in np.where(inside_map)[0]:
    # quick pre-cut
    pc = (np.abs(mdec.ravel() - dec[j]) < 0.55) & (np.abs((mra.ravel() - ra[j]) * np.cos(np.deg2rad(dec[j]))) < 0.55)
    s = SkyCoord(ra[j] * u.deg, dec[j] * u.deg).separation(mcoord[pc]).deg < 0.5
    vals = moiii.ravel()[pc][s]
    good = np.isfinite(vals)
    o_cov[j] = good.mean() if len(vals) else 0
    if good.sum() > 50:
        o_beam[j] = np.nanmean(np.clip(vals[good], -1, 20))   # clip residual stars/artefacts
print('beams with NSNS [O III] beam average:', np.isfinite(o_beam).sum())


# ------------------------------------------------------------------------------------------
# Multi-Gaussian decomposition of every spectrum (observed widths; intrinsic = sqrt(obs^2 - IP^2))
# ------------------------------------------------------------------------------------------
def gsum(x, *p):
    y = np.full_like(x, p[0])
    for k in range((len(p) - 1) // 3):
        a, x0, s = p[1 + 3 * k:4 + 3 * k]
        y = y + a * np.exp(-0.5 * ((x - x0) / s) ** 2)
    return y


def fit_spec(v, f, e, nmax=3):
    m = np.isfinite(f) & np.isfinite(e) & (v > -85) & (v < 75)
    v, f, sig = v[m], f[m], np.sqrt(np.clip(e[m], 1e-8, None))
    best = None
    for n in range(1, nmax + 1):
        p0 = [0.0]
        lo, hi = [-0.1], [0.1]
        guesses = [(f.max(), v[np.argmax(f)], 12.0), (0.3 * f.max(), -30.0, 15.0), (0.3 * f.max(), 25.0, 12.0)]
        for k in range(n):
            a, x0, s = guesses[k]
            p0 += [max(a, 1e-3), x0, s]
            lo += [0.0, -80.0, SIG_IP]
            hi += [5.0, 70.0, 45.0]
        try:
            p, cov = curve_fit(gsum, v, f, p0=p0, sigma=sig, bounds=(lo, hi), maxfev=20000)
        except Exception:
            continue
        chi2 = float(np.sum(((f - gsum(v, *p)) / sig) ** 2))
        bic = chi2 + len(p) * np.log(len(v))
        perr = np.sqrt(np.clip(np.diag(cov), 0, None))
        if best is None or bic < best['bic'] - 10:
            best = dict(n=n, p=p, perr=perr, chi2=chi2, dof=len(v) - len(p), bic=bic)
    return best


rows = []
for j in range(N):
    fr = fit_spec(VGRID, spec[j], svar[j])
    comps = []
    if fr is not None:
        for k in range(fr['n']):
            a, x0, s = fr['p'][1 + 3 * k:4 + 3 * k]
            ea, ex0, es = fr['perr'][1 + 3 * k:4 + 3 * k]
            area = a * s * np.sqrt(2 * np.pi)
            earea = area * np.hypot(ea / max(a, 1e-9), es / s)
            fw = 2.3548 * s
            comps.append(dict(I_R=area, eI_R=earea, v=x0, ev=ex0, fwhm_obs=fw,
                              fwhm_int=np.sqrt(max(fw ** 2 - FWHM_IP ** 2, 0.0))))
        comps.sort(key=lambda c: -c['I_R'])
    mm = np.isfinite(spec[j]) & (VGRID > -40) & (VGRID < 50)
    vbar = np.nansum(spec[j, mm] * VGRID[mm]) / np.nansum(spec[j, mm])

    def win(v1, v2):
        mw = np.isfinite(spec[j]) & (VGRID >= v1) & (VGRID < v2)
        return float(np.sum(spec[j, mw]) * 2.0), float(np.sqrt(np.sum(svar[j, mw])) * 2.0)

    rows.append(dict(ra=ra[j], dec=dec[j], l=l[j], b=b[j], dperp=dperp[j], dec_spine=dec_sp[j], beyond=bool(beyond[j]),
                     I_tot=inten[j], eI_tot=ierr[j], star=bool(star[j]), oiii_beam_R=o_beam[j], vbar=vbar,
                     ncomp=fr['n'] if fr else 0, chi2nu=fr['chi2'] / fr['dof'] if fr else np.nan, comps=comps,
                     w_m60_m20=win(-60, -20), w_m20_m5=win(-20, -5), w_m5_p9=win(-5, 9), w_p9_p20=win(9, 20),
                     w_p20_p60=win(20, 60)))
print('fitted', len(rows))

# ------------------------------------------------------------------------------------------
# Beams nearest P1-P4 and the paper's velocities in the LSR frame
# ------------------------------------------------------------------------------------------
results = {'input': 'WHAM-SS DR1 v161116 (Haffner et al. 2003, 2016-11-16 build, 2017-09-12 FITS)',
           'n_beams_region': int(N), 'ip_fwhm_kms': FWHM_IP}
dv = {k: float(vlsr_minus_vhelio(*REGIONS[k])[0]) for k in REGIONS}
results['vlsr_minus_vhelio_std_kms'] = dv
paper_vhelio = {'P1': {'OIII': 26.8, 'SII': 11.8, 'NII': 11.8, 'Ha': 16.3},
                'P2': {'OIII': 26.8, 'SII': 12.0, 'NII': 12.0, 'Ha': 17.0},
                'P3': {'OIII': 23.8, 'SII': 11.7, 'NII': 11.7, 'Ha': 18.4},
                'P4': {'SII': 12.7, 'NII': 13.3, 'Ha': 19.7}}
paper_vlsr_table = {'P1': {'OIII': 15.5, 'SII': -1.0, 'NII': -0.9, 'Ha': 4.7},
                    'P2': {'OIII': 14.5, 'SII': 0.5, 'NII': -0.5, 'Ha': 4.0},
                    'P3': {'OIII': 11.2, 'SII': 1.7, 'NII': 0.4, 'Ha': 6.1},
                    'P4': {'SII': 1.1, 'NII': 1.0, 'Ha': 4.5}}
results['paper_vhelio_table4'] = paper_vhelio
results['paper_vlsr_table4'] = paper_vlsr_table
results['vlsr_recomputed_from_vhelio'] = {k: {ln: round(v + dv[k], 1) for ln, v in d.items()} for k, d in paper_vhelio.items()}

cw = SkyCoord(ra * u.deg, dec * u.deg)
near = {}
for k, (r0, d0) in REGIONS.items():
    sep = cw.separation(SkyCoord(r0 * u.deg, d0 * u.deg)).deg
    idx = np.argsort(sep)[:3]
    near[k] = []
    for j in idx:
        rr = rows[j]
        near[k].append({'sep_deg': round(float(sep[j]), 3), 'l': round(rr['l'], 2), 'b': round(rr['b'], 2),
                        'ra': round(rr['ra'], 3), 'dec': round(rr['dec'], 3), 'dperp_deg': round(float(rr['dperp']), 2),
                        'I_tot_R': round(rr['I_tot'], 3), 'eI_tot_R': round(rr['eI_tot'], 3),
                        'vbar_lsr': round(float(rr['vbar']), 2), 'ncomp': rr['ncomp'],
                        'chi2nu': round(float(rr['chi2nu']), 2),
                        'oiii_beam_R_nsns': None if not np.isfinite(rr['oiii_beam_R']) else round(float(rr['oiii_beam_R']), 3),
                        'components': [{kk: round(float(vv), 2) for kk, vv in c.items()} for c in rr['comps']],
                        'I_vlsr_p9_p20_R': [round(x, 3) for x in rr['w_p9_p20']],
                        'I_vlsr_m5_p9_R': [round(x, 3) for x in rr['w_m5_p9']]})
results['nearest_beams'] = near

# Mean spectrum over the beams around each region (within 0.8 deg) for plotting
regspec = {}
for k, (r0, d0) in REGIONS.items():
    sep = cw.separation(SkyCoord(r0 * u.deg, d0 * u.deg)).deg
    j = np.argmin(sep)
    regspec[k] = (spec[j], np.sqrt(svar[j]))

# ------------------------------------------------------------------------------------------
# On-filament minus off-filament stacks (matched in Dec), and [O III]-weighted regression
# ------------------------------------------------------------------------------------------
good = (~beyond) & (dec_sp > 12.0) & (dec_sp < 28.0) & np.all(np.isfinite(spec[:, (VGRID > -70) & (VGRID < 60)]), axis=1)
ON_MAX = 0.5
on = good & (np.abs(dperp) < ON_MAX)
off = good & (np.abs(dperp) > 1.0) & (np.abs(dperp) < 3.0)
print('on beams', on.sum(), 'off beams', off.sum())
vm = (VGRID > -70) & (VGRID < 60)
diffs, ons, offs, decs_on, dO = [], [], [], [], []
for j in np.where(on)[0]:
    nb = off & (np.abs(dec - dec[j]) < 0.8)
    east = nb & (dperp > 0)
    west = nb & (dperp < 0)
    if east.sum() < 1 or west.sum() < 1:
        continue
    ref = 0.5 * (np.nanmean(spec[east], 0) + np.nanmean(spec[west], 0))  # symmetric E/W background
    diffs.append(spec[j] - ref)
    ons.append(spec[j])
    offs.append(ref)
    decs_on.append(dec[j])
    dO.append(o_beam[j] - 0.5 * (np.nanmean(o_beam[east]) + np.nanmean(o_beam[west])))
diffs, ons, offs = np.array(diffs), np.array(ons), np.array(offs)
print('matched on/off pairs', len(diffs))
d_mean = np.nanmean(diffs, 0)
# error from bootstrap over pairs
boots = np.array([np.nanmean(diffs[RNG.integers(0, len(diffs), len(diffs))], 0) for _ in range(2000)])
d_err = np.nanstd(boots, 0)


def gauss_fixed(v, a, v0, s):
    return a * np.exp(-0.5 * ((v - v0) / s) ** 2)


dO = np.array(dO)
stack = {'on_max_deg': ON_MAX, 'off_range_deg': [1.0, 3.0], 'n_pairs': int(len(diffs)),
         'delta_oiii_beam_R_nsns': [float(np.nanmean(dO)), float(np.nanstd(dO) / np.sqrt(np.isfinite(dO).sum())),
                                    int(np.isfinite(dO).sum())], 'dec_on_range': [float(np.min(decs_on)), float(np.max(decs_on))],
         'on_mean_I_R': float(np.nansum(np.nanmean(ons, 0)[vm]) * 2), 'off_mean_I_R': float(np.nansum(np.nanmean(offs, 0)[vm]) * 2)}
# window integrals of the difference spectrum with bootstrap errors
for nm_, (v1, v2) in {'m70_p60': (-70, 60), 'm5_p9': (-5, 9), 'p9_p20': (9, 20), 'p5_p25': (5, 25), 'm20_m5': (-20, -5),
                      'p20_p40': (20, 40)}.items():
    mw = (VGRID >= v1) & (VGRID < v2)
    stack[f'diff_I_{nm_}_R'] = [float(np.nansum(d_mean[mw]) * 2), float(np.nanstd(np.nansum(boots[:, mw], 1) * 2))]
# fixed-velocity Gaussian at the [O III] V_LSR (+13.6, mean of P1-P3 recomputed) with sigma from FWHM_obs
v_oiii = float(np.mean([results['vlsr_recomputed_from_vhelio'][k]['OIII'] for k in ('P1', 'P2', 'P3')]))
fw_oiii_obs = np.sqrt(46.3 ** 2 + FWHM_IP ** 2)    # [O III] intrinsic FWHM ~46 km/s (mean P1-P3) as seen by WHAM
# H-alpha from the same gas is thermally broader than [O III]: FWHM_Ha^2 = FWHM_O^2 + 21.4^2 (T/1e4)(1-1/16)
fw_ha_obs = np.sqrt(46.3 ** 2 + 21.4 ** 2 * (1 - 1 / 16.) + FWHM_IP ** 2)
for tag, fw in (('fwhm_oiii_like', fw_oiii_obs), ('fwhm_ha_from_same_gas_T1e4', fw_ha_obs), ('fwhm_25', np.hypot(25, FWHM_IP))):
    s = fw / 2.3548
    prof = np.exp(-0.5 * ((VGRID - v_oiii) / s) ** 2)
    mfit = vm & np.isfinite(d_mean)
    w = 1 / d_err[mfit] ** 2
    a = np.sum(w * prof[mfit] * d_mean[mfit]) / np.sum(w * prof[mfit] ** 2)
    ea = 1 / np.sqrt(np.sum(w * prof[mfit] ** 2))
    # bootstrap the amplitude too (channels are correlated)
    ab = np.array([np.sum(w * prof[mfit] * bb[mfit]) / np.sum(w * prof[mfit] ** 2) for bb in boots])
    area = a * s * np.sqrt(2 * np.pi)
    stack[f'fixed_v{v_oiii:.1f}_{tag}'] = {'fwhm_obs': round(float(fw), 1), 'I_R': round(float(area), 4),
                                           'eI_R_formal': round(float(ea * s * np.sqrt(2 * np.pi)), 4),
                                           'eI_R_boot': round(float(np.std(ab) * s * np.sqrt(2 * np.pi)), 4)}
# free single-Gaussian fit to the difference
try:
    p, cov = curve_fit(gauss_fixed, VGRID[vm], d_mean[vm], p0=[0.01, 10, 15], sigma=d_err[vm],
                       bounds=([-1, -60, SIG_IP], [1, 60, 40]), maxfev=20000)
    pe = np.sqrt(np.diag(cov))
    stack['free_gauss'] = {'amp': float(p[0]), 'v': float(p[1]), 'ev': float(pe[1]), 'fwhm_obs': float(2.3548 * p[2]),
                           'I_R': float(p[0] * p[2] * np.sqrt(2 * np.pi)),
                           'eI_R': float(np.hypot(pe[0] / max(abs(p[0]), 1e-9), pe[2] / p[2]) * abs(p[0] * p[2]) * np.sqrt(2 * np.pi))}
except Exception as ex:
    stack['free_gauss'] = repr(ex)
results['onoff_stack'] = stack
results['v_oiii_lsr_mean_P1P3'] = v_oiii

# Cross-filament profiles: residual window intensity (beam minus local E+W background at similar Dec)
prof = {}
bins = np.arange(-3.0, 3.01, 0.5)
wins = {'m60_m20': (-60, -20), 'm20_m5': (-20, -5), 'm5_p9': (-5, 9), 'p9_p20': (9, 20), 'p20_p60': (20, 60), 'tot': (-70, 60)}
for nm_, (v1, v2) in wins.items():
    mw = (VGRID >= v1) & (VGRID < v2)
    Iw = np.sum(spec[:, mw], 1) * 2.0
    resid = np.full(N, np.nan)
    for j in np.where(good)[0]:
        nb = off & (np.abs(dec - dec[j]) < 0.8)
        if (nb & (dperp > 0)).sum() and (nb & (dperp < 0)).sum():
            resid[j] = Iw[j] - 0.5 * (np.mean(Iw[nb & (dperp > 0)]) + np.mean(Iw[nb & (dperp < 0)]))
    mu, er, nn = [], [], []
    for b1, b2 in zip(bins[:-1], bins[1:]):
        m_ = good & np.isfinite(resid) & (dperp >= b1) & (dperp < b2)
        mu.append(float(np.mean(resid[m_])) if m_.sum() else np.nan)
        er.append(float(np.std(resid[m_], ddof=1) / np.sqrt(m_.sum())) if m_.sum() > 1 else np.nan)
        nn.append(int(m_.sum()))
    prof[nm_] = {'mean_R': mu, 'sem_R': er, 'n': nn}
Ob = np.where(np.isfinite(o_beam), o_beam, np.nan)
mu, nn = [], []
for b1, b2 in zip(bins[:-1], bins[1:]):
    m_ = good & np.isfinite(Ob) & (dperp >= b1) & (dperp < b2)
    mu.append(float(np.mean(Ob[m_])) if m_.sum() else np.nan)
    nn.append(int(m_.sum()))
prof['oiii_beam_nsns'] = {'mean_R': mu, 'n': nn}
prof['bin_edges_deg'] = bins.tolist()
results['cross_profiles'] = prof

# Regression of each velocity channel on beam-averaged NSNS [O III], with a quadratic surface in
# (RA, Dec) for the smooth DIG background; region restricted to the NSNS map footprint.
reg = np.isfinite(o_beam) & np.all(np.isfinite(spec[:, vm]), axis=1) & (o_cov > 0.9)
# exclude beams overlapping the northern H-alpha patch (teammate polygon) to avoid a spurious correlation
_, _, tr = spine_arrays()
patch = np.asarray(tr['halpha_patch']) if tr is not None and 'halpha_patch' in tr else None
if patch is not None:
    pc_ = SkyCoord(patch[:, 0] * u.deg, patch[:, 1] * u.deg)
    near_patch = np.array([SkyCoord(ra[j] * u.deg, dec[j] * u.deg).separation(pc_).deg.min() < 0.6 for j in range(N)])
else:
    near_patch = np.zeros(N, bool)
regres = {}
for tag, msk in (('all', reg), ('no_Ha_patch', reg & ~near_patch)):
    x = (ra[msk] - 116.5) * np.cos(np.deg2rad(dec[msk]))
    y = dec[msk] - 18.5
    O = o_beam[msk]
    A = np.c_[np.ones_like(x), x, y, x * x, x * y, y * y, O]
    Y = spec[msk][:, vm]
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    # bootstrap over beams
    bc = []
    for _ in range(1000):
        ii = RNG.integers(0, msk.sum(), msk.sum())
        cc, *_ = np.linalg.lstsq(A[ii], Y[ii], rcond=None)
        bc.append(cc[-1])
    bc = np.array(bc)
    slope = coef[-1]                       # R/(km/s) of H-alpha per R of beam-averaged [O III]
    eslope = np.std(bc, 0)
    vv = VGRID[vm]
    out = {'n_beams': int(msk.sum()), 'O_beam_range_R': [float(O.min()), float(O.max())],
           'O_beam_p50_p95': [float(np.percentile(O, 50)), float(np.percentile(O, 95))]}
    for nm_, (v1, v2) in {'p5_p25': (5, 25), 'm5_p9': (-5, 9), 'p9_p20': (9, 20), 'm70_p60': (-70, 60), 'm20_m5': (-20, -5)}.items():
        mw = (vv >= v1) & (vv < v2)
        out[f'Ha_per_OIII_{nm_}'] = [float(np.sum(slope[mw]) * 2), float(np.std(np.sum(bc[:, mw], 1) * 2))]
    regres[tag] = out
    regres[tag + '_slope_spectrum'] = None
    if tag == 'no_Ha_patch':
        slope_keep, eslope_keep = slope, eslope
results['oiii_regression'] = regres
results['n_beams_near_Ha_patch'] = int(near_patch.sum())

# Maps of mean V_LSR and channel intensities for the region grid are in the CSV

# ------------------------------------------------------------------------------------------
# LAMOST relative red-line intensities -> rough absolute (R) using WHAM H-alpha
# ------------------------------------------------------------------------------------------
T2 = {'P1': dict(Ha=(0.30, 0.02), N6583=(0.24, 0.02), S6716=(0.14, 0.02), S6731=(0.13, 0.02)),
      'P2': dict(Ha=(0.37, 0.02), N6583=(0.26, 0.03), S6716=(0.17, 0.02), S6731=(0.15, 0.02)),
      'P3': dict(Ha=(0.26, 0.02), N6583=(0.20, 0.02), S6716=(0.11, 0.02), S6731=(0.10, 0.02)),
      'P4': dict(Ha=(0.25, 0.02), N6583=(0.20, 0.02), S6716=(0.10, 0.02), S6731=(0.10, 0.02))}


def interp_I(r0, d0, sigma=0.6):
    """Gaussian-weighted (sigma deg) WHAM total intensity at a position, from beams within 1.5 deg."""
    sep = cw.separation(SkyCoord(r0 * u.deg, d0 * u.deg)).deg
    m = (sep < 1.5) & ~star
    w = np.exp(-0.5 * (sep[m] / sigma) ** 2)
    return float(np.sum(w * inten[m]) / np.sum(w)), float(np.sqrt(np.sum(w ** 2 * ierr[m] ** 2)) / np.sum(w)), \
        float(np.std(inten[m]))


absol = {}
kvals = []
for k in REGIONS:
    I, eI, scat = interp_I(*REGIONS[k])
    absol[k] = {'WHAM_Ha_interp_R': round(I, 3), 'err_R': round(eI, 3), 'scatter_within_1.5deg_R': round(scat, 3)}
    kvals.append(I / T2[k]['Ha'][0])
kvals = np.array(kvals)
kc, ke = float(np.mean(kvals)), float(np.std(kvals, ddof=1))
results['lamost_to_R'] = {'per_region_R_per_unit': kvals.round(3).tolist(), 'mean': round(kc, 3), 'std': round(ke, 3),
                          'note': 'k = WHAM H-alpha (beam-interpolated, all velocities) / LAMOST relative H-alpha (sky-line units); '
                                  'assumes the 6553 A sky-line normalisation is the same for all fibres of the plate and that '
                                  'H-alpha has little structure inside the 1-deg WHAM beam.'}
for k in REGIONS:
    d = {}
    for ln, (val, e) in T2[k].items():
        d[ln + '_R'] = [round(val * kc, 3), round(float(np.hypot(e * kc, val * ke)), 3)]
    d['SII_sum_R'] = [round((T2[k]['S6716'][0] + T2[k]['S6731'][0]) * kc, 3), None]
    d['WHAM'] = absol[k]
    results.setdefault('lamost_absolute', {})[k] = d

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------
import csv  # noqa: E402
with open(DER / 'spectra_ism_wham_beams.csv', 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['ra', 'dec', 'l', 'b', 'dperp_deg', 'dec_spine', 'beyond_spine', 'I_tot_R', 'eI_tot_R', 'star', 'oiii_beam_R_nsns',
                 'vbar_lsr_m40_p50', 'ncomp', 'chi2nu',
                 'c1_I', 'c1_eI', 'c1_v', 'c1_ev', 'c1_fwhm_int', 'c2_I', 'c2_eI', 'c2_v', 'c2_ev', 'c2_fwhm_int',
                 'c3_I', 'c3_eI', 'c3_v', 'c3_ev', 'c3_fwhm_int',
                 'I_m60_m20', 'I_m20_m5', 'I_m5_p9', 'I_p9_p20', 'I_p20_p60', 'e_I_p9_p20'])
    for r in rows:
        cc = []
        for k in range(3):
            if k < len(r['comps']):
                c = r['comps'][k]
                cc += [round(c['I_R'], 4), round(c['eI_R'], 4), round(c['v'], 2), round(c['ev'], 2), round(c['fwhm_int'], 1)]
            else:
                cc += [''] * 5
        wr.writerow([round(r['ra'], 4), round(r['dec'], 4), round(r['l'], 3), round(r['b'], 3), round(float(r['dperp']), 3),
                     round(float(r['dec_spine']), 3), int(r['beyond']), round(r['I_tot'], 4), round(r['eI_tot'], 4), int(r['star']),
                     '' if not np.isfinite(r['oiii_beam_R']) else round(float(r['oiii_beam_R']), 4),
                     round(float(r['vbar']), 2), r['ncomp'], round(float(r['chi2nu']), 2), *cc,
                     *[round(r[w][0], 4) for w in ('w_m60_m20', 'w_m20_m5', 'w_m5_p9', 'w_p9_p20', 'w_p20_p60')],
                     round(r['w_p9_p20'][1], 4)])
np.savez(DER / 'spectra_ism_wham_stacks.npz', vgrid=VGRID, diff_mean=d_mean, diff_err=d_err,
         on_mean=np.nanmean(ons, 0), off_mean=np.nanmean(offs, 0), slope=slope_keep, slope_err=eslope_keep,
         slope_v=VGRID[vm], **{f'spec_{k}': regspec[k][0] for k in regspec}, **{f'err_{k}': regspec[k][1] for k in regspec},
         ra=ra, dec=dec, spec=spec, dperp=dperp, oiii_beam=o_beam)
json.dump(results, open(DER / 'spectra_ism_wham_results.json', 'w'), indent=1, default=float)
print(json.dumps({k: results[k] for k in ('vlsr_minus_vhelio_std_kms', 'vlsr_recomputed_from_vhelio', 'onoff_stack',
                                          'oiii_regression', 'lamost_to_R')}, indent=1, default=float))
for k in near:
    print(k, near[k][0])
