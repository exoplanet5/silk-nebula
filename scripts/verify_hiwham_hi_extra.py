"""verify_hiwham: HI along the spine (velocity structure), at the W strand, near PGCC G199.09+19.95 / the
H-alpha patch, and the catalogued HI shells GS202+17+019 (Suad+2014) and GS202+19-019 (Ehlerova & Palous 2005).

Data: GALFA-HI DR2 2' mosaic (verify_hiwham_common.galfa_cube), LAB HEALPix 10-km/s slices (Kalberla+2005, Land &
Slosar repackaging; 'cutXXX' = XXX..XXX+10 km/s), shell parameters from data/raw/context_vizier_hi/*.tsv.

A. Along-spine HI spectra: in 0.5-deg along-spine bins, mean T_B(v) on the spine (|d|<0.15) and in the local
   off band (0.4<|d|<1.2, both sides); difference spectra; HI moment-1 and peak velocity of the on-spine
   spectrum vs spine Dec; Gaussian decomposition of the per-segment on-spine spectrum.
B. W strand: cross-strand profiles (0.05-deg bins) around the W strand (Dec 15-19.75) in the heliocentric
   W-strand window (V_helio +16..+21), the E-strand window (+10.5..+16.5) and LSR slices, with random
   placement nulls of the W polyline.
C. PGCC G199.09+19.95 and the H-alpha patch: on (r<0.1 deg / inside the patch polygon) minus annulus
   (0.3-0.8 deg) spectra; narrow-line (CNM) search.
D. HI shells: velocity-radius diagrams around the catalogued centres; rim-minus-centre contrast vs velocity and
   against nulls (same ellipse placed at random LAB positions with similar latitude); Silk spine radius in
   units of the shell radius.
Outputs: data/derived/verify_hiwham_hi_extra.json, verify_hiwham_hi_extra.npz
"""
import json
import sys

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
from matplotlib.path import Path as MPath

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_hiwham_common import (RAW, DER, TRACE, strand, Trace, galfa_cube, nhi_slice, vlsr_minus_vhelio,
                                  tp_fwd)

rng = np.random.default_rng(7)
cube, v, gra, gdec = galfa_cube()
RA, DEC = np.meshgrid(gra, gdec)
pra, pdec = RA.ravel(), DEC.ravel()
C2 = cube.reshape(len(v), -1)
cover = np.isfinite(C2[40])
corr = vlsr_minus_vhelio(pra, pdec)
out = {}
npz = {}

# ======================================================================== A. along-spine spectra
ra_E, dec_E = strand('E')
E = Trace(ra_E, dec_E)
d, s, ins = E.coords(pra, pdec, dmax=1.5)
sb = np.arange(0, E.L + 1e-9, 0.5)
on_sp, off_sp, dec_b, nb_on = [], [], [], []
for i in range(len(sb) - 1):
    m = ins & (s >= sb[i]) & (s < sb[i + 1]) & cover
    on = m & (np.abs(d) < 0.15)
    off = m & (np.abs(d) > 0.4) & (np.abs(d) < 1.2)
    on_sp.append(np.nanmean(C2[:, on], axis=1) if on.sum() > 3 else np.full(len(v), np.nan))
    off_sp.append(np.nanmean(C2[:, off], axis=1) if off.sum() > 3 else np.full(len(v), np.nan))
    dec_b.append(float(np.interp(0.5 * (sb[i] + sb[i + 1]), E.s, E.dec))); nb_on.append(int(on.sum()))
on_sp, off_sp, dec_b = np.array(on_sp), np.array(off_sp), np.array(dec_b)
vlsr_oiii = 13.5 + vlsr_minus_vhelio(np.interp(dec_b, E.dec, E.ra), dec_b)


def moments(spec, lo=-15, hi=30):
    k = (v >= lo) & (v <= hi)
    w = np.clip(spec[:, k], 0, None)
    m1 = np.sum(w * v[k], axis=1) / np.sum(w, axis=1)
    vp = v[k][np.nanargmax(np.nan_to_num(spec[:, k], nan=-9), axis=1)]
    return m1, vp


m1_on, vp_on = moments(on_sp)
m1_off, vp_off = moments(off_sp)
diff = on_sp - off_sp
# significance of difference spectra: per-channel scatter across bins (empirical)
dsd = np.nanstd(diff, axis=0)
# linear fit of on-spine HI moment-1 and peak velocity vs Dec
ok = np.isfinite(m1_on)
pm1 = np.polyfit(dec_b[ok], m1_on[ok], 1, cov=True)
pvp = np.polyfit(dec_b[ok], vp_on[ok], 1, cov=True)
out['A_along_spine'] = {
    'dec_bins': dec_b.round(2).tolist(), 'n_on_pix': nb_on,
    'm1_on_kms': m1_on.round(2).tolist(), 'm1_off_kms': m1_off.round(2).tolist(),
    'vpeak_on_kms': vp_on.round(2).tolist(), 'vpeak_off_kms': vp_off.round(2).tolist(),
    'vlsr_oiii_expected': vlsr_oiii.round(2).tolist(),
    'm1_slope_kms_per_deg': [float(pm1[0][0]), float(np.sqrt(pm1[1][0, 0]))],
    'vpeak_slope_kms_per_deg': [float(pvp[0][0]), float(np.sqrt(pvp[1][0, 0]))],
    'm1_on_minus_off_mean': float(np.nanmean(m1_on - m1_off)), 'm1_on_minus_off_std': float(np.nanstd(m1_on - m1_off)),
}
# Gaussian decomposition of segment (1.5 deg) on-spine spectra (-30..+40), up to 3 comps by BIC


def gs(x, *p):
    y = np.zeros_like(x)
    for k in range(len(p) // 3):
        a, x0, sg = p[3 * k:3 * k + 3]
        y = y + a * np.exp(-0.5 * ((x - x0) / sg) ** 2)
    return y


segs = []
for i in range(0, len(sb) - 1, 3):
    sp = np.nanmean(on_sp[i:i + 3], axis=0)
    if not np.isfinite(sp).all():
        continue
    best = None
    for n in (1, 2, 3):
        p0, lo, hi = [], [], []
        guesses = [(sp.max(), v[np.argmax(sp)], 4.0), (0.3 * sp.max(), 15.0, 5.0), (0.2 * sp.max(), -8.0, 8.0)]
        for k in range(n):
            p0 += list(guesses[k]); lo += [0, -30, 0.8]; hi += [50, 40, 25]
        try:
            p, cv = curve_fit(gs, v, sp, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception:
            continue
        rss = np.sum((sp - gs(v, *p)) ** 2)
        bic = len(v) * np.log(rss / len(v)) + len(p) * np.log(len(v))
        if best is None or bic < best[0] - 6:
            best = (bic, p)
    comps = sorted([dict(T=float(best[1][3 * k]), v=float(best[1][3 * k + 1]), fwhm=float(2.3548 * best[1][3 * k + 2]),
                         NHI=float(1.823e18 * best[1][3 * k] * best[1][3 * k + 2] * np.sqrt(2 * np.pi)))
                    for k in range(len(best[1]) // 3)], key=lambda c: c['v'])
    segs.append({'dec_mid': float(np.mean(dec_b[i:i + 3])), 'comps': comps})
out['A_along_spine']['gauss_segments'] = segs
npz.update(dict(A_v=v, A_on=on_sp, A_off=off_sp, A_dec=dec_b, A_vlsr_oiii=vlsr_oiii))
print('A done: m1 slope', out['A_along_spine']['m1_slope_kms_per_deg'], 'vpeak slope', out['A_along_spine']['vpeak_slope_kms_per_deg'])

# ======================================================================== B. W strand
ra_W, dec_W = strand('W')
Wt = Trace(ra_W, dec_W)


def helio_window(lo_h, hi_h):
    dv = abs(v[1] - v[0]); e0, e1 = v - dv / 2, v + dv / 2
    lo = lo_h + corr; hi = hi_h + corr
    o = np.zeros(len(pra)); bad = np.zeros(len(pra), bool)
    for k in range(len(v)):
        w = np.clip(np.minimum(e1[k], hi) - np.maximum(e0[k], lo), 0, None)
        if np.any(w > 0):
            o += w * np.nan_to_num(C2[k]); bad |= (w > 0) & ~np.isfinite(C2[k])
    o *= 1.823e18; o[bad] = np.nan
    return o


WIN = {'W_helio(+16..+21)': helio_window(16, 21), 'E_helio(+10.5..+16.5)': helio_window(10.5, 16.5),
       'lsr+4..+9': nhi_slice(cube, v, 4, 9).ravel(), 'lsr+11..+16': nhi_slice(cube, v, 11, 16).ravel(),
       'all': nhi_slice(cube, v, -29.9, 40.2).ravel()}
EDG = np.arange(-2.0, 2.0001, 0.05); CEN = 0.5 * (EDG[1:] + EDG[:-1])
SG = np.hypot(0.276, 4 / 60) / 2.3548


def prof(tr, val, pix):
    dd, ss, ii = tr.coords(pra[pix], pdec[pix], dmax=2.0)
    m = ii & np.isfinite(dd) & np.isfinite(val[pix]) & (np.abs(dd) < 2)
    b = ((dd[m] + 2) / 0.05).astype(int)
    n = np.bincount(b, minlength=len(CEN))
    p = np.bincount(b, weights=val[pix][m], minlength=len(CEN)) / np.maximum(n, 1)
    p[n < 4] = np.nan
    return p


def famp(p):
    w = (np.abs(CEN) < 1.0) & np.isfinite(p)
    if w.sum() < 15:
        return np.nan
    X = np.c_[np.ones(w.sum()), CEN[w], np.exp(-0.5 * (CEN[w] / SG) ** 2)]
    c, *_ = np.linalg.lstsq(X, p[w], rcond=None)
    return c[2]


allpix = np.arange(len(pra))
dEr, sEr, iEr = E.coords(pra, pdec, dmax=2.5)
dWr, sWr, iWr = Wt.coords(pra, pdec, dmax=2.5)
realband = (iEr & (np.abs(dEr) < 1.2)) | (iWr & (np.abs(dWr) < 1.2))
resW = {}
nullsW = {k: [] for k in WIN}
cen_ok = np.where(cover)[0]
tries = 0
while len(nullsW['all']) < 600 and tries < 20000:
    tries += 1
    j = cen_ok[rng.integers(len(cen_ok))]
    t = Trace(ra_W, dec_W, rot=rng.uniform(0, 360), mirror=bool(rng.integers(2)), centre=(pra[j], pdec[j]))
    near = (np.abs(pra - pra[j]) * np.cos(np.radians(pdec[j])) < 4) & (np.abs(pdec - pdec[j]) < 4)
    pix = np.where(near)[0]
    dd, ss, ii = t.coords(pra[pix], pdec[pix], dmax=2.0)
    band = ii & (np.abs(dd) < 1.5)
    if band.sum() < 100 or cover[pix][band].mean() < 0.9 or realband[pix][band].mean() > 0.05:
        continue
    for k, val in WIN.items():
        nullsW[k].append(famp(prof(t, val, pix)))
profW = {}
for k, val in WIN.items():
    p = prof(Wt, val, allpix)
    a = famp(p); nl = np.array(nullsW[k]); nl = nl[np.isfinite(nl)]
    sr = 0.5 * (np.percentile(nl, 84) - np.percentile(nl, 16))
    resW[k] = {'A_fil_cm2': float(a), 'null_mean': float(nl.mean()), 'null_std': float(nl.std()), 'null_sig_robust': float(sr),
               'z_robust': float((a - np.median(nl)) / sr), 'UL2_cm2': float(max(a - np.median(nl), 0) + 2 * sr), 'n_null': int(len(nl))}
    profW[k] = p
    print(f'W strand {k:24s} A={a:+.2e}  null {np.median(nl):+.2e} +- {sr:.2e}  z={resW[k]["z_robust"]:+.2f}')
out['B_W_strand'] = resW
npz.update({'B_cen': CEN, **{f'B_prof_{k}': p for k, p in profW.items()}})

# ======================================================================== C. PGCC / H-alpha patch
pg = SkyCoord(l=199.09 * u.deg, b=19.95 * u.deg, frame='galactic').icrs
rpg = SkyCoord(pra * u.deg, pdec * u.deg).separation(pg).deg
poly = np.array(TRACE['halpha_patch'])
inpatch = MPath(poly).contains_points(np.c_[pra, pdec])
pc = TRACE['halpha_patch_info']['moments_halfmax']
rpatch = SkyCoord(pra * u.deg, pdec * u.deg).separation(SkyCoord(pc['centroid_ra'] * u.deg, pc['centroid_dec'] * u.deg)).deg
resC = {'pgcc_radec': [float(pg.ra.deg), float(pg.dec.deg)]}
for lab, on, off in (('PGCC_r<0.1', (rpg < 0.1) & cover, (rpg > 0.3) & (rpg < 0.8) & cover),
                     ('patch_polygon', inpatch & cover, (~inpatch) & (rpatch < 1.2) & (rpatch > 0.5) & cover)):
    s_on = np.nanmean(C2[:, on], axis=1); s_off = np.nanmean(C2[:, off], axis=1)
    sd_off = np.nanstd(C2[:, off], axis=1) / np.sqrt(max(on.sum() / 4.0, 1))   # ~independent 4' beams in on region
    dif = s_on - s_off
    k = (v > -25) & (v < 35)
    ipk = np.argmax(np.where(k, dif, -9)); imn = np.argmin(np.where(k, dif, 9))
    # narrow component fit to difference spectrum
    try:
        p, cv = curve_fit(lambda x, a, x0, sg, c0: a * np.exp(-0.5 * ((x - x0) / sg) ** 2) + c0, v[k], dif[k],
                          p0=[dif[ipk], v[ipk], 2.0, 0], bounds=([-20, -25, 0.6, -5], [20, 35, 15, 5]))
        fit = {'T_K': float(p[0]), 'v': float(p[1]), 'fwhm': float(2.3548 * p[2]), 'err': np.sqrt(np.diag(cv)).tolist()}
    except Exception:
        fit = None
    # N(HI) differences in windows
    win = {}
    for wl, (lo, hi) in {'all': (-29.9, 40.2), 'lsr-2..+4': (-2, 4), 'lsr+4..+9': (4, 9), 'lsr+11..+16': (11, 16),
                         'lsr-15..-5': (-15, -5)}.items():
        kk = (v >= lo) & (v <= hi)
        win[wl] = float(1.823e18 * np.sum(dif[kk]) * abs(v[1] - v[0]))
    resC[lab] = {'n_on': int(on.sum()), 'n_off': int(off.sum()), 'dT_peak_K': float(dif[ipk]), 'v_peak': float(v[ipk]),
                 'dT_min_K': float(dif[imn]), 'v_min': float(v[imn]), 'fit_narrow': fit, 'dNHI_windows': win,
                 'on_peak_T': float(np.nanmax(s_on[k])), 'on_peak_v': float(v[k][np.nanargmax(s_on[k])]),
                 'on_m1': float(np.sum(np.clip(s_on[k], 0, None) * v[k]) / np.sum(np.clip(s_on[k], 0, None))),
                 'off_m1': float(np.sum(np.clip(s_off[k], 0, None) * v[k]) / np.sum(np.clip(s_off[k], 0, None)))}
    npz[f'C_{lab}_on'] = s_on; npz[f'C_{lab}_off'] = s_off
    print('C', lab, json.dumps(resC[lab])[:400])
out['C_pgcc_patch'] = resC

# ======================================================================== D. catalogued HI shells (LAB + GALFA)
labmaps = {}
for lo in (-20, -10, 0, 10, 20, 30):
    fn = RAW / 'lambda_imaging_mw/lab_extract/LAB_Healpix' / (f'LAB_cut{lo:+d}.fits' if lo != 0 else 'LAB_cut+0.fits')
    with fits.open(fn) as h:
        tab = h[1].data
        col = tab.columns.names[0]
        arr = np.asarray(tab[col]).ravel().astype(float)
        nside = h[1].header.get('NSIDE'); order = h[1].header.get('ORDERING', 'RING')
    labmaps[lo] = (arr, nside, order)
print('LAB', {k: (len(a), n, o) for k, (a, n, o) in labmaps.items()})
shells = {
    'GS202+17+019': dict(l=202.4, b=16.9, a=5.7, bb=4.5, phi=29.6, v0=18.6, dv=15.5, lmin=203.3, bmin=19.1, src='Suad+2014 J/A+A/564/A116'),
    'GS202+19-019_sizes_as_diameters': dict(l=202.3, b=18.8, a=3.25, bb=2.75, phi=0.0, v0=-18.8, dv=21.6, src='Ehlerova & Palous 2005 (Delta l, Delta b taken as full sizes)'),
    'GS202+19-019_sizes_as_radii': dict(l=202.3, b=18.8, a=6.5, bb=5.5, phi=0.0, v0=-18.8, dv=21.6, src='Ehlerova & Palous 2005 (Delta l, Delta b taken as semi-axes)'),
}


def ell_rho(lv, bv, sh, phi_from='l'):
    x = (lv - sh['l'] + 180) % 360 - 180
    x = x * np.cos(np.radians(sh['b']))
    y = bv - sh['b']
    ph = np.radians(sh['phi'])
    xr = x * np.cos(ph) + y * np.sin(ph); yr = -x * np.sin(ph) + y * np.cos(ph)
    return np.sqrt((xr / sh['a']) ** 2 + (yr / sh['bb']) ** 2)


spc = SkyCoord(E.ra * u.deg, E.dec * u.deg).galactic
resD = {}
for nm, sh in shells.items():
    rec = {'params': sh}
    rho_sp = ell_rho(spc.l.deg, spc.b.deg, sh)
    rec['spine_rho'] = {'min': float(rho_sp.min()), 'max': float(rho_sp.max()), 'median': float(np.median(rho_sp)),
                        'frac_0.8_1.2': float(np.mean((rho_sp > 0.8) & (rho_sp < 1.2)))}
    # LAB contrast rim(0.8-1.2) - centre(<0.5) per slice, plus nulls with the ellipse shifted in l (same b)
    rimc = {}
    for lo, (arr, nside, order) in labmaps.items():
        ipx = hp.query_disc(nside, hp.ang2vec(sh['l'], sh['b'], lonlat=True), np.radians(2.2 * max(sh['a'], sh['bb'])),
                            nest=(order == 'NESTED'))
        lv, bv = hp.pix2ang(nside, ipx, nest=(order == 'NESTED'), lonlat=True)
        val = arr[ipx]
        rho = ell_rho(lv, bv, sh)
        K = np.nanmean(val[(rho > 0.8) & (rho < 1.2)]) - np.nanmean(val[rho < 0.5])
        nul = []
        for dl in list(range(-80, -14, 6)) + list(range(15, 81, 6)):
            for dbb in (-4, 0, 4):
                sh2 = dict(sh); sh2['l'] = sh['l'] + dl; sh2['b'] = sh['b'] + dbb
                ip2 = hp.query_disc(nside, hp.ang2vec(sh2['l'], sh2['b'], lonlat=True), np.radians(1.3 * max(sh['a'], sh['bb'])),
                                    nest=(order == 'NESTED'))
                l2, b2 = hp.pix2ang(nside, ip2, nest=(order == 'NESTED'), lonlat=True)
                r2 = ell_rho(l2, b2, sh2); v2 = arr[ip2]
                nul.append(np.nanmean(v2[(r2 > 0.8) & (r2 < 1.2)]) - np.nanmean(v2[r2 < 0.5]))
        nul = np.array(nul)
        # normalise by the mean level so that slices are comparable: express as fraction of mean T in rho<1.2
        lvl = np.nanmean(val[rho < 1.2])
        rimc[f'{lo:+d}..{lo + 10:+d}'] = {'K': float(K), 'K_over_mean': float(K / lvl), 'null_median': float(np.median(nul)),
                                           'null_sig_robust': float(0.5 * (np.percentile(nul, 84) - np.percentile(nul, 16))),
                                           'z_robust': float((K - np.median(nul)) / (0.5 * (np.percentile(nul, 84) - np.percentile(nul, 16)))),
                                           'n_null': int(len(nul))}
    rec['LAB_rim_minus_centre'] = rimc
    # GALFA velocity-radius diagram (where covered)
    gc = SkyCoord(l=sh['l'] * u.deg, b=sh['b'] * u.deg, frame='galactic').icrs
    gl_ = SkyCoord(pra * u.deg, pdec * u.deg).galactic
    rho_g = ell_rho(gl_.l.deg, gl_.b.deg, sh)
    redg = np.arange(0, 1.61, 0.1)
    vr = np.full((len(redg) - 1, len(v)), np.nan); fcov = []
    for i in range(len(redg) - 1):
        m = (rho_g >= redg[i]) & (rho_g < redg[i + 1])
        fcov.append(float(cover[m].mean()) if m.any() else 0.0)
        mm = m & cover
        if mm.sum() > 20:
            vr[i] = np.nanmean(C2[:, mm], axis=1)
    rec['galfa_rho_coverage'] = fcov
    rec['centre_radec'] = [float(gc.ra.deg), float(gc.dec.deg)]
    # GALFA contrast at v0 +- 4 and at the revised [O III] velocity
    for wl, (lo, hi) in {'v0+-4': (sh['v0'] - 4, sh['v0'] + 4), 'lsr-2..+4': (-2, 4)}.items():
        kk = (v >= lo) & (v <= hi)
        prof_r = 1.823e18 * np.nansum(vr[:, kk], axis=1) * abs(v[1] - v[0])
        rec[f'galfa_NHI_vs_rho_{wl}'] = prof_r.tolist()
    npz[f'D_{nm}_vr'] = vr
    resD[nm] = rec
    print('D', nm, 'spine rho', rec['spine_rho'], {k: round(x['z_robust'], 2) for k, x in rimc.items()})
npz['D_redges'] = np.arange(0, 1.61, 0.1)
out['D_shells'] = resD
(DER / 'verify_hiwham_hi_extra.json').write_text(json.dumps(out, indent=1, default=float))
np.savez_compressed(DER / 'verify_hiwham_hi_extra.npz', **npz)
print('saved')
