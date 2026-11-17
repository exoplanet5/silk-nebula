#!/usr/bin/env python
"""spectra_ism: measure interstellar Ca II K/H and Na I D1/D2 absorption in archival high-resolution spectra
of early-type stars around the Silk Nebula (files listed in data/raw/spectra_ism_spectra/manifest.csv).

Method (per spectrum and line):
  * heliocentric velocity scale (UVES/ELODIE topocentric + BERV; HARPS/SOPHIE barycentric; frames checked with
    telluric O2 lines, spectra_ism_frame_test.py);
  * resample to a uniform velocity grid; continuum = iterated running median (70 km/s) + 15 km/s boxcar, with
    narrow absorption dips masked, so stellar lines broader than ~50 km/s are divided out;
  * noise = robust (MAD) scatter of the normalised flux at 80-250 km/s from the line;
  * 0-3 Gaussian absorption components (sigma >= instrumental) fitted in V_helio -60..+80 km/s, chosen by BIC;
  * forced fit at V_helio = +25 km/s (the [O III] velocity) with b = 5 km/s to give EW and a 3-sigma limit;
  * column densities in the optically thin (AOD) limit, N = 1.130e17 W(mA) / (f lambda^2) cm^-2.
Outputs: data/derived/spectra_ism_abs_components.csv, spectra_ism_abs_stars.csv, spectra_ism_abs_profiles.npz
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_absorption_measure.py
"""
import csv
import io
import json
import sys

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.optimize import curve_fit
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, filament_coords, vlsr_minus_vhelio  # noqa: E402
from spectra_ism_speclib import load, window, C, LINES, FVAL  # noqa: E402
from spectra_ism_absorption_catalogs_lib import spine_sep  # noqa: E402
import spectra_ism_telluric as tel  # noqa: E402

SP = RAW / 'spectra_ism_spectra'
DV = 0.5
VG = np.arange(-300, 300 + DV / 2, DV)
V_OIII = 25.0          # heliocentric [O III] velocity (Wu+2026 Table 4: +23.8..+26.8)

# star metadata from the SIMBAD-based archive inventory
cand = {r['name'].replace('  ', ' ').strip(): r for r in csv.DictReader(open(DER / 'spectra_ism_archive_candidates.csv'))}
simb = list(csv.DictReader(io.StringIO((RAW / 'spectra_ism_archives/simbad_OBA_V10_box.csv').read_text(errors='replace'))))
srv = {}
for r in simb:
    srv[' '.join(r['main_id'].split())] = r


def meta_for(star, sp):
    key = ' '.join(star.split())
    hit = None
    for k, r in cand.items():
        if ' '.join(k.split()).replace('* ', '') == key or ' '.join(k.split()) == key:
            hit = r
    if hit is None:        # fall back to the closest candidate to the FITS pointing
        c0 = SkyCoord(float(sp['ra']) * u.deg, float(sp['dec']) * u.deg)
        best = min(cand.values(), key=lambda r: c0.separation(SkyCoord(float(r['ra']) * u.deg, float(r['dec']) * u.deg)).deg)
        hit = best
    rv = np.nan
    for k, r in srv.items():
        if ' '.join(k.split()) == ' '.join(hit['name'].split()):
            try:
                rv = float(r['rvz_radvel'])
            except ValueError:
                pass
    return hit, rv


def gneg(v, *p):
    y = np.ones_like(v)
    for k in range(len(p) // 3):
        a, v0, s = p[3 * k:3 * k + 3]
        y = y - a * np.exp(-0.5 * ((v - v0) / s) ** 2)
    return y


def normalise(v, f):
    fv = np.interp(VG, v, f, left=np.nan, right=np.nan)
    good = np.isfinite(fv)
    x = np.where(good, fv, np.nanmedian(fv))
    nwin = int(70 / DV) | 1
    cont = uniform_filter1d(median_filter(x, size=nwin, mode='nearest'), int(15 / DV), mode='nearest')
    for _ in range(3):
        r = x / cont - 1
        sig = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        y = np.where(r < -2.5 * sig, cont, x)
        cont = uniform_filter1d(median_filter(y, size=nwin, mode='nearest'), int(15 / DV), mode='nearest')
    nf = np.where(good, fv / cont, np.nan)
    return nf


def measure(v, nf, lam0, fwhm_ip, line):
    sig_ip = fwhm_ip / 2.3548
    off = (np.abs(VG) > 90) & (np.abs(VG) < 250) & np.isfinite(nf)
    rms = 1.4826 * np.nanmedian(np.abs(nf[off] - np.nanmedian(nf[off])))
    m = (VG > -60) & (VG < 80) & np.isfinite(nf)
    vv, yy = VG[m], nf[m]
    if len(vv) < 100 or not np.isfinite(rms):
        return None
    sig = np.full_like(yy, rms)
    best = dict(n=0, p=[], perr=[], chi2=float(np.sum(((yy - 1) / sig) ** 2)))
    best['bic'] = best['chi2']
    for n in range(1, 4):
        resid = 1 - yy
        p0, lo, hi = [], [], []
        used = []
        for k in range(n):
            r_ = resid.copy()
            for u0 in used:
                r_[np.abs(vv - u0) < 6] = -9
            j = int(np.argmax(uniform_filter1d(r_, 5)))
            used.append(vv[j])
            p0 += [max(resid[j], 3 * rms), vv[j], max(sig_ip * 1.2, 2.5)]
            lo += [0.0, -60.0, sig_ip * 0.8]
            hi += [1.0, 80.0, 25.0]
        try:
            p, cov = curve_fit(gneg, vv, yy, p0=p0, sigma=sig, bounds=(lo, hi), maxfev=20000)
        except Exception:
            continue
        chi2 = float(np.sum(((yy - gneg(vv, *p)) / sig) ** 2))
        bic = chi2 + 3 * n * np.log(len(vv) / (fwhm_ip / DV))   # effective independent pixels ~ per resolution element
        if bic < best['bic'] - 10:
            best = dict(n=n, p=p, perr=np.sqrt(np.clip(np.diag(cov), 0, None)), chi2=chi2, bic=bic)
    comps = []
    for k in range(best['n']):
        a, v0, s = best['p'][3 * k:3 * k + 3]
        ea, ev, es = best['perr'][3 * k:3 * k + 3]
        ew = a * s * np.sqrt(2 * np.pi) * lam0 / C * 1000
        eew = ew * np.hypot(ea / max(a, 1e-9), es / s)
        # noise-based EW error floor: rms * sqrt(Npix_eff) * dlam
        ew_floor = rms * np.sqrt(2.5 * s / DV) * DV * lam0 / C * 1000 * np.sqrt(fwhm_ip / DV / 2.0)
        eew = float(np.hypot(eew, ew_floor))
        b = np.sqrt(max(2 * (s ** 2 - sig_ip ** 2), 0.0))
        comps.append(dict(v_helio=float(v0), ev=float(ev), depth=float(a), fwhm_obs=float(2.3548 * s), b_kms=float(b),
                          EW_mA=float(ew), eEW_mA=eew, snr=float(ew / eew) if eew > 0 else np.nan,
                          broad_flag=bool(s > 9.0)))
    comps.sort(key=lambda c: c['v_helio'])
    # forced component at the [O III] velocity
    s_f = np.hypot(sig_ip, 5.0 / np.sqrt(2))
    prof = np.exp(-0.5 * ((vv - V_OIII) / s_f) ** 2)
    a_f = np.sum(prof * (1 - yy)) / np.sum(prof ** 2)
    ea_f = rms / np.sqrt(np.sum(prof ** 2)) * np.sqrt(fwhm_ip / DV / 2.0)   # correct for correlated pixels
    k_ew = s_f * np.sqrt(2 * np.pi) * lam0 / C * 1000
    # EW integrated in windows
    def wint(v1, v2):
        mm = (VG >= v1) & (VG < v2) & np.isfinite(nf)
        return float(np.sum(1 - nf[mm]) * DV * lam0 / C * 1000), float(rms * np.sqrt(mm.sum()) * DV * lam0 / C * 1000 * np.sqrt(fwhm_ip / DV / 2.0))
    return dict(rms=float(rms), ncomp=best['n'], comps=comps, forced_EW=float(a_f * k_ew), forced_eEW=float(ea_f * k_ew),
                EW_tot_m40_p70=wint(-40, 70), EW_p15_p30=wint(15, 30), EW_hv_neg=wint(-60, -5), EW_hv_pos=wint(45, 80))


def ncol(ew_ma, line):
    return 1.130e17 * ew_ma / (FVAL[line] * LINES[line] ** 2)


man = [m for m in csv.DictReader(open(SP / 'manifest.csv')) if m['ok'] == 'True']
SPEC = {m['file']: load(SP / m['file']) for m in man}
POOL = [sp for sp in SPEC.values() if tel.tau_on_grid(sp) is not None]
tel_rows = []
comp_rows, star_rows, prof_store = [], [], {}
for mrow in man:
    sp = SPEC[mrow['file']]
    hit, rv = meta_for(mrow['star'], sp)
    fwhm_ip = C / sp['R'] if np.isfinite(sp['R']) and sp['R'] > 0 else 7.0
    sp_nat = sp
    corr = tel.correct(sp, POOL) if tel.tau_on_grid(sp) is not None else None
    if corr is not None:
        sp_tc = tel.apply(sp, corr)
        tel_rows.append(dict(file=mrow['file'], alpha=round(float(corr[2]), 3), max_tau_template=round(float(corr[3].max()), 3)))
    for line in LINES:
        spx = sp_tc if (line.startswith('Na') and corr is not None) else sp_nat
        v, f, e = window(spx, line, -320, 320)
        if len(v) < 200 or (v.max() - v.min()) < 500:
            continue
        nf = normalise(v, f)
        res = measure(v, nf, LINES[line], fwhm_ip, line)
        if res is None:
            continue
        key = f"{mrow['star']}|{sp['inst']}|{line}"
        prof_store[key] = nf.astype(np.float32)
        base = dict(star=mrow['star'], simbad=hit['name'].strip(), inst=sp['inst'], file=mrow['file'], date=sp['date'],
                    R=sp['R'], frame=sp['frame'], berv=round(sp['berv'], 2), line=line,
                    telluric_corr=bool(line.startswith('Na') and corr is not None), d_pc=hit['dist_pc'], dperp_deg=hit['dperp_deg'], sep_spine_deg=hit['dmin_spine_deg'],
                    dec_spine=hit['dec_spine'], sp=hit['sp'], rv_star=rv, rms=round(res['rms'], 4))
        for c in res['comps']:
            comp_rows.append(dict(base, **{k: (round(v_, 3) if isinstance(v_, float) else v_) for k, v_ in c.items()},
                                  logN_AOD=round(float(np.log10(max(ncol(c['EW_mA'], line), 1))), 2)))
        star_rows.append(dict(base, ncomp=res['ncomp'],
                              comps=';'.join(f"{c['v_helio']:+.1f}({c['EW_mA']:.1f}+-{c['eEW_mA']:.1f}{'B' if c['broad_flag'] else ''})"
                                             for c in res['comps']),
                              EW_tot_mA=round(res['EW_tot_m40_p70'][0], 2), eEW_tot_mA=round(res['EW_tot_m40_p70'][1], 2),
                              EW_forced_v25_mA=round(res['forced_EW'], 2), eEW_forced_v25_mA=round(res['forced_eEW'], 2),
                              EW_p15_p30_mA=round(res['EW_p15_p30'][0], 2), eEW_p15_p30_mA=round(res['EW_p15_p30'][1], 2),
                              EW_neg_m60_m5_mA=round(res['EW_hv_neg'][0], 2), EW_pos_p45_p80_mA=round(res['EW_hv_pos'][0], 2),
                              eEW_win_mA=round(res['EW_hv_neg'][1], 2)))
        print(f"{mrow['star']:11s} {sp['inst']:6s} {line:5s} d={hit['dist_pc']:>6s} sep={hit['dmin_spine_deg']:>6s} rms={res['rms']:.4f} "
              f"n={res['ncomp']} {star_rows[-1]['comps']:60s} EW25={res['forced_EW']:+.2f}+-{res['forced_eEW']:.2f}")

with open(DER / 'spectra_ism_abs_components.csv', 'w', newline='') as fh:
    keys = list(comp_rows[0].keys())
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(comp_rows)
with open(DER / 'spectra_ism_abs_stars.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(star_rows[0].keys()))
    w.writeheader()
    w.writerows(star_rows)
np.savez_compressed(DER / 'spectra_ism_abs_profiles.npz', vgrid=VG, keys=np.array(list(prof_store)),
                    **{f'p{i}': v for i, v in enumerate(prof_store.values())})
with open(DER / 'spectra_ism_abs_telluric.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(tel_rows[0].keys()))
    w.writeheader()
    w.writerows(tel_rows)
print('components', len(comp_rows), 'star-lines', len(star_rows))
