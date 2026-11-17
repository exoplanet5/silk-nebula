"""verify_sdss: shared helpers (own code; independent of the wave-1 spectra_sdss_lamost_* scripts).

Conventions: vacuum wavelengths; intensities in Rayleigh (photon units); 1 R = 1e6/(4 pi) ph s^-1 cm^-2 sr^-1.
Fibre diameters: SDSS-I/II 3.0", BOSS (incl. SDSS-V plates) 2.0".
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
C = 2.99792458e5           # km/s
H_ERG_S = 6.62607015e-27
C_AS = 2.99792458e18      # A/s
ARCSEC2_SR = (np.pi / 648000.0) ** 2


def air_to_vac(lam_air):
    """IAU standard (Morton 2000) air->vacuum, Angstrom."""
    lam = np.asarray(lam_air, float)
    s2 = (1e4 / lam) ** 2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return lam * n


# vacuum rest wavelengths (air values from NIST / standard nebular lists, converted)
AIR = {'OII3726': 3726.032, 'OII3729': 3728.815, 'Hg4358': 4358.328, 'OIII4363': 4363.210, 'Hgam': 4340.471,
       'HeII4686': 4685.710, 'Hb': 4861.333, 'OIII4959': 4958.911, 'OIII5007': 5006.843,
       'NI5198': 5197.902, 'NI5200': 5200.257, 'OI5577': 5577.339, 'OI6300': 6300.304, 'OI6364': 6363.776,
       'NII6548': 6548.050, 'Ha': 6562.819, 'NII6583': 6583.460, 'SII6716': 6716.440, 'SII6731': 6730.810,
       'Hg5461': 5460.735, 'NaD1': 5895.924}
VAC = {k: float(air_to_vac(v)) for k, v in AIR.items()}
# OH Meinel lines near Halpha (air wavelengths from the UVES sky atlas list, Hanuschik 2003; strongest ones only)
OH_AIR_HA = [6533.04, 6544.03, 6553.62, 6568.78, 6577.18, 6596.64, 6603.99]
OH_VAC_HA = [float(air_to_vac(x)) for x in OH_AIR_HA]


def r_per_flux(lam, diam):
    """Line flux per fibre (erg/s/cm^2) from 1 R filling a fibre of diameter diam arcsec."""
    omega = np.pi * (diam / 2.0) ** 2 * ARCSEC2_SR
    return 1e6 / (4 * np.pi) * omega * H_ERG_S * C_AS / lam


def flam_to_R_per_A(flam_1e17, lam, diam):
    """f_lambda per fibre in 1e-17 erg/s/cm^2/A -> R per A (as calibrated, no diffuse correction)."""
    return np.asarray(flam_1e17) * 1e-17 / r_per_flux(lam, diam)


def load_windows():
    z = np.load(ROOT / 'data/derived/verify_sdss_windows.npz', allow_pickle=False)
    D = {k: z[k] for k in z.files}
    D['lam'] = 10 ** (float(D['ll0']) + float(D['dll']) * D['kidx'])
    return D


# ------------------------------------------------------------------ geometry
def unitvec(ra, dec):
    ra = np.radians(ra); dec = np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def densify(pts, step=0.005):
    pts = np.asarray(pts, float)
    out = []
    for (r1, d1), (r2, d2) in zip(pts[:-1], pts[1:]):
        n = max(2, int(np.hypot((r2 - r1) * np.cos(np.radians(0.5 * (d1 + d2))), d2 - d1) / step) + 1)
        t = np.linspace(0, 1, n, endpoint=False)
        out.append(np.c_[r1 + t * (r2 - r1), d1 + t * (d2 - d1)])
    out.append(pts[-1:])
    return np.vstack(out)


def load_strands():
    tr = json.loads((ROOT / 'data/derived/silk_trace.json').read_text())
    S = {}
    for k, v in tr['strands'].items():
        S[k] = densify([(p['ra'], p['dec']) for p in v['points']])
    patch = np.array(tr['halpha_patch'], float)
    return S, patch, tr


def strand_distance(ra, dec, poly):
    """Angular distance (deg) from points to a densified polyline; signed + east (larger RA) of the line.
    Also returns whether the nearest point is an end point (beyond the strand end)."""
    v = unitvec(ra, dec)
    P = unitvec(poly[:, 0], poly[:, 1])
    dots = np.clip(v @ P.T, -1, 1)
    j = np.argmax(dots, axis=1)
    d = np.degrees(np.arccos(dots[np.arange(len(ra)), j]))
    sgn = np.sign((np.asarray(ra) - poly[j, 0]))
    end = (j == 0) | (j == len(poly) - 1)
    return d * np.where(sgn == 0, 1, sgn), end, poly[j, 1]


def in_polygon(ra, dec, poly):
    from matplotlib.path import Path as MPath
    return MPath(poly).contains_points(np.c_[ra, dec])


# ------------------------------------------------------------------ fitting
def gauss_unit(lam, mu, sig):
    return np.exp(-0.5 * ((lam - mu) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)


def linfit(A, y, w):
    """Weighted linear least squares. Returns coefficients, covariance, chi2."""
    Aw = A * np.sqrt(w)[:, None]
    yw = y * np.sqrt(w)
    cov = np.linalg.pinv(Aw.T @ Aw)
    c = cov @ (Aw.T @ yw)
    chi2 = float(np.sum(w * (y - A @ c) ** 2))
    return c, cov, chi2


def fit_lines(lam, y, w, lines, v_kms, sig_A, npoly=2, ties=None, extra=None):
    """Fit y = poly(lam) + sum_k a_k G(lam; lam_k(1+v/c), sig_A(lam_k)) with linear amplitudes.
    lines: list of names (in VAC) or (name, lam_vac). ties: dict name -> (base_name, ratio) for tied amplitudes.
    sig_A: callable lam -> sigma in A (LSF), or float.
    Returns dict name -> (amp, err) in units of y*A, plus chi2/dof, model."""
    ties = ties or {}
    ok = (w > 0) & np.isfinite(y)
    lamc = 0.5 * (lam[ok].min() + lam[ok].max()) if ok.any() else lam.mean()
    cols, names = [], []
    for p in range(npoly):
        cols.append(((lam - lamc) / 50.0) ** p); names.append(f'poly{p}')
    lv = {}
    for L in lines:
        nm, l0 = (L, VAC[L]) if isinstance(L, str) else L
        lv[nm] = l0
    free = [n for n in lv if n not in ties]
    for nm in free:
        g = 0
        for n2, l0 in lv.items():
            if n2 == nm:
                fac = 1.0
            elif n2 in ties and ties[n2][0] == nm:
                fac = ties[n2][1]
            else:
                continue
            mu = l0 * (1 + v_kms / C)
            s = sig_A(mu) if callable(sig_A) else sig_A
            g = g + fac * gauss_unit(lam, mu, s)
        cols.append(g); names.append(nm)
    if extra is not None:
        for nm, col in extra.items():
            cols.append(col); names.append(nm)
    A = np.array(cols).T
    c, cov, chi2 = linfit(A[ok], y[ok], w[ok])
    res = {nm: (float(c[i]), float(np.sqrt(max(cov[i, i], 0)))) for i, nm in enumerate(names)}
    res['_chi2'] = chi2; res['_dof'] = int(ok.sum() - len(names))
    res['_model'] = A @ c
    res['_cov'] = cov; res['_names'] = names
    return res
