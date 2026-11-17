"""Shared helpers for the 'context' topic: filament curvature, LSR-frame back-tracing, blackbody ionising rates.

Import: sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts'); from context_common import *
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import quad
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import CENTER, REGIONS, ROOT, load_trace, gal

KMS_MYR_TO_PC = 1.0227121650537077          # 1 km/s for 1 Myr = 1.0227 pc
K_PM = 4.740470463533348                    # km/s per (mas/yr * kpc)
SUN_UVW = np.array([11.1, 12.24, 7.25])      # km/s, Schoenrich, Binney & Dehnen (2010)

# Literature X-ray/optical outlines of the Monogem complex (centre RA, Dec, radius deg)
MONOGEM_OUTLINES = {
    'Monogem single (Knies+24 eROSITA)': (107.08, 13.50, 11.97),
    'Monogem two-SNR (Knies+24)': (104.15, 15.30, 11.00),
    'G205.6+12.4 (Knies+24)': (110.33, 12.94, 8.82),
    'Monogem ring on PSR (Thorsett+03)': None,   # 9.2 deg circle centred on PSR B0656+14 (l,b)=(201.1,+8.3)
}
PSR_B0656_RADEC = (104.9506, 14.2393)

# Provisional arc constraints (used only if no refined trace exists):
#   P3, P1, P2 from Wu+2026 Table 1 (sigma 0.15 deg) and the two ends of the NSNS 'O200.0' arc
#   (Ziegenbalg NSNS gem-01 annotation: centre 115.4516,+19.9427, 1015' x 182', PA 173 deg) -> ends at
#   (116.50,+11.54) and (114.28,+28.34), sigma 0.8 deg (bow/width of the 2.5 deg E-W box).
PROVISIONAL_POINTS = [(116.45, 14.72, 0.15, 'P3'), (116.65, 16.60, 0.15, 'P1'), (116.21, 17.74, 0.15, 'P2'),
                      (116.50, 11.54, 0.8, 'NSNS south end'), (114.28, 28.34, 0.8, 'NSNS north end')]


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def radec_of(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360, np.degrees(np.arcsin(v[..., 2]))


def angsep(ra1, dec1, ra2, dec2):
    a, b = unitvec(ra1, dec1), unitvec(ra2, dec2)
    return np.degrees(np.arccos(np.clip(np.sum(a * b, -1), -1, 1)))


def trace_points():
    """Return (ra, dec, sigma, label, source_tag) arrays for curvature fitting."""
    tr = load_trace()
    if tr and not tr.get('preliminary', True) and tr.get('spine'):
        sp = tr['spine']
        return (np.array([p['ra'] for p in sp]), np.array([p['dec'] for p in sp]),
                np.array([p.get('sigma_deg', 0.2) for p in sp]), [p.get('source', '') for p in sp],
                f"silk_trace.json v{tr.get('version')}")
    a = np.array([p[:3] for p in PROVISIONAL_POINTS])
    return a[:, 0], a[:, 1], a[:, 2], [p[3] for p in PROVISIONAL_POINTS], 'provisional (P1-P3 + NSNS arc ends)'


def fit_small_circle(ra, dec, sig, west=True, x0=None):
    """Least-squares small circle on the sphere: minimise (angdist(c, p_i) - rho)/sig_i.
    Returns (ra_c, dec_c, rho). The solution with the centre on the concave (west) side is chosen."""
    P = unitvec(ra, dec)
    m = P.mean(0); m /= np.linalg.norm(m)
    if x0 is None:
        ra0, dec0 = radec_of(m)
        x0 = [ra0 - 20 if west else ra0 + 20, dec0, 20.0]

    def res(x):
        c = unitvec(x[0], x[1])
        d = np.degrees(np.arccos(np.clip(P @ c, -1, 1)))
        return (d - x[2]) / sig
    best = None
    for start in ([x0] + [[x0[0] + dr, x0[1] + dd, rr] for dr in (-30, -10, 10) for dd in (-10, 10) for rr in (8, 30, 60)]):
        try:
            s = least_squares(res, start, bounds=([0, -89, 0.5], [360, 89, 89.5]))
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    return best.x[0], best.x[1], best.x[2], best.cost


def curvature_mc(n=2000, seed=1):
    """Monte-Carlo the small-circle fit; returns dict with nominal and sample centres/radii."""
    ra, dec, sig, lab, src = trace_points()
    rng = np.random.default_rng(seed)
    rc, dc, rho, _ = fit_small_circle(ra, dec, sig)
    samp = []
    cosd = np.cos(np.radians(dec))
    for k in range(n):
        rr = ra + rng.normal(0, sig) / cosd
        dd = dec + rng.normal(0, sig)
        try:
            samp.append(fit_small_circle(rr, dd, sig, x0=[rc, dc, rho])[:3])
        except Exception:
            pass
    samp = np.array(samp)
    return dict(source=src, labels=lab, ra=ra, dec=dec, sig=sig, ra_c=rc, dec_c=dc, rho=rho, samples=samp)


def centre_region_grid(ra_rng=(30, 125), dec_rng=(-15, 40), step=0.5, strand='E'):
    """Delta-chi2 map of small-circle centres for the trace (single-circle model).
    Intrinsic scatter sigma_int is added in quadrature so that chi2_min/dof = 1 at the best centre.
    Returns dict(ra, dec, dchi2 (2-D), radius (2-D), sigma_int, best=(ra,dec,rho), n)."""
    tr = load_trace()
    if tr and not tr.get('preliminary', True):
        pts = tr['strands'][strand]['points'] if strand in tr.get('strands', {}) else tr['spine']
        ra = np.array([p['ra'] for p in pts]); dec = np.array([p['dec'] for p in pts])
        sig = np.array([p.get('sigma_deg', 0.05) for p in pts])
    else:
        ra, dec, sig, _, _ = trace_points()
    P = unitvec(ra, dec)
    gra = np.arange(ra_rng[0], ra_rng[1] + 1e-9, step)
    gdec = np.arange(dec_rng[0], dec_rng[1] + 1e-9, step)
    GR, GD = np.meshgrid(gra, gdec)
    C = unitvec(GR, GD)                                     # (nd, nr, 3)
    D = np.degrees(np.arccos(np.clip(np.einsum('ijk,nk->ijn', C, P), -1, 1)))  # (nd,nr,npts)

    def chi2_for(sint):
        w = 1 / (sig ** 2 + sint ** 2)
        rho = (D * w).sum(-1) / w.sum()
        return (((D - rho[..., None]) ** 2) * w).sum(-1), rho
    sint = 0.0
    for _ in range(30):                                      # iterate sigma_int so chi2_min/dof = 1
        chi2, rho = chi2_for(sint)
        red = chi2.min() / (len(ra) - 3)
        if abs(red - 1) < 0.02:
            break
        resid = D[np.unravel_index(np.argmin(chi2), chi2.shape)] - rho.flat[np.argmin(chi2)]
        sint = max(np.sqrt(max(np.mean(resid ** 2) - np.mean(sig ** 2), 1e-6)), 1e-3) * (red ** 0.5 if sint > 0 else 1)
        sint = float(np.sqrt(max(np.var(resid) * len(ra) / (len(ra) - 3) - np.median(sig ** 2), 1e-6)))
    chi2, rho = chi2_for(sint)
    k = np.unravel_index(np.argmin(chi2), chi2.shape)
    return dict(ra=gra, dec=gdec, dchi2=chi2 - chi2.min(), radius=rho, sigma_int=sint,
                best=(GR[k], GD[k], rho[k]), n=len(ra), chi2min=float(chi2.min()))


def in_centre_region(reg, ra, dec, level=6.18):
    """True where (ra, dec) falls inside the Delta-chi2 <= level region (nearest grid cell)."""
    ra = np.atleast_1d(ra); dec = np.atleast_1d(dec)
    i = np.clip(np.round((dec - reg['dec'][0]) / (reg['dec'][1] - reg['dec'][0])).astype(int), 0, len(reg['dec']) - 1)
    j = np.clip(np.round(((ra % 360) - reg['ra'][0]) / (reg['ra'][1] - reg['ra'][0])).astype(int), 0, len(reg['ra']) - 1)
    inside = reg['dchi2'][i, j] <= level
    outside_grid = (ra % 360 < reg['ra'][0]) | (ra % 360 > reg['ra'][-1]) | (dec < reg['dec'][0]) | (dec > reg['dec'][-1])
    return inside & ~outside_grid


# ---------------------------------------------------------------- kinematics
def galactic_xyz(l, b, d_pc):
    l, b = np.radians(l), np.radians(b)
    return np.stack([d_pc * np.cos(b) * np.cos(l), d_pc * np.cos(b) * np.sin(l), d_pc * np.sin(b)], -1)


def lb_of_xyz(X):
    d = np.linalg.norm(X, axis=-1)
    l = np.degrees(np.arctan2(X[..., 1], X[..., 0])) % 360
    b = np.degrees(np.arcsin(X[..., 2] / d))
    return l, b, d


def helio_velocity(ra, dec, d_pc, pmra, pmdec, vr):
    """Heliocentric Galactic-Cartesian velocity (km/s), arrays broadcast."""
    c = SkyCoord(ra=np.atleast_1d(ra) * u.deg, dec=np.atleast_1d(dec) * u.deg,
                 pm_ra_cosdec=np.atleast_1d(pmra) * u.mas / u.yr, pm_dec=np.atleast_1d(pmdec) * u.mas / u.yr,
                 distance=np.atleast_1d(d_pc) * u.pc, radial_velocity=np.atleast_1d(vr) * u.km / u.s, frame='icrs')
    g = c.galactic
    v = g.velocity
    return np.stack([v.d_x.to_value(u.km / u.s), v.d_y.to_value(u.km / u.s), v.d_z.to_value(u.km / u.s)], -1)


def backtrack_lsr(ra, dec, d_pc, pmra, pmdec, vr, t_myr):
    """Position (relative to the present Sun, in a frame co-moving with the LSR) of the ISM parcel the object
    occupied t_myr ago, i.e. X_now - v_LSR * t. Returns (l, b, d) with shape (n_obj, n_t)."""
    l0, b0 = gal(ra, dec)
    X0 = galactic_xyz(l0, b0, np.atleast_1d(d_pc))
    v = helio_velocity(ra, dec, d_pc, pmra, pmdec, vr) + SUN_UVW
    t = np.atleast_1d(t_myr)
    Xt = X0[:, None, :] - v[:, None, :] * t[None, :, None] * KMS_MYR_TO_PC
    return lb_of_xyz(Xt)


def lsr_pm(ra, dec, d_pc, pmra, pmdec, vr=0.0):
    """Proper motion relative to the LSR (mas/yr, along RA*cosDec, Dec) -- removes the solar reflex."""
    v = helio_velocity(ra, dec, d_pc, pmra, pmdec, vr) + SUN_UVW
    l0, b0 = gal(ra, dec)
    X0 = galactic_xyz(l0, b0, np.atleast_1d(d_pc))
    X1 = X0 + v * 1e-3 * KMS_MYR_TO_PC        # 1 kyr step
    l1, b1, _ = lb_of_xyz(X1)
    c0 = SkyCoord(l0 * u.deg, b0 * u.deg, frame='galactic').icrs
    c1 = SkyCoord(l1 * u.deg, b1 * u.deg, frame='galactic').icrs
    dra = ((c1.ra.deg - c0.ra.deg + 180) % 360 - 180) * np.cos(np.radians(c0.dec.deg)) * 3.6e6 / 1e3
    ddec = (c1.dec.deg - c0.dec.deg) * 3.6e6 / 1e3
    return dra, ddec


# ---------------------------------------------------------------- ionising photons
H, KB, C = 6.62607015e-27, 1.380649e-16, 2.99792458e10
RSUN = 6.957e10
ALPHA_B = 2.59e-13            # cm^3/s, case B, T_e = 1e4 K (Osterbrock & Ferland 2006)


def q_bb(teff, r_rsun, e_ev=13.598):
    """Blackbody photon rate above e_ev (photons/s) for radius r_rsun."""
    x0 = e_ev * 1.602176634e-12 / (KB * teff)
    integ = quad(lambda x: x * x / np.expm1(x), x0, max(60, x0 + 60))[0]
    flux = 2 * np.pi * (KB * teff / H) ** 3 / C ** 2 * integ          # photons/cm^2/s at surface
    return 4 * np.pi * (r_rsun * RSUN) ** 2 * flux


def stromgren_pc(q, n):
    return (3 * q / (4 * np.pi * n * n * ALPHA_B)) ** (1 / 3) / 3.0857e18
