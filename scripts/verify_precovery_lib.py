"""verify_precovery: independent library for the SDSS imaging-sky precovery check (wave-2 verification).

Written from scratch for the `precovery` verification task; it does NOT import any wave-1 imaging_optical code.
Inputs are read-only:
  data/raw/verify_precovery_sdss/fields_box_full.csv   SkyServer DR18 Field table, all bands and 3 sky estimators
  data/derived/silk_trace.json                         geometry (positions/widths only, no NSNS intensities)
  data/derived/spectra_desi_results.json               DESI petal forward-model ridge amplitudes (R)
  data/derived/spectra_sdss_lamost_sdss_onoff.json     SDSS sky-fibre strand amplitudes (R)
  data/raw/verify_precovery_sdss/filter_curves/*, data/raw/filter_curves/*   SDSS 2010 response curves (speclite)
"""
from pathlib import Path
import json
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw'
DER = ROOT / 'data/derived'
FIG = ROOT / 'figures'
MYRAW = RAW / 'verify_precovery_sdss'

# ----------------------------------------------------------------------------------------------- unit vectors
def uvec(ra, dec):
    ra = np.radians(np.asarray(ra, float)); dec = np.radians(np.asarray(dec, float))
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def radec(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360, np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1)))


# ----------------------------------------------------------------------------------------------- Field table
def load_fields(path=MYRAW / 'fields_box_full.csv'):
    return np.genfromtxt(path, delimiter=',', names=True, skip_header=1)


# ----------------------------------------------------------------------------------------------- strands
class Strand:
    """Polyline on the sky, densified along great circles; signed perpendicular distances (+ = east)."""

    def __init__(self, ra, dec, hwhm=None, name='', step_deg=0.005):
        o = np.argsort(dec)
        ra, dec = np.asarray(ra, float)[o], np.asarray(dec, float)[o]
        self.name = name
        v = uvec(ra, dec)
        pts, dd, hw = [], [], []
        hwn = np.full(len(ra), np.nan) if hwhm is None else np.asarray(hwhm, float)[o]
        for i in range(len(v) - 1):
            ang = np.degrees(np.arccos(np.clip(np.dot(v[i], v[i + 1]), -1, 1)))
            n = max(int(np.ceil(ang / step_deg)), 1)
            t = np.linspace(0, 1, n, endpoint=False)[:, None]
            p = v[i] * (1 - t) + v[i + 1] * t
            pts.append(p / np.linalg.norm(p, axis=1, keepdims=True))
            dd.append(dec[i] + (dec[i + 1] - dec[i]) * t[:, 0])
            hw.append(hwn[i] + (hwn[i + 1] - hwn[i]) * t[:, 0])
        pts.append(v[-1:]); dd.append(dec[-1:]); hw.append(hwn[-1:])
        self.p = np.vstack(pts); self.pdec = np.concatenate(dd); self.phw = np.concatenate(hw)
        tan = np.gradient(self.p, axis=0)
        tan -= (tan * self.p).sum(1, keepdims=True) * self.p
        self.t = tan / np.linalg.norm(tan, axis=1, keepdims=True)      # direction of increasing Dec (~north)
        self.tree = cKDTree(self.p)
        self.dec_min, self.dec_max = dec.min(), dec.max()
        self.ra_nodes, self.dec_nodes = ra, dec

    def geometry(self, ra, dec):
        """Return signed perpendicular distance d (deg, + east), Dec of the nearest strand point,
        along-track overshoot beyond the strand ends (deg, 0 inside), and the local HWHM."""
        q = uvec(ra, dec)
        shp = q.shape[:-1]
        q = q.reshape(-1, 3)
        _, k = self.tree.query(q)
        p = self.p[k]; t = self.t[k]
        rel = q - p
        along = (rel * t).sum(1)
        e = np.cross(t, p)                       # with t ~ north and p outward, t x p points ~ east? check sign
        across = (rel * e).sum(1)
        # orientation check: east = increasing RA; ensure sign convention by testing a point due east
        d = np.degrees(np.arctan2(across, (q * p).sum(1)))
        over = np.zeros(len(q))
        lo = k == 0; hi = k == len(self.p) - 1
        over[lo] = np.clip(-np.degrees(along[lo]), 0, None)
        over[hi] = np.clip(np.degrees(along[hi]), 0, None)
        return (d * self._sign).reshape(shp), self.pdec[k].reshape(shp), over.reshape(shp), self.phw[k].reshape(shp)

    @property
    def _sign(self):
        if not hasattr(self, '_s'):
            self._s = 1.0
            i = len(self.p) // 2
            r0, d0 = radec(self.p[i])
            dd, *_ = self.geometry(np.array([r0 + 0.5]), np.array([d0]))
            self._s = 1.0 if dd[0] > 0 else -1.0
        return self._s


def load_strands():
    t = json.loads((DER / 'silk_trace.json').read_text())
    sp = t['spine']
    out = {'E': Strand([p['ra'] for p in sp], [p['dec'] for p in sp], [p['hwhm_deg'] for p in sp], 'E')}
    for k in ('W', 'W2'):
        pts = t['strands'][k]['points']
        out[k] = Strand([p['ra'] for p in pts], [p['dec'] for p in pts],
                        [p.get('hwhm_deg', t['half_width_deg']) for p in pts], k)
    return out, t


# ----------------------------------------------------------------------------------------------- ridge models
SIG_028 = 0.28 / 2.3548            # FWHM 0.28 deg


def desi_nodes():
    r = json.loads((DER / 'spectra_desi_results.json').read_text())['forward_model']
    nodes = {}
    for name, a, e in zip(r['basis'], r['amp_R'], r['err_R']):
        k, d = name.split('@')
        nodes.setdefault(k, []).append((float(d), a, e))
    return {k: np.array(sorted(v)) for k, v in nodes.items()}


def fibre_nodes():
    """SDSS sky-fibre strand amplitudes (R) per Dec bin (spectra_sdss_lamost_sdss_onoff.json), at bin centres."""
    d = json.loads((DER / 'spectra_sdss_lamost_sdss_onoff.json').read_text())['strand_fit']['OIII']['dec_bins']
    E, W = [], []
    for b in d:
        c = 0.5 * (b['dec_lo'] + b['dec_hi'])
        if np.isfinite(b['A_E']):
            E.append((c, b['A_E'], b['A_E_err']))
        if np.isfinite(b['A_W']):
            W.append((c, b['A_W'], b['A_W_err']))
    return {'E': np.array(E), 'W': np.array(W)}


class RidgeModel:
    """Sum of Gaussian ridges along strands.  amp: dict strand -> callable(dec_on_strand) or float.
    width: 'fixed' (FWHM fwhm) or 'trace' (trace HWHM clipped to 0.05-0.4 deg, as the DESI model).
    end_taper: Gaussian roll-off sigma (deg) beyond strand ends."""

    def __init__(self, strands, amp, width='fixed', fwhm=0.28, end_taper=0.15, ra_shift=0.0):
        self.strands = strands; self.amp = amp; self.width = width; self.fwhm = fwhm
        self.end_taper = end_taper; self.ra_shift = ra_shift

    def __call__(self, ra, dec):
        ra = np.asarray(ra, float) - self.ra_shift
        out = np.zeros(np.shape(ra))
        for k, a in self.amp.items():
            d, pd, over, hw = self.strands[k].geometry(ra, dec)
            if self.width == 'trace':
                sig = np.clip(hw, 0.05, 0.4) / 1.1774
            else:
                sig = self.fwhm / 2.3548
            amp = a(pd) if callable(a) else a
            out += amp * np.exp(-0.5 * (d / sig) ** 2) * np.exp(-0.5 * (over / self.end_taper) ** 2)
        return out


def make_model(kind, strands, ra_shift=0.0):
    """Named ridge models used in the verification."""
    if kind == 'gauss028':          # unit-amplitude E+W strands, FWHM 0.28 deg (model-free amplitudes)
        return RidgeModel(strands, {'E': 1.0, 'W': 1.0}, 'fixed', 0.28, ra_shift=ra_shift)
    if kind == 'gauss028_EWW2':
        return RidgeModel(strands, {'E': 1.0, 'W': 1.0, 'W2': 1.0}, 'fixed', 0.28, ra_shift=ra_shift)
    if kind == 'gaussE028':
        return RidgeModel(strands, {'E': 1.0}, 'fixed', 0.28, ra_shift=ra_shift)
    if kind.startswith('gaussfw'):  # e.g. gaussfw0.60
        fw = float(kind[7:])
        return RidgeModel(strands, {'E': 1.0, 'W': 1.0}, 'fixed', fw, ra_shift=ra_shift)
    if kind == 'desi':              # DESI forward-model amplitudes (R of [OIII]5007 peak), trace widths
        n = desi_nodes()
        amp = {k: (lambda dd, v=v: np.interp(dd, v[:, 0], v[:, 1])) for k, v in n.items()}
        return RidgeModel(strands, amp, 'trace', ra_shift=ra_shift)
    if kind == 'fibre':             # SDSS sky-fibre strand amplitudes (R), sigma 0.12 deg as in that fit
        n = fibre_nodes()
        amp = {k: (lambda dd, v=v: np.interp(dd, v[:, 0], v[:, 1])) for k, v in n.items()}
        return RidgeModel(strands, amp, 'fixed', 0.12 * 2.3548, ra_shift=ra_shift)
    raise ValueError(kind)


# ----------------------------------------------------------------------------------------------- field footprints
PIX = 0.396 / 3600.0
ACROSS_HALF = 1024 * PIX            # 13.5' / 2 (columns, perpendicular to scan)
ALONG_HALF = 744.5 * PIX            # 9.8' / 2 (1489 rows incl. 128-row overlap)


def scan_frames(a):
    """Per field: unit vectors (centre n, along-scan e_a, across-scan e_c) from neighbouring field centres."""
    n = uvec(a['ra'], a['dec'])
    ea = np.full_like(n, np.nan)
    key = a['run'].astype(int) * 10 + a['camcol'].astype(int)
    for k in np.unique(key):
        idx = np.where(key == k)[0]
        idx = idx[np.argsort(a['field'][idx])]
        if len(idx) < 2:
            continue
        p = n[idx]
        g = np.gradient(p, axis=0)
        ea[idx] = g
    ea -= (ea * n).sum(1, keepdims=True) * n
    ea /= np.linalg.norm(ea, axis=1, keepdims=True)
    ec = np.cross(n, ea)
    return n, ea, ec


def footprint_points(a, nx=9, ny=7, frames=None):
    """[N, nx*ny] RA/Dec samples inside each field footprint (13.5' across x 9.8' along scan)."""
    n, ea, ec = frames if frames is not None else scan_frames(a)
    xs = np.radians(np.linspace(-ACROSS_HALF, ACROSS_HALF, nx))
    ys = np.radians(np.linspace(-ALONG_HALF, ALONG_HALF, ny))
    X, Y = np.meshgrid(xs, ys)
    X = X.ravel(); Y = Y.ravel()
    P = n[:, None, :] + X[None, :, None] * ec[:, None, :] + Y[None, :, None] * ea[:, None, :]
    ok = np.isfinite(P).all(-1)
    P[~ok] = np.repeat(n[:, None, :], P.shape[1], 1)[~ok]
    return radec(P)


def field_template(model, a, fp=None, stat='mean'):
    ra, dec = fp if fp is not None else footprint_points(a)
    v = model(ra, dec)
    return v.mean(1) if stat == 'mean' else np.median(v, 1)


# ----------------------------------------------------------------------------------------------- photometry
NMGY_CGS = 3.631e-29            # erg s-1 cm-2 Hz-1
SR_PER_ARCSEC2 = (np.pi / 180 / 3600) ** 2
C_ANG = 2.99792458e18
H = 6.62607015e-27


def load_filter(b):
    for p in (RAW / f'filter_curves/sdss2010-{b}.ecsv', MYRAW / f'filter_curves/sdss2010-{b}.ecsv'):
        if p.exists():
            rows = [l.split() for l in p.read_text().splitlines() if l and not l.startswith('#')]
            arr = np.array([[float(x) for x in r] for r in rows[1:]])
            return arr[:, 0], arr[:, 1]
    raise FileNotFoundError(b)


def nmgy_per_R(b, lam_vac):
    """AB surface brightness (nanomaggies arcsec^-2) produced in SDSS band b by 1 R of a line at lam (vacuum A).
    Photon-counting response T(lam) (Doi+2010, airmass 1.3); response curves are in air wavelengths."""
    w, t = load_filter(b)
    lam_air = lam_vac / (1.0 + 2.735182e-4 + 131.4182 / lam_vac ** 2 + 2.76249e8 / lam_vac ** 4)
    T = np.interp(lam_air, w, t, left=0, right=0)
    norm = np.trapezoid(t / w, w)
    photons = 1e6 / (4 * np.pi)                     # photons s-1 cm-2 sr-1
    F = photons * H * C_ANG / lam_vac               # erg s-1 cm-2 sr-1
    fnu = F * T * lam_air / (C_ANG * norm)          # erg s-1 cm-2 Hz-1 sr-1
    return fnu * SR_PER_ARCSEC2 / NMGY_CGS


# filament line list: vacuum wavelengths and ENERGY ratios relative to [OIII]5007 (DESI differential, notes/spectra_desi.md;
# Halpha from Hbeta x 2.86; [NII],[SII] from the DESI positive-residual ratios; 4959 = 5007/2.98)
LINES = {
    'OIII5007': (5008.24, 1.0), 'OIII4959': (4960.30, 1 / 2.98), 'Hb': (4862.69, 1 / 20.0),
    'OII3727': (3728.48, 1.49), 'NeIII3869': (3869.86, 0.10), 'NeIII3968': (3968.59, 0.03),
    'Ha': (6564.61, 2.86 / 20.0), 'NII6583': (6585.27, 0.53 * 2.86 / 20.0), 'NII6548': (6549.86, 0.53 * 2.86 / 20.0 / 3.05),
    'SII6716': (6718.29, 0.40 * 2.86 / 20 * 0.57), 'SII6731': (6732.67, 0.40 * 2.86 / 20 * 0.43),
    'ArIII7136': (7137.77, 0.018), 'ArIII7751': (7753.2, 0.018 * 0.24), 'OII7325': (7325.0, 0.037),
    'SIII9069': (9071.1, 0.031), 'SIII9531': (9533.2, 0.031 * 2.47),
}


def band_response_per_R5007(b, lines=LINES, oii_ratio=None):
    """nmgy arcsec^-2 in band b per 1 R (photons) of [OIII]5007 for the filament spectrum (energy ratios -> photons)."""
    lam0 = lines['OIII5007'][0]
    tot = 0.0; parts = {}
    for k, (lam, eratio) in lines.items():
        if k == 'OII3727' and oii_ratio is not None:
            eratio = oii_ratio
        nphot = eratio * lam / lam0          # photon ratio
        v = nphot * nmgy_per_R(b, lam)
        parts[k] = v; tot += v
    return tot, parts


# ----------------------------------------------------------------------------------------------- fitting
def design_nuisance(f, c, drift_order=2, common=True):
    """Columns: one per field number (common mode, simultaneous in all camcols), camcol offsets (first dropped),
    per-camcol polynomial drifts in field number up to drift_order."""
    fu, fi = np.unique(f, return_inverse=True)
    cu, ci = np.unique(c, return_inverse=True)
    n = len(f)
    cols = []
    if common:
        S = np.zeros((n, len(fu))); S[np.arange(n), fi] = 1.0; cols.append(S)
    else:
        cols.append(np.ones((n, 1)))
    G = np.zeros((n, len(cu))); G[np.arange(n), ci] = 1.0
    cols.append(G[:, 1:])
    t = (f - np.median(f)) / 50.0
    for p in range(1, drift_order + 1):
        D = np.zeros((n, len(cu))); D[np.arange(n), ci] = t ** p
        if not common:
            cols.append(D)
        else:
            cols.append(D[:, 1:] if p >= 1 else D)
    return np.hstack(cols)


def smooth_time_basis(f, knot_fields=6):
    """Cubic B-spline basis in field number (for the 'no common-mode' model variant)."""
    from scipy.interpolate import BSpline
    lo, hi = f.min() - 0.5, f.max() + 0.5
    nk = max(int(np.ceil((hi - lo) / knot_fields)), 1)
    inner = np.linspace(lo, hi, nk + 1)
    k = 3
    t = np.r_[[lo] * k, inner, [hi] * k]
    nb = len(t) - k - 1
    B = np.zeros((len(f), nb))
    for j in range(nb):
        c = np.zeros(nb); c[j] = 1
        B[:, j] = BSpline(t, c, k, extrapolate=False)(f)
    B = np.nan_to_num(B)
    return B[:, B.sum(0) > 0]


def robust_fit(y, X, ntempl, clip=4.0, niter=10, w=None):
    """Least squares with iterative MAD clipping. The last `ntempl` columns of X are the templates.
    Returns dict(coef[ntempl], cov[ntempl,ntempl], rms, n, mask)."""
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    for _ in range(niter):
        Xo = X[ok]; yo = y[ok]
        keep = np.abs(Xo).sum(0) > 0
        coef, *_ = np.linalg.lstsq(Xo[:, keep], yo, rcond=None)
        res = y - X[:, keep] @ coef
        s = 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok])))
        new = ok & (np.abs(res) <= clip * max(s, 1e-9))
        new |= ~np.isfinite(y) & False
        if (new == ok).all():
            break
        ok = new
    Xo = X[ok][:, keep]; yo = y[ok]
    coef, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    res = yo - Xo @ coef
    dof = max(len(yo) - Xo.shape[1], 1)
    s2 = float(res @ res / dof)
    cov = s2 * np.linalg.pinv(Xo.T @ Xo)
    nt = ntempl
    full = np.zeros(X.shape[1]); full[keep] = coef
    return dict(coef=full[-nt:], cov=cov[-nt:, -nt:], rms=np.sqrt(s2), n=int(ok.sum()), dof=dof, mask=ok,
                resid=y - X @ full)
