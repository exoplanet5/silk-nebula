"""Small image-registration helpers (star blobs -> Tycho-2 matching) for the geometry work."""
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = '/Users/mickey/astronomy/silk-nebula'


def load_tycho(vtmax=11.0):
    p = f'{ROOT}/data/raw/vizier_tycho2/tyc2_VT11_box.tsv'
    rows = []
    for line in open(p):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        try:
            ra, de = float(f[0]), float(f[1])
        except ValueError:
            continue
        vt = float(f[3]) if f[3].strip() else (float(f[2]) if f[2].strip() else 99)
        rows.append((ra, de, vt))
    a = np.array(rows)
    return a[a[:, 2] < vtmax]


def blobs(mask, min_area=3, weight=None):
    lab, n = ndimage.label(mask)
    idx = np.arange(1, n + 1)
    area = ndimage.sum(np.ones_like(lab), lab, idx)
    w = weight if weight is not None else np.ones_like(lab, dtype=float)
    cy, cx = np.array(ndimage.center_of_mass(w, lab, idx)).T
    flux = ndimage.sum(w, lab, idx)
    m = area >= min_area
    return np.c_[cx[m], cy[m], area[m], flux[m]]


def gnomonic(ra, de, ra0, de0):
    ra, de, ra0, de0 = map(np.radians, (ra, de, ra0, de0))
    c = np.sin(de0) * np.sin(de) + np.cos(de0) * np.cos(de) * np.cos(ra - ra0)
    x = np.cos(de) * np.sin(ra - ra0) / c
    y = (np.cos(de0) * np.sin(de) - np.sin(de0) * np.cos(de) * np.cos(ra - ra0)) / c
    return np.degrees(x), np.degrees(y)


def inv_gnomonic(x, y, ra0, de0):
    x, y = np.radians(x), np.radians(y)
    ra0, de0 = np.radians(ra0), np.radians(de0)
    rho = np.hypot(x, y)
    c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        de = np.arcsin(np.cos(c) * np.sin(de0) + np.where(rho > 0, y * np.sin(c) * np.cos(de0) / rho, 0))
        ra = ra0 + np.arctan2(x * np.sin(c), rho * np.cos(de0) * np.cos(c) - y * np.sin(de0) * np.sin(c))
    return np.degrees(ra) % 360, np.degrees(de)


def poly_terms(x, y, order):
    cols = []
    for i in range(order + 1):
        for j in range(order + 1 - i):
            cols.append(x ** i * y ** j)
    return np.vstack(cols).T


class PolyMap:
    """pixel (x,y) -> (xi, eta) gnomonic plane about (ra0,de0), polynomial of given order."""

    def __init__(self, ra0, de0, order=1):
        self.ra0, self.de0, self.order = ra0, de0, order
        self.cx = self.cy = None
        self.x0 = self.y0 = 0.0
        self.s = 1.0

    def fit(self, px, py, ra, de, w=None):
        xi, eta = gnomonic(ra, de, self.ra0, self.de0)
        self.x0, self.y0, self.s = np.mean(px), np.mean(py), max(np.std(px), np.std(py))
        A = poly_terms((px - self.x0) / self.s, (py - self.y0) / self.s, self.order)
        W = np.ones(len(px)) if w is None else w
        self.cx, *_ = np.linalg.lstsq(A * W[:, None], xi * W, rcond=None)
        self.cy, *_ = np.linalg.lstsq(A * W[:, None], eta * W, rcond=None)
        return self

    def __call__(self, px, py):
        A = poly_terms((np.asarray(px, float) - self.x0) / self.s, (np.asarray(py, float) - self.y0) / self.s, self.order)
        return inv_gnomonic(A @ self.cx, A @ self.cy, self.ra0, self.de0)

    def inverse(self, ra, de, guess_xy=None, n_iter=30):
        """numerical inverse via Newton iterations."""
        ra = np.atleast_1d(ra).astype(float); de = np.atleast_1d(de).astype(float)
        xi_t, eta_t = gnomonic(ra, de, self.ra0, self.de0)
        if guess_xy is None:
            px = np.full_like(ra, self.x0); py = np.full_like(ra, self.y0)
        else:
            px, py = [np.array(g, float) for g in guess_xy]
        h = 0.5
        for _ in range(n_iter):
            def f(x, y):
                A = poly_terms((x - self.x0) / self.s, (y - self.y0) / self.s, self.order)
                return A @ self.cx, A @ self.cy
            fx, fy = f(px, py)
            fxa, fya = f(px + h, py); fxb, fyb = f(px, py + h)
            j11, j21 = (fxa - fx) / h, (fya - fy) / h
            j12, j22 = (fxb - fx) / h, (fyb - fy) / h
            det = j11 * j22 - j12 * j21
            dx, dy = xi_t - fx, eta_t - fy
            px = px + (j22 * dx - j12 * dy) / det
            py = py + (-j21 * dx + j11 * dy) / det
        return px, py


def match_fit(det_xy, cat, mapping, tols=(0.3, 0.15, 0.08, 0.05, 0.03), order=1, ra0=None, de0=None,
              ndet=None, ncat=None):
    """Iteratively match detections (sorted bright-first) to catalog (ra,dec,mag) and fit PolyMap."""
    d = det_xy if ndet is None else det_xy[:ndet]
    c = cat[np.argsort(cat[:, 2])]
    c = c if ncat is None else c[:ncat]
    ra0 = np.median(c[:, 0]) if ra0 is None else ra0
    de0 = np.median(c[:, 1]) if de0 is None else de0
    tree = cKDTree(np.c_[c[:, 0] * np.cos(np.radians(de0)), c[:, 1]])
    pm = None
    for k, tol in enumerate(tols):
        ra, de = (mapping if pm is None else pm)(d[:, 0], d[:, 1])
        dist, idx = tree.query(np.c_[ra * np.cos(np.radians(de0)), de])
        m = dist < tol
        # unique matches only
        u, cnt = np.unique(idx[m], return_counts=True)
        dup = set(u[cnt > 1])
        m &= np.array([i not in dup for i in idx])
        o = order if (k >= 2 and m.sum() > 6 * (order + 1) * (order + 2) / 2) else 1
        pm = PolyMap(ra0, de0, o).fit(d[m, 0], d[m, 1], c[idx[m], 0], c[idx[m], 1])
    ra, de = pm(d[:, 0], d[:, 1])
    dist, idx = tree.query(np.c_[ra * np.cos(np.radians(de0)), de])
    m = dist < tols[-1]
    res = dist[m] * 3600
    return pm, dict(n_match=int(m.sum()), n_det=len(d), n_cat=len(c), rms_arcsec=float(np.sqrt(np.mean(res ** 2))),
                    median_arcsec=float(np.median(res)), order=pm.order), (d[m], c[idx[m]])
