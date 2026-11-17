"""verify_desi_common.py -- shared, independently written helpers for the wave-2 DESI verification (tag desi).

Only the on-disk FORMAT of the wave-1 windows is reused (the +-40 A pixel index WIN_IDX of the fixed DESI coadd
grid 3600 + 0.8 i A, from spectra_desi_lines.py); everything analytic below is re-implemented:
  * air->vacuum with the Edlen (1966) / IAU (Morton 1991) formula (wave 1 used Morton 2000 / Ciddor);
  * per-coadd tile-petal contributions from desi_dr1.fiberassign, weights = summed EFFTIME_SPEC of each tile;
  * per-coadd barycentric correction (astropy, KPNO, every exposure, EFFTIME-weighted);
  * per-coadd TRUE point-source fibre acceptance from the GFA (desi_dr1.exposure.fiber_fracflux_gfa),
    EFFTIME-weighted harmonic mean over exposures, with a Moffat(beta=3.5) + lambda^-0.2 chromatic model;
  * a geometry-only ridge template on the trace (nearest-sample distance on the sphere; works for any rotation)
    and the per-tile-petal mean subtraction that DESI's per-petal sky model imposes.
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import WIN_WAVE  # noqa: E402  (pixel grid of the saved windows = data format)

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
C = 299792.458
WAVE = np.asarray(WIN_WAVE, float)


# ----------------------------------------------------------------------------- wavelengths
def air2vac(lam):
    """Edlen 1966 / IAU (Morton 1991, ApJS 77, 119) standard air -> vacuum (A)."""
    lam = np.asarray(lam, float)
    vac = lam.copy()
    for _ in range(4):
        s2 = (1e4 / vac) ** 2
        n = 1 + 6.4328e-5 + 2.94981e-2 / (146 - s2) + 2.5540e-4 / (41 - s2)
        vac = lam * n
    return vac


_AIR = {  # air wavelengths (A): NIST ASD / Osterbrock & Ferland 2006
    'OII3726': 3726.032, 'OII3729': 3728.815, 'NeIII3869': 3868.760, 'Hd': 4101.734, 'Hg': 4340.464,
    'OIII4363': 4363.209, 'HeII4686': 4685.710, 'Hb': 4861.325, 'OIII4959': 4958.911, 'OIII5007': 5006.843,
'NII5755': 5754.59, 'HeI5876': 5875.62, 'OI6300': 6300.304, 'OI6364': 6363.776,
    'NII6548': 6548.04, 'Ha': 6562.80, 'NII6583': 6583.46, 'SII6716': 6716.44, 'SII6731': 6730.81,
    'ArIII7136': 7135.79, 'OII7320': 7319.99, 'OII7330': 7330.73, 'SIII9069': 9068.6,
    'FeII4359': 4359.34, 'FeII7155': 7155.16, 'CaII7291': 7291.47, 'CaII7324': 7323.89,
    'NI5198': 5197.90, 'NI5200b': 5200.26,
    # night-sky lines at Kitt Peak (topocentric): mercury street-lamp lines, Na D, [O I] airglow
    'HgI4047': 4046.563, 'HgI4358': 4358.328, 'HgI5461': 5460.735, 'HgI5770': 5769.598, 'HgI5791': 5790.663,
    'NaD1': 5895.924, 'NaD2': 5889.950, 'sky5577': 5577.339, 'sky6300': 6300.304, 'sky6364': 6363.776,
}
L = {k: float(air2vac(v)) for k, v in _AIR.items()}
SKY_LINES = [L[k] for k in ('HgI4047', 'HgI4358', 'HgI5461', 'HgI5770', 'HgI5791', 'NaD1', 'NaD2', 'sky5577')]
# generic continuum masks: nebular lines, except the weak hypothetical shock lines [Fe II]4359 and [Ca II]7324 that
# sit 1-4 A from Hg I 4358 / [O III]4363 / [O II]7320,7330 (masking them would blank those lines); they are only fitted
NEB_LINES = [v for k, v in L.items() if not (k.startswith('Hg') or k.startswith('Na') or k.startswith('sky'))
             and k not in ('FeII4359', 'CaII7324')]
GAL_REST = np.array([1215.67, 1240.14, 1397.6, 1549.06, 1640.4, 1908.73, 2799.12, 3346.8, 3426.9, 3727.1,
                     3729.9, 3869.86, 3934.78, 3969.59, 4102.89, 4305.6, 4341.69, 4862.68, 4960.30, 5008.24,
                     5176.7, 5895.6, 6302.05, 6549.86, 6564.61, 6585.27, 6718.29, 6732.67])


# ----------------------------------------------------------------------------- LSF
def lsf_sigma():
    """Median per-pixel LSF Gaussian sigma (A) on WAVE from the SPARCL wave_sigma arrays of the pilot download."""
    S = []
    for f in sorted((D / 'spectra_desi_windows').glob('pilot_*.npz')):
        with np.load(f, allow_pickle=False) as d:
            if 'wave_sigma' in d.files:
                S.append(d['wave_sigma'])
    S = np.concatenate(S).astype(float)
    S[S <= 0] = np.nan
    return np.nanmedian(S, 0), S.shape[0]


# ----------------------------------------------------------------------------- spectra
def load_spectra(include_ext=True):
    files = sorted((D / 'spectra_desi_windows').glob('*_[0-9][0-9][0-9][0-9][0-9].npz'))
    if include_ext:
        files += sorted((D / 'verify_desi_windows').glob('ext_[0-9][0-9][0-9][0-9][0-9].npz'))
    ids, F, I, M, src = [], [], [], [], []
    for f in files:
        with np.load(f, allow_pickle=False) as d:
            ids.append(d['sparcl_id'])
            F.append(d['flux'])
            I.append(d['ivar'])
            M.append(d['model'])
            src += ['ext' if f.parent.name == 'verify_desi_windows' else 'w1'] * len(d['sparcl_id'])
    ids = np.concatenate(ids)
    _, first = np.unique(ids, return_index=True)
    first = np.sort(first)
    return (ids[first], np.concatenate(F)[first], np.concatenate(I)[first], np.concatenate(M)[first],
            np.array(src)[first])


def load_meta(ids):
    cat = {}
    for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv')):
        cat[r['sparcl_id']] = r
    keep = np.array([i in cat for i in ids])
    rows = [cat[i] for i in ids[keep]]
    m = {k: np.array([r[k] for r in rows]) for k in ('targetid', 'survey', 'program', 'spectype', 'healpix')}
    for k in ('ra', 'dec', 'redshift', 'mean_psf_to_fiber_specflux', 'd_any', 'd_spine', 'dRA_spine', 'min_mjd',
              'max_mjd', 'coadd_exptime', 'deltachi2'):
        m[k] = np.array([float(r[k]) if r[k] not in ('', 'nan') else np.nan for r in rows])
    m['targetid'] = m['targetid'].astype(np.int64)
    m['healpix'] = m['healpix'].astype(int)
    return keep, m


# ----------------------------------------------------------------------------- exposures / tiles
def moffat_frac(fwhm, r=0.76, beta=3.5):
    alpha = fwhm / (2 * np.sqrt(2 ** (1 / beta) - 1))
    return 1 - (1 + (r / alpha) ** 2) ** (1 - beta)


def _fwhm_from_frac(f, r=0.76, beta=3.5):
    g = np.linspace(0.3, 4.0, 3701)
    fr = moffat_frac(g, r, beta)
    return np.interp(-np.asarray(f), -fr, g)


LAM_GFA = 6400.0   # effective wavelength of the GFA r-like band (approximate; +-300 A tested)


class Exposures:
    """Per-exposure GFA fibre fraction, seeing, EFFTIME, barycentric correction; per-tile aggregates."""

    def __init__(self):
        from astropy.coordinates import SkyCoord, EarthLocation
        from astropy.time import Time
        import astropy.units as u
        rows = list(csv.DictReader(open(ROOT / 'data/raw/verify_desi_tap/exposure_gfa_box.csv')))

        def g(r, k):
            try:
                v = float(r[k])
                return v if np.isfinite(v) else np.nan
            except (ValueError, KeyError):
                return np.nan
        self.rows = rows
        tile = np.array([int(r['tileid']) for r in rows])
        mjd = np.array([g(r, 'mjd') for r in rows])
        eff = np.array([g(r, 'efftime_spec') for r in rows])
        ff = np.array([g(r, 'fiber_fracflux_gfa') for r in rows])
        see = np.array([g(r, 'seeing_gfa') for r in rows])
        see_etc = np.array([g(r, 'seeing_etc') for r in rows])
        prog = np.array([r['program'] for r in rows])
        # fill missing GFA fibre fractions from seeing (Moffat model rescaled to the GFA values) or program median
        okf = np.isfinite(ff) & (ff > 0.05)
        see = np.where(see > 0.3, see, np.nan)          # 0.0 is used for 'not measured' in the table
        see_etc = np.where(see_etc > 0.3, see_etc, np.nan)
        s_any = np.where(np.isfinite(see), see, see_etc)
        okfit = okf & np.isfinite(s_any)
        k = np.median(ff[okfit] / moffat_frac(s_any[okfit])) if okfit.any() else 1.0
        self.ff_scale_vs_moffat = float(k)
        fill = ~okf
        ff_f = ff.copy()
        m1 = fill & np.isfinite(s_any)
        ff_f[m1] = k * moffat_frac(s_any[m1])
        for p in set(prog):
            m2 = fill & ~np.isfinite(s_any) & (prog == p)
            ff_f[m2] = np.median(ff[okf & (prog == p)]) if (okf & (prog == p)).any() else np.median(ff[okf])
        self.ff_filled = fill
        tra = {int(r['tileid']): None for r in rows}
        tiles = {int(r['tileid']): (r['survey'], r['program']) for r in
                 csv.DictReader(open(ROOT / 'data/raw/datalab_tap/tiles_box.csv'))}
        tilepos = {int(r['tileid']): (float(r['tilera']), float(r['tiledec'])) for r in
                   csv.DictReader(open(ROOT / 'data/raw/datalab_tap/tiles_box.csv'))}
        c = SkyCoord([tilepos[t][0] for t in tile] * u.deg, [tilepos[t][1] for t in tile] * u.deg)
        kpno = EarthLocation.from_geodetic(-111.59989 * u.deg, 31.96403 * u.deg, 2097 * u.m)
        vb = c.radial_velocity_correction('barycentric', obstime=Time(mjd, format='mjd'),
                                          location=kpno).to(u.km / u.s).value
        self.tile, self.mjd, self.eff, self.ff, self.see, self.prog, self.vb = tile, mjd, eff, ff_f, s_any, prog, vb
        self.tile_prog = tiles
        self.by_tile = defaultdict(list)
        for j, t in enumerate(tile):
            if np.isfinite(eff[j]) and eff[j] > 0:
                self.by_tile[t].append(j)

    def chromatic(self, j, lam, lam_ref=LAM_GFA, atm_frac=1.0):
        """Point-source fibre fraction of exposure j at wavelength lam (A): Moffat FWHM matched to the GFA value
        at lam_ref, the atmospheric share (atm_frac) of the FWHM scaling as lam^-0.2."""
        f0 = self.ff[j] / self.ff_scale_vs_moffat
        fw = _fwhm_from_frac(f0)
        fw_l = fw * (atm_frac * (lam / lam_ref) ** -0.2 + (1 - atm_frac))
        return self.ff_scale_vs_moffat * moffat_frac(fw_l)


class Coadds:
    """Tile-petal contributions, weights, f_true and v_corr for each coadd (targetid, survey, program)."""

    def __init__(self, meta, ex):
        bt = defaultdict(list)
        self.pet_pos = defaultdict(lambda: ([], []))
        for r in csv.DictReader(open(ROOT / 'data/raw/datalab_tap/fiberassign_box.csv')):
            t, p = int(r['tileid']), int(r['petal_loc'])
            bt[int(r['targetid'])].append((t, p))
            a, b = self.pet_pos[(t, p)]
            a.append(float(r['target_ra']))
            b.append(float(r['target_dec']))
        self.pet_keys = list(self.pet_pos.keys())
        self.pet_ra = {k: np.array(v[0]) for k, v in self.pet_pos.items()}
        self.pet_dec = {k: np.array(v[1]) for k, v in self.pet_pos.items()}
        N = len(meta['ra'])
        self.contrib = []
        self.ex = ex
        for i in range(N):
            out = []
            for t, p in bt.get(int(meta['targetid'][i]), []):
                if ex.tile_prog.get(t) == (meta['survey'][i], meta['program'][i]) and ex.by_tile.get(t):
                    w = float(np.sum(ex.eff[ex.by_tile[t]]))
                    out.append((t, p, w))
            s = sum(w for _, _, w in out)
            self.contrib.append([(t, p, w / s) for t, p, w in out] if s > 0 else [])
        self.ok = np.array([len(c) > 0 for c in self.contrib])
        self.n_tiles = np.array([len(c) for c in self.contrib])
        # group label = dominant tile-petal
        self.group = np.array([max(c, key=lambda x: x[2])[:2] if c else (-1, -1) for c in self.contrib])

    def exposures_of(self, i):
        js, ws = [], []
        for t, p, _ in self.contrib[i]:
            for j in self.ex.by_tile[t]:
                js.append(j)
                ws.append(self.ex.eff[j])
        return np.array(js, int), np.array(ws)

    def f_true(self, lam=5008.0, atm_frac=1.0, lam_ref=LAM_GFA):
        """EFFTIME-weighted harmonic-mean point-source acceptance at lam (coadd of F_e = SB*Omega/f_e)."""
        out = np.full(len(self.contrib), np.nan)
        cache = {}
        for i in range(len(self.contrib)):
            js, ws = self.exposures_of(i)
            if js.size == 0:
                continue
            fs = []
            for j in js:
                key = (j, lam, atm_frac, lam_ref)
                if key not in cache:
                    cache[key] = self.ex.chromatic(j, lam, lam_ref, atm_frac)
                fs.append(cache[key])
            fs = np.array(fs)
            out[i] = ws.sum() / np.sum(ws / fs)
        return out

    def vcorr(self):
        out = np.full(len(self.contrib), np.nan)
        spread = np.full(len(self.contrib), np.nan)
        for i in range(len(self.contrib)):
            js, ws = self.exposures_of(i)
            if js.size:
                out[i] = np.sum(ws * self.ex.vb[js]) / ws.sum()
                spread[i] = np.ptp(self.ex.vb[js])
        return out, spread

    def tiles_of(self, i):
        return sorted({t for t, _, _ in self.contrib[i]})


# ----------------------------------------------------------------------------- geometry template
def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]


def vec2radec(v):
    v = v / np.linalg.norm(v, axis=1)[:, None]
    return np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360, np.degrees(np.arcsin(v[:, 2]))


def rotmat(axis, ang_deg):
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    t = np.radians(ang_deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def trace_strands():
    t = json.loads((D / 'silk_trace.json').read_text())
    out = {'E': [(p['ra'], p['dec'], p['hwhm_deg']) for p in t['spine']]}
    for k in ('W', 'W2'):
        out[k] = [(p['ra'], p['dec'], p.get('hwhm_deg', t['half_width_deg'])) for p in t['strands'][k]['points']]
    return {k: np.array(v, float) for k, v in out.items()}


class Template:
    """Sum over strands of Gaussian cross-sections (HWHM from the trace, clipped 0.05-0.4 deg) about densely
    resampled polylines; R = optional 3x3 rotation applied to the trace (null tests)."""

    def __init__(self, R=None, strands=('E', 'W', 'W2'), amps=None, step=0.01):
        self.parts = []
        for s, arr in trace_strands().items():
            if s not in strands:
                continue
            ra, dec, hw = arr.T
            v = unitvec(ra, dec)
            seg = np.degrees(np.arccos(np.clip(np.sum(v[1:] * v[:-1], 1), -1, 1)))
            pts, hws = [], []
            for k in range(len(ra) - 1):
                n = max(int(np.ceil(seg[k] / step)), 1)
                fr = np.arange(n) / n
                pk = v[k][None, :] * (1 - fr[:, None]) + v[k + 1][None, :] * fr[:, None]
                pts.append(pk)
                hws.append(hw[k] * (1 - fr) + hw[k + 1] * fr)
            pts.append(v[-1:])
            hws.append(hw[-1:])
            P = np.concatenate(pts)
            P /= np.linalg.norm(P, axis=1)[:, None]
            if R is not None:
                P = P @ R.T
            self.parts.append((s, cKDTree(P), np.clip(np.concatenate(hws), 0.05, 0.4) / 1.17741,
                               1.0 if amps is None else amps.get(s, 1.0), P))

    def eval_parts(self, ra, dec):
        X = unitvec(ra, dec)
        out = {}
        for s, tree, sig, a, _ in self.parts:
            dist, idx = tree.query(X, distance_upper_bound=np.radians(2.5))
            ok = np.isfinite(dist)
            val = np.zeros(len(ra))
            d = np.degrees(dist[ok])
            val[ok] = a * np.exp(-0.5 * (d / sig[idx[ok]]) ** 2)
            out[s] = val
        return out

    def __call__(self, ra, dec):
        return sum(self.eval_parts(ra, dec).values())


def petal_subtracted(template, meta, co, parts=False):
    """Predicted residual Z_i = T(x_i) - sum_t w_t <T>_{tile t, petal p}; mean over all fibres of the tile-petal."""
    keys = co.pet_keys
    allra = np.concatenate([co.pet_ra[k] for k in keys])
    alldec = np.concatenate([co.pet_dec[k] for k in keys])
    lens = np.array([co.pet_ra[k].size for k in keys])
    ends = np.cumsum(lens)
    starts = ends - lens
    vals_all = template.eval_parts(allra, alldec)
    vals_fib = template.eval_parts(meta['ra'], meta['dec'])
    names = list(vals_all.keys()) if parts else ['all']
    if not parts:
        vals_all = {'all': sum(vals_all.values())}
        vals_fib = {'all': sum(vals_fib.values())}
    Z = {}
    for nm in names:
        cs = np.concatenate([[0.0], np.cumsum(vals_all[nm])])
        pm = {k: (cs[e] - cs[s]) / (e - s) for k, s, e in zip(keys, starts, ends)}
        sub = np.array([sum(w * pm[(t, p)] for t, p, w in c) if c else np.nan for c in co.contrib])
        Z[nm] = (vals_fib[nm] - sub, vals_fib[nm], sub)
    return Z if parts else Z['all'][0]


def good_fit(lf, g, col=0, scale=1.0, errmult=4.0, npix_frac=0.6):
    """Quality mask for a fitted group: finite, err > 0, err < errmult x median err, >= npix_frac of the
    window's maximum unmasked pixels, |amp/err| < 30 (gross fit failures)."""
    a, e, n = lf[g + '__amp'][:, col], lf[g + '__err'][:, col] * scale, lf[g + '__npix']
    m = np.isfinite(a) & np.isfinite(e) & (e > 0)
    m &= e < errmult * np.median(e[m])
    m &= n >= npix_frac * np.max(n)
    m &= np.abs(a / np.where(e > 0, e, np.inf)) < 30
    return m
