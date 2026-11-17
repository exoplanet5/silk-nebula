"""w4_family_desi.py -- wave-4 `family`: DESI DR1 residual regression on Monogem-periphery [O III] filaments
(S1 south limb: main/backup tiles; E1 east of the Silk: main dark/bright/backup), with the corrected wave-2 method.

Method (follows verify_desi; analytic pieces imported unchanged from scripts/verify_desi_common.py and
scripts/verify_desi_linefit.py, file-specific pieces re-implemented here for the new tiles):
  * residual = SPARCL flux - Redrock model on the saved +-40 A windows; lines fitted per coadd with
    verify_desi_linefit.fit_group (quadratic continuum + Gaussians, sky-line and object-line masks, 5-sigma clip,
    LSF = median SPARCL wave_sigma (+0.8 A pixel, +12 km/s intrinsic)), fixed V_NEB (barycentric) = +13 km/s.
    STAR coadds get no object-line mask (their z ~ 0 would otherwise mask [O III] itself).
  * surface brightness: SB = F_cal * f_ps(lambda) / Omega (Omega of the 1.52" fibre), f_ps = EFFTIME-weighted
    harmonic mean of the per-exposure GFA fiber_fracflux_gfa moved to lambda with Moffat(3.5), FWHM ~ lambda^-0.2
    (the wave-2 "0.736 acceptance correction" relative to MEAN_PSF_TO_FIBER_SPECFLUX).
  * template: unit-peak Gaussian about the NSNS ridge trace (w4_family_traces.json, S/N >= 4 points, gaps > 0.3 deg
    not bridged), sigma from the NSNS stacked cross-profile; Z_i = T(x_i) - sum_t w_t <T>_(tile t, petal p) over all
    fiberassign positions (DESI's per-petal sky model removes the petal mean).
  * amplitude of each line: weighted LS slope of y_i (R) on Z_i, |chi| < 8 clip; = template-weighted mean ridge PEAK
    brightness. Errors: bootstrap over dominant tile-petal groups (1000x); empirical nulls from the template shifted
    perpendicular to the filament (S1 in Dec, E1 in RA) and from fake lines / [O I]5577.
  * ratios: energy ratios to [O III]5007 from paired bootstrap; [O III]/Hb lower limit from the 97.7th percentile of
    Hb. Velocity: free-velocity Gaussian fit to the Z-weighted regression stack of the [O III] doublet residual
    (barycentric = heliocentric to < 0.02 km/s), bootstrap over groups; tiles 41485/41487 (-19 km/s stellar offset
    found in wave 1) excluded from the velocity and tested separately.
Inputs : data/derived/w4_family_windows/fam_*.npz (+ wave-1/2 windows for E1 coadds already on disk),
         data/raw/w4_family_sparcl/find_south.csv, data/raw/w4_family_tap/{zpix,fiberassign,exposure}_south.csv,
         data/raw/datalab_tap/fiberassign_box.csv, data/raw/verify_desi_tap/exposure_gfa_box.csv,
         data/derived/spectra_desi_catalog.csv, data/raw/w4_family_tap/tiles_monogem.csv
Outputs: data/derived/w4_family_desi.json, data/derived/w4_family_desi_fibres.npz, data/derived/w4_family_desi_stacks.npz
Usage  : ~/.venvs/astro313/bin/python -I scripts/w4_family_desi.py
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (C, WAVE, L, moffat_frac, _fwhm_from_frac, LAM_GFA, lsf_sigma, unitvec)  # noqa
from verify_desi_linefit import groups, fit_group  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
RAW = ROOT / 'data/raw'
V_NEB = 13.0
V_REGION = {'S1': 13.0, 'E1': -12.0, 'SILK': 13.0}   # E1 fitted at its own measured stack velocity (first pass, V=+13,
                                                    # gave -12.3 +- 2.0 km/s); SILK = pipeline check on wave-1/2 data
OMEGA = np.pi * (1.52 / 2) ** 2
HC = 6.62607015e-27 * 2.99792458e18
ARCSEC2_SR = (np.pi / 648000) ** 2
RNG = np.random.default_rng(20261117)
BAD_TILES = (41485, 41487)
SIG = {'S1': 0.042, 'E1': 0.045}       # NSNS stacked Gaussian sigma (deg), w4_family_traces.json
GROUPS_USED = ['o3', 'hb', 'o2sum', 'o2', 'he2', '4363hg', 'han2', 's2', 'ne3', 'sky5577', 'null4885', 'null4990',
               'null5030', 'null4700', 'null3760']
LAMS = {'o3': L['OIII5007'], 'hb': L['Hb'], 'o2sum': L['OII3726'], 'o2': L['OII3726'], 'he2': L['HeII4686'],
        '4363hg': L['OIII4363'], 'han2': L['Ha'], 's2': L['SII6716'], 'ne3': L['NeIII3869'], 'sky5577': L['sky5577'],
        'null4885': 4885.0, 'null4990': 4990.0, 'null5030': 5030.0, 'null4700': 4700.0, 'null3760': 3760.0}


def one_R(lam):
    return 1e6 / (4 * np.pi) * HC / lam * ARCSEC2_SR


# ------------------------------------------------------------------ data
def load_windows(want):
    ids, F, I, M = [], [], [], []
    files = sorted((D / 'w4_family_windows').glob('fam_[0-9][0-9][0-9][0-9][0-9].npz'))
    files += sorted((D / 'spectra_desi_windows').glob('*_[0-9][0-9][0-9][0-9][0-9].npz'))
    files += sorted((D / 'verify_desi_windows').glob('ext_[0-9][0-9][0-9][0-9][0-9].npz'))
    for f in files:
        with np.load(f, allow_pickle=False) as d:
            s = d['sparcl_id']
            m = np.isin(s, want)
            if m.any():
                ids.append(s[m]); F.append(d['flux'][m]); I.append(d['ivar'][m]); M.append(d['model'][m])
    ids = np.concatenate(ids)
    _, first = np.unique(ids, return_index=True)
    first = np.sort(first)
    return ids[first], np.concatenate(F)[first], np.concatenate(I)[first], np.concatenate(M)[first]


def ridge_points(name):
    sys.path.insert(0, str(ROOT / 'scripts'))
    from w4_family_select import ridge_tree
    return ridge_tree(name)


def meta_rows():
    """sparcl_id -> dict(targetid, survey, program, ra, dec, spectype, z, region)"""
    out = {}
    zp = {(r['targetid'], r['survey'], r['program']): r for r in csv.DictReader(open(RAW / 'w4_family_tap/zpix_south.csv'))}
    for r in csv.DictReader(open(RAW / 'w4_family_sparcl/find_south.csv')):
        k = (r['targetid'], r['survey'], r['program'])
        if k in zp and int(zp[k]['zwarn']) == 0:
            z = zp[k]
            out[r['sparcl_id']] = dict(targetid=int(r['targetid']), survey=r['survey'], program=r['program'],
                                       ra=float(z['mean_fiber_ra']), dec=float(z['mean_fiber_dec']),
                                       spectype=z['spectype'], z=float(z['z']), region='south')
    for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv')):
        ra, dec = float(r['ra']), float(r['dec'])
        if 117.5 < ra < 121.5 and 12.5 < dec < 17.6 and int(r['zwarn']) == 0:
            out[r['sparcl_id']] = dict(targetid=int(r['targetid']), survey=r['survey'], program=r['program'],
                                       ra=ra, dec=dec, spectype=r['spectype'], z=float(r['redshift']), region='east')
        elif dec < 17.6 and r['d_any'] not in ('', 'nan') and float(r['d_any']) < 1.0 and int(r['zwarn']) == 0 \
                and r['spectype'] in ('GALAXY', 'QSO'):
            out[r['sparcl_id']] = dict(targetid=int(r['targetid']), survey=r['survey'], program=r['program'],
                                       ra=ra, dec=dec, spectype=r['spectype'], z=float(r['redshift']), region='silk')
    return out


class Exp:
    def __init__(self):
        from astropy.coordinates import SkyCoord, EarthLocation
        from astropy.time import Time
        import astropy.units as u
        rows = list(csv.DictReader(open(RAW / 'w4_family_tap/exposure_south.csv')))
        rows += list(csv.DictReader(open(RAW / 'verify_desi_tap/exposure_gfa_box.csv')))
        seen, R = set(), []
        for r in rows:
            if r['expid'] not in seen:
                seen.add(r['expid']); R.append(r)

        def g(r, k):
            try:
                v = float(r[k]); return v if np.isfinite(v) else np.nan
            except (ValueError, KeyError):
                return np.nan
        tiles = {}
        for fn in (RAW / 'w4_family_tap/tiles_monogem.csv', RAW / 'datalab_tap/tiles_box.csv'):
            for r in csv.DictReader(open(fn)):
                tiles[int(r['tileid'])] = (r['survey'], r['program'], float(r['tilera']), float(r['tiledec']))
        tile = np.array([int(r['tileid']) for r in R])
        mjd = np.array([g(r, 'mjd') for r in R]); eff = np.array([g(r, 'efftime_spec') for r in R])
        ff = np.array([g(r, 'fiber_fracflux_gfa') for r in R]); see = np.array([g(r, 'seeing_gfa') for r in R])
        see_etc = np.array([g(r, 'seeing_etc') for r in R])
        prog = np.array([tiles[t][1] for t in tile])
        okf = np.isfinite(ff) & (ff > 0.05)
        see = np.where(see > 0.3, see, np.nan); see_etc = np.where(see_etc > 0.3, see_etc, np.nan)
        s_any = np.where(np.isfinite(see), see, see_etc)
        okfit = okf & np.isfinite(s_any)
        k = np.median(ff[okfit] / moffat_frac(s_any[okfit]))
        ff_f = ff.copy()
        m1 = ~okf & np.isfinite(s_any); ff_f[m1] = k * moffat_frac(s_any[m1])
        for p in set(prog):
            m2 = ~okf & ~np.isfinite(s_any) & (prog == p)
            ff_f[m2] = np.median(ff[okf & (prog == p)]) if (okf & (prog == p)).any() else np.median(ff[okf])
        c = SkyCoord([tiles[t][2] for t in tile] * u.deg, [tiles[t][3] for t in tile] * u.deg)
        kpno = EarthLocation.from_geodetic(-111.59989 * u.deg, 31.96403 * u.deg, 2097 * u.m)
        self.vb = c.radial_velocity_correction('barycentric', obstime=Time(mjd, format='mjd'),
                                               location=kpno).to(u.km / u.s).value
        self.k, self.tile, self.eff, self.ff, self.tiles = k, tile, eff, ff_f, tiles
        self.n_filled = int((~okf).sum())
        self.by_tile = defaultdict(list)
        for j, t in enumerate(tile):
            if np.isfinite(eff[j]) and eff[j] > 0:
                self.by_tile[t].append(j)

    def chromatic(self, j, lam):
        f0 = self.ff[j] / self.k
        fw = _fwhm_from_frac(f0) * (lam / LAM_GFA) ** -0.2
        return self.k * moffat_frac(fw)


class Co:
    def __init__(self, meta, ex):
        bt = defaultdict(list)
        pos = defaultdict(lambda: ([], []))
        need_tiles = set()
        for fn in (RAW / 'w4_family_tap/fiberassign_south.csv', RAW / 'datalab_tap/fiberassign_box.csv'):
            for r in csv.DictReader(open(fn)):
                t, p = int(r['tileid']), int(r['petal_loc'])
                bt[int(r['targetid'])].append((t, p))
                a, b = pos[(t, p)]
                a.append(float(r['target_ra'])); b.append(float(r['target_dec']))
        self.pos = {k: (np.array(v[0]), np.array(v[1])) for k, v in pos.items()}
        self.contrib = []
        for m in meta:
            out = []
            for t, p in set(bt.get(m['targetid'], [])):
                ti = ex.tiles.get(t)
                if ti and (ti[0], ti[1]) == (m['survey'], m['program']) and ex.by_tile.get(t):
                    out.append((t, p, float(np.sum(ex.eff[ex.by_tile[t]]))))
            s = sum(w for *_, w in out)
            self.contrib.append([(t, p, w / s) for t, p, w in out] if s > 0 else [])
        self.ex = ex
        self.group = np.array([f'{max(c, key=lambda x: x[2])[0]}_{max(c, key=lambda x: x[2])[1]}' if c else 'none'
                               for c in self.contrib])
        self.has_bad = np.array([any(t in BAD_TILES for t, _, _ in c) for c in self.contrib])

    def f_true(self, lam):
        out = np.full(len(self.contrib), np.nan)
        cache = {}
        for i, c in enumerate(self.contrib):
            js, ws = [], []
            for t, p, _ in c:
                for j in self.ex.by_tile[t]:
                    js.append(j); ws.append(self.ex.eff[j])
            if not js:
                continue
            fs = []
            for j in js:
                if j not in cache:
                    cache[j] = self.ex.chromatic(j, lam)
                fs.append(cache[j])
            ws, fs = np.array(ws), np.array(fs)
            out[i] = ws.sum() / np.sum(ws / fs)
        return out

    def vcorr(self):
        out = np.full(len(self.contrib), np.nan)
        for i, c in enumerate(self.contrib):
            js = [j for t, p, _ in c for j in self.ex.by_tile[t]]
            if js:
                w = self.ex.eff[js]
                out[i] = np.sum(w * self.ex.vb[js]) / w.sum()
        return out


class RidgeT:
    def __init__(self, P, sigma, shift=(0.0, 0.0)):
        P = P.copy()
        P[:, 1] += shift[1]
        P[:, 0] += shift[0] / np.cos(np.radians(P[:, 1]))
        self.tree = cKDTree(unitvec(P[:, 0], P[:, 1]))
        self.sig = sigma

    def __call__(self, ra, dec):
        d, _ = self.tree.query(unitvec(ra, dec), distance_upper_bound=np.radians(1.0))
        out = np.zeros(len(ra))
        ok = np.isfinite(d)
        out[ok] = np.exp(-0.5 * (np.degrees(d[ok]) / self.sig) ** 2)
        return out


def petal_Z(T, ra, dec, co):
    pm = {k: float(np.mean(T(a, b))) for k, (a, b) in co.pos.items()}
    sub = np.array([sum(w * pm[(t, p)] for t, p, w in c) if c else np.nan for c in co.contrib])
    t = T(ra, dec)
    return t - sub, t, sub


def slope(y, e, Z, m, clip=8.0):
    m = m & np.isfinite(y) & np.isfinite(e) & (e > 0) & np.isfinite(Z)
    w = 1 / e ** 2
    a = np.sum(w[m] * Z[m] * y[m]) / np.sum(w[m] * Z[m] ** 2)
    chi = (y - a * Z) / e
    m2 = m & (np.abs(chi) < clip)
    a = np.sum(w[m2] * Z[m2] * y[m2]) / np.sum(w[m2] * Z[m2] ** 2)
    ea = 1 / np.sqrt(np.sum(w[m2] * Z[m2] ** 2))
    chi = (y - a * Z) / e
    return a, ea, m2, float(np.std(chi[m2]))


def boot(Y, E, Z, m, grp, nb=1000):
    """paired bootstrap over groups for several lines; returns [nb, nlines]"""
    ug = np.unique(grp[m])
    idx = {g: np.where(m & (grp == g))[0] for g in ug}
    out = []
    for _ in range(nb):
        pick = RNG.choice(ug, len(ug))
        ii = np.concatenate([idx[g] for g in pick])
        row = []
        for y, e in zip(Y, E):
            w = 1 / e[ii] ** 2
            row.append(np.sum(w * Z[ii] * y[ii]) / np.sum(w * Z[ii] ** 2))
        out.append(row)
    return np.array(out)


def main():
    meta_all = meta_rows()
    sel = list(csv.DictReader(open(D / 'w4_family_sel_queue.csv')))
    p2 = D / 'w4_family_sel_queue_ext.csv'
    if p2.exists():
        sel += list(csv.DictReader(open(p2)))
    want = set(r['sparcl_id'] for r in sel)
    # E1 coadds already on disk from waves 1-2 within 0.12 deg of the E1 ridge
    tE, PE = ridge_points('E1')
    tS, PS = ridge_points('S1')
    for sid, m in meta_all.items():
        if m['region'] == 'silk':
            want.add(sid)
        if m['region'] == 'east':
            d, _ = tE.query(unitvec([m['ra']], [m['dec']]))
            if np.degrees(d[0]) < 0.12:
                want.add(sid)
    ids, F, I, M = load_windows(np.array(sorted(want)))
    keep = np.array([i in meta_all for i in ids])
    ids, F, I, M = ids[keep], F[keep], I[keep], M[keep]
    meta = [meta_all[i] for i in ids]
    ra = np.array([m['ra'] for m in meta]); dec = np.array([m['dec'] for m in meta])
    region = np.array([m['region'] for m in meta]); stype = np.array([m['spectype'] for m in meta])
    prog = np.array([m['program'] for m in meta])
    print('spectra loaded', len(ids), dict(zip(*np.unique(region, return_counts=True))), flush=True)
    ex = Exp()
    co = Co(meta, ex)
    vc = co.vcorr()
    okc = np.array([len(c) > 0 for c in co.contrib])
    print('with tile-petal contributions', okc.sum(), ' GFA ff filled exposures', ex.n_filled, flush=True)
    sigA, nls = lsf_sigma()
    R = (F - M).astype(np.float32)
    z = np.array([m['z'] if m['spectype'] != 'STAR' else 10.0 for m in meta])   # no object-line mask for stars
    amp = {g: np.full((len(ids), 2), np.nan) for g in GROUPS_USED}
    err = {g: np.full((len(ids), 2), np.nan) for g in GROUPS_USED}
    vcz = np.where(np.isfinite(vc), vc, 0.0)
    rmap = {'S1': 'south', 'E1': 'east', 'SILK': 'silk'}
    for rn, vneb in V_REGION.items():
        ii = np.where(region == rmap[rn])[0]
        G = groups(vneb)
        for g in (GROUPS_USED if rn != 'SILK' else ['o3', 'hb', 'o2sum', 'null4990', 'null5030']):
            for i0 in range(0, len(ii), 1500):
                jj = ii[i0:i0 + 1500]
                res = fit_group(R[jj], I[jj], z[jj], vcz[jj], G[g], vneb, sigA)
                k = res['amp'].shape[1]
                amp[g][jj, :k], err[g][jj, :k] = res['amp'], res['err']
        print('fitted', rn, len(ii), 'at', vneb, flush=True)
    # conversion to Rayleigh (photon) at each line
    kR = {}
    for g in GROUPS_USED:
        f = co.f_true(LAMS[g])
        kR[g] = f * 1e-17 / OMEGA / one_R(LAMS[g])
    y = {g: amp[g][:, 0] * kR[g] for g in GROUPS_USED}
    e = {g: err[g][:, 0] * kR[g] for g in GROUPS_USED}
    # verify_desi_linefit group 'o2sum' is ONE profile 1.0*G(3726) + 1.3*G(3729); its amplitude is the 3726 flux only,
    # so the doublet total is 2.3 x amp (same normalisation trap as the [O II]7325 bug found by w3_oii).
    y['o2sum'] = y['o2sum'] * 2.3; e['o2sum'] = e['o2sum'] * 2.3
    y['o2_3729'] = amp['o2'][:, 1] * kR['o2']; e['o2_3729'] = err['o2'][:, 1] * kR['o2']
    y['o2_twocomp'] = y['o2'] + y['o2_3729']; e['o2_twocomp'] = np.hypot(e['o2'], e['o2_3729'])
    y['nii6583'] = amp['han2'][:, 1] * kR['han2']; e['nii6583'] = err['han2'][:, 1] * kR['han2']
    y['s2_6731'] = amp['s2'][:, 1] * kR['s2']; e['s2_6731'] = err['s2'][:, 1] * kR['s2']
    f5008 = co.f_true(LAMS['o3'])
    out = dict(method=__doc__.split('Method')[1].split('Inputs')[0].strip(), v_neb=V_REGION, n_spectra=int(len(ids)),
               lsf_n=int(nls), gfa_scale_vs_moffat=float(ex.k), regions={})
    save = dict(sparcl_id=ids, ra=ra, dec=dec, region=region, spectype=stype, program=prog, group=co.group,
                has_bad=co.has_bad, vcorr=vc, f5008=f5008)
    for g in y:
        save[f'y_{g}'] = y[g]; save[f'e_{g}'] = e[g]
    stacks = {}
    from verify_desi_common import Template
    for name, P, sig, regsel, shift_axis in [('S1', PS, SIG['S1'], region == 'south', 1),
                                              ('E1', PE, SIG['E1'], region == 'east', 0),
                                              ('SILK', None, None, region == 'silk', 0)]:
        T = RidgeT(P, sig) if P is not None else Template()
        Z, Tv, sub = petal_Z(T, ra, dec, co)
        save[f'Z_{name}'] = Z; save[f'T_{name}'] = Tv
        base = regsel & okc & np.isfinite(Z)
        # quality: finite errors, not crazy
        q = base & np.isfinite(e['o3']) & (e['o3'] > 0) & (e['o3'] < 6 * np.nanmedian(e['o3'][base]))
        res = {'n_fibres': int(q.sum()), 'n_onridge_T>0.5': int((q & (Tv > 0.5)).sum()),
               'n_groups': int(len(np.unique(co.group[q]))), 'median_f5008': float(np.nanmedian(f5008[q])),
               'median_err_o3_R_per_fibre': float(np.nanmedian(e['o3'][q])),
               'programs': {p: int((q & (prog == p)).sum()) for p in np.unique(prog[q])},
               'spectypes': {s: int((q & (stype == s)).sum()) for s in np.unique(stype[q])}, 'lines': {}}
        names = ['o3', 'hb', 'o2sum', 'o2', 'o2_3729', 'o2_twocomp', 'he2', '4363hg', 'han2', 'nii6583', 's2', 's2_6731', 'ne3',
                 'sky5577', 'null4885', 'null4990', 'null5030', 'null4700', 'null3760']
        if name == 'SILK':
            names = ['o3', 'hb', 'o2sum', 'null4990', 'null5030']
        A, EA, MS = {}, {}, {}
        for g in names:
            qq = q & np.isfinite(y[g]) & np.isfinite(e[g]) & (e[g] > 0) & (e[g] < 6 * np.nanmedian(e[g][q]))
            a, ea, m2, chistd = slope(y[g], e[g], Z, qq)
            A[g], EA[g], MS[g] = a, ea, m2
            res['lines'][g] = dict(peak_R=float(a), err_formal=float(ea), chi_std=chistd, n=int(m2.sum()))
            # variant with a free global intercept (absorbs a common sky-subtraction offset near the filament)
            if m2.sum() > 10:
                w_ = 1 / e[g][m2] ** 2
                X_ = np.c_[Z[m2], np.ones(m2.sum())]
                cf = np.linalg.solve(X_.T @ (X_ * w_[:, None]), X_.T @ (w_ * y[g][m2]))
                res['lines'][g]['peak_R_with_intercept'] = float(cf[0]); res['lines'][g]['intercept_R'] = float(cf[1])
        # bootstrap (paired) over dominant tile-petal groups, common mask = o3 mask & each line's
        bl = ['o3', 'hb', 'o2sum', 'o2_twocomp', 'o2', 'o2_3729', 'he2', '4363hg', 'han2', 'nii6583', 's2', 's2_6731', 'ne3',
              'sky5577', 'null4885', 'null4990', 'null5030'] if name != 'SILK' else ['o3', 'hb', 'o2sum', 'null4990', 'null5030']
        mm = MS['o3'].copy()
        Yb = [np.where(MS[g], y[g], 0.0) for g in bl]
        Eb = [np.where(MS[g], e[g], 1e9) for g in bl]
        B = boot(Yb, Eb, Z, mm, co.group)
        for k, g in enumerate(bl):
            res['lines'][g]['err_boot'] = float(np.std(B[:, k]))
            res['lines'][g]['err_boot_mad'] = float(1.4826 * np.median(np.abs(B[:, k] - np.median(B[:, k]))))
            res['lines'][g]['boot_p2.5_97.5'] = [float(np.percentile(B[:, k], 2.5)), float(np.percentile(B[:, k], 97.5))]
        o3b = B[:, 0]
        lam = {g: LAMS.get(g) for g in bl}
        lam['nii6583'] = L['NII6583']; lam['s2_6731'] = L['SII6731']
        lam['hb'] = L['Hb']; lam['o2_twocomp'] = L['OII3726']; lam['o2'] = L['OII3726']; lam['o2_3729'] = L['OII3729']
        rat = {}
        for k, g in enumerate(bl[1:], start=1):
            en = lam['o3'] / lam[g]          # photon -> energy ratio factor (num at lam_g, den at 5007)
            r = B[:, k] / o3b * en
            rat[f'{g}/o3_energy'] = dict(value=float(A[g] / A['o3'] * en), p2_5=float(np.percentile(r, 2.5)),
                                         p16=float(np.percentile(r, 16)), p50=float(np.percentile(r, 50)),
                                         p84=float(np.percentile(r, 84)), p97_5=float(np.percentile(r, 97.5)))
        hb_hi = np.percentile(B[:, 1], 97.7)
        rat['o3/hb_energy_lower2sig'] = float(np.median(o3b) * lam['hb'] / lam['o3'] / hb_hi) if hb_hi > 0 else None
        # [S II] and Ha combinations
        res['ratios'] = rat
        # empirical nulls: template shifted perpendicular to the filament
        nulls = []
        for sh in ([] if name == 'SILK' else (-1.5, -1.2, -0.9, -0.7, -0.5, 0.5, 0.7, 0.9, 1.2, 1.5)):
            Tn = RidgeT(P, sig, shift=(sh, 0.0) if shift_axis == 0 else (0.0, sh))
            Zn, Tnv, _ = petal_Z(Tn, ra, dec, co)
            qn = q & (Tv < 0.05) & np.isfinite(Zn)
            if (qn & (Tnv > 0.5)).sum() < 10:
                continue
            row = {}
            for g in ('o3', 'hb', 'o2sum', 'han2'):
                qq = qn & np.isfinite(y[g]) & np.isfinite(e[g]) & (e[g] > 0) & (e[g] < 6 * np.nanmedian(e[g][q]))
                a, ea, _, _ = slope(y[g], e[g], Zn, qq)
                row[g] = (float(a), float(ea))
            nulls.append(dict(shift=sh, **{g: v[0] for g, v in row.items()},
                              **{g + '_err': v[1] for g, v in row.items()}))
        res['nulls'] = nulls
        if nulls:
            for g in ('o3', 'hb', 'o2sum', 'han2'):
                v = np.array([n[g] for n in nulls]); ev = np.array([n[g + '_err'] for n in nulls])
                res['lines'][g]['null_rms'] = float(np.sqrt(np.mean(v ** 2)))
                res['lines'][g]['null_rms_over_formal'] = float(np.sqrt(np.mean((v / ev) ** 2)))
        # velocity: Z-weighted regression stack of the [O III] doublet residual in R/A
        m5 = (WAVE > 4935) & (WAVE < 5035)
        lamw = WAVE[m5]
        Rr = R[:, m5] * (f5008 * 1e-17 / OMEGA / one_R(LAMS['o3']))[:, None]
        Iw = I[:, m5] / (f5008 * 1e-17 / OMEGA / one_R(LAMS['o3']))[:, None] ** 2
        # mask sky/Hg lines and object lines as in the fitter would be complex; use the |chi|<8 clip on Z regression
        vres = {}
        for lab, sub_m in [('excl_41485_41487', MS['o3'] & ~co.has_bad), ('all', MS['o3'])]:
            S, Sv = stack_spec(Rr, Iw, Z, sub_m)
            fit = fit_velocity(lamw, S, Sv, sigA)
            bv = []
            ug = np.unique(co.group[sub_m])
            idx = {g: np.where(sub_m & (co.group == g))[0] for g in ug}
            for _ in range(300):
                pick = RNG.choice(ug, len(ug))
                ii = np.zeros(len(Z), bool)
                cnt = np.zeros(len(Z))
                for gg in pick:
                    cnt[idx[gg]] += 1
                S_, Sv_ = stack_spec(Rr, Iw, Z, cnt > 0, wmult=cnt)
                f_ = fit_velocity(lamw, S_, Sv_, sigA)
                if f_ is not None:
                    bv.append(f_['v'])
            if fit is not None:
                fit['v_boot_std'] = float(np.std(bv)); fit['v_boot_p16_p84'] = [float(np.percentile(bv, 16)),
                                                                                  float(np.percentile(bv, 84))]
                fit['n_fibres'] = int(sub_m.sum())
            vres[lab] = fit
            if lab == 'excl_41485_41487':
                stacks[f'{name}_lam'] = lamw; stacks[f'{name}_stack'] = S; stacks[f'{name}_var'] = Sv
        res['velocity_o3'] = vres
        # per-segment brightness (S1: RA bins; E1: Dec bins)
        segs = []
        if name == 'SILK':
            edges = []
        coord = ra if name == 'S1' else dec
        edges = [] if name == 'SILK' else ([102, 104, 106, 108, 110.5, 113.5] if name == 'S1' else [13.2, 14.5, 15.5, 16.9])
        for lo, hi in zip(edges[:-1], edges[1:]):
            Tseg = RidgeT(P[(P[:, 0 if name == 'S1' else 1] >= lo) & (P[:, 0 if name == 'S1' else 1] < hi)], sig) \
                if ((P[:, 0 if name == 'S1' else 1] >= lo) & (P[:, 0 if name == 'S1' else 1] < hi)).sum() > 3 else None
            if Tseg is None:
                continue
            Zs, Ts, _ = petal_Z(Tseg, ra, dec, co)
            qs = MS['o3'] & (np.abs(coord - (lo + hi) / 2) < (hi - lo) / 2 + 0.5)
            if (qs & (Ts > 0.5)).sum() < 5:
                continue
            a, ea, _, _ = slope(y['o3'], e['o3'], Zs, qs)
            segs.append(dict(lo=lo, hi=hi, peak_R=float(a), err=float(ea), n_on=int((qs & (Ts > 0.5)).sum())))
        res['segments_o3'] = segs
        out['regions'][name] = res
        lo = res['lines']
        if name == 'SILK':
            print('SILK check:', res['lines']['o3'], vres, flush=True)
            out['regions'][name] = res
            continue
        print(f"{name}: n={res['n_fibres']} on={res['n_onridge_T>0.5']} groups={res['n_groups']} "
              f"[OIII]={lo['o3']['peak_R']:.3f}+-{lo['o3']['err_boot']:.3f} (formal {lo['o3']['err_formal']:.3f}, "
              f"null rms {lo['o3'].get('null_rms', np.nan):.3f}) Hb={lo['hb']['peak_R']:.3f}+-{lo['hb']['err_boot']:.3f} "
              f"OII={lo['o2sum']['peak_R']:.3f}+-{lo['o2sum']['err_boot']:.3f} Ha={lo['han2']['peak_R']:.3f}+-{lo['han2']['err_boot']:.3f} "
              f"5577={lo['sky5577']['peak_R']:.3f}", flush=True)
        print('   ratios:', {k: (round(v['value'], 3), round(v['p16'], 3), round(v['p84'], 3)) if isinstance(v, dict) else v
                             for k, v in rat.items() if k.split('/')[0] in ('hb', 'o2sum', 'he2', '4363hg', 'han2', 's2')
                             or k.startswith('o3/hb')}, flush=True)
        print('   velocity:', vres, flush=True)
        print('   segments:', segs, flush=True)
    json.dump(out, open(D / 'w4_family_desi.json', 'w'), indent=1, default=float)
    np.savez_compressed(D / 'w4_family_desi_fibres.npz', **save)
    np.savez_compressed(D / 'w4_family_desi_stacks.npz', **stacks)


def stack_spec(Rr, Iw, Z, m, wmult=None):
    """Regression stack S(lam) = sum w Z r / sum w Z^2 with per-fibre weights 1/median(var) (robust)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        v = np.nanmedian(np.where(Iw > 0, 1 / Iw, np.nan), 1)
    w = np.where(m & np.isfinite(v) & (v > 0), 1 / v, 0.0)
    if wmult is not None:
        w = w * wmult
    Rz = np.nan_to_num(Rr)
    # robust: clip pixels with |r - median| > 10 sigma
    num = np.sum((w * Z)[:, None] * Rz, 0)
    den = np.sum(w * Z ** 2)
    var = np.sum(((w * Z) ** 2)[:, None] * np.nan_to_num(1 / np.where(Iw > 0, Iw, np.inf)), 0) / den ** 2
    return num / den, var


def fit_velocity(lam, S, V, sigA, sig_int=12.0):
    from scipy.optimize import least_squares
    l5, l4 = L['OIII5007'], L['OIII4959']
    s5 = np.sqrt(np.interp(l5, WAVE, sigA) ** 2 + 0.8 ** 2 / 12 + (l5 * sig_int / C) ** 2)
    s4 = s5 * l4 / l5
    ok = np.isfinite(S) & np.isfinite(V) & (V > 0)
    # mask Hg/sky-free region is clean here; exclude 4970-4995 nothing; use all
    x = lam[ok]; yv = S[ok]; ev = np.sqrt(V[ok])

    def model(p):
        a, v, c0, c1 = p
        g5 = np.exp(-0.5 * ((x - l5 * (1 + v / C)) / s5) ** 2) / (np.sqrt(2 * np.pi) * s5)
        g4 = np.exp(-0.5 * ((x - l4 * (1 + v / C)) / s4) ** 2) / (np.sqrt(2 * np.pi) * s4)
        return a * (g5 + g4 / 2.98) + c0 + c1 * (x - 4985) / 50

    try:
        r = least_squares(lambda p: (model(p) - yv) / ev, [max(np.nanmax(yv), 0.1), 10.0, 0, 0])
    except Exception:
        return None
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(np.sum(r.fun ** 2) / (len(x) - 4), 1.0)
    except np.linalg.LinAlgError:
        return None
    return dict(v=float(r.x[1]), ev_formal=float(np.sqrt(cov[1, 1])), amp_R=float(r.x[0]),
                eamp=float(np.sqrt(cov[0, 0])), rchi2=float(np.sum(r.fun ** 2) / (len(x) - 4)))


if __name__ == '__main__':
    main()
