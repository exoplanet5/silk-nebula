"""spectra_desi_analyze.py -- maps, forward-model fit, regression spectra, ratios and velocities.

Inputs : data/derived/spectra_desi_fits.npz (spectra_desi_measure.py), spectra_desi_catalog.csv,
         windows npz files (for regression/stacked spectra), petal forward model (spectra_desi_petalmodel.py).
Outputs: data/derived/spectra_desi_results.json, spectra_desi_fiber_table.csv, spectra_desi_regspec.npz
Usage  : python -I scripts/spectra_desi_analyze.py
"""
import csv
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_desi_lines import LINES, WIN_WAVE, C_KMS  # noqa: E402
import spectra_desi_measure as meas  # noqa: E402
from spectra_desi_petalmodel import PetalModel, RidgeBasis  # noqa: E402
from silk_geometry import REGIONS  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'data/derived'
FIB_D_ARCSEC = 1.52                      # mean DESI fibre diameter (desispec skymag.py, DESI-6043)
FIB_AREA = np.pi * (FIB_D_ARCSEC / 2) ** 2
SR_PER_ARCSEC2 = (np.pi / 648000) ** 2
HC = 1.98644586e-8                       # erg A


def rayleigh_per_unit(lam):
    """SB in Rayleigh produced by 1e-17 erg/s/cm2 of in-fibre flux at wavelength lam (A)."""
    one_R = 1e6 / (4 * np.pi) * (HC / lam) * SR_PER_ARCSEC2      # erg/s/cm2/arcsec2 per R
    return 1e-17 / FIB_AREA / one_R


def load():
    f = np.load(D / 'spectra_desi_fits.npz', allow_pickle=False)
    res = {k: f[k] for k in f.files}
    cat = {r['sparcl_id']: r for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv'))}
    ids = res['sparcl_id']
    meta = {k: np.array([cat[i][k] for i in ids]) for k in
            ('targetid', 'survey', 'program', 'spectype', 'ra', 'dec', 'redshift', 'zwarn', 'dRA_spine',
             'd_spine', 'd_any', 'mean_psf_to_fiber_specflux', 'coadd_exptime', 'mean_mjd', 'coadd_numtile',
             'glon', 'glat')}
    for k in ('ra', 'dec', 'redshift', 'dRA_spine', 'd_spine', 'd_any', 'mean_psf_to_fiber_specflux',
              'coadd_exptime', 'mean_mjd', 'glon', 'glat'):
        meta[k] = meta[k].astype(float)
    meta['coadd_numtile'] = meta['coadd_numtile'].astype(int)
    return res, meta


def good(res, g, j=0, mode='fix'):
    a, e = res[f'{g}__amp_{mode}'][:, j], res[f'{g}__err_{mode}'][:, j]
    return np.isfinite(a) & np.isfinite(e) & (e > 0) & ~res[f'{g}__flag_own'] & (res[f'{g}__npix_good'] > 40)


def wls(y, e, X):
    """Weighted least squares y = X b; returns b, cov, chi2, dof."""
    w = 1 / e ** 2
    A = X.T @ (X * w[:, None])
    b = np.linalg.solve(A, X.T @ (w * y))
    cov = np.linalg.inv(A)
    chi2 = float(np.sum(w * (y - X @ b) ** 2))
    return b, cov, chi2, len(y) - X.shape[1]


def neighbour_map(ra, dec, val, err, ok, radius=0.12):
    """Inverse-variance mean of val over OTHER fibres within radius (deg); independent-noise instrument."""
    from scipy.spatial import cKDTree
    x = np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y = np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z = np.sin(np.radians(dec))
    P = np.c_[x, y, z]
    tree = cKDTree(P[ok])
    idx_ok = np.where(ok)[0]
    r = 2 * np.sin(np.radians(radius) / 2)
    out = np.full(len(ra), np.nan)
    nn = np.zeros(len(ra), int)
    w = 1 / err[ok] ** 2
    v = val[ok]
    for i, nb in enumerate(tree.query_ball_point(P, r)):
        nb = [j for j in nb if idx_ok[j] != i]
        if len(nb) >= 3:
            out[i] = np.sum(w[nb] * v[nb]) / np.sum(w[nb])
            nn[i] = len(nb)
    return out, nn


def fibre_vcorr(pm, meta, info):
    """EFFTIME-weighted barycentric correction (km/s) of the exposures contributing to each coadd."""
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy.time import Time
    import astropy.units as u
    kpno = EarthLocation.from_geodetic(-111.5997 * u.deg, 31.9634 * u.deg, 2120 * u.m)
    ex = list(csv.DictReader(open(ROOT / 'data/raw/datalab_tap/exposure_box.csv')))
    c = SkyCoord([float(r['tilera']) for r in ex] * u.deg, [float(r['tiledec']) for r in ex] * u.deg)
    vb = c.radial_velocity_correction('barycentric', obstime=Time([float(r['mjd']) for r in ex], format='mjd'),
                                      location=kpno).to(u.km / u.s).value
    per_tile = {}
    for r, v in zip(ex, vb):
        try:
            w = float(r['efftime_spec'])
        except ValueError:
            continue
        per_tile.setdefault(int(r['tileid']), []).append((v, w))
    tile_v = {t: np.sum([v * w for v, w in L]) / max(np.sum([w for _, w in L]), 1e-9) for t, L in per_tile.items()}
    out = np.full(len(info), np.nan)
    for i, c_ in enumerate(info):
        if c_:
            vv = [(tile_v[t], w) for t, p, w in c_ if t in tile_v]
            if vv:
                out[i] = np.sum([a * b for a, b in vv]) / np.sum([b for _, b in vv])
    return out


def main():
    res, meta = load()
    N = len(res['sparcl_id'])
    out = {'n_spectra': N}
    ffrac = np.clip(meta['mean_psf_to_fiber_specflux'], 0.3, 1.0)
    # ---------------- [OIII] per-fibre table in Rayleigh ----------------
    g = 'OIII5007'
    okO = good(res, g)
    kR = rayleigh_per_unit(LINES['OIII5007'])
    O, eO = res[f'{g}__amp_fix'][:, 0], res[f'{g}__err_fix'][:, 0]
    OR, eOR = O * ffrac * kR, eO * ffrac * kR
    out['rayleigh_per_1e-17_infibre_5008'] = kR
    out['median_ffrac'] = float(np.median(ffrac))
    # null tests
    nulls = {}
    for nn in ('null4990', 'null5030', 'null6620', 'null4880'):
        m = good(res, nn)
        chi = res[f'{nn}__amp_fix'][m, 0] / res[f'{nn}__err_fix'][m, 0]
        nulls[nn] = dict(n=int(m.sum()), mean_chi=float(np.mean(chi)), std_chi=float(np.std(chi)),
                         mad_std_chi=float(1.4826 * np.median(np.abs(chi - np.median(chi)))),
                         frac_gt5=float(np.mean(np.abs(chi) > 5)))
    out['null_tests'] = nulls
    errscale = float(np.median([v['mad_std_chi'] for v in nulls.values()]))
    out['error_scale_from_nulls'] = errscale
    chiO = O[okO] / eO[okO]
    out['OIII_chi_stats'] = dict(n=int(okO.sum()), mean=float(chiO.mean()), std=float(chiO.std()),
                                 frac_gt5=float(np.mean(chiO > 5)), frac_lt_m5=float(np.mean(chiO < -5)),
                                 median_err_1e17=float(np.median(eO[okO])),
                                 median_err_R=float(np.median(eOR[okO])))
    # ---------------- forward model ----------------
    pm = PetalModel()
    basis = RidgeBasis(strands=('E', 'W', 'W2'), seg_deg=1.0)
    Dm, X, S, okc, info = pm.design(meta['ra'], meta['dec'], meta['targetid'], meta['survey'], meta['program'],
                                    basis)
    use = okO & okc
    active = np.abs(Dm[use]).sum(0) > 1.0          # basis functions constrained by the data
    names = [n for n, a in zip(basis.names, active) if a]
    Xd = Dm[use][:, active]
    y, e = OR[use], eOR[use] * errscale
    # robust: iterate WLS with 4-sigma clipping
    keep = np.ones(len(y), bool)
    for it in range(4):
        b, cov, chi2, dof = wls(y[keep], e[keep], Xd[keep])
        r = (y - Xd @ b) / e
        keep = np.abs(r) < 4
    eb = np.sqrt(np.diag(cov) * max(chi2 / dof, 1))
    chi2_null = float(np.sum((y[keep] / e[keep]) ** 2))
    det = dict(chi2_null=chi2_null, chi2_model=chi2, delta_chi2=chi2_null - chi2, n_par=int(Xd.shape[1]),
               n=int(keep.sum()))
    # null test: identical model with the whole trace displaced in RA (petal means recomputed)
    shifts = {}
    for sh in (-4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0):
        bs = RidgeBasis(strands=('E', 'W', 'W2'), seg_deg=1.0, ra_shift=sh)
        pm.petal_means(bs)
        Dm_s, _, _, okc_s, _ = pm.design(meta['ra'], meta['dec'], meta['targetid'], meta['survey'],
                                         meta['program'], bs)
        u_s = okO & okc_s
        act = np.abs(Dm_s[u_s]).sum(0) > 1.0
        if act.sum() == 0:
            continue
        Xs = Dm_s[u_s][:, act]
        ys, es = OR[u_s], eOR[u_s] * errscale
        ks = np.abs(ys / es) < 8
        bb_, cc_, chi2s, dofs = wls(ys[ks], es[ks], Xs[ks])
        shifts[f'{sh:+.0f}'] = dict(delta_chi2=float(np.sum((ys[ks] / es[ks]) ** 2) - chi2s), n_par=int(act.sum()),
                                     n=int(ks.sum()))
    # same clip (|chi|<8) for the true trace, for a like-for-like comparison
    k8 = np.abs(y / e) < 8
    b8, c8, chi28, dof8 = wls(y[k8], e[k8], Xd[k8])
    det['delta_chi2_clip8'] = float(np.sum((y[k8] / e[k8]) ** 2) - chi28)
    det['shifted_trace_null'] = shifts
    pm.petal_means(basis)
    out['detection'] = det
    out['forward_model'] = dict(basis=names, amp_R=b.tolist(), err_R=eb.tolist(), chi2=chi2, dof=dof,
                                n_used=int(keep.sum()), n_clipped=int((~keep).sum()),
                                note='peak [OIII]5007 surface brightness (R) of Gaussian ridges at Dec nodes')
    pred_full = np.full(N, np.nan)
    pred_full[use] = Xd @ b
    T_full = np.full(N, np.nan)
    T_full[use] = X[use][:, active] @ b                 # model sky brightness at fibre (before subtraction)
    # also a single global amplitude scale on the E spine alone and comparison of model vs data
    # ---------------- neighbour instrument ----------------
    nbO, nn = neighbour_map(meta['ra'], meta['dec'], OR, eOR, okO, radius=0.12)
    # ---------------- per-fibre table ----------------
    with open(D / 'spectra_desi_fiber_table.csv', 'w', newline='') as fh:
        wr = csv.writer(fh)
        hdr = ['sparcl_id', 'targetid', 'survey', 'program', 'spectype', 'ra', 'dec', 'dRA_spine', 'd_any',
               'ffrac', 'OIII_fix_1e17', 'OIII_fix_err', 'OIII_R', 'OIII_R_err', 'OIII_v_free', 'OIII_v_err',
               'OIII_dchi_free', 'model_resid_R', 'model_sky_R', 'nbr_OIII_R', 'n_tiles', 'ok']
        wr.writerow(hdr)
        for i in range(N):
            wr.writerow([res['sparcl_id'][i], meta['targetid'][i], meta['survey'][i], meta['program'][i],
                         meta['spectype'][i], f"{meta['ra'][i]:.5f}", f"{meta['dec'][i]:.5f}",
                         f"{meta['dRA_spine'][i]:.4f}", f"{meta['d_any'][i]:.4f}", f'{ffrac[i]:.3f}',
                         f'{O[i]:.3f}', f'{eO[i]:.3f}', f'{OR[i]:.3f}', f'{eOR[i]:.3f}',
                         f"{res[g + '__v_free'][i]:.1f}", f"{res[g + '__ev_free'][i]:.1f}",
                         f"{res[g + '__dchi_free'][i]:.1f}", f'{pred_full[i]:.3f}', f'{T_full[i]:.3f}',
                         f'{nbO[i]:.3f}', len(info[i]), int(okO[i])])
    # ---------------- velocities ----------------
    v, ev = res[g + '__v_free'], res[g + '__ev_free']
    sn = np.abs(O / eO)
    vel = {}
    for lab, m in (('positive_residual_SN>5', okO & (O / eO > 5)), ('negative_residual_SN>5', okO & (O / eO < -5)),
                   ('all_|SN|>5', okO & (sn > 5))):
        mm = m & np.isfinite(ev) & ~res[g + '__v_edge'] & (ev < 30)
        if mm.sum() > 2:
            w = 1 / ev[mm] ** 2
            vm = np.sum(w * v[mm]) / w.sum()
            vel[lab] = dict(n=int(mm.sum()), wmean=float(vm), err_formal=float(w.sum() ** -0.5),
                            median=float(np.median(v[mm])), std=float(np.std(v[mm])),
                            err_scatter=float(np.std(v[mm]) / np.sqrt(mm.sum())))
    out['OIII_velocity_bary'] = vel
    # ---- barycentric-correction test: v([OIII]) should not depend on the exposures' v_corr ----
    vcorr = fibre_vcorr(pm, meta, info)
    mm = okO & (sn > 5) & np.isfinite(ev) & ~res[g + '__v_edge'] & (ev < 30) & np.isfinite(vcorr)
    if mm.sum() > 10:
        X1 = np.c_[np.ones(mm.sum()), vcorr[mm]]
        evs = np.hypot(ev[mm], 5.0)       # 5 km/s allowance for intrinsic/positional velocity scatter
        bb, cc, chi2v, dofv = wls(v[mm], evs, X1)
        sc = max(chi2v / dofv, 1)
        out['OIII_v_vs_vcorr'] = dict(n=int(mm.sum()), intercept=float(bb[0]), slope=float(bb[1]),
                                      err_intercept=float(np.sqrt(cc[0, 0] * sc)),
                                      err_slope=float(np.sqrt(cc[1, 1] * sc)), red_chi2=float(chi2v / dofv),
                                      vcorr_range=[float(vcorr[mm].min()), float(vcorr[mm].max())],
                                      note='slope 0 = barycentric correction applied; -1 = missing')
    # ---- by program ----
    byprog = {}
    for p in ('dark', 'bright'):
        m = okO & (sn > 5) & np.isfinite(ev) & ~res[g + '__v_edge'] & (ev < 30) & (meta['program'] == p)
        if m.sum() > 3:
            w = 1 / np.hypot(ev[m], 5.0) ** 2
            byprog[p] = dict(n=int(m.sum()), wmean=float(np.sum(w * v[m]) / w.sum()),
                             err=float(np.std(v[m]) / np.sqrt(m.sum())))
    out['OIII_velocity_by_program'] = byprog
    # ---- velocity vs Dec along the filament (|S/N|>4 fibres within 1 deg of a strand) ----
    vd = []
    for lo in np.arange(13.0, 28.01, 1.0):
        m = okO & (sn > 4) & np.isfinite(ev) & ~res[g + '__v_edge'] & (ev < 40) & (meta['dec'] >= lo) & \
            (meta['dec'] < lo + 1) & (meta['d_any'] < 1.0)
        if m.sum() >= 3:
            w = 1 / np.hypot(ev[m], 5.0) ** 2
            vm = np.sum(w * v[m]) / w.sum()
            vd.append(dict(dec_lo=float(lo), n=int(m.sum()), v=float(vm),
                           err=float(max(w.sum() ** -0.5, np.std(v[m]) / np.sqrt(m.sum()))),
                           n_pos=int((O[m] > 0).sum()), n_neg=int((O[m] < 0).sum())))
    out['OIII_velocity_vs_dec'] = vd
    out['vcorr_note'] = 'v_corr = EFFTIME-weighted barycentric correction (astropy, KPNO) of contributing exposures'
    np.savez_compressed(D / 'spectra_desi_vcorr.npz', sparcl_id=res['sparcl_id'], vcorr=vcorr)
    np.savez_compressed(D / 'spectra_desi_model_arrays.npz', sparcl_id=res['sparcl_id'], pred=pred_full,
                        T=T_full, nbO=nbO, nn=nn, OR=OR, eOR=eOR, okO=okO, use=use, Dm=Dm, basis=np.array(basis.names))
    (D / 'spectra_desi_results.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ('n_spectra', 'null_tests', 'error_scale_from_nulls', 'OIII_chi_stats',
                                          'OIII_velocity_bary')}, indent=1))
    print('forward model:', [(n, round(a, 2), round(s, 2)) for n, a, s in zip(names, b, eb)], chi2, dof)


if __name__ == '__main__':
    main()
