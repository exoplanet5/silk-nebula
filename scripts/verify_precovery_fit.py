"""verify_precovery step 2: independent re-fit of the SDSS Field-table sky levels with NSNS-free templates.

Model per run (my own formulation):
    y[f,c] = S_f + g_c + sum_p h_{c,p} t^p (p=1..2) [+ extra regressors x per-run coefficients] + A * T[f,c]
S_f common mode per field number (simultaneous in the 6 camcols), g_c camcol offsets, quadratic per-camcol drifts,
T = footprint-averaged ridge template (verify_precovery_prep.py). Iterative 4-sigma MAD clipping.
Window: |d_E| <= 2.5 deg from the E spine, inside its Dec range; quality >= 2; no bright-star / planet fields.
A run is used if it crosses: >= 6 fields at d < -0.7, >= 6 at d > +0.7, >= 4 with |d| < 0.3.

Errors / nulls (all per run, then combined over runs):
  formal    : least-squares error from the clipped residual rms;
  null_in   : template shifted by s = +-0.8..2.0 deg in RA, fitted JOINTLY with the unshifted template;
  null_off  : the whole analysis repeated in windows centred on the shifted spine (|s| = 3..6 deg),
              excluding fields within 1.2 deg of the real E/W strands;  the same windows get injections.
Combined amplitude: inverse-variance mean with e_eff = max(formal, std(null_in)); error scaled by sqrt(chi2/dof)>=1.
Units: A [nmgy arcsec^-2] / CONV_b, CONV_b = response of band b to 1 R of [OIII]5007 for the filament spectrum
(lib.band_response_per_R5007).  For 'desi' / 'fibre' templates (already in R) A/CONV is a scale factor.
Output: data/derived/verify_precovery_fit.json
"""
import sys, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from verify_precovery_lib import DER, design_nuisance, robust_fit, band_response_per_R5007, nmgy_per_R

Z = dict(np.load(DER / 'verify_precovery_fields.npz'))
SH = Z['shifts']
CONV = {b: band_response_per_R5007(b)[0] for b in 'ugriz'}
CONV_OIII_ONLY = nmgy_per_R('g', 5008.24) + (4960.30 / 5008.24 / 2.98) * nmgy_per_R('g', 4960.30)
CONV_OII_U = nmgy_per_R('u', 3728.48)          # 1 R of the [OII] doublet in u (nmgy arcsec^-2)
W = 2.5
IN_SHIFTS = [s for s in SH if 0.79 <= abs(s) <= 2.01]
OFF_SHIFTS = [s for s in SH if abs(s) >= 2.99]


def est_key(band, est):
    return {'psp': f'f_sky_{band}', 'frames': f'f_skyFrames_{band}', 'framessub': f'f_skyFramesSub_{band}'}[est]


_BIGMASK = None


def big_star_mask():
    """Enlarged bright-star exclusion: 2.0 deg (V<1.5), 1.0 deg (V<4.5), 0.5 deg (V<6), 0.25 deg (V<7)."""
    global _BIGMASK
    if _BIGMASK is None:
        from verify_precovery_prep import tycho_bright
        from verify_precovery_lib import uvec
        from scipy.spatial import cKDTree
        st = tycho_bright(7.0)
        t = cKDTree(uvec(Z['ra'], Z['dec']))
        bad = np.zeros(len(Z['ra']), bool)
        for ra, de, v in st:
            r = 2.0 if v < 1.5 else 1.0 if v < 4.5 else 0.5 if v < 6 else 0.25
            bad[t.query_ball_point(uvec(ra, de), np.radians(r))] = True
        _BIGMASK = bad
    return _BIGMASK


STARMASK = {'mode': 'standard'}


def base_mask(extra_quality=2):
    m = (Z['quality'] >= extra_quality) & ~Z['star_bad'] & (Z['planet_sep'] > 1.0)
    if STARMASK['mode'] == 'big':
        m &= ~big_star_mask()
    return m


def select(run, ishift=None, quality=2, window=W):
    m = (Z['run'] == run) & base_mask(quality)
    if ishift is None:
        d = Z['dE']; o = Z['overE']
    else:
        d = Z['dE_shift'][ishift]; o = Z['overE_shift'][ishift]
        m &= (np.abs(Z['dE']) > 1.2) & (np.abs(Z['dW']) > 1.2)
    m &= (np.abs(d) <= window) & (o == 0)
    idx = np.where(m)[0]
    dd = d[idx]
    ok = (dd < -0.7).sum() >= 6 and (dd > 0.7).sum() >= 6 and (np.abs(dd) < 0.3).sum() >= 4
    ok &= len(np.unique(Z['field'][idx])) >= 8 and len(np.unique(Z['camcol'][idx])) >= 3
    return idx if ok else None


def moon_regressor(idx):
    alt = np.clip(Z['alt'][idx], 5, 90); X = 1 / np.sin(np.radians(alt))
    malt = Z['moon_alt'][idx]
    if np.nanmedian(malt) <= -2:
        return np.zeros(len(idx))
    rho = np.clip(Z['moon_sep'][idx], 2, 180)
    f = 10 ** 5.36 * (1.06 + np.cos(np.radians(rho)) ** 2) + 10 ** (6.15 - rho / 40)
    Xm = 1 / np.sin(np.radians(np.clip(malt, 2, 90)))
    k = 0.18
    return f * 10 ** (-0.4 * k * Xm) * (1 - 10 ** (-0.4 * k * X)) * Z['moon_illum'][idx]


def extra_columns(idx, which):
    cols = []
    for w in which:
        if w == 'airmass':
            v = Z['f_airmass_g'][idx]
        elif w == 'moon':
            v = moon_regressor(idx)
        elif w == 'zl5':
            v = np.exp(-np.abs(Z['ecl_beta'][idx]) / 5.0)
        elif w == 'zl15':
            v = np.exp(-np.abs(Z['ecl_beta'][idx]) / 15.0)
        elif w == 'gegen':
            v = np.exp(-0.5 * ((180 - np.abs(Z['ecl_dlam'][idx])) ** 2 + Z['ecl_beta'][idx] ** 2) / 10.0 ** 2)
        elif w == 'twi':        # twilight: grows ~0.4 dex per degree of solar altitude above -18
            v = 10 ** (0.4 * np.clip(Z['sun_alt'][idx] + 18.0, -10, 10))
        elif w == 'iris':
            v = Z['iris100'][idx]
        elif w == 'halpha':
            v = Z['halpha_f03'][idx]
        elif w == 'cubic':      # per-camcol cubic drift (extra flexibility)
            f = Z['field'][idx]; t = ((f - np.median(f)) / 50.0) ** 3
            cu, ci = np.unique(Z['camcol'][idx], return_inverse=True)
            D = np.zeros((len(idx), len(cu))); D[np.arange(len(idx)), ci] = t
            cols.append(D[:, 1:]); continue
        else:
            raise ValueError(w)
        v = np.nan_to_num(v - np.nanmean(v))
        if np.ptp(v) > 0:
            cols.append(v[:, None])
    return np.hstack(cols) if cols else np.zeros((len(idx), 0))


def fit(idx, y, templates, extra=(), drift=2, common=True, clip=4.0):
    f = Z['field'][idx]; c = Z['camcol'][idx]
    X = design_nuisance(f, c, drift, common)
    if not common:
        from verify_precovery_lib import smooth_time_basis
        X = np.hstack([X, smooth_time_basis(f.astype(float))])
    E = extra_columns(idx, extra)
    T = np.column_stack(templates)
    X = np.hstack([X, E, T])
    return robust_fit(y[idx], X, T.shape[1], clip=clip)


def tmpl(name, idx, stat='mean'):
    return Z[f'T_{name}_{stat}'][idx]


def run_list():
    return np.unique(Z['run'])


def per_run(band='g', est='psp', tname='gauss028', stat='mean', extra=(), drift=2, common=True, quality=2,
            do_nulls=True, inject=None, window=W, conv=None, starmask='standard'):
    STARMASK['mode'] = starmask
    y = Z[est_key(band, est)].astype(float)
    conv = CONV[band] if conv is None else conv
    shname = 'gauss028' if tname not in ('desi',) else 'desi'
    out = []
    for run in run_list():
        idx = select(run, None, quality, window)
        if idx is None:
            continue
        T0 = tmpl(tname, idx, stat)
        if np.ptp(T0) < 0.05 * np.nanmax(Z[f'T_{tname}_mean']):
            continue
        r = fit(idx, y, [T0], extra, drift, common)
        rec = dict(run=int(run), mjd=float(np.median(Z['mjd'][idx])),
                   dec_cross=float(np.median(Z['dec'][idx][np.abs(Z['dE'][idx]) < 0.3])),
                   ra_cross=float(np.median(Z['ra'][idx][np.abs(Z['dE'][idx]) < 0.3])),
                   n=r['n'], n_fields=int(len(np.unique(Z['field'][idx]))), rms=r['rms'] / conv,
                   A=float(r['coef'][0] / conv), eA_formal=float(np.sqrt(r['cov'][0, 0]) / conv),
                   sky=float(np.median(y[idx])), moon_alt=float(np.median(Z['moon_alt'][idx])),
                   moon_illum=float(np.median(Z['moon_illum'][idx])), moon_sep=float(np.median(Z['moon_sep'][idx])),
                   ecl_beta=float(np.median(Z['ecl_beta'][idx])), ecl_dlam=float(np.median(Z['ecl_dlam'][idx])),
                   Tptp=float(np.ptp(T0)))
        if do_nulls:
            # templates without their own shifted copies use the shifted unit Gaussian; its null amplitudes are
            # rescaled by the ratio of formal errors so that null/err keeps its meaning in the template's units
            scale = 1.0
            if tname not in ('gauss028', 'desi'):
                rg = fit(idx, y, [tmpl('gauss028', idx)], extra, drift, common)
                scale = rec['eA_formal'] / (np.sqrt(rg['cov'][0, 0]) / conv)
            rec['null_scale'] = float(scale)
            nin = []
            for s in IN_SHIFTS:
                i = int(np.where(np.isclose(SH, s))[0][0])
                Ts = Z[f'TS_{shname}'][i][idx]
                # skip shifts whose template is (nearly) absent from the window: unconstrained coefficient
                if np.ptp(Ts) < 0.3 * np.ptp(T0) or np.sum(Ts > 0.3 * np.nanmax(T0)) < 4:
                    nin.append(np.nan); continue
                rr = fit(idx, y, [T0, Ts], extra, drift, common)
                nin.append(scale * rr['coef'][1] / conv)
            nin = np.array(nin)
            if np.isfinite(nin).sum() < 5:
                nin[:] = np.nan
            rec['null_in'] = [float(v) for v in nin]
            rec['e_null_in'] = float(1.4826 * np.nanmedian(np.abs(nin - np.nanmedian(nin)))) if np.isfinite(nin).any() else None
            rec['e_null_in_std'] = float(np.nanstd(nin)) if np.isfinite(nin).any() else None
            noff = {}
            for s in OFF_SHIFTS:
                i = int(np.where(np.isclose(SH, s))[0][0])
                jdx = select(run, i, quality, window)
                if jdx is None:
                    continue
                Ts = Z[f'TS_{shname}'][i][jdx]
                rr = fit(jdx, y, [Ts], extra, drift, common)
                ent = dict(A=float(scale * rr['coef'][0] / conv), e=float(scale * np.sqrt(rr['cov'][0, 0]) / conv))
                if inject is not None:
                    yy = y.copy()
                    for amp in inject:
                        yy[jdx] = y[jdx] + amp * conv * Ts
                        ri = fit(jdx, yy, [Ts], extra, drift, common)
                        ent[f'inj{amp:g}'] = float(ri['coef'][0] / conv)
                noff[f'{s:+.1f}'] = ent
            rec['null_off'] = noff
            v = np.array([e['A'] for e in noff.values()])
            rec['e_null_off'] = float(np.std(v)) if len(v) >= 5 else None
        out.append(rec)
    return out


def combine(runs, key='A', ekey='e_eff', sel=None):
    rr = [r for r in runs if (sel is None or sel(r))]
    if len(rr) < 2:
        return None
    A = np.array([r[key] for r in rr]); e = np.array([r[ekey] for r in rr])
    w = 1 / e ** 2
    m = float(np.sum(w * A) / np.sum(w)); em = float(1 / np.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (A - m) ** 2)); dof = len(A) - 1
    return dict(A=m, eA=em * max(1.0, np.sqrt(chi2 / dof)), eA_raw=em, chi2_dof=chi2 / dof, n=len(A),
                n_pos=int((A > 0).sum()), median=float(np.median(A)))


def add_eff(runs, mode='max'):
    """max: max(formal, robust in-window null scatter); off: max(formal, off-window null std) when >= 5 off
    windows exist, else 'max'; formal: least-squares error only."""
    for r in runs:
        en = r.get('e_null_in') or 0.0
        eo = r.get('e_null_off')
        if mode == 'formal':
            r['e_eff'] = r['eA_formal']
        elif mode == 'off' and eo:
            r['e_eff'] = max(r['eA_formal'], eo)
        else:
            r['e_eff'] = max(r['eA_formal'], en)
    return runs


def null_ensembles(runs, sel):
    """Combined statistic for the same shift applied to all runs (in-window joint and off-window)."""
    rr = [r for r in runs if sel(r)]
    w = np.array([1 / r['e_eff'] ** 2 for r in rr])
    res = {}
    res['in'] = []
    for j in range(len(IN_SHIFTS)):
        v = np.array([r['null_in'][j] for r in rr])
        k = np.isfinite(v)
        if k.sum() >= len(rr) // 2:
            res['in'].append(float(np.sum(w[k] * v[k]) / np.sum(w[k])))
    off = {}
    for s in OFF_SHIFTS:
        k = f'{s:+.1f}'
        have = [(i, r['null_off'][k]['A']) for i, r in enumerate(rr) if k in r['null_off']]
        if len(have) >= max(3, len(rr) // 2):
            ii = np.array([h[0] for h in have]); v = np.array([h[1] for h in have])
            off[k] = dict(A=float(np.sum(w[ii] * v) / np.sum(w[ii])), n=len(have))
    res['off'] = off
    # random-shift bootstrap: each run draws an independent off-window null (or in-window if none)
    rng = np.random.default_rng(1)
    pools = []
    for r in rr:
        p = [v['A'] for v in r['null_off'].values()]
        q = np.array(r['null_in'], float); q = q[np.isfinite(q)]
        pools.append(np.array(p) if len(p) >= 3 else (q if len(q) else np.array([0.0])))
    draws = np.array([[rng.choice(p) for p in pools] for _ in range(4000)])
    comb = (draws * w).sum(1) / w.sum()
    res['random_off_std'] = float(np.std(comb)); res['random_off_mean'] = float(np.mean(comb))
    pin = [np.array(r['null_in'])[np.isfinite(r['null_in'])] for r in rr]
    pin = [p if len(p) else np.array([0.0]) for p in pin]
    draws = np.array([[rng.choice(p) for p in pin] for _ in range(4000)])
    comb = (draws * w).sum(1) / w.sum()
    res['random_in_std'] = float(np.std(comb)); res['random_in_mean'] = float(np.mean(comb))
    return res


CORE = lambda r: 12.0 <= r['dec_cross'] < 20.5
NORTH = lambda r: r['dec_cross'] >= 20.5


def summarize(runs, label):
    s = dict(label=label, n_runs=len(runs))
    add_eff(runs, 'formal'); s['core_formal'] = combine(runs, sel=CORE)
    add_eff(runs, 'off'); s['core_offerr'] = combine(runs, sel=CORE)
    if s['core_offerr'] and all('null_in' in r for r in runs):
        s['core_offerr_nulls'] = null_ensembles(runs, CORE)
    add_eff(runs, 'max')
    s['core'] = combine(runs, sel=CORE)
    s['north'] = combine(runs, sel=NORTH)
    s['all'] = combine(runs)
    if s['core'] and all('null_in' in r for r in runs):
        ne = null_ensembles(runs, CORE)
        s['core_nulls'] = ne
        offv = np.array([v['A'] for v in ne['off'].values()])
        s['core_sig_vs_random_off'] = s['core']['A'] / ne['random_off_std']
        s['core_sig_vs_random_in'] = s['core']['A'] / ne['random_in_std']
        if len(offv) >= 3:
            s['core_off_same_shift_std'] = float(np.std(offv)); s['core_off_same_shift_mean'] = float(np.mean(offv))
            s['core_sig_vs_off_same_shift'] = s['core']['A'] / float(np.std(offv))
    return s


def fmt(s):
    c = s.get('core')
    if not c:
        return f"{s['label']:44s} (no core runs)"
    t = (f"{s['label']:44s} core A={c['A']:+.3f}+-{c['eA']:.3f} (chi2/dof {c['chi2_dof']:.2f}, n={c['n']}, pos {c['n_pos']}, "
         f"med {c['median']:+.2f})")
    if 'core_sig_vs_random_off' in s:
        t += f"  S/N vs nulls: off-rand {s['core_sig_vs_random_off']:.1f}, in-rand {s['core_sig_vs_random_in']:.1f}"
        if 'core_sig_vs_off_same_shift' in s:
            t += f", off-same {s['core_sig_vs_off_same_shift']:.1f}"
    o = s.get('core_offerr')
    if o:
        t += f" || off-err A={o['A']:+.2f}+-{o['eA']:.2f} chi2 {o['chi2_dof']:.2f}"
        if 'core_offerr_nulls' in s:
            ne = s['core_offerr_nulls']; ov = np.array([v['A'] for v in ne['off'].values()])
            t += f" S/N rand-off {o['A'] / ne['random_off_std']:.1f} same-off {o['A'] / np.std(ov):.1f}" if len(ov) > 2 else ''
    n = s.get('north')
    if n:
        t += f" | north A={n['A']:+.2f}+-{n['eA']:.2f} (n={n['n']})"
    return t


def main():
    res = {'conv_nmgy_arcsec2_per_R5007': CONV, 'conv_g_OIII_only': CONV_OIII_ONLY, 'conv_u_OII_only': CONV_OII_U, 'window_deg': W,
           'in_shifts': IN_SHIFTS, 'off_shifts': OFF_SHIFTS, 'cases': {}}
    cases = []
    # (a) templates, PSP estimator, g
    for t in ['gauss028', 'gauss028_EWW2', 'gaussE028', 'gaussfw0.60', 'gaussfw1.00', 'desi', 'fibre']:
        cases.append((f'g_psp_{t}', dict(band='g', est='psp', tname=t)))
    cases.append(('g_psp_gauss028_median', dict(band='g', est='psp', tname='gauss028', stat='med')))
    cases.append(('g_psp_gauss028_centre', dict(band='g', est='psp', tname='gauss028', stat='centre')))
    # estimators
    for e in ['frames', 'framessub']:
        for t in ['gauss028', 'desi']:
            cases.append((f'g_{e}_{t}', dict(band='g', est=e, tname=t)))
    # model variants
    cases.append(('g_psp_gauss028_drift1', dict(band='g', est='psp', tname='gauss028', drift=1)))
    cases.append(('g_psp_gauss028_cubic', dict(band='g', est='psp', tname='gauss028', extra=('cubic',))))
    cases.append(('g_psp_gauss028_nocommon', dict(band='g', est='psp', tname='gauss028', common=False)))
    cases.append(('g_psp_gauss028_q3', dict(band='g', est='psp', tname='gauss028', quality=3)))
    cases.append(('g_psp_gauss028_win1.5', dict(band='g', est='psp', tname='gauss028', window=1.5)))
    cases.append(('g_psp_gauss028_win4', dict(band='g', est='psp', tname='gauss028', window=4.0)))
    # confounders
    cases.append(('g_psp_gauss028_conf_all', dict(band='g', est='psp', tname='gauss028',
                                                  extra=('airmass', 'twi', 'zl5', 'zl15', 'gegen', 'iris'))))
    for x in ['airmass', 'twi', 'zl5', 'gegen', 'iris']:
        cases.append((f'g_psp_gauss028_conf_{x}', dict(band='g', est='psp', tname='gauss028', extra=(x,))))
    for t in ['gauss028_EWW2', 'gaussE028', 'gaussfw0.60', 'gaussfw1.00', 'fibre']:
        cases.append((f'g_framessub_{t}', dict(band='g', est='framessub', tname=t)))
    cases.append(('g_framessub_gauss028_median', dict(band='g', est='framessub', tname='gauss028', stat='med')))
    cases.append(('g_framessub_gauss028_drift1', dict(band='g', est='framessub', tname='gauss028', drift=1)))
    cases.append(('g_framessub_gauss028_nocommon', dict(band='g', est='framessub', tname='gauss028', common=False)))
    cases.append(('g_framessub_gauss028_win4', dict(band='g', est='framessub', tname='gauss028', window=4.0)))
    cases.append(('g_framessub_gauss028_conf_all', dict(band='g', est='framessub', tname='gauss028',
                                                       extra=('airmass', 'twi', 'zl5', 'zl15', 'gegen', 'iris'))))
    cases.append(('u_psp_desi', dict(band='u', est='psp', tname='desi')))
    cases.append(('u_framessub_desi', dict(band='u', est='framessub', tname='desi')))
    cases.append(('u_psp_gauss028_conf_all', dict(band='u', est='psp', tname='gauss028',
                                                 extra=('airmass', 'twi', 'zl5', 'zl15', 'gegen', 'iris'))))
    cases.append(('u_frames_gauss028', dict(band='u', est='frames', tname='gauss028')))
    # u band in units of R of [OII]3726+3729 (photons) and g band in units of R of [OIII]5007 (+4959 at 1/2.98)
    for e in ('psp', 'framessub', 'frames'):
        cases.append((f'u_{e}_gauss028_OIIunits', dict(band='u', est=e, tname='gauss028', conv=CONV_OII_U)))
        cases.append((f'g_{e}_gauss028_OIIIunits', dict(band='g', est=e, tname='gauss028', conv=CONV_OIII_ONLY)))
    # enlarged bright-star exclusion (kappa Gem 0.25 deg E of the spine at Dec 24.4; V=5.0 star on the ridge at Dec 18.5)
    for e in ('psp', 'framessub'):
        cases.append((f'u_{e}_gauss028_OIIunits_bigmask', dict(band='u', est=e, tname='gauss028', conv=CONV_OII_U, starmask='big')))
        cases.append((f'g_{e}_gauss028_OIIIunits_bigmask', dict(band='g', est=e, tname='gauss028', conv=CONV_OIII_ONLY, starmask='big')))
        cases.append((f'r_{e}_gauss028_bigmask', dict(band='r', est=e, tname='gauss028', starmask='big')))
        cases.append((f'i_{e}_gauss028_bigmask', dict(band='i', est=e, tname='gauss028', starmask='big')))
    # other bands (controls / [OII] test)
    for b in 'urizg':
        for e in ['psp', 'framessub']:
            if b == 'g' and e == 'psp':
                continue
            cases.append((f'{b}_{e}_gauss028', dict(band=b, est=e, tname='gauss028')))
    cases.append(('r_psp_gauss028_halpha', dict(band='r', est='psp', tname='gauss028', extra=('halpha',))))
    # injection run (off windows) for the baseline
    cases.insert(0, ('g_psp_gauss028', dict(band='g', est='psp', tname='gauss028', inject=(1.0, 3.0))))
    cases.insert(1, ('g_psp_desi', dict(band='g', est='psp', tname='desi', inject=(1.0,))))
    seen = set()
    for label, kw in cases:
        if label in seen:
            continue
        seen.add(label)
        runs = per_run(**kw)
        s = summarize(runs, label)
        s['kwargs'] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in kw.items()}
        s['runs'] = runs
        res['cases'][label] = s
        print(fmt(s), flush=True)
    (DER / 'verify_precovery_fit.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
