"""w3_judge_skeptic: two quick falsification tests of the consensus (truncated shock on the Monogem limb).

A. How special is the "Monogem-limb registration" that carries the 252-318 pc distance?
   Observed: the southern 23 spine points (Dec 12.0-17.5) lie at 12.10 +- 0.12 deg (rms) from the Knies+24
   two-SNR centre K2 (R = 11.00 +- 1.1 deg). Null: keep the real trace fixed and place a circle of the same
   radius at random centres around the trace midpoint (equivalent to random rigid translation + spin of the
   trace about a fixed Monogem). Event E: some contiguous 23-point window has rms(distance) <= RMS_MAX and
   |mean distance - R| <= TOL. Reported unconditionally (for annuli of centre distance) and conditional on the
   trace entering the rim zone (some point within 1.5 deg of the circle), which is the honest comparison:
   "given a long [O III] filament that reaches Monogem's rim, how often is part of it tangent to the rim?"

B. Cross-filament ordering of ionisation stages in the southern Silk (single-strand Dec 12-15 and Dec 12-17.6).
   A shock moving outward (east, away from Monogem) puts the hot C IV / FUV zone EAST of [O III] (shockmodel
   prediction: FUV ridge ~0.1 deg outside) and [O II] WEST (downstream) of [O III]. A hot-gas interface
   (conduction front, mixing layer, NEI cooling of Monogem gas) has its hot side WEST, so FUV west of [O III]
   and [O II] east. A static heated layer gives no ordering. Measured with the w3_figs F6 machinery (NSNS [O III]
   1' map, GALEX FUV/NUV 1' skybg blocks, DESI fibre line fluxes), split-normal + linear-baseline fits, and
   along-spine jackknife (10 chunks) errors.

Inputs (all on disk, read-only): silk_trace.json, w3_figs_f1_maps.npz, GALEX blocks via
verify_xrayfuv_galex_measure.prep_band, verify_desi_linefluxes.npz, verify_desi_coadd_meta.npz,
verify_desi_detect.json.
Outputs: data/derived/w3_judge_skeptic_tests.json, figures/w3_judge_skeptic_tests.png
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_skeptic_tests.py
"""
import sys
import json
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import (CurveCoords, radec_to_vec, load_spine, CENTRES, sep_pa, offset_point, DER)
import w3_figs_f6_profiles_data as F6   # module import only (main() is guarded)

OUT_JSON = DER / 'w3_judge_skeptic_tests.json'
OUT_FIG = '/Users/mickey/astronomy/silk-nebula/figures/w3_judge_skeptic_tests.png'
rng = np.random.default_rng(20261117)

ra_sp, dec_sp, _ = load_spine()
CC0 = F6.CC0
S15 = float(np.interp(15.0, CC0.dec, CC0.s))
S176 = F6.S_SPLIT


# ----------------------------------------------------------------------------------------------- test A
def window_stats(d, w=23):
    c = np.cumsum(np.r_[0, d]); c2 = np.cumsum(np.r_[0, d * d])
    m = (c[w:] - c[:-w]) / w
    var = (c2[w:] - c2[:-w]) / w - m * m
    return m, np.sqrt(np.maximum(var, 0))


def test_A(n=40000, rms_max=0.15, tol=1.1, rim_zone=1.5):
    res = {}
    for key, R in (('K2', 11.00), ('K1', 11.97)):
        cra, cde = CENTRES[key]
        d_obs, _ = sep_pa(cra, cde, ra_sp, dec_sp)
        m_obs, r_obs = window_stats(d_obs)
        south = dec_sp <= 17.5 + 1e-6
        obs = dict(south_mean=float(d_obs[south].mean()), south_rms=float(d_obs[south].std()),
                   best_window_rms=float(r_obs.min()),
                   best_window_mean=float(m_obs[np.argmin(r_obs)]),
                   n_windows_event=int(np.sum((r_obs <= rms_max) & (np.abs(m_obs - R) <= tol))),
                   observed_event=bool(np.any((r_obs <= rms_max) & (np.abs(m_obs - R) <= tol))))
        # real midpoint-centre separation
        imid = len(ra_sp) // 2
        rho_obs, pa_obs = sep_pa(ra_sp[imid], dec_sp[imid], cra, cde)
        obs['rho_mid_to_centre'] = float(rho_obs)
        # random centres around the trace midpoint, uniform in area within [0, 26] deg
        rho = np.degrees(np.arccos(1 - rng.random(n) * (1 - np.cos(np.radians(26.0)))))
        pa = rng.random(n) * 360
        cr, cd = offset_point(ra_sp[imid], dec_sp[imid], rho, pa)
        ev = np.zeros(n, bool); touch = np.zeros(n, bool); cross = np.zeros(n, bool)
        for i in range(n):
            d, _ = sep_pa(cr[i], cd[i], ra_sp, dec_sp)
            m, r = window_stats(d)
            ev[i] = np.any((r <= rms_max) & (np.abs(m - R) <= tol))
            touch[i] = np.any(np.abs(d - R) <= rim_zone)
            cross[i] = (d.min() < R) & (d.max() > R)
        out = {}
        for lo, hi in ((0, 26), (6, 20), (10, 16)):
            sel = (rho >= lo) & (rho < hi)
            out[f'rho_{lo}_{hi}'] = dict(n=int(sel.sum()), p_event=float(ev[sel].mean()),
                                         p_touch=float(touch[sel].mean()),
                                         p_event_given_touch=float(ev[sel & touch].mean()) if (sel & touch).any() else None,
                                         p_event_given_cross=float(ev[sel & cross].mean()) if (sel & cross).any() else None)
        res[key] = dict(R=R, rms_max=rms_max, tol=tol, rim_zone=rim_zone, observed=obs, null=out, n=n)
        print(key, json.dumps(res[key], indent=None), flush=True)
    # tangency angle of the southern spine vs the K2 circle tangent
    cra, cde = CENTRES['K2']
    s = dec_sp <= 17.5 + 1e-6
    _, pa_c = sep_pa(cra, cde, ra_sp[s], dec_sp[s])           # PA of radius vector at each point
    dra = np.gradient(ra_sp[s]) * np.cos(np.radians(dec_sp[s])); dde = np.gradient(dec_sp[s])
    pa_t = np.degrees(np.arctan2(dra, dde)) % 180               # PA of the spine direction (mod 180)
    ang = ((pa_t - (pa_c % 180)) + 90) % 180 - 90               # 0 means spine along the radius
    res['south_angle_to_radius_deg'] = dict(median=float(np.median(np.abs(ang))), values=np.round(ang, 1).tolist())
    return res


# ----------------------------------------------------------------------------------------------- test B
def splitnorm(p, x):
    A, x0, sw, se, b0, b1 = p
    sig = np.where(x < x0, sw, se)
    return A * np.exp(-0.5 * ((x - x0) / sig) ** 2) + b0 + b1 * x


def fit_split(x, y, xmax=0.9, p0=None):
    m = np.isfinite(y) & (np.abs(x) <= xmax)
    if m.sum() < 12:
        return None
    if p0 is None:
        p0 = [np.nanmax(y[m]) - np.nanmedian(y[m]), 0.0, 0.12, 0.12, np.nanmedian(y[m]), 0.0]
    lb = [0, -0.3, 0.02, 0.02, -np.inf, -np.inf]; ub = [np.inf, 0.3, 0.8, 0.8, np.inf, np.inf]
    try:
        r = least_squares(lambda p: splitnorm(p, x[m]) - y[m], p0, bounds=(lb, ub))
    except Exception:
        return None
    A, x0, sw, se, b0, b1 = r.x
    # first moment of the baseline-subtracted profile within x0 +- 0.4
    mm = m & (np.abs(x - x0) <= 0.4)
    yy = y[mm] - (b0 + b1 * x[mm])
    mom = float(np.sum(x[mm] * yy) / np.sum(yy)) if np.sum(yy) > 0 else np.nan
    return dict(A=float(A), x0=float(x0), sig_w=float(sw), sig_e=float(se), asym_w_over_e=float(sw / se),
                fwhm=float(1.1774 * (sw + se)), moment=mom)


def shift_fit(x, y, xt, yt, xmax=0.9):
    """Fit y(x) = a * T(x - delta) + b0 + b1 x, with T the (interpolated, baseline-free) template."""
    m = np.isfinite(y) & (np.abs(x) <= xmax)
    good = np.isfinite(yt)
    T = lambda xx: np.interp(xx, xt[good], yt[good], left=0.0, right=0.0)
    best = None
    for dl in np.arange(-0.3, 0.3001, 0.005):
        X = np.vstack([T(x[m] - dl), np.ones(m.sum()), x[m]]).T
        coef, res, *_ = np.linalg.lstsq(X, y[m], rcond=None)
        chi = np.sum((X @ coef - y[m]) ** 2)
        if best is None or chi < best[0]:
            best = (chi, dl, coef)
    return float(best[1]), float(best[2][0])


def nsns_data():
    M = np.load(DER / 'w3_figs_f1_maps.npz')
    o3 = M['oiii'].astype(float)
    RA, DE = np.meshgrid(M['ra'], M['dec'])
    RA, DE, V = RA.ravel(), DE.ravel(), o3.ravel()
    box = (RA > 113.0) & (RA < 120.5) & (DE > 10.5) & (DE < 19.0)   # southern Silk +- ~2 deg
    return RA[box], DE[box], V[box]


def prof_generic(ra, dec, val, edges, smin, smax, stat='mean', tile=None, demean_out=0.3):
    x, s = CC0.coords(radec_to_vec(ra, dec), maxdist_deg=1.6)
    sel = np.isfinite(x) & (s >= smin) & (s <= smax) & np.isfinite(val)
    v = val.copy()
    if tile is not None:       # GALEX: per-tile demeaning against |x| > demean_out
        for t in np.unique(tile[sel]):
            mm = (tile == t) & sel
            ref = mm & (np.abs(x) > demean_out)
            if ref.sum() > 20:
                v[mm] -= np.mean(v[ref])
            else:
                sel &= ~mm
    bi = np.digitize(x[sel], edges) - 1
    ok = (bi >= 0) & (bi < len(edges) - 1)
    nb = len(edges) - 1
    cnt = np.bincount(bi[ok], minlength=nb)
    if stat == 'mean':
        pr = np.bincount(bi[ok], v[sel][ok], nb) / np.maximum(cnt, 1)
    else:
        pr = np.full(nb, np.nan)
        vv = v[sel][ok]; bb = bi[ok]
        for k in range(nb):
            q = vv[bb == k]
            if q.size >= 5:
                pr[k] = np.median(q)
    pr[cnt < (20 if tile is not None else 5)] = np.nan
    return pr


def jackknife(fun, smin, smax, K=10):
    edges = np.linspace(smin, smax, K + 1)
    full = fun(smin, smax, None)
    reps = []
    for k in range(K):
        reps.append(fun(smin, smax, (edges[k], edges[k + 1])))
    return full, reps


def test_B():
    import verify_xrayfuv_galex_measure as GM
    from verify_desi_common import good_fit
    out = {}
    ra_n, de_n, v_n = nsns_data()
    G = {b: GM.prep_band(b) for b in ('FUV', 'NUV')}
    lf = np.load(DER / 'verify_desi_linefluxes.npz'); cm = np.load(DER / 'verify_desi_coadd_meta.npz')
    det = json.loads((DER / 'verify_desi_detect.json').read_text())
    ok3 = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & cm['ok_contrib'] & np.isfinite(cm['kR_5008'])
    ok2 = good_fit(lf, 'o2sum', scale=det['err_scale_from_nulls']) & cm['ok_contrib'] & np.isfinite(cm['kR_3727'])
    okd = ok3 & ok2
    d_o3 = lf['o3__amp'][:, 0] * cm['kR_5008']
    d_o2 = lf['o2sum__amp'][:, 0] * cm['kR_3727']
    e025 = np.arange(-1.2, 1.2001, 0.025); x025 = 0.5 * (e025[1:] + e025[:-1])
    e10 = np.arange(-1.5, 1.5001, 0.1); x10 = 0.5 * (e10[1:] + e10[:-1])
    store = {}
    for seg, (smin, smax) in (('dec12_15', (0.0, S15)), ('dec12_17p6', (0.0, S176))):
        segout = {}

        def run(band, smin=smin, smax=smax):
            def fun(a, b, drop):
                if band == 'NSNS':
                    vals = v_n.copy()
                    if drop is not None:
                        x, s = CC0.coords(radec_to_vec(ra_n, de_n), maxdist_deg=1.6)
                        vals = np.where((s >= drop[0]) & (s < drop[1]), np.nan, vals)
                    return prof_generic(ra_n, de_n, vals, e025, a, b, stat='median')
                if band in ('FUV', 'NUV'):
                    D = G[band]
                    vals = D['v'].astype(float).copy()
                    if drop is not None:
                        x, s = CC0.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.6)
                        vals = np.where((s >= drop[0]) & (s < drop[1]), np.nan, vals)
                    return prof_generic(D['ra'], D['dec'], vals, e025, a, b, stat='mean', tile=D['tile'])
                vals = (d_o3 if band == 'DESI_O3' else d_o2).copy()
                vals = np.where(okd, vals, np.nan)
                if drop is not None:
                    x, s = CC0.coords(radec_to_vec(cm['ra'], cm['dec']), maxdist_deg=1.6)
                    vals = np.where((s >= drop[0]) & (s < drop[1]), np.nan, vals)
                return prof_generic(cm['ra'], cm['dec'], vals, e10, a, b, stat='mean')
            return jackknife(fun, smin, smax)
        profs = {}
        for band in ('NSNS', 'FUV', 'NUV', 'DESI_O3', 'DESI_O2'):
            profs[band] = run(band)
        store[seg] = profs
        # fits
        for band, (full, reps) in profs.items():
            xx = x10 if band.startswith('DESI') else x025
            xmax = 1.2 if band.startswith('DESI') else 0.9
            f = fit_split(xx, full, xmax=xmax)
            fr = [fit_split(xx, r, xmax=xmax) for r in reps]
            fr = [q for q in fr if q is not None]
            K = len(fr)
            je = {}
            for key in ('x0', 'sig_w', 'sig_e', 'asym_w_over_e', 'moment', 'fwhm'):
                vals = np.array([q[key] for q in fr])
                je[key] = float(np.sqrt((K - 1) / K * np.sum((vals - vals.mean()) ** 2))) if K > 2 else np.nan
            segout[band] = dict(fit=f, jk_err=je, n_jk=K)
        # shifts relative to the NSNS [O III] template (baseline-free split-normal model of NSNS)
        fN = segout['NSNS']['fit']
        tx = np.linspace(-1.2, 1.2, 961)
        tmpl = splitnorm([fN['A'], fN['x0'], fN['sig_w'], fN['sig_e'], 0, 0], tx)
        for band in ('FUV', 'NUV'):
            full, reps = profs[band]
            dl, a = shift_fit(x025, full, tx, tmpl)
            dls = np.array([shift_fit(x025, r, tx, tmpl)[0] for r in reps])
            K = len(dls)
            segout[f'shift_{band}_vs_NSNS'] = dict(delta=dl, jk_err=float(np.sqrt((K - 1) / K * np.sum((dls - dls.mean()) ** 2))),
                                                   amp=a)
        # DESI [O II] vs DESI [O III] (same fibres)
        fullO3, repsO3 = profs['DESI_O3']; fullO2, repsO2 = profs['DESI_O2']
        f3 = segout['DESI_O3']['fit']
        if f3 is not None:
            t3 = splitnorm([f3['A'], f3['x0'], f3['sig_w'], f3['sig_e'], 0, 0], tx)
            dl, a = shift_fit(x10, fullO2, tx, t3, xmax=1.2)
            dls = []
            for r3, r2 in zip(repsO3, repsO2):
                q3 = fit_split(x10, r3, xmax=1.2)
                if q3 is None:
                    continue
                tt = splitnorm([q3['A'], q3['x0'], q3['sig_w'], q3['sig_e'], 0, 0], tx)
                dls.append(shift_fit(x10, r2, tx, tt, xmax=1.2)[0])
            dls = np.array(dls); K = len(dls)
            segout['shift_DESI_O2_vs_O3'] = dict(delta=dl, jk_err=float(np.sqrt((K - 1) / K * np.sum((dls - dls.mean()) ** 2))) if K > 2 else None,
                                                 amp_ratio=a)
        out[seg] = segout
        print(seg, json.dumps(segout, default=float), flush=True)
    return out, store, x025, x10


def figure(A, B, store, x025, x10):
    from silk_style import new_sky_figure, finish, label, OIII_C, ACCENT, WARN, FG, FGDIM, OK, HB_C
    import matplotlib.pyplot as plt
    fig, axes = new_sky_figure(13, 8.0, 'Silk skeptic tests',
                               'A: chance of a Monogem-rim tangency  ·  B: cross-filament ordering (south)',
                               nrows=2, ncols=2)
    # A1: distance from K2 along the spine
    ax = axes[0, 0]
    cra, cde = CENTRES['K2']
    d, _ = sep_pa(cra, cde, ra_sp, dec_sp)
    d1, _ = sep_pa(*CENTRES['K1'], ra_sp, dec_sp)
    ax.plot(dec_sp, d, color=OIII_C, lw=1.6, label='from K2 (2-SNR)')
    ax.plot(dec_sp, d1, color=WARN, lw=1.2, label='from K1 (1-SNR)')
    ax.axhspan(11.0 - 1.1, 11.0 + 1.1, color=OIII_C, alpha=0.12)
    ax.axhline(11.0, color=OIII_C, lw=0.8, ls=':')
    ax.axhline(11.97, color=WARN, lw=0.8, ls=':')
    ax.axvline(17.5, color=FGDIM, lw=0.6)
    ax.set_xlabel('spine Dec (deg)'); ax.set_ylabel('distance from centre (deg)')
    ax.legend(loc='upper left', ncol=2, fontsize=7.5)
    pk = A['K2']['null']
    label(ax, 12.2, 15.6, f"P(tangent window | trace reaches rim) = {pk['rho_0_26']['p_event_given_touch']:.2f} (K2)",
          size=7.5)
    label(ax, 12.2, 15.1, f"P(tangent window | trace crosses rim) = {pk['rho_0_26']['p_event_given_cross']:.2f} (K2)",
          size=7.5)
    ax.set_ylim(8.5, 16.5)
    # B panels
    for j, seg in enumerate(('dec12_15', 'dec12_17p6')):
        ax = axes[0, 1] if j == 0 else axes[1, 0]
        pr = store[seg]
        def norm(p, xx):
            m = np.isfinite(p)
            base = np.nanmedian(p[m & (np.abs(xx) > 0.6) & (np.abs(xx) < 1.2)])
            q = p - base
            return q / np.nanmax(q[np.abs(xx) < 0.3])
        ax.plot(x025, norm(pr['NSNS'][0], x025), color=OIII_C, lw=1.8, label='NSNS [O III]')
        ax.plot(x025, norm(pr['FUV'][0], x025), color='#b388ff', lw=1.3, label='GALEX FUV')
        ax.plot(x025, norm(pr['NUV'][0], x025), color=OK, lw=1.0, alpha=0.8, label='GALEX NUV')
        ax.plot(x10, norm(pr['DESI_O3'][0], x10), 'o', color=OIII_C, ms=3.5, mfc='none', label='DESI [O III]')
        ax.plot(x10, norm(pr['DESI_O2'][0], x10), 's', color=WARN, ms=3.5, mfc='none', label='DESI [O II]')
        ax.axvline(0, color=FGDIM, lw=0.6)
        ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.5, 1.35)
        ax.set_xlabel('offset from E spine (deg; + = east, away from Monogem)'); ax.set_ylabel('normalised')
        r = B[seg]
        txt = (f"Δ(FUV−[O III]) = {r['shift_FUV_vs_NSNS']['delta']:+.3f} ± {r['shift_FUV_vs_NSNS']['jk_err']:.3f}°\n"
               f"Δ(NUV−[O III]) = {r['shift_NUV_vs_NSNS']['delta']:+.3f} ± {r['shift_NUV_vs_NSNS']['jk_err']:.3f}°\n"
               f"[O III] σW/σE = {r['NSNS']['fit']['asym_w_over_e']:.2f} ± {r['NSNS']['jk_err']['asym_w_over_e']:.2f}")
        if 'shift_DESI_O2_vs_O3' in r and r['shift_DESI_O2_vs_O3']['jk_err'] is not None:
            txt += f"\nΔ(DESI [O II]−[O III]) = {r['shift_DESI_O2_vs_O3']['delta']:+.3f} ± {r['shift_DESI_O2_vs_O3']['jk_err']:.3f}°"
        label(ax, -0.97, 1.30, ('Dec 12–15 (single strand)' if j == 0 else 'Dec 12–17.6 (incl. W strand)'), size=8, va='top')
        label(ax, 0.25, 1.30, txt, size=7, va='top')
        ax.legend(loc='lower left', ncol=2, fontsize=6.8)
    # expectation panel
    ax = axes[1, 1]
    ax.axis('off')
    lines = [
        'Expected ordering across the filament (east = outside Monogem):',
        '',
        'Outward shock (consensus):  front | C IV/FUV | [O III] | [O II] → west',
        '   FUV east of [O III] by ~0.1° (shockmodel); [O II] slightly west;',
        '   sharp east edge, 1/√h tail to the west (σW/σE > 1).',
        '',
        'Hot-gas interface (conduction, mixing layer, NEI cooling of',
        '   Monogem gas; hot side west): FUV WEST of [O III], [O II] east.',
        '',
        'Static heated layer / photoionised sheet: no ordering, symmetric.',
        '',
        'Caveats: GALEX skybg maps are smoothed; the W strand adds west',
        'emission at Dec 15–17.6; DESI values are petal-relative.',
    ]
    for i, t in enumerate(lines):
        ax.text(0.0, 0.95 - i * 0.07, t, color=FG if i == 0 else FGDIM, fontsize=8.2, transform=ax.transAxes)
    finish(fig, OUT_FIG)


def main():
    A = test_A()
    B, store, x025, x10 = test_B()
    res = dict(test_A_registration=A, test_B_ordering=B,
               notes=dict(S15=S15, S176=S176, bins_deg=[0.025, 0.1],
                          method='see module docstring; jackknife over 10 along-spine chunks'))
    OUT_JSON.write_text(json.dumps(res, indent=1, default=float))
    figure(A, B, store, x025, x10)
    print('saved', OUT_JSON, OUT_FIG)


if __name__ == '__main__':
    main()
