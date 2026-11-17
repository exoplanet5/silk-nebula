"""Model-light view of the SDSS sky-level detection + figure.

For each SDSS run crossing the ridge at Dec 12-20.5 (bright-star fields excluded) we remove only the nuisance model
(S[f] common mode per field number + per-camcol offset and drift; no template) from (a) the measured sky_g and
(b) the NSNS [OIII] 'ridge' template sampled at the same fields.  Both residual sets are binned in distance from the
NSNS ridge; the template residuals are scaled by the fitted A and by the expected-if-calibrated A.
Outputs: data/derived/imaging_optical_sdss_profile.json, figures/imaging_optical_sdss_sky.png
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import gaussian_filter
from imaging_optical_common import RAW, DER, FIG, save_json, along_across, overlay_v1
from imaging_optical_nsns import ridge_spine
from imaging_optical_sdss_sky import load_fields, template_sampler, bright_star_fields, EXPECT
import silk_style as ss
import silk_geometry as sg


def nuisance_resid(f, c, y, ok):
    fu = np.unique(f[ok]); cu = np.unique(c[ok])
    fi = {v: i for i, v in enumerate(fu)}; ci = {v: i for i, v in enumerate(cu)}
    rows = np.where(ok)[0]; nf = len(fu); nc = len(cu)
    X = np.zeros((len(rows), nf + 2 * nc)); f0 = np.median(f[ok])
    for k, i in enumerate(rows):
        X[k, fi[f[i]]] = 1
        if ci[c[i]] > 0:
            X[k, nf + ci[c[i]]] = 1
        X[k, nf + nc + ci[c[i]]] = (f[i] - f0) / 100
    X = np.delete(X, nf, axis=1)
    out = np.full(len(y), np.nan)
    coef, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
    out[rows] = y[rows] - X @ coef
    return out


def main():
    a = load_fields()
    res = json.loads((DER / 'imaging_optical_sdss_sky.json').read_text())
    case = res['cases']['g_vs_OIII_ridge_nostars']
    core_runs = [r['run'] for r in case['runs'] if 12 <= r['dec_cross'] < 20.5]
    A_fit = case['core_dec12_20p5']['A']; eA = case['core_dec12_20p5']['eA']
    sra, sdec = ridge_spine()
    bad, _ = bright_star_fields(a)
    across, along, beyond = along_across(a['ra'], a['dec'], sra, sdec)
    samp = template_sampler('oiii', variant='ridge')
    T_all = samp(a['ra'], a['dec'])
    R_y, R_t, X_a, D_a, RA_a = [], [], [], [], []
    for run in core_runs:
        s = (a['run'] == run) & (np.abs(across) <= 2.5) & ~beyond & (a['quality'] >= 2) & ~bad
        idx = np.where(s)[0]
        f = a['field'][idx].astype(int); c = a['camcol'][idx].astype(int)
        y = a['sky_g'][idx]; T = T_all[idx]
        ok = np.isfinite(y) & np.isfinite(T) & (y > 0)
        ry = nuisance_resid(f, c, y, ok); rt = nuisance_resid(f, c, T, ok)
        # clip 4 sigma outliers of the data residual
        sres = 1.4826 * np.nanmedian(np.abs(ry - np.nanmedian(ry)))
        good = ok & (np.abs(ry - A_fit * rt) < 4 * sres)
        R_y.append(ry[good]); R_t.append(rt[good]); X_a.append(across[idx][good]); D_a.append(a['dec'][idx][good]); RA_a.append(a['ra'][idx][good])
    ry = np.concatenate(R_y); rt = np.concatenate(R_t); x = np.concatenate(X_a); dd = np.concatenate(D_a); rr = np.concatenate(RA_a)
    edges = np.arange(-2.5, 2.501, 0.2); cen = 0.5 * (edges[1:] + edges[:-1])
    py, ey, pt, nn = [], [], [], []
    for i in range(len(cen)):
        m = (x >= edges[i]) & (x < edges[i + 1])
        nn.append(int(m.sum()))
        if m.sum() >= 5:
            py.append(np.mean(ry[m])); ey.append(np.std(ry[m]) / np.sqrt(m.sum() / 2.0)); pt.append(np.mean(rt[m]))
        else:
            py.append(np.nan); ey.append(np.nan); pt.append(np.nan)
    py, ey, pt = map(np.array, (py, ey, pt))
    out = dict(runs=core_runs, n_fields=int(len(ry)), across_deg=cen.tolist(), sky_resid_nmgy=py.tolist(),
               sky_resid_err=ey.tolist(), template_resid_R=pt.tolist(), n_per_bin=nn, A_fit=A_fit, A_err=eA,
               A_expected=EXPECT['g'],
               note='errors assume pairs of fields share one sky estimate (N_eff = N/2)')
    save_json(out, DER / 'imaging_optical_sdss_profile.json')

    # ---------------------------------------------------------------- figure
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    w = WCS(str(z['wcs_header']))
    fig = plt.figure(figsize=(14.5, 7.2))
    fig.text(0.012, 0.985, 'Silk Nebula · SDSS sky levels 2002–2007', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    fig.text(0.012, 0.948, 'Field.sky_g of 13.5′×9.8′ drift-scan fields vs NSNS [OIII] template', color=ss.FGDIM, fontsize=9.5, va='top')
    ax = fig.add_axes([0.035, 0.07, 0.30, 0.80], projection=w)
    tr = ax.get_transform('world')
    ax.imshow(gaussian_filter(np.nan_to_num(z['oiii']), 1.5), origin='lower', cmap='bone', vmin=-0.3, vmax=3.0)
    allsel = np.isin(a['run'], core_runs) & (np.abs(across) <= 2.5) & ~beyond
    ax.scatter(a['ra'][allsel], a['dec'][allsel], s=1.2, color=ss.WARN, alpha=0.35, transform=tr, lw=0)
    for k, (ra0, de0) in sg.REGIONS.items():
        c = ss.WARN if k == 'P4' else '#ffffff'
        ax.plot([ra0], [de0], 'o', ms=5, mfc='none', mec=c, transform=tr)
        ss.label(ax, ra0 + 0.45, de0, k, color=c, size=7.5, transform=tr, va='center', ha='right')
    ax.set_xlim(*w.wcs_world2pix([121.5, 110.5], [16, 16], 0)[0])
    ax.set_ylim(*w.wcs_world2pix([116, 116], [9.0, 23.5], 0)[1])
    ax.coords[0].set_axislabel('RA'); ax.coords[1].set_axislabel('Dec')
    ax.coords[0].set_major_formatter('d'); ax.coords[1].set_major_formatter('d')
    ax.coords[0].set_ticks(spacing=3 * u.deg); ax.coords[1].set_ticks(spacing=2 * u.deg)
    ax.coords.grid(color='#2a323c', alpha=0.6, ls=':', lw=0.5)
    ax.set_title(f'{len(core_runs)} runs, fields used (orange)', color=ss.FG, fontsize=9)

    ax2 = fig.add_axes([0.42, 0.56, 0.52, 0.34])
    ok = np.isfinite(py)
    ax2.errorbar(cen[ok], 1e3 * py[ok], yerr=1e3 * ey[ok], fmt='o', ms=3.5, color=ss.OIII_C, ecolor=ss.OIII_C, elinewidth=0.8,
                 label='SDSS g sky residual')
    ax2.plot(cen[ok], 1e3 * A_fit * pt[ok], '-', color='#ffffff', lw=1.2, label=f'NSNS template × fit ({A_fit / EXPECT["g"]:.2f}±{eA / EXPECT["g"]:.2f})')
    ax2.plot(cen[ok], 1e3 * EXPECT['g'] * pt[ok], '--', color=ss.FGDIM, lw=1.0, label='NSNS template × 1.0')
    ax2.plot(cen[ok], 1e3 * 2.7 * EXPECT['g'] * pt[ok], ':', color=ss.DANGER, lw=1.2, label='× 2.7 (paper 10–13 R, min.)')
    ax2.axhline(0, color='#2a323c', lw=0.8)
    ax2.set_xlim(-2.5, 2.5)
    ax2.set_ylabel('Δ sky_g [10⁻³ nmgy/″²]')
    ax2.set_xlabel('distance from [OIII] ridge [deg]  (+ east)')
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=7, ncol=2, loc='upper left')
    ax2b = ax2.secondary_yaxis('right', functions=(lambda v: v / 1.7, lambda v: v * 1.7))
    ax2b.set_ylabel('% of sky (1.7 nmgy/″²)', color=ss.FGDIM)
    ax2b.tick_params(colors=ss.FGDIM)

    ax3 = fig.add_axes([0.42, 0.08, 0.52, 0.36])
    marks = {'g_vs_OIII_ridge_nostars': (ss.OIII_C, 'o', 'g vs [OIII]'), 'r_vs_OIII_ridge': (ss.HA_C, 's', 'r vs [OIII] (null)'),
             'g_vs_Halpha': (ss.FGDIM, '^', 'g vs Hα (null)')}
    for j, (key, (col, mk, lab)) in enumerate(marks.items()):
        cs = res['cases'][key]; ex = EXPECT[cs['band']]
        d = np.array([r['dec_cross'] for r in cs['runs']]); A = np.array([r['A'] for r in cs['runs']]) / ex
        e = np.array([r['eA'] for r in cs['runs']]) / ex
        m = (d >= 12) & (d < 20.5)
        cc = cs['core_dec12_20p5']
        ax3.errorbar(d[m] + 0.06 * (j - 1), A[m], yerr=e[m], fmt=mk, ms=4, color=col, elinewidth=0.7, alpha=0.9,
                     label=f"{lab}: {cc['A_over_expected']:.2f}±{cc['err_over_expected']:.2f}")
    ax3.axhline(0, color='#2a323c', lw=0.8); ax3.axhline(1, color=ss.FGDIM, lw=0.8, ls='--')
    ax3.axhspan(2.7, 6.2, color=ss.DANGER, alpha=0.08)
    ss.label(ax3, 12.2, 5.6, 'paper 10–13 R', size=7, color=ss.DANGER, va='center')
    ax3.set_xlim(12, 20.6); ax3.set_ylim(-3, 7)
    ax3.set_xlabel('Dec of ridge crossing [deg]'); ax3.set_ylabel('A / A(NSNS calibration)')
    ax3.grid(alpha=0.3); ax3.legend(fontsize=7, loc='upper right', ncol=1)
    ss.credit(fig)
    fig.savefig(FIG / 'imaging_optical_sdss_sky.png', dpi=200)
    plt.close(fig)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}, indent=1))
    print('peak bins', [(round(c, 1), round(1e3 * p, 2), round(1e3 * e, 2), round(t, 2)) for c, p, e, t in zip(cen, py, ey, pt) if abs(c) < 0.8])


if __name__ == '__main__':
    main()
