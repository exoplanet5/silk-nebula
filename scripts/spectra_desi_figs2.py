"""spectra_desi_figs2.py -- final figures for the DESI DR1 residual-emission analysis of the Silk Nebula.

Usage: python -I scripts/spectra_desi_figs2.py [map|hamap|profile|regspec|velocity|stacks|ridge|valid|coverage|all]
Inputs: data/derived/spectra_desi_fiber_table.csv, spectra_desi_results.json, spectra_desi_model_arrays.npz,
        spectra_desi_fits.npz, spectra_desi_ratios.json, spectra_desi_regspec_full.npz, spectra_desi_rvcheck.json,
        spectra_desi_vcorr.npz, silk_trace.json; teammate imaging_optical_nsns_measurements.json (NSNS R values).
"""
import csv
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, FG, FGDIM, OIII_C, HA_C, SII_C, NII_C, HB_C, WARN,  # noqa
                        SUN, OK, ACCENT, BG2, DANGER)
from silk_geometry import ROOT, REGIONS, load_trace, V_HELIO  # noqa: E402
from spectra_desi_fig_coverage import draw_trace, draw_regions  # noqa: E402
from spectra_desi_figs import fibtab, sky_axes, binned_map  # noqa: E402
from spectra_desi_lines import LINES, WIN_WAVE, C_KMS  # noqa: E402
from spectra_desi_analyze import rayleigh_per_unit  # noqa: E402
import spectra_desi_measure as meas  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

D = ROOT / 'data/derived'
FIG = ROOT / 'figures'
XL, YL = (119.6, 113.4), (11.8, 28.6)
LAMOST_OIII = {'P1': (26.8, 4.3), 'P2': (26.8, 4.1), 'P3': (23.8, 3.3)}      # Wu+2026 Table 4, V_helio
def nf(w, h, title, subtitle, **kw):
    """new_sky_figure with the subtitle placed for any figure height (avoids title/subtitle overlap)."""
    fig, axes = new_sky_figure(w, h, title, '', **kw)
    fig.text(0.012, 0.985 - 0.34 / h, subtitle, color=FGDIM, fontsize=9.5, va='top')
    fig._silk_top = 1 - 0.66 / h
    return fig, axes


def fin(fig, path):
    fig.tight_layout(rect=(0, 0.025, 1, fig._silk_top))
    return finish(fig, path, tight=False)


LAMOST_HA = {'P1': (16.3, 1.2), 'P2': (17.0, 1.2), 'P3': (18.4, 1.4), 'P4': (19.7, 1.4)}


def patch_outline(ax, t, color=HA_C):
    p = np.array(t['halpha_patch'])
    ax.plot(np.r_[p[:, 0], p[0, 0]], np.r_[p[:, 1], p[0, 1]], '-', color=color, lw=0.8, alpha=0.9)


def tall_map(values_list, titles, fname, title, subtitle, cbl, vlim=2.0, binsz=0.2, extra=None, cmap='RdBu_r'):
    t = load_trace()
    n = len(values_list)
    fig, axes = new_sky_figure(3.3 * n + 1.6, 10.6, title, subtitle, ncols=n, sharey=True)
    axes = np.atleast_1d(axes)
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
    dec0 = np.mean(YL)
    xe = np.arange(XL[1], XL[0] + 1e-6, binsz / np.cos(np.radians(20)))
    ye = np.arange(YL[0], YL[1] + 1e-6, binsz)
    im = None
    for a, (kind, ra, dec, v, e), ttl in zip(axes, values_list, titles):
        if kind == 'scatter':
            o = np.argsort(np.abs(v))
            im = a.scatter(ra[o], dec[o], c=v[o], s=3.5, cmap=cmap, norm=norm, lw=0)
        else:
            m, _ = binned_map(ra, dec, v, e, xe, ye)
            im = a.imshow(m, origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]), cmap=cmap, norm=norm,
                          aspect='auto', interpolation='nearest')
        draw_trace(a, t, lw=0.7, alpha=0.55, envelope=False, color='#d8f3ff')
        if extra:
            extra(a, t)
        for k, (rr, dd) in REGIONS.items():
            a.plot(rr, dd, 'x' if k != 'P4' else '+', color=SUN if k != 'P4' else FGDIM, ms=5, mew=1.2)
        sky_axes(a, XL, YL, dec0)
        a.set_title(ttl, color=FG, fontsize=9.5)
        if a is not axes[0]:
            a.set_ylabel('')
    for k, (rr, dd) in REGIONS.items():
        label(axes[0], rr - 0.15, dd, k, color=SUN if k != 'P4' else FGDIM, size=7.5, va='center')
    fig.subplots_adjust(left=0.06, right=0.9, top=0.9, bottom=0.06, wspace=0.08)
    cax = fig.add_axes([0.915, 0.3, 0.012, 0.4])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbl, color=FG)
    finish(fig, FIG / fname, tight=False)


def fig_map():
    t = fibtab()
    ok = t['ok'] & np.isfinite(t['OIII_R'])
    ra, dec, v, e = t['ra'][ok], t['dec'][ok], t['OIII_R'][ok], t['OIII_R_err'][ok]
    p = t['model_resid_R'][ok]
    mp = np.isfinite(p)
    vals = [('scatter', ra, dec, v, e), ('bin', ra, dec, v, e),
            ('bin', ra[mp], dec[mp], p[mp], np.ones(mp.sum())),
            ('bin', ra[mp], dec[mp], (v - p)[mp], e[mp])]
    tall_map(vals, ['per fibre', 'binned 0.2°', 'ridge model', 'data − model'], 'spectra_desi_oiii_map.png',
             'DESI DR1 · [O III] 5007 residual', f'{ok.sum():,} galaxy/QSO fibres · fibre minus petal-mean sky',
             '[O III] residual (R)', vlim=2.0)


def fig_hamap():
    t = fibtab()
    f = np.load(D / 'spectra_desi_fits.npz', allow_pickle=False)
    ids = list(f['sparcl_id'])
    pos = {s: i for i, s in enumerate(ids)}
    k = np.array([pos[s] for s in t['sparcl_id']])
    ff = t['ffrac']
    vals, titles = [], []
    for g, j, lam, ttl in (('HaNII', 0, LINES['Ha'], 'Hα'), ('HaNII', 1, LINES['NII6583'], '[N II] 6583'),
                           ('SII', 0, LINES['SII6716'], '[S II] 6716+6731')):
        a = f[f'{g}__amp_fix'][k, j]
        e = f[f'{g}__err_fix'][k, j]
        if g == 'SII':
            a = a + f[f'{g}__amp_fix'][k, 1]
            e = np.hypot(e, f[f'{g}__err_fix'][k, 1])
        okk = np.isfinite(a) & np.isfinite(e) & (e > 0) & ~f[f'{g}__flag_own'][k] & (f[f'{g}__npix_good'][k] > 40)
        kr = rayleigh_per_unit(lam)
        vals.append(('bin', t['ra'][okk], t['dec'][okk], (a * ff * kr)[okk], (e * ff * kr)[okk]))
        titles.append(ttl)
    tall_map(vals, titles, 'spectra_desi_halpha_map.png', 'DESI DR1 · red-line residuals',
             'binned 0.2° · fibre minus petal-mean sky · red outline = NSNS Hα patch', 'residual (R)', vlim=1.5,
             extra=patch_outline)


def fig_profile():
    t = fibtab()
    ok = t['ok'] & np.isfinite(t['OIII_R']) & np.isfinite(t['dRA_spine'])
    bands = [(12.0, 15.0), (15.0, 16.0), (16.0, 16.8), (16.8, 17.6), (17.6, 19.0), (19.0, 21.0), (21.0, 24.0),
             (24.0, 28.6)]
    bands = [b for b in bands if (ok & (t['dec'] >= b[0]) & (t['dec'] < b[1])).sum() >= 40]
    nb = len(bands)
    ncol = 2
    nrow = int(np.ceil(nb / ncol))
    fig, axes = new_sky_figure(13.0, 2.3 * nrow + 1.3, 'DESI [O III] residual across the filament',
                               'ivar-weighted means in 0.06° bins of ΔRA from the E spine · '
                               'line = ridge model after petal subtraction', nrows=nrow, ncols=ncol, squeeze=False,
                               sharex=True)
    for (lo, hi), a in zip(bands, axes.ravel()):
        m = ok & (t['dec'] >= lo) & (t['dec'] < hi)
        x, y, e, p = t['dRA_spine'][m], t['OIII_R'][m], t['OIII_R_err'][m], t['model_resid_R'][m]
        be = np.arange(-3.0, 3.0001, 0.06)
        idx = np.digitize(x, be) - 1
        xc, ym, em, pm = [], [], [], []
        for kk in range(len(be) - 1):
            s = idx == kk
            if s.sum() >= 3:
                w = 1 / e[s] ** 2
                xc.append(0.5 * (be[kk] + be[kk + 1]))
                ym.append(np.sum(w * y[s]) / w.sum())
                em.append(w.sum() ** -0.5)
                pp = p[s][np.isfinite(p[s])]
                pm.append(pp.mean() if pp.size else np.nan)
        a.axhline(0, color=FGDIM, lw=0.6)
        a.axvline(0, color=OIII_C, lw=0.5, alpha=0.5, ls=':')
        a.errorbar(xc, ym, em, fmt='o', ms=2.6, color=OIII_C, ecolor=OIII_C, elinewidth=0.7, label='DESI')
        xc_, pm_ = np.array(xc, float), np.array(pm, float)
        gap = np.r_[False, np.diff(xc_) > 0.13]
        xc_ = np.insert(xc_, np.where(gap)[0], np.nan)
        pm_ = np.insert(pm_, np.where(gap)[0], np.nan)
        a.plot(xc_, pm_, '-', color=SUN, lw=1.3, label='ridge model')
        a.set_xlim(3.0, -3.0)
        a.set_ylim(-2.2, 2.8)
        label(a, 0.02, 0.86, f'Dec {lo:.1f}–{hi:.1f}  n={m.sum()}', transform=a.transAxes, size=8)
        a.grid(alpha=0.2)
    for a in axes.ravel()[nb:]:
        a.set_visible(False)
    for a in axes[:, 0]:
        a.set_ylabel('residual (R)')
    for a in axes[-1]:
        a.set_xlabel('ΔRA·cosδ from E spine (deg; east → left)')
    axes[0, 0].legend(loc='upper right', fontsize=7, ncol=2)
    finish(fig, FIG / 'spectra_desi_profile.png', tight=False)


PANELS = [('[O II] 3726,3729', 3712, 3745, ['OII3726', 'OII3729'], 'OII3727'),
          ('[Ne III] 3869', 3855, 3885, ['NeIII3869'], 'NeIII3869'),
          ('[O III] 4363', 4350, 4378, ['OIII4363'], 'OIII4363'),
          ('He II 4686', 4672, 4702, ['HeII4686'], 'HeII4686'),
          ('Hβ', 4848, 4878, ['Hb'], 'Hb'),
          ('[O III] 4959,5007', 4948, 5022, ['OIII4959', 'OIII5007'], 'OIII5007'),
          ('[N I] 5198,5200', 5186, 5214, ['NI5198', 'NI5200'], 'NI5200'),
          ('O I 5577 (sky)', 5564, 5594, ['sky5577'], 'sky5577'),
          ('[O I] 6300', 6288, 6316, ['OI6300'], 'OI6300'),
          ('Hα + [N II]', 6538, 6600, ['NII6548', 'Ha', 'NII6583'], 'HaNII'),
          ('[S II] 6716,6731', 6705, 6745, ['SII6716', 'SII6731'], 'SII'),
          ('[S III] 9069', 9055, 9085, ['SIII9069'], 'SIII9069')]


def model_curve(fit, g, lam):
    from spectra_desi_regspec import GROUPS_1D
    cols = GROUPS_1D[g]
    y = np.zeros_like(lam)
    for amp, c in zip(fit['amp'], cols):
        for l0, rel in c:
            s = np.hypot(meas.lsf_sigma_A(l0), l0 * 20.0 / C_KMS)
            mu = l0 * (1 + fit['v'] / C_KMS)
            y += amp * rel * np.exp(-0.5 * ((lam - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
    return y


def fig_regspec(run='model_all'):
    r = np.load(D / 'spectra_desi_regspec_full.npz', allow_pickle=False)
    J = json.load(open(D / 'spectra_desi_ratios.json'))[run]
    S, eS, lo_, hi_ = r[f'S_{run}'], r[f'eS_{run}'], r[f'p16_{run}'], r[f'p84_{run}']
    fig, axes = nf(14.0, 8.4, 'DESI regression spectrum of the [O III] filament',
                               f'per-pixel fit c = a + S·P, P = predicted [O III] residual · {J["n_fibres"]:,} fibres · '
                               f'band = healpix bootstrap 16–84% · per unit [O III] 5007', nrows=3, ncols=4)
    for (ttl, lo, hi, lines, g), a in zip(PANELS, axes.ravel()):
        s = (WIN_WAVE >= lo) & (WIN_WAVE <= hi)
        a.fill_between(WIN_WAVE[s], lo_[s], hi_[s], step='mid', color=OIII_C, alpha=0.25, lw=0)
        a.step(WIN_WAVE[s], S[s], where='mid', color=OIII_C, lw=0.9)
        fx = J['fixed'].get(f'{g}@vO')
        if fx:
            lam = np.linspace(lo, hi, 400)
            a.plot(lam, model_curve(fx, g, lam), '-', color=SUN, lw=1.0, alpha=0.9)
        a.axhline(0, color=FGDIM, lw=0.5)
        for ln in lines:
            a.axvline(LINES[ln] * (1 + J['v_ref_free']['v'] / C_KMS), color=FGDIM, lw=0.5, alpha=0.6, ls=':')
        a.set_title(ttl, color=FG, fontsize=9)
        a.set_xlim(lo, hi)
        a.tick_params(labelsize=7)
        a.grid(alpha=0.15)
        if g != 'OIII5007':
            yy = np.r_[lo_[s], hi_[s]]
            yy = yy[np.isfinite(yy)]
            if yy.size:
                lim = max(np.percentile(np.abs(yy), 98) * 1.3, 0.05)
                a.set_ylim(-lim, lim * 1.4)
        if fx:
            amp = fx['amp'][0] if g not in ('SII', 'OII3727') else fx['amp'][0] + fx['amp'][1]
            err = fx['boot_err'][0] if g not in ('SII', 'OII3727') else np.hypot(*fx['boot_err'][:2])
            label(a, 0.03, 0.86, f'{amp:.3f}±{err:.3f}', transform=a.transAxes, size=7.5, color=SUN)
    for a in axes[:, 0]:
        a.set_ylabel('R Å⁻¹ per R [O III]')
    for a in axes[-1]:
        a.set_xlabel('vacuum λ (Å, barycentric)')
    fin(fig, FIG / 'spectra_desi_regspec.png')


def fig_velocity():
    J = json.load(open(D / 'spectra_desi_ratios.json'))
    t = fibtab()
    vc = np.load(D / 'spectra_desi_vcorr.npz', allow_pickle=False)
    vcm = dict(zip(vc['sparcl_id'], vc['vcorr']))
    res = json.load(open(D / 'spectra_desi_results.json'))
    rv = json.load(open(D / 'spectra_desi_rvcheck.json'))
    fig, axes = nf(13.5, 5.6, 'DESI [O III] velocity along the Silk Nebula',
                               'heliocentric (DESI barycentric) · DESI regression bands vs LAMOST MRS-N (Wu+2026)',
                               ncols=2, gridspec_kw=dict(width_ratios=[1.6, 1]))
    a = axes[0]
    sn = np.abs(t['OIII_fix_1e17'] / t['OIII_fix_err'])
    m = t['ok'] & (sn > 5) & np.isfinite(t['OIII_v_err']) & (t['OIII_v_err'] < 30) & (t['d_any'] < 1.2)
    a.errorbar(t['dec'][m], t['OIII_v_free'][m], t['OIII_v_err'][m], fmt='.', ms=3, color=FGDIM, alpha=0.45,
               elinewidth=0.4, label='DESI fibre |S/N|>5')
    for b in J['v_by_dec']:
        x = 0.5 * (b['dec_lo'] + b['dec_hi'])
        lo, hi = max(b['v_boot_p16'], b['v_lo'] - 15), min(b['v_boot_p84'], b['v_hi'] + 15)
        a.errorbar([x], [b['v']], [[max(b['v'] - lo, 0.3)], [max(hi - b['v'], 0.3)]],
                   xerr=[[x - b['dec_lo']], [b['dec_hi'] - x]], fmt='o', ms=6, color=OIII_C, ecolor=OIII_C,
                   elinewidth=1.3, capsize=2, zorder=5)
    a.plot([], [], 'o', color=OIII_C, label='DESI regression (Dec band)')
    for k, (v, e) in LAMOST_OIII.items():
        d = REGIONS[k][1]
        a.errorbar([d], [v], [e], fmt='s', ms=6, color=SUN, ecolor=SUN, capsize=3, zorder=6)
        label(a, d + 0.15, v + e + 0.8, k, color=SUN, size=8)
    a.plot([], [], 's', color=SUN, label='LAMOST [O III]')
    a.axhspan(11.7, 19.7, color=HA_C, alpha=0.12, lw=0)
    a.plot([], [], 's', color=HA_C, alpha=0.4, ms=8, label='LAMOST Hα/[N II]/[S II] range')
    vO = J['model_all']['v_ref_free']['v']
    a.axhline(vO, color=OIII_C, lw=0.8, ls='--', alpha=0.7)
    a.axhspan(-0.5 + 12.1, 1.5 + 12.1, color=OK, alpha=0.18, lw=0)
    a.plot([], [], 's', color=OK, alpha=0.5, ms=8, label='WHAM Hα (V_LSR −0.5…+1.5)')
    if 'nbrHa_all' in J:
        vh = J['nbrHa_all']['v_ref_boot']
        a.errorbar([27.6], [J['nbrHa_all']['v_ref_free']['v']],
                   [[J['nbrHa_all']['v_ref_free']['v'] - vh['p16']], [vh['p84'] - J['nbrHa_all']['v_ref_free']['v']]],
                   fmt='^', ms=7, color=HA_C, ecolor=HA_C, capsize=3, zorder=6, label='DESI Hα structure (all Dec)')
    a.set_xlim(11.5, 28.8)
    a.set_ylim(-30, 60)
    a.set_xlabel('Dec (deg)')
    a.set_ylabel('V_helio (km/s)')
    a.grid(alpha=0.2)
    a.legend(loc='upper left', fontsize=6.8, ncol=3)
    # right: v vs vcorr + RV check
    a = axes[1]
    vcorr = np.array([vcm.get(s, np.nan) for s in t['sparcl_id']])
    mm = m & np.isfinite(vcorr)
    a.errorbar(vcorr[mm], t['OIII_v_free'][mm], t['OIII_v_err'][mm], fmt='.', ms=3, color=OIII_C, alpha=0.5,
               elinewidth=0.4)
    q = res.get('OIII_v_vs_vcorr')
    if q:
        xx = np.linspace(-35, 35, 10)
        a.plot(xx, q['intercept'] + q['slope'] * xx, '-', color=SUN, lw=1.2,
               label=f"slope {q['slope']:+.2f}±{q['err_slope']:.2f}")
        a.plot(xx, q['intercept'] - xx, ':', color=DANGER, lw=1.0, label='no bary. corr. (slope −1)')
    tw = rv.get('three_way', {})
    txt = []
    for key, lab in (('DESI_main/dark_minus_LRS', 'DESI dark−LRS'), ('DESI_main/bright_minus_LRS', 'DESI bright−LRS'),
                     ('LRS_minus_Gaia', 'LRS−Gaia'), ('DESI_main/backup_minus_Gaia', 'DESI backup−Gaia')):
        if key in tw:
            txt.append(f"{lab:16s} {tw[key]['median']:+6.1f}±{tw[key]['err_median']:.1f} (n={tw[key]['n']})")
    label(a, 0.03, 0.04, 'stellar RV checks (km/s)\n' + '\n'.join(txt), transform=a.transAxes, size=6.8,
          color=FG, va='bottom', family='Menlo')
    a.set_xlim(-35, 35)
    a.set_ylim(-30, 60)
    a.set_xlabel('barycentric correction of exposures (km/s)')
    a.set_ylabel('[O III] V_helio (km/s)')
    a.grid(alpha=0.2)
    a.legend(loc='upper right', fontsize=7)
    fin(fig, FIG / 'spectra_desi_velocity.png')


def fig_stacks():
    r = np.load(D / 'spectra_desi_regspec_full.npz', allow_pickle=False)
    J = json.load(open(D / 'spectra_desi_ratios.json'))
    wins = [('[O III] 4959,5007', 4945, 5025), ('Hβ', 4845, 4880), ('Hα + [N II]', 6535, 6600),
            ('[S II]', 6705, 6745)]
    fig, axes = nf(13.5, 4.4, 'DESI on/off-filament stacks',
                               'ivar-weighted mean residual (R Å⁻¹, in fibre) · on = top 10% predicted [O III], '
                               'neg = bottom 10%, ctrl = >5° from trace', ncols=4)
    cols = {'on': OIII_C, 'neg': WARN, 'ctrl': FGDIM, 'patch': HA_C}
    for (ttl, lo, hi), a in zip(wins, axes):
        s = (WIN_WAVE >= lo) & (WIN_WAVE <= hi)
        for lab in ('ctrl', 'neg', 'on', 'patch'):
            if f'stack_{lab}' not in r.files:
                continue
            n = J['stacks'][lab]['n']
            a.step(WIN_WAVE[s], r[f'stack_{lab}'][s], where='mid', color=cols[lab], lw=0.9,
                   label=f'{lab} (n={n})')
        a.axhline(0, color=FGDIM, lw=0.5)
        a.set_title(ttl, color=FG, fontsize=9)
        a.set_xlim(lo, hi)
        a.grid(alpha=0.15)
        a.set_xlabel('vacuum λ (Å)')
    axes[0].set_ylabel('R Å⁻¹')
    axes[0].legend(loc='upper left', fontsize=6.8)
    fin(fig, FIG / 'spectra_desi_stacks.png')


def fig_ridge():
    res = json.load(open(D / 'spectra_desi_results.json'))
    fm = res['forward_model']
    nsns = json.load(open(D / 'imaging_optical_nsns_measurements.json'))
    fig, a = nf(10.5, 4.8, 'Silk Nebula [O III] ridge brightness',
                            'DESI forward-model ridge peak (absolute, flux-calibrated) vs NSNS imaging and Wu+2026')
    for strand, c, mk in (('E', OIII_C, 'o'), ('W', ACCENT, 's'), ('W2', NII_C, 'D')):
        xs, ys, es = [], [], []
        for nm, amp, err in zip(fm['basis'], fm['amp_R'], fm['err_R']):
            st, d = nm.split('@')
            if st == strand and err < 0.8:
                xs.append(float(d))
                ys.append(amp)
                es.append(err)
        if xs:
            a.errorbar(xs, ys, es, fmt=mk, ms=5, color=c, ecolor=c, capsize=2, label=f'DESI ridge {strand}',
                       mfc=c if strand == 'E' else 'none')
    for k, v in nsns.items():
        if not isinstance(v, dict) or 'dec_center' not in v or 'oiii' not in v:
            continue
        o = v['oiii']
        a.plot([v['dec_center'] - v['half_length_deg'], v['dec_center'] + v['half_length_deg']],
               [o['peak_excess_R']] * 2, '-', color=OK, lw=2.2, alpha=0.8)
    a.plot([], [], '-', color=OK, lw=2.2, label='NSNS peak excess (imaging teammate)')
    for k in ('P1', 'P2', 'P3'):
        d = REGIONS[k][1]
        a.errorbar([d], [11.5], [[1.6], [1.4]], fmt='D', color=SUN, ms=7, capsize=2)
        label(a, d + 0.12, 12.6, k, color=SUN, size=8)
    a.plot([], [], 'D', color=SUN, label='Wu+2026 Table 3 (~10–13 R)')
    a.axhline(0, color=FGDIM, lw=0.5)
    a.set_xlim(11.5, 28.8)
    a.set_ylim(-1.5, 15)
    a.set_xlabel('Dec (deg)')
    a.set_ylabel('[O III] 5007 (R)')
    a.grid(alpha=0.2)
    a.legend(loc='upper right', fontsize=7.5, ncol=2)
    fin(fig, FIG / 'spectra_desi_ridge_sb.png')


def fig_valid():
    f = np.load(D / 'spectra_desi_fits.npz', allow_pickle=False)
    t = fibtab()
    cat = {r['sparcl_id']: r for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv'))}
    dany = np.array([float(cat[s]['d_any']) for s in f['sparcl_id']])
    arr = np.load(D / 'spectra_desi_model_arrays.npz', allow_pickle=False)
    fig, axes = nf(13.5, 4.6, 'DESI residual validation',
                               'line-free null windows, control fibres >5° from the trace, airglow 5577',
                               ncols=3)
    a = axes[0]
    bins = np.linspace(-10, 10, 81)

    def chi(g, m=None):
        ok = np.isfinite(f[f'{g}__amp_fix'][:, 0]) & (f[f'{g}__err_fix'][:, 0] > 0) & ~f[f'{g}__flag_own']
        if m is not None:
            ok &= m
        return f[f'{g}__amp_fix'][ok, 0] / f[f'{g}__err_fix'][ok, 0]

    nul = np.concatenate([chi(g) for g in ('null4990', 'null5030', 'null4880')])
    a.hist(nul, bins, histtype='step', color=FGDIM, density=True, label=f'null windows (σ={np.std(nul):.2f})')
    co = chi('OIII5007', dany > 5)
    if co.size > 20:
        a.hist(co, bins, histtype='step', color=OK, density=True, label=f'[O III] ctrl >5° (σ={np.std(co):.2f})')
    so = chi('OIII5007', dany < 1)
    a.hist(so, bins, histtype='step', color=OIII_C, density=True, label=f'[O III] <1° (σ={np.std(so):.2f})')
    xx = np.linspace(-10, 10, 200)
    a.plot(xx, np.exp(-0.5 * xx ** 2) / np.sqrt(2 * np.pi), ':', color=SUN, lw=0.8)
    a.set_yscale('log')
    a.set_ylim(1e-4, 1)
    a.set_xlabel('amplitude / error')
    a.set_title('per-fibre significance', color=FG, fontsize=9)
    a.legend(fontsize=6.8, loc='upper left')
    # [O III] residual vs predicted (binned)
    a = axes[1]
    P, OR, eOR, ok = arr['pred'], arr['OR'], arr['eOR'], arr['okO'] & np.isfinite(arr['pred'])
    be = np.linspace(-2, 3, 21)
    idx = np.digitize(P[ok], be) - 1
    xc, ym, em = [], [], []
    for k in range(len(be) - 1):
        s = idx == k
        if s.sum() >= 5:
            w = 1 / eOR[ok][s] ** 2
            xc.append(np.average(P[ok][s], weights=w))
            ym.append(np.sum(w * OR[ok][s]) / w.sum())
            em.append(w.sum() ** -0.5)
    a.errorbar(xc, ym, em, fmt='o', ms=4, color=OIII_C, ecolor=OIII_C, label='[O III] 5007')
    g5577 = f['sky5577__amp_fix'][:, 0] * np.clip(t['ffrac'], 0.3, 1) * rayleigh_per_unit(5578.9)
    e5577 = f['sky5577__err_fix'][:, 0] * np.clip(t['ffrac'], 0.3, 1) * rayleigh_per_unit(5578.9)
    ok2 = ok & np.isfinite(g5577) & (e5577 > 0)
    idx = np.digitize(P[ok2], be) - 1
    xc2, ym2, em2 = [], [], []
    for k in range(len(be) - 1):
        s = idx == k
        if s.sum() >= 5:
            v = g5577[ok2][s]
            xc2.append(np.mean(P[ok2][s]))
            ym2.append(np.median(v))
            em2.append(1.2533 * 1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(s.sum()))
    a.errorbar(xc2, ym2, em2, fmt='s', ms=3.5, color=WARN, ecolor=WARN, label='5577 airglow (median)')
    a.plot([-2, 3], [-2, 3], ':', color=FGDIM, lw=0.8)
    a.axhline(0, color=FGDIM, lw=0.5)
    a.set_xlabel('predicted [O III] residual (R)')
    a.set_ylabel('measured residual (R)')
    a.set_title('measured vs ridge-model prediction', color=FG, fontsize=9)
    a.legend(fontsize=7, loc='upper left')
    a.grid(alpha=0.2)
    # residual (data - model) chi vs distance from trace
    a = axes[2]
    r = (OR - P) / eOR
    m = ok & np.isfinite(r)
    be = np.arange(0, 6.01, 0.25)
    idx = np.digitize(dany[m], be) - 1
    xc, s1, s2 = [], [], []
    for k in range(len(be) - 1):
        s = idx == k
        if s.sum() >= 15:
            xc.append(0.5 * (be[k] + be[k + 1]))
            s1.append(np.std(OR[m][s] / eOR[m][s]))
            s2.append(np.std(r[m][s]))
    a.plot(xc, s1, 'o-', color=OIII_C, ms=3.5, lw=1, label='[O III] data / error')
    a.plot(xc, s2, 's-', color=SUN, ms=3.5, lw=1, label='(data − model) / error')
    a.axhline(1, color=FGDIM, lw=0.6, ls=':')
    a.set_xlabel('distance to nearest strand (deg)')
    a.set_ylabel('rms significance')
    a.set_title('excess scatter vs distance', color=FG, fontsize=9)
    a.legend(fontsize=7)
    a.grid(alpha=0.2)
    fin(fig, FIG / 'spectra_desi_validation.png')


def fig_coverage():
    import subprocess  # noqa: F401  (not used; coverage uses its own script)
    t = load_trace()
    tb = fibtab()
    ra, dec = [], []
    for r in csv.DictReader(open(D / 'spectra_desi_catalog.csv')):
        ra.append(float(r['ra']))
        dec.append(float(r['dec']))
    ra, dec = np.array(ra), np.array(dec)
    fig, ax = new_sky_figure(8.2, 11.0, 'DESI DR1 spectra · Silk Nebula field',
                             f'{len(ra):,} DR1 coadds (density) · {tb["ok"].sum():,} analysed galaxy/QSO fibres')
    H, xe, ye = np.histogram2d(ra, dec, bins=[np.arange(111, 122.01, 0.2), np.arange(6, 28.01, 0.2)])
    dens = H.T / (0.2 * 0.2 * np.cos(np.radians(0.5 * (ye[:-1] + ye[1:])))[:, None])
    im = ax.imshow(np.where(dens > 0, dens, np.nan), origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]),
                   cmap='magma', vmin=0, vmax=4000, aspect='auto', interpolation='nearest', alpha=0.85)
    draw_trace(ax, t, lw=1.4)
    draw_regions(ax)
    ax.plot(tb['ra'][tb['ok']], tb['dec'][tb['ok']], '.', ms=1.0, color=OK, alpha=0.8, label='analysed fibres')
    ax.set_xlim(122, 111)
    ax.set_ylim(6, 28)
    ax.set_aspect(1 / np.cos(np.radians(17)))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.grid(alpha=0.25, lw=0.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('DR1 spectra per deg²', color=FG)
    ax.legend(loc='lower left', ncol=2, fontsize=7, markerscale=6)
    finish(fig, FIG / 'spectra_desi_coverage.png')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    fns = dict(map=fig_map, hamap=fig_hamap, profile=fig_profile, regspec=fig_regspec, velocity=fig_velocity,
               stacks=fig_stacks, ridge=fig_ridge, valid=fig_valid, coverage=fig_coverage)
    for k, fn in fns.items():
        if what in (k, 'all'):
            fn()
            print('made', k)


def fig_ratios():
    """Forest plot of the filament-component line ratios (energy units, healpix-bootstrap intervals)."""
    J = json.load(open(D / 'spectra_desi_ratios.json'))
    pc = json.load(open(D / 'paper_checks.json'))['inputs']['table2']
    fig, axes = nf(13.5, 6.4, 'Silk Nebula · DESI differential line ratios',
                   'component spatially correlated with [O III] (left) or with Hα structure (right) · '
                   'bars: bootstrap 68% / 95%', ncols=2)
    runs = [('model_all', OIII_C, 'all fibres'), ('nbr_all', '#cfd6dd', 'neighbour instrument'),
            ('model_south', OK, 'Dec < 19'), ('model_north', NII_C, 'Dec ≥ 22')]
    left = [('OII3727/OIII5007', '[O II] 3727'), ('NeIII3869/OIII5007', '[Ne III] 3869'),
            ('OIII4363/OIII5007', '[O III] 4363'), ('HeII4686/OIII5007', 'He II 4686'),
            ('Hb/OIII5007', 'Hβ'), ('NI5200/OIII5007', '[N I] 5200'), ('HeI5876/OIII5007', 'He I 5876'),
            ('OI6300/OIII5007', '[O I] 6300'), ('Ha/OIII5007', 'Hα'), ('NII6583/OIII5007', '[N II] 6583'),
            ('SII6716+31/OIII5007', '[S II] 6716+31'), ('OII7325/OIII5007', '[O II] 7325'),
            ('SIII9069/OIII5007', '[S III] 9069'), ('sky5577/OIII5007', 'O I 5577 sky')]
    a = axes[0]
    for k, (rk, lab) in enumerate(left):
        for j, (run, c, rl) in enumerate(runs):
            if run not in J or rk not in J[run]['ratios']:
                continue
            q = J[run]['ratios'][rk]
            y = k + 0.18 * (j - 1.5)
            a.plot([q['p2p5'], q['p97p5']], [y, y], '-', color=c, lw=0.8, alpha=0.6)
            a.plot([q['p16'], q['p84']], [y, y], '-', color=c, lw=2.6)
            a.plot([q['energy']], [y], 'o', color=c, ms=3.5, mec='none')
    for j, (run, c, rl) in enumerate(runs):
        if run in J:
            a.plot([], [], '-', color=c, lw=2.6, label=f'{rl} (n={J[run]["n_fibres"]})')
    a.set_yticks(range(len(left)))
    a.set_yticklabels([l for _, l in left])
    a.invert_yaxis()
    a.set_xscale('symlog', linthresh=0.03, linscale=0.5)
    a.set_xlim(-0.4, 10)
    a.set_xticks([-0.3, -0.1, 0, 0.03, 0.1, 0.3, 1, 3, 10])
    a.set_xticklabels(['−0.3', '−0.1', '0', '0.03', '0.1', '0.3', '1', '3', '10'])
    a.axvline(0, color=FGDIM, lw=0.6)
    a.axvline(1, color=FGDIM, lw=0.5, ls=':')
    a.set_xlabel('line / [O III] 5007 (energy)')
    a.grid(alpha=0.15, axis='x')
    a.legend(loc='upper left', fontsize=6.6, framealpha=0.92)
    a = axes[1]
    right = [('NII6583/Ha', '[N II] 6583 / Hα'), ('SII6716+31/Ha', '[S II] / Hα'), ('OI6300/Ha', '[O I] 6300 / Hα'),
             ('OIII5007/Ha', '[O III] 5007 / Hα'), ('Hb/Ha', 'Hβ / Hα'), ('SII6716/6731', '[S II] 6716/6731')]
    runs2 = [('nbrHa_all', HA_C, 'Hα-structure, all'), ('nbrHa_nopatch', NII_C, 'Hα-structure, no N patch'),
             ('nbrHa_patch', SII_C, 'Hα patch region')]
    for k, (rk, lab) in enumerate(right):
        for j, (run, c, rl) in enumerate(runs2):
            if run not in J or rk not in J[run]['ratios']:
                continue
            q = J[run]['ratios'][rk]
            y = k + 0.2 * (j - 1)
            a.plot([q['p2p5'], q['p97p5']], [y, y], '-', color=c, lw=0.8, alpha=0.6)
            a.plot([q['p16'], q['p84']], [y, y], '-', color=c, lw=2.6)
            a.plot([q['energy']], [y], 'o', color=c, ms=3.5, mec='none')
        # LAMOST P1-P4 (Wu+2026 Table 2)
        for reg, mk in (('P1', 'x'), ('P2', 'x'), ('P3', 'x'), ('P4', '+')):
            t2 = pc[reg]
            v = {'NII6583/Ha': t2['N6583'][0] / t2['Ha'][0],
                 'SII6716+31/Ha': (t2['S6716'][0] + t2['S6731'][0]) / t2['Ha'][0],
                 'SII6716/6731': t2['S6716'][0] / t2['S6731'][0]}.get(rk)
            if v is not None:
                a.plot([v], [k + 0.32], mk, color=SUN if reg != 'P4' else FGDIM, ms=5, mew=1.2)
    for j, (run, c, rl) in enumerate(runs2):
        if run in J:
            a.plot([], [], '-', color=c, lw=2.6, label=f'{rl} (n={J[run]["n_fibres"]})')
    a.plot([], [], 'x', color=SUN, label='LAMOST P1–P3')
    a.plot([], [], '+', color=FGDIM, label='LAMOST P4 (bkg)')
    a.set_yticks(range(len(right)))
    a.set_yticklabels([l for _, l in right])
    a.invert_yaxis()
    a.set_xscale('symlog', linthresh=0.03, linscale=0.5)
    a.set_xlim(-0.4, 6)
    a.set_xticks([-0.3, -0.1, 0, 0.03, 0.1, 0.3, 1, 3])
    a.set_xticklabels(['−0.3', '−0.1', '0', '0.03', '0.1', '0.3', '1', '3'])
    a.axvline(0, color=FGDIM, lw=0.6)
    a.set_xlabel('ratio (energy)')
    a.grid(alpha=0.15, axis='x')
    a.legend(loc='lower left', fontsize=6.6, framealpha=0.92)
    fin(fig, FIG / 'spectra_desi_ratios.png')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] in ('ratios', 'all'):
    fig_ratios()
    print('made ratios')
