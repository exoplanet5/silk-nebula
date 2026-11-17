"""verify_desi_figs.py -- figures for the wave-2 DESI verification (tag desi). dpi 200, silk_style.

figures/verify_desi_detection.png    null distributions of the matched-filter S/N and split detections
figures/verify_desi_calibration.png  true fibre acceptance vs PSF_TO_FIBER_SPECFLUX; ridge peaks vs Dec (3 conversions)
figures/verify_desi_regspec.png      IV regression spectra: [O III]5007, Hb, [O II]3727, 4340-4385 (Hg vs 4363)
figures/verify_desi_ratios.png       forest plot of differential ratios (mine vs wave 1)
figures/verify_desi_velocity.png     sky-line frame test and [O III] velocity by Dec / strand
figures/verify_desi_extension_map.png  [O III] residual SB (GFA conversion) in the newly fetched Dec 12-15.5 region
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, FG, FGDIM, ACCENT, WARN, DANGER, OK, SUN, OIII_C,  # noqa
                        HA_C, HB_C, SII_C, NII_C)
import matplotlib.pyplot as plt  # noqa: E402
from verify_desi_common import D, ROOT, L, C, trace_strands, load_meta, good_fit  # noqa: E402

FIG = ROOT / 'figures'


def fig_detection():
    j = json.loads((D / 'verify_desi_detect.json').read_text())
    n = np.load(D / 'verify_desi_detect_nulls.npz')
    rot = n['rot']
    fig, ax = new_sky_figure(13, 5.2, 'DESI [O III] ridge · detection nulls',
                             'matched filter against a geometry-only template; true trace vs nulls', ncols=3)
    a = ax[0]
    bins = np.linspace(-30, 60, 91)
    for thr, col, lab in ((0.2, FGDIM, '|r|<0.2'), (0.05, ACCENT, '|r|<0.05')):
        m = np.abs(rot[:, 2]) < thr
        a.hist(rot[m, 0], bins, color=col, alpha=0.6, label=f'rigid rotations {lab} (n={m.sum()})')
    a.hist(n['shuffle'], bins, color=WARN, alpha=0.6, label='shuffle within tile-petal')
    a.hist(n['fake'], bins, color=OK, alpha=0.9, label=f"{len(n['fake'])} fake lines (true trace)")
    a.axvline(j['true']['snr'], color=DANGER, lw=2)
    label(a, j['true']['snr'] - 3, a.get_ylim()[1] * 0.5, f"true trace S/N {j['true']['snr']:.0f}", color=DANGER,
          ha='right')
    a.set_xlabel('matched-filter S/N')
    a.set_ylabel('N')
    a.legend(loc='upper center', fontsize=6.5)
    a.set_xlim(-30, max(160, j['true']['snr'] + 15))
    a = ax[1]
    tr = j['ra_translations']
    xs = [float(k[3:]) for k in tr]
    a.plot(xs, [v['snr'] for v in tr.values()], 'o-', color=ACCENT, label='all')
    a.plot(xs, [v['snr_clip8'] for v in tr.values()], 's--', color=WARN, ms=4, label='|chi|<8')
    a.axhline(0, color=FGDIM, lw=0.6)
    a.set_xlabel('RA shift of the whole trace (deg)')
    a.set_ylabel('S/N')
    a.legend(fontsize=7)
    a = ax[2]
    sp = j['splits']
    names = list(sp.keys())
    a.barh(range(len(names)), [sp[k]['snr'] for k in names], color=OIII_C, alpha=0.8)
    for i, k in enumerate(names):
        label(a, sp[k]['snr'] + 1, i, f"n={sp[k]['n']}", size=6.5, va='center')
    a.set_yticks(range(len(names)))
    a.set_yticklabels(names, fontsize=7)
    a.set_xlabel('S/N in subset (true trace)')
    a.invert_yaxis()
    return finish(fig, FIG / 'verify_desi_detection.png')


def fig_calibration():
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    b = json.loads((D / 'verify_desi_brightness.json').read_text())
    w1 = json.loads((D / 'spectra_desi_summary.json').read_text())['ridge_peak_R']['nodes']
    fig, ax = new_sky_figure(13, 5.4, 'DESI [O III] · point-source calibration vs surface brightness',
                             'SB = F_cal x f_ps / Omega_fibre; f_ps from the GFA, not PSF_TO_FIBER_SPECFLUX',
                             ncols=2, gridspec_kw=dict(width_ratios=[1, 1.6]))
    a = ax[0]
    for p, col in (('dark', ACCENT), ('bright', WARN)):
        m = cm['program'] == p
        a.scatter(cm['mean_psf_to_fiber_specflux'][m], cm['f_chrom_5008'][m], s=3, color=col, alpha=0.4,
                  label=f'{p} (n={m.sum()})')
    a.plot([0.3, 1], [0.3, 1], color=FGDIM, lw=0.8)
    a.axvline(0.789, color=DANGER, lw=0.8, ls=':')
    label(a, 0.792, 0.25, '0.789 = constant for PSF targets', color=DANGER, size=7, rotation=90)
    a.set_xlim(0.65, 0.82)
    a.set_ylim(0.1, 0.85)
    a.set_xlabel('MEAN_PSF_TO_FIBER_SPECFLUX (wave 1)')
    a.set_ylabel('true point-source acceptance f_ps(5008 A), GFA')
    a.legend(fontsize=7, loc='upper left')
    a = ax[1]
    sym = {'E': ('o', OIII_C), 'W': ('s', WARN), 'W2': ('^', OK)}
    fits = b['forward_model']
    for lab, off, alpha, fill in (('wave1_psf_to_fiber', -0.08, 0.5, 'none'), ('gfa_chromatic', 0.08, 1.0, None)):
        for d in fits[lab]['nodes']:
            s, dec = d['node'].split('@')
            if d['err_R'] > 1.0:
                continue
            mk, col = sym[s]
            a.errorbar(float(dec) + off, d['amp_R'], d['err_R'], fmt=mk, color=col, mfc=fill or col, alpha=alpha,
                       ms=5, lw=0.8)
    for d in w1:
        if d['err_R'] < 1.0:
            a.plot(d['dec'], d['amp_R'], 'x', color=FGDIM, ms=5)
    a.axhspan(10, 13, color=DANGER, alpha=0.12)
    label(a, 12.3, 10.6, 'Wu+2026: 10-13 R', color=DANGER, size=7.5)
    a.plot([], [], 'o', color=OIII_C, label='E strand, GFA conversion (this work)')
    a.plot([], [], 's', color=WARN, label='W strand, GFA conversion')
    a.plot([], [], '^', color=OK, label='W2 strand, GFA conversion')
    a.plot([], [], 'o', color=FGDIM, mfc='none', label='my model, wave-1 conversion')
    a.plot([], [], 'x', color=FGDIM, label='wave-1 forward model')
    a.set_xlabel('Dec node (deg)')
    a.set_ylabel('ridge peak [O III]5007 (R)')
    a.set_ylim(-0.5, 13.5)
    a.legend(fontsize=6.8, ncol=2, loc='upper right')
    return finish(fig, FIG / 'verify_desi_calibration.png')


def fig_regspec():
    r = np.load(D / 'verify_desi_regspec.npz')
    f = json.loads((D / 'verify_desi_regspec_fits.json').read_text())
    v = 14.0
    fig, ax = new_sky_figure(13, 7.6, 'DESI [O III]-correlated spectrum · IV regression',
                             'flux density per unit [O III]5007 flux (1/A); shaded: block-bootstrap 16-84%',
                             nrows=2, ncols=3)
    panels = [('5007', 'all', (4950, 5020), [('OIII4959', OIII_C), ('OIII5007', OIII_C)]),
              ('Hb', 'all', (4840, 4885), [('Hb', HB_C)]),
              ('3727', 'all', (3705, 3750), [('OII3726', SII_C), ('OII3729', SII_C)]),
              ('4363', 'all', (4335, 4385), [('Hg', HB_C), ('OIII4363', OIII_C), ('FeII4359', NII_C)]),
              ('4363', 'south', (4335, 4385), [('Hg', HB_C), ('OIII4363', OIII_C), ('FeII4359', NII_C)]),
              ('7300', 'all', (7285, 7345), [('CaII7291', NII_C), ('OII7320', SII_C), ('CaII7324', NII_C),
                                             ('OII7330', SII_C)])]
    for a, (wn, sn, xr, lines) in zip(ax.ravel(), panels):
        lam = r[f'lam_{wn}']
        S = r[f'S_{wn}_{sn}']
        bs = r[f'bs_{wn}_{sn}']
        lo, hi = np.nanpercentile(bs, [16, 84], 0)
        m = (lam >= xr[0]) & (lam <= xr[1])
        a.fill_between(lam[m], lo[m], hi[m], color=ACCENT, alpha=0.25, step='mid')
        a.step(lam[m], S[m], where='mid', color=FG, lw=0.9)
        a.axhline(0, color=FGDIM, lw=0.5)
        yl = a.get_ylim()
        for nm, col in lines:
            x = L[nm] * (1 + v / C)
            a.axvline(x, color=col, lw=0.8, ls='--')
            label(a, x + 0.2, yl[1] * 0.85, nm, color=col, size=6.5, rotation=90, va='top')
        if wn == '4363':
            xh = f[sn]['expected']['HgI4358_at_vcorr']
            a.axvline(xh, color=DANGER, lw=1.2)
            label(a, xh - 0.2, yl[1] * 0.85, 'Hg I 4358 sky', color=DANGER, size=6.5, rotation=90, va='top',
                  ha='right')
            fx = f[sn]['fixed_4363_Hg']
            label(a, xr[0] + 1, yl[0] + 0.08 * (yl[1] - yl[0]),
                  f"4363: {fx['amp'][0]:.3f}+-{fx['err'][0]:.3f}   Hg: {fx['amp'][1]:.3f}+-{fx['err'][1]:.3f}",
                  size=7)
        a.set_title(f'{wn} window · {sn}', fontsize=9, loc='left')
        a.set_xlabel('barycentric vacuum wavelength (A)')
    return finish(fig, FIG / 'verify_desi_regspec.png')


def fig_ratios():
    r = json.loads((D / 'verify_desi_ratios.json').read_text())
    w1 = json.loads((D / 'spectra_desi_ratios.json').read_text())['model_all']['ratios']
    rows = [('Hb', 'Hb/[O III]', None), ('OIII4959', '4959/5007', 'OIII5007/OIII4959'),
            ('OIII4363_jointHg', '4363/5007', None), ('HgI4358', 'Hg I 4358 (sky) / 5007', None),
            ('HeII4686', 'He II/[O III]', None), ('NeIII3869', '[Ne III]/[O III]', None),
            ('OII3726', '[O II]3726/[O III]', None), ('OII3729', '[O II]3729/[O III]', None),
            ('Ha', 'Ha/[O III]', None), ('NII6583', '[N II]6583/[O III]', None),
            ('SII6716', '[S II]6716/[O III]', None), ('SIII9069', '[S III]9069/[O III]', None),
            ('OII7325_jointCa', '[O II]7325/[O III]', None), ('CaII7291', '[Ca II]7291/[O III]', None),
            ('FeII7155', '[Fe II]7155/[O III]', None), ('sky5577', '[O I]5577 sky / 5007', None)]
    subs = [('all', FG, 0), ('petalFE', ACCENT, -0.25), ('QSO_or_z>0.5', WARN, -0.12),
            ('south_dec<19.5', OK, 0.12), ('nbr_instrument', NII_C, 0.25)]
    fig, ax = new_sky_figure(9.5, 8.5, 'DESI filament spectrum · differential line ratios',
                             'energy ratio to [O III]5007, IV estimator; bars: block-bootstrap 16-84 %')
    a = ax
    for i, (k, lab, _) in enumerate(rows):
        for sn, col, off in subs:
            d = r.get(sn, {}).get(k)
            if d is None:
                continue
            a.errorbar(d['ratio'], i + off, xerr=[[max(d['ratio'] - d['p16'], 0)], [max(d['p84'] - d['ratio'], 0)]],
                       fmt='o', color=col, ms=3.5, lw=0.9, label=sn if i == 0 else None)
    a.set_yticks(range(len(rows)))
    a.set_yticklabels([x[1] for x in rows], fontsize=7.5)
    a.axvline(0, color=FGDIM, lw=0.6)
    a.set_xscale('symlog', linthresh=0.02)
    a.set_xlim(-0.08, 3)
    a.invert_yaxis()
    a.set_xlabel('ratio (symlog)')
    a.legend(fontsize=7, loc='lower right')
    return finish(fig, FIG / 'verify_desi_ratios.png')


def fig_velocity():
    v = json.loads((D / 'verify_desi_velocity.json').read_text())
    sk = np.load(D / 'verify_desi_velocity_sky.npz')
    fig, ax = new_sky_figure(13, 5.4, 'DESI [O III] velocity · frame and strands',
                             'left: [O I]5577 sky-line variance centroid vs barycentric correction', ncols=2)
    a = ax[0]
    a.scatter(sk['OI5577_vcorr'], sk['OI5577_vcen'], s=2, color=ACCENT, alpha=0.3)
    xs = np.linspace(-25, 32, 10)
    t = v['sky_frame_test']['OI5577']
    a.plot(xs, t['intercept_kms'] + t['slope'] * xs, color=WARN, label=f"fit slope {t['slope']:.3f}")
    a.plot(xs, t['intercept_kms'] + t['sim_slope_stacked_width'] * xs, color=OK, ls='--',
           label=f"barycentric frame, simulated estimator (slope {t['sim_slope_stacked_width']:.2f})")
    st = np.array(v['sky_frame_test']['OI5577_stacked']['bins'])
    a.errorbar(st[:, 0], st[:, 1], st[:, 2], fmt='D', color=SUN, ms=6, zorder=5,
               label=f"stacked profiles by v_corr (slope {v['sky_frame_test']['OI5577_stacked']['slope']:.2f})")
    a.plot(xs, 0 * xs, color=DANGER, ls=':', label='no barycentric correction (slope 0)')
    a.set_xlabel('EFFTIME-weighted v_corr of the coadd (km/s)')
    a.set_ylabel('sky-line centroid (km/s, barycentric frame)')
    a.set_ylim(-40, 50)
    a.legend(fontsize=7, loc='upper left')
    a = ax[1]
    rs = v['regspec_centroids']
    bands = [('dec12-15.5', 13.75), ('dec15.5-17.5', 16.5), ('dec17.5-19.5', 18.5), ('dec19.5-22', 20.75),
             ('dec22-25', 23.5), ('dec25-28.6', 26.8)]
    for k, y in bands:
        if k in rs:
            d = rs[k]
            a.errorbar(y, d['v'], [[d['v'] - d['p16']], [d['p84'] - d['v']]], fmt='o', color=OIII_C, ms=5,
                       mfc='none')
            p = d.get('per_fibre_pos_SN4')
            if p and p['n'] >= 10:
                a.errorbar(y - 0.15, p['wmean'], p['err'], fmt='o', color=OIII_C, ms=5)
    for nm, x in (('P1', 16.60), ('P2', 17.74), ('P3', 14.72)):
        p = rs.get(f'within0.5deg_{nm}', {}).get('per_fibre_pos_SN4')
        if p:
            a.errorbar(x - 0.12, p['wmean'], p['err'], fmt='*', color=SUN, ms=9)
    w1 = json.loads((D / 'spectra_desi_ratios.json').read_text())['v_by_dec']
    for d in w1:
        a.errorbar(0.5 * (d['dec_lo'] + d['dec_hi']) + 0.15, d['v'], [[d['v'] - d['v_boot_p16']],
                   [d['v_boot_p84'] - d['v']]], fmt='x', color=FGDIM, ms=5, lw=0.7)
    for k, col, x in (('E_strand_dec15.5-17.5', OK, 16.2), ('W_strand_dec15.5-17.5', WARN, 16.8)):
        if k in rs and rs[k].get('per_fibre_pos_SN4'):
            p = rs[k]['per_fibre_pos_SN4']
            a.errorbar(x, p['wmean'], p['err'], fmt='s', color=col, ms=5)
    for nm, (dec, vv, ev) in {'P1': (16.60, 26.8, 4.3), 'P2': (17.74, 26.8, 4.1), 'P3': (14.72, 23.8, 3.3)}.items():
        a.errorbar(dec, vv, ev, fmt='D', color=DANGER, ms=5)
        label(a, dec + 0.12, vv + {'P1': 5.0, 'P2': 2.5, 'P3': 1.0}[nm], f'LAMOST {nm}', color=DANGER, size=7)
    a.axhline(rs['all']['v'], color=OIII_C, lw=0.8, ls='--')
    a.plot([], [], 'o', color=OIII_C, label='positive fibres S/N>4, by Dec band')
    a.plot([], [], 'o', color=OIII_C, mfc='none', label='IV regression, by Dec band')
    a.plot([], [], '*', color=SUN, ms=9, label='fibres within 0.5 deg of P1/P2/P3')
    a.plot([], [], 's', color=OK, label='E strand fibres, Dec 15.5-17.5')
    a.plot([], [], 's', color=WARN, label='W strand fibres, Dec 15.5-17.5')
    a.plot([], [], 'x', color=FGDIM, label='wave-1 bands')
    a.set_xlabel('Dec (deg)')
    a.set_ylabel('[O III] V_helio (km/s)')
    a.set_ylim(-5, 35)
    a.legend(fontsize=6.5, loc='lower right', ncol=2)
    return finish(fig, FIG / 'verify_desi_velocity.png')


def fig_extension():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    keep, meta = load_meta(lf['sparcl_id'])
    ok = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & np.isfinite(cm['kR_5008'])
    SB = lf['o3__amp'][:, 0] * cm['kR_5008']
    fig, ax = new_sky_figure(11, 6.0, 'DESI [O III] residual · Dec 12.5-18.5',
                             'per-fibre residual SB (R, GFA conversion); white rim: wave-2 fetch; lines: E/W/W2 strands',
                             ncols=2)
    fig.subplots_adjust(left=0.07, right=0.88, bottom=0.1, top=0.88, wspace=0.18)
    for a, (lo, hi) in zip(ax, ((12.5, 15.6), (15.4, 18.5))):
        m = ok & (meta['dec'] > lo) & (meta['dec'] < hi) & (meta['ra'] > 114.5) & (meta['ra'] < 118.5)
        sc = a.scatter(meta['ra'][m], meta['dec'][m], c=np.clip(SB[m], -2, 3), cmap='coolwarm', vmin=-2, vmax=3,
                       s=10, edgecolors=np.where(lf['src'][m] == 'ext', '#ffffff', 'none'), linewidths=0.3)
        for s, arr in trace_strands().items():
            a.plot(arr[:, 0], arr[:, 1], color=SUN if s == 'E' else WARN, lw=0.8, alpha=0.8)
        for nm, (x, y) in {'P1': (116.65, 16.60), 'P2': (116.21, 17.74), 'P3': (116.45, 14.72),
                           'P4': (115.20, 16.18)}.items():
            if lo < y < hi:
                a.plot(x, y, '+', color=FG, ms=9)
                label(a, x + 0.05, y + 0.05, nm, size=8)
        a.set_xlim(118.5, 114.5)
        a.set_ylim(lo, hi)
        a.set_aspect(1 / np.cos(np.radians(0.5 * (lo + hi))))
        a.set_xlabel('RA (deg)')
        a.set_ylabel('Dec (deg)')
    cb = fig.colorbar(sc, cax=fig.add_axes([0.9, 0.15, 0.015, 0.7]))
    cb.set_label('[O III]5007 residual (R)')
    return finish(fig, FIG / 'verify_desi_extension_map.png', tight=False)


if __name__ == '__main__':
    which = sys.argv[1:] or ['detection', 'calibration', 'regspec', 'ratios', 'velocity', 'extension']
    for w in which:
        print(globals()[f'fig_{w}']())
