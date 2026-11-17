"""w3_oii_figures.py -- (wave 3, tag oii) figures for the [O II] doublet forensic test.

figures/w3_oii_spectra.png      DESI regression spectra and SDSS stacks around [O II]3727 with fitted doublets.
figures/w3_oii_systematics.png  simulated bias vs velocity / LSF error (DESI, BOSS), BOSS blue-end calibration lines,
                                DESI extraction-phase replica.
figures/w3_oii_summary.png      all 3729/3726 and 6716/6731 measurements vs the low-density limits; 7325/3727 vs T.
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, FG, FGDIM, ACCENT, WARN, DANGER, OK, OIII_C, SII_C, HA_C  # noqa
from w3_oii_common import L, C, design, oii_atoms, oii_ratio  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula/'
D = ROOT + 'data/derived/'
O2 = [L['OII3726'], L['OII3729']]


def fit_curve(x, y, v, s, ratio=None):
    X = design(x, O2, v, [s, s], npoly=2)
    if ratio is not None:
        X = np.vstack([X[0] + ratio * X[1], X[2:]])
    p = np.linalg.lstsq(X.T, y, rcond=None)[0]
    xf = np.linspace(x.min(), x.max(), 600)
    Xf = design(xf, O2, v, [s, s], npoly=2)
    if ratio is not None:
        Xf = np.vstack([Xf[0] + ratio * Xf[1], Xf[2:]])
    return xf, p @ Xf, p


def fig_spectra():
    Z = np.load(D + 'w3_oii_desi_regspec.npz')
    F = json.load(open(D + 'w3_oii_desi_fits.json'))
    S = np.load(D + 'w3_oii_sdss_stacks.npz')
    SJ = json.load(open(D + 'w3_oii_sdss.json'))
    W = np.load(D + 'w3_oii_sdss_wim.npz')
    WJ = json.load(open(D + 'w3_oii_sdss_wim.json'))
    fig, axs = new_sky_figure(13, 8.2, '[O II] 3726/3729 in DESI and SDSS',
                              'filament-correlated spectra (DESI), ON−OFF stacks and the WIM control (SDSS/BOSS); '
                              'solid: best fit, dashed: ratio fixed at 1.46', nrows=2, ncols=3)
    panels = [
        ('DESI · all fibres', Z['lam_3727'], Z['S_geo_all_3727_F'], Z['bs_geo_all_3727_F'], F['geo|all|F']['o2']),
        ('DESI · south (Dec<19.5)', Z['lam_3727'], Z['S_geo_south_3727_F'], Z['bs_geo_south_3727_F'],
         F['geo|south|F']['o2']),
        ('DESI · Hα-correlated (WIM-like)', Z['lam_3727'], Z['S_ha_all_3727_F'], Z['bs_ha_all_3727_F'],
         F['ha|all|F']['o2']),
        ('SDSS · ON−OFF south', S['south_o2_lam'], S['south_o2_S'], S['south_o2_bs'], SJ['onoff']['south']['o2']),
        ('SDSS · ON−OFF Hα patch (30 %)', S['patch30_o2_lam'], S['patch30_o2_S'], S['patch30_o2_bs'],
         SJ['onoff']['patch30']['o2']),
        ('SDSS · WIM (WHAM-Hα regression)', W['lam_o2'], W['within_o2_S'], W['within_o2_bs'], WJ['within']['o2']),
    ]
    for ax, (title, lam, Sp, bs, fit) in zip(axs.ravel(), panels):
        m = (lam > 3712) & (lam < 3745)
        lo, hi = np.percentile(bs, 16, 0), np.percentile(bs, 84, 0)
        ax.fill_between(lam[m], lo[m], hi[m], color='#2a323c', step='mid', lw=0)
        ax.step(lam[m], Sp[m], where='mid', color=FG, lw=0.9)
        sdss = 'SDSS' in title
        if sdss:
            k1, k2 = 'fixed', 'vsfree'
        else:
            k1, k2 = 'kin', 'vsfree'
        f1, f2 = fit[k1], fit[k2]
        xf, yf, _ = fit_curve(lam[m], Sp[m], f1['v'], f1['sigma'])
        ax.plot(xf, yf, color=WARN, lw=1.2, label=f'{k1}: v={f1["v"]:+.0f}, R={f1["ratio"]:.2f}')
        if sdss:
            xf, yf, _ = fit_curve(lam[m], Sp[m], f2['v'], f2['sigma'])
            ax.plot(xf, yf, color=OK, lw=1.2, label=f'v,σ free: v={f2["v"]:+.0f}, R={f2["ratio"]:.2f}')
        xf, yf, _ = fit_curve(lam[m], Sp[m], f2['v'] if sdss else f1['v'], f2['sigma'] if sdss else f1['sigma'],
                              ratio=1.46)
        ax.plot(xf, yf, color=ACCENT, lw=1.0, ls='--', label='R = 1.46')
        for l0 in O2:
            ax.axvline(l0 * (1 + 11 / C), color=FGDIM, lw=0.5, ls=':')
        ax.set_title(title, color=FG, fontsize=9.5, loc='left')
        ax.set_xlabel('vacuum wavelength (Å)')
        ax.legend(loc='upper left', fontsize=7, frameon=False)
        ax.set_xlim(3712, 3745)
    axs[0, 0].set_ylabel('per unit [O III]5007 (Å⁻¹)')
    axs[1, 0].set_ylabel('R Å⁻¹ (as calibrated)')
    axs[1, 2].set_ylabel('R Å⁻¹ per R of WHAM Hα')
    return finish(fig, ROOT + 'figures/w3_oii_spectra.png')


def fig_systematics():
    sim = json.load(open(D + 'w3_oii_sim.json'))
    wc = json.load(open(D + 'w3_oii_sdss_wavecal.json'))['summary']
    SJ = json.load(open(D + 'w3_oii_sdss.json'))
    WJ = json.load(open(D + 'w3_oii_sdss_wim.json'))
    fig, axs = new_sky_figure(13, 5.6, '[O II] doublet fits: where a 1.46 doublet can turn into ~1.3',
                              'noiseless simulations at the measured LSFs; BOSS blue-end calibration from street-lamp Hg '
                              'lines and the [O II] centroids', nrows=1, ncols=3)
    ax = axs[0]
    for inst, col in (('DESI', ACCENT), ('BOSS', WARN)):
        d = sim[inst]['ratio_vs_v_offset_Rtrue1.46']
        x = np.array([float(k) for k in d])
        ax.plot(x, list(d.values()), color=col, lw=1.4, label=f'{inst}: velocity error')
        d2 = sim[inst]['ratio_vs_lsf_true_over_fit_Rtrue1.46']
        ax2x = (np.array([float(k) for k in d2]) - 1) * 100
        ax.plot(ax2x, list(d2.values()), color=col, lw=1.0, ls='--', label=f'{inst}: LSF error (%)')
    ax.axhline(1.46, color=FGDIM, lw=0.6)
    ax.axhline(1.30, color=DANGER, lw=0.6, ls=':')
    label(ax, 8, 1.305, 'reported 1.30', color=DANGER, size=7.5)
    ax.set_xlabel('fitted − true velocity (km/s)  |  true/fitted σ − 1 (%)')
    ax.set_ylabel('fitted 3729/3726 (true 1.46)')
    ax.legend(fontsize=7, frameon=False, loc='upper right')
    ax.set_title('bias vs velocity and LSF error', color=FG, fontsize=9.5, loc='left')
    ax = axs[1]
    # categorical axis ordered by wavelength: street-lamp Hg / airglow lines and the [O II] doublet centroids
    items = []
    for nm, lab in (('Hg3650', 'Hg 3650'), ('Hg3663', 'Hg 3663')):
        q = wc.get(nm)
        if q:
            items.append((lab, q['dv_median'], q['dv_err'], q['dv_err'], SII_C))
    pts = [('[O II] WIM within', WJ['within']['o2']['vfree'], 12.0),
           ('[O II] WIM between', WJ['between']['o2']['vfree'], 12.0),
           ('[O II] Hα patch', SJ['onoff']['patch30']['o2']['vfree'], 14.0),
           ('[O II] north', SJ['onoff']['north']['o2']['vfree'], 16.0),
           ('[O II] all', SJ['onoff']['all']['o2']['vfree'], 12.0),
           ('[O II] south', SJ['onoff']['south']['o2']['vfree'], 10.0)]
    for lab, f, expv in pts:
        items.append((lab, f['v'] - expv, max(f['v'] - f['v_p16'], 0.5), max(f['v_p84'] - f['v'], 0.5), OIII_C))
    for nm, lab in (('Hg4047', 'Hg 4047'), ('Hg4358', 'Hg 4358'), ('OI5577', '[O I] 5577')):
        q = wc.get(nm)
        if q:
            items.append((lab, q['dv_median'], q['dv_err'], q['dv_err'], SII_C))
    for i, (lab, y, elo, ehi, col) in enumerate(items):
        ax.errorbar(i, y, yerr=[[elo], [ehi]], fmt='o' if col == SII_C else 's', color=col, ms=4.5, lw=1)
    ax.set_xticks(range(len(items)))
    ax.set_xticklabels([it[0] for it in items], rotation=60, ha='right', fontsize=7)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.set_ylabel('BOSS centroid − expected (km/s)')
    label(ax, 0.02, 0.04, 'circles: Hg street lamps, [O I] airglow (topocentric)\nsquares: [O II] doublet, v free '
          '(expected: [S II]/DESI velocity)', color=FGDIM, size=6.5, transform=ax.transAxes)
    ax.set_title('BOSS wavelength zero point, blue end', color=FG, fontsize=9.5, loc='left')
    ax = axs[2]
    sp = sim['DESI_specter_replica_Rtrue1.46']
    for sc, col in (('0.55', ACCENT), ('0.62', OK)):
        ks = [k for k in sp if k.startswith(f's{sc}')]
        ph = [sp[k]['phase3729'] - sp[k]['phase3726'] for k in ks]
        ph29 = [sp[k]['phase3729'] for k in ks]
        rr = [sp[k]['ratio_fit_sigmaR'] for k in ks]
        o = np.argsort(ph29)
        ax.plot(np.array(ph29)[o], np.array(rr)[o], 'o-', color=col, ms=4, lw=1, label=f'CCD PSF σ = {sc} Å')
    ax.axhline(1.46, color=FGDIM, lw=0.6)
    ax.set_ylim(1.28, 1.50)
    ax.axhline(1.30, color=DANGER, lw=0.6, ls=':')
    ax.set_xlabel('sub-pixel phase of 3729 on the 0.8 Å grid')
    ax.set_ylabel('fitted 3729/3726 (true 1.46)')
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    ax.set_title('DESI extraction replica (specter ex2d)', color=FG, fontsize=9.5, loc='left')
    return finish(fig, ROOT + 'figures/w3_oii_systematics.png')


def fig_summary():
    A = json.load(open(D + 'w3_oii_atomic.json'))
    m = A['measured']
    fig, axs = new_sky_figure(13, 6.4, '[O II] and [S II] density diagnostics of the Silk',
                              'grey band: low-density limits (atomic-data spread); DESI errors: bootstrap (thin) and ×2 '
                              '(thick, subset scatter)', nrows=1, ncols=3, gridspec_kw={'width_ratios': [1.5, 1, 1]})
    ax = axs[0]
    rows = [('DESI all', 'DESI_all', ACCENT), ('DESI south', 'DESI_south', ACCENT), ('DESI north', 'DESI_north', ACCENT),
            ('DESI Z>0 (on ridge)', 'DESI_Zpos', ACCENT), ('DESI Z<0 (via sky model)', 'DESI_Zneg', ACCENT),
            ('DESI QSO/z>0.5', 'DESI_qso_hiz', ACCENT), ('DESI z<0.5 gal.', 'DESI_lowz_gal', ACCENT),
            ('DESI dark', 'DESI_dark', ACCENT), ('DESI bright', 'DESI_bright', ACCENT),
            ('SDSS south (v fixed)', 'SDSS_south_fixed', WARN), ('SDSS south (v,σ free)', 'SDSS_south_vsfree', OK),
            ('SDSS north (v fixed)', 'SDSS_north_fixed', WARN), ('SDSS north (v,σ free)', 'SDSS_north_vsfree', OK),
            ('SDSS Hα patch (v fixed)', 'SDSS_patch30_fixed', WARN),
            ('SDSS Hα patch (v,σ free)', 'SDSS_patch30_vsfree', OK),
            ('SDSS WIM control (v fixed)', 'SDSS_WIM_within_fixed', WARN),
            ('SDSS WIM control (v,σ free)', 'SDSS_WIM_within_vsfree', OK)]
    ldl = [v['12500'] for k, v in A['low_density_limits'].items() if k != 'SII']
    ax.axvspan(min(ldl), max(ldl), color='#2a323c')
    for i, (nm, key, col) in enumerate(rows):
        if key not in m:
            continue
        q = m[key]
        y = len(rows) - i
        ax.plot([q['p16'], q['p84']], [y, y], color=col, lw=1.2)
        if 'sig_inflated_x2' in q:
            ax.plot([q['R'] - q['sig_inflated_x2'], q['R'] + q['sig_inflated_x2']], [y, y], color=col, lw=3.5,
                    alpha=0.45)
        ax.plot(q['R'], y, 'o', color=col, ms=4.5)
    ax.set_xlim(0.8, 1.9)
    ax.set_yticks([len(rows) - i for i in range(len(rows))])
    ax.set_yticklabels([r_[0] for r_ in rows], fontsize=7.2, color=FG)
    ax.set_xlabel('[O II] 3729/3726')
    at = oii_atoms()['Kal09+FFT04']
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    nes = [1000, 300, 100, 30]
    top.set_xticks([oii_ratio(at, 12500.0, n) for n in nes])
    top.set_xticklabels([str(n) for n in nes], color=FGDIM, fontsize=7)
    top.set_xlabel('n_e (cm⁻³) at 1.25×10⁴ K', color=FGDIM, fontsize=8)
    ax = axs[1]
    s = A['sii']
    rows = [('SDSS WIM (within)', 'SDSS_WIM_within', OK), ('SDSS WIM (between)', 'SDSS_WIM_between', OK),
            ('SDSS Hα patch 30 %', 'SDSS_patch30', OK), ('SDSS patch half-max', 'SDSS_patch_halfmax', OK),
            ('DESI Hα-correlated', 'DESI_Ha_all', ACCENT), ('DESI filament (all)', 'DESI_geo_all', ACCENT),
            ('DESI filament (south)', 'DESI_geo_south', ACCENT), ('SDSS filament south', 'SDSS_south', WARN),
            ('SDSS filament north', 'SDSS_north', WARN)]
    lsii = A['low_density_limits']['SII']['10000']
    ax.axvspan(1.42, max(1.46, lsii), color='#2a323c')
    for i, (nm, key, col) in enumerate(rows):
        if key not in s:
            continue
        q = s[key]
        y = len(rows) - i
        ax.plot([q['p16'], q['p84']], [y, y], color=col, lw=1.2)
        ax.plot(q['R'], y, 'o', color=col, ms=4.5)
    ax.set_xlim(0.6, 3.0)
    ax.set_yticks([len(rows) - i for i in range(len(rows))])
    ax.set_yticklabels([r_[0] for r_ in rows], fontsize=7.2, color=FG)
    ax.set_xlabel('[S II] 6716/6731')
    ax = axs[2]
    cur = A['r7325_3727_curves']
    T = np.array(cur['T'])
    for ne, ls in (('ne1', '-'), ('ne100', '--'), ('ne300', ':')):
        ax.plot(T / 1e4, cur['Kal09+FFT04'][ne], color=OIII_C, ls=ls, lw=1.2, label=f'n_e = {ne[2:]} cm⁻³')
        ax.plot(T / 1e4, cur['T07+FFT04'][ne], color=FGDIM, ls=ls, lw=0.8)
    r7 = A['desi_7325_refit']
    for sn, col, yl in (('all', ACCENT, 0.108), ('south', OK, 0.100)):
        if sn in r7:
            q = r7[sn]['OH_free']
            ax.axhspan(q['ratio_p16'], q['ratio_p84'], color=col, alpha=0.15)
            ax.axhline(q['ratio_7325_3727'], color=col, lw=1.0)
            label(ax, 0.75, yl, f'DESI {sn} (4-line fit, OH terms): {q["ratio_7325_3727"]:.3f}', color=col, size=7)
    ax.axhline(0.026, color=DANGER, lw=0.8, ls=':')
    label(ax, 0.75, 0.0275, 'wave-2: 0.026 (λ7320 amplitude only)', color=DANGER, size=7)
    ax.axhline(0.026 * 1.8, color=DANGER, lw=0.8, ls='--')
    label(ax, 0.75, 0.0485, 'wave-2 × 1.8 (tied λ7330 added)', color=DANGER, size=7)
    ax.set_xlim(0.7, 3.0)
    ax.set_ylim(0, 0.12)
    ax.set_xlabel('T_e (10⁴ K)')
    ax.set_ylabel('[O II] 7320+7330 / 3726+3729')
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    ax.set_title('cyan: Kal09, grey: T07', color=FG, fontsize=9, loc='left')
    return finish(fig, ROOT + 'figures/w3_oii_summary.png')


if __name__ == '__main__':
    for f in (fig_spectra, fig_systematics, fig_summary):
        print(f())
