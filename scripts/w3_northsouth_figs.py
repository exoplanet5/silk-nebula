"""w3_northsouth_figs.py -- figures for wave-3 tag northsouth.
Reads data/derived/w3_northsouth_{desi,desi_velcheck,sdss,rm,geom,shock}.json (+ _desi_spec.npz, _rm_nulls.npz,
verify_shockmodel_tracks.npz, verify_hiwham_hi_extra.json, silk_trace.json).
Writes figures/w3_northsouth_{profiles,velocity,rm,shock}.png (dpi 200).
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, finish, label, FG, FGDIM, OIII_C, HA_C, SII_C, NII_C, HB_C,  # noqa
                        ACCENT, WARN, OK, DANGER, SUN)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula'
DER = f'{ROOT}/data/derived'
FIG = f'{ROOT}/figures'
dj = json.load(open(f'{DER}/w3_northsouth_desi.json'))
vc = json.load(open(f'{DER}/w3_northsouth_desi_velcheck.json'))
sj = json.load(open(f'{DER}/w3_northsouth_sdss.json'))
gj = json.load(open(f'{DER}/w3_northsouth_geom.json'))
rj = json.load(open(f'{DER}/w3_northsouth_rm.json'))
kj = json.load(open(f'{DER}/w3_northsouth_shock.json'))
hx = json.load(open(f'{DER}/verify_hiwham_hi_extra.json'))
MARKS = [(17.6, 'kink'), (19.5, 'N/S split'), (20.54, 'Hα patch'), (21.18, '11.97° rim'), (22.0, 'GALEX limit')]
DESI_C, SDSS_C = OIII_C, WARN


def vmarks(ax, labels=False):
    for x, t in MARKS:
        ax.axvline(x, color='#3a4450', lw=0.8, ls='-', zorder=0)
        if labels:
            label(ax, x + 0.08, 0.96, t, size=6.5, color=FGDIM, transform=ax.get_xaxis_transform(), va='top',
                  rotation=90)


def centres(edges):
    e = np.array(edges)
    return 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])


def sdss_bins():
    rows = []
    for k, v in sj['regions'].items():
        if not k.startswith('bin_') or v is None:
            continue
        lo, hi = [float(x) for x in k[4:].split('-')]
        rows.append((lo, hi, v))
    return sorted(rows, key=lambda r: r[0])


# ===================================================================== Figure 1: profiles along the spine
def fig_profiles():
    fig, axs = new_sky_figure(11, 15.5, 'Silk Nebula · north–south profiles along the spine',
                              'DESI joint Dec-bin templates; SDSS total-sky ON−OFF. Errors: DESI formal × placebo, '
                              'SDSS group bootstrap',
                              nrows=5, ncols=1, sharex=True)
    xc, xw = centres(dj['edges'])
    A = dj['A_all_strands']['quantities']
    Bq = dj['B_Estrand_bins']['quantities']
    sb = sdss_bins()
    # (a) brightness
    ax = axs[0]
    y, e = np.array(A['O3R']['slope']), np.array(A['O3R']['slope_err_infl'])
    ax.errorbar(xc, y, yerr=e, xerr=xw, fmt='o', color=DESI_C, ms=5, lw=1, label='DESI, all strands per bin')
    yE, eE = np.array(Bq['O3R']['slope'][:10]), np.array(Bq['O3R']['slope_err_infl'][:10])
    ax.errorbar(xc + 0.12, yE, yerr=eE, fmt='s', mfc='none', color=ACCENT, ms=4, lw=0.8, label='DESI, E strand only')
    for lo, hi, v in sb:
        o = v['OIII:OIII5007']
        ax.errorbar(0.5 * (lo + hi) - 0.12, o['I_R'], yerr=o['err_boot'], xerr=0.5 * (hi - lo), fmt='D', color=SDSS_C,
                    ms=4, lw=1, label='SDSS ON−OFF (r < 0.15°)' if lo == sb[0][0] else None)
    ax.set_ylabel('[O III] 5007 (R)')
    ax.set_ylim(0, 2.8)
    ax.legend(loc='lower left', fontsize=7, ncol=1)
    vmarks(ax, labels=True)
    # (b) [O II]/[O III]
    ax = axs[1]
    r, er = np.array(A['OII3727sum']['ratio']), np.array(A['OII3727sum']['ratio_err_infl'])
    good = er < 0.8 * np.abs(r) + 1.0
    good &= r > 0
    ax.errorbar(xc[good], r[good], yerr=er[good], xerr=xw[good], fmt='o', color=DESI_C, ms=5, lw=1, label='DESI')
    for k in np.where(~good)[0]:
        label(ax, xc[k], 0.5, 'undetermined', size=6.5, color=FGDIM, ha='center')
    for lo, hi, v in sb:
        q = v['ratios'].get('OII/OIII5007_energy_BOSS')
        if q:
            ax.errorbar(0.5 * (lo + hi) - 0.12, q['value'], yerr=[[q['value'] - q['p16']], [q['p84'] - q['value']]],
                        xerr=0.5 * (hi - lo), fmt='D', color=SDSS_C, ms=4, lw=1,
                        label='SDSS (BOSS groups)' if lo == sb[0][0] else None)
    ax.set_yscale('log')
    ax.set_ylim(0.4, 20)
    ax.set_ylabel('[O II] 3727 / [O III] 5007')
    ax.legend(loc='upper left', fontsize=7)
    vmarks(ax)
    # (c) [O III]/Hb
    ax = axs[2]
    for lo, hi, v in sb:
        q = v['ratios']['OIII5007/Hb_photon']
        ll = v['ratios']['OIII5007/Hb_lower2sig_boot'] * 5008.24 / 4862.68
        x = 0.5 * (lo + hi)
        if q['p16'] > 0:
            ax.errorbar(x - 0.12, q['value'] * 5008.24 / 4862.68,
                        yerr=[[(q['value'] - q['p16']) * 1.03], [(q['p84'] - q['value']) * 1.03]], xerr=0.5 * (hi - lo),
                        fmt='D', color=SDSS_C, ms=4, lw=1)
        ax.annotate('', xy=(x - 0.12, ll * 1.6), xytext=(x - 0.12, ll),
                    arrowprops=dict(arrowstyle='->', color=SDSS_C, lw=0.9))
        ax.plot([x - 0.12 - 0.5 * (hi - lo), x - 0.12 + 0.5 * (hi - lo)], [ll, ll], color=SDSS_C, lw=0.9)
    for zone, (lo, hi), j in (('E_three_zones', (12, 17.6), 0), ('E_three_zones', (17.6, 21), 1),
                              ('E_three_zones', (21, 28.6), 2)):
        h = dj[zone]['quantities']['Hb']
        lim = 1 / (h['ratio'][j] + 2 * h['ratio_err_infl'][j])
        ax.plot([lo, hi], [lim, lim], color=DESI_C, lw=1.4)
        ax.annotate('', xy=(0.5 * (lo + hi), lim * 1.6), xytext=(0.5 * (lo + hi), lim),
                    arrowprops=dict(arrowstyle='->', color=DESI_C, lw=1.2))
    ax.plot([], [], color=DESI_C, lw=1.4, label='DESI 2σ lower limit (zones)')
    ax.plot([], [], color=SDSS_C, lw=0.9, label='SDSS value (16–84 %) and 2σ limit')
    ax.axhspan(1.2, 5.0, color='#ffffff', alpha=0.06, lw=0)
    label(ax, 12.2, 2.6, 'complete shocks', size=6.5, color=FGDIM, ha='left')
    ax.set_yscale('log')
    ax.set_ylim(1, 300)
    ax.set_ylabel('[O III] / Hβ (energy)')
    ax.legend(loc='upper right', fontsize=7)
    vmarks(ax)
    # (d) velocity
    ax = axs[3]
    vb = gj['velocity_bins']
    xv = np.array([0.5 * (b['bin'][0] + b['bin'][1]) for b in vb])
    xwv = np.array([0.5 * (b['bin'][1] - b['bin'][0]) for b in vb])
    ax.errorbar(xv, [b['desi_v_lsr'] for b in vb], yerr=[b['desi_err_formal_infl'] for b in vb], xerr=xwv, fmt='o',
                color=DESI_C, ms=5, lw=1, label='DESI bins (formal × placebo)')
    for zn, z in gj['velocity_zones'].items():
        lo, hi = {'south_<19.5': (12, 19.5), 'north_>=19.5': (19.5, 28.6), 'south_<17.6': (12, 17.6),
                  'mid_17.6-21': (17.6, 21), 'north_>=21': (21, 28.6)}[zn]
        if zn in ('south_<19.5', 'north_>=19.5'):
            continue
        r0, r1 = z['desi_lsr_range']
        ax.fill_between([lo, hi], [r0 - 2, r0 - 2], [r1 + 2, r1 + 2], color=DESI_C, alpha=0.12, lw=0)
        ax.errorbar(0.5 * (lo + hi) + 0.3, z['sdss_lsr'], yerr=z['sdss_err'], fmt='D', color=SDSS_C, ms=6, lw=1.4,
                    mfc='none')
    ax.fill_between([], [], [], color=DESI_C, alpha=0.12, label='DESI zones, estimator range ± 2')
    for lo, hi, v in sb:
        vv = v['velocity_OIII']
        if 'v_corr' in vv:
            corr = np.mean([b['lsr_corr'] for b in vb if b['bin'][0] < hi and b['bin'][1] > lo])
            ax.errorbar(0.5 * (lo + hi) - 0.15, vv['v_corr'] + corr, yerr=vv['ev_boot'], xerr=0.5 * (hi - lo), fmt='D',
                        color=SDSS_C, ms=4, lw=0.8, label='SDSS bins / zones' if lo == sb[0][0] else None)
    for s in hx['A_along_spine']['gauss_segments']:
        for c in s['comps']:
            if c['NHI'] > 6e19 and -25 < c['v'] < 30:
                ax.scatter(s['dec_mid'], c['v'], s=8 + 25 * c['NHI'] / 3e20, color=OK, alpha=0.6, lw=0, zorder=1)
    ax.scatter([], [], s=20, color=OK, alpha=0.6, label='GALFA HI components on spine')
    ax.plot([20.84, 20.84], [-0.74, 6.46], 'x', color=NII_C, ms=7, mew=1.5, label='PGCC G199.09+19.95 CNM')
    ax.set_ylabel('V_LSR (km/s)')
    ax.set_ylim(-14, 26)
    ax.legend(loc='upper left', fontsize=6.6, ncol=3)
    vmarks(ax)
    # (e) distance outside the X-ray rims
    ax = axs[4]
    de = np.array(gj['spine_dec'])
    for nm, R, col, lab in (('Monogem_2SNR_Knies24', 12.10, ACCENT, 'Monogem 2-SNR centre − 12.10°'),
                            ('Monogem_1SNR_Knies24', 11.97, SUN, '1-SNR centre − 11.97°'),
                            ('G205.6+12.4_Knies24', 8.82, HA_C, 'G205.6+12.4 − 8.82°'),
                            ('Silk_curvature_centre', 21.33, OK, 'curvature centre − 21.33°')):
        th = np.array(gj['centres'][nm]['theta'])
        ax.plot(de, th - R, color=col, lw=1.4, label=lab)
    pp = gj['patch_pgcc']['halpha_patch_centroid']
    ax.plot(pp['dec'], pp['theta_from_centres']['Monogem_2SNR_Knies24'] - 12.10, '*', color=HA_C, ms=11,
            mec='k', label='Hα patch (2-SNR centre)')
    ax.axhline(0, color=FGDIM, lw=0.7)
    ax.set_ylabel('θ − R (deg)')
    ax.set_xlabel('Dec of trace (deg)')
    ax.set_xlim(12, 28.6)
    ax.set_ylim(-3, 6)
    ax.legend(loc='upper left', fontsize=7, ncol=2)
    vmarks(ax)
    finish(fig, f'{FIG}/w3_northsouth_profiles.png')


# ===================================================================== Figure 2: velocity
def fig_velocity():
    fig, axs = new_sky_figure(15, 5.6, 'Silk Nebula · north–south [O III] velocity',
                              'DESI joint-template regression spectra (halves); estimator spread; shell projection models',
                              nrows=1, ncols=3, gridspec_kw={'width_ratios': [1.15, 1.0, 1.1]})
    sp = np.load(f'{DER}/w3_northsouth_desi_spec.npz')
    ax = axs[0]
    C = 299792.458
    l0 = 5008.24
    for key, col, lab, vfit in (('velD_0', DESI_C, 'south (trace Dec < 19.5)', dj['vel_D_halves'][0]['v_helio']),
                                ('velD_1', HA_C, 'north (trace Dec ≥ 19.5)', dj['vel_D_halves'][1]['v_helio'])):
        lam, b, s = sp[key]
        m = np.abs(lam - l0) < 16
        x = (lam[m] / l0 - 1) * C
        base = np.median(b[m][np.abs(x) > 300])
        pk = np.max(b[m] - base)
        ax.step(x, (b[m] - base) / pk, where='mid', color=col, lw=1.2, label=f'{lab}: fit {vfit:+.1f} km/s')
        ax.fill_between(x, (b[m] - base - s[m]) / pk, (b[m] - base + s[m]) / pk, step='mid', color=col, alpha=0.18,
                        lw=0)
        ax.axvline(vfit, color=col, lw=0.8, ls='--')
    ax.set_xlim(-300, 300)
    ax.set_xlabel('V_helio relative to [O III] 5008.24 Å (km/s)')
    ax.set_ylabel('regression coefficient (normalised)')
    ax.legend(loc='upper left', fontsize=7)
    label(ax, 0.98, 0.04, 'DESI LSF σ ≈ 42 km/s', transform=ax.transAxes, ha='right', size=7, color=FGDIM)
    # (b) estimators
    ax = axs[1]
    names = {'W1_boot_weights': 'bootstrap σ weights', 'W2_formal_weights': 'formal σ weights',
             'W3_uniform': 'uniform weights', 'L5_5008_only': '5008 only', 'L4_4960_only': '4960 only',
             'MOM_moment': 'model-free moment'}
    for j, (k, nm) in enumerate(names.items()):
        for zn, off, col in (('halves', -0.15, DESI_C), ('three', 0.15, OK)):
            d = vc[zn]['diff_last_minus_first'][k]
            ax.errorbar(d['p50'], j + off, xerr=[[d['p50'] - d['p16']], [d['p84'] - d['p50']]], fmt='o', color=col,
                        ms=4, lw=1)
    sS = sj['regions']['zone_dec<19.5']['velocity_OIII']
    sN = sj['regions']['zone_dec>=19.5']['velocity_OIII']
    dS = sN['v_corr'] - sS['v_corr']
    ax.errorbar(dS, len(names), xerr=np.hypot(sS['ev_boot'], sN['ev_boot']), fmt='D', color=SDSS_C, ms=5, lw=1)
    ax.set_yticks(range(len(names) + 1))
    ax.set_yticklabels(list(names.values()) + ['SDSS (≥19.5 − <19.5)'], fontsize=7.5)
    ax.plot([], [], 'o', color=DESI_C, label='DESI halves (≥19.5 − <19.5)')
    ax.plot([], [], 'o', color=OK, label='DESI zones (≥21 − <17.6)')
    ax.axvline(0, color=FGDIM, lw=0.7)
    ax.axvline(20.9 - 10.6, color=WARN, lw=0.8, ls=':')
    label(ax, 20.9 - 10.6 + 0.3, 6.35, 'wave-2 audit', size=6.5, color=WARN)
    ax.set_xlabel('ΔV_helio north − south (km/s)')
    ax.set_xlim(-12, 26)
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.93), fontsize=6.3)
    # (c) models
    ax = axs[2]
    m = gj['models']
    lo, hi = m['observed_dv_lsr_range']
    ax.axhspan(lo - 2.5, hi + 2.5, color=HA_C, alpha=0.12, lw=0)
    ax.axhspan(lo, hi, color=HA_C, alpha=0.25, lw=0, label='observed ΔV_LSR (N − S)')
    ii = np.linspace(0, 30, 200)
    for vs, col in ((75, SUN), (100, OK), (125, ACCENT), (150, NII_C)):
        ax.plot(ii, 0.85 * vs * np.sin(np.radians(ii)), color=col, lw=1.3, label=f'tilted sheet, v_s = {vs}')
    bn = m['T_sky_plane_bend_beta_north_deg']
    ax.axvline(bn, color=FGDIM, lw=0.9, ls='--')
    label(ax, bn + 0.4, 46, f'sky-plane bend β ≈ {bn:.0f}°', size=7, color=FGDIM)
    ch = [c['centroid_one_sided'] for c in m['CH_one_sided_chord']]
    ax.plot([27.5, 27.5], [min(ch), max(ch)], color=WARN, lw=4, solid_capstyle='butt')
    label(ax, 27.0, max(ch) + 1.5, 'one-sided\nchord', size=7, color=WARN, ha='right')
    hv = gj['velocity_zones']
    dHI = 0.5 * (8.5 + 11.3) - 2.1
    ax.plot([29.2, 29.2], [8.5 - 2.6, 11.3 - 1.8], color=OK, lw=4, solid_capstyle='butt')
    label(ax, 29.6, 2.5, 'HI shift', size=7, color=OK, rotation=90)
    ax.set_xlim(0, 30.5)
    ax.set_ylim(0, 55)
    ax.set_xlabel('tilt of shock normal out of the sky plane, i (deg)')
    ax.set_ylabel('ΔV relative to pre-shock gas (km/s)')
    ax.legend(loc='upper left', fontsize=6.8)
    finish(fig, f'{FIG}/w3_northsouth_velocity.png')


# ===================================================================== Figure 3: RM
def fig_rm():
    z = np.load(f'{DER}/w3_northsouth_rm_nulls.npz')
    fig, axs = new_sky_figure(15, 5.2, 'Silk Nebula · NVSS rotation-measure step: placement null',
                              'Taylor+2009 RMs; statistic = median RM(+2…+6°) − median RM(−6…−2°) about the spine; '
                              'whole trace moved rigidly', nrows=1, ncols=3)
    for ax, h in zip(axs[:2], ('north', 'south')):
        t = rj['true'][h]['rm']['S']
        for mode, col in (('translate', ACCENT), ('spin_translate', SUN)):
            x = z[f'{mode}_{h}_rm']
            p = rj['nulls'][mode][h]['p_abs_ge_true']
            ax.hist(x, bins=np.arange(-110, 111, 6), histtype='step', color=col, lw=1.4,
                    label=f'{mode.replace("_", " + ")}: p(|S| ≥ |S_true|) = {p:.2f}')
        ax.axvline(t, color=HA_C, lw=2, label=f'true {h}: {t:+.1f} ± {rj["true"][h]["rm"]["err_boot"]:.1f} (boot)')
        ax.set_xlabel('E − W median RM step (rad m⁻²)')
        ax.set_ylabel('placements')
        ax.set_title(f'{h} half (trace Dec {"20–28" if h == "north" else "12–20"})', loc='left', fontsize=9, color=FG)
        ax.legend(loc='upper left', fontsize=6.8)
    ax = axs[2]
    for h, col in (('north', HA_C), ('south', OIII_C)):
        p = rj['profiles'][h]
        e = np.array(p['edges'])
        x = 0.5 * (e[:-1] + e[1:])
        y = np.array([np.nan if v is None else v for v in p['median_rm']])
        er = np.array([np.nan if v is None else v for v in p['err']])
        ax.errorbar(x, y, yerr=er, fmt='o-', color=col, ms=4, lw=1.1, label=h)
    g = rj['smooth_gradient']['quadratic']['predicted_step']
    label(ax, 0.02, 0.04, f'smooth quadratic RM gradient predicts {g:+.1f} rad m⁻² for the north step\n'
          f'L1 step location (north): {rj["true"]["north"]["l1_step_location_deg"]:+.1f}° from the spine',
          transform=ax.transAxes, size=7, color=FG)
    ax.axhline(0, color=FGDIM, lw=0.6)
    ax.axvspan(-0.15, 0.15, color='#ffffff', alpha=0.12, lw=0)
    ax.set_xlabel('offset from spine (deg, + = east)')
    ax.set_ylabel('median RM (rad m⁻²)')
    ax.legend(loc='upper right', fontsize=7)
    finish(fig, f'{FIG}/w3_northsouth_rm.png')


# ===================================================================== Figure 4: geometry + shock scenarios
def fig_shock():
    fig, axs = new_sky_figure(16, 6.6, 'Silk Nebula · what changes in the north',
                              'left: geometry; middle: shock-model ratio tracks (verify_shockmodel code); '
                              'right: north/south predictions for a slower shock', nrows=1, ncols=3,
                              gridspec_kw={'width_ratios': [1.0, 1.1, 0.9]})
    tr = json.load(open(f'{DER}/silk_trace.json'))
    ax = axs[0]
    ra = np.array(gj['spine_ra'])
    de = np.array(gj['spine_dec'])
    sc = ax.scatter(ra, de, c=np.array(gj['centres']['Monogem_2SNR_Knies24']['theta']), cmap='viridis', s=12,
                    vmin=11.8, vmax=16, zorder=3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label('θ from Monogem 2-SNR centre (deg)', color=FGDIM, fontsize=7)
    cb.ax.tick_params(colors=FGDIM, labelsize=6.5)
    for k in ('W', 'W2'):
        p = np.array([(q['ra'], q['dec']) for q in tr['strands'][k]['points']])
        ax.plot(p[:, 0], p[:, 1], color=OIII_C, lw=0.8, alpha=0.7)
    pa = np.array(tr['halpha_patch'])
    ax.fill(pa[:, 0], pa[:, 1], color=HA_C, alpha=0.45, lw=0)
    pg = gj['patch_pgcc']['PGCC_G199.09+19.95']
    ax.plot(pg['ra'], pg['dec'], 'x', color=NII_C, ms=8, mew=1.6)
    label(ax, pa[:, 0].min() - 0.15, 20.3, 'Hα patch', size=7, color=HA_C, ha='right')
    label(ax, pg['ra'] - 0.15, pg['dec'] + 0.35, 'PGCC', size=7, color=NII_C, ha='right')
    tt = np.linspace(0, 2 * np.pi, 720)
    for (cra, cde), R, col, lab in (((104.15, 15.30), 12.10, ACCENT, '2-SNR 12.1°'),
                                    ((107.08, 13.50), 11.97, SUN, '1-SNR 11.97°'),
                                    ((110.33, 12.94), 8.82, HA_C, 'G205.6 8.82°'),
                                    ((94.5, 16.76), 21.33, OK, 'curvature 21.3°')):
        c0 = np.array([np.cos(np.radians(cde)) * np.cos(np.radians(cra)),
                       np.cos(np.radians(cde)) * np.sin(np.radians(cra)), np.sin(np.radians(cde))])
        a = np.cross(c0, [0, 0, 1.0])
        a /= np.linalg.norm(a)
        b_ = np.cross(c0, a)
        r = np.radians(R)
        P = np.cos(r) * c0[None] + np.sin(r) * (np.cos(tt)[:, None] * a[None] + np.sin(tt)[:, None] * b_[None])
        cr = np.degrees(np.arctan2(P[:, 1], P[:, 0])) % 360
        cd = np.degrees(np.arcsin(P[:, 2]))
        mm = (cr > 108) & (cr < 124) & (cd > 10) & (cd < 30)
        cr2, cd2 = np.where(mm, cr, np.nan), np.where(mm, cd, np.nan)
        ax.plot(cr2, cd2, color=col, lw=1.0, ls='--', label=lab)
    ax.set_xlim(123, 109)
    ax.set_ylim(11, 29.5)
    ax.set_aspect(1 / np.cos(np.radians(20)))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.legend(loc='lower left', fontsize=6.5)
    # (b) ratio tracks
    ax = axs[1]
    t = np.load(f'{DER}/verify_shockmodel_tracks.npz')
    fit = json.load(open(f'{DER}/verify_shockmodel_fit.json'))
    summ = json.load(open(f'{DER}/verify_shockmodel_summary.json'))
    runs = {(r['v'], r['xn'], r['b'], r['ab'], r['complete']): r['tag'] for r in fit['runs']}
    done = set()
    for sm in summ['consistent']:
        tag = runs.get((sm['v'], sm['xn'], sm['b'], 'asplund09', False))
        o2 = t[tag + '__F_O2'] / t[tag + '__F_O3']
        o3 = t[tag + '__F_O3'] / t[tag + '__F_Hb']
        ax.plot(o3, o2, color=ACCENT, lw=0.6, alpha=0.5)
        i = int(np.argmin(np.abs(np.log10(t[tag + '__N']) - sm['logN'])))
        ax.plot(o3[i], o2[i], 'o', color=ACCENT, ms=3)
        for vN in (70, 80, 90):
            key = (vN, sm['xn'], sm['b'])
            if key in done:
                continue
            done.add(key)
            tg = runs.get((vN, sm['xn'], sm['b'], 'asplund09', False))
            if tg:
                ax.plot(t[tg + '__F_O3'] / t[tg + '__F_Hb'], t[tg + '__F_O2'] / t[tg + '__F_O3'],
                        color={70: SUN, 80: WARN, 90: HA_C}[vN], lw=0.7, alpha=0.6)
    for c in kj['patch_complete_runs']:
        r = c['ratios']
        if r['O3_Hb'] > 0.1:
            ax.plot(r['O3_Hb'], r['O2_Hb'] / r['O3_Hb'], 's', color=FG, ms=4, mfc='none')
    ax.plot([], [], color=ACCENT, label='25 southern models (100–150 km/s), dot = N*')
    for vN, col in ((70, SUN), (80, WARN), (90, HA_C)):
        ax.plot([], [], color=col, label=f'{vN} km/s, same b and X')
    ax.plot([], [], 's', color=FG, mfc='none', label='complete shocks (80–200)')
    ax.add_patch(plt.Rectangle((15, 0.9), 285, 0.6, fc=OIII_C, alpha=0.18, ec=OIII_C))
    ax.add_patch(plt.Rectangle((4, 3.5), 26, 5.5, fc=HA_C, alpha=0.18, ec=HA_C))
    label(ax, 60, 1.15, 'south', size=8, color=OIII_C)
    label(ax, 4.5, 7.6, 'north', size=8, color=HA_C)
    ax.plot(1.8, 13 / 1.8, '*', color=HA_C, ms=12, mec='k')
    label(ax, 1.9, 9.5, 'Hα patch', size=7, color=HA_C)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 300)
    ax.set_ylim(0.3, 20)
    ax.set_xlabel('[O III] / Hβ (cumulative, energy)')
    ax.set_ylabel('[O II] / [O III]')
    ax.legend(loc='upper right', fontsize=6.4)
    # (c) slower-shock predictions
    ax = axs[2]
    V = [r for r in kj['V'] if r['found']]
    for r in V:
        x = r['v_N'] + np.random.default_rng(int(r['v'] * 10 + r['b'])).uniform(-1.5, 1.5)
        ax.plot(x, r['sb_ratio_equalP'], 'o', color=OIII_C, ms=4)
        ax.plot(x, r['fuv_ratio'], '^', color=NII_C, ms=4)
    ax.axhspan(0.36, 0.45, color=OIII_C, alpha=0.2, lw=0)
    ax.axhspan(0.61 - 0.15, 0.61 + 0.15, color=NII_C, alpha=0.15, lw=0)
    ax.plot([], [], 'o', color=OIII_C, label='[O III] N/S at equal pressure')
    ax.plot([], [], '^', color=NII_C, label='FUV/[O III] N/S')
    ax.fill_between([], [], [], color=OIII_C, alpha=0.2, label='observed [O III] N/S')
    ax.fill_between([], [], [], color=NII_C, alpha=0.15, label='observed FUV/[O III] N/S (GALEX)')
    label(ax, 0.03, 0.96, 'later truncation of the 100–150 km/s\nsouthern shocks: no match (0/25)',
          transform=ax.transAxes, size=7, color=FG, va='top')
    ax.set_xlim(62, 106)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel('northern shock speed v_N (km/s)')
    ax.set_ylabel('north / south')
    ax.legend(loc='lower right', fontsize=6.6)
    finish(fig, f'{FIG}/w3_northsouth_shock.png')


if __name__ == '__main__':
    fig_profiles()
    fig_velocity()
    fig_rm()
    fig_shock()
    print('figures written')
