"""verify_desi_audit_fig.py -- figure + extra numbers for the resume audit (tag desi).

Reads data/derived/verify_desi_audit.json, verify_desi_audit_vel.json, verify_desi_audit_regspec.npz,
verify_desi_audit_boxfluxes.npz, verify_desi_coadd_meta.npz and verify_desi_ratios.json (first-pass Gaussian fitter).
Writes data/derived/verify_desi_audit_extra.json and figures/verify_desi_audit.png.
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, finish, label, FG, FGDIM, ACCENT, WARN, DANGER, OK, OIII_C, HA_C,  # noqa
                        NII_C, HB_C, SII_C)
from verify_desi_common import D, ROOT, L, C  # noqa: E402

HC = 6.62607015e-27 * 2.99792458e18
ARCSEC2_SR = (np.pi / 648000) ** 2
OMEGA = np.pi * 0.76 ** 2


def one_R(lam):
    return 1e6 / (4 * np.pi) * HC / lam * ARCSEC2_SR


def main():
    a = json.load(open(D / 'verify_desi_audit.json'))
    v = json.load(open(D / 'verify_desi_audit_vel.json'))
    g = json.load(open(D / 'verify_desi_ratios.json'))
    rs = np.load(D / 'verify_desi_audit_regspec.npz')
    bf = np.load(D / 'verify_desi_audit_boxfluxes.npz')
    meta = np.load(D / 'verify_desi_coadd_meta.npz')
    assert (bf['sparcl_id'] == meta['sparcl_id']).all()

    # --- global matched-filter amplitude of the box [O III] in R (GFA vs wave-1 conversion)
    O, U, Z, dec = bf['E_OIII5007'], bf['U_OIII5007'], meta['Z'], meta['dec']
    w1fac = np.clip(meta['mean_psf_to_fiber_specflux'], 0.3, 1.0) / meta['f_chrom_5008']
    extra = {}
    k = 1e-17 / OMEGA / one_R(L['OIII5007'])
    for nm, sel in (('all', np.ones(Z.size, bool)), ('south_dec<19.5', dec < 19.5), ('north_dec>=19.5', dec >= 19.5)):
        m = sel & np.isfinite(O) & np.isfinite(U) & (U > 0) & np.isfinite(Z)
        e = 1 / np.sqrt(U[m])
        m2 = np.where(m)[0][(e < 4 * np.median(e)) & (np.abs(O[m] / e) < 30)]
        A = np.sum(U[m2] * Z[m2] * O[m2]) / np.sum(U[m2] * Z[m2] ** 2)
        sA = 1 / np.sqrt(np.sum(U[m2] * Z[m2] ** 2))
        Ow = O[m2] * w1fac[m2]
        Uw = U[m2] / w1fac[m2] ** 2
        Aw = np.sum(Uw * Z[m2] * Ow) / np.sum(Uw * Z[m2] ** 2)
        extra[nm] = dict(n=int(m2.size), amp_R_gfa=float(A * k), amp_R_gfa_err_formal=float(sA * k),
                         amp_R_wave1conv=float(Aw * k), ratio=float(A / Aw))
    print(extra)

    # --- inflated ratio limits (box method, nside16 calibration of the bootstrap with null boxes)
    infl = a['null_summary']['rms_over_sig16']
    infl32 = a['null_summary']['rms_over_sig32']
    lim = {}
    for meth, src in (('box', a['ratios']['all']), ('gauss', g['all'])):
        hb = src['Hb']
        s = hb['nside16']['sig_boot'] * infl if meth == 'box' else hb['sig_boot'] * 1.5
        up = hb['ratio'] + 2 * s
        o2 = src['OII3727sum']['ratio']
        o3 = 1 + src['OIII4959']['ratio']
        lim[meth] = dict(Hb=hb['ratio'], Hb_sig_inflated=s, Hb_upper_2sig=up, OIII_over_Hb_lower_2sig=1 / up,
                         R23_lower_2sig=(o2 * 1.0 + o3) / up,
                         OIII4363=src['OIII4363']['ratio'] if meth == 'box' else src['OIII4363_jointHgFe']['ratio'],
                         OIII4363_sig_inflated=(src['OIII4363']['nside16']['sig_boot'] * infl if meth == 'box'
                                                else src['OIII4363_jointHgFe']['sig_boot'] * 1.8))
    extra['limits'] = lim
    extra['inflation_box_nside16'] = infl
    extra['inflation_box_nside32'] = infl32
    print(json.dumps(lim, indent=1))
    json.dump(extra, open(D / 'verify_desi_audit_extra.json', 'w'), indent=1)

    # ------------------------------------------------------------------ figure
    fig, ax = new_sky_figure(13, 8.4, 'DESI DR1 audit · box fluxes vs Gaussian fits',
                             'IV ratios to [O III]5007 (energy), HEALPix bootstrap 68 %; 4363 regression spectrum; '
                             '[O III] velocity', nrows=2, ncols=2)
    # (1) ratios forest
    a1 = ax[0, 0]
    lines = [('Hb', 'Hβ', HB_C), ('Ha', 'Hα', HA_C), ('NII6583', '[N II]6583', NII_C), ('OIII4363', '[O III]4363',
                                                                                         OIII_C),
             ('HgI4358', 'Hg I 4358 (sky)', WARN), ('HeII4686', 'He II 4686', FGDIM), ('sky5577', '[O I]5577 (sky)',
                                                                                        DANGER)]
    sets = [('all', 'all', 0.0), ('Zpos', 'Z>0 (on ridge)', 0.18), ('Zneg', 'Z<0 (petal-sky side)', 0.36),
            ('south_dec<19.5', 'Dec<19.5', 0.54)]
    for j, (key, nm, col) in enumerate(lines):
        for sname, slab, off in sets:
            r = a['ratios'][sname].get(key)
            if r is None:
                continue
            y = j + off
            a1.errorbar(r['ratio'], y, xerr=[[r['ratio'] - r['p16']], [r['p84'] - r['ratio']]], fmt='o', ms=4,
                        color=col, alpha=1.0 - 0.18 * sets.index((sname, slab, off)), capsize=0, lw=1.2)
        gk = {'OIII4363': 'OIII4363_jointHgFe', 'HgI4358': 'HgI4358_jointFe'}.get(key, key)
        if gk in g['all']:
            r = g['all'][gk]
            a1.plot(r['ratio'], j - 0.16, marker='D', ms=4, color='white', mfc='none')
    a1.axvline(0, color=FGDIM, lw=0.6)
    a1.set_yticks(range(len(lines)))
    a1.set_yticklabels([nm for _, nm, _ in lines])
    a1.set_xlim(-0.15, 0.32)
    a1.invert_yaxis()
    a1.set_xlabel('line / [O III]5007 (energy)')
    label(a1, 0.15, 0.35, 'rows: all · Z>0 · Z<0 · Dec<19.5\n◇ first-pass Gaussian fitter', size=7.5, color=FGDIM)
    a1.set_title('Box-integral IV ratios: Hα asymmetric, sky lines non-zero', loc='left', color=FG)

    # (2) 4363 regression spectrum
    a2 = ax[0, 1]
    lam, sp, sg, mod = rs['lam4363'], rs['spec4363'], rs['sig4363'], rs['model4363']
    a2.fill_between(lam, sp - sg, sp + sg, step='mid', color=ACCENT, alpha=0.25, lw=0)
    a2.step(lam, sp, where='mid', color=ACCENT, lw=1.0)
    a2.plot(lam, mod, color=WARN, lw=1.0)
    r4 = a['regspec_4363']
    for x, t, c in ((r4['expected_OIII4363'], '[O III]4363', OIII_C), (r4['expected_Hg'], 'Hg I 4358', WARN),
                    (L['Hg'] * (1 + 13 / C), 'Hγ', HB_C), (L['FeII4359'] * (1 + 13 / C), '[Fe II]', FGDIM)):
        a2.axvline(x, color=c, lw=0.7, ls='-')
        label(a2, x + 0.3, 0.92 * np.nanmax(sp + sg), t, color=c, size=7.5, rotation=90, va='top')
    a2.axhline(0, color=FGDIM, lw=0.6)
    a2.set_xlim(4340, 4390)
    a2.set_xlabel('vacuum barycentric wavelength (Å)')
    a2.set_ylabel('per unit [O III]5007 (Å⁻¹)')
    a2.set_title('Pixel IV spectrum: peak %.2f Å ([O III] %.2f)' %
                 (r4['centre_4363'], r4['expected_OIII4363']), loc='left', color=FG)

    # (3) velocities
    a3 = ax[1, 0]
    rows = [('total_all', 'Z', 'all, total Z'), ('joint_all', 'Z_E', 'E strand (joint)'),
            ('joint_all', 'Z_W', 'W strand (joint)'), ('total_dec<19.5', 'Z', 'Dec < 19.5'),
            ('total_dec>=19.5', 'Z', 'Dec ≥ 19.5'), ('joint_dec15.5-17.5', 'Z_E', 'E, Dec 15.5–17.5'),
            ('joint_dec15.5-17.5', 'Z_W', 'W, Dec 15.5–17.5')]
    for j, (s, kk, nm) in enumerate(rows):
        r = v[s][kk]
        lo, hi = r['v_boot_p16_p84']
        a3.errorbar(r['v'], j, xerr=[[max(r['v'] - lo, 0)], [max(hi - r['v'], 0)]], fmt='o', color=OIII_C, ms=4)
    a3.axvspan(23.8 - 3.3, 26.8 + 4.3, color=DANGER, alpha=0.15, lw=0)
    label(a3, 25.3, -0.45, 'LAMOST P1–P3', color=DANGER, size=7.5, ha='center')
    a3.axvspan(10.0 - 1.8, 10.0 + 1.8, color=OK, alpha=0.15, lw=0)
    label(a3, 10.0, -0.45, 'SDSS', color=OK, size=7.5, ha='center')
    a3.set_ylim(len(rows) - 0.5, -0.8)
    a3.set_yticks(range(len(rows)))
    a3.set_yticklabels([nm for _, _, nm in rows])
    a3.set_xlim(-5, 35)
    a3.set_xlabel('[O III] V_helio (km/s), regression spectrum, both [O III] lines')
    a3.set_title('Velocity: no robust E/W difference; north redder', loc='left', color=FG)

    # (4) null boxes
    a4 = ax[1, 1]
    ns = a['null_summary']
    rr = np.array(ns['ratios'])
    s16 = np.array(ns['sig_boot16'])
    s32 = np.array(ns['sig_boot32'])
    xx = np.arange(rr.size)
    a4.errorbar(xx - 0.1, rr, yerr=s32, fmt='o', color=ACCENT, ms=4, label='nside-32 bootstrap')
    a4.errorbar(xx + 0.1, rr, yerr=s16, fmt='s', color=WARN, ms=3, mfc='none', label='nside-16 bootstrap')
    for key, col, nm in (('sky5577', DANGER, '[O I]5577'), ('HgI4358', WARN, 'Hg 4358')):
        r = a['ratios']['all'][key]
        a4.errorbar(rr.size + (0 if key == 'sky5577' else 1), r['ratio'], yerr=r['sig_boot'], fmt='D', color=col,
                    ms=4)
    a4.axhline(0, color=FGDIM, lw=0.6)
    a4.set_xticks(list(xx) + [rr.size, rr.size + 1])
    a4.set_xticklabels([n.replace('null_', '').replace('OIII', 'O3 ').replace('_', ' ') for n in ns['names']] +
                       ['5577', 'Hg'], rotation=60, fontsize=7)
    a4.set_ylabel('IV ratio of a line-free box')
    a4.legend(loc='upper left', ncol=2, fontsize=7)
    a4.set_title('Null boxes: rms/σ = %.2f (nside 32), %.2f (nside 16)' % (ns['rms_over_sig32'],
                                                                          ns['rms_over_sig16']),
                 loc='left', color=FG)
    finish(fig, ROOT / 'figures/verify_desi_audit.png')


if __name__ == '__main__':
    main()
