"""Put the Wu et al. (2026) LAMOST MRS-N relative intensities and velocities next to the SDSS / LAMOST-DR12
measurements of this topic (tag spectra_sdss_lamost).

(1) Absolute scale of the paper's relative intensities.  Table 2 of arXiv:2610.13329 gives every line as the
    ratio to a sky line in the same arm: [N I] 5198 (blue; the 5197.9 component is resolved from 5200.3 at
    R ~ 7500) and OH 6553 (red).  Two anchors are combined:
      a. Halpha: paper relative Halpha vs WHAM-SS total Halpha of the nearest beams (read from the ISM
         teammate's data/derived/spectra_ism_wham_results.json) -> R per red relative unit = I(OH 6553) at
         Xinglong on 2024-02-01.
      b. the ratio I([N I] 5197.9) / I(OH 6553) from (i) the UVES sky atlas (Hanuschik 2003, VizieR
         J/A+A/407/1157 tables 6-7; 'Flux' = peak flux density, integrated = Flux * FWHM * 1.0645) and
         (ii) the SDSS sky fibres of this project (APO; [N I] doublet fitted with 0.6/0.4 profile weights).
    -> [O III] (R) = rel_OIII * I(OH 6553) * I(NI5198)/I(OH6553);  [O III]/Halpha = rel_OIII/rel_Ha * same ratio.
(2) SDSS absolute [O III], Hbeta, Halpha near P1-P4 (sky fibres within 0.30 deg; diffuse-corrected).
(3) Velocities: paper Table 4 vs SDSS strand stacks (substacks JSON) vs LAMOST DR12 MRS super-sky [O III]
    dips (stackfit JSON) vs WHAM Halpha.  V_LSR - V_helio from the ISM teammate (standard solar motion).

Outputs: data/derived/spectra_sdss_lamost_paper_compare.json, figures/spectra_sdss_lamost_paper_compare.png
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import air_to_vac
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, NII_C, SII_C, FGDIM, FG, WARN, ACCENT, OK
import silk_geometry as sg

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
H, C, ARCSEC2 = 6.62607015e-27, 2.99792458e10, (np.pi / 180 / 3600) ** 2

# ---- paper Table 2 (relative intensities) and Table 4 (V_helio), arXiv:2610.13329v1
PAPER_REL = {'P1': dict(OIII=1.67, Ha=0.30, NII=0.24, SII=0.27),
             'P2': dict(OIII=1.83, Ha=0.37, NII=0.26, SII=0.32),
             'P3': dict(OIII=1.90, Ha=0.26, NII=0.20, SII=0.21),
             'P4': dict(OIII=np.nan, Ha=0.25, NII=0.20, SII=0.20)}
PAPER_REL_ERR_HA = 0.02
PAPER_VH = {'P1': dict(OIII=26.8, Ha=16.3, NII=11.8, SII=11.8), 'P2': dict(OIII=26.8, Ha=17.0, NII=12.0, SII=12.0),
            'P3': dict(OIII=23.8, Ha=18.4, NII=11.7, SII=11.7), 'P4': dict(Ha=19.7, NII=13.3, SII=12.7)}
PAPER_V_ERR = {'OIII': 3.5, 'Ha': 1.4}      # paper text: [O III] 3-4 km/s, Halpha ~1.4 km/s
IMAGING_R = {'P1': 12.8, 'P2': 9.9, 'P3': 10.1}   # paper Table 3 SB -> R (paper_review teammate, AB, 60 A width)


def rayleigh_arcsec2(lam_A):
    """erg s-1 cm-2 arcsec-2 of 1 R at lam_A."""
    return 1e6 / (4 * np.pi) * H * C / (lam_A * 1e-8) * ARCSEC2


def hanuschik():
    out = {}
    for tab, lam_air, key in [(6, 5197.9282, 'NI5198'), (6, 5200.2856, 'NI5200'), (6, 5577.3467, 'OI5577'),
                              (7, 6553.6260, 'OH6553')]:
        T = Table.read(ROOT / f'data/raw/vizier_hanuschik2003/table{tab}.tsv', format='ascii.tab', comment='#',
                       data_start=3, names=['Index', 'lambda', 'FWHM', 'Flux'])
        lam = np.array(T['lambda'], float)
        i = int(np.argmin(np.abs(lam - lam_air)))
        integ = float(T['Flux'][i]) * float(T['FWHM'][i]) * 1.0645 * 1e-16
        out[key] = dict(lambda_air=float(lam[i]), peak=float(T['Flux'][i]), fwhm=float(T['FWHM'][i]),
                        integrated_cgs_arcsec2=integ, I_R=integ / rayleigh_arcsec2(lam[i]))
    out['ratio_NI5198_OH6553'] = out['NI5198']['I_R'] / out['OH6553']['I_R']
    out['ratio_NI5198_over_NI5200'] = out['NI5198']['I_R'] / out['NI5200']['I_R']
    return out


def sdss_airglow(G):
    bad = (G['sb5040_ABmag'] < 20.4) | (G['I_OIII_cal_err'] > 0.6)
    g = G[~bad]
    f = np.array(g['f_psf'])
    ni = np.array(g['I_NI5199_cal']) * f * 0.6          # 5197.9 component share of the fitted doublet
    oh = np.array(g['I_OH6555_cal']) * f
    o1 = np.array(g['I_OI5577_cal']) * f
    r = ni / oh
    ok = np.isfinite(r) & (oh > 0.5)
    return dict(n_groups=int(ok.sum()),
                NI5198_R_median=float(np.nanmedian(ni[ok])), NI5198_R_p16_p84=[float(x) for x in np.nanpercentile(ni[ok], [16, 84])],
                OH6553_R_median=float(np.nanmedian(oh[ok])), OH6553_R_p16_p84=[float(x) for x in np.nanpercentile(oh[ok], [16, 84])],
                OI5577_R_median=float(np.nanmedian(o1[ok])),
                ratio_NI5198_OH6553_median=float(np.nanmedian(r[ok])),
                ratio_NI5198_OH6553_p16_p84=[float(x) for x in np.nanpercentile(r[ok], [16, 84])],
                sb5040_ABmag_corrected_median=float(np.nanmedian(np.array(g['sb5040_ABmag']) - 2.5 * np.log10(f))))


def sdss_near_regions(F, G, radius=0.30):
    bad = {(int(r['plate']), int(r['mjd']), int(r['sp'])) for r in G
           if (r['sb5040_ABmag'] < 20.4) or (r['I_OIII_cal_err'] > 0.6)}
    keep = np.array([(int(p), int(m), int(s)) not in bad for p, m, s in zip(F['plate'], F['mjd'], F['sp'])])
    F = F[keep]
    out = {}
    for k, (ra, dec) in sg.REGIONS.items():
        dd = np.hypot((np.array(F['ra']) - ra) * np.cos(np.radians(dec)), np.array(F['dec']) - dec)
        m = dd < radius
        r = dict(n_fibres=int(m.sum()), n_groups=int(len(set(zip(np.array(F['plate'])[m], np.array(F['mjd'])[m])))))
        for ln in ['OIII', 'Hb', 'Ha', 'NII']:
            v = np.array(F[f'I_{ln}_cal'])[m] * np.array(F['f_psf'])[m]
            e = np.array(F[f'I_{ln}_cal_err'])[m] * np.array(F['f_psf'])[m]
            ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
            if ok.sum() < 2:
                r[ln] = r[ln + '_err'] = np.nan
                continue
            w = 1 / e[ok] ** 2
            mu = np.sum(w * v[ok]) / np.sum(w)
            chi = np.sum(w * (v[ok] - mu) ** 2) / max(ok.sum() - 1, 1)
            r[ln] = float(mu)
            r[ln + '_err'] = float(np.sqrt(1 / np.sum(w) * max(chi, 1)))
            r[ln + '_max'] = float(np.nanmax(v[ok]))
        out[k] = r
    return out


def main():
    G = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_groups.csv')
    F = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_fibres.csv')
    W = json.loads((ROOT / 'data/derived/spectra_ism_wham_results.json').read_text())
    SS = json.loads((ROOT / 'data/derived/spectra_sdss_lamost_sdss_substacks.json').read_text())
    LS = json.loads((ROOT / 'data/derived/spectra_sdss_lamost_lamost_stackfit.json').read_text())
    MV = json.loads((ROOT / 'data/derived/spectra_sdss_lamost_mrs_onoff_velocity.json').read_text())
    hn = hanuschik()
    ap = sdss_airglow(G)
    # (a) red-arm scale from WHAM
    wham = {p: W['nearest_beams'][p][0] for p in PAPER_REL}
    scale = {p: wham[p]['I_tot_R'] / PAPER_REL[p]['Ha'] for p in PAPER_REL}
    s_mean = float(np.mean(list(scale.values())))
    s_sd = float(np.std(list(scale.values()), ddof=1))
    conv = {}
    for tag, ratio in [('Hanuschik2003', hn['ratio_NI5198_OH6553']), ('APO_SDSS_median', ap['ratio_NI5198_OH6553_median']),
                       ('APO_SDSS_p16', ap['ratio_NI5198_OH6553_p16_p84'][0]), ('APO_SDSS_p84', ap['ratio_NI5198_OH6553_p16_p84'][1])]:
        c = {'ratio_NI5198_OH6553': ratio, 'I_NI5198_R': ratio * s_mean}
        for p in ('P1', 'P2', 'P3'):
            c[p] = dict(OIII_R=PAPER_REL[p]['OIII'] * ratio * s_mean,
                        OIII_over_Ha_total=PAPER_REL[p]['OIII'] / PAPER_REL[p]['Ha'] * ratio,
                        OIII_over_Ha_excess_vs_P4=(PAPER_REL[p]['OIII'] / (PAPER_REL[p]['Ha'] - PAPER_REL['P4']['Ha']) * ratio
                                                   if PAPER_REL[p]['Ha'] - PAPER_REL['P4']['Ha'] > 2 * PAPER_REL_ERR_HA * 1.414 else None))
        conv[tag] = c
    near = sdss_near_regions(F, G)
    # (3) velocities
    dl = W['vlsr_minus_vhelio_std_kms']
    dlm = float(np.mean(list(dl.values())))
    vel = {'paper_vhelio': PAPER_VH, 'vlsr_minus_vhelio': dl,
           'sdss_OIII_south': dict(v=SS['regions']['south_dec_lt_20']['lines']['OIII']['v_helio'],
                                   err=float(np.hypot(SS['regions']['south_dec_lt_20']['lines']['OIII']['v_helio_err_stat'],
                                                      max(SS['regions']['south_dec_lt_20']['lines']['OIII']['v_helio_err_zp'], 1.0)))),
           'sdss_OIII_Wstrand': dict(v=SS['regions']['W_strand_dec15_20']['lines']['OIII']['v_helio'],
                                     err=float(np.hypot(SS['regions']['W_strand_dec15_20']['lines']['OIII']['v_helio_err_stat'],
                                                        max(SS['regions']['W_strand_dec15_20']['lines']['OIII']['v_helio_err_zp'], 1.0)))),
           'sdss_OIII_north': dict(v=SS['regions']['north_dec_ge_20']['lines']['OIII']['v_helio'],
                                   err=float(np.hypot(SS['regions']['north_dec_ge_20']['lines']['OIII']['v_helio_err_stat'],
                                                      max(SS['regions']['north_dec_ge_20']['lines']['OIII']['v_helio_err_zp'], 1.0)))),
           'sdss_HaNII_all': dict(v=SS['regions']['all']['lines']['HaNII']['v_helio'],
                                  err=float(np.hypot(SS['regions']['all']['lines']['HaNII']['v_helio_err_stat'],
                                                     max(SS['regions']['all']['lines']['HaNII']['v_helio_err_zp'], 1.0)))),
           'lamost_mrs_supersky_OIII_freewidth': {k.replace('OIII_off_', ''): dict(v=v['v_helio'], err=v['v_err'], snr=v['snr'], sigma_kms=v['sigma_kms'])
                                                  for k, v in LS['mrs']['oiii_dip_fits'].items() if k.startswith('OIII_off_')},
           'lamost_mrs_supersky_OIII': {k: dict(v=v['off']['v_helio'], err=v['off']['v_err'], snr=v['off']['snr'])
                                        for k, v in MV.items() if isinstance(v, dict) and 'off' in v and k != 'all_plates'},
           'lamost_mrs_onoff_all_plates': dict(v=MV['all_plates']['on_minus_off']['v_helio'], err=MV['all_plates']['on_minus_off']['v_err'],
                                               F_R=MV['all_plates']['on_minus_off']['F_R']),
           'wham_Ha_vhelio': {p: float(wham[p]['vbar_lsr'] - dl[p]) for p in wham}}
    out = dict(hanuschik2003=hn, sdss_apo_airglow=ap,
               wham_scale_R_per_red_unit=dict(per_region=scale, mean=s_mean, sd=s_sd,
                                              note='= I(OH 6553) at Xinglong 2024-02-01 if LAMOST Halpha = WHAM total Halpha'),
               paper_OIII_converted=conv, sdss_near_regions=near, imaging_R_paper_review=IMAGING_R, velocities=vel,
               vlsr_minus_vhelio_mean=dlm)
    (ROOT / 'data/derived/spectra_sdss_lamost_paper_compare.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: out[k] for k in ['wham_scale_R_per_red_unit', 'sdss_near_regions']}, indent=1, default=float))
    print('Hanuschik', {k: (round(v['I_R'], 2) if isinstance(v, dict) else round(v, 3)) for k, v in hn.items()})
    print('APO', ap)
    for t, c in conv.items():
        print(t, 'ratio %.3f' % c['ratio_NI5198_OH6553'], ' '.join(f"{p}: {c[p]['OIII_R']:.1f} R, OIII/Ha {c[p]['OIII_over_Ha_total']:.2f}" for p in ('P1', 'P2', 'P3')))
    print(json.dumps(vel, indent=0, default=float)[:2500])

    # ------------------------------------------------------------------ figure
    fig, (a1, a2) = new_sky_figure(13, 5.8, 'Silk Nebula · [O III] brightness and velocity: paper vs SDSS / LAMOST DR12',
                                   'left: [O III] 5007 at P1–P3 by method (R); right: heliocentric velocities (km/s)',
                                   ncols=2, gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    xs = np.arange(3)
    regs = ['P1', 'P2', 'P3']
    sets = [('imaging SB (paper Tab. 3)', [IMAGING_R[p] for p in regs], None, WARN, 's'),
            ('MRS-N rel. × UVES sky ratio', [conv['Hanuschik2003'][p]['OIII_R'] for p in regs], None, '#c792ea', 'D'),
            ('MRS-N rel. × APO sky ratio', [conv['APO_SDSS_median'][p]['OIII_R'] for p in regs],
             [[conv['APO_SDSS_median'][p]['OIII_R'] - conv['APO_SDSS_p16'][p]['OIII_R'] for p in regs],
              [conv['APO_SDSS_p84'][p]['OIII_R'] - conv['APO_SDSS_median'][p]['OIII_R'] for p in regs]], '#ce93d8', 'd'),
            ('SDSS sky fibres ≤0.3° (mean)', [near[p]['OIII'] for p in regs], [near[p]['OIII_err'] for p in regs], OIII_C, 'o'),
            ('SDSS sky fibres ≤0.3° (max)', [near[p]['OIII_max'] for p in regs], None, ACCENT, '^')]
    for i, (lab, y, e, col, mk) in enumerate(sets):
        xx = xs + (i - 2) * 0.13
        a1.errorbar(xx, y, yerr=e, fmt=mk, color=col, ms=6, capsize=0, elinewidth=1, label=lab,
                    mfc=col if mk in 'o^' else 'none')
    a1.set_xticks(xs)
    a1.set_xticklabels(regs)
    a1.set_yscale('log')
    a1.set_ylim(0.3, 50)
    a1.set_ylabel('[O III] 5007 (R)')
    a1.grid(alpha=0.25, which='both')
    a1.legend(fontsize=7, loc='upper left', ncol=1)
    # velocities
    rows = []
    for p in ('P1', 'P2', 'P3'):
        rows.append((f'paper [O III] {p}', PAPER_VH[p]['OIII'], PAPER_V_ERR['OIII'], OIII_C, 's'))
    for p in ('P1', 'P2', 'P3', 'P4'):
        rows.append((f'paper Hα {p}', PAPER_VH[p]['Ha'], PAPER_V_ERR['Ha'], HA_C, 's'))
    rows.append(('paper [N II]/[S II] P1–P4', np.mean([PAPER_VH[p]['NII'] for p in PAPER_VH]), 1.0, NII_C, 's'))
    rows.append(('WHAM Hα (4 beams)', np.mean(list(vel['wham_Ha_vhelio'].values())), np.std(list(vel['wham_Ha_vhelio'].values())), HA_C, 'o'))
    rows.append(('SDSS [O III] strand, Dec<20', vel['sdss_OIII_south']['v'], vel['sdss_OIII_south']['err'], OIII_C, 'o'))
    rows.append(('SDSS [O III] W strand', vel['sdss_OIII_Wstrand']['v'], vel['sdss_OIII_Wstrand']['err'], OIII_C, 'o'))
    rows.append(('SDSS [O III] strand, Dec≥20', vel['sdss_OIII_north']['v'], vel['sdss_OIII_north']['err'], OIII_C, 'o'))
    rows.append(('SDSS Hα+[N II] strand', vel['sdss_HaNII_all']['v'], vel['sdss_HaNII_all']['err'], NII_C, 'o'))
    for k, v in vel['lamost_mrs_supersky_OIII'].items():
        rows.append((f'LAMOST MRS super-sky [O III] {k[:8]}', v['v'], v['err'], ACCENT, 'D'))
    rows.append(('LAMOST MRS ON−OFF, 3 plates', vel['lamost_mrs_onoff_all_plates']['v'], vel['lamost_mrs_onoff_all_plates']['err'], ACCENT, 'D'))
    yy = np.arange(len(rows))[::-1]
    for y, (lab, v, e, col, mk) in zip(yy, rows):
        a2.errorbar([v], [y], xerr=[e], fmt=mk, color=col, ms=5.5, capsize=0, elinewidth=1.1,
                    mfc=col if mk == 'o' else 'none')
    a2.set_yticks(yy)
    a2.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    a2.axvspan(23.8 - 3.5, 26.8 + 3.5, color=OIII_C, alpha=0.08, lw=0)
    a2.axvline(-dlm, color=FGDIM, lw=0.7, ls=':')
    label(a2, -dlm + 0.5, -0.9, 'V_LSR = 0', color=FGDIM, size=7.5)
    a2.set_ylim(-1.4, len(rows) - 0.4)
    a2.set_xlabel('V_helio (km/s)')
    a2.set_xlim(-10, 45)
    a2.grid(alpha=0.25)
    finish(fig, ROOT / 'figures/spectra_sdss_lamost_paper_compare.png')


if __name__ == '__main__':
    main()
