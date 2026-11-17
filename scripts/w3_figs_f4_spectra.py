"""F4 spectral evidence: DESI DR1 regression spectra and SDSS/BOSS ON-OFF total-sky stacks (south).

Inputs:
  data/derived/verify_desi_regspec.npz   (wave 2; windows 3727, Hb, 5007, 4363; subset 'south' = Dec < 19.5)
  data/derived/w3_figs_desi_regspec_ha.npz (wave-2 method, Ha + [N II] window; w3_figs_desi_regspec_ha.py)
  data/derived/w3_figs_sdss_stacks.npz  (verify_sdss method, south = Dec < 20; w3_figs_sdss_stacks.py)
SDSS panels get a linear baseline fitted away from the lines (ON-OFF continuum offsets). Both rows are shown as flux density per A per unit [O III]5007 flux (energy units). The SDSS photon stacks are
converted with x (5008/lambda) and divided by the integrated 5007 of the same stack ([O II]: BOSS-only groups).
Measured ratios quoted in the panels are the wave-2 values (verify_desi.md rev. 2, verify_sdss.md).
Output: figures/w3_figs_f4_spectra.png, data/derived/w3_figs_f4_norm.json
"""
import sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, C_KMS, json_dump, panel_tag

V = 13.0
AIR = {'OII3726': 3726.032, 'OII3729': 3728.815, 'Hg': 4340.471, 'OIII4363': 4363.210, 'Hb': 4861.333,
       'OIII4959': 4958.911, 'OIII5007': 5006.843, 'NII6548': 6548.05, 'Ha': 6562.819, 'NII6583': 6583.454,
       'HgI4358': 4358.328}


def a2v(l):
    s2 = (1e4 / l) ** 2
    return l * (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


VAC = {k: a2v(v) for k, v in AIR.items()}
OBS = {k: v * (1 + V / C_KMS) for k, v in VAC.items()}


def line_int(lam, s, c, half=6.0, side=(8.0, 20.0)):
    """Integrated line flux about centre c, minus a linear baseline from side bands."""
    on = np.abs(lam - c) <= half
    sb = (np.abs(lam - c) >= side[0]) & (np.abs(lam - c) <= side[1]) & np.isfinite(s)
    p = np.polyfit(lam[sb], s[sb], 1)
    r = s - np.polyval(p, lam)
    dl = np.gradient(lam)
    from scipy.optimize import curve_fit
    fm = (np.abs(lam - c) <= 8.0) & np.isfinite(r)
    p0 = [np.nansum(r[on] * dl[on]), c, 1.0]
    (fl, cc, sg), _ = curve_fit(lambda x, f, m, s_: gauss(x, m, s_, f), lam[fm], r[fm], p0=p0)
    return float(np.nansum(r[on] * dl[on])), r, float(abs(sg))


def baseline_sub(lam, s, centres, half=4.5):
    """Subtract a linear baseline fitted to pixels farther than `half` A from every line in the panel."""
    m = np.isfinite(s)
    for c in centres:
        m &= np.abs(lam - c) > half
    p = np.polyfit(lam[m], s[m], 1)
    return s - np.polyval(p, lam)


def gauss(lam, c, sig, flux):
    return flux * np.exp(-0.5 * ((lam - c) / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)


def main():
    Dd = np.load(DER / 'verify_desi_regspec.npz')
    Dh = np.load(DER / 'w3_figs_desi_regspec_ha.npz')
    Sd = np.load(DER / 'w3_figs_sdss_stacks.npz')
    rows = {}
    # ---- DESI (already per unit 5007 flux)
    def dband(lamk, sk, bk, src):
        lam = src[lamk]; s = src[sk]; b = src[bk]
        return lam, s, np.nanpercentile(b, 16, axis=0), np.nanpercentile(b, 84, axis=0)
    rows['DESI'] = {'OII': dband('lam_3727', 'S_3727_south', 'bs_3727_south', Dd),
                    'Hb': dband('lam_Hb', 'S_Hb_south', 'bs_Hb_south', Dd),
                    'OIII': dband('lam_5007', 'S_5007_south', 'bs_5007_south', Dd),
                    '4363': dband('lam_4363', 'S_4363_south', 'bs_4363_south', Dd),
                    'Ha': dband('lam_Ha', 'S_Ha_south', 'bs_Ha_south', Dh)}
    I5d, _, sigd = line_int(Dd['lam_5007'], Dd['S_5007_south'], OBS['OIII5007'])
    # ---- SDSS: photon R/A -> energy per unit 5007 energy
    lam = Sd['lam']
    norm = {}
    for key in ('all_south', 'boss_south', 'all_south176'):
        I5, _, sg = line_int(lam, Sd[key], OBS['OIII5007'])
        norm[key] = (I5, sg)
    sigs = norm['all_south'][1]

    def sband(key, lo, hi):
        m = (lam > lo) & (lam < hi)
        I5 = norm[key][0]
        s = Sd[key][m] * (5008.24 / lam[m]) / I5
        return lam[m], s, None, None
    rows['SDSS'] = {'OII': sband('boss_south', 3700, 3760), 'Hb': sband('all_south', 4830, 4895),
                    'OIII': sband('all_south', 4935, 5040), '4363': sband('all_south', 4325, 4395),
                    'Ha': sband('all_south', 6530, 6605)}
    json_dump({'desi_5007_integral_check': I5d, 'desi_5007_sigma_A': sigd,
               'sdss_I5007_R_all_south': norm['all_south'][0], 'sdss_I5007_R_boss_south': norm['boss_south'][0],
               'sdss_5007_sigma_A': sigs}, DER / 'w3_figs_f4_norm.json')
    print('DESI 5007 integral %.3f sigma %.2f A; SDSS I5007 all %.3f boss %.3f R, sigma %.2f A'
          % (I5d, sigd, norm['all_south'][0], norm['boss_south'][0], sigs))

    cols = [('OII', '[O II] 3726/29', (3712, 3746), [('OII3726', ''), ('OII3729', '3726+29')]),
            ('Hb', 'Hβ', (4838, 4888), [('Hb', 'Hβ')]),
            ('OIII', '[O III] 4959/5007', (4945, 5022), [('OIII4959', '4959'), ('OIII5007', '5007')]),
            ('4363', 'Hγ · Hg I · [O III] 4363', (4330, 4380), [('Hg', 'Hγ'), ('HgI4358', 'Hg'), ('OIII4363', '4363')]),
            ('Ha', 'Hα · [N II]', (6540, 6596), [('NII6548', '6548'), ('Ha', 'Hα'), ('NII6583', '6583')])]
    notes = {
        'DESI': {'OII': '[O II]/[O III] = 1.11–1.16', 'Hb': '[O III]/Hβ > 11–18 (2σ)',
                 'OIII': '5007/4959 = 3.04', '4363': '4363/5007 = 0.06–0.09\n(2–3σ; Hg excluded)',
                 'Ha': '[N II]/Hα ≈ 1.0–1.1\n(systematics-limited)'},
        'SDSS': {'OII': '[O II]/[O III] = 1.84\n(BOSS fibres)', 'Hb': '[O III]/Hβ ≈ 23 (13–48)',
                 'OIII': 'south 1.39 ± 0.10 R\n(ON−OFF, r_on 0.15°)', '4363': '4363/5007 = 0.04 ± 0.04\n(not constraining)',
                 'Ha': '[N II]/Hα uncertain ×2'}}
    ylims = {'DESI': {'OII': (-0.15, 0.42), 'Hb': (-0.07, 0.20), 'OIII': (-0.2, 0.6), '4363': (-0.06, 0.075),
                      'Ha': (-0.045, 0.07)},
             'SDSS': {'OII': (-0.15, 0.42), 'Hb': (-0.07, 0.20), 'OIII': (-0.2, 0.6), '4363': (-0.06, 0.075),
                      'Ha': (-0.05, 0.13)}}
    fig, axs = plt.subplots(2, 5, figsize=(15.5, 7.4), gridspec_kw=dict(wspace=0.28, hspace=0.42))
    fig.text(0.012, 0.985, 'Silk Nebula · spectral evidence (southern filament)', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.952, 'flux density per Å per unit [O III] 5007 flux (energy); DESI: template regression, Dec < 19.5, '
             'band = 16–84 % block bootstrap · SDSS/BOSS: total-sky ON−OFF, Dec < 20',
             color=FGDIM, fontsize=9.2, va='top')
    for i, inst in enumerate(('DESI', 'SDSS')):
        sigma = sigd if inst == 'DESI' else sigs
        for j, (key, title, xl, lines) in enumerate(cols):
            ax = axs[i, j]
            lam_, s, lo, hi = rows[inst][key]
            m = (lam_ >= xl[0] - 2) & (lam_ <= xl[1] + 2)
            cen = [OBS[lk] for lk, _ in lines]
            if inst == 'SDSS':
                s = s.copy(); s[m] = baseline_sub(lam_[m], s[m], cen)
            if lo is not None:
                ax.fill_between(lam_[m], lo[m], hi[m], step='mid', color=OIII_C, alpha=0.22, lw=0)
            ax.step(lam_[m], s[m], where='mid', color=FG, lw=0.9)
            if inst == 'SDSS' and key == 'Ha':
                l2, s2, _, _ = sband('all_south176', 6530, 6605)
                m2 = (l2 >= xl[0] - 2) & (l2 <= xl[1] + 2)
                ax.step(l2[m2], baseline_sub(l2[m2], s2[m2], cen), where='mid', color=WARN, lw=0.9, alpha=0.9)
                ax.text(0.03, 0.86, 'white: Dec < 20 (incl. Hα\nstructure near the patch)', transform=ax.transAxes,
                        color=FG, fontsize=6.8, ha='left', va='top', path_effects=[HALO])
                ax.text(0.03, 0.72, 'orange: Dec < 17.6', transform=ax.transAxes, color=WARN, fontsize=6.8,
                        ha='left', va='top', path_effects=[HALO])
            ax.axhline(0, color=FGDIM, lw=0.5)
            for lk, lt in lines:
                c = WARN if lk == 'HgI4358' else (HA_C if lk in ('Ha', 'Hb', 'Hg') else
                                                  (NII_C if lk.startswith('NII') else OIII_C))
                ax.axvline(OBS[lk], color=c, lw=0.7, ls=':', alpha=0.9)
                if lt:
                    ax.text(OBS[lk], 0.985, lt, transform=ax.get_xaxis_transform(), ha='center', va='top',
                            fontsize=7, color=c, path_effects=[HALO])
            if key == 'Hb':
                xx = np.linspace(xl[0], xl[1], 400)
                ax.plot(xx, gauss(xx, OBS['Hb'], sigma, 1 / 3.0), color=HA_C, lw=1.0, ls='--')
                ax.plot(xx, gauss(xx, OBS['Hb'], sigma, 1 / 11.0), color=HA_C, lw=0.9, ls=':')
                label(ax, OBS['Hb'] + 2.2 * sigma, 0.75 * gauss(OBS['Hb'], OBS['Hb'], sigma, 1 / 3.0),
                      '[O III]/Hβ = 3\n(complete shock)', color=HA_C, size=6.8, ha='left', va='center')
                label(ax, OBS['Hb'] - 1.8 * sigma, 0.9 * gauss(OBS['Hb'], OBS['Hb'], sigma, 1 / 11.0),
                      '= 11', color=HA_C, size=6.8, ha='right', va='center')
            if key == 'Ha':
                if inst == 'DESI':
                    label(ax, 0.03, 0.86, 'complete-shock Hα ≈ 0.95\n(off scale, ×13)', color=HA_C, size=6.6,
                          ha='left', va='top', transform=ax.transAxes)
            ax.set_xlim(*xl); ax.set_ylim(*ylims[inst][key])
            ax.tick_params(labelsize=7.5)
            ax.text(0.97, 0.035, notes[inst][key], transform=ax.transAxes, ha='right', va='bottom', fontsize=7.3,
                    color=FG, path_effects=[HALO])
            if i == 1:
                ax.set_xlabel('vacuum λ (Å)', fontsize=8)
            if j == 0:
                ax.set_ylabel(('DESI DR1' if inst == 'DESI' else 'SDSS/BOSS') + '\nper Å per unit [O III]', fontsize=8.5)
            ax.set_title(title, fontsize=9, loc='left', color='#cfd6dd')
    ss.credit(fig)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.08)
    fig.savefig(FIG / 'w3_figs_f4_spectra.png', dpi=200)
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    main()
