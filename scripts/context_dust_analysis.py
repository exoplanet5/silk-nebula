"""Local Bubble wall and 3-D dust toward the Silk Nebula (P1-P4 and along the spine).

Inputs (data/raw/context_dust/, accessed 2026-11-16):
  edenhofer24_los_profiles.npz            Edenhofer+2024 (A&A 685 A82) 3-D dust, Zenodo 10658339, via HTTP Range
                                          (scripts/context_dust_fetch.py); units E (Zhang, Green & Rix 2023) per pc
  lallement22_hpx7_A0.05_dist.vot         distance where Lallement+2022 (A&A 661 A147) cumulative A0 = 0.05 mag,
                                          HEALPix level 7 (Zenodo 22129167, van Wolde, de Bruijne & Brown 2026).
                                          NOTE: this third-party table is NULL for essentially all tiles with
                                          l = 120-210 deg (even Taurus, l=170 b=-15), so it is reported but not used.
  pelgrims2020_L19_map-inner_final.fits   Pelgrims+2020 (A&A 636 A17) inner surface of the Local Bubble shell
                                          (Harvard Dataverse doi:10.7910/DVN/RHPVNC), nside 128 RING
  oneill2024_LB_ShellProperties_MeanDraws.fits  O'Neill+2024 (ApJ 973 136) Local Bubble shell (peak distance d,
                                          d_inner, d_outer, A_G), nside 256 (doi:10.7910/DVN/INB1RB)
Outputs: data/derived/context_dust_los.csv, figures/context_dust_profiles.png
"""
import sys, json
import numpy as np
import pandas as pd
import healpy as hp
from astropy.io import fits
from astropy.io.votable import parse_single_table
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *
from silk_style import *

RAW = ROOT / 'data/raw/context_dust'
AV_PER_E = 2.8   # approximate A_V per unit ZGR23 E (see note)


def main():
    z = np.load(RAW / 'edenhofer24_los_profiles.npz', allow_pickle=False)
    names = list(z['names']); rb = z['r_bounds']; rc = z['r_centers']
    dr = np.diff(rb)
    mean = z['mean'].astype(float).mean(2)            # average the 4 sub-pixels -> (nlos, 516)
    std = z['std'].astype(float).mean(2)
    cum = np.cumsum(mean * dr, 1)
    cum_err = np.sqrt(np.cumsum((std * dr) ** 2, 1))   # upper bound-ish (ignores radial correlations)
    lal = parse_single_table(str(RAW / 'lallement22_hpx7_A0.05_dist.vot')).to_table().to_pandas()
    pel = fits.open(RAW / 'pelgrims2020_L19_map-inner_final.fits')[1].data
    one = fits.open(RAW / 'oneill2024_LB_ShellProperties_MeanDraws.fits', memmap=True)[1].data
    one_l, one_b = np.array(one['l']), np.array(one['b'])
    lal_l, lal_b = lal['l'].values, lal['b'].values
    rows = []
    for i, nm in enumerate(names):
        ra, dec = float(z['ra'][i]), float(z['dec'][i])
        l, b = [x[0] for x in gal(ra, dec)]
        # Lallement: nearest level-7 tile centre (angsep works for any lon/lat pair)
        k = int(np.argmin(angsep(lal_l, lal_b, l, b)))
        d_lal = float(lal['dist'].values[k]) if np.isfinite(lal['dist'].values[k]) else np.nan
        pix = hp.ang2pix(128, l, b, lonlat=True)
        pel_raw, pel_10 = float(pel['r_inner_raw'][pix]), float(pel['r_in_lmax-10'][pix])
        j = np.argmin((one_l - l) ** 2 * np.cos(np.radians(b)) ** 2 + (one_b - b) ** 2)
        tot = cum[i, -1]
        d50 = float(np.interp(0.5 * tot, cum[i], rc))
        m400 = rc < 400
        d_peak = float(rc[m400][np.argmax(mean[i][m400])])
        E = lambda d: float(np.interp(d, rc, cum[i]))
        rows.append(dict(los=nm, ra=ra, dec=dec, l=l, b=b, E_lt_100=E(100), E_lt_160=E(160), E_lt_300=E(300),
                         E_lt_600=E(600), E_lt_1250=tot, e_E_lt_1250=float(cum_err[i, -1]),
                         AV_lt_1250_approx=tot * AV_PER_E, d_half_E1250_pc=d50, d_peak_density_lt400_pc=d_peak,
                         lallement22_d_A0_0p05_pc=d_lal, pelgrims20_r_inner_raw_pc=pel_raw,
                         pelgrims20_r_inner_lmax10_pc=pel_10,
                         oneill24_d_peak_pc=float(one['d'][j]), oneill24_sig_d=float(one['sig_d'][j]),
                         oneill24_d_inner_pc=float(one['d_inner'][j]), oneill24_d_outer_pc=float(one['d_outer'][j]),
                         oneill24_AG_shell=float(one['Ag'][j]), oneill24_l=float(one_l[j]), oneill24_b=float(one_b[j])))
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / 'data/derived/context_dust_los.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print(df[['los', 'l', 'b', 'E_lt_160', 'E_lt_300', 'E_lt_1250', 'd_half_E1250_pc', 'd_peak_density_lt400_pc',
              'lallement22_d_A0_0p05_pc', 'pelgrims20_r_inner_raw_pc', 'pelgrims20_r_inner_lmax10_pc',
              'oneill24_d_inner_pc', 'oneill24_d_peak_pc', 'oneill24_d_outer_pc', 'oneill24_AG_shell']].to_string())

    fig, axes = new_sky_figure(12, 5.4, 'Dust toward the Silk Nebula',
                               'Edenhofer+24 3-D dust density along P1–P4 · Local Bubble wall along the spine', ncols=2)
    ax = axes[0]
    cols = {'P1': OIII_C, 'P2': OK, 'P3': WARN, 'P4': FGDIM, 'east4_P1': '#c792ea', 'west4_P1': DANGER}
    for nm, c in cols.items():
        i = names.index(nm)
        ax.plot(rc, 1e3 * mean[i], color=c, lw=1.5, label=nm.replace('east4_P1', '4° E of P1').replace('west4_P1', '4° W of P1'))
        ax.fill_between(rc, 1e3 * (mean[i] - std[i]), 1e3 * (mean[i] + std[i]), color=c, alpha=0.12, lw=0)
    ax.axvspan(250, 320, color=ACCENT, alpha=0.07, lw=0)
    for d, t, yy in [(160, 'Wu+26 160 pc', 2.0), (288, 'PSR B0656+14 288 pc', 1.1)]:
        ax.axvline(d, color=FGDIM, lw=1, ls=':')
        label(ax, d * 1.02, yy, t, color=FGDIM, size=7.5, va='center')
    label(ax, 255, 0.05, 'Monogem rim\n250–320 pc', color=ACCENT, size=7, va='bottom')
    ax.set_xscale('log'); ax.set_xlim(69, 1250); ax.set_yscale('log'); ax.set_ylim(1e-3, 3)
    ax.set_xlabel('distance (pc)'); ax.set_ylabel('dust density (10$^{-3}$ E per pc, ZGR23 ≈ E(B−V))')
    ax.legend(fontsize=7.5, loc='lower left', ncol=3); ax.grid(alpha=0.3, which='both')
    ax = axes[1]
    sp = df[df.los.str.startswith('arc_dec')]
    ax.plot(sp.dec, sp.pelgrims20_r_inner_lmax10_pc, 's-', color=OK, label='Pelgrims+20 inner surface')
    ax.errorbar(sp.dec, sp.oneill24_d_peak_pc, yerr=sp.oneill24_sig_d, fmt='d-', color=ACCENT, ms=4, label="O'Neill+24 shell peak")
    ax.fill_between(sp.dec, sp.oneill24_d_inner_pc, sp.oneill24_d_outer_pc, color=ACCENT, alpha=0.15, lw=0,
                    label="O'Neill+24 shell inner–outer")
    ax.plot(sp.dec, sp.d_peak_density_lt400_pc, '^-', color='#c792ea', ms=4, label='Edenhofer+24 density peak < 400 pc')
    ax.axhspan(250, 320, color=ACCENT, alpha=0.07, lw=0)
    for d, t in [(160, 'Wu+26 160 pc'), (288, 'PSR B0656+14 288 pc')]:
        ax.axhline(d, color=FGDIM, ls=':', lw=1)
        label(ax, 27.9, d - 4, t, color=FGDIM, size=7.5, ha='right', va='top')
    ax.set_xlabel('Dec along the spine (deg)'); ax.set_ylabel('distance (pc)')
    ax.set_ylim(0, 450); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc='upper left')
    finish(fig, ROOT / 'figures/context_dust_profiles.png')


if __name__ == '__main__':
    main()
