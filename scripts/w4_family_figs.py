"""w4_family_figs.py -- wave-4 `family`: comparison figure (Silk vs Monogem-periphery filaments) and summary JSON.

Inputs : data/derived/w4_family_desi.json, w4_family_desi_stacks.npz, w4_family_sdss_e1.json, w4_family_traces.json,
         w4_family_nsns_wide.npz; Silk reference values quoted from notes/w3_canon.md (data/derived/w3_canon_values.json)
Outputs: figures/w4_family_compare.png, data/derived/w4_family_summary.json
Usage  : ~/.venvs/astro313/bin/python -I scripts/w4_family_figs.py
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from scipy.ndimage import median_filter, gaussian_filter
from silk_style import *

ROOT = '/Users/mickey/astronomy/silk-nebula'
C = 299792.458
L5007 = 5008.240
LSR = {'E1': -12.11, 'SILK': -12.37}            # V_LSR - V_helio (standard solar motion), derived here (astropy)


def circle(ra0, dec0, r, n=721):
    t = np.radians(np.linspace(0, 360, n))
    d0, a0, rr = np.radians(dec0), np.radians(ra0), np.radians(r)
    dec = np.arcsin(np.sin(d0) * np.cos(rr) + np.cos(d0) * np.sin(rr) * np.cos(t))
    ra = a0 + np.arctan2(np.sin(t) * np.sin(rr) * np.cos(d0), np.cos(rr) - np.sin(d0) * np.sin(dec))
    return np.degrees(ra), np.degrees(dec)


def main():
    dj = json.load(open(f'{ROOT}/data/derived/w4_family_desi.json'))
    st = np.load(f'{ROOT}/data/derived/w4_family_desi_stacks.npz')
    sj = json.load(open(f'{ROOT}/data/derived/w4_family_sdss_e1.json'))
    tj = json.load(open(f'{ROOT}/data/derived/w4_family_traces.json'))
    E, S, K = dj['regions']['E1'], dj['regions']['S1'], dj['regions']['SILK']
    vE = E['velocity_o3']['excl_41485_41487']; vK = K['velocity_o3']['excl_41485_41487']
    summ = {
        'E1': dict(
            position='RA 119.2-119.9, Dec 13.3-16.9 (NSNS ridge, 3.0 deg E of the Silk spine)',
            nsns_peak_R_stack=tj['filaments']['E1']['stack_fit']['peak_R'],
            nsns_fwhm_deg=tj['filaments']['E1']['stack_fit']['fwhm_deg'],
            desi_oiii_peak_R=E['lines']['o3']['peak_R'], desi_oiii_err_boot=E['lines']['o3']['err_boot'],
            desi_oiii_peak_R_with_intercept=E['lines']['o3']['peak_R_with_intercept'],
            v_helio=vE['v'], v_helio_err_boot=vE['v_boot_std'], v_lsr=vE['v'] + LSR['E1'],
            oii_oiii_energy=E['ratios']['o2sum/o3_energy'], oii_twocomp_oiii_energy=E['ratios']['o2_twocomp/o3_energy'],
            hb_oiii_energy=E['ratios']['hb/o3_energy'], oiii_hb_lower2sig=E['ratios']['o3/hb_energy_lower2sig'],
            heii_oiii_energy=E['ratios']['he2/o3_energy'], ha_oiii_energy=E['ratios']['han2/o3_energy'],
            nii_oiii_energy=E['ratios']['nii6583/o3_energy'],
            sii_total_oiii_energy_value=E['ratios']['s2/o3_energy']['value'] + E['ratios']['s2_6731/o3_energy']['value'],
            oiii4363_oiii=E['ratios']['4363hg/o3_energy'],
            sdss_on_off_r010=dict(OIII_R=sj['results']['r0.10']['OIII:OIII5007']['I_R'],
                                  err=sj['results']['r0.10']['OIII:OIII5007']['err_boot'],
                                  placebo_rms=sj['placebo']['r0.10']['OIII:OIII5007']['rms'])),
        'S1': dict(position='RA 102.1-113.4, Dec 3.4-5.0 (Monogem south limb, Fesen+26 Fig 32)',
                   nsns_peak_R_stack=tj['filaments']['S1']['stack_fit']['peak_R'],
                   nsns_fwhm_deg=tj['filaments']['S1']['stack_fit']['fwhm_deg'],
                   desi_oiii_peak_R=S['lines']['o3']['peak_R'], desi_err_formal=S['lines']['o3']['err_formal'],
                   desi_null_rms=S['lines']['o3']['null_rms'], n_onridge=S['n_onridge_T>0.5'],
                   per_fibre_err_R=S['median_err_o3_R_per_fibre'], verdict='not measurable: backup-program stars'),
        'SILK_pipeline_check': dict(oiii_peak_R=K['lines']['o3']['peak_R'], err_boot=K['lines']['o3']['err_boot'],
                                    v_helio=vK['v'], v_helio_err_boot=vK['v_boot_std'],
                                    oii_oiii_energy=K['ratios']['o2sum/o3_energy'],
                                    canon=dict(oiii_R=1.62, v_helio_dec_lt19p5='+10.6 (10.2-12.0)',
                                               oii_oiii_south='1.11-1.16'))}
    json.dump(summ, open(f'{ROOT}/data/derived/w4_family_summary.json', 'w'), indent=1, default=float)

    fig = plt.figure(figsize=(16, 9.4))
    fig.text(0.012, 0.985, 'Monogem family test · DESI DR1 spectra of other periphery [O III] filaments',
             color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    fig.text(0.012, 0.955, 'Same residual-regression pipeline for all; Silk south re-measured as a check',
             color=FGDIM, fontsize=9.5, va='top')
    # (a) map
    d = np.load(f'{ROOT}/data/derived/w4_family_nsns_wide.npz')
    W = WCS(fits.Header.fromstring(str(d['wcs_header'])))
    m = d['oiii']
    bg = median_filter(np.nan_to_num(m, nan=np.nanmedian(m)), size=21)
    img = gaussian_filter(np.nan_to_num(m - bg), 0.7); img[~np.isfinite(m)] = np.nan
    ax = fig.add_axes([0.03, 0.08, 0.38, 0.82], projection=W)
    ax.imshow(img, origin='lower', cmap='gray', vmin=-0.2, vmax=1.2)
    tx = ax.get_transform('world')
    for (ra0, dec0, r, c) in [(104.15, 15.30, 11.0, WARN), (107.08, 13.50, 11.97, SUN)]:
        cr, cd = circle(ra0, dec0, r); ax.plot(cr, cd, color=c, lw=0.8, ls='--', transform=tx)
    fb = np.load(f'{ROOT}/data/derived/w4_family_desi_fibres.npz')
    for reg, col in [('south', DANGER), ('east', OK), ('silk', ACCENT)]:
        s = fb['region'] == reg
        ax.scatter(fb['ra'][s], fb['dec'][s], s=0.6, color=col, alpha=0.5, transform=tx, lw=0)
    for nm, (ra0, dec0, txt) in {'S1': (107.5, 2.6, 'S1 south limb'), 'E1': (120.0, 17.4, 'E1'),
                                 'Silk': (113.6, 22.5, 'Silk'), 'G206': (105.6, 10.6, 'G206.6+6.1'),
                                 'NW': (100.5, 21.0, 'NW group')}.items():
        label(ax, ra0, dec0, txt, size=8, transform=tx, ha='center')
    ax.set_xlim(60, 780); ax.set_ylim(150, 720)
    ax.coords[0].set_axislabel('RA'); ax.coords[1].set_axislabel('Dec')
    ax.coords.grid(color='#3a4450', alpha=0.4, lw=0.4)
    ax.set_title('NSNS [O III]; DESI coadds used (red S1, green E1, blue Silk)', color=FG, fontsize=9)
    # (b) stacks in velocity
    ax2 = fig.add_axes([0.47, 0.56, 0.24, 0.34])
    for nm, col, sc in [('SILK', ACCENT, 1.0), ('E1', OK, 4.0), ('S1', DANGER, 1.0)]:
        lam, sp = st[f'{nm}_lam'], st[f'{nm}_stack']
        v = (lam / L5007 - 1) * C
        mm = np.abs(v) < 260
        ax2.step(v[mm], sp[mm] * sc, where='mid', color=col, lw=1.0,
                 label={'SILK': 'Silk south', 'E1': 'E1 ×4', 'S1': 'S1'}[nm])
    ax2.axvline(vK['v'], color=ACCENT, lw=0.6, ls=':'); ax2.axvline(vE['v'], color=OK, lw=0.6, ls=':')
    ax2.set_xlabel('V_helio of [O III] 5007 (km/s)'); ax2.set_ylabel('regression stack (R/Å at peak)')
    ax2.legend(fontsize=7, ncol=3, loc='upper left'); ax2.grid(alpha=0.25)
    ax2.set_title('[O III] 5007 regression stacks', color=FG, fontsize=9)
    # (c) velocities
    ax3 = fig.add_axes([0.76, 0.56, 0.22, 0.34])
    rows = [('Silk south (canon)', 1.2, 1.0, ACCENT), ('Silk south (this check)', vK['v'] + LSR['SILK'], vK['v_boot_std'], ACCENT),
            ('Silk north (canon)', 10.6, 2.5, '#9ad4f5'), ('E1 (this work)', vE['v'] + LSR['E1'], np.hypot(vE['v_boot_std'], 1.5), OK)]
    for k, (lab, v, ev, col) in enumerate(rows):
        ax3.errorbar(v, k, xerr=ev, fmt='o', color=col, ms=5, capsize=2)
        ax3.text(-38, k + 0.25, lab, fontsize=7, color=FG)
    ax3.axvspan(-2, 5, color=FGDIM, alpha=0.15)
    ax3.text(-1.5, 3.55, 'local HI / WIM', fontsize=6.5, color=FGDIM)
    ax3.set_xlim(-40, 20); ax3.set_ylim(-0.6, 3.9); ax3.set_yticks([])
    ax3.set_xlabel('V_LSR (km/s)'); ax3.grid(alpha=0.25)
    ax3.set_title('[O III] velocity', color=FG, fontsize=9)
    # (d) ratios
    ax4 = fig.add_axes([0.47, 0.08, 0.51, 0.37])
    names = ['[O II]/[O III]', 'Hβ/[O III]', 'He II/[O III]', 'Hα/[O III]', '[N II]/[O III]', '[S II]6716/[O III]']
    keys = ['o2sum/o3_energy', 'hb/o3_energy', 'he2/o3_energy', 'han2/o3_energy', 'nii6583/o3_energy', 's2/o3_energy']
    x = np.arange(len(names))
    for k, kk in enumerate(keys):
        r = E['ratios'][kk]
        ax4.errorbar(k + 0.12, r['value'], yerr=[[r['value'] - r['p16']], [r['p84'] - r['value']]], fmt='o', color=OK,
                     ms=5, capsize=2, label='E1 (DESI)' if k == 0 else None)
        ax4.plot([k + 0.12] * 2, [r['p2_5'], r['p97_5']], color=OK, lw=0.6, alpha=0.6)
    r = K['ratios']['o2sum/o3_energy']
    ax4.errorbar(-0.12, r['value'], yerr=[[r['value'] - r['p16']], [r['p84'] - r['value']]], fmt='s', color=ACCENT,
                 ms=5, capsize=2, label='Silk south (DESI, this check)')
    r = K['ratios']['hb/o3_energy']
    ax4.errorbar(0.88, r['value'], yerr=[[r['value'] - r['p16']], [r['p84'] - r['value']]], fmt='s', color=ACCENT, ms=5, capsize=2)
    # canon references (Silk): north [O II]/[O III] 4.2-4.4 (DESI); He II < 0.026; Hb/[O III] < 0.09 (2 sigma)
    ax4.errorbar(-0.24, 4.3, yerr=[[0.8], [0.9]], fmt='D', color='#9ad4f5', ms=5, capsize=2, label='Silk north (canon)')
    ax4.plot([1.88, 2.12], [0.026, 0.026], color=ACCENT, lw=1.5); ax4.plot([0.75, 1.0], [0.09, 0.09], color=ACCENT, lw=1.5)
    ax4.text(1.6, 0.13, 'Silk 2σ limits', fontsize=6.5, color=ACCENT)
    ax4.axhline(0, color=FGDIM, lw=0.5)
    ax4.set_xticks(x); ax4.set_xticklabels(names, fontsize=7.5)
    ax4.set_yscale('symlog', linthresh=0.1); ax4.set_ylim(-0.3, 12); ax4.set_ylabel('energy ratio (16–84 %, thin 95 %)')
    ax4.legend(fontsize=7, ncol=3, loc='lower right'); ax4.grid(alpha=0.25)
    ax4.set_title('Line ratios relative to [O III] 5007', color=FG, fontsize=9)
    credit(fig)
    fig.savefig(f'{ROOT}/figures/w4_family_compare.png', dpi=200)
    print(json.dumps(summ, indent=1, default=float)[:3000])


if __name__ == '__main__':
    main()
