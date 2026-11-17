"""Figure + summary for the Pan-STARRS1 DR2 stack test (reads imaging_optical_ps1_profiles.json).

(a) cross-filament profiles (star-masked, 1' medians, 2' bins, quadratic background) in g and r at P3, P1, P2 and the
    northern arc (Dec 22.3), with the SDSS-scaled [OIII] expectation (1.35 x NSNS template x 1.75e-3 nmgy/arcsec^2 per R
    in PS1 g, imaging_optical_expected.py);
(b) Abell 21 (Medusa PN, ~10' diameter) sky-subtraction test: PS1 radial profile vs the profile predicted from the NSNS
    DR0.2 [OIII] and Halpha maps (PS1 g <- [OIII]5007+4959, r <- Halpha+[NII]; [NII]/Halpha assumed 1.0 here).
Flux unit: 10^-3 nmgy/arcsec^2 (1 nmgy/arcsec^2 = 22.5 AB mag/arcsec^2).
Outputs: figures/imaging_optical_ps1_profiles.png, data/derived/imaging_optical_ps1_summary.json
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from imaging_optical_common import DER, FIG, save_json
import silk_style as ss

NM = 10 ** (-0.4 * 22.5)          # 1 nmgy/arcsec^2 in 10^(-0.4 m) units
PER_R = {'g': 1.75e-3, 'r': 1.88e-3}  # nmgy/arcsec^2 per R ([OIII]5007+4959 in g; Halpha in r)


def medusa_prediction():
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    w = WCS(str(z['wcs_header']))
    ra, dec = z['ra'], z['dec']
    r = np.hypot((ra - 112.2625) * np.cos(np.radians(13.2464)), dec - 13.2464) * 60
    out = {}
    for ch in ('oiii', 'halpha'):
        img = z[ch]
        bg = np.nanmedian(img[(r > 15) & (r < 25)])
        edges = np.arange(0, 13.01, 1.0)
        out[ch] = [float(np.nanmedian(img[(r >= edges[i]) & (r < edges[i + 1])]) - bg) for i in range(len(edges) - 1)]
        out['r_arcmin'] = (0.5 * (edges[1:] + edges[:-1])).tolist()
    return out


def main():
    d = json.loads((DER / 'imaging_optical_ps1_profiles.json').read_text())
    sdss = json.loads((DER / 'imaging_optical_sdss_sky.json').read_text())['cases']['g_vs_OIII_ridge_nostars']['core_dec12_20p5']
    scale = sdss['A_over_expected']
    med = medusa_prediction()
    fig = plt.figure(figsize=(14.5, 7.4))
    fig.text(0.012, 0.985, 'Silk Nebula · Pan-STARRS1 stacks', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    fig.text(0.012, 0.948, 'star-masked 1′ medians across the ridge · DR2 stack.unconv, sky-subtracted', color=ss.FGDIM, fontsize=9.5, va='top')
    strips = [('P3', 14.72), ('P1', 16.60), ('P2', 17.74), ('N22', 22.3)]
    summ = {'strips': {}}
    for i, (nm, dec) in enumerate(strips):
        ax = fig.add_axes([0.05 + (i % 2) * 0.33, 0.55 - (i // 2) * 0.45, 0.29, 0.33])
        for b, col in (('g', ss.OIII_C), ('r', ss.HA_C)):
            v = d.get(f'{nm}_{b}')
            if v is None:
                continue
            x = np.array(v['across_deg']); y = np.array(v['profile_bgsub']) / NM * 1e3
            ax.plot(x, y, '-', color=col, lw=1.1, label=f'PS1 {b}')
            A = v['A_per_R'] / NM; eA = v['A_err_placebo'] / NM
            summ['strips'][f'{nm}_{b}'] = dict(A_nmgy_per_NSNS_R=A, err=eA, A_over_expected=A / PER_R[b], err_over_expected=eA / PER_R[b],
                                               bin_scatter_1e3nmgy=float(10 ** (-0.4 * (v['bin_scatter_mag'] - 22.5)) * 1e3),
                                               exptime_median_s=v['exptime_median'], n_cutouts=v['n_cutouts'])
            if b == 'g':
                t = np.array(v['template_R'])
                bgm = np.isfinite(t) & (np.abs(x) >= 0.7)
                t3 = t - np.polyval(np.polyfit(x[bgm], t[bgm], 2), x)
                ax.plot(x, scale * PER_R['g'] * t3 * 1e3, '--', color='#ffffff', lw=1.0, label=f'SDSS-scaled [OIII] (×{scale:.2f})')
        ax.axhline(0, color='#2a323c', lw=0.8)
        ax.set_xlim(-1.7, 1.7); ax.set_ylim(-4, 8)
        ax.set_title(f'{nm}  Dec {dec:+.2f}', color=ss.FG, fontsize=9)
        ax.set_xlabel('distance from [OIII] ridge [deg]', fontsize=8); ax.set_ylabel('Δ [10⁻³ nmgy/″²]', fontsize=8)
        ax.grid(alpha=0.3); ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc='upper left')
    # Medusa test
    ax = fig.add_axes([0.72, 0.10, 0.26, 0.78])
    rm = np.array(med['r_arcmin'])
    for b, col, pred in (('g', ss.OIII_C, PER_R['g'] * np.array(med['oiii'])),
                         ('r', ss.HA_C, PER_R['r'] * 2.0 * np.array(med['halpha']))):
        v = d['medusa_A21'][b]
        r = np.array(v['r_arcmin']); y = np.array(v['profile_flux_arcsec2']) / NM
        y = y - np.nanmedian(y[r > 10])
        ax.plot(r, y * 1e3, '-', color=col, lw=1.3, label=f'PS1 {b} observed')
        ax.plot(rm, pred * 1e3, '--', color=col, lw=1.1, label=f'{b} predicted from NSNS')
        summ[f'medusa_{b}'] = dict(r_arcmin=r.tolist(), observed_1e3nmgy=(y * 1e3).tolist(),
                                   predicted_r_arcmin=rm.tolist(), predicted_1e3nmgy=(pred * 1e3).tolist())
    ax.axhline(0, color='#2a323c', lw=0.8)
    ax.set_xlabel('radius from Abell 21 [arcmin]', fontsize=8); ax.set_ylabel('Δ [10⁻³ nmgy/″²]', fontsize=8)
    ax.set_title('Medusa PN: extended-emission response', color=ss.FG, fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, loc='upper right'); ax.tick_params(labelsize=7)
    ax.set_xlim(0, 13)
    ss.credit(fig)
    fig.savefig(FIG / 'imaging_optical_ps1_profiles.png', dpi=200)
    plt.close(fig)
    summ['medusa_nsns_profile_R'] = med
    summ['notes'] = ('A = fitted amplitude of the NSNS template (g: [OIII], r: Halpha) per R, in nmgy/arcsec^2; '
                     'A_over_expected = A / (response to 1 R of the line). PS1 DR2 stacks are background-subtracted on '
                     'arcminute scales, so these amplitudes are not sensitive to the degree-scale filament (see Medusa test).')
    save_json(summ, DER / 'imaging_optical_ps1_summary.json')
    for k, v in summ['strips'].items():
        print(k, round(v['A_over_expected'], 2), '+-', round(v['err_over_expected'], 2), 'binsc', round(v['bin_scatter_1e3nmgy'], 2))
    print('medusa g obs', [round(x, 1) for x in summ['medusa_g']['observed_1e3nmgy'][:10]])
    print('medusa g pred', [round(x, 1) for x in summ['medusa_g']['predicted_1e3nmgy'][:8]])
    print('medusa r obs', [round(x, 1) for x in summ['medusa_r']['observed_1e3nmgy'][:10]])
    print('medusa r pred', [round(x, 1) for x in summ['medusa_r']['predicted_1e3nmgy'][:8]])


if __name__ == '__main__':
    main()
