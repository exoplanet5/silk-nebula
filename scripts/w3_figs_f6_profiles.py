"""F6 cross-spine profiles of the southern Silk (spine Dec < 17.6) from data/derived/w3_figs_f6_profiles.npz/.json.

Row 1: [O III] in Rayleigh (NSNS DR0.2 median profile with its null band; DESI DR1 fibre residuals, petal-relative,
baseline-subtracted, mean +- s.e.). Other rows: profile minus its mean over the wings 0.5 <= |x| <= x_max, divided by the null
scatter sigma_null (median over x of the robust std of the mean-subtracted null profiles), i.e. in units of the
placement-null sigma; bands = 16-84 % and 2.5-97.5 % of the null profiles on the same scale.
Output: figures/w3_figs_f6_profiles.png
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, panel_tag

ROWS = [('fuv', 'GALEX FUV', OIII_C, 1.2, 'CU', '%.0f'),
        ('nuv', 'GALEX NUV', '#c792ea', 1.2, 'CU', '%.0f'),
        ('shassa', 'SHASSA Hα+[N II] · Dec 12–15.6', HA_C, 2.0, 'R', '%.2f'),
        ('erass', 'eRASS1 0.2–0.6 keV', '#ff9e80', 2.5, 'cts s⁻¹ deg⁻²', '%.2f'),
        ('hi', 'GALFA HI · V_LSR −2…+6', HB_C, 2.0, '10¹⁹ cm⁻²', '%.1f'),
        ('iris', 'IRIS 100 µm', SII_C, 2.0, 'MJy sr⁻¹', '%.2f')]
NULLDESC = {'fuv': 'RA shifts', 'nuv': 'RA shifts', 'shassa': 'random placements', 'erass': 'rotations about PSR',
            'hi': 'RA/Dec shifts', 'iris': 'RA/Dec shifts'}


def norm_rows(x, p, nl, xmax):
    m = np.abs(x) <= xmax + 1e-9
    off = m & (np.abs(x) >= 0.5)                 # reference level = mean of the off-ridge wings
    pm = p - np.nanmean(p[off])
    nm = nl - np.nanmean(nl[:, off], axis=1, keepdims=True)
    med = np.nanmedian(nm, 0)
    sig = 1.4826 * np.nanmedian(np.abs(nm - med), 0)
    s = float(np.nanmedian(sig[m]))
    return m, pm / s, nm / s, s


def main():
    P = np.load(DER / 'w3_figs_f6_profiles.npz')
    J = json.loads((DER / 'w3_figs_f6_profiles.json').read_text())['bands']
    fig, axs = plt.subplots(len(ROWS) + 1, 1, figsize=(9.2, 13.2), sharex=True,
                            gridspec_kw=dict(hspace=0.10, height_ratios=[1.45] + [1] * len(ROWS)))
    fig.text(0.012, 0.992, 'Silk Nebula · cross-spine profiles (south, spine Dec < 17.6)', color='#cfd6dd',
             fontsize=13, fontweight='bold', va='top')
    fig.text(0.012, 0.975, 'offset from the E spine, + = east (outside, away from Monogem); grey = null placements '
             '(16–84 %, 2.5–97.5 %)', color=FGDIM, fontsize=9, va='top')
    # ---- [O III] row in R
    ax = axs[0]
    x = P['nsns_x']; p = P['nsns']; nl = P['nsns_null']
    m = np.abs(x) <= 2.0
    base = np.nanmedian(p[m & (np.abs(x) > 0.6)])
    nb = nl - np.nanmedian(np.where((np.abs(x) > 0.6) & m, nl, np.nan), axis=1, keepdims=True)
    ax.fill_between(x[m], np.nanpercentile(nb, 2.5, 0)[m], np.nanpercentile(nb, 97.5, 0)[m], color=FGDIM, alpha=0.15, lw=0)
    ax.fill_between(x[m], np.nanpercentile(nb, 16, 0)[m], np.nanpercentile(nb, 84, 0)[m], color=FGDIM, alpha=0.25, lw=0)
    ax.plot(x[m], p[m] - base, color=OIII_C, lw=1.8, label='NSNS DR0.2 (2′ blocks, median)')
    xd = P['desi_x']; md = np.abs(xd) <= 2.0
    ax.errorbar(xd[md], P['desi'][md], yerr=P['desi_err'][md], fmt='o', ms=3.5, color='#ffffff', mec='k', mew=0.4,
                lw=0.8, capsize=0, label='DESI DR1 fibre residuals (GFA-corrected)')
    ax.set_ylabel('[O III] (R)')
    ax.legend(loc='upper right', fontsize=7.6, framealpha=0.85)
    panel_tag(ax, '[O III] 5007 · NSNS ridge %.2f R (z = %.0f) · DESI %.2f R' %
              (J['nsns']['A'], J['nsns']['z'], J['desi']['A']), size=8.6)
    ax.set_ylim(-0.6, 3.2)
    # ---- normalised rows
    for ax, (key, name, c, xmax, unit, fmt) in zip(axs[1:], ROWS):
        x = P[f'{key}_x']; p = P[key]; nl = np.asarray(P[f'{key}_null'], float)
        m, pz, nz, s = norm_rows(x, p, nl, xmax)
        ax.fill_between(x[m], np.nanpercentile(nz, 2.5, 0)[m], np.nanpercentile(nz, 97.5, 0)[m], color=FGDIM, alpha=0.15,
                        lw=0)
        ax.fill_between(x[m], np.nanpercentile(nz, 16, 0)[m], np.nanpercentile(nz, 84, 0)[m], color=FGDIM, alpha=0.25,
                        lw=0)
        ax.plot(x[m], pz[m], color=c, lw=1.6)
        b = J[key]
        txt = f'{name} · ridge {fmt % b["A"]} ± {fmt % b["A_null_sig"]} {unit} (z = {b["z"]:+.1f}, {b["n_null"]} {NULLDESC[key]})'
        panel_tag(ax, txt, size=8.0)
        ax.text(0.995, 0.06, f'σ_null = {fmt % s} {unit}', transform=ax.transAxes, ha='right', va='bottom',
                fontsize=7.3, color=FGDIM, path_effects=[HALO])
        ax.set_ylabel('(I − ⟨I⟩_wings)/σ_null', fontsize=7.5)
        lo = np.nanmin([np.nanpercentile(nz, 2.5, 0)[m].min(), np.nanmin(pz[m])])
        hi = np.nanmax([np.nanpercentile(nz, 97.5, 0)[m].max(), np.nanmax(pz[m])])
        ax.set_ylim(lo - 0.6, hi + 0.35 * (hi - lo) + 0.8)
    for ax in axs:
        ax.axvline(0, color=FGDIM, lw=0.7, ls=':')
        ax.axhline(0, color=FGDIM, lw=0.5, alpha=0.6)
        ax.grid(alpha=0.18)
        ax.tick_params(labelsize=8)
    axs[-1].set_xlim(-2.0, 2.0)
    axs[-1].set_xlabel('offset from E spine (deg)')
    axs[-1].text(0.0, -0.33, '← west · inside, toward Monogem', transform=axs[-1].transAxes, fontsize=8, color=FGDIM,
                 ha='left', va='top')
    axs[-1].text(1.0, -0.33, 'east · outside →', transform=axs[-1].transAxes, fontsize=8, color=FGDIM, ha='right',
                 va='top')
    ss.credit(fig)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.962, bottom=0.065)
    fig.savefig(FIG / 'w3_figs_f6_profiles.png', dpi=200)
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    main()
