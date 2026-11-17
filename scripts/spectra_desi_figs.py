"""spectra_desi_figs.py -- figures for the DESI residual-emission analysis.

Usage: python -I scripts/spectra_desi_figs.py [map|profile|regspec|velocity|brightness|ratios|controls|all]
Reads data/derived/spectra_desi_fiber_table.csv, spectra_desi_results.json, spectra_desi_model_arrays.npz,
spectra_desi_regspec.npz, spectra_desi_regspec_lines.json, spectra_desi_fits.npz.
"""
import csv
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, FG, FGDIM, OIII_C, HA_C, SII_C, NII_C, HB_C, WARN,  # noqa
                        SUN, OK, ACCENT, BG2)
from silk_geometry import ROOT, REGIONS, load_trace, V_HELIO  # noqa: E402
from spectra_desi_fig_coverage import draw_trace, draw_regions  # noqa: E402
from spectra_desi_lines import LINES, WIN_WAVE  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

D = ROOT / 'data/derived'
FIG = ROOT / 'figures'


def fibtab():
    rows = list(csv.DictReader(open(D / 'spectra_desi_fiber_table.csv')))
    t = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    for k in rows[0]:
        if k not in ('sparcl_id', 'survey', 'program', 'spectype'):
            t[k] = t[k].astype(float)
    t['ok'] = t['ok'].astype(bool)
    return t


def sky_axes(ax, xlim, ylim, dec0):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect(1 / np.cos(np.radians(dec0)))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.grid(alpha=0.2, lw=0.5)


def binned_map(ra, dec, v, e, xe, ye):
    w = 1 / e ** 2
    sw, _, _ = np.histogram2d(ra, dec, [xe, ye], weights=w)
    swv, _, _ = np.histogram2d(ra, dec, [xe, ye], weights=w * v)
    n, _, _ = np.histogram2d(ra, dec, [xe, ye])
    with np.errstate(invalid='ignore', divide='ignore'):
        m = swv / sw
    m[n < 3] = np.nan
    return m.T, (1 / np.sqrt(sw)).T


def fig_map(xlim=(119.2, 113.6), ylim=(13.3, 28.4), tag='', binsz=0.15, vlim=2.0, w=14.5, h=11.5):
    t = fibtab()
    tr = load_trace()
    ok = t['ok'] & np.isfinite(t['OIII_R'])
    ra, dec, v, e = t['ra'][ok], t['dec'][ok], t['OIII_R'][ok], t['OIII_R_err'][ok]
    dec0 = np.mean(ylim)
    fig, axes = new_sky_figure(w, h, 'DESI DR1 · [O III] 5007 sky-subtraction residual',
                               f'{ok.sum():,} galaxy/QSO coadds · in-fibre surface brightness (R) · '
                               'fibre minus petal-mean sky', ncols=3)
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
    a = axes[0]
    o = np.argsort(np.abs(v))
    sc = a.scatter(ra[o], dec[o], c=v[o], s=4, cmap='RdBu_r', norm=norm, lw=0)
    draw_trace(a, tr, lw=0.8, alpha=0.7, envelope=False, color='#d8f3ff')
    draw_regions(a, size=6)
    sky_axes(a, xlim, ylim, dec0)
    a.set_title('per fibre', color=FG)
    xe = np.arange(xlim[1], xlim[0] + 1e-6, binsz / np.cos(np.radians(dec0)))
    ye = np.arange(ylim[0], ylim[1] + 1e-6, binsz)
    m, _ = binned_map(ra, dec, v, e, xe, ye)
    a = axes[1]
    a.imshow(m, origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]), cmap='RdBu_r', norm=norm,
             aspect='auto', interpolation='nearest')
    draw_trace(a, tr, lw=0.8, alpha=0.7, envelope=False, color='#d8f3ff')
    sky_axes(a, xlim, ylim, dec0)
    a.set_title(f'binned {binsz}° (ivar mean)', color=FG)
    a = axes[2]
    p = t['model_resid_R'][ok]
    mm = np.isfinite(p)
    mp, _ = binned_map(ra[mm], dec[mm], p[mm], np.ones(mm.sum()), xe, ye)
    a.imshow(mp, origin='lower', extent=(xe[0], xe[-1], ye[0], ye[-1]), cmap='RdBu_r', norm=norm,
             aspect='auto', interpolation='nearest')
    draw_trace(a, tr, lw=0.8, alpha=0.7, envelope=False, color='#d8f3ff')
    sky_axes(a, xlim, ylim, dec0)
    a.set_title('ridge model after petal subtraction', color=FG)
    cb = fig.colorbar(sc, ax=list(axes), fraction=0.02, pad=0.01)
    cb.set_label('[O III] residual (R)', color=FG)
    axes[0].legend(loc='lower left', fontsize=7)
    finish(fig, FIG / f'spectra_desi_oiii_map{tag}.png', tight=False)


def fig_profile(dec_bins=((16.1, 16.7),), tag=''):
    t = fibtab()
    ok = t['ok'] & np.isfinite(t['OIII_R']) & np.isfinite(t['dRA_spine'])
    nb = len(dec_bins)
    fig, axes = new_sky_figure(9.5, 2.6 * nb + 1.2, 'DESI [O III] residual across the filament',
                               'ivar-weighted means in 0.05° bins of RA offset from the E spine · '
                               'line = ridge model after petal subtraction', nrows=nb, squeeze=False)
    for (lo, hi), a in zip(dec_bins, axes[:, 0]):
        m = ok & (t['dec'] >= lo) & (t['dec'] < hi)
        x = t['dRA_spine'][m]
        y, e, p = t['OIII_R'][m], t['OIII_R_err'][m], t['model_resid_R'][m]
        be = np.arange(-2.5, 2.5001, 0.05)
        idx = np.digitize(x, be) - 1
        xc, ym, em, pm = [], [], [], []
        for k in range(len(be) - 1):
            s = idx == k
            if s.sum() >= 3:
                w = 1 / e[s] ** 2
                xc.append(0.5 * (be[k] + be[k + 1]))
                ym.append(np.sum(w * y[s]) / w.sum())
                em.append(w.sum() ** -0.5)
                pp = p[s][np.isfinite(p[s])]
                pm.append(pp.mean() if pp.size else np.nan)
        a.axhline(0, color=FGDIM, lw=0.6)
        a.errorbar(xc, ym, em, fmt='o', ms=3, color=OIII_C, ecolor=OIII_C, elinewidth=0.8, label='DESI fibres')
        a.plot(xc, pm, '-', color=SUN, lw=1.4, label='ridge model')
        a.set_xlim(2.5, -2.5)
        a.set_ylabel('[O III] residual (R)')
        label(a, 0.02, 0.88, f'Dec {lo:.1f}–{hi:.1f}  n={m.sum()}', transform=a.transAxes, size=8.5)
        a.grid(alpha=0.2)
    axes[-1, 0].set_xlabel('ΔRA·cosδ from E spine (deg; east → left)')
    axes[0, 0].legend(loc='upper right', fontsize=7.5)
    finish(fig, FIG / f'spectra_desi_profile{tag}.png')


PANELS = [('[O II] 3727', 3712, 3745, ['OII3726', 'OII3729']), ('Hβ', 4842, 4884, ['Hb']),
          ('[O III] 4959,5007', 4945, 5025, ['OIII4959', 'OIII5007']), ('He II 4686', 4672, 4702, ['HeII4686']),
          ('[N I] 5200', 5185, 5215, ['NI5198', 'NI5200']), ('[O I] 6300', 6287, 6318, ['OI6300']),
          ('Hα + [N II]', 6535, 6600, ['NII6548', 'Ha', 'NII6583']), ('[S II]', 6705, 6745, ['SII6716', 'SII6731']),
          ('He I 5876', 5862, 5892, ['HeI5876']), ('[O I] 5577 sky', 5562, 5595, ['sky5577'])]


def fig_regspec(inst='model', tag=''):
    r = np.load(D / 'spectra_desi_regspec.npz', allow_pickle=False)
    fits = json.load(open(D / 'spectra_desi_regspec_lines.json'))
    S, eS = r[f'S_{inst}'], r[f'eS_{inst}']
    fig, axes = new_sky_figure(13.5, 7.6, 'DESI regression spectrum of the [O III] filament',
                               f'Σ w P c / Σ w P² with P = expected [O III] residual ({inst}); '
                               f'{fits[inst]["n_fibres"]:,} fibres; units: per unit [O III] 5007', nrows=2, ncols=5)
    for (ttl, lo, hi, lines), a in zip(PANELS, axes.ravel()):
        s = (WIN_WAVE >= lo) & (WIN_WAVE <= hi)
        a.fill_between(WIN_WAVE[s], (S - eS)[s], (S + eS)[s], step='mid', color=OIII_C, alpha=0.25, lw=0)
        a.step(WIN_WAVE[s], S[s], where='mid', color=OIII_C, lw=0.9)
        a.axhline(0, color=FGDIM, lw=0.5)
        for ln in lines:
            l0 = LINES[ln] * (1 + 20 / 299792.458)
            a.axvline(l0, color=SUN, lw=0.5, alpha=0.6, ls=':')
        a.set_title(ttl, color=FG, fontsize=9)
        a.set_xlim(lo, hi)
        a.tick_params(labelsize=7)
        a.grid(alpha=0.15)
    for a in axes[:, 0]:
        a.set_ylabel('per Å / [O III]')
    for a in axes[1]:
        a.set_xlabel('vacuum λ (Å, barycentric)')
    finish(fig, FIG / f'spectra_desi_regspec_{inst}{tag}.png')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('map', 'all'):
        fig_map()
    if what in ('pilotmap',):
        fig_map(xlim=(118.8, 114.4), ylim=(13.8, 19.2), tag='_zoom', binsz=0.1, vlim=2.0, w=14, h=7.0)
    if what in ('profile', 'all', 'pilotprofile'):
        fig_profile(tag='_pilot' if what == 'pilotprofile' else '')
    if what in ('regspec', 'all'):
        fig_regspec('model')
        fig_regspec('nbr')
