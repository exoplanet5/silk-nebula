"""Maps, cross-filament profiles and on/off summaries of the SDSS sky-fibre nebular-line intensities
(topic tag spectra_sdss_lamost).  Reads the outputs of spectra_sdss_lamost_sdss_fit.py.

Outputs
  figures/spectra_sdss_lamost_sdss_oiii_map.png      per-fibre absolute [OIII]5007 over the field
  figures/spectra_sdss_lamost_sdss_profiles.png      line intensity vs distance from the nearest [OIII] strand
  data/derived/spectra_sdss_lamost_sdss_profiles.csv binned profiles
  data/derived/spectra_sdss_lamost_sdss_onoff.json   on/off-filament weighted means, excesses and ratios
"""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, NII_C, SII_C, FGDIM, FG, ACCENT, WARN
import spectra_sdss_lamost_trace as TR
import silk_geometry as sg

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
LINES = [('OIII', '[O III] 5007', OIII_C), ('Hb', 'Hβ', HB_C), ('Ha', 'Hα', HA_C), ('NII', '[N II] 6583', NII_C),
         ('SII', '[S II] 6716+6731', SII_C)]
ON, OFF = (0.0, 0.25), (0.6, 2.5)        # |distance| ranges (deg) for on / off strand samples


def wmean(v, e):
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if ok.sum() < 2:
        return np.nan, np.nan, int(ok.sum())
    w = 1 / e[ok] ** 2
    m = np.sum(w * v[ok]) / np.sum(w)
    # error: max of formal and scatter-based
    se = np.sqrt(1 / np.sum(w))
    chi = np.sum(w * (v[ok] - m) ** 2) / max(ok.sum() - 1, 1)
    return m, se * np.sqrt(max(chi, 1)), int(ok.sum())


def load():
    F = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_fibres.csv')
    G = Table.read(ROOT / 'data/derived/spectra_sdss_lamost_sdss_groups.csv')
    F['I_SII_cal'] = F['I_SII6716_cal'] + F['I_SII6731_cal']
    F['I_SII_cal_err'] = np.hypot(F['I_SII6716_cal_err'], F['I_SII6731_cal_err'])
    F['dI_SII_cal'] = F['dI_SII6716_cal'] + F['dI_SII6731_cal']
    geo = TR.load()
    d, inside, which = TR.nearest_strand_distance(np.array(F['ra']), np.array(F['dec']), geo)
    F['dist'] = d
    F['in_span'] = inside
    F['dE'], _ = TR.signed_distance(np.array(F['ra']), np.array(F['dec']), geo['E'])
    if geo.get('W') is not None:
        F['dW'], F['inW'] = TR.signed_distance(np.array(F['ra']), np.array(F['dec']), geo['W'])
    # group quality: drop bright (moon/twilight) skies and groups with poor [OIII] fits
    G['bad'] = (G['sb5040_ABmag'] < 20.4) | (G['I_OIII_cal_err'] > 0.6) | (G['n_sky'] < 6)
    bad = {(int(r['plate']), int(r['mjd']), int(r['sp'])) for r in G if r['bad']}
    keep = np.array([(int(p), int(m), int(s)) not in bad for p, m, s in zip(F['plate'], F['mjd'], F['sp'])])
    print('groups flagged bad:', sorted(bad))
    F = F[keep]
    return F, G, geo


def strand_fit(F, key, dec_bins=None, sigma=0.12, step=False):
    """Within-group fit  I_f = b_group + A_E G(dE) + A_W G(dW) [+ S * west-step of E]  (diffuse-corrected R).
    Group offsets absorb all plate-level systematics, so A measures the strand contrast at the ridge."""
    v = np.array(F[f'I_{key}_cal']) * np.array(F['f_psf'])
    e = np.array(F[f'I_{key}_cal_err']) * np.array(F['f_psf'])
    e = np.sqrt(e ** 2 + 0.05 ** 2)
    gid = np.array([f'{p}-{m}-{s}' for p, m, s in zip(F['plate'], F['mjd'], F['sp'])])
    dec = np.array(F['dec'])
    dE = np.array(F['dE'])
    dW = np.array(F['dW']) if 'dW' in F.colnames else np.full(len(F), 9.0)
    inW = np.array(F['inW']) if 'inW' in F.colnames else np.zeros(len(F), bool)
    ok = np.isfinite(v) & np.isfinite(e) & (np.array(F[f'I_{key}_cal_err']) < 5 * np.nanmedian(np.array(F[f'I_{key}_cal_err'])))
    ok &= (dec > 11.0) & (dec < 29.0)
    bins = dec_bins or [(11.0, 29.0)]
    gE = np.exp(-0.5 * (dE / sigma) ** 2)
    gW = np.where(inW, np.exp(-0.5 * (dW / sigma) ** 2), 0.0)
    stp = 1 / (1 + np.exp(dE / 0.15))           # ~1 west of strand E, 0 east
    out = []
    for lo, hi in bins:
        m = ok & (dec >= lo) & (dec < hi)
        nE0 = int(np.sum(m & (np.abs(dE) < 2 * sigma)))
        nW0 = int(np.sum(m & (np.abs(dW) < 2 * sigma) & inW))
        if m.sum() < 10 or (nE0 < 5 and nW0 < 5):
            out.append(dict(dec_lo=lo, dec_hi=hi, n=int(m.sum()), n_groups=0, rchi2=np.nan, n_near_E=nE0,
                            n_near_W=nW0, A_E=np.nan, A_E_err=np.nan, A_W=np.nan, A_W_err=np.nan,
                            S_west=np.nan, S_west_err=np.nan))
            continue
        for it in range(3):                      # iterative 5-sigma clipping
            ug = np.unique(gid[m])
            cols = [(gid[m] == u).astype(float) for u in ug]
            names = list(ug)
            cols.append(gE[m]); names.append('A_E')
            if gW[m].max() > 0.05:
                cols.append(gW[m]); names.append('A_W')
            if step:
                cols.append(stp[m]); names.append('S_west')
            A = np.vstack(cols).T / e[m][:, None]
            y = v[m] / e[m]
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            res = y - A @ coef
            rchi = np.sum(res ** 2) / max(len(y) - len(coef), 1)
            idx = np.where(m)[0]
            clipm = np.abs(res) > 5 * np.sqrt(rchi)
            if not clipm.any():
                break
            m[idx[clipm]] = False
        cov = np.linalg.pinv(A.T @ A) * max(rchi, 1)
        r = dict(dec_lo=lo, dec_hi=hi, n=int(m.sum()), n_groups=len(ug), rchi2=float(rchi),
                 n_near_E=int(np.sum(m & (np.abs(dE) < 2 * sigma))), n_near_W=int(np.sum(m & (np.abs(dW) < 2 * sigma) & inW)))
        for nm in ['A_E', 'A_W', 'S_west']:
            if nm in names:
                i = names.index(nm)
                r[nm], r[nm + '_err'] = float(coef[i]), float(np.sqrt(cov[i, i]))
            else:
                r[nm], r[nm + '_err'] = np.nan, np.nan
        out.append(r)
    return out


def onoff(F, key, use='I'):
    col = f'{use}_{key}_cal' if use == 'I' else f'dI_{key}_cal'
    v = np.array(F[col]) * np.array(F['f_psf'])            # diffuse-light corrected
    e = np.array(F[f'I_{key}_cal_err']) * np.array(F['f_psf'])
    good = np.isfinite(v) & np.isfinite(e) & (e < 5 * np.nanmedian(e)) & np.array(F['in_span'])
    ad = np.abs(np.array(F['dist']))
    on = good & (ad >= ON[0]) & (ad < ON[1])
    off = good & (ad >= OFF[0]) & (ad < OFF[1])
    m1, e1, n1 = wmean(v[on], e[on])
    m0, e0, n0 = wmean(v[off], e[off])
    return dict(on=m1, on_err=e1, n_on=n1, off=m0, off_err=e0, n_off=n0, excess=m1 - m0, excess_err=np.hypot(e1, e0))


def profiles(F, edges):
    rows = []
    ad = np.array(F['dist'])
    ins = np.array(F['in_span'])
    for a, b in zip(edges[:-1], edges[1:]):
        m = ins & (ad >= a) & (ad < b)
        r = dict(d_lo=a, d_hi=b, d_mid=0.5 * (a + b), n=int(m.sum()))
        for key, *_ in LINES:
            v = np.array(F[f'I_{key}_cal'])[m] * np.array(F['f_psf'])[m]
            e = np.array(F[f'I_{key}_cal_err'])[m] * np.array(F['f_psf'])[m]
            ok = np.isfinite(e) & (e < 5 * np.nanmedian(np.array(F[f'I_{key}_cal_err'])))
            mm, ee, nn = wmean(v[ok], e[ok])
            r[key], r[key + '_err'] = mm, ee
            vd = np.array(F[f'dI_{key}_cal'])[m] * np.array(F['f_psf'])[m]
            mm, ee, nn = wmean(vd[ok], e[ok])
            r['d' + key], r['d' + key + '_err'] = mm, ee
        rows.append(r)
    return Table(rows)


def fig_profiles(P, path, title_extra=''):
    fig, axes = new_sky_figure(11, 8.2, 'SDSS sky fibres · line intensity across the Silk filament',
                               'distance from nearest [O III] strand (E or W; + = east) · diffuse-corrected, Rayleigh',
                               nrows=3, ncols=2, sharex=True)
    axes = axes.ravel()
    x = np.array(P['d_mid'])
    for ax, (key, name, col) in zip(axes, LINES):
        y, e = np.array(P[key]), np.array(P[key + '_err'])
        ok = np.isfinite(y)
        ax.errorbar(x[ok], y[ok], e[ok], fmt='o', ms=3.5, color=col, ecolor=col, elinewidth=0.8, capsize=0)
        ax.axvspan(-ON[1], ON[1], color=ACCENT, alpha=0.08, lw=0)
        ax.axhline(0, color=FGDIM, lw=0.5)
        label(ax, 0.02, 0.88, name, color=col, size=9, transform=ax.transAxes)
        ax.set_ylabel('I (R)')
        ax.grid(alpha=0.25)
    # panel 6: number of fibres per bin
    ax = axes[5]
    ax.bar(x, np.array(P['n']), width=np.diff(np.array(P['d_lo']))[0] * 0.9, color=FGDIM, alpha=0.6)
    ax.set_ylabel('sky fibres / bin')
    label(ax, 0.02, 0.88, 'sampling', color=FG, size=9, transform=ax.transAxes)
    for ax in axes[4:]:
        ax.set_xlabel('distance from strand (deg)')
    finish(fig, path)


def _decorate(ax, geo, labels=True):
    if geo.get('envelope') is not None:
        ax.plot(geo['envelope'][:, 0], geo['envelope'][:, 1], color=ACCENT, lw=0.5, alpha=0.45)
    ax.plot(geo['E'][:, 0], geo['E'][:, 1], color=ACCENT, lw=1.1, alpha=0.9, label='[O III] strand E')
    if geo.get('W') is not None:
        ax.plot(geo['W'][:, 0], geo['W'][:, 1], color=ACCENT, lw=1.0, ls='--', alpha=0.9, label='strand W')
    if geo.get('halpha') is not None:
        hp = geo['halpha']
        ax.fill(hp[:, 0], hp[:, 1], facecolor='none', edgecolor=HA_C, lw=1.0, label='Hα patch')
    for k, (ra, dec) in sg.REGIONS.items():
        ax.plot(ra, dec, marker='s', mfc='none', mec=WARN, ms=6, mew=1.1)
        if labels:
            label(ax, ra - 0.1, dec + 0.05, k, color=WARN, size=8)
    ax.add_patch(Ellipse((116.2525, 15.895), 2 * 2.5 / np.cos(np.radians(15.9)), 5.0, fill=False, ec=FGDIM,
                         lw=0.7, ls=':'))
    ax.set_xlabel('RA (deg)')
    ax.grid(alpha=0.2)


def fig_map(F, G, geo, path, key='OIII', vmax=3.5):
    fig, axes = new_sky_figure(14, 9.6, 'SDSS sky fibres · absolute [O III] 5007',
                               'one dot per SDSS/BOSS sky fibre (DR17 + SDSS-V DR19), diffuse-corrected Rayleigh; fibres with σ < 0.6 R',
                               ncols=2, gridspec_kw=dict(width_ratios=[1.0, 1.25]))
    v = np.array(F[f'I_{key}_cal']) * np.array(F['f_psf'])
    e = np.array(F[f'I_{key}_cal_err']) * np.array(F['f_psf'])
    good = np.isfinite(v) & (e < 0.6)
    o = np.argsort(v[good])
    ra, dec = np.array(F['ra'])[good][o], np.array(F['dec'])[good][o]
    for ax, ms, lim in [(axes[0], 7, (123.3, 109.7, 4.7, 29.4)), (axes[1], 16, (119.4, 113.2, 13.0, 21.6))]:
        sc = ax.scatter(ra, dec, c=v[good][o], s=ms, cmap='magma', vmin=-0.3, vmax=vmax, lw=0)
        _decorate(ax, geo, labels=(ax is axes[1]))
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[2], lim[3])
        ax.set_aspect(1 / np.cos(np.radians(0.5 * (lim[2] + lim[3]))))
    axes[0].set_ylabel('Dec (deg)')
    axes[0].plot([119.4, 113.2, 113.2, 119.4, 119.4], [13.0, 13.0, 21.6, 21.6, 13.0], color=FG, lw=0.6, alpha=0.6)
    label(axes[1], 119.2, 13.25, 'LAMOST MRS-N field', color=FGDIM, size=7)
    axes[1].legend(loc='upper right', fontsize=7)
    fig.subplots_adjust(left=0.07, right=0.86, top=0.91, bottom=0.07, wspace=0.12)
    cb = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.015)
    cb.set_label('[O III] 5007 (R)')
    finish(fig, path, tight=False)


def fig_strand(SF, path):
    """Strand amplitude A_E (and A_W) per Dec bin for each line."""
    fig, ax = new_sky_figure(10, 5.6, 'SDSS sky fibres · strand contrast along the filament',
                             'within-plate Gaussian amplitude at the ridge (σ = 0.12°), diffuse-corrected')
    off = {'OIII': 0.0, 'Hb': 0.12, 'Ha': -0.12, 'NII': 0.24, 'SII': -0.24}
    for key, name, col in LINES:
        rows = SF[key]['dec_bins']
        x = np.array([0.5 * (r['dec_lo'] + r['dec_hi']) for r in rows]) + off[key]
        y = np.array([r['A_E'] for r in rows]); e = np.array([r['A_E_err'] for r in rows])
        ok = np.isfinite(y) & (e < 3)
        ax.errorbar(x[ok], y[ok], e[ok], fmt='o', ms=4, color=col, ecolor=col, elinewidth=0.9, capsize=0, label=name)
        w = SF[key]['W_strand'][0]
        if np.isfinite(w['A_W']):
            ax.errorbar([17.6 + off[key]], [w['A_W']], [w['A_W_err']], fmt='D', ms=5, mfc='none', color=col, ecolor=col,
                        elinewidth=0.9)
    ax.axhline(0, color=FGDIM, lw=0.6)
    label(ax, 17.25, ax.get_ylim()[1] * 0.9 if ax.get_ylim()[1] > 0 else 1, '◇ = W strand', color=FGDIM, size=8)
    ax.set_xlabel('Dec (deg)')
    ax.set_ylabel('strand amplitude (R)')
    ax.grid(alpha=0.25)
    ax.legend(ncol=5, fontsize=7.5, loc='lower center', bbox_to_anchor=(0.5, -0.32))
    finish(fig, path)


def main():
    F, G, geo = load()
    edges = np.arange(-3.0, 3.01, 0.2)
    P = profiles(F, edges)
    P.write(ROOT / 'data/derived/spectra_sdss_lamost_sdss_profiles.csv', overwrite=True)
    res = {'on_range_deg': ON, 'off_range_deg': OFF, 'trace_version': geo['version'],
           'units': 'Rayleigh, diffuse-light corrected (x f_psf); sky fibres with Dec inside the strand span',
           'lines': {}}
    for key, *_ in LINES:
        res['lines'][key] = onoff(F, key)
        res['lines'][key + '_residual'] = onoff(F, key, use='d')
    L = res['lines']
    for a, b in [('OIII', 'Hb'), ('OIII', 'Ha'), ('SII', 'Ha'), ('NII', 'Ha')]:
        ea, eb = L[a]['excess'], L[b]['excess']
        res[f'excess_ratio_{a}/{b}'] = ea / eb if eb != 0 else np.nan
        res[f'excess_ratio_{a}/{b}_note'] = f'{a} excess {ea:.2f}+-{L[a]["excess_err"]:.2f} R, {b} excess {eb:.2f}+-{L[b]["excess_err"]:.2f} R'
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_onoff.json').write_text(json.dumps(res, indent=1, default=float))
    bins = [(11.0, 14.0), (14.0, 15.5), (15.5, 17.0), (17.0, 18.5), (18.5, 20.0), (20.0, 23.0), (23.0, 29.0)]
    SF = {}
    for key, *_ in LINES:
        SF[key] = dict(all=strand_fit(F, key, step=key in ('Ha', 'NII', 'SII', 'Hb')),
                       dec_bins=strand_fit(F, key, dec_bins=bins),
                       W_strand=strand_fit(F, key, dec_bins=[(15.0, 20.0)]))
    res['strand_fit'] = SF
    res['strand_fit_note'] = ('I_f = b_group + A_E*G(dE) + A_W*G(dW) [+ S_west*step] with G Gaussian sigma 0.12 deg; '
                              'per plate-spectrograph offsets b_group; diffuse-corrected R; 5-sigma clipping')
    (ROOT / 'data/derived/spectra_sdss_lamost_sdss_onoff.json').write_text(json.dumps(res, indent=1, default=float))
    for key in SF:
        a = SF[key]['all'][0]
        print(key, 'A_E %.2f+-%.2f A_W %.2f+-%.2f S_west %.2f+-%.2f n %d' % (a['A_E'], a['A_E_err'], a['A_W'], a['A_W_err'], a['S_west'], a['S_west_err'], a['n']))
        print('   ', ['%.1f-%.1f: %.2f+-%.2f (n=%d,nE=%d)' % (r['dec_lo'], r['dec_hi'], r['A_E'], r['A_E_err'], r['n'], r['n_near_E']) for r in SF[key]['dec_bins']])
    fig_strand(SF, ROOT / 'figures/spectra_sdss_lamost_sdss_strand_dec.png')
    fig_profiles(P, ROOT / 'figures/spectra_sdss_lamost_sdss_profiles.png')
    fig_map(F, G, geo, ROOT / 'figures/spectra_sdss_lamost_sdss_oiii_map.png')


if __name__ == '__main__':
    main()
