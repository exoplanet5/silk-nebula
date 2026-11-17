"""verify_xrayfuv: GALEX FUV/NUV on the Silk vs shock / interface UV-line predictions."""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *
from verify_xrayfuv_common import *
import verify_xrayfuv_galex_measure as M
from scipy import ndimage

meas = json.loads((DER / 'verify_xrayfuv_galex_measure.json').read_text())
pred = json.loads((DER / 'verify_xrayfuv_uv_predictions.json').read_text())
dust = json.loads((DER / 'verify_xrayfuv_galex_dust.json').read_text())
ra_sp, dec_sp, _ = load_spine()
fig = plt.figure(figsize=(15, 10.5))
fig.text(0.012, 0.985, 'Silk Nebula · GALEX FUV/NUV excess vs UV-line predictions', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
fig.text(0.012, 0.958, 'skybg 1′ blocks, tile offsets + mean flat pattern removed, Tycho-2/NSNS star masks · predictions per Rayleigh of [O III] 5007, E(B−V)=0.03 (CCM89)', color=FGDIM, fontsize=9.5, va='top')
gs = fig.add_gridspec(1, 3, left=0.05, right=0.955, top=0.92, bottom=0.555, wspace=0.30)
gs2 = fig.add_gridspec(1, 3, left=0.20, right=0.985, top=0.45, bottom=0.07, wspace=0.22, width_ratios=[1, 0.8, 1.1])
Dd = {}
for k, band in enumerate(['FUV', 'NUV']):
    D = M.prep_band(band); Dd[band] = D
    ax = fig.add_subplot(gs[0, k])
    rae = np.arange(114.0, 119.01, 1 / 30); dee = np.arange(11.5, 20.01, 1 / 30)
    S, _, _ = np.histogram2d(D['dec'], D['ra'], bins=[dee, rae], weights=D['v'])
    N, _, _ = np.histogram2d(D['dec'], D['ra'], bins=[dee, rae])
    img = np.where(N > 0, S / np.maximum(N, 1), np.nan)
    f = np.where(np.isfinite(img), img, 0); w = np.isfinite(img).astype(float)
    sm = ndimage.gaussian_filter(f, 2) / np.maximum(ndimage.gaussian_filter(w, 2), 1e-3); sm[w == 0] = np.nan
    im = ax.imshow(sm, origin='lower', extent=[rae[0], rae[-1], dee[0], dee[-1]], cmap='RdBu_r', vmin=-120, vmax=120, aspect='auto')
    ax.invert_xaxis()
    T = np.load(DER / 'imaging_optical_nsns_maps.npz')
    o3 = ndimage.gaussian_filter(np.nan_to_num(T['oiii']), 3)
    ax.contour(T['ra'], T['dec'], o3, levels=[0.8, 1.6], colors=[ACCENT, '#b8ecff'], linewidths=0.7)
    ax.set_xlim(119, 114); ax.set_ylim(11.5, 20)
    ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01); cb.set_label('CU (relative)', fontsize=8)
    ax.set_title(f'({"ab"[k]}) GALEX {band} residual, NSNS [O III] contours 0.8/1.6 R', fontsize=9.5, loc='left')

# (c) cross-spine profiles (south)
ax = fig.add_subplot(gs[0, 2])
ns = np.load(DER / 'imaging_optical_nsns_maps.npz')
cc = CurveCoords(ra_sp, dec_sp)
edges = np.arange(-1.0, 1.0001, 0.05); xc = 0.5 * (edges[1:] + edges[:-1])
for band, col in [('FUV', OIII_C), ('NUV', '#c792ea')]:
    D = Dd[band]
    x, s = cc.coords(radec_to_vec(D['ra'], D['dec']), maxdist_deg=1.2)
    sel = np.isfinite(x) & (D['dec'] < 17.6)
    v = D['v'].copy()
    for t in np.unique(D['tile'][sel]):     # demean per tile using |x|>0.3 blocks
        m = (D['tile'] == t) & sel
        ref = m & (np.abs(x) > 0.3)
        if ref.sum() > 20:
            v[m] -= np.mean(v[ref])
        else:
            sel &= ~m
    bi = np.digitize(x[sel], edges) - 1; ok = (bi >= 0) & (bi < len(xc))
    prof = np.bincount(bi[ok], v[sel][ok], len(xc)) / np.maximum(np.bincount(bi[ok], minlength=len(xc)), 1)
    ax.plot(xc, prof, color=col, lw=1.6, label=f'GALEX {band}')
    print(band, 'profile mean |x|<0.1:', np.mean(prof[np.abs(xc) < 0.1]), 'mean |x|<0.25:', np.mean(prof[np.abs(xc) < 0.25]))
    if band == 'FUV':
        t3 = M.sample_template(np.nan_to_num(ns['oiii']), D['ra'], D['dec'])
        tt = t3.copy()
        for t in np.unique(D['tile'][sel]):
            m = (D['tile'] == t) & sel; ref = m & (np.abs(x) > 0.3)
            tt[m] -= np.nanmean(t3[ref])
        pt = np.bincount(bi[ok], np.nan_to_num(tt[sel][ok]), len(xc)) / np.maximum(np.bincount(bi[ok], minlength=len(xc)), 1)
print('NSNS profile mean |x|<0.1:', np.mean(pt[np.abs(xc) < 0.1]), '|x|<0.25:', np.mean(pt[np.abs(xc) < 0.25]))
ax2 = ax.twinx(); ax2.plot(xc, pt, color=WARN, lw=1.2, ls='--', label='NSNS [O III] (R, right)')
ax2.set_ylabel('NSNS [O III] excess (R)', color=WARN); ax2.tick_params(colors=WARN)
ax.axhline(0, color=FGDIM, lw=0.6); ax.axvline(0, color=FGDIM, lw=0.6, ls=':')
ax.set_xlabel('offset from E spine (deg, + east)'); ax.set_ylabel('UV excess (CU)')
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=7.5)
ax.set_title('(c) stacked cross-spine profiles, Dec < 17.6', fontsize=9.5, loc='left')

# (d,e) predicted vs measured per Rayleigh
groups = [('R79 complete', ['R79_E', 'R79_F', 'R79_G', 'R79_H', 'R79_I']),
          ('CR85 complete', ['CR85_run1', 'CR85_run3', 'CR85_run5']),
          ('DS96 150/200', ['DS96_s15010', 'DS96_s20010']),
          ('R88 incomplete', ['R88_SE', 'R88_M'])]
for k, band in enumerate(['FUV', 'NUV']):
    ax = fig.add_subplot(gs2[0, k])
    y = 0; yt = []; yl = []
    for gname, names in groups:
        for nm in names:
            m = pred['models'][nm]; p = m['per_R']
            tot = p['ebv0.03_civsupp1.0'][band]['total']; s35 = p['ebv0.03_civsupp3.5'][band]['total']
            s20 = p['ebv0.03_civsupp20.0'][band]['total']; nr = p['ebv0.03_civsupp1.0'][band]['non_resonant']
            ax.plot([nr, tot], [y, y], color=FGDIM, lw=1)
            ax.plot(tot, y, 'o', color=DANGER, ms=5); ax.plot(s35, y, 'o', color=WARN, ms=4); ax.plot(s20, y, 'o', color=SUN, ms=4)
            ax.plot(nr, y, 's', color=OK, ms=4)
            yt.append(y); yl.append(f"{nm.replace('_', ' ')} {m['vs']} km/s [O III]/Hβ {m['5007/Hb']:.0f}"); y += 1
        y += 0.6
    for ep in ['2.2e6yr', '7.4e6yr']:
        r = pred['interfaces'][f'BBF90_{ep}']
        ax.plot([r[f'{band}_nonres_CU_per_R5007_min'], r[f'{band}_CU_per_R5007_min']], [y, y], color=FGDIM, lw=1)
        ax.plot(r[f'{band}_CU_per_R5007_min'], y, '>', color='#c792ea', ms=6)
        yt.append(y); yl.append(f'BBF90 interface {ep} (≥, needs ×{r["limb_brightening_needed_for_1R"]:.0f})'); y += 1
    jb = meas[band]['joint_south_dec<17.6']; sb = meas[band]['beta_OIII_south_dec<17.6']; ab = meas[band]['beta_OIII_all']
    for val, err, col, lab in [(sb['beta_CU_per_R'], sb['null_std'], ACCENT, 'measured β south (per NSNS R)'),
                               (jb['beta_OIII'], jb['null_std_OIII'], '#7fdcff', 'joint fit south'),
                               (ab['beta_CU_per_R'] / 1.35, ab['null_std'] / 1.35, FG, 'all, per SDSS-scaled R (÷1.35)')]:
        ax.axvspan(val - err, val + err, color=col, alpha=0.18); ax.axvline(val, color=col, lw=1.3, label=lab)
    ax.set_yticks(yt); ax.set_yticklabels(yl if k == 0 else [], fontsize=6.8); ax.invert_yaxis()
    ax.set_xscale('log'); ax.set_xlim(20, 5000 if band == 'FUV' else 1500)
    ax.set_xlabel(f'GALEX {band} per 1 R of [O III] 5007 (CU)')
    ax.legend(loc='upper right', fontsize=6.6)
    ax.set_title(f'({"de"[k]}) {band}: ● C IV ×1 ● ÷3.5 ● ÷20 ■ non-res.' if k == 0 else f'(e) {band}', fontsize=9.0, loc='left')

# (f) FUV vs NUV plane
ax = fig.add_subplot(gs2[0, 2])
for nm, m in pred['models'].items():
    p = m['per_R']
    xs = [p[f'ebv0.03_civsupp{s}']['NUV']['total'] for s in ('1.0', '3.5', '20.0')]
    ys = [p[f'ebv0.03_civsupp{s}']['FUV']['total'] for s in ('1.0', '3.5', '20.0')]
    col = OK if nm.startswith('R88') else FGDIM
    ax.plot(xs, ys, color=col, lw=0.8, alpha=0.9); ax.plot(xs[0], ys[0], 'o', color=col, ms=3)
    if nm.startswith('R88'):
        label(ax, xs[0] * 1.03, ys[0], nm.replace('_', ' '), color=OK, size=7)
for band_reg, col in [('south_dec<17.6', ACCENT), ('all', FG), ('north_dec>=17.6', WARN)]:
    f = meas['FUV'][f'beta_OIII_{band_reg}']; n = meas['NUV'][f'beta_OIII_{band_reg}']
    ax.errorbar(n['beta_CU_per_R'], f['beta_CU_per_R'], xerr=n['null_std'], yerr=f['null_std'], fmt='D', color=col, ms=5, label=f'measured {band_reg.split("_")[0]}')
xx = np.array([20, 400])
ax.plot(xx, xx * dust['FUV']['slope_CU_per_MJysr'] / dust['NUV']['slope_CU_per_MJysr'], color=SII_C, ls='--', lw=1, label=f"dust-scattered colour (FUV/NUV {dust['FUV']['slope_CU_per_MJysr']/dust['NUV']['slope_CU_per_MJysr']:.1f})")
ax.plot(xx, xx * dust['diffuse_FUV_over_NUV_median'], color=FGDIM, ls=':', lw=1, label=f"diffuse background colour ({dust['diffuse_FUV_over_NUV_median']:.2f})")
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(20, 400); ax.set_ylim(20, 1500)
ax.set_xlabel('NUV per R [O III] (CU)'); ax.set_ylabel('FUV per R [O III] (CU)')
label(ax, 22, 1100, 'model tracks: C IV ×1 → ÷3.5 → ÷20 (dot = ×1)', color=FGDIM, size=7)
ax.legend(loc='lower right', fontsize=6.8)
ax.set_title('(f) FUV vs NUV per Rayleigh', fontsize=9.5, loc='left')
finish(fig, str(FIG / 'verify_xrayfuv_uv.png'), tight=False)
print('ok')
