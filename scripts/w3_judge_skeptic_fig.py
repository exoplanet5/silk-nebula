"""w3_judge_skeptic summary figure: (a) Monogem-rim registration and its chance probability, (b) southern
single-strand cross profiles vs the consensus shock model seen as a curved limb or a flat sheet, (c) observed vs
predicted ordering/asymmetry statistics, (d) the skeptic judge's posterior over origin hypotheses.

Inputs: data/derived/w3_judge_skeptic_tests.json, w3_judge_skeptic_robust.json, w3_judge_skeptic_modelprof.json,
NSNS / GALEX profiles recomputed with w3_judge_skeptic_tests functions, shock tracks for one representative model.
Output: figures/w3_judge_skeptic_summary.png (dpi 200)
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_skeptic_fig.py
"""
import sys
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import DER, CENTRES, sep_pa
import w3_judge_skeptic_tests as T
import w3_judge_skeptic_modelprof as MP
from silk_style import new_sky_figure, finish, label, OIII_C, WARN, FG, FGDIM, OK, DANGER, ACCENT

OUT = '/Users/mickey/astronomy/silk-nebula/figures/w3_judge_skeptic_summary.png'
POST = [('H1 Monogem limb shock', 0.60, OIII_C), ('H2 other shock, d≲1 kpc', 0.27, '#b388ff'),
        ('H3 paper: SNR at 160 pc', 0.005, DANGER), ('H4 superbubble 1–2 kpc', 0.035, WARN),
        ('H5 photoionised / DIG', 0.01, '#9aa4ae'), ('H6 other non-shock', 0.08, OK)]


def norm(x, p):
    m = np.isfinite(p)
    b = m & (np.abs(x) > 0.6) & (np.abs(x) < 1.1)
    c = np.polyfit(x[b], p[b], 1)
    q = p - np.polyval(c, x)
    return q / np.nanmax(q[np.abs(x) < 0.3])


def main():
    A = json.loads((DER / 'w3_judge_skeptic_tests.json').read_text())
    R = json.loads((DER / 'w3_judge_skeptic_robust.json').read_text())
    fig, ax = new_sky_figure(13, 8.6, 'Silk skeptic judge',
                             'registration chance, limb test, ordering test, posterior', nrows=2, ncols=2)
    # (a)
    a = ax[0, 0]
    ra, de, _ = T.load_spine()
    for key, R0, col, nm in (('K2', 11.00, OIII_C, 'K2 (2-SNR, R 11.0°)'), ('K1', 11.97, WARN, 'K1 (1-SNR, R 11.97°)')):
        d, _ = sep_pa(*CENTRES[key], ra, de)
        a.plot(de, d, color=col, lw=1.6, label=nm)
        a.axhline(R0, color=col, lw=0.8, ls=':')
    a.axhspan(9.9, 12.1, color=OIII_C, alpha=0.10)
    a.axvline(17.5, color=FGDIM, lw=0.6)
    a.set_xlabel('spine Dec (deg)'); a.set_ylabel('distance from Monogem centre (deg)'); a.set_ylim(8.5, 16.8)
    a.legend(loc='upper left', fontsize=7.5, ncol=1)
    nk = A['test_A_registration']['K2']['null']
    label(a, 21.2, 9.6, 'chance of a tangent 5.5° window at R±1.1°\n'
          f"given the trace reaches the rim: {nk['rho_0_26']['p_event_given_touch']:.2f}–{nk['rho_10_16']['p_event_given_touch']:.2f} (K2)\n"
          'K1: observed trace does not register', size=7.2)
    a.set_title('(a) Monogem registration', loc='left', fontsize=9, color=FG)
    # (b) profiles Dec 12-15
    b = ax[0, 1]
    ra_n, de_n, v_n = T.nsns_data()
    e = np.arange(-1.2, 1.2001, 0.025); xc = 0.5 * (e[1:] + e[:-1])
    pn = T.prof_generic(ra_n, de_n, v_n, e, 0.0, T.S15, stat='median')
    import verify_xrayfuv_galex_measure as GM
    D = GM.prep_band('FUV')
    pf = T.prof_generic(D['ra'], D['dec'], D['v'].astype(float), e, 0.0, T.S15, stat='mean', tile=D['tile'])
    b.plot(xc, norm(xc, pn), color=OIII_C, lw=2.0, label='NSNS [O III] (obs)')
    b.plot(xc, norm(xc, pf), color='#b388ff', lw=1.3, label='GALEX FUV (obs)')
    # representative model: median-chi2 consistent model
    summ = json.loads((DER / 'verify_shockmodel_summary.json').read_text())
    TR = np.load(DER / 'verify_shockmodel_tracks.npz')
    m = summ['consistent'][0]
    tag = f"v{m['v']}_x{m['xn']}_b{m['b']}_asplund09_trn"
    N = TR[tag + '__N']; istar = int(np.argmin(np.abs(np.log10(np.maximum(N, 1)) - m['logN'])))
    xg = np.arange(-1.5, 1.5001, 0.005)
    hs, I3, _ = MP.proj(TR[tag + '__s'], TR[tag + '__F_O3'], istar, m['n0'])
    _, I4, _ = MP.proj(TR[tag + '__s'], TR[tag + '__F_CIV'], istar, m['n0'])
    hd = np.degrees(hs / (MP.D_PC * MP.PC)); hpk = hd[np.argmax(I3)]
    c3 = np.interp(hpk - xg, hd, I3, left=0, right=0); c4 = np.interp(hpk - xg, hd, I4, left=0, right=0)
    f0 = T.fit_split(xg, c3, xmax=1.2)
    sk = np.sqrt(max(0.406 ** 2 - f0['fwhm'] ** 2, 1e-6)) / 2.3548 / 0.005
    c3c = gaussian_filter1d(c3, sk, mode='constant'); c4c = gaussian_filter1d(c4, sk, mode='constant')
    b.plot(xg, c3c / c3c.max(), color=OIII_C, lw=1.0, ls='--', label='model [O III], curved limb')
    b.plot(xg, c4c / c4c.max(), color='#b388ff', lw=1.0, ls='--', label='model C IV, curved limb')
    s = TR[tag + '__s'][:istar + 1] * MP.N0REF / m['n0']; Fc = TR[tag + '__F_O3'][:istar + 1]
    hf = np.degrees(0.5 * (s[1:] + s[:-1]) / (MP.D_PC * MP.PC)); jf = np.diff(Fc) / np.maximum(np.diff(s), 1e-30)
    fl = np.interp(hf[np.argmax(jf)] - xg, hf, jf, left=0, right=0)
    f1 = T.fit_split(xg, fl, xmax=1.2)
    sk1 = np.sqrt(max(0.406 ** 2 - f1['fwhm'] ** 2, 1e-6)) / 2.3548 / 0.005
    flc = gaussian_filter1d(fl, sk1, mode='constant')
    b.plot(xg, flc / flc.max(), color=OK, lw=1.0, ls=':', label='model [O III], flat sheet')
    b.axvline(0, color=FGDIM, lw=0.6)
    b.set_xlim(-1.0, 1.0); b.set_ylim(-0.25, 1.15)
    b.set_xlabel('offset from E spine (deg; + = east, outside Monogem)'); b.set_ylabel('normalised')
    b.legend(loc='upper left', fontsize=6.8)
    b.set_title(f'(b) Dec 12–15 single strand; model {tag.split("_asplund")[0]}', loc='left', fontsize=9, color=FG)
    # (c) scoreboard
    c = ax[1, 0]
    obs = R['observed']['dec12_15']
    asy = [q['asym'] for q in obs['nsns']]; asyj = np.median([q['asym_jk'] for q in obs['nsns']])
    fuv = [q['shift'] for q in obs['galex'] if q['band'] == 'FUV']; fuvj = np.median([q['shift_jk'] for q in obs['galex'] if q['band'] == 'FUV'])
    o2 = [q['shift'] for q in obs['desi'] if np.isfinite(q['shift'])]; o2j = np.nanmedian([q['shift_jk'] for q in obs['desi'] if np.isfinite(q['shift_jk']) and q['shift_jk'] > 0])
    msum = R['model']['summary_min_median_max']
    rows = [('[O III] σW/σE', np.median(asy), asyj, (min(asy), max(asy)),
             {'curved limb': msum['curved_dec12_15']['asym'], 'flat sheet': msum['flat_dec12_15']['asym']}, (0, 4)),
            ('ΔFUV−[O III] (°)', np.median(fuv), fuvj, (min(fuv), max(fuv)),
             {'shock (C IV)': msum['curved_dec12_15']['shift_CIV'],
              'shock (50 % C IV)': [0.035, 0.045, 0.065], 'hot-side-west interface': [-0.12, -0.06, -0.02]}, (-0.15, 0.2)),
            ('Δ[O II]−[O III] (°), DESI', np.median(o2), o2j, (min(o2), max(o2)),
             {'shock': msum['curved_dec12_15']['shift_O2'], 'hot-side-west interface': [0.02, 0.05, 0.1]}, (-0.15, 0.2))]
    c.axis('off')
    for k, (name, val, err, rng, preds, xl) in enumerate(rows):
        sub = c.inset_axes([0.36, 0.72 - k * 0.33, 0.58, 0.22])
        sub.set_xlim(*xl); sub.set_yticks([])
        y = 0
        for j, (pn_, (lo, md, hi)) in enumerate(preds.items()):
            yy = 1 + j
            sub.plot([lo, hi], [yy, yy], color=WARN if 'interface' not in pn_ else DANGER, lw=5, alpha=0.6, solid_capstyle='butt')
            sub.plot([md], [yy], '|', color=FG, ms=8)
            sub.text(xl[1], yy, ' ' + pn_, color=FGDIM, fontsize=6.5, va='center')
        sub.errorbar([val], [0], xerr=[[err], [err]], fmt='o', color=OIII_C, ms=5, capsize=3)
        sub.plot(rng, [0, 0], color=OIII_C, lw=1, alpha=0.6)
        sub.set_ylim(-0.7, len(preds) + 0.7)
        sub.axvline(0 if 'Δ' in name else 1, color=FGDIM, lw=0.5)
        sub.tick_params(labelsize=7)
        c.text(0.0, 0.83 - k * 0.33, name, color=FG, fontsize=7.8, transform=c.transAxes)
        c.text(0.0, 0.77 - k * 0.33, f'obs {val:+.3f} ± {err:.3f}' if 'Δ' in name else f'obs {val:.2f} ± {err:.2f}',
               color=OIII_C, fontsize=7.5, transform=c.transAxes)
    c.set_title('(c) observed (blue) vs predicted (bands), Dec 12–15', loc='left', fontsize=9, color=FG)
    # (d) posterior
    d = ax[1, 1]
    names = [p[0] for p in POST][::-1]; vals = [p[1] for p in POST][::-1]; cols = [p[2] for p in POST][::-1]
    d.barh(names, vals, color=cols, alpha=0.85)
    for i, v in enumerate(vals):
        d.text(v + 0.01, i, f'{v:.3f}' if v < 0.05 else f'{v:.2f}', color=FG, fontsize=8, va='center')
    d.set_xlim(0, 0.75); d.set_xlabel('posterior probability (skeptic judge)')
    d.tick_params(axis='y', labelsize=8)
    d.set_title('(d) origin hypotheses', loc='left', fontsize=9, color=FG)
    finish(fig, OUT)
    print('saved', OUT)


if __name__ == '__main__':
    main()
