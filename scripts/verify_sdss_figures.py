"""verify_sdss figures (dpi 200, silk_style).
  verify_sdss_skymodel_bias.png  per-group [O III] ON-OFF from total sky vs from pipeline residual (flux) and sky model
  verify_sdss_abscal.png         absolute-calibration checks (imaging sky, 5577, Halpha) by instrument
  verify_sdss_stacks.png         ON-OFF stacked spectra, south / north / Halpha patch
  verify_sdss_velocity.png       zero points from sky lines and [O III] velocities
  verify_sdss_precovery.png      per-plate [O III] detections vs date and on the sky
Usage: python -I verify_sdss_figures.py <diff_tag> <f_sdss> <f_boss>
"""
import sys, csv, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import silk_style as ss
import matplotlib.pyplot as plt
from verify_sdss_lib import ROOT, VAC, load_windows, flam_to_R_per_A, load_strands
from verify_sdss_results import load, good_groups

FIG = ROOT / 'figures'


def fcsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def col(rows, k):
    o = []
    for r in rows:
        try:
            o.append(float(r[k]))
        except (ValueError, TypeError, KeyError):
            o.append(np.nan)
    return np.array(o)


def stacks(D, fib, grp, gok, on, off, fac, w_by_g):
    """Weighted mean over groups of (mean ON - mean OFF) total sky, R/A (diffuse-corrected)."""
    lam = D['lam']
    num = np.zeros(lam.size); den = np.zeros(lam.size)
    for g in np.where(gok)[0]:
        if not np.isfinite(w_by_g[g]) or w_by_g[g] <= 0:
            continue
        s = np.where(fib['gid'] == g)[0]
        so = s[on[s]]; sf = s[off[s]]
        if len(so) < 1 or len(sf) < 3:
            continue
        diam = fib['diam'][s[0]]
        X = flam_to_R_per_A(D['tot'][s].astype(float), lam, diam)
        X[D['ivar'][s] <= 0] = np.nan
        Xf = X[np.isin(s, sf)]
        med = np.nanmedian(Xf, 0)
        mad = 1.4826 * np.nanmedian(np.abs(Xf - med), 0)
        Xf = np.where(np.abs(Xf - med) < 4 * np.maximum(mad, 1e-6), Xf, np.nan)
        d = (np.nanmean(X[np.isin(s, so)], 0) - np.nanmean(Xf, 0)) * fac[g]
        ok = np.isfinite(d)
        num[ok] += w_by_g[g] * d[ok]; den[ok] += w_by_g[g]
    with np.errstate(all='ignore'):
        return num / den


def main(tag, f_sdss, f_boss):
    D = load_windows(); lam = D['lam']
    fib, grp, amp, Z = load(tag)
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    R = json.load(open(ROOT / f'data/derived/verify_sdss_results_{tag}_r0.15.json'))
    dec = fib['dec']
    ON = fib['dmin'] < 0.15
    OFF = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    # group weights from [O III] OFF scatter
    a = amp['OIII:OIII5007']
    ng = len(grp['gid'])
    wg = np.full(ng, np.nan)
    for g in np.where(gok)[0]:
        af = a[(fib['gid'] == g) & OFF & np.isfinite(a)]
        if len(af) >= 3:
            sd = 1.4826 * np.median(np.abs(af - np.median(af))) * fac[g]
            wg[g] = 1 / max(sd, 1e-3) ** 2
    # ------------------------------------------------ 1. sky-model bias
    m = (lam > 4985) & (lam < 5032)
    rows = []
    for g in np.where(gok)[0]:
        s = np.where(fib['gid'] == g)[0]
        so = s[(ON & (dec < 20))[s]]; sf = s[OFF[s]]
        if len(so) < 1 or len(sf) < 3:
            continue
        diam = fib['diam'][s[0]]
        cv = lambda X: flam_to_R_per_A(X[:, m].astype(float), lam[m], diam)
        tot = cv(D['tot'][so]).mean(0) - cv(D['tot'][sf]).mean(0)
        sky = cv(D['sky'][so]).mean(0) - cv(D['sky'][sf]).mean(0)
        x = lam[m]; mu = VAC['OIII5007'] * (1 + 13 / 2.998e5)
        G = np.exp(-0.5 * ((x - mu) / 1.1) ** 2) / (np.sqrt(2 * np.pi) * 1.1)
        A = np.c_[np.ones_like(x), x - x.mean(), G]
        at = np.linalg.lstsq(A, tot, rcond=None)[0][2]; ak = np.linalg.lstsq(A, sky, rcond=None)[0][2]
        rows.append((grp['instrument'][g], at * fac[g], ak * fac[g], (at - ak) * fac[g], len(so)))
    rows = np.array(rows, dtype=object)
    fig, ax = ss.new_sky_figure(10, 5.2, 'SDSS sky fibres · pipeline sky model absorbs [O III] structure',
                                'south strand ON−OFF per plate–spectrograph group, diffuse-corrected R', ncols=2)
    for ins, c in (('SDSS', ss.WARN), ('BOSS', ss.ACCENT)):
        s = rows[:, 0] == ins
        ax[0].scatter(rows[s, 1].astype(float), rows[s, 3].astype(float), s=12 + 3 * rows[s, 4].astype(float), color=c, alpha=0.8, label=f'{ins} residual (flux)')
        ax[0].scatter(rows[s, 1].astype(float), rows[s, 2].astype(float), s=12 + 3 * rows[s, 4].astype(float), facecolor='none', edgecolor=c, label=f'{ins} sky model')
    lim = [0, max(3.0, float(np.nanmax(rows[:, 1].astype(float))) * 1.05)]
    ax[0].plot(lim, lim, color=ss.FGDIM, lw=0.8); ax[0].plot(lim, [0.5 * l for l in lim], color=ss.FGDIM, lw=0.6, ls=':')
    ax[0].set_xlabel('ON−OFF from total sky (flux+sky) [R]'); ax[0].set_ylabel('ON−OFF from one component [R]')
    ax[0].legend(loc='upper left'); ax[0].set_xlim(lim); ax[0].set_ylim(-0.3, lim[1])
    frac = rows[:, 3].astype(float) / rows[:, 1].astype(float)
    for ins, c in (('SDSS', ss.WARN), ('BOSS', ss.ACCENT)):
        s = rows[:, 0] == ins
        ax[1].hist(np.clip(frac[s], -0.5, 1.5), bins=np.linspace(-0.5, 1.5, 21), color=c, alpha=0.6, label=f'{ins} median {np.median(frac[s]):.2f}')
    ax[1].axvline(1, color=ss.FGDIM, lw=0.8)
    ax[1].set_xlabel('fraction of ON−OFF [O III] left in the residual'); ax[1].set_ylabel('groups'); ax[1].legend()
    ss.finish(fig, FIG / 'verify_sdss_skymodel_bias.png')
    sm = {'n_groups': int(len(rows)), 'median_frac_residual_SDSS': float(np.median(frac[rows[:, 0] == 'SDSS'])),
          'median_frac_residual_BOSS': float(np.median(frac[rows[:, 0] == 'BOSS'])),
          'mean_total_R': float(np.mean(rows[:, 1].astype(float))), 'mean_resid_R': float(np.mean(rows[:, 3].astype(float))),
          'mean_skymodel_R': float(np.mean(rows[:, 2].astype(float)))}
    json.dump(sm, open(ROOT / 'data/derived/verify_sdss_skymodel_bias.json', 'w'), indent=1)
    print('skymodel bias', sm)
    # ------------------------------------------------ 2. abscal
    A = np.load(ROOT / 'data/derived/verify_sdss_abscal_groups.npz')
    AJ = json.load(open(ROOT / 'data/derived/verify_sdss_abscal.json'))
    fig, ax = ss.new_sky_figure(13, 5.4, 'SDSS sky fibres · absolute diffuse-light calibration',
                                'dark-sky plate–spectrograph groups; as calibrated (point-source convention)', ncols=3)
    for ins, c in (('SDSS', ss.WARN), ('BOSS', ss.ACCENT)):
        s = A['dark'] & (A['inst'] == ins) & np.isfinite(A['img_g'])
        d = A['sb_g'][s] - A['img_g'][s]
        ax[0].hist(d, bins=np.linspace(-1.5, 1.0, 26), color=c, alpha=0.6, label=f'{ins}: median {np.median(d):+.2f} → f={10 ** (0.4 * np.median(d)):.2f}')
    ax[0].axvline(0, color=ss.FGDIM, lw=0.8)
    ax[0].set_xlabel('g sky: spectro (as calibrated) − imaging [mag/″²]'); ax[0].set_ylabel('groups'); ax[0].legend(fontsize=7.5)
    for ins, c, f in (('SDSS', ss.WARN, f_sdss), ('BOSS', ss.ACCENT, f_boss)):
        s = A['dark'] & (A['inst'] == ins)
        ax[1].hist(A['I5577'][s], bins=np.linspace(0, 900, 31), histtype='step', color=c, lw=1.0, ls=':', label=f'{ins} as calibrated')
        ax[1].hist(A['I5577'][s] * f, bins=np.linspace(0, 900, 31), histtype='step', color=c, lw=1.6, label=f'{ins} ×{f:.2f}: median {np.median(A["I5577"][s] * f):.0f} R')
        ax[2].hist(A['IHa'][s] * f, bins=np.linspace(0, 20, 26), histtype='step', color=c, lw=1.6, label=f'{ins} ×{f:.2f}: median {np.nanmedian(A["IHa"][s] * f):.1f} R')
    ax[1].set_xlabel('[O I] 5577 [R]'); ax[1].legend(fontsize=7.5)
    ax[2].axvspan(2.9, 4.4, color=ss.HA_C, alpha=0.15, label='WHAM Galactic Hα 2.9–4.4 R')
    ax[2].set_xlabel('total-sky Hα (Galactic + geocoronal) [R]'); ax[2].legend(fontsize=7.5)
    ss.finish(fig, FIG / 'verify_sdss_abscal.png')
    # ------------------------------------------------ 3. stacks
    regs = [('south (Dec<20)', ON & (dec < 20), OFF), ('north (Dec≥20)', ON & (dec >= 20), OFF),
            ('Hα patch', fib['inpatch30'], (fib['dmin'] > 0.45) & (fib['dmin'] < 2.5) & ~fib['inpatch30'])]
    wins = [('[O II] 3726/29 (BOSS)', 3705, 3752, ['OII3726', 'OII3729']), ('Hβ', 4845, 4880, ['Hb']),
            ('[O III] 4959/5007', 4945, 5022, ['OIII4959', 'OIII5007']), ('Hγ · Hg · [O III]4363', 4330, 4380, ['Hgam', 'OIII4363', 'Hg4358']),
            ('Hα · [N II]', 6540, 6600, ['NII6548', 'Ha', 'NII6583']), ('[S II]', 6705, 6745, ['SII6716', 'SII6731'])]
    fig, ax = ss.new_sky_figure(14, 8.4, 'SDSS sky fibres · ON−OFF strand spectra',
                                f'total sky (flux+sky), same-group OFF 0.45–2°, diffuse-corrected (×{f_sdss}/×{f_boss}); R per Å', nrows=3, ncols=6)
    for i, (nm, on, off) in enumerate(regs):
        sp = stacks(D, fib, grp, gok, on, off & ~on, fac, wg)
        for j, (wn, lo, hi, lines) in enumerate(wins):
            x = ax[i, j]
            mm = (lam > lo) & (lam < hi)
            x.step(lam[mm], sp[mm], where='mid', color=ss.FG, lw=0.8)
            x.axhline(0, color=ss.FGDIM, lw=0.5)
            for L in lines:
                x.axvline(VAC[L] * (1 + 13 / 2.998e5), color=ss.OIII_C if 'OIII' in L else (ss.WARN if 'Hg' in L else ss.HA_C), lw=0.6, ls=':')
            if i == 0:
                x.set_title(wn, fontsize=8.5)
            if j == 0:
                x.set_ylabel(nm + '\nR/Å')
            x.tick_params(labelsize=7)
    ss.finish(fig, FIG / 'verify_sdss_stacks.png')
    # ------------------------------------------------ 4. velocity / zero points
    Zp = fcsv(ROOT / 'data/derived/verify_sdss_zeropoints.csv')
    inst = np.array([r['instrument'] for r in Zp]); hr = col(Zp, 'helio_rv')
    fig, ax = ss.new_sky_figure(12, 5.6, 'SDSS sky fibres · wavelength zero point and [O III] velocity',
                                'sky-line offsets from the topocentric position (v_exp = −HELIO_RV); heliocentric [O III]', ncols=2)
    Gd = fcsv(ROOT / 'data/derived/verify_sdss_groups.csv')
    sbg = col(Gd, 'sb5070_AB_cal'); nfg = col(Gd, 'n_fib'); ig = np.array([r['instrument'] for r in Gd])
    medg = {i_: np.nanmedian(sbg[(ig == i_) & (nfg >= 6)]) for i_ in ('SDSS', 'BOSS')}
    darkg = (nfg >= 6) & (sbg > np.array([medg[i_] for i_ in ig]) - 0.75)
    for k, c, lab in (('dv_OI5577', ss.OK, '[O I]5577'), ('dv_Hg4358', ss.WARN, 'Hg I 4358'), ('dv_NI', ss.NII_C, '[N I]5200 (ratio 1.6)'), ('dv_OI6300', ss.HA_C, '[O I]6300')):
        y = col(Zp, k); s = np.isfinite(y) & (np.abs(y) < 40) & darkg
        ax[0].scatter(hr[s], y[s], s=6, color=c, alpha=0.6, label=f'{lab}: median {np.median(y[s]):+.1f}')
    ax[0].axhline(0, color=ss.FGDIM, lw=0.6)
    ax[0].set_xlabel('HELIO_RV [km/s]'); ax[0].set_ylabel('measured − expected [km/s]'); ax[0].set_ylim(-30, 30); ax[0].legend(fontsize=7.5)
    items = []
    for nm, lab in (('south_dec_lt20', 'SDSS south'), ('W_strand_15_20', 'SDSS W strand'), ('E_strand_15_20', 'SDSS E strand'),
                    ('north_dec_ge20', 'SDSS north'), ('all', 'SDSS all')):
        v = R['regions'][nm].get('velocity_OIII', {})
        if 'v_corr' in v:
            items.append((lab, v['v_corr'], v['ev_boot'], ss.ACCENT))
    items += [('wave-1 SDSS south', 12.5, 2.6, ss.FGDIM), ('DESI all (wave-1)', 14.5, 1.5, ss.OK), ('DESI E strand', 13.0, 0.5, ss.OK),
              ('DESI W strand', 18.5, 1.25, ss.OK), ('Wu+26 LAMOST P1', 26.8, 4.3, ss.DANGER), ('Wu+26 P2', 26.8, 4.1, ss.DANGER),
              ('Wu+26 P3', 23.8, 3.3, ss.DANGER), ('WHAM Hα (V_helio)', 12.3, 0.5, ss.HA_C)]
    for i, (lab, v, e, c) in enumerate(items):
        ax[1].errorbar(v, -i, xerr=e, fmt='o', color=c, ms=4, capsize=2)
        ss.label(ax[1], 33, -i, lab, size=7.5, va='center', color=c)
    ax[1].set_xlim(-2, 52); ax[1].set_yticks([]); ax[1].set_xlabel('V_helio [O III] (km/s)')
    ax[1].axvline(12.2, color=ss.FGDIM, lw=0.6, ls=':')
    ss.finish(fig, FIG / 'verify_sdss_velocity.png')
    # ------------------------------------------------ 5. precovery
    P = fcsv(ROOT / 'data/derived/verify_sdss_precovery.csv')
    PP = {(r['plate'], r['mjd']): r for r in fcsv(ROOT / 'data/derived/verify_sdss_precovery_placebo.csv')}
    mj = col(P, 'mjdfirst'); snr = col(P, 'onoff_snr'); psnr = col(P, 'prof_snr'); ra = col(P, 'ra'); de = col(P, 'dec')
    ins = np.array([r['instrument'] for r in P])
    best = np.array([float(PP[(r['plate'], r['mjd'])]['snr_placebo']) if (r['plate'], r['mjd']) in PP and PP[(r['plate'], r['mjd'])]['snr_placebo'] not in ('nan', '') else np.nan for r in P])
    from astropy.time import Time
    yr = Time(mj, format='mjd').decimalyear
    fig, ax = ss.new_sky_figure(12, 5.0, 'SDSS sky fibres · spectroscopic precovery of the Silk [O III]',
                                'per plate–MJD ON−OFF [O III] significance including the per-plate placebo scatter (trace shifted ±0.6–1.2°)', ncols=2)
    for i_, c in (('SDSS', ss.WARN), ('BOSS', ss.ACCENT)):
        s = ins == i_
        ax[0].scatter(yr[s], best[s], s=14, color=c, label=i_)
    ax[0].axhline(3, color=ss.DANGER, lw=0.8, ls='--'); ax[0].axhline(0, color=ss.FGDIM, lw=0.5)
    ax[0].set_xlabel('year of first exposure'); ax[0].set_ylabel('S/N (stat ⊕ placebo)'); ax[0].legend()
    for r_, b_, y_ in zip(P, best, yr):
        if r_['plate'] in ('1203', '857', '2054', '2074'):
            ss.label(ax[0], y_ + 0.15, (b_ if np.isfinite(b_) else 0) + 0.4, r_['plate'], size=7)
    S, patch, tr = load_strands()
    for k in ('E', 'W', 'W2'):
        ax[1].plot(S[k][:, 0], S[k][:, 1], color=ss.OIII_C, lw=0.8)
    det = np.nan_to_num(best, nan=-9) >= 3
    ax[1].scatter(ra[~det], de[~det], s=10, color=ss.FGDIM, label='plate, < 3σ or untested')
    sc = ax[1].scatter(ra[det], de[det], s=24, c=yr[det], cmap='plasma', label='plate, ≥ 3σ')
    plt.colorbar(sc, ax=ax[1], label='year')
    ax[1].invert_xaxis(); ax[1].set_xlabel('RA [deg]'); ax[1].set_ylabel('Dec [deg]'); ax[1].legend(fontsize=7.5, loc='lower left')
    ax[1].set_xlim(122, 109); ax[1].set_ylim(8, 30)
    ss.finish(fig, FIG / 'verify_sdss_precovery.png')
    print('figures done')


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), float(a[2]))
