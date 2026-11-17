"""verify_sdss stage 3: region ON-OFF line intensities, ratios, [O III] velocity with airglow zero point, precovery.

Inputs: verify_sdss_fibres.csv, verify_sdss_groups.csv, verify_sdss_diff_<tag>.npz, verify_sdss_windows.npz.
Usage: python -I verify_sdss_results.py <diff_tag> <r_on> [f_sdss f_boss]
  r_on = ON radius (deg) around the nearest strand.  f_* = diffuse-light factors applied to the final intensities
  (defaults 1 = as calibrated).  Errors: bootstrap over plate-spectrograph groups (2000 resamples).
Writes data/derived/verify_sdss_results_<diff_tag>_r<r_on>.json
"""
import sys, csv, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, VAC, C, load_windows, flam_to_R_per_A, fit_lines

RNG = np.random.default_rng(20261116)
NB = 2000


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def fcol(rows, k, typ=float):
    out = []
    for r in rows:
        try:
            out.append(typ(r[k]))
        except (ValueError, TypeError):
            out.append(np.nan if typ is float else -1)
    return np.array(out)


def load(tag):
    F = read_csv(ROOT / 'data/derived/verify_sdss_fibres.csv')
    G = read_csv(ROOT / 'data/derived/verify_sdss_groups.csv')
    Z = np.load(ROOT / f'data/derived/verify_sdss_diff_{tag}.npz')
    fib = {k: fcol(F, k) for k in ('ra', 'dec', 'dE', 'dW', 'dW2', 'dmin', 'diam', 'mjd', 'mjdfirst', 'plate')}
    fib['gid'] = fcol(F, 'gid', int)
    fib['nearest'] = np.array([r['nearest'] for r in F])
    fib['inpatch'] = fcol(F, 'in_patch', int).astype(bool)
    fib['inpatch30'] = fcol(F, 'in_patch30', int).astype(bool)
    fib['instrument'] = np.array([r['instrument'] for r in F])
    fib['release'] = np.array([r['release'] for r in F])
    grp = {k: fcol(G, k) for k in G[0].keys() if k not in ('release', 'instrument')}
    grp['instrument'] = np.array([r['instrument'] for r in G]); grp['release'] = np.array([r['release'] for r in G])
    amp = {k.replace('__', ':'): Z[k] for k in Z.files if '__' in k}
    return fib, grp, amp, Z


def good_groups(grp):
    """Reject bright-sky groups (moon/twilight) and tiny groups."""
    sb = grp['sb5070_AB_cal']
    ok = (grp['n_fib'] >= 6) & np.isfinite(sb)
    med = {ins: np.nanmedian(sb[ok & (grp['instrument'] == ins)]) for ins in ('SDSS', 'BOSS')}
    ref = np.array([med[i] for i in grp['instrument']])
    ok &= sb > ref - 0.75           # no more than 0.75 mag brighter than the instrument median
    return ok, med


def group_stats(fib, amp, key, on, off, gok):
    """Per group: mean ON amplitude, error from OFF scatter. Returns arrays over groups."""
    a = amp[key]
    ng = int(fib['gid'].max()) + 1
    m_on = np.full(ng, np.nan); e = np.full(ng, np.nan); n_on = np.zeros(ng, int); n_off = np.zeros(ng, int)
    for g in np.where(gok)[0]:
        s = fib['gid'] == g
        ao = a[s & on & np.isfinite(a)]; af = a[s & off & np.isfinite(a)]
        if len(ao) < 1 or len(af) < 3:
            continue
        sd = 1.4826 * np.median(np.abs(af - np.median(af)))
        sd = max(sd, np.std(af) * 0.5)
        m_on[g] = np.mean(ao) - np.mean(af)
        e[g] = sd * np.sqrt(1 / len(ao) + 1 / len(af))
        n_on[g] = len(ao); n_off[g] = len(af)
    return m_on, e, n_on, n_off


def combine(vals, w, idx=None):
    if idx is None:
        idx = np.arange(len(vals))
    v = vals[idx]; ww = w[idx]
    m = np.isfinite(v) & np.isfinite(ww) & (ww > 0)
    return np.sum(v[m] * ww[m]) / np.sum(ww[m]) if m.any() else np.nan


def region_analysis(fib, grp, amp, gok, on, off, fac, lines):
    """Combine groups with common weights from the [O III] error; bootstrap ratios."""
    st = {L: group_stats(fib, amp, L, on, off, gok) for L in lines}
    mo, eo, non, noff = st['OIII:OIII5007']
    use = np.isfinite(mo) & np.isfinite(eo) & (eo > 0)
    gids = np.where(use)[0]
    if len(gids) < 2:
        return None
    w = np.zeros_like(mo); w[use] = 1 / eo[use] ** 2
    out = {'n_groups': int(len(gids)), 'n_on_fibres': int(non[use].sum()), 'n_off_fibres': int(noff[use].sum())}
    boots = {}
    for L in lines:
        m, e, _, _ = st[L]
        # per-line, common weights; groups lacking this line (e.g. [O II] in SDSS-I/II) drop out
        val = combine(m * fac, w, gids)
        bb = np.array([combine(m * fac, w, RNG.choice(gids, len(gids))) for _ in range(NB)])
        # formal error with common weights
        mm = np.isfinite(m[gids])
        ef = np.sqrt(np.sum((w[gids][mm] * (e[gids][mm] * fac[gids][mm])) ** 2)) / np.sum(w[gids][mm]) if mm.any() else np.nan
        out[L] = {'I_R': float(val), 'err_boot': float(np.nanstd(bb)), 'err_formal': float(ef),
                  'n_groups': int(mm.sum())}
        boots[L] = bb
    return out, boots, st, w, gids


def ratio_summary(bnum, bden, val_num, val_den):
    r = bnum / bden
    q = np.nanpercentile(r, [2.5, 16, 50, 84, 97.5])
    return {'value': float(val_num / val_den), 'p2.5': float(q[0]), 'p16': float(q[1]), 'p50': float(q[2]),
            'p84': float(q[3]), 'p97.5': float(q[4])}


def lower_limit_ratio(bnum, bden):
    """For a denominator consistent with zero: 2-sigma (97.7%) lower limit of num/den from num/(den_upper)."""
    den_hi = np.nanpercentile(bden, 97.7)
    return float(np.nanmedian(bnum) / den_hi) if den_hi > 0 else np.inf


# -------------------------------------------------------------- velocity stack of [O III]
def oiii_velocity(fib, grp, gok, on, off, w_g, gids, D, zp, nboot=400):
    lam = D['lam']
    m = (lam > 4935) & (lam < 5032)
    lm = lam[m]
    specs = {}; sigs = {}
    for g in gids:
        s = np.where(fib['gid'] == g)[0]
        so = s[on[s]]; sf = s[off[s]]
        if len(so) < 1 or len(sf) < 3:
            continue
        diam = fib['diam'][s[0]]
        X = flam_to_R_per_A(D['tot'][:, m][s].astype(float), lm, diam)
        iv = D['ivar'][:, m][s]
        X[iv <= 0] = np.nan
        med = np.nanmedian(X[np.isin(s, sf)], 0)
        Xo = X[np.isin(s, so)]; Xf = X[np.isin(s, sf)]
        with np.errstate(all='ignore'):
            mad = 1.4826 * np.nanmedian(np.abs(Xf - med), 0)
            Xf = np.where(np.abs(Xf - med) < 4 * np.maximum(mad, 1e-6), Xf, np.nan)
        specs[g] = np.nanmean(Xo, 0) - np.nanmean(Xf, 0)
        sigs[g] = float(np.nanmedian(D['wdisp'][:, m][s])) * 5008.0 * np.log(10) * 1e-4
    keys = list(specs)

    def stack(ks):
        W = np.array([w_g[k] for k in ks])
        S = np.array([specs[k] for k in ks])
        good = np.isfinite(S)
        num = np.nansum(np.where(good, S * W[:, None], 0), 0); den = np.sum(np.where(good, W[:, None], 0), 0)
        sp = num / den
        sigma = np.sum(W * np.array([sigs[k] for k in ks])) / np.sum(W)
        z = np.sum(W * np.array([zp[k] for k in ks])) / np.sum(W)
        return sp, sigma, z

    def fitv(sp, sigma):
        ok = np.isfinite(sp)
        # error spectrum from line-free pixels
        rms = np.nanstd(sp[ok & ((lm < 4950) | ((lm > 4970) & (lm < 4998)) | (lm > 5018))])
        wgt = np.full(ok.sum(), 1 / rms ** 2)
        vs = np.arange(-40, 70.01, 0.5)
        ch = np.array([fit_lines(lm[ok], sp[ok], wgt, ['OIII4959', 'OIII5007'], v, sigma, npoly=2,
                                 ties={'OIII4959': ('OIII5007', 1 / 3.01)})['_chi2'] for v in vs])
        j = np.argmin(ch)
        jj = slice(max(j - 6, 0), min(j + 7, len(vs)))
        pc = np.polyfit(vs[jj], ch[jj], 2)
        vb = -pc[1] / (2 * pc[0]) if pc[0] > 0 else vs[j]
        f = fit_lines(lm[ok], sp[ok], wgt, ['OIII4959', 'OIII5007'], vb, sigma, npoly=2, ties={'OIII4959': ('OIII5007', 1 / 3.01)})
        # free-ratio check
        f2 = fit_lines(lm[ok], sp[ok], wgt, ['OIII4959', 'OIII5007'], vb, sigma, npoly=2)
        return vb, float(np.sqrt(2 / pc[0])) if pc[0] > 0 else np.nan, f['OIII5007'], f2

    sp, sigma, z = stack(keys)
    vb, ev, a5007, f2 = fitv(sp, sigma)
    bv = []
    for _ in range(nboot):
        ks = list(RNG.choice(keys, len(keys)))
        s_, sg_, z_ = stack(ks)
        try:
            v_, _, _, _ = fitv(s_, sg_)
            bv.append(v_ - z_)
        except Exception:
            pass
    bv = np.array(bv)
    return {'v_raw': float(vb), 'ev_formal': ev, 'zp_mean': float(z), 'v_corr': float(vb - z),
            'ev_boot': float(np.nanstd(bv)), 'v_boot_p16_p84': [float(np.nanpercentile(bv, 16)), float(np.nanpercentile(bv, 84))],
            'I5007_stack': a5007[0], 'ratio5007_4959_free': float(f2['OIII5007'][0] / f2['OIII4959'][0]) if f2['OIII4959'][0] != 0 else np.nan,
            'n_groups': len(keys), 'lsf_sigma_A': sigma}, (lm, sp)


def main(tag, r_on, f_sdss=1.0, f_boss=1.0):
    fib, grp, amp, Z = load(tag)
    D = load_windows()
    gok, med = good_groups(grp)
    off_lo, off_hi = float(Z['off_lo']), float(Z['off_hi'])
    excl = int(Z['excl_patch'])
    off = (fib['dmin'] > off_lo) & (fib['dmin'] < off_hi)
    if excl:
        off &= ~fib['inpatch30']
    fac = np.where(grp['instrument'] == 'SDSS', f_sdss, f_boss)
    lines = ['OIII:OIII5007', 'Hb:Hb', 'HaNII:Ha', 'HaNII:NII6583', 'SII:SII6716', 'SII:SII6731', 'OII:OII3726',
             'OII:OII3729', 'B4363:OIII4363', 'B4363:Hg4358', 'B4363:Hgam', 'HeII:HeII4686']
    ON = fib['dmin'] < r_on
    dec = fib['dec']
    regions = {
        'all': ON,
        'south_dec_lt20': ON & (dec < 20),
        'north_dec_ge20': ON & (dec >= 20),
        'W_strand_15_20': ON & (fib['nearest'] == 'W') & (dec >= 15) & (dec < 20),
        'E_strand_15_20': ON & (fib['nearest'] == 'E') & (dec >= 15) & (dec < 20),
        'south_dec_lt17p6': ON & (dec < 17.6),
        'north_dec_ge20_excl_patch30': ON & (dec >= 20) & ~fib['inpatch30'],
        'north_dec_20_24': ON & (dec >= 20) & (dec < 24),
        'north_dec_ge24': ON & (dec >= 24),
        'halpha_patch_30pct': fib['inpatch30'],
        'halpha_patch_halfmax': fib['inpatch'],
    }
    # velocity zero point per group from [O I]5577: measured minus expected topocentric (expected = -HELIO_RV)
    zp = grp['v_OI5577'] + grp['helio_rv']
    zp6300 = grp['v_OI6300'] + grp['helio_rv']
    out = {'tag': tag, 'r_on': r_on, 'off': [off_lo, off_hi], 'exclude_patch_from_off': bool(excl),
           'pairing': str(Z['pairing']), 'f_sdss': f_sdss, 'f_boss': f_boss,
           'n_groups_good': int(gok.sum()), 'sb5070_median_cal': med, 'regions': {}}
    # helio sign test
    ok5577 = gok & np.isfinite(grp['v_OI5577'])
    A = np.c_[np.ones(ok5577.sum()), grp['helio_rv'][ok5577]]
    cf = np.linalg.lstsq(A, grp['v_OI5577'][ok5577], rcond=None)[0]
    tai = np.isfinite(grp['vcorr_astropy'])
    out['helio_sign_test'] = {'slope_v5577_vs_HELIO_RV': float(cf[1]), 'intercept': float(cf[0]),
                              'median_HELIO_RV_plus_vcorr_astropy': float(np.nanmedian(grp['helio_rv'][tai] + grp['vcorr_astropy'][tai])),
                              'max_abs_HELIO_RV_plus_vcorr': float(np.nanmax(np.abs(grp['helio_rv'][tai] + grp['vcorr_astropy'][tai]))),
                              'n_tai': int(tai.sum()),
                              'zp5577_median': float(np.nanmedian(zp[gok])), 'zp5577_p16_p84': [float(np.nanpercentile(zp[gok], 16)), float(np.nanpercentile(zp[gok], 84))],
                              'zp6300_median': float(np.nanmedian(zp6300[gok]))}
    for nm, sel in regions.items():
        r = region_analysis(fib, grp, amp, gok, sel, off & ~sel, fac, lines)
        if r is None:
            out['regions'][nm] = None
            continue
        res, boots, st, w, gids = r
        # ratios (photon units)
        rat = {}
        def RS(a, b):
            return ratio_summary(boots[a], boots[b], res[a]['I_R'], res[b]['I_R'])
        rat['OIII5007/Hb'] = RS('OIII:OIII5007', 'Hb:Hb')
        rat['OIII5007/Hb_lower2sig'] = lower_limit_ratio(boots['OIII:OIII5007'], boots['Hb:Hb'])
        rat['OIII5007/Ha'] = RS('OIII:OIII5007', 'HaNII:Ha')
        rat['NII6583/Ha'] = RS('HaNII:NII6583', 'HaNII:Ha')
        bS = boots['SII:SII6716'] + boots['SII:SII6731']
        rat['SII/Ha'] = ratio_summary(bS, boots['HaNII:Ha'], res['SII:SII6716']['I_R'] + res['SII:SII6731']['I_R'], res['HaNII:Ha']['I_R'])
        rat['Ha/Hb'] = RS('HaNII:Ha', 'Hb:Hb')
        if res['OII:OII3726']['n_groups'] >= 2:
            bO = boots['OII:OII3726'] + boots['OII:OII3729']
            # [O II] only from BOSS/SDSS-V groups: compare with [O III] from the same groups
            m_o2, _, _, _ = st['OII:OII3726']
            gb = gids[np.isfinite(m_o2[gids])]
            o3b = combine(st['OIII:OIII5007'][0] * fac, w, gb)
            o3bb = np.array([combine(st['OIII:OIII5007'][0] * fac, w, RNG.choice(gb, len(gb))) for _ in range(NB)])
            # paired bootstrap for the [O II]/[O III] ratio on the BOSS subset
            pb = []
            for _ in range(NB):
                ii = RNG.choice(gb, len(gb))
                pb.append((combine(st['OII:OII3726'][0] * fac, w, ii) + combine(st['OII:OII3729'][0] * fac, w, ii)) /
                          combine(st['OIII:OIII5007'][0] * fac, w, ii))
            pb = np.array(pb)
            o2v = res['OII:OII3726']['I_R'] + res['OII:OII3729']['I_R']
            rat['OII/OIII5007_photon_BOSS'] = {'value': float(o2v / o3b), 'p2.5': float(np.nanpercentile(pb, 2.5)),
                                              'p16': float(np.nanpercentile(pb, 16)), 'p84': float(np.nanpercentile(pb, 84)),
                                              'p97.5': float(np.nanpercentile(pb, 97.5)),
                                              'OIII5007_same_groups_R': float(o3b), 'n_groups_boss': int(len(gb))}
            rat['OII/OIII5007_energy_BOSS'] = {k: (v * 5008.24 / 3728.5 if isinstance(v, float) and k != 'OIII5007_same_groups_R' else v)
                                              for k, v in rat['OII/OIII5007_photon_BOSS'].items()}
            rat['OII3729/3726'] = RS('OII:OII3729', 'OII:OII3726')
        rat['OIII4363/5007'] = RS('B4363:OIII4363', 'OIII:OIII5007')
        rat['HeII4686/OIII5007'] = RS('HeII:HeII4686', 'OIII:OIII5007')
        res['ratios_photon'] = rat
        # velocity
        if nm in ('all', 'south_dec_lt20', 'north_dec_ge20', 'W_strand_15_20', 'E_strand_15_20', 'south_dec_lt17p6'):
            try:
                vres, (lm, sp) = oiii_velocity(fib, grp, gok, sel, off & ~sel, w, gids, D, zp)
                res['velocity_OIII'] = vres
                np.save(ROOT / f'data/derived/verify_sdss_stack_{tag}_r{r_on}_{nm}.npy', np.vstack([lm, sp]))
                for ins in ('SDSS', 'BOSS'):
                    gi = gids[grp['instrument'][gids] == ins]
                    if len(gi) >= 3:
                        v2, _ = oiii_velocity(fib, grp, gok, sel, off & ~sel, w, gi, D, zp, nboot=200)
                        res[f'velocity_OIII_{ins}'] = v2
            except Exception as ex:
                res['velocity_OIII'] = {'error': repr(ex)}
        out['regions'][nm] = res
        o = res['OIII:OIII5007']; hb = res['Hb:Hb']; ha = res['HaNII:Ha']
        print(f"{nm:28s} ng={res['n_groups']:3d} non={res['n_on_fibres']:4d} [OIII]={o['I_R']:.3f}±{o['err_boot']:.3f} "
              f"Hb={hb['I_R']:.3f}±{hb['err_boot']:.3f} Ha={ha['I_R']:.3f}±{ha['err_boot']:.3f} "
              f"NII={res['HaNII:NII6583']['I_R']:.3f}±{res['HaNII:NII6583']['err_boot']:.3f} "
              f"SII={res['SII:SII6716']['I_R'] + res['SII:SII6731']['I_R']:.3f}  "
              f"O3/Hb>{rat['OIII5007/Hb_lower2sig']:.1f} N2/Ha={rat['NII6583/Ha']['value']:.2f} S2/Ha={rat['SII/Ha']['value']:.2f}"
              + (f" v={res['velocity_OIII']['v_corr']:+.1f}±{res['velocity_OIII']['ev_boot']:.1f}" if 'velocity_OIII' in res and 'v_corr' in res['velocity_OIII'] else ''),
              flush=True)
    json.dump(out, open(ROOT / f'data/derived/verify_sdss_results_{tag}_r{r_on}.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], float(a[1]), *(float(x) for x in a[2:4]))
