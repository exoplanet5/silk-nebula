"""w3_oii_sdss.py -- (wave 3, tag oii) SDSS/BOSS [O II] 3729/3726 instrumental control and ON-OFF re-fits.

(A) ABSOLUTE total sky (flux + sky) of BOSS/SDSS-V 2" sky fibres, per dark plate-spectrograph group (wave-2 groups):
    the Galactic WIM [O II] doublet that every sightline carries must sit near the low-density limit (1.45-1.50).
    Model, 3700-3760 A: CALSPEC solar model broadened to the group's wdisp LSF (zodiacal light + scattered starlight;
    free Doppler shift on a grid) x (a + b x) + quadratic additive pseudo-continuum (Herzberg O2 airglow etc.)
    + two point-sampled Gaussians at V_helio = +12 km/s (variants: free v; sigma = wdisp x 0.90..1.30 or free).
    Per-group ratios are combined with inverse-variance weights; errors from a bootstrap over groups.
    Airglow test: [O II] intensity vs solar depression / UT hour (groups with TAI).
(B) ON-OFF stacks rebuilt from the total sky in the [O II], [S II] and [O III] windows (wave-2 ON/OFF definitions:
    ON within 0.15 deg of the nearest strand, OFF 0.45-2.0 deg in the same group with the Halpha-patch 30 % contour
    excluded; regions south/north/all/patch), weighted by the wave-2 [O III] group errors, then fitted with the
    same model variants (sigma from the stack's mean wdisp).
(C) ABSOLUTE total-sky [S II] 6716/6731 (WIM) per group: red-arm control.
NOTE (result): the direct absolute fits (A, C) FAIL as a control -- the CALSPEC solar template does not match the
total-sky Fraunhofer structure at R~1600 (median chi2/dof ~ 23, solar Doppler shift runs to the grid edges), and OH
lines were not modelled in (C). They are kept only as a diagnostic (DO_ABS = True). The low-density control is done
by regression on WHAM Halpha instead: scripts/w3_oii_sdss_wim.py.
Fits of the ON-OFF stacks use uniform pixel weights (a GLS with the bootstrap covariance would be fitted to the very
replicates that define the intervals). Variant 'eps': an extra regressor eps x (stacked OFF total sky), i.e. an ON/OFF
throughput or continuum mismatch that carries the WIM [O II] doublet and the solar Fe I 3727.6 line.
Output: data/derived/w3_oii_sdss.json, data/derived/w3_oii_sdss_stacks.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w3_oii_sdss.py
"""
import csv
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import load_windows, flam_to_R_per_A  # noqa: E402
from verify_sdss_groups import sun_on, clipped_mean  # noqa: E402  (solar template convolution, clipped mean)
from w3_oii_common import L, C, design, grid_fit, chi2_curve, shrink_cov  # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula/'
DO_ABS = False
RNG = np.random.default_rng(20261118)
O2 = [L['OII3726'], L['OII3729']]
S2 = [L['SII6731'], L['SII6716']]   # order chosen so that ratio = line2/line1 = 6716/6731
V_WIM = 12.0
FR = np.arange(0.80, 1.50001, 0.01)
VG = np.arange(-30.0, 50.01, 1.0)


def read_csv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))


def fcol(rows, k):
    out = []
    for r in rows:
        try:
            out.append(float(r[k]))
        except (ValueError, TypeError):
            out.append(np.nan)
    return np.array(out)


def good_groups(G):
    sb = fcol(G, 'sb5070_AB_cal')
    nf = fcol(G, 'n_fib')
    ins = np.array([r['instrument'] for r in G])
    ok = (nf >= 6) & np.isfinite(sb)
    med = {i: np.nanmedian(sb[ok & (ins == i)]) for i in ('SDSS', 'BOSS')}
    ref = np.array([med[i] for i in ins])
    return ok & (sb > ref - 0.75), ins


def abs_fit(lam, y, var, sig, lines, lo, hi, vwim=V_WIM, extra_free=True):
    """Absolute total-sky fit for one group. Returns dict of variants."""
    m = (lam > lo) & (lam < hi) & np.isfinite(y) & np.isfinite(var) & (var > 0) & np.isfinite(sig)
    x, yy, w = lam[m], y[m], 1 / var[m]
    s0 = float(np.nanmedian(sig[m]))
    Wi = np.diag(w)
    xc = (x - x.mean()) / 30.0
    # solar template on a Doppler grid (zodiacal light); keep the best per data set
    suns = {vs: sun_on(x, sig[m], v=vs) for vs in np.arange(-60, 61, 10.0)}

    def builders(vgrid, sgrid):
        for vs, sv in suns.items():
            for v in vgrid:
                for s in sgrid:
                    X = design(x, lines, v, [s, s], npoly=0)
                    B = np.vstack([X, sv, sv * xc, np.ones_like(x), xc, xc ** 2])
                    yield (vs, v, s), B
    out = {}
    Ys = yy[None, :]
    for tag, vg, sg in (('fixed', [vwim], [s0]), ('s0.90', [vwim], [0.9 * s0]), ('s1.10', [vwim], [1.1 * s0]),
                        ('s1.20', [vwim], [1.2 * s0]), ('vfree', VG[::2], [s0]),
                        ('sfree', [vwim], list(FR * s0)), ('vsfree', VG[::4], list(FR[::2] * s0))):
        if not extra_free and tag in ('vfree', 'vsfree', 'sfree'):
            continue
        c2, par, p = grid_fit(Ys, Wi, builders(vg, sg))
        cov = None
        # amplitude errors at the best parameters
        vs, v, s = par[0]
        X = design(x, lines, v, [s, s], npoly=0)
        B = np.vstack([X, suns[vs], suns[vs] * xc, np.ones_like(x), xc, xc ** 2])
        A = B @ Wi @ B.T
        cov = np.linalg.inv(A)
        out[tag] = dict(I1=float(p[0][0]), I2=float(p[0][1]), e1=float(np.sqrt(cov[0, 0])), e2=float(np.sqrt(cov[1, 1])),
                        c12=float(cov[0, 1]), v=float(v), sigma=float(s), s_over_wdisp=float(s / s0), v_sun=float(vs),
                        chi2=float(c2[0]), dof=int(x.size - B.shape[0]))
    return out, s0


def combine_ratio(rows, key, nboot=2000):
    I1 = np.array([r[key]['I1'] for r in rows])
    I2 = np.array([r[key]['I2'] for r in rows])
    e1 = np.array([r[key]['e1'] for r in rows])
    e2 = np.array([r[key]['e2'] for r in rows])
    w = 1 / (e1 ** 2 + e2 ** 2)
    ok = np.isfinite(I1) & np.isfinite(I2) & np.isfinite(w)
    I1, I2, w = I1[ok], I2[ok], w[ok]
    R = np.sum(w * I2) / np.sum(w * I1)
    bs = []
    for _ in range(nboot):
        i = RNG.integers(0, I1.size, I1.size)
        bs.append(np.sum(w[i] * I2[i]) / np.sum(w[i] * I1[i]))
    q = np.percentile(bs, [2.5, 16, 50, 84, 97.5])
    per = I2 / I1
    return dict(ratio=float(R), p2p5=float(q[0]), p16=float(q[1]), p84=float(q[3]), p97p5=float(q[4]), n=int(I1.size),
                median_group_ratio=float(np.median(per)), I_sum_mean_cal=float(np.sum(w * (I1 + I2)) / np.sum(w)))


def main():
    D = load_windows()
    lam = D['lam']
    F = read_csv(ROOT + 'data/derived/verify_sdss_fibres.csv')
    G = read_csv(ROOT + 'data/derived/verify_sdss_groups.csv')
    gid = np.array([int(r['gid']) for r in F])
    dmin = fcol(F, 'dmin')
    dec = fcol(F, 'dec')
    in30 = np.array([int(r['in_patch30']) for r in F]).astype(bool)
    inhm = np.array([int(r['in_patch']) for r in F]).astype(bool)
    diam = fcol(F, 'diam')
    gok, ins = good_groups(G)
    boss = gok & (ins == 'BOSS')
    TOT, IV, WD = D['tot'].astype(float), D['ivar'].astype(float), D['wdisp'].astype(float)
    res = {'n_boss_dark_groups': int(boss.sum())}
    # ------------------------------------------------------------- (A)+(C) absolute total sky per group
    # LSF statistics of BOSS dark-group sky fibres
    selb = np.isin(gid, np.where(boss)[0])
    i3 = int(np.argmin(np.abs(lam - 3728.0)))
    i6 = int(np.argmin(np.abs(lam - 6724.0)))
    s3 = WD[selb, i3] * lam[i3] * np.log(10) * 1e-4
    s6 = WD[selb, i6] * lam[i6] * np.log(10) * 1e-4
    s3, s6 = s3[np.isfinite(s3) & (s3 > 0)], s6[np.isfinite(s6) & (s6 > 0)]
    lsf = dict(boss_wdisp_sigma_3728_median=float(np.median(s3)), p16_3728=float(np.percentile(s3, 16)),
               p84_3728=float(np.percentile(s3, 84)), boss_wdisp_sigma_6724_median=float(np.median(s6)),
               pixel_A_3728=float(lam[i3 + 1] - lam[i3]), n=int(s3.size))
    json.dump(lsf, open(ROOT + 'data/derived/w3_oii_sdss_lsf.json', 'w'), indent=1)
    print('BOSS LSF', lsf)
    rows_o2, rows_s2, meta = [], [], []
    gstack = {'o2_resid': [], 'o2_w': [], 's2_resid': []}
    tai = fcol(G, 'tai')
    for g in (np.where(boss)[0] if DO_ABS else []):
        sel = np.where(gid == g)[0]
        dm = diam[sel[0]]
        X = flam_to_R_per_A(TOT[sel], lam, dm)
        V = np.where(IV[sel] > 0, flam_to_R_per_A(1.0, lam, dm) ** 2 / np.where(IV[sel] > 0, IV[sel], 1), np.inf)
        M, MV, n = clipped_mean(X, V)
        sig = np.nanmedian(WD[sel], 0) * lam * np.log(10) * 1e-4
        try:
            fo, so = abs_fit(lam, M, MV, sig, O2, 3700.0, 3760.0)
            fs, ss = abs_fit(lam, M, MV, sig, S2, 6700.0, 6752.0, extra_free=False)
        except np.linalg.LinAlgError:
            continue
        rows_o2.append(fo)
        rows_s2.append(fs)
        meta.append(dict(gid=int(g), plate=G[g]['plate'], mjd=G[g]['mjd'], sp=G[g]['sp'], n_fib=len(sel),
                         wdisp_sigma_3728=so, wdisp_sigma_6724=ss, tai=float(tai[g]), dec=float(fcol([G[g]], 'dec')[0]),
                         ra=float(fcol([G[g]], 'ra')[0]), I5577=float(fcol([G[g]], 'I_OI5577_cal')[0])))
    if DO_ABS:
        res['abs_groups'] = [dict(meta[i], o2=rows_o2[i], s2=rows_s2[i]) for i in range(len(meta))]
        res['abs_o2'] = {k: combine_ratio(rows_o2, k) for k in rows_o2[0]}
        res['abs_s2'] = {k: combine_ratio(rows_s2, k) for k in rows_s2[0]}
        for k, v in res['abs_o2'].items():
            print('ABS [O II] %-7s 3729/3726 = %.3f [%.3f, %.3f] (95%% %.3f-%.3f) n=%d  median group %.3f  <I>=%.2f R cal' % (
                k, v['ratio'], v['p16'], v['p84'], v['p2p5'], v['p97p5'], v['n'], v['median_group_ratio'],
                v['I_sum_mean_cal']))
        sfr = np.array([r['sfree']['s_over_wdisp'] for r in rows_o2])
        vfr = np.array([r['vfree']['v'] for r in rows_o2])
        res['abs_o2_sfree_sigma_over_wdisp'] = dict(median=float(np.median(sfr)), p16=float(np.percentile(sfr, 16)),
                                                    p84=float(np.percentile(sfr, 84)))
        res['abs_o2_vfree_v'] = dict(median=float(np.median(vfr)), p16=float(np.percentile(vfr, 16)),
                                     p84=float(np.percentile(vfr, 84)))
        print('sfree sigma/wdisp', res['abs_o2_sfree_sigma_over_wdisp'], 'vfree v', res['abs_o2_vfree_v'])
        for k, v in res['abs_s2'].items():
            print('ABS [S II] %-7s 6716/6731 = %.3f [%.3f, %.3f] n=%d  <I>=%.2f R cal' % (
                k, v['ratio'], v['p16'], v['p84'], v['n'], v['I_sum_mean_cal']))
        # airglow test: solar depression from TAI where available
        try:
            from astropy.time import Time
            from astropy.coordinates import EarthLocation, AltAz, get_sun
            import astropy.units as u
            apo = EarthLocation.from_geodetic(lon=-105.8202 * u.deg, lat=32.7802 * u.deg, height=2788 * u.m)
            tt = np.array([m['tai'] for m in meta])
            okt = np.isfinite(tt) & (tt > 0)
            dep = np.full(len(meta), np.nan)
            if okt.any():
                t = Time(tt[okt] / 86400.0, format='mjd', scale='tai')
                alt = get_sun(t).transform_to(AltAz(obstime=t, location=apo)).alt.deg
                dep[okt] = -alt
            Iso = np.array([r['fixed']['I1'] + r['fixed']['I2'] for r in rows_o2])
            res['airglow_test'] = dict(n_with_tai=int(okt.sum()),
                                       corr_I_vs_depression=float(np.corrcoef(Iso[okt], dep[okt])[0, 1]) if okt.sum() > 3 else None)
            for i, m in enumerate(meta):
                m['solar_depression'] = float(dep[i])
        except Exception as ex:  # noqa: BLE001
            res['airglow_test'] = {'error': repr(ex)}
        # time-of-night proxy for all groups: correlation of [O II] with [O I]5577 (airglow) and with the 3729/3726 ratio
        Iso = np.array([r['fixed']['I1'] + r['fixed']['I2'] for r in rows_o2])
        rat = np.array([r['fixed']['I2'] / r['fixed']['I1'] for r in rows_o2])
        I5577 = np.array([m['I5577'] for m in meta])
        ok = np.isfinite(I5577) & np.isfinite(Iso)
        res['airglow_test']['corr_IOII_vs_I5577'] = float(np.corrcoef(Iso[ok], I5577[ok])[0, 1])
        res['airglow_test']['corr_ratio_vs_I5577'] = float(np.corrcoef(rat[ok], I5577[ok])[0, 1])
        print('airglow test', res['airglow_test'])
    # ------------------------------------------------------------- (B) ON-OFF stacks
    Zd = np.load(ROOT + 'data/derived/verify_sdss_diff_base.npz')
    a3, e3 = Zd['OIII__OIII5007'], Zd['OIII__OIII5007__err']
    off = (dmin > 0.45) & (dmin < 2.0) & ~in30
    ON = dmin < 0.15
    regions = {'all': ON, 'south': ON & (dec < 20), 'north': ON & (dec >= 20), 'south_lt17p6': ON & (dec < 17.6),
               'patch30': in30, 'patch_halfmax': inhm, 'north_excl_patch30': ON & (dec >= 20) & ~in30}
    wins = {'o2': (3690.0, 3770.0), 's2': (6690.0, 6770.0), 'o3': (4935.0, 5032.0)}
    stacks = {}
    res['onoff'] = {}
    for rn, onm in regions.items():
        spec = {k: [] for k in wins}
        offm = {k: [] for k in wins}
        sigs = {k: [] for k in wins}
        wts, gl = [], []
        for g in np.where(boss)[0]:
            s = np.where(gid == g)[0]
            so, sf = s[onm[s]], s[off[s] & ~onm[s]]
            if len(so) < 1 or len(sf) < 3:
                continue
            # group weight from the wave-2 [O III] ON-OFF error (OFF scatter)
            af = a3[sf]
            sd = max(1.4826 * np.nanmedian(np.abs(af - np.nanmedian(af))), 0.5 * np.nanstd(af))
            eg = sd * np.sqrt(1 / len(so) + 1 / len(sf))
            if not np.isfinite(eg) or eg <= 0:
                continue
            dm = diam[s[0]]
            for k, (lo, hi) in wins.items():
                mm = (lam > lo) & (lam < hi)
                Xg = flam_to_R_per_A(TOT[s][:, mm], lam[mm], dm)
                Vg = np.where(IV[s][:, mm] > 0, flam_to_R_per_A(1.0, lam[mm], dm) ** 2 /
                              np.where(IV[s][:, mm] > 0, IV[s][:, mm], 1), np.inf)
                Mo, _ = clipped_mean(Xg[np.isin(s, sf)], Vg[np.isin(s, sf)])[:2]
                Xon = np.where(np.isfinite(Vg[np.isin(s, so)]), Xg[np.isin(s, so)], np.nan)
                spec[k].append(np.nanmean(Xon, 0) - Mo)
                offm[k].append(Mo)
                sigs[k].append(np.nanmedian(WD[so][:, mm], 0) * lam[mm] * np.log(10) * 1e-4)
            wts.append(1 / eg ** 2)
            gl.append(g)
        if len(wts) < 3:
            continue
        wts = np.array(wts)
        out = {'n_groups': len(wts)}
        for k, (lo, hi) in wins.items():
            mm = (lam > lo) & (lam < hi)
            A = np.array(spec[k])
            Sg = np.array(sigs[k])
            Mo_all = np.array(offm[k])
            good = np.all(np.isfinite(A), 1) & np.all(np.isfinite(Mo_all), 1)
            A, Sg, ww, Mo_all = A[good], Sg[good], wts[good], Mo_all[good]
            stacks[f'{rn}_{k}_off'] = np.sum(ww[:, None] * Mo_all, 0) / np.sum(ww)
            st = np.sum(ww[:, None] * A, 0) / np.sum(ww)
            bs = []
            for _ in range(400):
                i = RNG.integers(0, len(ww), len(ww))
                bs.append(np.sum(ww[i, None] * A[i], 0) / np.sum(ww[i]))
            bs = np.array(bs)
            stacks[f'{rn}_{k}_lam'] = lam[mm]
            stacks[f'{rn}_{k}_S'] = st
            stacks[f'{rn}_{k}_bs'] = bs.astype(np.float32)
            stacks[f'{rn}_{k}_sig'] = np.sum(ww[:, None] * Sg, 0) / np.sum(ww)
        res['onoff'][rn] = out
    np.savez_compressed(ROOT + 'data/derived/w3_oii_sdss_stacks.npz', **stacks)
    # fits of the ON-OFF stacks
    for rn in res['onoff']:
        for k, lines, lo, hi, l5 in (('o2', O2, 3708.0, 3750.0, 3728.0), ('s2', S2, 6702.0, 6750.0, 6724.0)):
            x = stacks[f'{rn}_{k}_lam']
            m = (x > lo) & (x < hi)
            x = x[m]
            S = stacks[f'{rn}_{k}_S'][m]
            bs = stacks[f'{rn}_{k}_bs'][:, m].astype(float)
            s0 = float(np.interp(l5, stacks[f'{rn}_{k}_lam'], stacks[f'{rn}_{k}_sig']))
            Wi = np.eye(x.size)
            Ys = np.vstack([S[None], bs])
            offt = stacks[f'{rn}_{k}_off'][m]
            offt = (offt - offt.mean()) / np.abs(offt).max()
            fits = {}
            # eps variant: extra regressor = stacked OFF total sky (throughput / continuum mismatch)
            for tag, vg in (('eps_fixed', [13.0]), ('eps_vfree', VG)):
                c2, par, p = grid_fit(Ys, Wi, (((v, s0), np.vstack([design(x, lines, v, [s0, s0], npoly=2), offt]))
                                               for v in vg))
                r = np.array([q[1] / q[0] for q in p])
                qq = np.nanpercentile(r[1:], [2.5, 16, 50, 84, 97.5])
                fits[tag] = dict(ratio=float(r[0]), p2p5=float(qq[0]), p16=float(qq[1]), p84=float(qq[3]),
                                 p97p5=float(qq[4]), v=float(par[0][0]), eps_coef=float(p[0][-1]))
            for tag, vg, sg in (('fixed', [13.0], [s0]), ('s0.90', [13.0], [0.9 * s0]), ('s1.10', [13.0], [1.1 * s0]),
                                ('s1.20', [13.0], [1.2 * s0]), ('vfree', VG, [s0]), ('sfree', [13.0], list(FR * s0)),
                                ('vsfree', VG[::2], list(FR[::2] * s0))):
                c2, par, p = grid_fit(Ys, Wi, (((v, s), design(x, lines, v, [s, s], npoly=2)) for v in vg for s in sg))
                r = np.array([q[1] / q[0] for q in p])
                qq = np.nanpercentile(r[1:], [2.5, 16, 50, 84, 97.5])
                fits[tag] = dict(ratio=float(r[0]), p2p5=float(qq[0]), p16=float(qq[1]), p84=float(qq[3]),
                                 p97p5=float(qq[4]), v=float(par[0][0]), sigma=float(par[0][1]),
                                 s_over_wdisp=float(par[0][1] / s0), I1=float(p[0][0]), I2=float(p[0][1]),
                                 chi2=float(c2[0]), npix=int(x.size),
                                 v_p16=float(np.percentile([q[0] for q in par[1:]], 16)),
                                 v_p84=float(np.percentile([q[0] for q in par[1:]], 84)),
                                 s_p16=float(np.percentile([q[1] for q in par[1:]], 16) / s0),
                                 s_p84=float(np.percentile([q[1] for q in par[1:]], 84) / s0))
            # profile likelihood (v, sigma free)
            rg = np.arange(0.6, 2.01, 0.02)
            prof = []
            for rr in rg:
                Xs = []
                for v in VG[::4]:
                    for f in FR[::5]:
                        X = design(x, lines, v, [f * s0, f * s0], npoly=2)
                        Xs.append(np.vstack([X[0] + rr * X[1], X[2:]]))
                prof.append(chi2_curve(S[None], Wi, Xs).min())
            prof = np.array(prof)
            prof = (prof - prof.min()) / max(prof.min() / max(x.size - 5, 1), 1e-9)   # chi2 scaled to chi2/dof = 1
            fits['profile'] = dict(r_best=float(rg[np.argmin(prof)]), r68=[float(rg[prof <= 1].min()), float(rg[prof <= 1].max())],
                                   r95=[float(rg[prof <= 4].min()), float(rg[prof <= 4].max())],
                                   dchi2_at_1p46=float(np.interp(1.46, rg, prof)))
            res['onoff'][rn][k] = fits
            print(f'ON-OFF {rn:18s} {k} ng={res["onoff"][rn]["n_groups"]:3d} ' + ' '.join(
                f'{t}={fits[t]["ratio"]:.2f}[{fits[t]["p16"]:.2f},{fits[t]["p84"]:.2f}]' for t in
                ('fixed', 's1.10', 's1.20', 'vfree', 'sfree', 'vsfree', 'eps_fixed')) +
                f' (vfree v={fits["vfree"]["v"]:+.0f} [{fits["vfree"]["v_p16"]:+.0f},{fits["vfree"]["v_p84"]:+.0f}];'
                f' sfree s/wd={fits["sfree"]["s_over_wdisp"]:.2f}; vsfree v={fits["vsfree"]["v"]:+.0f} s/wd={fits["vsfree"]["s_over_wdisp"]:.2f}) '
                f'prof {fits["profile"]["r68"]} d146={fits["profile"]["dchi2_at_1p46"]:.1f}', flush=True)
    json.dump(res, open(ROOT + 'data/derived/w3_oii_sdss.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
