"""verify_desi_ratios.py -- (c) independent differential line ratios of the [O III]-correlated component.

Estimator (instrumental variables, my own): for line L with per-fibre flux L_i (my fitter, fixed V_NEB) and the
chromatic fibre-acceptance correction c_i = f_i(lambda_L)/f_i(5008) (verify_desi_setup.py), the energy ratio of the
filament component is
    rho_L = sum u (Zf - <Zf>) (c L - <c L>) / sum u (Zf - <Zf>) (O - <O>),   u_i = 1/sigma(cL)_i^2,
with Zf_i = Z_i / kR_i the geometry-only template prediction in flux units (no data used) and <.> weighted means
(intercept). Variants: petal fixed effects (de-meaning within the dominant tile-petal); neighbour instrument (mean O
of fibres within 0.15 deg in OTHER tile-petals); subsets (QSO+z>0.5, smooth continua, dark, south/north, no Ha patch,
wave-1 fibres only). Errors: bootstrap over HEALPix nside=32 blocks (~1.8 deg), 400 resamples.
Hg tests: Hg I 4358 (topocentric, at v_corr) and [O I]5577 treated as 'lines' (should give rho = 0); 4363 fitted
alone and jointly with Hg; fake lines next to 4363 give the empirical null of rho at that wavelength; a pixel-level
IV regression spectrum around 4300-4400 A shows where the 4363 feature sits.
Output: data/derived/verify_desi_ratios.json, verify_desi_regspec.npz
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, C, WAVE, L, NEB_LINES, SKY_LINES, GAL_REST, load_meta, load_spectra, good_fit  # noqa

RNG = np.random.default_rng(11)
NB = 400
LINEDEF = {   # name: (group, column, chromatic key, rest vac A)
    'Hb': ('hb', 0, '4861'), 'OIII4959': ('o3_4959', 0, '5008'), 'OIII5007_single': ('o3_5007', 0, '5008'),
    'Hg': ('hg', 0, '4363'), 'Hd': ('hd', 0, '4363'), 'OIII4363': ('4363', 0, '4363'),
    'OIII4363_jointHg': ('4363hg', 0, '4363'), 'HgI4358_joint': ('4363hg', 1, '4363'),
    'HgI4358': ('hg4358', 0, '4363'), 'HeII4686': ('he2', 0, '4686'), 'OII3726': ('o2', 0, '3727'),
    'OII3729': ('o2', 1, '3727'), 'OII3727sum': ('o2sum', 0, '3727'), 'NeIII3869': ('ne3', 0, '3869'),
    'Ha': ('han2', 0, '6563'), 'NII6583': ('han2', 1, '6583'), 'SII6716': ('s2', 0, '6724'),
    'SII6731': ('s2', 1, '6724'), 'HeI5876': ('he1', 0, '5008'), 'NII5755': ('n2_5755', 0, '5008'),
    'ArIII7136': ('ar3', 0, '7136'), 'OII7325': ('o2_7325', 0, '7325'), 'SIII9069': ('s3', 0, '9069'),
    'sky5577': ('sky5577', 0, '5008'),
    'CaII7291': ('ca2_7291', 0, '7325'), 'FeII7155': ('fe2_7155', 0, '7136'), 'FeII4359': ('fe2_4359', 0, '4363'),
    'OII7325_jointCa': ('o2ca2', 0, '7325'), 'CaII7324_joint': ('o2ca2', 1, '7325'),
    'OIII4363_jointHgFe': ('hg_fe2_4363', 0, '4363'), 'HgI4358_jointFe': ('hg_fe2_4363', 1, '4363'),
    'FeII4359_jointHg': ('hg_fe2_4363', 2, '4363'), 'NI5200': ('n1', 0, '5008'),
}


def fake_defs(lf):
    out = {}
    for k in lf.files:
        if k.startswith('null') and k.endswith('__amp'):
            g = k[:-5]
            lam = float(g[4:])
            ck = min(['3727', '3869', '4363', '4686', '4861', '5008', '6563', '7136', '9069'],
                     key=lambda s: abs(float(s) - lam))
            out[g] = (g, 0, ck)
    return out


def iv(y, Zi, O, u, grp=None):
    """IV slope of y on O with instrument Zi, weights u; intercept or group fixed effects."""
    def dm(v):
        if grp is None:
            return v - np.sum(u * v) / np.sum(u)
        s = np.bincount(grp, u * v) / np.maximum(np.bincount(grp, u), 1e-300)
        return v - s[grp]
    zz = dm(Zi)
    return np.sum(u * zz * dm(y)), np.sum(u * zz * dm(O))


def boot_ratio(y, Zi, O, u, blocks, grp=None, CNT=None, UB=None):
    """CNT: [NB, n_blocks] resample counts shared by all lines of a subset (joint bootstrap); UB: block labels."""
    inv = np.searchsorted(UB, blocks)
    num, den = iv(y, Zi, O, u, grp)
    out = np.empty(NB)
    nums = np.empty(NB)
    dens = np.empty(NB)
    for b in range(NB):
        cnt = CNT[b][inv].astype(float)
        uu = u * cnt
        n_, d_ = iv(y, Zi, O, uu, grp)
        nums[b], dens[b] = n_, d_
        out[b] = n_ / d_
    return num / den, out, nums, dens


def summarize(r, bs):
    p = np.percentile(bs, [2.5, 16, 50, 84, 97.5])
    return dict(ratio=float(r), p2p5=float(p[0]), p16=float(p[1]), p50=float(p[2]), p84=float(p[3]),
                p97p5=float(p[4]), sig_boot=float(0.5 * (p[3] - p[1])))


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids = lf['sparcl_id']
    keep, meta = load_meta(ids)
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    scale = det['err_scale_from_nulls']
    import healpy as hp
    blocks = hp.ang2pix(32, meta['ra'], meta['dec'], lonlat=True)
    blocks64 = hp.ang2pix(64, meta['ra'], meta['dec'], lonlat=True)
    O, EO = lf['o3__amp'][:, 0], lf['o3__err'][:, 0] * scale
    Zf = cm['Z'] / cm['kR_5008']
    base = good_fit(lf, 'o3') & np.isfinite(Zf) & cm['ok_contrib']
    grp_lab = np.array([f'{t}_{p}' for t, p in zip(cm['group_tile'], cm['group_petal'])])
    _, grp = np.unique(grp_lab, return_inverse=True)
    # neighbour instrument: mean O (in SB units) of fibres within 0.15 deg in other tile-petals
    from scipy.spatial import cKDTree
    xyz = np.c_[np.cos(np.radians(meta['dec'])) * np.cos(np.radians(meta['ra'])),
                np.cos(np.radians(meta['dec'])) * np.sin(np.radians(meta['ra'])), np.sin(np.radians(meta['dec']))]
    tree = cKDTree(xyz[base])
    bi = np.where(base)[0]
    OSB, wSB = O * cm['kR_5008'], 1 / (EO * cm['kR_5008']) ** 2
    nbr = np.full(len(O), np.nan)
    for i, nb in enumerate(tree.query_ball_point(xyz, 2 * np.sin(np.radians(0.15) / 2))):
        js = bi[nb]
        js = js[grp_lab[js] != grp_lab[i]]
        if js.size >= 3:
            nbr[i] = np.sum(wSB[js] * OSB[js]) / np.sum(wSB[js])
    nbr_f = nbr / cm['kR_5008']
    ra, dec = meta['ra'], meta['dec']
    patch = (ra > 115.0) & (ra < 116.2) & (dec > 19.6) & (dec < 21.6)
    hiz = (meta['redshift'] > 0.5) | (meta['spectype'] == 'QSO')
    # smooth-continuum proxy: low continuum-fit reduced chi2 in the Hb and 4363 windows
    smooth = (lf['hb__rchi2'] < np.nanpercentile(lf['hb__rchi2'], 60)) & \
             (lf['4363__rchi2'] < np.nanpercentile(lf['4363__rchi2'], 60))
    subsets = {'all': base, 'petalFE': base, 'nbr_instrument': base & np.isfinite(nbr_f),
               'QSO_or_z>0.5': base & hiz, 'smooth_continuum': base & smooth,
               'dark': base & (meta['program'] == 'dark'), 'south_dec<19.5': base & (dec < 19.5),
               'north_dec>=19.5': base & (dec >= 19.5), 'no_Ha_patch': base & ~patch,
               'wave1_fibres_only': base & (lf['src'] == 'w1'), 'blocks_nside64': base}
    defs = dict(LINEDEF)
    defs.update(fake_defs(lf))
    res = {'err_scale': scale, 'n_base': int(base.sum()), 'v_neb': float(lf['v_neb'])}
    boots = {}
    for sname, m0 in subsets.items():
        rs = {}
        bl_all = blocks64 if sname == 'blocks_nside64' else blocks
        UB = np.unique(bl_all[m0])
        CNT = np.array([np.bincount(RNG.integers(0, UB.size, UB.size), minlength=UB.size) for _ in range(NB)])
        bsub = {}
        for name, (g, col, ck) in defs.items():
            if sname not in ('all', 'petalFE') and name.startswith('null') and sname != 'QSO_or_z>0.5':
                continue
            a = lf[g + '__amp'][:, col]
            e = lf[g + '__err'][:, col] * scale
            c = cm[f'f_chrom_{ck}'] / cm['f_chrom_5008']
            y, ey = a * c, e * c
            m = m0 & good_fit(lf, g, col) & np.isfinite(c)
            if m.sum() < 50:
                continue
            u = 1 / ey[m] ** 2
            Zi = nbr_f[m] if sname == 'nbr_instrument' else Zf[m]
            gg = grp[m] if sname == 'petalFE' else None
            if gg is not None:
                _, gg = np.unique(gg, return_inverse=True)
            bl = blocks64[m] if sname == 'blocks_nside64' else blocks[m]
            r, bs, nums, dens = boot_ratio(y[m], Zi, O[m], u, bl, gg, CNT, UB)
            rs[name] = summarize(r, bs)
            bsub[name] = bs
            rs[name]['n'] = int(m.sum())
            # achromatic variant (no f(lambda) correction) for the main run
            if sname == 'all' and not name.startswith('null'):
                num0, den0 = iv(a[m], Zi, O[m], 1 / e[m] ** 2)
                rs[name]['ratio_achromatic'] = float(num0 / den0)
                c2 = cm[f'f_half_{ck}'] / cm['f_half_5008']
                num1, den1 = iv(a[m] * c2[m], Zi, O[m], 1 / (e[m] * c2[m]) ** 2)
                rs[name]['ratio_halfchrom'] = float(num1 / den1)
        boots[sname] = bsub
        res[sname] = rs
        print(sname, {k: round(v['ratio'], 4) for k, v in rs.items() if not k.startswith('null')}, flush=True)
    # empirical null of rho at 4300-4400 A (fake lines) vs 4363
    fk = [k for k in res['all'] if k.startswith('null43') or k.startswith('null439') or k.startswith('null438')]
    res['fake_4363_neighbourhood'] = {k: res['all'][k] for k in fk}
    vals = np.array([res['all'][k]['ratio'] for k in fk])
    sigs = np.array([res['all'][k]['sig_boot'] for k in fk])
    res['fake_4363_summary'] = dict(rms_ratio=float(np.sqrt(np.mean(vals ** 2))), mean_sig_boot=float(sigs.mean()),
                                    chi_rms=float(np.sqrt(np.mean((vals / sigs) ** 2))), n=len(fk))
    allfk = [k for k in res['all'] if k.startswith('null')]
    v2 = np.array([res['all'][k]['ratio'] / res['all'][k]['sig_boot'] for k in allfk])
    res['fake_all_chi'] = dict(n=len(allfk), chi_rms=float(np.sqrt(np.mean(v2 ** 2))), chi=v2.tolist(),
                               names=allfk)
    res['derived'] = {sn: derived(res[sn], boots[sn]) for sn in ('all', 'south_dec<19.5', 'petalFE', 'QSO_or_z>0.5',
                                                                'north_dec>=19.5', 'nbr_instrument')}
    np.savez_compressed(D / 'verify_desi_ratio_boots.npz',
                        **{f'{sn}|{k}': v for sn, d in boots.items() for k, v in d.items()})
    (D / 'verify_desi_ratios.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res['derived'], indent=1))
    print('fake 4363 neighbourhood', json.dumps(res['fake_4363_summary']))
    print('all fake chi', res['fake_all_chi']['chi_rms'], np.round(v2, 2))


def derived(A, b):
    """Derived diagnostics from one subset's ratios A and its JOINT bootstrap arrays b (energy units)."""
    def pct(x, q):
        return float(np.percentile(x, q))
    der = {}
    hb = b['Hb']
    der['OIII5007/Hb'] = dict(value=float(1 / A['Hb']['ratio']),
                              lower_97p5=float(1 / pct(hb, 97.5)) if pct(hb, 97.5) > 0 else None,
                              lower_84=float(1 / pct(hb, 84)) if pct(hb, 84) > 0 else None,
                              hb_over_o3_p16_p84=[pct(hb, 16), pct(hb, 84)])
    o2 = b['OII3726'] + b['OII3729']
    der['OII3727total/OIII5007'] = dict(value=float(A['OII3726']['ratio'] + A['OII3729']['ratio']),
                                        p16=pct(o2, 16), p84=pct(o2, 84), p2p5=pct(o2, 2.5), p97p5=pct(o2, 97.5))
    q = hb / (o2 + 1 + b['OIII4959'])
    der['R23'] = dict(value=float((A['OII3726']['ratio'] + A['OII3729']['ratio'] + 1 + A['OIII4959']['ratio']) /
                                  A['Hb']['ratio']),
                      lower_84=float(1 / pct(q, 84)) if pct(q, 84) > 0 else None,
                      lower_97p5=float(1 / pct(q, 97.5)) if pct(q, 97.5) > 0 else None)
    for nm in ('OIII4363', 'OIII4363_jointHg', 'OIII4363_jointHgFe'):
        if nm not in A:
            continue
        rr = A[nm]['ratio']
        t = (1 + A['OIII4959']['ratio']) / rr if rr > 0 else np.inf
        te = 32900 / np.log(t / 7.90) if t > 7.9 else np.inf
        tb = (1 + b['OIII4959']) / np.where(b[nm] > 0, b[nm], np.nan)
        teb = np.where(tb > 7.9, 32900 / np.log(tb / 7.90), np.inf)
        teb = np.where(np.isfinite(tb), teb, 0.0)       # 4363 <= 0 -> no temperature (cold)
        der[f'Te_{nm}'] = dict(R_O3=float(t), Te_K=float(te), ratio=rr, sig_boot=A[nm]['sig_boot'],
                               Te_p16=pct(teb, 16), Te_p84=pct(teb, 84), Te_p2p5=pct(teb, 2.5),
                               frac_boot_ratio_le0=float(np.mean(b[nm] <= 0)))
    r3 = b['OII3729'] / b['OII3726']
    der['OII3729/3726'] = dict(value=float(A['OII3729']['ratio'] / A['OII3726']['ratio']), p16=pct(r3, 16),
                               p84=pct(r3, 84), p2p5=pct(r3, 2.5), p97p5=pct(r3, 97.5))
    r4 = b['OIII5007_single'] / b['OIII4959']
    der['OIII5007/4959'] = dict(value=float(A['OIII5007_single']['ratio'] / A['OIII4959']['ratio']),
                                p16=pct(r4, 16), p84=pct(r4, 84))
    der['HeII/OIII_upper97p5'] = A['HeII4686']['p97p5']
    der['HeII/Hb'] = dict(value=float(A['HeII4686']['ratio'] / A['Hb']['ratio']))
    ha = b['Ha']
    der['Ha/Hb'] = dict(value=float(A['Ha']['ratio'] / A['Hb']['ratio']), p16=pct(ha / hb, 16), p84=pct(ha / hb, 84))
    for ln in ('NII6583', 'SII'):
        num = b['NII6583'] if ln == 'NII6583' else b['SII6716'] + b['SII6731']
        val = A['NII6583']['ratio'] if ln == 'NII6583' else A['SII6716']['ratio'] + A['SII6731']['ratio']
        der[f'{ln}/Ha'] = dict(value=float(val / A['Ha']['ratio']), p16=pct(num / ha, 16), p84=pct(num / ha, 84))
        hbased = 2.86 * hb
        der[f'{ln}/Ha_from_Hb'] = dict(value=float(val / (2.86 * A['Hb']['ratio'])), p16=pct(num / hbased, 16),
                                       p84=pct(num / hbased, 84))
    if 'CaII7291' in A:
        der['CaII7291/OII7325'] = dict(value=float(A['CaII7291']['ratio'] / A['OII7325_jointCa']['ratio']))
    return der


if __name__ == '__main__':
    main()
