"""verify_desi_detect.py -- (a) independent detection test of the Silk [O III] ridge in DESI DR1 residuals.

Statistic: matched filter of the per-fibre [O III] residual O_i (my fitter, doublet tied, fixed V_NEB) against a
GEOMETRY-ONLY template Z_i = T(x_i) - sum_t w_t <T>_{tile-petal} (unit-amplitude Gaussian ridges on the E/W/W2
trace, verify_desi_common.Template; no amplitude fitted to the data):
    A = sum w Z O / sum w Z^2,   SNR = A * sqrt(sum w Z^2),   w = 1/sigma^2,
with sigma scaled by the MAD of the null-line fluxes. Null distributions:
  1. rigid random rotations + translations of the trace on the sphere (300 draws; draws whose Z correlates with
     the true Z by |r| > 0.2 are set aside);
  2. label shuffling of O within the dominant tile-petal group (1000 permutations);
  3. the same statistic for 12 fake lines at line-free wavelengths (true template);
  4. RA/Dec translations on a grid (-4..+4 deg);
  5. splits: by Dec band, program, survey, single- vs multi-tile coadds, spectype; and jackknife over nside-4 HEALPix.
Output: data/derived/verify_desi_detect.json, verify_desi_detect_nulls.npz
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (D, load_meta, Exposures, Coadds, Template, petal_subtracted, rotmat, unitvec,  # noqa
                                good_fit)

RNG = np.random.default_rng(2026)
CENTRE = (115.8, 19.5)


def mf(O, E, Z, m):
    w = 1 / E[m] ** 2
    s2 = np.sum(w * Z[m] ** 2)
    A = np.sum(w * Z[m] * O[m]) / s2
    return float(A), float(A * np.sqrt(s2))


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids = lf['sparcl_id']
    assert (cm['sparcl_id'] == ids).all()
    keep, meta = load_meta(ids)
    O, E = lf['o3__amp'][:, 0], lf['o3__err'][:, 0]
    nulls = [k[:-5] for k in lf.files if k.startswith('null') and k.endswith('__amp')]
    chis = []
    for g in nulls:
        a, e = lf[g + '__amp'][:, 0], lf[g + '__err'][:, 0]
        k = np.isfinite(a) & (e > 0)
        chis.append(a[k] / e[k])
    scale = float(np.median([1.4826 * np.median(np.abs(c - np.median(c))) for c in chis]))
    kR = cm['kR_5008']
    Z = cm['Z']
    ok = np.isfinite(O) & np.isfinite(E) & (E > 0) & np.isfinite(Z) & np.isfinite(kR) & cm['ok_contrib']
    ok &= good_fit(lf, 'o3')                 # gross fit failures (masked pixels, huge errors)
    Es = E * scale
    out = dict(n=int(ok.sum()), err_scale_from_nulls=scale, n_null_lines=len(nulls))
    A0, S0 = mf(O, Es, Z, ok)
    clip = ok & (np.abs(O / Es) < 8)
    A0c, S0c = mf(O, Es, Z, clip)
    out['true'] = dict(A_1e17=A0, snr=S0, snr_clip8=S0c, A_clip8=A0c,
                       A_R_median_kR=float(A0 * np.nanmedian(kR[ok])))
    # ---- 3. fake lines with the true template
    fl = {}
    for g in nulls:
        a, e = lf[g + '__amp'][:, 0], lf[g + '__err'][:, 0] * scale
        m = ok & good_fit(lf, g)
        fl[g] = mf(a, e, Z, m)[1]
    out['fake_line_snr'] = fl
    # ---- 2. shuffle within dominant tile-petal
    grp = np.array([f'{t}_{p}' for t, p in zip(cm['group_tile'], cm['group_petal'])])
    idx = np.where(ok)[0]
    ug, inv = np.unique(grp[idx], return_inverse=True)
    members = [idx[inv == k] for k in range(ug.size)]
    sh = []
    for it in range(1000):
        perm = np.arange(len(O))
        for mem in members:
            if mem.size > 1:
                perm[mem] = mem[RNG.permutation(mem.size)]
        sh.append(mf(O[perm], Es[perm], Z, ok)[1])
    sh = np.array(sh)
    out['shuffle_within_petal'] = dict(n=1000, mean=float(sh.mean()), std=float(sh.std()), max=float(sh.max()),
                                       min=float(sh.min()), n_groups=int(ug.size),
                                       z_true=float((S0 - sh.mean()) / sh.std()))
    # ---- 1. rigid random rotations + translations
    ex = Exposures()
    co = Coadds(meta, ex)
    c0 = unitvec(*CENTRE)[0]
    rot, rr = [], []
    for it in range(300):
        R1 = rotmat(c0, RNG.uniform(0, 360))
        # translate: rotate the centre towards a random direction by up to 4 deg
        ang = RNG.uniform(0, 2 * np.pi)
        amp = RNG.uniform(0, 4.0)
        east = np.cross([0, 0, 1], c0)
        east /= np.linalg.norm(east)
        north = np.cross(c0, east)
        axis = np.cross(c0, np.cos(ang) * east + np.sin(ang) * north)
        R = rotmat(axis, amp) @ R1
        Zr = petal_subtracted(Template(R=R), meta, co)
        m = ok & np.isfinite(Zr)
        if np.sum(Zr[m] ** 2) <= 0:
            continue
        r = np.corrcoef(Zr[m], Z[m])[0, 1]
        A, S = mf(O, Es, Zr, m)
        Ac, Sc = mf(O, Es, Zr, m & clip)
        rot.append((S, Sc, r, np.sqrt(np.sum(Zr[m] ** 2 / Es[m] ** 2)), amp))
        rr.append(R.ravel())
    rot = np.array(rot)
    clean = np.abs(rot[:, 2]) < 0.2
    out['rigid_rotations'] = dict(n=int(len(rot)), n_clean=int(clean.sum()),
                                  snr_clean_mean=float(rot[clean, 0].mean()),
                                  snr_clean_std=float(rot[clean, 0].std()),
                                  snr_clean_max=float(rot[clean, 0].max()),
                                  snr_clean_p99=float(np.percentile(rot[clean, 0], 99)),
                                  snr_clip_clean_std=float(rot[clean, 1].std()),
                                  snr_clip_clean_max=float(rot[clean, 1].max()),
                                  z_true_vs_clean=float((S0 - rot[clean, 0].mean()) / rot[clean, 0].std()),
                                  z_true_clip_vs_clean=float((S0c - rot[clean, 1].mean()) / rot[clean, 1].std()),
                                  frac_clean_ge_true=float(np.mean(rot[clean, 0] >= S0)),
                                  snr_vs_corr_slope=float(np.polyfit(rot[:, 2], rot[:, 0], 1)[0]))
    # ---- 4. translations on a grid
    tr = {}
    for dra in (-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4):
        R = rotmat([0, 0, 1], dra)
        Zr = petal_subtracted(Template(R=R), meta, co)
        m = ok & np.isfinite(Zr)
        tr[f'dRA{dra:+g}'] = dict(snr=mf(O, Es, Zr, m)[1], snr_clip8=mf(O, Es, Zr, m & clip)[1],
                                  corr_with_true=float(np.corrcoef(Zr[m], Z[m])[0, 1]))
    out['ra_translations'] = tr
    # ---- 5. splits and jackknife
    dec = meta['dec']
    spl = {}
    for lab, m in [('dec<15.5', dec < 15.5), ('15.5-17.5', (dec >= 15.5) & (dec < 17.5)),
                   ('17.5-19.5', (dec >= 17.5) & (dec < 19.5)), ('19.5-22', (dec >= 19.5) & (dec < 22)),
                   ('22-25', (dec >= 22) & (dec < 25)), ('>=25', dec >= 25),
                   ('dark', meta['program'] == 'dark'), ('bright', meta['program'] == 'bright'),
                   ('main', meta['survey'] == 'main'), ('sv', meta['survey'] != 'main'),
                   ('single_tile', cm['n_tiles'] == 1), ('multi_tile', cm['n_tiles'] > 1),
                   ('QSO', meta['spectype'] == 'QSO'), ('GALAXY', meta['spectype'] == 'GALAXY'),
                   ('z>0.5_or_QSO', (meta['redshift'] > 0.5) | (meta['spectype'] == 'QSO')),
                   ('wave1_fibres', lf['src'] == 'w1'), ('ext_fibres', lf['src'] == 'ext')]:
        mm = ok & m
        if mm.sum() > 30 and np.sum(Z[mm] ** 2) > 0:
            A, S = mf(O, Es, Z, mm)
            spl[lab] = dict(n=int(mm.sum()), A_1e17=A, snr=S, snr_clip8=mf(O, Es, Z, mm & clip)[1])
    out['splits'] = spl
    import healpy as hp
    hpx = hp.ang2pix(4, meta['ra'], meta['dec'], lonlat=True)
    jk = {}
    for h in np.unique(hpx[ok]):
        mm = ok & (hpx != h)
        jk[int(h)] = dict(n_dropped=int((ok & (hpx == h)).sum()), snr=mf(O, Es, Z, mm)[1])
    out['jackknife_nside4'] = jk
    # tile jackknife: drop every coadd that includes a given tile
    tl = np.array([set(map(int, s.split(','))) if s else set() for s in cm['tiles']], dtype=object)
    alltiles = sorted(set().union(*tl[ok]))
    tj = {}
    for t in alltiles:
        has = np.array([t in s for s in tl])
        if (ok & has).sum() == 0:
            continue
        tj[int(t)] = dict(n_dropped=int((ok & has).sum()), snr=mf(O, Es, Z, ok & ~has)[1])
    out['jackknife_tiles'] = dict(n_tiles=len(tj), min_snr=min(v['snr'] for v in tj.values()),
                                  max_snr=max(v['snr'] for v in tj.values()),
                                  worst=sorted(tj.items(), key=lambda kv: kv[1]['snr'])[:5])
    # dominant tile-petal jackknife
    pj = []
    for g in np.unique(grp[ok]):
        pj.append(mf(O, Es, Z, ok & (grp != g))[1])
    out['jackknife_tilepetal'] = dict(n=len(pj), min_snr=float(min(pj)), max_snr=float(max(pj)))
    np.savez_compressed(D / 'verify_desi_detect_nulls.npz', rot=rot, rotR=np.array(rr), shuffle=sh,
                        fake=np.array(list(fl.values())), fake_names=np.array(list(fl.keys())))
    (D / 'verify_desi_detect.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ('n', 'err_scale_from_nulls', 'true', 'shuffle_within_petal',
                                          'rigid_rotations')}, indent=1))
    print('fake', {k: round(v, 2) for k, v in fl.items()})
    print('translations', {k: round(v['snr'], 1) for k, v in tr.items()})
    print('splits', {k: (v['n'], round(v['snr'], 1)) for k, v in spl.items()})
    print('jackknife min snr', min(v['snr'] for v in jk.values()))
    print('tile jackknife', out['jackknife_tiles']['min_snr'], 'tile-petal', out['jackknife_tilepetal'])


if __name__ == '__main__':
    main()
