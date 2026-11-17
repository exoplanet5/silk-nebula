"""w3_oii_desi_regress.py -- (wave 3, tag oii) DESI DR1 pixel-level IV regression spectra for the [O II] doublet test.

For every usable coadd (wave-2 'base' set: clean [O III] fit, finite geometry instrument, tile-petal contributions)
the sky-subtraction residual r = flux - Redrock model is continuum-subtracted (quadratic, all nebular/sky lines and the
object's own redshifted lines masked, 2 passes with 4-sigma clipping) in five windows:
    3727 (3690-3770), 5007 (4930-5045), SII (6690-6770), 7325 (7282-7375), Ha (6520-6610).
Instruments:
  * 'geo'  : Zf = Z / kR_5008, the geometry-only petal-subtracted template (wave-2 verify_desi_coadd_meta.npz).
             Denominator = per-fibre [O III]5007 Gaussian amplitude (wave-2 fixed-velocity fit, o3).
             S(lambda) = flux density per unit [O III]5007 flux of the filament-correlated component.
  * 'ha'   : mean Halpha amplitude of OTHER-tile-petal fibres within 0.15 deg (neighbour-Halpha instrument), with the
             fibre's own Halpha amplitude as denominator -> spectrum of whatever Halpha-bright structure exists
             (WIM / Halpha patch), per unit Halpha.
Two weightings: 'P' per-pixel ivar (as wave 2); 'F' one weight per fibre (median ivar over the window core),
so every pixel uses the same fibre mix. Chromatic acceptance f(lambda)/f(5008) applied as in wave 2.
Errors: HEALPix nside=32 block bootstrap, 400 resamples, same resamples for all windows of a subset.

Also writes per-fibre LSF statistics from the SPARCL wave_sigma arrays of the pilot download.
Robust variant ('robust' argument): pixels with |r| sqrt(ivar) > 5 after continuum subtraction get zero weight, and
fibres whose off-line chi^2/pixel in the window exceeds 3x the median are dropped; geometry instrument only,
windows 3727/5007/SII/7325 -> data/derived/w3_oii_desi_regspec_robust.npz.
Output: data/derived/w3_oii_desi_regspec.npz, data/derived/w3_oii_desi_lsf.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_oii_desi_regress.py [robust]
"""
import json
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import (D, WAVE, NEB_LINES, SKY_LINES, GAL_REST, load_spectra, load_meta,  # noqa: E402
                                good_fit)

RNG = np.random.default_rng(20261117)
NB = 400
WINS = {'3727': (3690.0, 3770.0, '3727', (3722.0, 3735.0)),
        '5007': (4930.0, 5045.0, '5008', (4955.0, 5013.0)),
        'SII': (6690.0, 6770.0, '6724', (6714.0, 6737.0)),
        '7325': (7282.0, 7375.0, '7325', (7316.0, 7336.0)),
        'Ha': (6520.0, 6610.0, '6563', (6560.0, 6588.0))}


def contsub(R, W, z, lo, hi):
    s = (WAVE >= lo) & (WAVE <= hi)
    lam = WAVE[s]
    r = R[:, s].astype(float)
    w = W[:, s].astype(float)
    bad = ~np.isfinite(r) | ~np.isfinite(w) | (w <= 0)
    r[bad], w[bad] = 0.0, 0.0
    m = np.ones(lam.size, bool)
    for l0 in NEB_LINES:
        m &= np.abs(lam - l0) > 3.0
    for l0 in SKY_LINES:
        m &= np.abs(lam - l0) > 2.5
    own = np.zeros(r.shape, bool)
    for gl in GAL_REST:
        own |= np.abs(lam[None, :] - gl * (1 + z)[:, None]) < 8.0
    w[own] = 0.0
    x = (lam - lam.mean()) / (hi - lo)
    B = np.stack([np.ones_like(x), x, x * x])
    wc = w * m[None, :]
    for it in range(2):
        A = np.einsum('kp,np,jp->nkj', B, wc, B)
        b = np.einsum('kp,np,np->nk', B, wc, r)
        ok = np.linalg.cond(A) < 1e12
        A[~ok] = np.eye(3)
        p = np.linalg.solve(A, b[..., None])[..., 0]
        cont = p @ B
        if it == 0:
            zz = (r - cont) * np.sqrt(wc)
            wc = wc * (np.abs(zz) < 4)
    w[~ok] = 0.0
    return lam, r - cont, w


def regspec_boot(c, w, Zi, Y, m, inv, nblk, CNT):
    """IV slope per pixel: sum u dz (c - cb) / sum u dz (Y - yb). c, w: [N, P]; Zi, Y: [N]. CNT: [NB, nblk] counts.
    Returns S [P] and bootstrap [NB, P]. Uses moment expansion so all resamples are matrix products."""
    c, w, Zi, Y = c[m], w[m], Zi[m], Y[m]
    inv = inv[m]
    Wt = np.vstack([np.ones(nblk), CNT]).astype(float)[:, inv]          # [1+NB, N] fibre multiplicities
    wu = w
    s0 = Wt @ wu
    sz = Wt @ (wu * Zi[:, None])
    sy = Wt @ (wu * Y[:, None])
    sc = Wt @ (wu * c)
    szc = Wt @ (wu * c * Zi[:, None])
    szy = Wt @ (wu * (Zi * Y)[:, None])
    num = szc - sz * sc / s0
    den = szy - sz * sy / s0
    S = num / den
    return S[0], S[1:]


def robustify(lam, c, w, lines_core):
    """Zero-weight 5-sigma pixel outliers and fibres with off-line chi2/pix > 3x median (robust variant)."""
    z = c * np.sqrt(w)
    w = np.where(np.abs(z) > 5, 0.0, w)
    off = np.ones(lam.size, bool)
    for lo, hi in lines_core:
        off &= ~((lam > lo) & (lam < hi))
    npx = np.maximum((w[:, off] > 0).sum(1), 1)
    chi = np.sum(np.where(w[:, off] > 0, z[:, off] ** 2, 0), 1) / npx
    bad = chi > 3 * np.median(chi[np.isfinite(chi) & (chi > 0)])
    w[bad] = 0.0
    return w, int(bad.sum())


def main():
    ROBUST = len(sys.argv) > 1 and sys.argv[1] == 'robust'
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    ids_lf = lf['sparcl_id']
    ids, F, I, M, src = load_spectra()
    pos = {s: k for k, s in enumerate(ids)}
    k = np.array([pos[s] for s in ids_lf])
    R = (F - M)[k]
    I = I[k]
    del F, M
    keep, meta = load_meta(ids_lf)
    assert keep.all()
    import healpy as hp
    blocks = hp.ang2pix(32, meta['ra'], meta['dec'], lonlat=True)
    O = lf['o3__amp'][:, 0]
    Zf = cm['Z'] / cm['kR_5008']
    base = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & np.isfinite(Zf) & cm['ok_contrib']
    dec, ra = meta['dec'], meta['ra']
    hiz = (meta['redshift'] > 0.5) | (meta['spectype'] == 'QSO')
    # neighbour-Halpha instrument (other tile-petals within 0.15 deg)
    Ha = lf['han2__amp'][:, 0]
    eHa = lf['han2__err'][:, 0]
    okHa = good_fit(lf, 'han2', 0) & cm['ok_contrib']
    grp_lab = np.array([f'{t}_{p}' for t, p in zip(cm['group_tile'], cm['group_petal'])])
    from scipy.spatial import cKDTree
    xyz = np.c_[np.cos(np.radians(dec)) * np.cos(np.radians(ra)), np.cos(np.radians(dec)) * np.sin(np.radians(ra)),
                np.sin(np.radians(dec))]
    bi = np.where(okHa)[0]
    tree = cKDTree(xyz[bi])
    wHa = 1 / eHa ** 2
    nbrHa = np.full(len(O), np.nan)
    for i, nb in enumerate(tree.query_ball_point(xyz, 2 * np.sin(np.radians(0.15) / 2))):
        js = bi[nb]
        js = js[grp_lab[js] != grp_lab[i]]
        if js.size >= 3:
            nbrHa[i] = np.sum(wHa[js] * Ha[js]) / np.sum(wHa[js])
    patch = (ra > 115.0) & (ra < 116.2) & (dec > 19.6) & (dec < 21.6)
    baseHa = okHa & np.isfinite(nbrHa) & np.isfinite(Ha)
    subsets = {
        'geo': {
            'all': base, 'south': base & (dec < 19.5), 'north': base & (dec >= 19.5),
            'Zpos': base & (Zf > 0), 'Zneg': base & (Zf < 0),
            'south_Zpos': base & (dec < 19.5) & (Zf > 0), 'south_Zneg': base & (dec < 19.5) & (Zf < 0),
            'north_Zpos': base & (dec >= 19.5) & (Zf > 0), 'north_Zneg': base & (dec >= 19.5) & (Zf < 0),
            'qso_hiz': base & hiz, 'lowz_gal': base & ~hiz,
            'dark': base & (meta['program'] == 'dark'), 'bright': base & (meta['program'] == 'bright'),
            'south_qso_hiz': base & hiz & (dec < 19.5), 'no_patch': base & ~patch,
            'dec12_15p5': base & (dec < 15.5), 'dec15p5_17p5': base & (dec >= 15.5) & (dec < 17.5),
            'dec17p5_19p5': base & (dec >= 17.5) & (dec < 19.5), 'dec19p5_22': base & (dec >= 19.5) & (dec < 22),
            'dec22_25': base & (dec >= 22) & (dec < 25), 'dec25_29': base & (dec >= 25),
            'w1_only': base & (lf['src'] == 'w1'), 'ext_only': base & (lf['src'] == 'ext'),
        },
        'ha': {'all': baseHa, 'patch': baseHa & patch, 'no_patch': baseHa & ~patch,
               'south': baseHa & (dec < 19.5), 'north': baseHa & (dec >= 19.5)},
    }
    out = {'v_neb_o3_fit': float(lf['v_neb']), 'n_base': int(base.sum()), 'n_baseHa': int(baseHa.sum())}
    cont = {}
    for wn, (lo, hi, ck, core) in WINS.items():
        lam, c, w = contsub(R, I, meta['redshift'], lo, hi)
        chrom = (cm[f'f_chrom_{ck}'] / cm['f_chrom_5008'])
        c = c * chrom[:, None]
        w = w / chrom[:, None] ** 2
        if ROBUST:
            w, nbad = robustify(lam, c, w, [WINS[wn][3]])
            out[f'n_dropped_{wn}'] = nbad
        cont[wn] = (lam, c, w)
        out[f'lam_{wn}'] = lam
        print(wn, 'contsub done', flush=True)
    if ROBUST:
        keep_s = ('all', 'south', 'north', 'Zpos', 'Zneg', 'qso_hiz', 'lowz_gal', 'dark', 'bright', 'w1_only',
                  'ext_only', 'south_qso_hiz')
        subsets = {'geo': {k: v for k, v in subsets['geo'].items() if k in keep_s}}
        cont.pop('Ha', None)
    for inst, subs in subsets.items():
        for sn, m in subs.items():
            if m.sum() < 80:
                print('skip', inst, sn, m.sum())
                continue
            ub, inv = np.unique(blocks, return_inverse=True)
            nblk = ub.size
            # resample only blocks present in the subset; others get zero weight anyway
            present = np.unique(inv[m])
            CNT = np.zeros((NB, nblk), np.int32)
            for b in range(NB):
                draw = RNG.choice(present, present.size, replace=True)
                CNT[b] = np.bincount(draw, minlength=nblk)
            if inst == 'geo':
                Zi, Y = Zf, O
            else:
                Zi, Y = nbrHa, Ha
            out[f'n_{inst}_{sn}'] = int(m.sum())
            out[f'nblk_{inst}_{sn}'] = int(present.size)
            for wn, (lam, c, w) in cont.items():
                if inst == 'ha' and wn in ('5007', '7325'):
                    continue
                for wt in ('P', 'F'):
                    if wt == 'P':
                        ww = w
                    else:
                        lo_c, hi_c = WINS[wn][3]
                        cs = (lam >= lo_c) & (lam <= hi_c)
                        wcore = np.where(w[:, cs] > 0, w[:, cs], np.nan)
                        wf = np.nan_to_num(np.nanmedian(wcore, 1), nan=0.0)
                        ww = (w > 0) * wf[:, None]
                    S, bs = regspec_boot(c, ww, Zi, Y, m, inv, nblk, CNT)
                    out[f'S_{inst}_{sn}_{wn}_{wt}'] = S
                    out[f'bs_{inst}_{sn}_{wn}_{wt}'] = bs.astype(np.float32)
            print(inst, sn, int(m.sum()), 'blocks', present.size, flush=True)
    np.savez_compressed(D / ('w3_oii_desi_regspec_robust.npz' if ROBUST else 'w3_oii_desi_regspec.npz'), **out)
    if ROBUST:
        return
    # ---------------------------------------------------------------- LSF statistics (pilot wave_sigma)
    S_ = []
    for f in sorted((D / 'spectra_desi_windows').glob('pilot_*.npz')):
        with np.load(f, allow_pickle=False) as d:
            if 'wave_sigma' in d.files:
                S_.append(d['wave_sigma'])
    S_ = np.concatenate(S_).astype(float)
    S_[S_ <= 0] = np.nan
    lsf = {'n_spectra': int(S_.shape[0]), 'source': 'SPARCL DESI-DR1 wave_sigma (camcoadd BRZ wavedisp), pilot download'}
    for lam0 in (3727.1, 3729.9, 4960.3, 5008.2, 6718.3, 6732.7, 7321.0, 7331.4):
        v = S_[:, int(np.argmin(np.abs(WAVE - lam0)))]
        q = np.nanpercentile(v, [2.5, 16, 50, 84, 97.5])
        lsf[f'{lam0:.1f}'] = dict(median=float(q[2]), p2p5=float(q[0]), p16=float(q[1]), p84=float(q[3]),
                                  p97p5=float(q[4]), mean=float(np.nanmean(v)), std=float(np.nanstd(v)))
    # within-fibre difference between the two [O II] lines and between 3727 and 5008
    i1, i2, i5 = (int(np.argmin(np.abs(WAVE - x))) for x in (3727.1, 3729.9, 5008.2))
    lsf['ratio_3729_over_3726_per_fibre'] = dict(median=float(np.nanmedian(S_[:, i2] / S_[:, i1])),
                                                 std=float(np.nanstd(S_[:, i2] / S_[:, i1])))
    lsf['ratio_3727_over_5008_per_fibre'] = dict(median=float(np.nanmedian(S_[:, i1] / S_[:, i5])),
                                                 std=float(np.nanstd(S_[:, i1] / S_[:, i5])))
    (D / 'w3_oii_desi_lsf.json').write_text(json.dumps(lsf, indent=1))
    print(json.dumps(lsf, indent=1))


if __name__ == '__main__':
    main()
