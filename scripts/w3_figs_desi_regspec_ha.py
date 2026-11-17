"""DESI DR1 pixel-level IV regression spectrum around Ha + [N II] (and [S II]) with the wave-2 method.

verify_desi_regspec.npz (wave 2) holds 3727, 3869, 4363, 4686, Hb, 5007, 7150, 7300 windows but no Ha window.
This script calls the unchanged wave-2 functions verify_desi_regspec.contsub / regspec (continuum-subtracted
residuals x chromatic acceptance, IV slope on the geometry template, HEALPix nside-32 block bootstrap) for the
6525-6605 A and 6700-6750 A windows, subsets 'all' and 'south' (Dec < 19.5) exactly as in verify_desi_regspec.main.
S(lambda) is flux density per unit [O III]5007 flux of the filament component (energy, per A).
Output: data/derived/w3_figs_desi_regspec_ha.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_desi_regspec_ha.py
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_desi_regspec as VR
from verify_desi_common import D, load_spectra, load_meta, good_fit

WINS = {'Ha': (6525.0, 6605.0, '6563'), 'SII': (6700.0, 6752.0, '6724')}


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    cm = np.load(D / 'verify_desi_coadd_meta.npz')
    ids_lf = lf['sparcl_id']
    ids, F, I, M, src = load_spectra()
    pos = {s: k for k, s in enumerate(ids)}
    k = np.array([pos[s] for s in ids_lf])
    R = (F - M)[k]
    I = I[k]
    del F, M
    keep, meta = load_meta(ids_lf)
    import healpy as hp
    blocks = hp.ang2pix(32, meta['ra'], meta['dec'], lonlat=True)
    det = json.loads((D / 'verify_desi_detect.json').read_text())
    O = lf['o3__amp'][:, 0]
    Zf = cm['Z'] / cm['kR_5008']
    base = good_fit(lf, 'o3', scale=det['err_scale_from_nulls']) & np.isfinite(Zf) & cm['ok_contrib']
    subsets = {'all': base, 'south': base & (meta['dec'] < 19.5)}
    out = {}
    for wn, (lo, hi, ch) in WINS.items():
        lam, c, w = VR.contsub(R, I, meta['redshift'], lo, hi)
        c = c * (cm[f'f_chrom_{ch}'] / cm['f_chrom_5008'])[:, None]
        w = w / (cm[f'f_chrom_{ch}'] / cm['f_chrom_5008'])[:, None] ** 2
        out[f'lam_{wn}'] = lam
        for sn, m in subsets.items():
            S, bs = VR.regspec(c, w, Zf, O, m, blocks)
            out[f'S_{wn}_{sn}'] = S
            out[f'bs_{wn}_{sn}'] = bs.astype(np.float32)
            print(wn, sn, 'n=%d' % m.sum(), 'max %.3f' % np.nanmax(S), flush=True)
    out['note'] = np.array('wave-2 method (verify_desi_regspec.contsub/regspec), run by w3_figs; per unit [O III]5007 flux')
    np.savez_compressed(D / 'w3_figs_desi_regspec_ha.npz', **out)
    print('saved')


if __name__ == '__main__':
    main()
