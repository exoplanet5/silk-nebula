"""verify_precovery step 3b: four-band (u, g, r, i) frame-level test of run 4458 and robustness variants.

Uses process_frame / backfit from verify_precovery_frames.py.  Blocks for u and r are added to a separate cache
(data/derived/verify_precovery_frames_blocks_ur.npz).  For every band the ridge amplitude is fitted with the unit
E+W template (FWHM 0.28 deg) and reported in nmgy arcsec^-2 so that the bands can be compared as an SED:
  [OIII]+[OII] line spectrum (DESI ratios, south)  -> u/g ~ 1.1-1.4, r/g ~ 0.37, i/g ~ 0.10
  flat f_nu continuum                              -> 1 : 1 : 1 : 1 (redder for dust-scattered starlight)
Variants: per-camcol drift order 2 and 4; asymmetric |d| windows (west-heavy / east-heavy) to break the symmetry of
the selection about the spine; Dec halves. Nulls: template shifted by +-0.6..1.2 deg with a full re-fit.
Output: data/derived/verify_precovery_frames_sed.json
"""
import sys, glob, json
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from verify_precovery_frames import process_frame, backfit, RUN, B
from verify_precovery_lib import MYRAW, DER, load_strands, make_model, band_response_per_R5007, nmgy_per_R

CACHE = DER / 'verify_precovery_frames_blocks_ur.npz'


def blocks_ur():
    if CACHE.exists():
        return dict(np.load(CACHE))
    out = {k: [] for k in ('band', 'camcol', 'field', 'ix', 'iy', 'tot', 'skymodel', 'img', 'fgood', 'ra', 'dec')}
    files = sorted(glob.glob(str(MYRAW / f'frames_{RUN}/frame-[ur]-{RUN:06d}-*.fits.bz2')))
    for k, p in enumerate(files):
        name = p.split('/')[-1]
        b = name.split('-')[1]; c = int(name.split('-')[3]); f = int(name.split('-')[4][:4])
        try:
            d = process_frame(p)
        except Exception as e:
            print('skip', name, e); continue
        n = len(d['ra'])
        out['band'] += [b] * n; out['camcol'] += [c] * n; out['field'] += [f] * n
        for kk in ('ix', 'iy', 'tot', 'skymodel', 'img', 'fgood', 'ra', 'dec'):
            out[kk] += list(d[kk])
        if k % 50 == 0:
            print(k, len(files), name, flush=True)
    arr = {k: np.array(v) for k, v in out.items()}
    np.savez_compressed(CACHE, **arr)
    return arr


def main():
    gi = dict(np.load(DER / 'verify_precovery_frames_blocks.npz'))
    ur = blocks_ur()
    arr = {k: np.concatenate([gi[k], ur[k]]) for k in gi}
    strands, _ = load_strands()
    model = make_model('gauss028', strands)
    res = {'run': RUN, 'units': 'A in 1e-3 nmgy arcsec^-2 (peak of unit E+W ridges, FWHM 0.28 deg)', 'bands': {}}
    unit_R = {'u': nmgy_per_R('u', 3728.48), 'g': nmgy_per_R('g', 5008.24) + nmgy_per_R('g', 4960.30) * 4960.30 / 5008.24 / 2.98,
              'r': nmgy_per_R('r', 6564.61), 'i': band_response_per_R5007('i')[0]}
    for band in ('u', 'g', 'r', 'i'):
        m = (arr['band'] == band) & np.isfinite(arr['tot'])
        if m.sum() < 1000:
            continue
        cam = arr['camcol'][m]; fld = arr['field'][m]; iy = arr['iy'][m]; ix = arr['ix'][m]
        tau = ((fld * 1361 + iy * B) // B) * B; col = cam * 100 + ix
        ra, dec = arr['ra'][m], arr['dec'][m]; y = arr['tot'][m]
        T = model(ra, dec)
        dE, _, _, _ = strands['E'].geometry(ra, dec)
        out = {}
        variants = {'base': (np.ones(m.sum(), bool), 2), 'drift4': (np.ones(m.sum(), bool), 4),
                    'west_heavy': (dE < 0.9, 2), 'east_heavy': (dE > -0.9, 2),
                    'dec_lt_17': (dec < 17.0, 2), 'dec_ge_17': (dec >= 17.0, 2),
                    'no_star_18p5': (np.hypot((ra - 116.531) * np.cos(np.radians(18.5)), dec - 18.510) > 0.5, 2)}
        for name, (sel, drift) in variants.items():
            f = backfit(y[sel], tau[sel], col[sel], cam[sel], T[sel], drift=drift)
            out[name] = dict(A=float(f['A'] * 1e3), n=int(f['mask'].sum()))
        nulls = []
        for s in (-1.2, -1.0, -0.8, -0.6, 0.6, 0.8, 1.0, 1.2):
            Ts = make_model('gauss028', strands, ra_shift=s)(ra, dec)
            # joint: real template removed with its fitted amplitude, then the shifted one fitted with a full backfit
            f = backfit(y - out['base']['A'] * 1e-3 * T, tau, col, cam, Ts)
            nulls.append(float(f['A'] * 1e3))
        out['nulls'] = nulls; out['null_std'] = float(np.std(nulls)); out['null_mean'] = float(np.mean(nulls))
        out['unit_nmgy_per_R'] = unit_R[band]
        out['A_R_lineunits'] = out['base']['A'] * 1e-3 / unit_R[band]
        res['bands'][band] = out
        print(band, {k: (round(v['A'], 2) if isinstance(v, dict) else v) for k, v in out.items() if k not in ('nulls',)}, flush=True)
    g = res['bands']['g']['base']['A']
    res['sed_over_g'] = {b: dict(ratio=v['base']['A'] / g, err=v['null_std'] / g) for b, v in res['bands'].items()}
    res['expected_over_g'] = {'line_spectrum_DESI': {b: band_response_per_R5007(b)[0] / band_response_per_R5007('g')[0] for b in 'ugri'},
                              'flat_fnu': {b: 1.0 for b in 'ugri'}}
    print(json.dumps(res['sed_over_g'], indent=1)); print(json.dumps(res['expected_over_g'], indent=1))
    (DER / 'verify_precovery_frames_sed.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
