"""verify_desi_setup.py -- per-coadd metadata for the wave-2 DESI verification (tag desi).

For every spectrum in data/derived/verify_desi_linefluxes.npz computes and saves
  * true point-source fibre acceptance f_true(lambda) from the GFA (EFFTIME-weighted harmonic mean over exposures),
    chromatic (Moffat beta=3.5, FWHM ~ lambda^-0.2 for the whole FWHM) and achromatic (value at the GFA band);
  * conversion k_R(lambda): Rayleigh of uniform surface brightness per 1e-17 erg/s/cm^2 of calibrated flux,
    k_R = f_true / (Omega_fibre * 1 R in erg/s/cm^2/arcsec^2), Omega from the 1.52" mean fibre (DESI-6043);
  * wave-1 equivalent (f = MEAN_PSF_TO_FIBER_SPECFLUX) for comparison;
  * EFFTIME-weighted barycentric correction, tile list, dominant tile-petal, template prediction Z (geometry only).
Output: data/derived/verify_desi_coadd_meta.npz
"""
import sys
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_desi_common import D, load_meta, Exposures, Coadds, Template, petal_subtracted, L  # noqa: E402

HC = 6.62607015e-27 * 2.99792458e18     # erg A
ARCSEC2_SR = (np.pi / 648000) ** 2
OMEGA = np.pi * (1.52 / 2) ** 2         # arcsec^2


def one_R(lam):
    """1 Rayleigh in erg s^-1 cm^-2 arcsec^-2 at vacuum wavelength lam (A)."""
    return 1e6 / (4 * np.pi) * HC / lam * ARCSEC2_SR


LAMS = {'3727': L['OII3726'], '3869': L['NeIII3869'], '4363': L['OIII4363'], '4686': L['HeII4686'],
        '4861': L['Hb'], '5008': L['OIII5007'], '6563': L['Ha'], '6583': L['NII6583'], '6724': L['SII6716'],
        '7136': L['ArIII7136'], '7325': L['OII7320'], '9069': L['SIII9069']}


def main():
    lf = np.load(D / 'verify_desi_linefluxes.npz')
    ids = lf['sparcl_id']
    keep, meta = load_meta(ids)
    assert keep.all()
    ex = Exposures()
    co = Coadds(meta, ex)
    out = dict(sparcl_id=ids, **{k: v for k, v in meta.items()})
    for nm, lam in LAMS.items():
        out[f'f_chrom_{nm}'] = co.f_true(lam, atm_frac=1.0)
        out[f'f_half_{nm}'] = co.f_true(lam, atm_frac=0.5)
        out[f'kR_{nm}'] = out[f'f_chrom_{nm}'] * 1e-17 / OMEGA / one_R(lam)
    out['f_achrom'] = co.f_true(6400.0, atm_frac=1.0)
    out['f_achrom_lamref6100'] = co.f_true(6100.0, lam_ref=6100.0)
    out['f_chrom_5008_lamref6100'] = co.f_true(L['OIII5007'], lam_ref=6100.0)
    out['f_chrom_5008_lamref6700'] = co.f_true(L['OIII5007'], lam_ref=6700.0)
    out['kR_w1_5008'] = np.clip(meta['mean_psf_to_fiber_specflux'], 0.3, 1.0) * 1e-17 / OMEGA / one_R(L['OIII5007'])
    out['one_R_5008'] = one_R(L['OIII5007'])
    vc, vsp = co.vcorr()
    out['vcorr'], out['vspread'] = vc, vsp
    out['n_tiles'] = co.n_tiles
    out['ok_contrib'] = co.ok
    out['group_tile'] = co.group[:, 0]
    out['group_petal'] = co.group[:, 1]
    out['tiles'] = np.array([','.join(map(str, co.tiles_of(i))) for i in range(len(ids))])
    out['has_41485_41487'] = np.array([any(t in (41485, 41487) for t in co.tiles_of(i)) for i in range(len(ids))])
    # GFA conditions of the coadd
    see, effs, ffs, filled = [], [], [], []
    for i in range(len(ids)):
        js, ws = co.exposures_of(i)
        if js.size:
            see.append(np.sum(ws * ex.see[js]) / ws.sum())
            effs.append(ws.sum())
            filled.append(np.mean(ex.ff_filled[js]))
        else:
            see.append(np.nan)
            effs.append(np.nan)
            filled.append(np.nan)
    out['seeing_gfa'], out['efftime'], out['frac_ff_filled'] = np.array(see), np.array(effs), np.array(filled)
    T = Template()
    Zp = petal_subtracted(T, meta, co, parts=True)
    for s, (z, t, sub) in Zp.items():
        out[f'Z_{s}'], out[f'T_{s}'], out[f'sub_{s}'] = z, t, sub
    out['Z'] = sum(Zp[s][0] for s in Zp)
    out['T'] = sum(Zp[s][1] for s in Zp)
    out['ff_scale_vs_moffat'] = ex.ff_scale_vs_moffat
    np.savez_compressed(D / 'verify_desi_coadd_meta.npz', **out)
    r = out['f_chrom_5008'] / np.clip(meta['mean_psf_to_fiber_specflux'], 0.3, 1.0)
    ok = np.isfinite(r)
    print('n', len(ids), 'ok contrib', co.ok.sum(), 'ff_scale_vs_moffat', round(ex.ff_scale_vs_moffat, 3))
    print('f_true(5008) median %.3f p16 %.3f p84 %.3f' % tuple(np.nanpercentile(out['f_chrom_5008'], [50, 16, 84])))
    print('wave-1 ffrac median %.3f' % np.median(meta['mean_psf_to_fiber_specflux']))
    print('ratio f_true/ffrac median %.3f p16 %.3f p84 %.3f' % tuple(np.percentile(r[ok], [50, 16, 84])))
    for p in ('dark', 'bright', 'backup'):
        m = ok & (meta['program'] == p)
        if m.any():
            print(p, m.sum(), 'ratio median %.3f' % np.median(r[m]))
    for nm in LAMS:
        q = out[f'f_chrom_{nm}'] / out['f_chrom_5008']
        q2 = out[f'f_half_{nm}'] / out['f_half_5008']
        print(f'f({nm})/f(5008): full-atm {np.nanmedian(q):.3f}  half-atm {np.nanmedian(q2):.3f}')


if __name__ == '__main__':
    main()
