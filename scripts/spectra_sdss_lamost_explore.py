"""Quick look: mean sky-fibre spectra (flux+sky) of a few SDSS plates around [OIII], Hbeta, Halpha.
usage: python -I spectra_sdss_lamost_explore.py PLATE [PLATE ...]
"""
import sys
import glob
import numpy as np
from astropy.io import fits
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_sdss_lamost_lines import VAC, sun_at, gauss, linfit, flux_to_rayleigh, C_KMS
from silk_style import new_sky_figure, finish, label, OIII_C, HB_C, HA_C, FGDIM

ROOT = '/Users/mickey/astronomy/silk-nebula'


def mean_sky(plate):
    fs = sorted(glob.glob(f'{ROOT}/data/raw/sdss_dr17_speclite/{plate:04d}/spec-*.fits'))
    out = {}
    for f in fs:
        h = fits.open(f)
        d = h[1].data
        mjd = h[0].header['MJD']
        sp = 1 if h[2].data['FIBERID'][0] <= (320 if len(d) < 4200 else 500) else 2
        lam = 10 ** d['loglam']
        out.setdefault((mjd, sp), []).append((lam, d['flux'] + d['sky'], d['ivar'], d['sky'], d['wdisp']))
    res = {}
    for k, v in out.items():
        lam0 = v[0][0]
        tot = np.array([np.interp(lam0, x[0], x[1]) for x in v])
        sky = np.array([np.interp(lam0, x[0], x[3]) for x in v])
        wd = np.median([np.interp(lam0, x[0], x[4]) for x in v], axis=0)
        res[k] = dict(lam=lam0, tot=np.nanmean(tot, 0), totsd=np.nanstd(tot, 0) / np.sqrt(len(v)),
                      sky=np.median(sky, 0), n=len(v), wdisp=wd)
    return res


if __name__ == '__main__':
    plates = [int(p) for p in sys.argv[1:]]
    fig, axes = new_sky_figure(13, 3.2 * len(plates), 'SDSS sky-fibre spectra · quick look', 'mean total sky of SKY fibres per spectrograph (DR17 spec-lite); Hβ, [O III], Hα+[N II] windows; dotted = rest λ',
                               nrows=len(plates), ncols=3, squeeze=False)
    for i, p in enumerate(plates):
        r = mean_sky(p)
        for (mjd, sp), s in sorted(r.items()):
            for j, (lo, hi, lines) in enumerate([(4840, 4885, ['Hb']), (4940, 5030, ['OIII4959', 'OIII5007']), (6530, 6600, ['NII6548', 'Ha', 'NII6583'])]):
                m = (s['lam'] > lo) & (s['lam'] < hi)
                ax = axes[i, j]
                ax.plot(s['lam'][m], s['tot'][m], lw=0.9, label=f'{p}-{mjd} sp{sp} n={s["n"]}')
                for ln in lines:
                    ax.axvline(VAC[ln], color=FGDIM, lw=0.5, ls=':')
                if j == 0:
                    ax.legend(fontsize=6)
                ax.set_xlabel('vacuum wavelength, heliocentric (Å)')
                if j == 0:
                    ax.set_ylabel('flux+sky (1e-17 erg/s/cm²/Å/fibre)', fontsize=7.5)
            print(p, mjd, sp, s['n'], 'median wdisp', np.median(s['wdisp']))
    finish(fig, f'{ROOT}/figures/spectra_sdss_lamost_quicklook.png')
