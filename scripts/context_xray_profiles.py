"""Where does the Silk filament lie relative to the Monogem soft-X-ray emission?  ROSAT diffuse-background maps.

Input : data/raw/context_skyview/rass_bkg{1,2,4,5,6}_l205b12_car.fits  (SkyView 'RASS Background N' = Snowden+1997
        ROSAT All-Sky Survey diffuse maps R1..R7, Galactic CAR 60x50 deg, 0.2 deg px, centre (205, 12.5);
        fetched 2026-11-16 from https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl)
Output: data/derived/context_rosat_profiles.csv   spine-parallel profiles (R12=R1+R2, R45=R4+R5) in 4 Dec bands
        data/derived/context_rosat_spine.csv      inner (W, -5..-1 deg) / outer (E, +1..+5 deg) means and ratios
        figures/context_rosat_profiles.png
Method: maps smoothed with a Gaussian of FWHM 1.0 deg; for each Dec band the map is averaged along curves parallel to
        the spine (offset dRA*cos(Dec) = const, 0.25-deg Dec steps); errors = std of the samples / sqrt(N_indep),
        N_indep = band length / 1 deg.
"""
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, map_coordinates
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *
from silk_style import *

RAW = ROOT / 'data/raw/context_skyview'
BANDS = [(12.0, 16.0), (16.0, 20.0), (20.0, 24.0), (24.0, 28.25)]


def load(band):
    h = fits.open(RAW / f'rass_bkg{band}_l205b12_car.fits')[0]
    return np.nan_to_num(h.data.astype(float)), WCS(h.header)


def sampler(img, w, fwhm_deg=1.0):
    sm = gaussian_filter(img, fwhm_deg / 2.355 / abs(w.wcs.cdelt[0]))

    def f(ra, dec):
        l, b = gal(ra, dec)
        x, y = w.wcs_world2pix(l, b, 0)
        return map_coordinates(sm, [y, x], order=1, mode='nearest')
    return f


def main():
    imgs = {b: load(b) for b in (1, 2, 4, 5)}
    w = imgs[4][1]
    S = {'R12': sampler(imgs[1][0] + imgs[2][0], w), 'R45': sampler(imgs[4][0] + imgs[5][0], w)}
    tr = load_trace(); sp = tr['spine']
    sra = np.array([q['ra'] for q in sp]); sde = np.array([q['dec'] for q in sp])
    dras = np.arange(-15, 8.01, 0.5)
    rows, summ = [], []
    for d0, d1 in BANDS:
        decs = np.arange(d0, d1 + 1e-6, 0.25)
        ras_s = np.interp(decs, sde, sra)
        nind = max((d1 - d0) / 1.0, 1)
        for key, f in S.items():
            prof, err = [], []
            for a in dras:
                v = f(ras_s + a / np.cos(np.radians(decs)), decs)
                prof.append(v.mean()); err.append(v.std() / np.sqrt(nind))
                rows.append(dict(band=f'{d0:g}-{d1:g}', map=key, dra_deg=a, mean=v.mean(), err=v.std() / np.sqrt(nind)))
            prof, err = np.array(prof), np.array(err)
            win_in, win_out = (dras >= -5) & (dras <= -1), (dras >= 1) & (dras <= 5)
            mi, mo = prof[win_in].mean(), prof[win_out].mean()
            ei, eo = err[win_in].mean() / np.sqrt(2), err[win_out].mean() / np.sqrt(2)
            # peak west of the spine and the half-level crossing east of it (X-ray 'edge')
            wwin = (dras >= -8) & (dras <= 0)
            ipk = np.where(wwin)[0][np.argmax(prof[wwin])]
            far_e = prof[(dras >= 4) & (dras <= 8)].mean()
            half = 0.5 * (prof[ipk] + far_e)
            east = np.arange(ipk, len(dras))
            cross = east[np.where(prof[east] < half)[0]]
            edge = np.nan
            if len(cross):
                k = cross[0]
                edge = float(np.interp(half, [prof[k], prof[k - 1]], [dras[k], dras[k - 1]]))
            summ.append(dict(band=f'{d0:g}-{d1:g}', map=key, inner_W=mi, e_inner=ei, outer_E=mo, e_outer=eo,
                             ratio_W_over_E=mi / mo, e_ratio=(mi / mo) * np.hypot(ei / mi, eo / mo),
                             at_spine=float(prof[dras == 0][0]),
                             far_W_m12_m8=prof[(dras >= -12) & (dras <= -8)].mean(),
                             peak_W=float(prof[ipk]), peak_offset_deg=float(dras[ipk]), far_E_4_8=far_e,
                             peak_over_farE=float(prof[ipk] / far_e), half_level_edge_offset_deg=edge))
    pd.DataFrame(rows).to_csv(ROOT / 'data/derived/context_rosat_profiles.csv', index=False, float_format='%.4g')
    sm = pd.DataFrame(summ)
    sm.to_csv(ROOT / 'data/derived/context_rosat_spine.csv', index=False, float_format='%.4g')
    print(sm.to_string())

    prof = pd.DataFrame(rows)
    fig, axes = new_sky_figure(11, 5.6, 'Silk vs Monogem soft X-rays',
                               'ROSAT RASS diffuse maps, profiles parallel to the [O III] spine (FWHM 1 deg)', ncols=2)
    cols = [ACCENT, OK, WARN, DANGER]
    for ax, key, ttl in zip(axes, ['R12', 'R45'], ['1/4 keV (R1+R2)', '3/4 keV (R4+R5)']):
        for c, (d0, d1) in zip(cols, BANDS):
            g = prof[(prof['band'] == f'{d0:g}-{d1:g}') & (prof['map'] == key)]
            ax.plot(g.dra_deg, g['mean'], color=c, lw=1.6, label=f'Dec {d0:g} to {d1:g}')
            ax.fill_between(g.dra_deg, g['mean'] - g.err, g['mean'] + g.err, color=c, alpha=0.18, lw=0)
        ax.axvline(0, color=OIII_C, lw=1.2, ls='-', alpha=0.9)
        yl = ax.get_ylim()
        ax.set_ylim(yl[0], yl[1] + 0.30 * (yl[1] - yl[0]))          # headroom for the legend
        yl = ax.get_ylim()
        label(ax, 0.3, yl[0] + 0.03 * (yl[1] - yl[0]), 'Silk spine', color=OIII_C, size=8, va='bottom', ha='left')
        label(ax, -14.5, yl[1] - 0.03 * (yl[1] - yl[0]), 'west: Monogem interior', color=FGDIM, size=8, va='top')
        ax.set_xlabel('offset from spine, dRA cos(Dec) (deg; + = east)')
        ax.set_ylabel('1e-6 counts s$^{-1}$ arcmin$^{-2}$')
        ax.set_title(ttl, loc='left', color=FG)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5, loc='upper right', ncol=2)
    finish(fig, ROOT / 'figures/context_rosat_profiles.png')


if __name__ == '__main__':
    main()
