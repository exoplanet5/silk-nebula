"""verify_sdss stage 4: independent checks of the absolute (diffuse-light) calibration of SDSS sky-fibre spectra.

SDSS/BOSS spectrophotometry is tied to point-source (PSF) photometry of standard stars, so a uniform diffuse
source is reported too bright by 1/f, f = fraction of a standard star's light entering the fibre during the
exposure (seeing + guiding + ADR).  Three independent anchors for f (= true SB / as-calibrated SB):
  T1  SDSS imaging sky (Field.sky_g/sky_r, nmgy/arcsec^2, DR18 SkyServer table already on disk) near each plate
      versus the synthetic g/r SB of the same-area total sky spectrum (as calibrated).  Imaging sky is a correctly
      calibrated surface brightness; epochs differ (airglow varies), so only ensemble medians are meaningful.
  T2  The solar-spectrum (Fraunhofer-line) component of the continuum (fit coefficient 'sun_coef_blue', R/A at
      ~4965 A) versus the zodiacal-light prediction (Leinert+1998 via STScI ETC table) for the date and position:
      regression across plates observed at different solar elongations; slope = 1/f (x atmospheric diffuse-
      transmission factor), intercept = starlight (ISL) / f.  Insensitive to the airglow continuum.
  T3  BOSS (2") vs SDSS-I/II (3") as-calibrated SB ratio of the same components (relative aperture factor).
Also reported: [O I]5577 intensities (zenith-reduced with a van Rhijn factor for a 97 km layer) and the absolute
Halpha of the total sky (Galactic + geocoronal) compared with WHAM Galactic Halpha.
Writes data/derived/verify_sdss_abscal.json
"""
import sys, csv, json
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_sun, GeocentricTrueEcliptic
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT
from astropy.io import fits


def read_csv(p):
    with open(p) as fh:
        return [r for r in csv.DictReader(l for l in fh if not l.startswith('#'))]


def fl(rows, k):
    out = []
    for r in rows:
        try:
            out.append(float(r[k]))
        except (ValueError, TypeError, KeyError):
            out.append(np.nan)
    return np.array(out)


def zodi_interp():
    Z = json.load(open(ROOT / 'data/derived/verify_sdss_zodi_table_stsci_etc.json'))
    lp = np.array(Z['rows_lambda_prime_deg'], float); be = np.array(Z['cols_beta_deg'], float)
    V = np.array([[np.nan if v is None else v for v in row] for row in Z['V']], float)
    return RegularGridInterpolator((lp, be), V, bounds_error=False)


def solar_V_to_RperA(V, lam=4965.0):
    """Zodiacal light of V mag/arcsec^2 with a solar spectrum -> photon SB in R/A at lam (above atmosphere)."""
    with fits.open(ROOT / 'data/raw/calspec_sun/sun_mod_001.fits') as h:
        w = np.asarray(h[1].data['WAVELENGTH'], float); f = np.asarray(h[1].data['FLUX'], float)
    # Johnson V approximated by Bessell (1990) V tabulated response (coarse), Vega f_lambda(V=0) = 3.63e-9 at 5450 A
    bw = np.array([4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200,
                   6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000])
    bt = np.array([0.000, 0.030, 0.163, 0.458, 0.780, 0.967, 1.000, 0.973, 0.898, 0.792, 0.684, 0.574, 0.461, 0.359,
                   0.270, 0.197, 0.135, 0.081, 0.045, 0.025, 0.017, 0.013, 0.009, 0.000])
    T = np.interp(w, bw, bt, left=0, right=0)
    fV = np.trapezoid(f * T * w, w) / np.trapezoid(T * w, w)          # photon-weighted mean f_lambda in V
    m = (w > lam - 50) & (w < lam + 50)
    fl_lam = np.mean(f[m])
    ratio = fl_lam / fV
    flam = 3.63e-9 * 10 ** (-0.4 * V) * ratio                         # erg/s/cm2/A/arcsec2 at lam
    ph = flam / (6.62607015e-27 * 2.99792458e18 / lam)
    return ph / (1e6 / (4 * np.pi) * (np.pi / 648000.0) ** 2), float(ratio)


def main():
    G = read_csv(ROOT / 'data/derived/verify_sdss_groups.csv')
    g = {k: fl(G, k) for k in G[0].keys() if k not in ('release', 'instrument')}
    inst = np.array([r['instrument'] for r in G]); rel = np.array([r['release'] for r in G])
    ok = (g['n_fib'] >= 6) & np.isfinite(g['sb5070_AB_cal'])
    med = {i: np.nanmedian(g['sb5070_AB_cal'][ok & (inst == i)]) for i in ('SDSS', 'BOSS')}
    dark = ok & (g['sb5070_AB_cal'] > np.array([med[i] for i in inst]) - 0.75)
    out = {'n_groups': int(len(G)), 'n_dark': int(dark.sum()), 'sb5070_median_cal': med}
    # ---------------- T1 imaging sky
    Fd = read_csv(ROOT / 'data/raw/sdss_dr18/fields_box.csv')
    fra, fdec = fl(Fd, 'ra'), fl(Fd, 'dec'); sg, sr = fl(Fd, 'sky_g'), fl(Fd, 'sky_r'); fmjd = fl(Fd, 'mjd_g')
    fq = fl(Fd, 'quality')
    fc = SkyCoord(fra * u.deg, fdec * u.deg)
    gc = SkyCoord(g['ra'] * u.deg, g['dec'] * u.deg)
    img_g = np.full(len(G), np.nan); img_r = np.full(len(G), np.nan)
    for i in range(len(G)):
        s = (fc.separation(gc[i]).deg < 1.0) & (sg > 0) & (sr > 0)
        if s.sum() >= 5:
            img_g[i] = 22.5 - 2.5 * np.log10(np.median(sg[s])); img_r[i] = 22.5 - 2.5 * np.log10(np.median(sr[s]))
    t1 = {}
    for ins in ('SDSS', 'BOSS'):
        for b, im in (('g', img_g), ('r', img_r)):
            s = dark & (inst == ins) & np.isfinite(im) & np.isfinite(g[f'sb_{b}_cal'])
            d = g[f'sb_{b}_cal'][s] - im[s]                 # negative = spectro brighter
            t1[f'{ins}_{b}'] = {'n': int(s.sum()), 'median_dmag_spec_minus_img': float(np.median(d)),
                                'p16': float(np.percentile(d, 16)), 'p84': float(np.percentile(d, 84)),
                                'f_implied_median': float(10 ** (0.4 * np.median(d))),
                                'spec_median_mag': float(np.median(g[f'sb_{b}_cal'][s])), 'img_median_mag': float(np.median(im[s]))}
    # epoch-matched: SDSS-I/II spectra and imaging both in MJD 51900-54600
    allimg = {'g_median_mag_box': float(22.5 - 2.5 * np.log10(np.median(sg[(sg > 0)]))),
              'r_median_mag_box': float(22.5 - 2.5 * np.log10(np.median(sr[(sr > 0)]))),
              'mjd_range_img': [float(np.nanmin(fmjd)), float(np.nanmax(fmjd))],
              'mjd_p16_p84_img': [float(np.nanpercentile(fmjd, 16)), float(np.nanpercentile(fmjd, 84))]}
    t1['imaging_box'] = allimg
    out['T1_imaging_sky'] = t1
    # ---------------- T2 zodiacal regression
    Zi = zodi_interp()
    t = Time(g['mjd'] + 0.30, format='mjd')                     # ~ local midnight at APO (UT 07h)
    sun = get_sun(t).transform_to(GeocentricTrueEcliptic(equinox=t))
    fe = gc.transform_to(GeocentricTrueEcliptic(equinox=t))
    lp = np.abs(((fe.lon.deg - sun.lon.deg) + 180) % 360 - 180)
    beta = np.abs(fe.lat.deg)
    Vz = Zi(np.c_[lp, beta])
    Rz, ratio_sun = solar_V_to_RperA(Vz)
    sunc = g['sun_coef_blue']
    t2 = {'solar_flam4965_over_V': ratio_sun, 'lambda_prime_range': [float(np.nanmin(lp[dark])), float(np.nanmax(lp[dark]))],
          'beta_range': [float(np.nanmin(beta[dark])), float(np.nanmax(beta[dark]))],
          'Vzodi_range': [float(np.nanmin(Vz[dark])), float(np.nanmax(Vz[dark]))]}
    for ins in ('SDSS', 'BOSS'):
        s = dark & (inst == ins) & np.isfinite(sunc) & np.isfinite(Rz) & (sunc > 0)
        x = Rz[s]; y = sunc[s]
        A = np.c_[x, np.ones_like(x)]
        cf = np.linalg.lstsq(A, y, rcond=None)[0]
        # robust: bootstrap over groups
        rng = np.random.default_rng(1)
        bs = []
        for _ in range(1000):
            ii = rng.integers(0, len(x), len(x))
            bs.append(np.linalg.lstsq(A[ii], y[ii], rcond=None)[0])
        bs = np.array(bs)
        t2[ins] = {'n': int(s.sum()), 'slope_meas_per_zodi': float(cf[0]), 'slope_err': float(np.std(bs[:, 0])),
                   'intercept_RperA': float(cf[1]), 'intercept_err': float(np.std(bs[:, 1])),
                   'median_ratio_meas_over_zodi': float(np.median(y / x)),
                   'median_zodi_RperA': float(np.median(x)), 'median_sun_component_RperA_cal': float(np.median(y)),
                   'corr_coef': float(np.corrcoef(x, y)[0, 1])}
    out['T2_zodiacal'] = t2
    # ---------------- T3 BOSS vs SDSS (dark groups, medians)
    t3 = {}
    for k in ('sb5070_AB_cal', 'sb_g_cal', 'sb_r_cal'):
        a = g[k][dark & (inst == 'SDSS')]; b = g[k][dark & (inst == 'BOSS')]
        t3[k] = {'SDSS_median': float(np.nanmedian(a)), 'BOSS_median': float(np.nanmedian(b)),
                 'BOSS_minus_SDSS_mag': float(np.nanmedian(b) - np.nanmedian(a))}
    for k in ('sun_coef_blue', 'I_OI5577_cal', 'I_Ha_abs_cal', 'add_coef_blue'):
        a = g[k][dark & (inst == 'SDSS')]; b = g[k][dark & (inst == 'BOSS')]
        t3[k] = {'SDSS_median': float(np.nanmedian(a)), 'BOSS_median': float(np.nanmedian(b)), 'BOSS_over_SDSS': float(np.nanmedian(b) / np.nanmedian(a))}
    # zodi-normalised solar component ratio (removes the elongation mix)
    rz = sunc / Rz
    t3['sun_over_zodi'] = {'SDSS_median': float(np.nanmedian(rz[dark & (inst == 'SDSS')])),
                           'BOSS_median': float(np.nanmedian(rz[dark & (inst == 'BOSS')])),
                           'BOSS_over_SDSS': float(np.nanmedian(rz[dark & (inst == 'BOSS')]) / np.nanmedian(rz[dark & (inst == 'SDSS')]))}
    out['T3_boss_vs_sdss'] = t3
    # ---------------- airglow / Halpha
    X = g['airmass']
    # van Rhijn factor for emission height h=97 km: V(z) = 1/sqrt(1 - (R/(R+h))^2 sin^2 z)
    Re = 6371.0
    def vr(Xa, h=97.0):
        z = np.arccos(1 / Xa)
        return 1 / np.sqrt(1 - (Re / (Re + h)) ** 2 * np.sin(z) ** 2)
    Xf = np.where(np.isfinite(X) & (X >= 1), X, 1.2)
    I5577 = g['I_OI5577_cal']
    o = {}
    for ins in ('SDSS', 'BOSS'):
        s = dark & (inst == ins)
        o[ins] = {'I5577_cal_median': float(np.nanmedian(I5577[s])), 'I5577_cal_p16_p84': [float(np.nanpercentile(I5577[s], 16)), float(np.nanpercentile(I5577[s], 84))],
                  'I5577_cal_zenith_median': float(np.nanmedian(I5577[s] / vr(Xf[s]))),
                  'airmass_median_known': float(np.nanmedian(X[s])) if np.isfinite(X[s]).any() else None,
                  'IHa_abs_cal_median': float(np.nanmedian(g['I_Ha_abs_cal'][s])),
                  'IHa_abs_cal_p16_p84': [float(np.nanpercentile(g['I_Ha_abs_cal'][s], 16)), float(np.nanpercentile(g['I_Ha_abs_cal'][s], 84))],
                  'I_OH6555_cal_median': float(np.nanmedian(g['I_OH6555_cal'][s])),
                  'I_NI5198_cal_median': float(np.nanmedian(g['I_NI5198_cal'][s])),
                  'I_NI5200_cal_median': float(np.nanmedian(g['I_NI5200_cal'][s]))}
    out['airglow_halpha'] = o
    # per-group table for figures
    np.savez_compressed(ROOT / 'data/derived/verify_sdss_abscal_groups.npz', inst=inst, dark=dark, img_g=img_g, img_r=img_r,
                        sb_g=g['sb_g_cal'], sb_r=g['sb_r_cal'], sb5070=g['sb5070_AB_cal'], Rz=Rz, Vz=Vz, lp=lp, beta=beta,
                        sunc=sunc, mjd=g['mjd'], I5577=I5577, IHa=g['I_Ha_abs_cal'], airmass=X, seeing=g['seeing50'])
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_abscal.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
