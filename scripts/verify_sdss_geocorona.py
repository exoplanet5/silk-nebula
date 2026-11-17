"""verify_sdss: absolute total-sky Halpha vs solar depression (geocoronal Halpha decreases as the Earth's shadow rises).

Groups with a TAI timestamp (SDSS-I/II and SDSS-V DR19).  I(Halpha)_total (diffuse-corrected) = Galactic + geocoronal;
the Galactic part is compared with WHAM (2.9-4.4 R near the strand).  Writes data/derived/verify_sdss_geocorona.json
Usage: python -I verify_sdss_geocorona.py f_sdss f_boss
"""
import sys, csv, json
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun, SkyCoord
import astropy.units as u
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT

APO = EarthLocation.from_geodetic(lon=-105.8202 * u.deg, lat=32.7802 * u.deg, height=2788 * u.m)


def main(f_sdss, f_boss):
    G = list(csv.DictReader(open(ROOT / 'data/derived/verify_sdss_groups.csv')))
    def c(k):
        o = []
        for r in G:
            try:
                o.append(float(r[k]))
            except (ValueError, TypeError):
                o.append(np.nan)
        return np.array(o)
    tai = c('tai'); ha = c('I_Ha_abs_cal'); inst = np.array([r['instrument'] for r in G]); sb = c('sb5070_AB_cal'); nf = c('n_fib')
    ra, dec = c('ra'), c('dec')
    f = np.where(inst == 'SDSS', f_sdss, f_boss)
    ok = np.isfinite(tai) & (tai > 0) & np.isfinite(ha) & (nf >= 6)
    med = {i: np.nanmedian(sb[(inst == i) & (nf >= 6)]) for i in ('SDSS', 'BOSS')}
    ok &= sb > np.array([med[i] for i in inst]) - 0.75
    t = Time(tai[ok] / 86400.0, format='mjd', scale='tai')
    frame = AltAz(obstime=t, location=APO)
    sun_alt = get_sun(t).transform_to(frame).alt.deg
    # shadow height along the line of sight is approximated by the solar depression only
    y = ha[ok] * f[ok]
    dep = -sun_alt
    A = np.c_[np.ones(ok.sum()), dep]
    cf = np.linalg.lstsq(A, y, rcond=None)[0]
    bins = [(18, 30), (30, 45), (45, 60), (60, 90)]
    tab = []
    for lo, hi in bins:
        s = (dep >= lo) & (dep < hi)
        if s.sum() >= 3:
            tab.append({'dep_range': [lo, hi], 'n': int(s.sum()), 'median_Ha_R': float(np.median(y[s])),
                        'p16': float(np.percentile(y[s], 16)), 'p84': float(np.percentile(y[s], 84))})
    out = {'n': int(ok.sum()), 'f_sdss': f_sdss, 'f_boss': f_boss, 'linear_fit_Ha_vs_depression': {'intercept': float(cf[0]), 'slope_per_deg': float(cf[1])},
           'by_depression': tab, 'WHAM_galactic_R': [2.9, 4.4],
           'note': 'total = Galactic + geocoronal; geocoronal Halpha at large solar depression is of order a few R (approximate expectation from WHAM geocoronal studies quoted from memory, not verified here)'}
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_geocorona.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]))
