"""Retry the GALFA-HI cutout tiles that were rate-limited (HTTP 429): wait, then fetch slowly."""
import sys, time
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_fetch import fetch, RAW
API = 'https://galfa-dr2-cutout-dptnqjtpkq-uc.a.run.app/cutout'
time.sleep(300)
ra = 108 + 17 / 6 * 5
for dec in [13.75, 20.25, 26.75]:
    p = {'product': 'wide', 'ra': round(ra, 4), 'dec': dec, 'width': round(17 / 3, 4), 'height': 6.5, 'v_lo': -30, 'v_hi': 40}
    for attempt in range(4):
        r = fetch(API, RAW / 'galfa_hi_dr2' / f'galfa_dr2_wide_ra{ra:.2f}_dec{dec:.2f}_vm30p40.fits', 'galfa_hi_dr2',
                  'GALFA-HI DR2 (Peek et al. 2018, ApJS 234, 2) wide product via www.galfa-hi.space cutout service',
                  params=p, timeout=1800, retries=1)
        if r is not None:
            break
        time.sleep(240)
    time.sleep(90)
