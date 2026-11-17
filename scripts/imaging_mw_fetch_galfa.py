"""Fetch GALFA-HI DR2 wide-product cutouts (V_LSR -30..+40 km/s) tiling RA 108-125, Dec 4-30.

Full-field single request returned HTTP 500; tiles of ~5.67 x 6.5 deg work. Same 1' RA-Dec grid, no resampling.
"""
import sys, time
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from imaging_mw_fetch import fetch, RAW
API = 'https://galfa-dr2-cutout-dptnqjtpkq-uc.a.run.app/cutout'
ras = [108 + 17 / 6 * (2 * i + 1) for i in range(3)]       # 110.83, 116.5, 122.17 ; width 5.667 deg
decs = [4 + 6.5 * (j + 0.5) for j in range(4)]              # 7.25, 13.75, 20.25, 26.75 ; height 6.5 deg
for ra in ras:
    for dec in decs:
        p = {'product': 'wide', 'ra': round(ra, 4), 'dec': round(dec, 4), 'width': round(17 / 3, 4), 'height': 6.5,
             'v_lo': -30, 'v_hi': 40}
        fetch(API, RAW / 'galfa_hi_dr2' / f'galfa_dr2_wide_ra{ra:.2f}_dec{dec:.2f}_vm30p40.fits', 'galfa_hi_dr2',
              'GALFA-HI DR2 (Peek et al. 2018, ApJS 234, 2) wide product via www.galfa-hi.space cutout service',
              params=p, timeout=1800)
        time.sleep(3)
