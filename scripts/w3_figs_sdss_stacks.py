"""SDSS/BOSS ON-OFF total-sky stacks for F4, with the wave-2 (verify_sdss) method and calibration.

Uses the unchanged verify_sdss functions: verify_sdss_results.load / good_groups and verify_sdss_figures.stacks
(group-weighted mean of [mean ON - clipped mean OFF] of flux+sky per plate-spectrograph group, diffuse-corrected with
f = 0.84 (3" SDSS) / 0.60 (2" BOSS); ON = within 0.15 deg of a strand, OFF = 0.45-2 deg, Ha-patch 30 % contour
excluded from OFF; group weights from the OFF [O III] scatter), as in verify_sdss_figures.main.
Regions: south (Dec < 20), south176 (Dec < 17.6) and north (Dec >= 20). A BOSS-only stack is also made because [O II]3727 is covered
only by 2" BOSS/SDSS-V fibres, so [O II] must be normalised by the [O III] of the same groups.
Output: data/derived/w3_figs_sdss_stacks.npz  (R per A, photon units, per fibre)
Run: ~/.venvs/astro313/bin/python -I scripts/w3_figs_sdss_stacks.py
"""
import sys
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_sdss_lib import ROOT, load_windows
from verify_sdss_results import load, good_groups
from verify_sdss_figures import stacks

F_SDSS, F_BOSS = 0.84, 0.60


def main():
    D = load_windows(); lam = D['lam']
    fib, grp, amp, Z = load('base')
    gok, _ = good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', F_SDSS, F_BOSS)
    dec = fib['dec']
    ON = fib['dmin'] < 0.15
    OFF = (fib['dmin'] > 0.45) & (fib['dmin'] < 2.0) & ~fib['inpatch30']
    a = amp['OIII:OIII5007']
    wg = np.full(len(grp['gid']), np.nan)
    for g in np.where(gok)[0]:
        af = a[(fib['gid'] == g) & OFF & np.isfinite(a)]
        if len(af) >= 3:
            sd = 1.4826 * np.median(np.abs(af - np.median(af))) * fac[g]
            wg[g] = 1 / max(sd, 1e-3) ** 2
    boss = gok & (grp['instrument'] == 'BOSS')
    out = {'lam': lam}
    for nm, on in (('south', ON & (dec < 20)), ('south176', ON & (dec < 17.6)), ('north', ON & (dec >= 20))):
        out[f'all_{nm}'] = stacks(D, fib, grp, gok, on, OFF & ~on, fac, wg)
        out[f'boss_{nm}'] = stacks(D, fib, grp, boss, on, OFF & ~on, fac, wg)
        print(nm, 'done', flush=True)
    out['note'] = np.array('verify_sdss method (stacks of verify_sdss_figures), tag base, f 0.84/0.60, R/A photon')
    np.savez_compressed(ROOT / 'data/derived/w3_figs_sdss_stacks.npz', **out)
    print('saved')


if __name__ == '__main__':
    main()
