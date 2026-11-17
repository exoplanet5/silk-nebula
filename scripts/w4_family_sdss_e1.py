"""w4_family_sdss_e1.py -- wave-4 `family`: SDSS/BOSS sky-fibre ON-OFF measurement of the E1 [O III] filament
(~3 deg east of the Silk; NSNS ridge in data/derived/w4_family_traces.json) with the corrected wave-2 method.

Re-uses, unchanged and by import, the wave-2 per-fibre TOTAL-sky differential line fits
(data/derived/verify_sdss_diff_base.npz: total sky minus the clipped mean of the group's Silk-OFF fibres, fixed
velocity +13 km/s, continuum template) and the group/bootstrap/velocity code of scripts/verify_sdss_results.py.
Because those per-fibre amplitudes are linear in the spectrum, ON - OFF with a new OFF set is exact:
  ON  = d_E1 < r_on (0.10 and 0.15 deg)
  OFF = 0.45 < d_E1 < 2.0 deg, > 0.45 deg from every Silk strand (dmin), outside the Halpha-patch 30 % contour.
Diffuse-light factors f_eff = 0.84 (SDSS 3") / 0.60 (BOSS 2") (verify_sdss). Ratios photon; energy where stated.
Placebo: the E1 ridge shifted in RA by +-0.6..+-1.6 deg (ON/OFF redefined, same code) -> scatter of null ON-OFF.
Velocity: verify_sdss oiii_velocity (free-velocity fit of the weighted ON-OFF group stack, [O I]5577 zero point).
Output: data/derived/w4_family_sdss_e1.json
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_family_sdss_e1.py
"""
import sys
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_sdss_results as vr                      # noqa: E402
from w4_family_select import ridge_tree, dist          # noqa: E402

ROOT = '/Users/mickey/astronomy/silk-nebula'
LINES = ['OIII:OIII5007', 'Hb:Hb', 'HaNII:Ha', 'HaNII:NII6583', 'SII:SII6716', 'SII:SII6731', 'OII:OII3726',
         'OII:OII3729', 'B4363:OIII4363', 'HeII:HeII4686']


def shifted_tree(dra):
    from scipy.spatial import cKDTree
    from w4_family_select import uv
    _, P = ridge_tree('E1')
    P = P.copy(); P[:, 0] += dra / np.cos(np.radians(P[:, 1]))
    return cKDTree(uv(P[:, 0], P[:, 1]))


def analyse(fib, grp, amp, gok, fac, dE, r_on, D=None, zp=None, velocity=False):
    on = dE < r_on
    off = (dE > 0.45) & (dE < 2.0) & (fib['dmin'] > 0.45) & ~fib['inpatch30']
    r = vr.region_analysis(fib, grp, amp, gok, on, off & ~on, fac, LINES)
    if r is None:
        return None
    res, boots, st, w, gids = r
    rat = {}
    rat['OIII5007/Hb'] = vr.ratio_summary(boots['OIII:OIII5007'], boots['Hb:Hb'], res['OIII:OIII5007']['I_R'],
                                          res['Hb:Hb']['I_R'])
    rat['OIII5007/Hb_lower2sig'] = vr.lower_limit_ratio(boots['OIII:OIII5007'], boots['Hb:Hb'])
    rat['OIII5007/Ha'] = vr.ratio_summary(boots['OIII:OIII5007'], boots['HaNII:Ha'], res['OIII:OIII5007']['I_R'],
                                          res['HaNII:Ha']['I_R'])
    m_o2 = st['OII:OII3726'][0]
    gb = gids[np.isfinite(m_o2[gids])]
    if len(gb) >= 2:
        pb = []
        for _ in range(vr.NB):
            ii = vr.RNG.choice(gb, len(gb))
            pb.append((vr.combine(st['OII:OII3726'][0] * fac, w, ii) + vr.combine(st['OII:OII3729'][0] * fac, w, ii)) /
                      vr.combine(st['OIII:OIII5007'][0] * fac, w, ii))
        pb = np.array(pb) * 5008.24 / 3728.5
        o3b = vr.combine(st['OIII:OIII5007'][0] * fac, w, gb)
        o2b = vr.combine(st['OII:OII3726'][0] * fac, w, gb) + vr.combine(st['OII:OII3729'][0] * fac, w, gb)
        rat['OII/OIII5007_energy_BOSS'] = dict(value=float(o2b / o3b * 5008.24 / 3728.5),
                                               p2_5=float(np.nanpercentile(pb, 2.5)), p16=float(np.nanpercentile(pb, 16)),
                                               p84=float(np.nanpercentile(pb, 84)), p97_5=float(np.nanpercentile(pb, 97.5)),
                                               OIII_same_groups_R=float(o3b), OII_R=float(o2b), n_groups_boss=int(len(gb)))
    res['ratios_photon'] = rat
    if velocity:
        try:
            v, _ = vr.oiii_velocity(fib, grp, gok, on, off & ~on, w, gids, D, zp, nboot=300)
            res['velocity_OIII'] = v
        except Exception as ex:
            res['velocity_OIII'] = {'error': repr(ex)}
    return res


def main():
    fib, grp, amp, Z = vr.load('base')
    gok, med = vr.good_groups(grp)
    fac = np.where(grp['instrument'] == 'SDSS', 0.84, 0.60)
    tE, _ = ridge_tree('E1')
    dE = dist(tE, fib['ra'], fib['dec'])
    D = vr.load_windows()
    zp = grp['v_OI5577'] + grp['helio_rv']
    out = dict(method=__doc__.split('Re-uses')[1].split('Output:')[0].strip(), results={}, placebo={})
    for r_on in (0.10, 0.15):
        res = analyse(fib, grp, amp, gok, fac, dE, r_on, D, zp, velocity=True)
        out['results'][f'r{r_on:.2f}'] = res
        o = res['OIII:OIII5007']; hb = res['Hb:Hb']
        print(f'r_on {r_on}: ng {res["n_groups"]} non {res["n_on_fibres"]} [OIII] {o["I_R"]:.3f}+-{o["err_boot"]:.3f} '
              f'Hb {hb["I_R"]:.3f}+-{hb["err_boot"]:.3f} Ha {res["HaNII:Ha"]["I_R"]:.3f}+-{res["HaNII:Ha"]["err_boot"]:.3f} '
              f'OII {res["OII:OII3726"]["I_R"] + res["OII:OII3729"]["I_R"]:.3f} '
              f'O3/Hb>{res["ratios_photon"]["OIII5007/Hb_lower2sig"]:.1f} '
              f'{res["ratios_photon"].get("OII/OIII5007_energy_BOSS", {})} v={res.get("velocity_OIII")}', flush=True)
    # placebo: shifted ridges
    for r_on in (0.10, 0.15):
        vals = {L: [] for L in LINES}
        for dra in (-1.6, -1.3, -1.0, -0.8, -0.6, 0.6, 0.8, 1.0, 1.3, 1.6):
            dP = dist(shifted_tree(dra), fib['ra'], fib['dec'])
            dP = np.where(dE < 0.45, 99.0, dP)            # never let the real E1 enter a placebo ON/OFF
            res = analyse(fib, grp, amp, gok, fac, dP, r_on)
            if res is None or res['n_on_fibres'] < 10:
                continue
            for L in LINES:
                vals[L].append(res[L]['I_R'])
        out['placebo'][f'r{r_on:.2f}'] = {L: dict(n=len(v), mean=float(np.mean(v)) if v else None,
                                                  rms=float(np.sqrt(np.mean(np.square(v)))) if v else None,
                                                  values=[float(x) for x in v]) for L, v in vals.items()}
        print('placebo', r_on, {L: (round(out['placebo'][f'r{r_on:.2f}'][L]['rms'] or np.nan, 3)) for L in LINES})
    json.dump(out, open(f'{ROOT}/data/derived/w4_family_sdss_e1.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
