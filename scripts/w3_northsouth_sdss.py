"""w3_northsouth_sdss.py -- wave-3 (tag northsouth): SDSS/BOSS sky-fibre [O III] brightness, [O II]/[O III],
[O III]/Hb and [O III] velocity in Dec bins along the Silk, with the wave-2-corrected method (verify_sdss):
TOTAL sky (flux + sky), per-fibre differential fits against same plate-spectrograph OFF fibres (0.45-2.0 deg,
Halpha-patch 30 % contour excluded from OFF), diffuse-light factors f = 0.84 (3") / 0.60 (2"), group bootstrap,
[O III] velocity from group ON-OFF stacks with the per-group [O I]5577 airglow zero point.

The functions load / good_groups / region_analysis / oiii_velocity / ratio_summary / lower_limit_ratio are imported
read-only from scripts/verify_sdss_results.py (wave-2 code); the Dec-bin regions and the summary are new.
ON = within 0.15 deg of the nearest strand (E, W, W2) and outside the Halpha-patch 30 % contour.

Outputs: data/derived/w3_northsouth_sdss.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_northsouth_sdss.py
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_sdss_results as vr  # noqa: E402
from verify_sdss_lib import ROOT, load_windows  # noqa: E402

F_SDSS, F_BOSS = 0.84, 0.60
R_ON = 0.15
EDGES = [12.0, 16.5, 18.0, 19.5, 21.0, 22.5, 25.0, 28.6]
ZONES = {'dec<17.6': (12.0, 17.6), '17.6-21': (17.6, 21.0), '>=21': (21.0, 28.6),
         'dec<19.5': (12.0, 19.5), 'dec>=19.5': (19.5, 28.6), 'dec<20': (12.0, 20.0), 'dec>=20': (20.0, 28.6)}
LINES = ['OIII:OIII5007', 'Hb:Hb', 'HaNII:Ha', 'HaNII:NII6583', 'SII:SII6716', 'SII:SII6731', 'OII:OII3726',
         'OII:OII3729', 'B4363:OIII4363', 'HeII:HeII4686']


def ratios_for(res, boots, st, w, gids, fac):
    rat = {}
    rat['OIII5007/Hb_photon'] = vr.ratio_summary(boots['OIII:OIII5007'], boots['Hb:Hb'], res['OIII:OIII5007']['I_R'],
                                                 res['Hb:Hb']['I_R'])
    rat['OIII5007/Hb_lower2sig_boot'] = vr.lower_limit_ratio(boots['OIII:OIII5007'], boots['Hb:Hb'])
    rat['NII6583/Ha'] = vr.ratio_summary(boots['HaNII:NII6583'], boots['HaNII:Ha'], res['HaNII:NII6583']['I_R'],
                                         res['HaNII:Ha']['I_R'])
    bS = boots['SII:SII6716'] + boots['SII:SII6731']
    rat['SII/Ha'] = vr.ratio_summary(bS, boots['HaNII:Ha'], res['SII:SII6716']['I_R'] + res['SII:SII6731']['I_R'],
                                     res['HaNII:Ha']['I_R'])
    if res['OII:OII3726']['n_groups'] >= 2:
        m_o2 = st['OII:OII3726'][0]
        gb = gids[np.isfinite(m_o2[gids])]
        o3 = st['OIII:OIII5007'][0] * fac
        o3b = vr.combine(o3, w, gb)
        pb = []
        for _ in range(vr.NB):
            ii = vr.RNG.choice(gb, len(gb))
            pb.append((vr.combine(st['OII:OII3726'][0] * fac, w, ii) + vr.combine(st['OII:OII3729'][0] * fac, w, ii))
                      / vr.combine(o3, w, ii))
        pb = np.array(pb) * 5008.24 / 3728.5
        o2v = res['OII:OII3726']['I_R'] + res['OII:OII3729']['I_R']
        rat['OII/OIII5007_energy_BOSS'] = dict(value=float(o2v / o3b * 5008.24 / 3728.5),
                                               p2p5=float(np.nanpercentile(pb, 2.5)),
                                               p16=float(np.nanpercentile(pb, 16)), p84=float(np.nanpercentile(pb, 84)),
                                               p97p5=float(np.nanpercentile(pb, 97.5)),
                                               OIII5007_same_groups_R=float(o3b), OII_R=float(o2v),
                                               n_groups_boss=int(len(gb)))
        rat['OII3729/3726'] = vr.ratio_summary(boots['OII:OII3729'], boots['OII:OII3726'], res['OII:OII3729']['I_R'],
                                               res['OII:OII3726']['I_R'])
    return rat


def main():
    fib, grp, amp, Z = vr.load('base')
    D = load_windows()
    gok, med = vr.good_groups(grp)
    off = (fib['dmin'] > float(Z['off_lo'])) & (fib['dmin'] < float(Z['off_hi']))
    if int(Z['excl_patch']):
        off &= ~fib['inpatch30']
    fac = np.where(grp['instrument'] == 'SDSS', F_SDSS, F_BOSS)
    zp = grp['v_OI5577'] + grp['helio_rv']
    dec = fib['dec']
    ON = (fib['dmin'] < R_ON) & ~fib['inpatch30']
    regions = {}
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        regions[f'bin_{lo:g}-{hi:g}'] = ON & (dec >= lo) & (dec < hi)
    for nm, (lo, hi) in ZONES.items():
        regions[f'zone_{nm}'] = ON & (dec >= lo) & (dec < hi)
    regions['halpha_patch_30pct'] = fib['inpatch30']
    out = dict(method='verify_sdss total-sky ON-OFF (diff tag base), r_on 0.15 deg, patch excluded from ON',
               f_sdss=F_SDSS, f_boss=F_BOSS, edges=EDGES, zones=ZONES, regions={})
    for nm, sel in regions.items():
        r = vr.region_analysis(fib, grp, amp, gok, sel, off & ~sel, fac, LINES)
        if r is None:
            out['regions'][nm] = None
            print(nm, 'no data')
            continue
        res, boots, st, w, gids = r
        res['ratios'] = ratios_for(res, boots, st, w, gids, fac)
        res['dec_mean_on'] = float(np.mean(dec[sel]))
        try:
            vres, _ = vr.oiii_velocity(fib, grp, gok, sel, off & ~sel, w, gids, D, zp, nboot=300)
            res['velocity_OIII'] = vres
        except Exception as e:  # noqa: BLE001
            res['velocity_OIII'] = {'error': repr(e)}
        out['regions'][nm] = res
        o = res['OIII:OIII5007']
        rr = res['ratios'].get('OII/OIII5007_energy_BOSS', {})
        v = res['velocity_OIII']
        print(f"{nm:22s} ng={res['n_groups']:3d} non={res['n_on_fibres']:4d} [OIII]={o['I_R']:.2f}+-{o['err_boot']:.2f} "
              f"OII/OIII={rr.get('value', np.nan):.2f} ({rr.get('p16', np.nan):.2f}-{rr.get('p84', np.nan):.2f}) "
              f"O3/Hb={res['ratios']['OIII5007/Hb_photon']['value']:.1f} >{res['ratios']['OIII5007/Hb_lower2sig_boot']:.1f} "
              + (f"v={v['v_corr']:+.1f}+-{v['ev_boot']:.1f} (ng {v['n_groups']})" if 'v_corr' in v else str(v)),
              flush=True)
    json.dump(out, open(ROOT / 'data/derived/w3_northsouth_sdss.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
