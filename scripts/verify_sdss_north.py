"""verify_sdss: selection-variant summary of strand line ratios (south / north / patch) and forest figure.

Reads data/derived/verify_sdss_results_<tag>_r<r>.json for several variants and verify_sdss_robust_base.json
(gradient-robust and placebo), writes data/derived/verify_sdss_variants.json and figures/verify_sdss_north.png.
"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import silk_style as ss
from verify_sdss_lib import ROOT

VARIANTS = [('base', 0.15, 'base: ON<0.15°, OFF 0.45–2°'), ('base', 0.10, 'ON<0.10°'), ('base', 0.25, 'ON<0.25°'),
            ('offnear', 0.15, 'OFF 0.3–1.0°'), ('offfar', 0.15, 'OFF 1.0–2.5°'), ('plate', 0.15, 'OFF from whole plate'),
            ('inclpatch', 0.15, 'Hα patch kept in OFF'), ('notmpl', 0.15, 'no continuum/OH templates')]
REGS = ['south_dec_lt20', 'W_strand_15_20', 'north_dec_ge20', 'north_dec_20_24', 'north_dec_ge24', 'north_dec_ge20_excl_patch30',
        'halpha_patch_30pct']
RAT = ['OIII5007/Hb', 'OIII5007/Ha', 'NII6583/Ha', 'SII/Ha', 'OII/OIII5007_energy_BOSS']


def main():
    out = {}
    for tag, r, lab in VARIANTS:
        p = ROOT / f'data/derived/verify_sdss_results_{tag}_r{r}.json'
        if not p.exists():
            continue
        d = json.load(open(p))
        o = {}
        for reg in REGS:
            x = d['regions'].get(reg)
            if not x:
                continue
            o[reg] = {'OIII': x['OIII:OIII5007']['I_R'], 'OIII_err': x['OIII:OIII5007']['err_boot'],
                      'Hb': x['Hb:Hb']['I_R'], 'Hb_err': x['Hb:Hb']['err_boot'],
                      'Ha': x['HaNII:Ha']['I_R'], 'Ha_err': x['HaNII:Ha']['err_boot'],
                      'NII': x['HaNII:NII6583']['I_R'], 'NII_err': x['HaNII:NII6583']['err_boot'],
                      'SII': x['SII:SII6716']['I_R'] + x['SII:SII6731']['I_R'],
                      'ratios': {k: x['ratios_photon'].get(k) for k in RAT + ['OIII5007/Hb_lower2sig']},
                      'v': x.get('velocity_OIII', {}).get('v_corr'), 'v_err': x.get('velocity_OIII', {}).get('ev_boot')}
        out[f'{tag}_r{r}'] = {'label': lab, 'regions': o}
    rb = ROOT / 'data/derived/verify_sdss_robust_base.json'
    if rb.exists():
        R = json.load(open(rb))
        out['gradient_robust'] = R['gradient_robust']
        out['placebo_summary'] = R['placebo']['summary']
    json.dump(out, open(ROOT / 'data/derived/verify_sdss_variants.json', 'w'), indent=1, default=float)
    # forest plot
    fig, ax = ss.new_sky_figure(13, 6.4, 'SDSS sky fibres · strand line ratios under selection variants',
                                'photon ratios, ON−OFF; bars = group bootstrap 16–84 % (placebo systematics not included); DESI → photon units', ncols=3)
    panels = [('NII6583/Ha', '[N II]6583 / Hα', (0, 3.5)), ('SII/Ha', '[S II] / Hα', (0, 3.5)), ('OIII5007/Hb', '[O III]5007 / Hβ', (0, 60))]
    cols = {'south_dec_lt20': ss.OIII_C, 'north_dec_ge20': ss.DANGER, 'halpha_patch_30pct': ss.WARN}
    for k, (rk, title, xl) in enumerate(panels):
        a = ax[k]; y = 0; yt = []; yl = []
        for reg in ('south_dec_lt20', 'north_dec_ge20', 'halpha_patch_30pct'):
            for vk, v in out.items():
                if vk in ('gradient_robust', 'placebo_summary') or reg not in v['regions']:
                    continue
                q = v['regions'][reg]['ratios'][rk]
                if q is None:
                    continue
                lo, hi = q.get('p16', np.nan), q.get('p84', np.nan)
                a.errorbar(q['value'], y, xerr=[[max(q['value'] - lo, 0)], [max(hi - q['value'], 0)]], fmt='o', ms=3.5, color=cols[reg], capsize=1.5)
                yt.append(y); yl.append(f"{reg.split('_')[0]} · {v['label']}" if k == 0 else '')
                y -= 1
            y -= 0.8
        # DESI references (energy -> photon: x lambda_num/lambda_den)
        if rk == 'NII6583/Ha':
            refs = [(0.65 * 6585.3 / 6564.6, 'DESI patch'), (0.53 * 6585.3 / 6564.6, 'DESI +resid fibres'), (0.99 * 1.003, 'DESI north arc (Dec≥22)')]
        elif rk == 'SII/Ha':
            refs = [(0.83 * 6725 / 6564.6, 'DESI patch'), (0.40 * 6725 / 6564.6, 'DESI +resid fibres'), (1.03 * 1.024, 'DESI north arc (Dec≥22)')]
        else:
            refs = [(20 * 5008.2 / 4862.7, 'DESI all')]
        for xv, lab in refs:
            a.axvline(xv, color=ss.OK, lw=0.8, ls='--')
            ss.label(a, xv, y + 0.3, ' ' + lab, size=6.5, rotation=90, va='bottom', ha='right', color=ss.OK)
        a.set_title(title); a.set_xlim(*xl)
        a.set_yticks(yt); a.set_yticklabels(yl, fontsize=6.5)
        a.set_ylim(y, 1.0)
    ss.finish(fig, ROOT / 'figures/verify_sdss_north.png')
    for vk, v in out.items():
        if vk in ('gradient_robust', 'placebo_summary'):
            continue
        s = []
        for reg in ('south_dec_lt20', 'north_dec_ge20', 'halpha_patch_30pct'):
            if reg in v['regions']:
                q = v['regions'][reg]['ratios']
                s.append(f"{reg[:5]}: O3/Hb>{q['OIII5007/Hb_lower2sig']:.1f} N2/Ha={q['NII6583/Ha']['value']:.2f} S2/Ha={q['SII/Ha']['value']:.2f} "
                         f"O3={v['regions'][reg]['OIII']:.2f} Ha={v['regions'][reg]['Ha']:.2f}")
        print(f"{v['label']:34s}", ' | '.join(s))


if __name__ == '__main__':
    main()
