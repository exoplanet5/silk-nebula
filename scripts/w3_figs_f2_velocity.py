"""F2 velocity reconciliation: heliocentric [O III] and ISM velocities toward the Silk.

Numbers (all wave-2 unless stated; see notes/w3_figs.md for the full provenance list):
  paper Table 4 (sources/arxiv/2610.13329v1.txt): [O III] P1/P2/P3 26.8+-4.3 / 26.8+-4.1 / 23.8+-3.3;
      [N II] 11.8+-1.6 / 12.0+-1.6 / 11.7+-1.9 / 13.3+-2.3 (P4); [S II] 11.8+-2.5 / 12.0+-1.6 / 11.7+-1.9 / 12.7+-2.1;
      Ha 16.3+-1.2 / 17.0+-1.2 / 18.4+-1.4 / 19.7+-1.4
  correction -10.1+-1.0 and Fig. 3a vs [N I] airglow +13.6+-1.2: verify_velocity_summary.json
  DESI: verify_desi.md rev. 2 (all +11.3 [10.2,13.9]; Dec<19.5 +10.6 [10.2,12.0]; Dec>=19.5 +20.9 [18.5,23.1]);
        verify_velocity ridge stack +14.4+-1.9
  SDSS: verify_sdss_results_base_r0.15.json (south Dec<20 +10.0+-1.8 stat +-0.8 zp; north +16.4+-5.9)
  WHAM: verify_hiwham_wham_beams.csv (main-component Gaussian v_g, V_LSR, per beam; on = |d|<0.5 deg)
  HI:   verify_hiwham_hi_extra.json A_along_spine (on-spine GALFA peak velocity per Dec bin, Dec 12-22)
  adopted: verify_velocity combined V_helio +13.3+-1.0
Output: figures/w3_figs_f2_velocity.png, data/derived/w3_figs_f2_velocity.json
"""
import sys
import csv
import json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from w3_figs_common import *  # noqa: F401,F403
from w3_figs_common import DER, FIG, json_dump

LSR_OFF = -12.2   # V_LSR - V_helio toward the southern Silk (standard solar motion; -11.9...-12.6 across P1-P4)


def wham_hi():
    rows = list(csv.DictReader(open(DER / 'verify_hiwham_wham_beams.csv')))
    W = {k: np.array([float(r[k]) for r in rows]) for k in ('d_deg', 's_deg', 'v_c_lsr', 'v_g', 'dec')}
    res = json.loads((DER / 'verify_hiwham_wham_results.json').read_text())
    vh_adopt, s_split = res['v_helio_oiii_adopted'], res['s_split']
    vhel = W['v_g'] - W['v_c_lsr'] + vh_adopt          # beam LSR -> helio using the per-beam LSR offset
    out = {}
    for nm, sel in (('south_on', (np.abs(W['d_deg']) < 0.5) & (W['s_deg'] <= s_split)),
                    ('north_on', (np.abs(W['d_deg']) < 0.5) & (W['s_deg'] > s_split))):
        v = vhel[sel]
        out['wham_' + nm] = dict(n=int(sel.sum()), mean=float(np.mean(v)), std=float(np.std(v, ddof=1)),
                                 sem=float(np.std(v, ddof=1) / np.sqrt(sel.sum())),
                                 v_lsr_mean=float(np.mean(W['v_g'][sel])))
    hx = json.loads((DER / 'verify_hiwham_hi_extra.json').read_text())['A_along_spine']
    db = np.array(hx['dec_bins']); vp = np.array(hx['vpeak_on_kms']); vexp = np.array(hx['vlsr_oiii_expected'])
    s = (db >= 12) & (db <= 22)
    vh = vp[s] - vexp[s] + 13.5            # vlsr_oiii_expected is V_LSR of V_helio = +13.5 at that Dec
    out['hi_peak_dec12_22'] = dict(n=int(s.sum()), median=float(np.median(vh)), p16=float(np.percentile(vh, 16)),
                                   p84=float(np.percentile(vh, 84)), vlsr_median=float(np.median(vp[s])))
    return out


def main():
    I = wham_hi()
    sd = json.loads((DER / 'verify_sdss_results_base_r0.15.json').read_text())['regions']
    vs = sd['south_dec_lt20']['velocity_OIII']; vn = sd['north_dec_ge20']['velocity_OIII']
    vsum = json.loads((DER / 'verify_velocity_summary.json').read_text())
    corr = vsum['table4_offset']; adopt, eadopt = vsum['combined']['V_helio'], vsum['combined']['err']
    f3 = vsum['fig3a']
    T4 = {'P1': (26.8, 4.3), 'P2': (26.8, 4.1), 'P3': (23.8, 3.3)}
    NII = {'P1': (11.8, 1.6), 'P2': (12.0, 1.6), 'P3': (11.7, 1.9), 'P4': (13.3, 2.3)}
    SII = {'P1': (11.8, 2.5), 'P2': (12.0, 1.6), 'P3': (11.7, 1.9), 'P4': (12.7, 2.1)}
    HA = {'P1': (16.3, 1.2), 'P2': (17.0, 1.2), 'P3': (18.4, 1.4), 'P4': (19.7, 1.4)}

    # rows: (group, label, kind, payload, colour, marker)
    R = []
    G1 = 'Wu+26 LAMOST MRS-N'
    for k, (v, e) in T4.items():
        R.append((G1, f'[O III] {k} · Table 4 as published', 'pe', (v, e), DANGER, 'o'))
    for k, (v, e) in T4.items():
        R.append((G1, f'[O III] {k} · Table 4 − {corr:.1f}', 'pe', (v - corr, np.hypot(e, 1.0)), OIII_C, 'o'))
    R.append((G1, '[O III] Fig. 3a vs [N I] airglow', 'pe', (f3['V_helio_corr'], f3['err']), OIII_C, 'D'))
    R.append((G1, '[N II] P1–P4 · Table 4', 'multi', [NII[k] for k in NII], NII_C, 's'))
    R.append((G1, '[S II] P1–P4 · Table 4', 'multi', [SII[k] for k in SII], SII_C, 's'))
    R.append((G1, 'Hα P1–P4 · Table 4 (geocoronal bias)', 'multi', [HA[k] for k in HA], HA_C, 's'))
    G2 = 'DESI DR1 [O III] (sky residuals)'
    R.append((G2, 'on-ridge stack (verify_velocity)', 'pe', (14.4, 1.9), OIII_C, 'o'))
    R.append((G2, 'all fibres · joint 4959+5007', 'pi', (11.3, 10.2, 13.9), OIII_C, 'o'))
    R.append((G2, 'south · Dec < 19.5', 'pi', (10.6, 10.2, 12.0), OIII_C, 'o'))
    R.append((G2, 'north · Dec ≥ 19.5', 'pi', (20.9, 18.5, 23.1), '#b8ecff', 'o'))
    G3 = 'SDSS/BOSS [O III] (total sky ON−OFF)'
    R.append((G3, 'south · Dec < 20', 'pe', (vs['v_corr'], float(np.hypot(vs['ev_boot'], 0.8))), OIII_C, 'o'))
    R.append((G3, 'north · Dec ≥ 20', 'pe', (vn['v_corr'], float(np.hypot(vn['ev_boot'], 0.8))), '#b8ecff', 'o'))
    G4 = 'Local ISM'
    w = I['wham_south_on']
    R.append((G4, f'WHAM Hα main comp. · S on-filament ({w["n"]} beams)', 'pe', (w['mean'], w['std']), HA_C, 'D'))
    wn = I['wham_north_on']
    R.append((G4, f'WHAM Hα main comp. · N on-filament ({wn["n"]} beams)', 'pe', (wn['mean'], wn['std']), HA_C, 'D'))
    h = I['hi_peak_dec12_22']
    R.append((G4, 'GALFA HI on-spine peak · Dec 12–22', 'pi', (h['median'], h['p16'], h['p84']), HB_C, 'D'))

    n = len(R)
    fig = plt.figure(figsize=(11.5, 9.2))
    fig.text(0.012, 0.988, 'Silk Nebula · [O III] velocity reconciliation', color='#cfd6dd', fontsize=13,
             fontweight='bold', va='top')
    fig.text(0.012, 0.962, 'heliocentric velocities; bars = 68 % intervals (Table 4: quoted ±1σ)',
             color=FGDIM, fontsize=9.5, va='top')
    ax = fig.add_axes([0.385, 0.075, 0.585, 0.80])
    y = np.arange(n)[::-1].astype(float)
    # group spacing
    gy = []
    off = 0.0; prev = None
    for i, r in enumerate(R):
        if prev is not None and r[0] != prev:
            off += 0.8
        y[i] = -(i + off)
        prev = r[0]
    ax.axvspan(adopt - eadopt, adopt + eadopt, color=OIII_C, alpha=0.16, lw=0, zorder=0)
    ax.axvline(adopt, color=OIII_C, lw=0.8, alpha=0.7, zorder=0)
    for i, (g, lab, kind, pl, c, mk) in enumerate(R):
        yy = y[i]
        if kind == 'pe':
            v, e = pl
            ax.errorbar(v, yy, xerr=e, fmt=mk, color=c, ms=6, lw=1.4, capsize=2.5, mec='k', mew=0.4, zorder=4)
        elif kind == 'pi':
            v, lo, hi = pl
            ax.errorbar(v, yy, xerr=[[v - lo], [hi - v]], fmt=mk, color=c, ms=6, lw=1.4, capsize=2.5, mec='k',
                        mew=0.4, zorder=4)
        else:
            for j, (v, e) in enumerate(pl):
                ax.errorbar(v, yy + (j - 1.5) * 0.13, xerr=e, fmt=mk, color=c, ms=3.8, lw=0.9, capsize=1.5,
                            alpha=0.9, zorder=4)
        ax.text(-0.012, yy, lab, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=8.6,
                color=FG)
    # arrows from published to corrected
    for k_i in range(3):
        v0 = R[k_i][3][0]; v1 = R[3 + k_i][3][0]
    # group headers
    prev = None
    for i, r in enumerate(R):
        if r[0] != prev:
            ax.text(-0.012, y[i] + 0.72, r[0], transform=ax.get_yaxis_transform(), ha='right', va='center',
                    fontsize=9.2, color=WARN, fontweight='bold')
            prev = r[0]
    ax.set_ylim(y.min() - 0.8, y.max() + 1.3)
    ax.set_yticks([])
    ax.set_xlim(2, 33)
    ax.set_xlabel('V_helio (km s⁻¹)')
    ax.grid(axis='x', alpha=0.25)
    sec = ax.secondary_xaxis('top', functions=(lambda v: v + LSR_OFF, lambda v: v - LSR_OFF))
    sec.set_xlabel('V_LSR (km s⁻¹; −12.2 offset)', fontsize=8.5, color=FGDIM)
    sec.tick_params(labelsize=8, colors=FGDIM)
    label(ax, adopt + eadopt + 0.3, y.max() + 0.85, f'adopted +{adopt:.1f} ± {eadopt:.1f}', color=OIII_C, size=8.5,
          ha='left', va='center')
    label(ax, 29.5, y[4] + 0.5, f'airglow scale + missing\nv_corr: −{corr:.1f} km s⁻¹', color=FGDIM, size=7.5,
          ha='center', va='center')
    ax.annotate('', xy=(R[3][3][0] + 1.0, y[3] - 0.0), xytext=(R[0][3][0] - 1.2, y[0] - 0.0),
                arrowprops=dict(arrowstyle='-|>', color=FGDIM, lw=0.9, connectionstyle='arc3,rad=0.25',
                                mutation_scale=9), zorder=2)
    ss.credit(fig)
    fig.savefig(FIG / 'w3_figs_f2_velocity.png', dpi=200)
    plt.close(fig)
    out = {'adopted_V_helio': adopt, 'adopted_err': eadopt, 'table4_correction': corr, 'LSR_offset_used': LSR_OFF,
           'rows': [dict(group=r[0], label=r[1], kind=r[2], values=r[3]) for r in R], 'ism_derived': I}
    json_dump(out, DER / 'w3_figs_f2_velocity.json')
    print(json.dumps(I, indent=1))


if __name__ == '__main__':
    main()
