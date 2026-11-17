"""verify_distance: figures (silk_style, dpi 200).

  figures/verify_distance_constraints.png  distance ladder of all constraints + sky geometry + kinematics
  figures/verify_distance_hd61717.png       HD 61717: Ca II K/H fit, Na I D after telluric correction, SB2 lines
  figures/verify_distance_iue.png           IUE SWP-HI co-adds toward HD 60848 (and HD 58050) with doublet fits
Run: ~/.venvs/astro313/bin/python -I scripts/verify_distance_figures.py
"""
import json
import sys
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import (new_sky_figure, label, finish, plt, BG2, FG, FGDIM, ACCENT, WARN, DANGER, OK, SUN,  # noqa: E402
                        OIII_C, HA_C, SII_C, NII_C, HB_C)
from silk_geometry import load_trace, ROOT  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

S = json.loads((ROOT / 'data/derived/verify_distance_synthesis.json').read_text())
A = np.load(ROOT / 'data/derived/verify_distance_synthesis_arrays.npz')
H61 = json.loads((ROOT / 'data/derived/verify_distance_hd61717.json').read_text())
P61 = np.load(ROOT / 'data/derived/verify_distance_hd61717_profiles.npz')
IUE = json.loads((ROOT / 'data/derived/verify_distance_iue.json').read_text())
PI = np.load(ROOT / 'data/derived/verify_distance_iue_profiles.npz')
tr = load_trace()
sp_ra = np.array([p['ra'] for p in tr['spine']])
sp_dec = np.array([p['dec'] for p in tr['spine']])

# =============================================================== figure 1: constraints
fig = plt.figure(figsize=(15, 9.3))
fig.text(0.012, 0.985, 'Silk Nebula · distance constraints', color='#cfd6dd', fontsize=13, fontweight='bold', va='top')
fig.text(0.012, 0.955, 'favoured 252–318 pc (Monogem limb, 68 %) · permitted ~0.1–1 kpc · no absorption bracket on the ionized gas',
         color=FGDIM, fontsize=9.5, va='top')
ax = fig.add_axes([0.215, 0.08, 0.40, 0.82])
rows = [
    ('Monogem limb (S spine, PSR parallax)', 'range'),
    ('Monogem shell near / far side', 'range'),
    ('Local Bubble wall (O\'Neill+24)', 'ref'),
    ('HD 61717 −17 km/s cloud, if related', 'upper'),
    ('WD J0747+2106 COS, no Si IV/N V', 'lower'),
    ('Ca II / Na I at [O III] velocity', 'none'),
    ('Kinematic, DESI/SDSS V_LSR≈+1…+3', 'kin'),
    ('Kinematic, LAMOST V_LSR +11…+16', 'kin'),
    ('PSR J0751+1807 (0.95° E) DM', 'none'),
    ('Extinction / 3D dust', 'none'),
    ('IUE HD 60848 high ions (1.9 kpc)', 'none'),
    ('Size / z-height plausibility', 'soft'),
]
ny = len(rows)
ml = S['monogem_limb_mc']['d_limb_pc_percentiles_2p5_16_50_84_97p5']
nf = S['monogem_shell_near_far_mc']
lb = S['local_bubble_oneill24_spine']
kr = S['kinematic']['ranges_kpc']
y = {nm: ny - 1 - i for i, (nm, _) in enumerate(rows)}
XL, XR = 40, 4000
ax.set_xscale('log')
ax.set_xlim(XL, XR)
ax.set_ylim(-0.8, ny - 0.3)
for i in range(ny):
    ax.axhline(i, color='#1e242c', lw=0.6, zorder=0)
# combined bands
ax.axvspan(ml[1], ml[3], color=OK, alpha=0.13, zorder=0)
ax.axvspan(100, 1000, color=FGDIM, alpha=0.05, zorder=0)


def bar(yy, lo, hi, color, hatch=None, alpha=0.85, lw=0):
    ax.add_patch(Rectangle((lo, yy - 0.28), hi - lo, 0.56, facecolor=color if not hatch else 'none', edgecolor=color,
                           hatch=hatch, alpha=alpha, lw=lw if not hatch else 1.0))


def arrow(yy, x0, direction, color):
    x1 = x0 * (2.2 if direction > 0 else 1 / 2.2)
    ax.annotate('', xy=(x1, yy), xytext=(x0, yy), arrowprops=dict(arrowstyle='-|>', color=color, lw=1.8))
    ax.plot([x0, x0], [yy - 0.3, yy + 0.3], color=color, lw=2)


r = y['Monogem limb (S spine, PSR parallax)']
bar(r, ml[0], ml[4], OK, alpha=0.35)
bar(r, ml[1], ml[3], OK)
ax.plot(ml[2], r, 'o', color='white', ms=4, zorder=5)
label(ax, ml[4] * 1.06, r, f'{ml[1]}–{ml[3]} pc (95 %: {ml[0]}–{ml[4]})', size=8, va='center')
r = y['Monogem shell near / far side']
bar(r, nf['near_pc_16_50_84'][0], nf['far_pc_16_50_84'][2], OK, hatch='///')
label(ax, nf['far_pc_16_50_84'][2] * 1.06, r, f"{nf['near_pc_16_50_84'][0]}–{nf['far_pc_16_50_84'][2]} pc", size=8, va='center')
r = y["Local Bubble wall (O'Neill+24)"]
bar(r, lb['inner_pc_median'], lb['outer_pc_median'], SUN, alpha=0.35)
ax.plot([lb['peak_pc_median']] * 2, [r - 0.3, r + 0.3], color=SUN, lw=2.2)
label(ax, lb['outer_pc_median'] * 1.06, r, f"peak {lb['peak_pc_median']} pc · reference only", size=8, va='center')
r = y['HD 61717 −17 km/s cloud, if related']
arrow(r, S['hd61717_distance']['d_pc'], -1, WARN)
label(ax, S['hd61717_distance']['d_pc'] * 1.08, r, '< 272 pc (conditional)', size=8, va='center', color=WARN)
r = y['WD J0747+2106 COS, no Si IV/N V']
arrow(r, 99.5, +1, WARN)
label(ax, 99.5 * 2.4, r, '> 100 pc (conditional)', size=8, va='center', ha='left', color=WARN)
r = y['Ca II / Na I at [O III] velocity']
label(ax, 330, r, 'none: [O III] V = local gas; ionized gas', size=8, va='center', ha='center', color=FGDIM)
for key, nm in (('[O III] DESI (+14.5 helio)', 'Kinematic, DESI/SDSS V_LSR≈+1…+3'),
                ('[O III] LAMOST Wu+2026 (+23.8..+26.8)', 'Kinematic, LAMOST V_LSR +11…+16')):
    lo, hi = kr[key]
    lo = max(lo * 1000, XL)
    bar(y[nm], lo, hi * 1000, FGDIM, hatch='xx', alpha=0.7)
    label(ax, hi * 1000 * 1.05 if hi * 1000 < 2000 else lo / 1.05, y[nm], 'not physical (σ_pec 7 km/s)', size=8,
          va='center', ha='left' if hi * 1000 < 2000 else 'right', color=FGDIM)
r = y['PSR J0751+1807 (0.95° E) DM']
pj = S['pulsars']['J0751_1807']
ax.errorbar(1000 * pj['d_kpc'], r, xerr=[[1000 * (pj['d_kpc'] - pj['d_kpc_68'][0])], [1000 * (pj['d_kpc_68'][1] - pj['d_kpc'])]],
            fmt='D', color=ACCENT, ms=5)
label(ax, 1000 * pj['d_kpc'] / 1.1, r, 'pulsar here; Silk ΔDM ≲ 1 → no bracket', size=8, va='center', ha='right', color=FGDIM)
r = y['Extinction / 3D dust']
label(ax, 330, r, 'none: no dust in the Silk; E(B−V) < 0.006 to 300 pc', size=8, va='center', ha='center', color=FGDIM)
r = y['IUE HD 60848 high ions (1.9 kpc)']
ax.plot(1000 / 0.5227, r, 's', color=ACCENT, ms=5)
label(ax, 1000 / 0.5227 / 1.1, r, 'total Galactic C IV/Si IV only', size=8, va='center', ha='right', color=FGDIM)
r = y['Size / z-height plausibility']
arrow(r, 1000, -1, FGDIM)
label(ax, 1000 * 1.1, r, '1 kpc → 295 pc long, z 330 pc', size=8, va='center', color=FGDIM)
for xv, txt, col in ((160, 'Wu+26 "160 pc"', DANGER), (1600, 'kd output ≈1.6 kpc', DANGER)):
    ax.axvline(xv, color=col, lw=0.9, ls='-', alpha=0.8)
    label(ax, xv * 1.03, ny - 0.55, txt, size=7.5, color=col, va='center')
ax.set_yticks([y[nm] for nm, _ in rows])
ax.set_yticklabels([nm for nm, _ in rows], fontsize=8)
ax.set_xlabel('distance (pc)')
ax.set_xticks([50, 100, 200, 300, 500, 1000, 2000])
ax.set_xticklabels(['50', '100', '200', '300', '500', '1000', '2000'])

# sky geometry
ax2 = fig.add_axes([0.665, 0.42, 0.32, 0.50])
th = np.linspace(0, 2 * np.pi, 400)


def circle(ra0, dec0, rdeg):
    c0 = SkyCoord(ra0 * u.deg, dec0 * u.deg)
    pts = c0.directional_offset_by(th * u.rad, rdeg * u.deg)
    return pts.ra.deg, pts.dec.deg


for (ra0, dec0, R), ls, nm in (((104.15, 15.30, 11.0), '-', 'Monogem (2-SNR)'), ((107.08, 13.50, 11.97), '--', 'single SNR')):
    x, yy = circle(ra0, dec0, R)
    ax2.plot(x, yy, ls, color=ACCENT, lw=0.9, alpha=0.8, label=nm)
ax2.plot(sp_ra, sp_dec, '-', color=OIII_C, lw=3, alpha=0.8, label='Silk spine')
pts = [('HD 61717', 115.394, 22.117, '272 pc', WARN), ('HD 60848', 114.274, 16.904, '1.9 kpc', ACCENT),
       ('WD J0747', 116.900, 21.110, '99.5 pc', WARN), ('J0751+1807', 117.788, 18.127, '1.18 kpc', SUN),
       ('B0656+14', 104.951, 14.239, '288 pc', SUN), ('HD 64705', 119.027, 23.624, '280 pc', FG),
       ('HD 61038', 114.462, 14.042, '149 pc', FG), ('HD 63432', 117.488, 28.199, '250 pc', FG),
       ('HD 55869', 108.5, 21.0, None, None), ('NGC 2392', 112.295, 20.912, '1.8 kpc', ACCENT),
       ('HD 66665', 121.199, 6.186, '1.8 kpc', ACCENT), ('11 CMi', 116.568, 10.768, '92 pc', FG),
       ('HD 52124', 105.163, 16.959, '264 pc', FG), ('HD 53588', 106.573, 17.744, '240 pc', FG)]
for nm, ra, dec, dist, col in pts:
    if col is None:
        continue
    ax2.plot(ra, dec, '*' if 'J' in nm[:2] or 'B0' in nm else 'o', color=col, ms=5 if 'HD' in nm or 'CMi' in nm else 6)
    label(ax2, ra - 0.3, dec + 0.25, f'{nm}' + (f' {dist}' if dist else ''), size=6.5, color=col)
ax2.set_xlim(124, 100)
ax2.set_ylim(4, 31)
ax2.set_aspect(1 / np.cos(np.deg2rad(17)))
ax2.set_xlabel('RA (deg)')
ax2.set_ylabel('Dec (deg)')
ax2.legend(loc='lower right', fontsize=6.5, ncol=1, framealpha=0.8)
ax2.grid(alpha=0.25)

# kinematics
ax3 = fig.add_axes([0.695, 0.08, 0.29, 0.25])
dg, vr = A['dgrid'], A['vrot']
ax3.fill_between(dg * 1000, vr - 7, vr + 7, color=FGDIM, alpha=0.2, lw=0)
ax3.plot(dg * 1000, vr, color=FG, lw=1.2, label='flat rotation ± 7 km/s')
for k, col, nm in (('[O III] DESI (+14.5 helio)', OIII_C, 'DESI [O III]'), ('[O III] SDSS sky fibres (+12.5)', HB_C, 'SDSS [O III]'),
                   ('[O III] LAMOST Wu+2026 (+23.8..+26.8)', DANGER, 'LAMOST [O III]'), ('HD 61717 Ca II (-16.88)', WARN, 'HD 61717 Ca II')):
    ax3.axhline(S['vlsr'][k], color=col, lw=1.1, ls='-')
    dy = {'SDSS [O III]': -4.5}.get(nm, 1.0)
    label(ax3, 2900, S['vlsr'][k] + dy, nm, size=7, color=col, ha='right')
ax3.set_xlim(0, 3000)
ax3.set_ylim(-35, 30)
ax3.set_xlabel('distance (pc)')
ax3.set_ylabel('V_LSR (km/s)')
ax3.grid(alpha=0.25)
ax3.legend(loc='upper left', fontsize=6.5)
finish(fig, ROOT / 'figures/verify_distance_constraints.png', tight=False)

# =============================================================== figure 2: HD 61717
fig, axs = new_sky_figure(15, 9.3, 'HD 61717 · the −17 km/s Ca II component',
                          'UVES 2023-11-12 · interstellar-like (b≈1.5 km/s, both K and H), SB2 star, not at either stellar velocity',
                          nrows=2, ncols=2)
ax = axs[0, 0]
vK, fK, mK, cK = P61['vK'], P61['fK'], P61['mK'], P61['cK']
vH, fH, mH, cH = P61['vH'], P61['fH'], P61['mH'], P61['cH']
ax.plot(vK, fK / cK, color=FG, lw=0.8, label='Ca II K')
ax.plot(vK, mK / cK, color=OIII_C, lw=1.0)
ax.plot(vH, fH / cH - 0.35, color=FG, lw=0.8, label='Ca II H − 0.35')
ax.plot(vH, mH / cH - 0.35, color=OIII_C, lw=1.0, label='3-component τ fit')
lism = {'LIC': 19.33, 'Gem': 33.26, 'Aur': 19.45, 'Hyades': 4.79, 'Blue': 10.15}
for nm, vv in lism.items():
    ax.axvline(vv, color=SUN, lw=0.5, alpha=0.6)
    label(ax, vv + 0.5, 0.42, nm, size=6.5, color=SUN, rotation=90)
c0 = H61['caII_fit']['components'][0]
ax.axvline(c0['v'], color=WARN, lw=0.8)
label(ax, -80, 0.52, f"v = {c0['v']:.2f}±{c0['e_v']:.2f}  b = {c0['b']:.2f}±{c0['e_b']:.2f} km/s\n"
      f"log N(Ca II) = {c0['logN']:.2f}±{c0['e_logN']:.2f}", size=8, color=WARN)
ax.set_xlim(-85, 60)
ax.set_ylim(0.3, 1.1)
ax.set_xlabel('V_helio (km/s)')
ax.set_ylabel('normalised flux')
ax.legend(loc='lower right', fontsize=7)
ax.set_title('Ca II K and H; LISM cloud velocities (Redfield & Linsky 2008 calculator)', fontsize=9)
ax = axs[0, 1]
v2, v1, r17 = P61['vD2'], P61['vD1'], P61['resid17']
ax.plot(v2, r17, color=FG, lw=0.8, label='Na I D2 (telluric + local removed)')
ax.plot(v1, r17 - 0.04, color=SII_C, lw=0.8, label='Na I D1 − 0.04')
ax.axvline(c0['v'], color=WARN, lw=0.8)
na = H61['naI']['per_line_ew_m17']
label(ax, -78, 0.925, f"W(D2) = {na['NaD2']['W_mA']:.1f}±{na['NaD2']['eW_mA']:.1f} mÅ, W(D1) = {na['NaD1']['W_mA']:.1f}±{na['NaD1']['eW_mA']:.1f} mÅ\n"
      f"N(Ca II)/N(Na I) = {H61['ratio_m17']['CaII_over_NaI']:.0f} (68 %: {H61['ratio_m17']['CaII_over_NaI_68'][0]:.0f}–{H61['ratio_m17']['CaII_over_NaI_68'][1]:.0f});"
      f" local gas {H61['ratio_local']['CaII_over_NaI']:.1f}", size=7.5, color=WARN)
ax.set_xlim(-80, 50)
ax.set_ylim(0.9, 1.03)
ax.set_xlabel('V_helio (km/s)')
ax.legend(loc='lower right', fontsize=7)
ax.set_title('Na I D after HARPS telluric template', fontsize=9)
ax = axs[1, 0]
VG = P61['VG']
for i, nm in enumerate(P61['ph_names']):
    off = 0.07 * i
    ax.plot(VG, P61['ph_f'][i] - off, color=FG, lw=0.5)
    ax.plot(VG, P61['ph_m1'][i] - off, color=FGDIM, lw=0.9, ls='--')
    ax.plot(VG, P61['ph_m2'][i] - off, color=HA_C, lw=1.0)
    label(ax, 330, 1.0 - off, str(nm), size=7)
st = H61['stellar']
for vv in (st['two_component']['rv1'], st['two_component']['rv2']):
    ax.axvline(vv, color=HA_C, lw=0.6)
ax.axvline(c0['v'], color=WARN, lw=0.8)
ax.axvline(st['apogee_visits']['vhelio_mean'], color=NII_C, lw=0.8)
label(ax, -440, 0.595, f"2 components: {st['two_component']['rv1']:.0f} and +{st['two_component']['rv2']:.0f} km/s, "
      f"v sin i ≈ {st['two_component']['vsini1']:.0f}/{st['two_component']['vsini2']:.0f}\n"
      f"ΔBIC(1−2) = {st['delta_BIC_1minus2']:.0f} · APOGEE 22 visits mean +{st['apogee_visits']['vhelio_mean']:.0f} "
      f"(σ {st['apogee_visits']['vhelio_std']:.1f}) · Gaia RUWE 2.14", size=7.5, color=HA_C)
ax.set_xlim(-450, 450)
ax.set_ylim(0.58, 1.04)
ax.set_xlabel('V_helio (km/s)')
ax.set_title('photospheric lines: 1-component (dashed) vs SB2 (red) rotational fits', fontsize=9)
ax = axs[1, 1]
names = ['Ca II K', 'Ca II H', 'Na I D2', 'Na I D1', 'Ti II 3384', 'Ca I 4227', 'CH+ 4232', 'Fe II 4233', 'Ti II 3761']
o = H61['other_lines']
sh = H61['shell_test_m17']
vals = [(o['CaK']['W_m17_mA'], o['CaK']['eW_m17_mA']), (o['CaH']['W_m17_mA'], o['CaH']['eW_m17_mA']),
        (na['NaD2']['W_mA'], na['NaD2']['eW_mA']), (na['NaD1']['W_mA'], na['NaD1']['eW_mA']),
        (o['TiII3384']['W_m17_mA'], o['TiII3384']['eW_m17_mA']), (o['CaI4227']['W_m17_mA'], o['CaI4227']['eW_m17_mA']),
        (o['CHp4232']['W_m17_mA'], o['CHp4232']['eW_m17_mA']), (sh['FeII4233']['W_m17_mA'], sh['FeII4233']['eW_mA']),
        (sh['TiII3761']['W_m17_mA'], sh['TiII3761']['eW_mA'])]
xs = np.arange(len(names))
ax.errorbar(xs, [v for v, _ in vals], yerr=[e for _, e in vals], fmt='o', color=WARN, ms=5, capsize=2)
ax.axhline(0, color=FGDIM, lw=0.6)
ax.set_xticks(xs)
ax.set_xticklabels(names, rotation=35, fontsize=7.5, ha='right')
ax.set_ylabel('EW at −17 km/s (mÅ)')
ax.set_yscale('symlog', linthresh=2)
ax.set_ylim(-3, 45)
u_ = H61['programme'] if 'programme' in H61 else None
label(ax, 3.6, 12, f"Ti II/Ca II < {H61['TiII_over_CaII']['m17_3sig_upper']} (3σ) vs {H61['TiII_over_CaII']['local']} in local gas\n"
      "no narrow −17 km/s line in 5 other stars of the programme", size=7.5, color=FG)
ax.set_title('matched-filter EWs at −17 km/s (shell-star lines absent)', fontsize=9)
finish(fig, ROOT / 'figures/verify_distance_hd61717.png')

# =============================================================== figure 3: IUE
fig, axs = new_sky_figure(15, 8.5, 'HD 60848 · IUE SWP high dispersion',
                          '13 co-added spectra · O8:V:pe at ~1.9 kpc, 2.3° W of the spine · NEWSIPS velocity zero point',
                          nrows=2, ncols=4)
VGI = PI['VG']
pan = [('SII1253', 'S II 1253'), ('SiII1304', 'Si II 1304'), ('AlIII1854', 'Al III 1854'), ('AlIII1862', 'Al III 1862'),
       ('SiIV1393', 'Si IV 1393'), ('SiIV1402', 'Si IV 1402'), ('CIV1548', 'C IV 1548'), ('NV1238', 'N V 1238')]
ref = IUE['refined_doublet_fits']['HD 60848']
for ax, (ln, tt) in zip(axs.flat, pan):
    f = PI[f'HD_60848|{ln}|f']
    e = PI[f'HD_60848|{ln}|e']
    ax.fill_between(VGI, f - e, f + e, color=FGDIM, alpha=0.25, lw=0)
    ax.plot(VGI, f, color=FG, lw=0.8)
    if 'HD_58050|' + ln + '|f' in PI:
        ax.plot(VGI, PI['HD_58050|' + ln + '|f'] - 0.0, color=SII_C, lw=0.6, alpha=0.6)
    ax.axhline(1, color=FGDIM, lw=0.5)
    ax.axvspan(5, 30, color=OIII_C, alpha=0.12, lw=0)
    ax.axvspan(-200, -70, color=DANGER, alpha=0.06, lw=0)
    ax.axvspan(70, 200, color=DANGER, alpha=0.06, lw=0)
    key = ln[:-4]
    core = IUE['stars']['HD 60848']['lines'][ln]['core_m50_p60']
    if key == 'NV':
        label(ax, -290, 0.18, f"log N < {core['logN_3sig_upper']:.2f} (3σ, AOD)", size=7, color=OIII_C)
    elif key == 'CIV':
        c50 = IUE['stars']['HD 60848']['lines']['CIV1550']['core_m50_p60']
        label(ax, -290, 0.12, f"joint fit log N = {ref[key]['logN']:.2f}±{ref[key]['e_logN']:.2f}\n"
              f"1550 AOD {c50['logN']:.2f}±{c50['elogN']:.2f}\n1548 wings = stellar wind", size=7, color=OIII_C)
    elif key in ref and isinstance(ref[key], dict):
        rr = ref[key]
        label(ax, -290, 0.18, f"log N = {rr['logN']:.2f}±{rr['e_logN']:.2f}\nv = {rr['v_helio']:+.1f}±{rr['e_v']:.1f}  b = {rr['b_kms']:.0f}",
              size=7, color=OIII_C)
    ax.set_xlim(-300, 300)
    ax.set_ylim(0.0, 1.45)
    ax.set_title(tt, fontsize=9)
    ax.set_xlabel('V_helio (km/s)')
axs[0, 0].set_ylabel('normalised flux')
axs[1, 0].set_ylabel('normalised flux')
label(axs[1, 3], -290, 1.3, 'blue band: Silk-like V (+5…+30)\nred bands: |v| 70–200 km/s\norange: HD 58050 (B2Ve, 452 pc)',
      size=7, color=FGDIM)
finish(fig, ROOT / 'figures/verify_distance_iue.png')
print('figures written')
