#!/usr/bin/env python
"""spectra_ism: Faraday rotation (Taylor et al. 2009 NVSS RMs) and pulsar DM/RM around the Silk Nebula.

Inputs (VizieR CfA mirror, fetched by spectra_ism_fetch_vizier.py and a B/psr re-query with RM columns):
  data/raw/spectra_ism_vizier/taylor2009_nvss_rm_cone.tsv   (J/ApJ/702/1230; 20-deg cone on 116.0,+20.0)
  data/raw/spectra_ism_vizier/psrcat_bpsr_cone35_rm.tsv     (B/psr, ATNF catalogue snapshot on VizieR)
Outputs: data/derived/spectra_ism_rm_results.json, spectra_ism_rm_sources.csv, figures/spectra_ism_rm_psr.png
"""
import csv
import json
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, FIG, filament_coords, draw_trace, gal  # noqa: E402
from silk_style import new_sky_figure, finish, label, FG, FGDIM, OIII_C, HA_C, WARN, ACCENT  # noqa: E402

RNG = np.random.default_rng(7)


def read_viz(path):
    lines = [ln for ln in open(path).read().splitlines() if ln and not ln.startswith('#')]
    hdr = lines[0].split('\t')
    return [dict(zip(hdr, [x.strip() for x in ln.split('\t')])) for ln in lines[3:]]


T = read_viz(RAW / 'spectra_ism_vizier/taylor2009_nvss_rm_cone.tsv')
c = SkyCoord([r['RAJ2000'] for r in T], [r['DEJ2000'] for r in T], unit=(u.hourangle, u.deg))
ra, dec = c.ra.deg, c.dec.deg
rm = np.array([float(r['RM']) for r in T])
erm = np.array([float(r['e_RM']) for r in T])
dperp, dsp, beyond = filament_coords(ra, dec)
l, b = gal(ra, dec)
print('Taylor+09 RMs in 20-deg cone:', len(rm))

# smooth foreground: 2-D quadratic in (x, y) about the filament centre, robust (iterative 3-sigma clip),
# fitted to sources more than 1 deg from the spine
x = (ra - 116.5) * np.cos(np.deg2rad(dec))
y = dec - 20.0
A = np.c_[np.ones_like(x), x, y, x * x, x * y, y * y]
use = (np.abs(dperp) > 1.0) & (np.hypot(x, y) < 15)
for _ in range(5):
    coef, *_ = np.linalg.lstsq(A[use], rm[use], rcond=None)
    res = rm - A @ coef
    s = 1.4826 * np.median(np.abs(res[use] - np.median(res[use])))
    use = use & (np.abs(res) < 3 * s)
resid = rm - A @ coef
print('foreground fit: coef', coef.round(2), 'robust sigma', round(s, 1))

alongs = (~beyond) & (dsp > 12) & (dsp < 28.2)
out = {'n_sources_cone20': int(len(rm)), 'foreground_quadratic_coef': coef.tolist(), 'robust_sigma_resid': float(s),
       'median_RM_cone': float(np.median(rm)), 'median_eRM': float(np.median(erm))}
edges = np.arange(-8, 8.01, 1.0)
prof = {'edges': edges.tolist(), 'median_resid': [], 'err': [], 'n': [], 'median_RM': []}
for b1, b2 in zip(edges[:-1], edges[1:]):
    m = alongs & (dperp >= b1) & (dperp < b2)
    if m.sum() >= 3:
        boots = [np.median(resid[m][RNG.integers(0, m.sum(), m.sum())]) for _ in range(2000)]
        prof['median_resid'].append(float(np.median(resid[m])))
        prof['err'].append(float(np.std(boots)))
        prof['median_RM'].append(float(np.median(rm[m])))
    else:
        prof['median_resid'].append(np.nan)
        prof['err'].append(np.nan)
        prof['median_RM'].append(np.nan)
    prof['n'].append(int(m.sum()))
out['cross_profile'] = prof


def jump(w_in, w_out, key):
    e = alongs & (dperp > w_in) & (dperp < w_out)
    w = alongs & (dperp < -w_in) & (dperp > -w_out)
    d = np.median(resid[e]) - np.median(resid[w])
    bs = [np.median(resid[e][RNG.integers(0, e.sum(), e.sum())]) - np.median(resid[w][RNG.integers(0, w.sum(), w.sum())])
          for _ in range(4000)]
    dr = np.median(rm[e]) - np.median(rm[w])
    out[key] = {'east_minus_west_median_resid': float(d), 'err_boot': float(np.std(bs)), 'n_east': int(e.sum()),
                'n_west': int(w.sum()), 'east_minus_west_median_RM_raw': float(dr)}


jump(0.0, 2.0, 'jump_0_2deg')
jump(0.0, 4.0, 'jump_0_4deg')
jump(0.5, 3.0, 'jump_0p5_3deg')
on = alongs & (np.abs(dperp) < 0.5)
off = alongs & (np.abs(dperp) > 1.0) & (np.abs(dperp) < 4.0)
bs = [np.median(resid[on][RNG.integers(0, on.sum(), on.sum())]) - np.median(resid[off][RNG.integers(0, off.sum(), off.sum())])
      for _ in range(4000)]
out['on_minus_off'] = {'on_deg': 0.5, 'off_deg': [1, 4], 'n_on': int(on.sum()), 'n_off': int(off.sum()),
                       'median_resid_on': float(np.median(resid[on])), 'median_resid_off': float(np.median(resid[off])),
                       'diff': float(np.median(resid[on]) - np.median(resid[off])), 'err_boot': float(np.std(bs)),
                       'on_sources': [{'ra': round(float(ra[i]), 3), 'dec': round(float(dec[i]), 3), 'RM': float(rm[i]),
                                       'eRM': float(erm[i]), 'dperp': round(float(dperp[i]), 2)} for i in np.where(on)[0]]}

# North / south halves: raw-RM profiles and an L1 step-function fit for the location of any RM step
halves = {}
for tag, (d1, d2) in {'south_dec12_20': (12, 20), 'north_dec20_28': (20, 28.2)}.items():
    mh = (~beyond) & (dsp >= d1) & (dsp < d2) & (np.abs(dperp) < 6)
    hp = {'bins': [], 'median_RM': [], 'err': [], 'n': [], 'frac_neg': []}
    for b1, b2 in ((-6, -4), (-4, -2), (-2, -1), (-1, 0), (0, 1), (1, 2), (2, 4), (4, 6)):
        mm = mh & (dperp >= b1) & (dperp < b2)
        bsm = [np.median(rm[mm][RNG.integers(0, mm.sum(), mm.sum())]) for _ in range(1000)] if mm.sum() > 1 else [np.nan]
        hp['bins'].append([b1, b2])
        hp['median_RM'].append(float(np.median(rm[mm])) if mm.sum() else np.nan)
        hp['err'].append(float(np.std(bsm)))
        hp['n'].append(int(mm.sum()))
        hp['frac_neg'].append(float(np.mean(rm[mm] < 0)) if mm.sum() else np.nan)

    def l1_step(dd, rr):
        grid = np.arange(-3.0, 3.01, 0.05)
        cost = []
        for d0 in grid:
            a, bb = rr[dd < d0], rr[dd >= d0]
            if len(a) < 5 or len(bb) < 5:
                cost.append(np.inf)
                continue
            cost.append(np.sum(np.abs(a - np.median(a))) + np.sum(np.abs(bb - np.median(bb))))
        k = int(np.argmin(cost))
        d0 = grid[k]
        return d0, float(np.median(rr[dd >= d0]) - np.median(rr[dd < d0]))
    d0, step = l1_step(dperp[mh], rm[mh])
    bs = []
    ii_all = np.where(mh)[0]
    for _ in range(500):
        ii = ii_all[RNG.integers(0, len(ii_all), len(ii_all))]
        bs.append(l1_step(dperp[ii], rm[ii]))
    bs = np.array(bs)
    e = mh & (dperp > 2) & (dperp < 6)
    w = mh & (dperp < -2) & (dperp > -6)
    bsd = [np.median(rm[e][RNG.integers(0, e.sum(), e.sum())]) - np.median(rm[w][RNG.integers(0, w.sum(), w.sum())]) for _ in range(2000)]
    halves[tag] = {'profile': hp, 'step_location_deg': float(d0), 'step_location_16_84': np.percentile(bs[:, 0], [16, 84]).tolist(),
                   'step_size': step, 'step_size_16_84': np.percentile(bs[:, 1], [16, 84]).tolist(),
                   'E2_6_minus_W2_6_median_RM': float(np.median(rm[e]) - np.median(rm[w])), 'err_boot': float(np.std(bsd)),
                   'frac_neg_W2_6': float(np.mean(rm[w] < 0)), 'frac_neg_E2_6': float(np.mean(rm[e] < 0)),
                   'n_W': int(w.sum()), 'n_E': int(e.sum())}
    print(tag, {k: v for k, v in halves[tag].items() if k != 'profile'})
out['halves'] = halves

# Smoothed (median within 1.5 deg) RM map on a grid, for the figure
gra, gde = np.meshgrid(np.arange(100, 132.01, 0.5), np.arange(4, 36.01, 0.5))
gsm = np.full(gra.shape, np.nan)
cc = SkyCoord(ra * u.deg, dec * u.deg)
for i in range(gra.shape[0]):
    gg = SkyCoord(gra[i] * u.deg, gde[i] * u.deg)
    for j in range(gra.shape[1]):
        mm = gg[j].separation(cc).deg < 1.5
        if mm.sum() >= 4:
            gsm[i, j] = np.median(rm[mm])

# Pulsars
P = read_viz(RAW / 'spectra_ism_vizier/psrcat_bpsr_cone35_rm.tsv')
psr = []
for r in P:
    try:
        dm = float(r['DM'])
    except ValueError:
        continue
    pra, pde = float(r['RA2000']), float(r['DE2000'])
    dp, ds, bey = filament_coords(pra, pde)
    plx = float(r['Plx']) if r['Plx'] else np.nan
    eplx = float(r['e_Plx']) if r['e_Plx'] else np.nan
    rmv = float(r['RM']) if r['RM'] else np.nan
    ermv = float(r['e_RM']) if r['e_RM'] else np.nan
    d_ind = 1.0 / plx if (np.isfinite(plx) and plx / eplx > 3) else np.nan     # kpc, only >3-sigma parallaxes
    psr.append({'psrj': r['PSRJ'], 'name': r['Name'], 'ra': pra, 'dec': pde, 'l': float(r['GLON']), 'b': float(r['GLAT']),
                'r_deg': float(r['_r']), 'dperp': float(dp[0]), 'DM': dm, 'RM': rmv, 'eRM': ermv, 'plx': plx, 'eplx': eplx,
                'd_plx_kpc': d_ind, 'd_cat_kpc': float(r['Dist']) if r['Dist'] else np.nan,
                'ne_mean': dm / (d_ind * 1000) if np.isfinite(d_ind) else np.nan,
                'Bpar_uG': 1.232 * rmv / dm if np.isfinite(rmv) else np.nan})
psr.sort(key=lambda p: p['r_deg'])
out['pulsars_within_20deg'] = [p for p in psr if p['r_deg'] < 20]
out['pulsars_with_parallax'] = [p for p in psr if np.isfinite(p['d_plx_kpc'])]

json.dump(out, open(DER / 'spectra_ism_rm_results.json', 'w'), indent=1, default=lambda v: None if v != v else float(v))
with open(DER / 'spectra_ism_rm_sources.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['ra', 'dec', 'l', 'b', 'RM', 'eRM', 'dperp_deg', 'dec_spine', 'beyond', 'resid_RM'])
    for i in range(len(rm)):
        w.writerow([round(ra[i], 4), round(dec[i], 4), round(float(l[i]), 3), round(float(b[i]), 3), rm[i], erm[i],
                    round(float(dperp[i]), 3), round(float(dsp[i]), 2), int(beyond[i]), round(float(resid[i]), 1)])
for k in ('jump_0_2deg', 'jump_0_4deg', 'jump_0p5_3deg'):
    print(k, out[k])
print('on-off', {k: v for k, v in out['on_minus_off'].items() if k != 'on_sources'})
print('profile', np.round(prof['median_resid'], 1), prof['n'])

# ------------------------------------------------------------------------------------------
fig, axs = new_sky_figure(19, 7.2, 'Silk Nebula · Faraday rotation and pulsars',
                          'NVSS RMs (Taylor+2009); ATNF pulsars (B/psr); red = RM > 0, blue = RM < 0', nrows=1, ncols=3,
                          gridspec_kw={'width_ratios': [1.0, 1.0, 0.95]})
ax = axs[0]
m = np.hypot(x, y) < 16
sz = np.clip(np.abs(rm[m]), 2, 150) * 0.9
ax.scatter(ra[m][rm[m] > 0], dec[m][rm[m] > 0], s=sz[rm[m] > 0], facecolors='none', edgecolors=HA_C, lw=0.7, alpha=0.85)
ax.scatter(ra[m][rm[m] <= 0], dec[m][rm[m] <= 0], s=sz[rm[m] <= 0], facecolors='none', edgecolors=ACCENT, lw=0.7, alpha=0.85)
for v in (10, 30, 100):
    ax.scatter([], [], s=v * 0.9, facecolors='none', edgecolors=FGDIM, label=f'|RM| = {v} rad m⁻²')
draw_trace(ax, color='#ffffff', lw=1.2, label=False)
for p in psr:
    if p['r_deg'] < 16:
        ax.plot(p['ra'], p['dec'], '*', color=WARN, ms=9 if np.isfinite(p['RM']) else 6, mec='k', mew=0.4)
        if p['r_deg'] < 13:
            txt = f"{p['name']}\nDM {p['DM']:.1f}" + (f", RM {p['RM']:+.1f}" if np.isfinite(p['RM']) else '')
            label(ax, p['ra'] - 0.4, p['dec'] - 1.6 if p['name'] == 'J0751+1807' else p['dec'] + 0.3, txt, size=6.8, color=WARN)
ax.set_xlim(132, 100)
ax.set_ylim(4, 36)
ax.set_aspect(1 / np.cos(np.deg2rad(20)))
ax.set_xlabel('RA (deg)')
ax.set_ylabel('Dec (deg)')
ax.legend(loc='lower left', fontsize=6.5, ncol=3)
ax = axs[1]
import matplotlib.colors as mcolors  # noqa: E402
im = ax.pcolormesh(gra, gde, gsm, cmap='RdBu_r', norm=mcolors.TwoSlopeNorm(0, -40, 80), shading='auto')
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.set_label('median RM within 1.5° (rad m⁻²)', color=FGDIM)
cb.ax.tick_params(colors=FGDIM, labelsize=7)
draw_trace(ax, color='#111111', lw=1.6, label=False, regions=False)
draw_trace(ax, color='#ffffff', lw=0.8, label=False)
ax.set_xlim(132, 100)
ax.set_ylim(4, 36)
ax.set_aspect(1 / np.cos(np.deg2rad(20)))
ax.set_xlabel('RA (deg)')
ax.set_title('Smoothed RM', loc='left', color=FG, fontsize=9)
ax = axs[2]
for tag, cc_, lab in (('north_dec20_28', HA_C, 'north (spine Dec 20–28)'), ('south_dec12_20', OIII_C, 'south (spine Dec 12–20)')):
    hp = halves[tag]['profile']
    xb = [0.5 * (a + b_) for a, b_ in hp['bins']]
    ax.errorbar(xb, hp['median_RM'], yerr=hp['err'], fmt='o-', color=cc_, ms=4, lw=1.2, label=lab, capsize=0)
ax.axhline(0, color=FGDIM, lw=0.6)
ax.axvspan(-0.22, 0.22, color='#ffffff', alpha=0.10, lw=0)
ax.set_xlim(-6.5, 6.5)
ax.set_ylim(-40, 70)
ax.set_xlabel('offset from spine (deg, + = east)')
ax.set_ylabel('median RM (rad m⁻²)')
hn = halves['north_dec20_28']
label(ax, 0.02, 0.04, f"north: E(2–6°) − W(2–6°) = {hn['E2_6_minus_W2_6_median_RM']:+.0f} ± {hn['err_boot']:.0f} rad m⁻²\n"
      f"RM<0: {100*hn['frac_neg_W2_6']:.0f}% W vs {100*hn['frac_neg_E2_6']:.0f}% E; step at {hn['step_location_deg']:+.1f}° "
      f"({hn['step_location_16_84'][0]:+.1f}…{hn['step_location_16_84'][1]:+.1f})", transform=ax.transAxes, size=7.2, color=FG)
ax.legend(loc='upper left', fontsize=7)
ax.set_title('Raw RM profiles, north vs south', loc='left', color=FG, fontsize=9)
finish(fig, FIG / 'spectra_ism_rm_psr.png')
