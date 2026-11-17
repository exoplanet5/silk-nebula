"""verify_xrayfuv: X-ray edge figure (eRASS1 0.2-0.6 keV map around Monogem with null placements;
cross-spine profiles; null distributions; radial scan)."""
import sys, json
import numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_style import *
from verify_xrayfuv_common import *
import verify_xrayfuv_xray_edge as X
from scipy import ndimage

res = json.loads((DER / 'verify_xrayfuv_xray_edge.json').read_text())
runs = {r['key']: r for r in res['runs']}
nul = np.load(DER / 'verify_xrayfuv_xray_nulls.npz')
ra, dec, _ = load_spine()
E = X.load_erass()

fig = plt.figure(figsize=(15, 10.5))
fig.text(0.012, 0.985, 'Silk Nebula · Monogem soft X-ray edge test (eRASS1 0.2–0.6 keV, ROSAT)', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
fig.text(0.012, 0.958, 'x < 0 = west / toward Monogem centre · nulls: rotations about PSR B0656+14 (same radius) and random sky translations', color=FGDIM, fontsize=9.5, va='top')
gs = fig.add_gridspec(2, 3, left=0.05, right=0.985, top=0.92, bottom=0.07, wspace=0.28, hspace=0.32)

# (a) map: bin eRASS soft (full mask) onto RA/Dec grid 0.25 deg, smooth
ax = fig.add_subplot(gs[0, 0])
m = E['mask_full']
rae, dee = np.arange(88, 128.01, 0.25), np.arange(-6, 36.01, 0.25)
S, _, _ = np.histogram2d(E['dec'][m], E['ra'][m], bins=[dee, rae], weights=E['soft'][m])
N, _, _ = np.histogram2d(E['dec'][m], E['ra'][m], bins=[dee, rae])
img = np.where(N > 3, S / np.maximum(N, 1), np.nan)
fill = np.where(np.isfinite(img), img, 0); w8 = np.isfinite(img).astype(float)
sm = ndimage.gaussian_filter(fill, 1.5) / np.maximum(ndimage.gaussian_filter(w8, 1.5), 1e-3)
sm[w8 == 0] = np.nan
im = ax.imshow(sm, origin='lower', extent=[rae[0], rae[-1], dee[0], dee[-1]], cmap='magma', vmin=2.5, vmax=8.5, aspect='auto')
ax.invert_xaxis()
rc, dc = CENTRES['PSR']
for th in [-120, -80, 60, 100, 140]:
    r2, d2 = rotate_about(ra, dec, rc, dc, th)
    ax.plot(r2, d2, color=FGDIM, lw=0.8, ls='-', alpha=0.8)
ax.plot(ra, dec, color=ACCENT, lw=1.6)
ax.plot(*CENTRES['PSR'], marker='*', color=SUN, ms=10, ls='none')
ax.plot(*CENTRES['K2'], marker='+', color=OK, ms=10, ls='none')
t = np.linspace(0, 360, 361)
cr, cd = offset_point(rc, dc, 11.36, t); ax.plot(cr, cd, color=SUN, lw=0.6, ls=':')
label(ax, rc - 1, dc - 1.5, 'PSR B0656+14', color=SUN, size=7.5)
label(ax, 117.5, 29, 'Silk spine', color=ACCENT, size=8)
label(ax, 124, -4.5, 'grey: example rotated nulls', color=FGDIM, size=7)
ax.set_xlabel('RA (deg)'); ax.set_ylabel('Dec (deg)')
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.01); cb.set_label('cts s$^{-1}$ deg$^{-2}$', fontsize=8)
ax.set_title('(a) eRASS1 0.2–0.6 keV (sources masked, 0.25° bins)', fontsize=9.5, loc='left')

# (b) profiles
ax = fig.add_subplot(gs[0, 1])
cols = {'all': FG, 'south': OIII_C, 'north': WARN}
for seg in ['all', 'south', 'north']:
    k = f'eRASS|soft|full|raw|PSR|{seg}'
    x, p = map(np.array, runs[k]['real_profile'])
    p = np.array([np.nan if v is None else v for v in p], float)
    ax.plot(x, p, color=cols[seg], lw=1.4, label=f"{seg}: W/E {runs[k]['real']['WE']:.2f}, step {100*runs[k]['real']['cfrac']:+.0f}%")
ax.axvline(0, color=ACCENT, lw=0.8, ls=':')
ax.set_xlabel('offset x from spine (deg; + = east/outside)'); ax.set_ylabel('cts s$^{-1}$ deg$^{-2}$')
ax.legend(loc='upper right', fontsize=7.5, ncol=1)
ax.set_title('(b) cross-spine profiles, eRASS1 0.2–0.6 keV', fontsize=9.5, loc='left')

# (c) null distributions of step fraction (south) and W/E (all)
ax = fig.add_subplot(gs[0, 2])
k = 'eRASS|soft|full|raw|PSR|south'
for nm, col, lab in [('rot', SUN, 'rotations'), ('tr', ACCENT, 'translations')]:
    v = nul[f'{k}|{nm}|cfrac']
    ax.hist(100 * v, bins=np.arange(-30, 31, 2.5), histtype='step', color=col, lw=1.4, label=f'{lab} (n={v.size})')
obs = 100 * runs[k]['real']['cfrac']
ax.axvline(obs, color=DANGER, lw=1.8)
zr = runs[k]['cfrac']['rot']; zt = runs[k]['cfrac']['tr']
label(ax, obs + 1, ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] > 0 else 5, f'Silk south\n{obs:+.1f}%\nrot p={zr["p_ge"]:.3f}\ntr p={zt["p_ge"]:.3f}', color=DANGER, size=7.5)
ax.set_xlabel('step beyond linear gradient at x=0 (% of level)'); ax.set_ylabel('placements')
ax.legend(loc='upper left', fontsize=7.5)
ax.set_title('(c) local step, south (Dec < 17.6): nulls', fontsize=9.5, loc='left')

# (d) W/E null distributions (all)
ax = fig.add_subplot(gs[1, 0])
k = 'eRASS|soft|full|raw|PSR|all'
for nm, col, lab in [('rot', SUN, 'rotations about PSR'), ('tr', ACCENT, 'sky translations')]:
    v = nul[f'{k}|{nm}|WE']
    ax.hist(v, bins=np.arange(0.6, 1.6, 0.025), histtype='step', color=col, lw=1.4, label=f'{lab} (n={v.size})')
obs = runs[k]['real']['WE']
ax.axvline(obs, color=DANGER, lw=1.8)
label(ax, 1.42, 18, f"W/E {obs:.2f}\nrot p={runs[k]['WE']['rot']['p_ge']:.2f}\ntr p={runs[k]['WE']['tr']['p_ge']:.3f}", color=DANGER, size=7.5)
ax.set_xlabel('W/E ratio (−2.5…−0.5 vs +0.5…+2.5 deg)'); ax.set_ylabel('placements')
ax.legend(loc='upper left', fontsize=7.5)
ax.set_title('(d) wave-1 W/E statistic: nulls', fontsize=9.5, loc='left')

# (e) radial scan
ax = fig.add_subplot(gs[1, 1])
for k, col, lab in [('eRASS|soft|full|raw|PSR|south', OIII_C, 'eRASS south'), ('eRASS|soft|full|raw|PSR|north', WARN, 'eRASS north'),
                    ('ROSAT|r12|none|raw|PSR|south', OK, 'ROSAT R1+R2 south')]:
    sc = runs[k]['scan']
    ax.plot([s['dr'] for s in sc], [100 * s['cfrac'] for s in sc], color=col, lw=1.4, marker='o', ms=3, label=lab)
ax.axvline(0, color=FGDIM, lw=0.8, ls=':'); ax.axhline(0, color=FGDIM, lw=0.6)
ax.set_xlabel('trace moved radially from PSR by Δr (deg)'); ax.set_ylabel('step at trace (% of level)')
ax.legend(loc='upper right', fontsize=7.5)
ax.set_title('(e) radial scan of the local step', fontsize=9.5, loc='left')

# (f) Monogem radial profiles around PSR: Silk sectors vs others
ax = fig.add_subplot(gs[1, 2])
sm_ = res['shell_model_PSR']
P = np.array(json.loads(json.dumps(res['shell_model_PSR'].get('P', [])))) if 'P' in sm_ else None
rr, pa = sep_pa(rc, dc, E['ra'].astype(float), E['dec'].astype(float))
srad, spa = sep_pa(rc, dc, ra, dec)
ok = E['mask_full']
redg = np.arange(0, 20.01, 0.5); rcen = 0.5 * (redg[1:] + redg[:-1])
def prof(sel):
    return np.array([np.mean(E['soft'][sel & (rr >= redg[i]) & (rr < redg[i + 1])]) if np.sum(sel & (rr >= redg[i]) & (rr < redg[i + 1])) > 30 else np.nan for i in range(len(rcen))])
for (p0, p1), col in [((75, 100), OIII_C), ((30, 60), WARN), ((100, 160), OK), ((160, 220), '#c792ea'), ((300, 360), FGDIM)]:
    ax.plot(rcen, prof(ok & (pa >= p0) & (pa < p1)), color=col, lw=1.3, label=f'PA {p0}–{p1}°')
ax.axvspan(srad[dec < 17.6].min(), srad[dec < 17.6].max(), color=OIII_C, alpha=0.18)
ax.axvspan(srad[dec >= 17.6].min(), srad[dec >= 17.6].max(), color=WARN, alpha=0.10)
ax.set_xlabel('distance from PSR B0656+14 (deg)'); ax.set_ylabel('cts s$^{-1}$ deg$^{-2}$')
ax.legend(loc='lower left', fontsize=7, ncol=2)
label(ax, 11.3, 8.6, 'Silk S', color=OIII_C, size=7.5); label(ax, 14.3, 8.6, 'Silk N', color=WARN, size=7.5)
ax.set_title('(f) Monogem radial profiles by PA (Silk: PA 30–100°)', fontsize=9.5, loc='left')
finish(fig, str(FIG / 'verify_xrayfuv_xray_edge.png'), tight=False)
print('ok')
