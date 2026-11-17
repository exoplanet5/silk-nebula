"""Register Fig 2 and Fig 4 RA axes (whose tick labels are truncated h-m strings) by cross-matching their
vector markers to the Fig 5 points (axes labelled in decimal degrees), then write the digitized LAMOST
[OIII] fiber table data/derived/geometry_lamost_oiii_points_digitized.csv.
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_register_figs.py
"""
import json
import numpy as np
from scipy.spatial import cKDTree

ROOT = '/Users/mickey/astronomy/silk-nebula'
D = f'{ROOT}/data/derived'
ax = json.load(open(f'{D}/geometry_pdf_axes.json'))
f2 = np.loadtxt(f'{D}/geometry_fig2_fibers_pagecoords.csv', delimiter=',', skiprows=1)
f4 = np.loadtxt(f'{D}/geometry_fig4_triangles_pagecoords.csv', delimiter=',', skiprows=1)
f5 = np.loadtxt(f'{D}/geometry_fig5_points_pagecoords.csv', delimiter=',', skiprows=1)

ra5 = ax['fig5']['ra_coef'][0] + ax['fig5']['ra_coef'][1] * f5[:, 0]
de5 = ax['fig5']['dec_coef'][0] + ax['fig5']['dec_coef'][1] * f5[:, 1]


def register(xp, yp, ref_ra, ref_de, a0, b0, dec_coef, tol=0.03, fit_dec=True, niter=8):
    """Fit RA = a + b*x and (optionally) Dec = c + d*y by iterative nearest-neighbour matching."""
    c, d = dec_coef
    a, b = a0, b0
    tree = cKDTree(np.c_[ref_ra, ref_de])
    for it in range(niter):
        ra = a + b * xp; de = c + d * yp
        dist, idx = tree.query(np.c_[ra, de])
        m = dist < (tol if it > 1 else 0.1)
        A = np.c_[np.ones(m.sum()), xp[m]]
        (a, b), *_ = np.linalg.lstsq(A, ref_ra[idx[m]], rcond=None)
        if fit_dec:
            A = np.c_[np.ones(m.sum()), yp[m]]
            (c, d), *_ = np.linalg.lstsq(A, ref_de[idx[m]], rcond=None)
    ra = a + b * xp; de = c + d * yp
    dist, idx = tree.query(np.c_[ra, de])
    m = dist < tol
    rr = (ra[m] - ref_ra[idx[m]]) * np.cos(np.radians(de[m])) * 3600
    rd = (de[m] - ref_de[idx[m]]) * 3600
    return dict(a=float(a), b=float(b), c=float(c), d=float(d), n_match=int(m.sum()), n=len(xp),
                rms_ra_arcsec=float(np.sqrt(np.mean(rr ** 2))), rms_dec_arcsec=float(np.sqrt(np.mean(rd ** 2))),
                median_sep_arcsec=float(np.median(dist[m]) * 3600)), ra, de, dist, idx


# --- Fig 2: initial guess assumes equal aspect with left edge 119.53 (from truncated labels)
b0 = -abs(ax['fig2']['dec_coef'][1]); a0 = 119.53 - b0 * 108.4
r2, ra2, de2, dist2, idx2 = register(f2[:, 0], f2[:, 1], ra5, de5, a0, b0, ax['fig2']['dec_coef'])
print('Fig2 registration', r2)
_d = dist2 * 3600
r2['frac_within_2arcsec'] = float(np.mean(_d[_d < 108] < 2)); r2['n_within_2arcsec'] = int((_d < 2).sum())
print('Fig2 fibres matched within 2":', r2['n_within_2arcsec'], r2['frac_within_2arcsec'])
# check tick labels consistency
for xt, lab in zip(ax['fig2']['ra_ticks_page'], ax['fig2']['ra_tick_labels']):
    ra = r2['a'] + r2['b'] * xt; h = ra / 15; mm = (h - int(h)) * 60
    print(f'  fig2 tick {lab}: RA {ra:.3f} deg = {int(h)}h{mm:05.2f}m')
# --- Fig 4: register triangles against Fig 2 fibres (complete fibre set) in data coordinates
b0 = -abs(ax['fig4']['dec_coef'][1]); a0 = 118.83 - b0 * 112.78
r4, ra4, de4, dist4, idx4 = register(f4[:, 0], f4[:, 1], ra2, de2, a0, b0, ax['fig4']['dec_coef'], tol=0.02)
print('Fig4 registration (vs Fig2 fibres)', r4)
for xt, lab in zip(ax['fig4']['ra_ticks_page'], ax['fig4']['ra_tick_labels']):
    ra = r4['a'] + r4['b'] * xt; h = ra / 15; mm = (h - int(h)) * 60
    print(f'  fig4 tick {lab}: RA {ra:.3f} deg = {int(h)}h{mm:05.2f}m')
# Fig 5 points vs Fig 2 fibres (consistency)
t2 = cKDTree(np.c_[ra2, de2]); d52, i52 = t2.query(np.c_[ra5, de5])
print('Fig5 points matched to Fig2 fibres within 0.02 deg: %d/%d, median sep %.1f"' %
      ((d52 < 0.02).sum(), len(ra5), np.median(d52[d52 < 0.02]) * 3600))

# V_LSR (Halpha) at each [OIII] fibre: Fig 5 point nearest to the matched fibre position
t5 = cKDTree(np.c_[ra5, de5]); d45, i45 = t5.query(np.c_[ra4, de4])
reg = {'fig2': r2, 'fig4': r4,
       'fig5_vs_fig2': {'n_match': int((d52 < 0.02).sum()), 'n': len(ra5),
                        'median_sep_arcsec': float(np.median(d52[d52 < 0.02]) * 3600)}}
json.dump(reg, open(f'{D}/geometry_fig_registration.json', 'w'), indent=1)

hdr = ('# DIGITIZED from Wu et al. 2026 (arXiv:2610.13329v1) Fig 4 (vector triangles, colour -> relative [OIII] '
       'intensity via jet LUT, colorbar %.2f-%.2f) and Fig 5 (Halpha V_LSR colour, colorbar 0-15 km/s).\n'
       '# NOT original measurements. RA axis of Fig 4 registered by matching triangles to Fig 2 fibres '
       '(Fig 2 registered to Fig 5 degree axes): rms %.1f" (RA) %.1f" (Dec), %d/%d matched.\n'
       '# Colour quantisation: intensity +-%.2f, V_LSR +-%.3f km/s (half a jet LUT step). halpha_vlsr is blank when no '
       'Fig 5 point lies within 0.02 deg; halpha_vlsr_censored=-1/+1 marks colours at the colorbar ends (true value <=0 or >=15 km/s).\n' % (*ax['fig4']['cbar_range'], r4['rms_ra_arcsec'], r4['rms_dec_arcsec'],
                                               r4['n_match'], r4['n'], (ax['fig4']['cbar_range'][1] - ax['fig4']['cbar_range'][0]) / 512, 15 / 512))
lines = ['id,ra_deg,dec_deg,oiii_rel_intensity,sep_to_fig2_fibre_arcsec,fibre_ra_deg,fibre_dec_deg,halpha_vlsr_kms,halpha_vlsr_censored,sep_to_fig5_point_arcsec']
for k in range(len(f4)):
    fra, fde = ra2[idx4[k]], de2[idx4[k]]
    v = f5[i45[k], 7] if d45[k] < 0.02 else np.nan
    cens = '' if np.isnan(v) else ('-1' if f5[i45[k], 6] == 0 else ('1' if f5[i45[k], 6] == 255 else '0'))
    lines.append(f'{k},{ra4[k]:.4f},{de4[k]:.4f},{f4[k, 7]:.2f},{dist4[k] * 3600:.1f},{fra:.4f},{fde:.4f},'
                 f'{"" if np.isnan(v) else f"{v:.2f}"},{cens},{d45[k] * 3600:.1f}')
open(f'{D}/geometry_lamost_oiii_points_digitized.csv', 'w').write(hdr + '\n'.join(lines) + '\n')
# full fibre and V_LSR tables in data coords
np.savetxt(f'{D}/geometry_lamost_fibres_fig2_digitized.csv', np.c_[ra2, de2], delimiter=',', fmt='%.4f',
           header='ra_deg,dec_deg  # DIGITIZED Fig 2 fibre positions (2684 usable spectra), Wu+2026', comments='')
cens5 = np.where(f5[:, 6] == 0, -1, np.where(f5[:, 6] == 255, 1, 0))
np.savetxt(f'{D}/geometry_lamost_halpha_vlsr_fig5_digitized.csv', np.c_[ra5, de5, f5[:, 7], f5[:, 10], cens5], delimiter=',', fmt='%.4f',
           header='ra_deg,dec_deg,halpha_vlsr_kms,layer2,censored  # censored=-1/+1: at colorbar end (<=0 / >=15 km/s); layer2=1: second scatter layer (drawn last). DIGITIZED Fig 5 (jet colour -> V_LSR, 0-15 km/s, clipped at ends), Wu+2026', comments='')
print('n oiii', len(f4), 'with vlsr', np.isfinite([float(l.split(',')[7]) if l.split(',')[7] else np.nan for l in lines[1:]]).sum())
