"""Extract vector markers (fibers, [OIII] triangles, V_LSR points) from Wu et al. 2026 (arXiv:2610.13329v1)
Figures 2, 4 and 5 using PyMuPDF, and register them to RA/Dec from the axis tick marks.

Outputs (data/derived/):
  geometry_fig2_fibers_pagecoords.csv, geometry_fig4_triangles_pagecoords.csv, geometry_fig5_points_pagecoords.csv
  geometry_pdf_axes.json  (tick-based axis calibrations)
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_extract_pdf_vectors.py
"""
import json
import numpy as np
import fitz
import matplotlib
matplotlib.use('Agg')
from matplotlib import colormaps

ROOT = '/Users/mickey/astronomy/silk-nebula'
PDF = f'{ROOT}/sources/arxiv/2610.13329v1.pdf'
doc = fitz.open(PDF)


def lin(p1, v1, p2, v2):
    """value(page) = a + b*page from two (page, value) pairs."""
    b = (v2 - v1) / (p2 - p1)
    return v1 - b * p1, b


def fit_ticks(pages, values):
    A = np.vstack([np.ones(len(pages)), pages]).T
    coef, res, *_ = np.linalg.lstsq(A, values, rcond=None)
    resid = values - A @ coef
    return coef, float(np.max(np.abs(resid)))


axes = {}
# ---------------- Fig 5 (page 11): RA/Dec in degrees, labelled ticks -----------------
# x ticks labelled 119..113 (left->right), y ticks 12..20 (bottom->top), colorbar 0..14 by 2 (bottom->top)
f5x = np.array([94.73, 130.32, 165.9, 201.49, 237.08, 272.66, 308.25]); f5xv = np.array([119, 118, 117, 116, 115, 114, 113.])
f5y = np.array([391.26, 319.92, 248.59, 177.26, 105.93]); f5yv = np.array([12, 14, 16, 18, 20.])
f5c = np.array([415.76, 369.88, 324.01, 278.13, 232.26, 186.38, 140.5, 94.63]); f5cv = np.arange(0, 15, 2.)
(ax5, rx5) = fit_ticks(f5x, f5xv); (ay5, ry5) = fit_ticks(f5y, f5yv); (ac5, rc5) = fit_ticks(f5c, f5cv)
axes['fig5'] = {'ra_coef': ax5.tolist(), 'dec_coef': ay5.tolist(), 'cbar_coef': ac5.tolist(),
                'tick_resid_deg': [rx5, ry5], 'axis_box_page': [76.36, 71.69, 308.25, 415.76],
                'ra_range': [float(ax5[0] + ax5[1] * 76.36), float(ax5[0] + ax5[1] * 308.25)],
                'dec_range': [float(ay5[0] + ay5[1] * 415.76), float(ay5[0] + ay5[1] * 71.69)],
                'cbar_range': [float(ac5[0] + ac5[1] * 415.76), float(ac5[0] + ac5[1] * 71.69)]}
# ---------------- Fig 2 (page 5): Dec ticks labelled; RA ticks labelled in truncated h m ---------------
f2y = np.array([634.97, 508.36, 381.74, 255.13, 128.52]); f2yv = np.array([12, 14, 16, 18, 20.])
(ay2, ry2) = fit_ticks(f2y, f2yv)
axes['fig2'] = {'dec_coef': ay2.tolist(), 'tick_resid_deg': [None, ry2],
                'ra_ticks_page': [108.4, 190.85, 273.3, 355.75, 438.19, 520.64],
                'ra_tick_labels': ['7h58m', '7h52m', '7h47m', '7h42m', '7h37m', '7h32m'],
                'axis_box_page': [108.4, 67.75, 520.64, 678.46]}
# ---------------- Fig 4 (page 9) ----------------
f4y = np.array([652.24, 569.01, 485.78, 402.55, 319.32, 236.09, 152.86, 69.64]); f4yv = np.arange(13, 21, 1.)
(ay4, ry4) = fit_ticks(f4y, f4yv)
f4c = np.array([655.49, 568.48, 481.48, 394.47, 307.46, 220.45, 133.44]); f4cv = np.arange(10, 80, 10.)
(ac4, rc4) = fit_ticks(f4c, f4cv)
axes['fig4'] = {'dec_coef': ay4.tolist(), 'cbar_coef': ac4.tolist(), 'tick_resid_deg': [None, ry4],
                'ra_ticks_page': [112.78, 200.27, 287.75, 375.23, 462.72],
                'ra_tick_labels': ['7h55m', '7h51m', '7h46m', '7h42m', '7h38m'],
                'axis_box_page': [112.78, 69.64, 462.72, 668.88],
                'cbar_box_page': [471.6, 69.6, 501.1, 668.9],
                'cbar_range': [float(ac4[0] + ac4[1] * 668.88), float(ac4[0] + ac4[1] * 69.64)]}

jet = colormaps['jet'].resampled(256)
lut = jet(np.arange(256))[:, :3]


def inv_jet(rgb, vmin, vmax):
    d = np.sum((lut - np.asarray(rgb)[None, :]) ** 2, axis=1)
    i = int(np.argmin(d))
    v = vmin + (i + 0.5) / 256 * (vmax - vmin)
    return v, i, float(np.sqrt(d[i]))


# Fig 2 fibers: red stroked circles, 8 bezier segments
rows2 = []
for x in doc[4].get_drawings():
    if x['type'] == 's' and len(x['items']) == 8 and x['color'] and x['color'][0] > 0.9 and x['color'][1] < 0.1:
        r = x['rect']
        rows2.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2, r.width, r.height))
rows2 = np.array(rows2)
# Fig 4 triangles
rows4 = []
for x in doc[8].get_drawings():
    if x['type'] == 'fs' and len(x['items']) == 3:
        pts = np.array([[it[1].x, it[1].y] for it in x['items']])
        ys = pts[:, 1]
        xc = pts[np.argmin(ys), 0]               # apex x
        yc = 0.5 * (ys.min() + ys.max())        # matplotlib '^' marker centre = midpoint of apex & base
        v, i, dist = inv_jet(x['fill'], *axes['fig4']['cbar_range'])
        rows4.append((xc, yc, ys.max() - ys.min(), *x['fill'], i, v, dist))
rows4 = np.array(rows4)
# Fig 5 points: filled 8-segment circles inside the map axes
rows5 = []
bx = axes['fig5']['axis_box_page']
for x in doc[10].get_drawings():
    # two scatter layers: 'f' (fill-only path, stroke drawn separately) and 'fs' (fill+stroke in one path)
    if x['type'] in ('f', 'fs') and len(x['items']) == 8 and x['fill'] is not None:
        r = x['rect']
        xc, yc = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        if bx[0] <= xc <= bx[2] and bx[1] <= yc <= bx[3]:
            v, i, dist = inv_jet(x['fill'], *axes['fig5']['cbar_range'])
            rows5.append((xc, yc, r.width, *x['fill'], i, v, dist, x['seqno'], 1 if x['type'] == 'fs' else 0))
rows5 = np.array(rows5)

np.savetxt(f'{ROOT}/data/derived/geometry_fig2_fibers_pagecoords.csv', rows2, delimiter=',', fmt='%.4f',
           header='x_pt,y_pt,w_pt,h_pt  # Fig 2 (page 5) red fiber circles, PDF page points', comments='')
np.savetxt(f'{ROOT}/data/derived/geometry_fig4_triangles_pagecoords.csv', rows4, delimiter=',', fmt='%.5f',
           header='x_pt,y_pt,h_pt,r,g,b,jet_index,intensity_from_jet,color_dist  # Fig 4 (page 9) [OIII] triangles', comments='')
np.savetxt(f'{ROOT}/data/derived/geometry_fig5_points_pagecoords.csv', rows5, delimiter=',', fmt='%.5f',
           header='x_pt,y_pt,w_pt,r,g,b,jet_index,vlsr_from_jet,color_dist,seqno,layer_fs  # Fig 5 (page 11) Halpha V_LSR points', comments='')
json.dump(axes, open(f'{ROOT}/data/derived/geometry_pdf_axes.json', 'w'), indent=1)
print('fig2 fibers', len(rows2), 'fig4 tri', len(rows4), 'fig5 pts', len(rows5))
print(json.dumps(axes, indent=1))
print('max jet colour mismatch fig4 %.4f fig5 %.4f' % (rows4[:, -1].max(), rows5[:, 8].max()))
