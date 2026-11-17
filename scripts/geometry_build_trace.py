"""Combine the per-source strand fits into the Silk Nebula spine, fit circles (small circles on the sphere),
derive the Halpha patch polygon and the [OIII] envelope, and write data/derived/silk_trace.json
(preliminary=false) plus data/derived/geometry_spine_points.csv.
Run: ~/.venvs/astro313/bin/python -I scripts/geometry_build_trace.py
"""
import csv, json, sys
import numpy as np
from scipy import ndimage
from scipy.optimize import least_squares
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import REGIONS, gal
from geometry_imreg import PolyMap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/Users/mickey/astronomy/silk-nebula'
D = f'{ROOT}/data/derived'
rng = np.random.default_rng(20261116)
HIP37797 = (116.2524833, 15.8950042)     # SIMBAD (CfA mirror) ICRS J2000, Gaia DR3
REG_FLOOR = 0.005                        # deg: registration + resampling floor (star-fit rms 3-19", 0.02 deg grid)
SINGLE_SRC_SCATTER = 0.02                # deg: assumed source-to-source scatter where only one image covers

fits = list(csv.DictReader(open(f'{D}/geometry_strand_fits.csv')))
for q in fits:
    for k in ('dec', 'ra', 'ra_err', 'amp_snr', 'fwhm_sky_deg'):
        q[k] = float(q[k])
    q['ok'] = (q['amp_snr'] > 5) and (q['at_bound'] == 'False') and np.isfinite(q['ra_err']) and q['ra_err'] < 0.1


def combine(strand, step=0.25):
    q = [x for x in fits if x['strand'] == strand and x['ok']]
    dmin, dmax = min(x['dec'] for x in q), max(x['dec'] for x in q)
    out = []
    for d0 in np.arange(np.floor(dmin * 4) / 4, dmax + 1e-6, step):
        b = [x for x in q if d0 - step / 2 <= x['dec'] < d0 + step / 2]
        if not b:
            continue
        srcs = sorted(set(x['source'] for x in b))
        sm = {s: np.median([x['ra'] for x in b if x['source'] == s]) for s in srcs}
        serr = {s: np.median([x['ra_err'] for x in b if x['source'] == s]) for s in srcs}
        fw = np.median([x['fwhm_sky_deg'] for x in b])
        snr = np.median([x['amp_snr'] for x in b])
        vals = np.array(list(sm.values()))
        ra_c = float(np.mean(vals))
        cosd = np.cos(np.radians(d0))
        # conservative per-point uncertainty (great-circle deg, ~perpendicular to the N-S filament):
        # full source-to-source scatter (not divided by sqrt(n): differences are systematic), median Gaussian-fit
        # error, and a registration/resampling floor
        scat = float(np.std(vals, ddof=1)) * cosd if len(vals) > 1 else SINGLE_SRC_SCATTER
        fiterr = float(np.median(list(serr.values()))) * cosd
        sig_sky = float(np.sqrt(scat ** 2 + fiterr ** 2 + REG_FLOOR ** 2))
        out.append(dict(ra=round(ra_c, 4), dec=round(float(d0), 3), sigma_deg=round(sig_sky, 4),
                        hwhm_deg=round(float(fw) / 2, 3), n_sources=len(vals), sources='+'.join(srcs),
                        source_scatter_deg=round(scat, 4) if len(vals) > 1 else None, median_snr=round(float(snr), 1),
                        per_source={s: round(float(v), 4) for s, v in sm.items()}))
    return out


strands = {s: combine(s) for s in ('E', 'W', 'W2')}

# ---------- trim ends where only marginal detections remain (require snr>8 in >=1 source, contiguous)
for s in strands:
    pts = strands[s]
    good = [p['median_snr'] >= 8 for p in pts]
    if any(good):
        i0 = good.index(True); i1 = len(good) - 1 - good[::-1].index(True)
        strands[s] = pts[i0:i1 + 1]
# W fades into the Halpha patch: drop its northern points where the four images disagree by >0.1 deg
while strands['W'] and (strands['W'][-1]['source_scatter_deg'] or 0) > 0.1:
    strands['W'].pop()
for s in strands:
    for p in strands[s]:
        p['flag'] = ''
        if s == 'E' and 17.6 <= p['dec'] <= 18.8:
            p['flag'] = 'junction: E ridge shifts ~0.2 deg east between Dec 17.5 and 17.75 and is faint (SNR 6-12) to Dec ~19'
        if p['n_sources'] == 1:
            p['flag'] = (p['flag'] + '; ' if p['flag'] else '') + 'single source (NSNS wide)'


# ---------- small-circle fit on the sphere
def unit(ra, de):
    ra, de = np.radians(ra), np.radians(de)
    return np.c_[np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)]


def circ_fit(ra, de, sig, x0=(80.0, 10.0, 35.0)):
    P = unit(ra, de)

    def res(q):
        c = unit(q[0], q[1])[0]
        ang = np.degrees(np.arccos(np.clip(P @ c, -1, 1)))
        return (ang - q[2]) / sig
    r = least_squares(res, x0, x_scale=[5, 5, 5])
    J = r.jac
    chi2 = float(np.sum(r.fun ** 2)); dof = max(1, len(ra) - 3)
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, chi2 / dof)
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(3, np.nan)
    return r.x, err, chi2, dof, r.fun * sig


def circ_with_boot(pts, nboot=400, label=''):
    ra = np.array([p['ra'] for p in pts]); de = np.array([p['dec'] for p in pts]); sg = np.array([p['sigma_deg'] for p in pts])
    # initial guess: tangent-plane circle through three points
    best = None
    for x0 in [(80, 10, 35), (95, 15, 20), (60, 0, 55), (100, 18, 15), (140, 25, 25)]:
        f = circ_fit(ra, de, sg, x0)
        if best is None or f[2] < best[2]:
            best = f
    p, e, chi2, dof, resid = best
    if p[2] > 90:                                    # same small circle described about the antipode
        p = np.array([(p[0] + 180) % 360, -p[1], 180 - p[2]])
    bs = []
    for _ in range(nboot):
        k = rng.integers(0, len(ra), len(ra))
        try:
            bp = circ_fit(ra[k], de[k], sg[k], tuple(p))[0]
            if bp[2] > 90:
                bp = np.array([(bp[0] + 180) % 360, -bp[1], 180 - bp[2]])
            bp[0] = p[0] + ((bp[0] - p[0] + 180) % 360 - 180)   # unwrap RA near the best value
            bs.append(bp)
        except Exception:
            pass
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [16, 84], axis=0)
    cl, cb = gal(p[0] % 360, p[1])
    return dict(label=label, n_points=len(ra), dec_range=[float(de.min()), float(de.max())],
                center_ra=round(float(p[0] % 360), 2), center_dec=round(float(p[1]), 2), radius_deg=round(float(p[2]), 2),
                center_l=round(float(cl[0]), 2), center_b=round(float(cb[0]), 2),
                err_formal=[round(float(x), 2) for x in e],
                boot_16_84={'center_ra': [round(float(lo[0] % 360), 2), round(float(hi[0] % 360), 2)],
                            'center_dec': [round(float(lo[1]), 2), round(float(hi[1]), 2)],
                            'radius_deg': [round(float(lo[2]), 2), round(float(hi[2]), 2)]},
                chi2=round(chi2, 1), dof=dof, rms_resid_deg=round(float(np.sqrt(np.mean(resid ** 2))), 4))


circles = {
    'E_full': circ_with_boot(strands['E'], label='E (main) strand, full length'),
    'E_south_le17p5': circ_with_boot([p for p in strands['E'] if p['dec'] <= 17.5], label='E strand south of the junction (Dec <= 17.5)'),
    'E_north_ge17p75': circ_with_boot([p for p in strands['E'] if p['dec'] >= 17.75], label='E strand north of the junction (Dec >= 17.75)'),
    'W': circ_with_boot(strands['W'], label='W (bright, ends at Halpha patch) strand'),
    'E_full_plus_W': circ_with_boot(strands['E'] + strands['W'], label='E and W strands together'),
}
# radius consistency of candidate centres (literature positions) with the E spine
cands = {'Geminga (now)': (98.4756, 17.7703), 'Geminga birthplace (Pellizza+2019 via NSNS)': (84.7088, 7.6303),
         'PSR B0656+14': (104.9506, 14.2393), 'Monogem single-SNR centre (Knies+2024 via NSNS)': (107.08, 13.50),
         'Monogem two-SNR centre (Knies+2024 via NSNS)': (104.15, 15.30), 'G205.6+12.4 centre (Knies+2024 via NSNS)': (110.33, 12.94),
         'NSNS C0 (sphere centre if arc on limb)': (73.1139, 6.5749)}
PE = unit([p['ra'] for p in strands['E']], [p['dec'] for p in strands['E']])
cand_out = {}
for k, (cr, cd) in cands.items():
    a = np.degrees(np.arccos(np.clip(PE @ unit(cr, cd)[0], -1, 1)))
    cand_out[k] = {'ra': cr, 'dec': cd, 'mean_dist_deg': round(float(a.mean()), 2), 'rms_about_mean_deg': round(float(a.std()), 3),
                   'range_deg': [round(float(a.min()), 2), round(float(a.max()), 2)]}
bc = unit(circles['E_full']['center_ra'], circles['E_full']['center_dec'])[0]
for k in cand_out:
    cand_out[k]['sep_from_fitted_centre_deg'] = round(float(np.degrees(np.arccos(np.clip(unit(*cands[k])[0] @ bc, -1, 1)))), 2)
# straight-line (great-circle) fit for comparison: chi2 of best great circle (radius fixed 90)


def gc_fit(pts):
    ra = np.array([p['ra'] for p in pts]); de = np.array([p['dec'] for p in pts]); sg = np.array([p['sigma_deg'] for p in pts])
    P = unit(ra, de)

    def res(q):
        c = unit(q[0], q[1])[0]
        return (np.degrees(np.arccos(np.clip(P @ c, -1, 1))) - 90) / sg
    best = min((least_squares(res, x0) for x0 in [(20, 0), (200, 0), (30, 30)]), key=lambda r: r.cost)
    return float(2 * best.cost), len(ra) - 2


gchi, gdof = gc_fit(strands['E'])
circles['E_full']['great_circle_chi2_dof'] = [round(gchi, 1), gdof]

# ---------- [OIII] envelope polygon (east edge / west edge at HWHM of outermost strands)
env = []
for d0 in np.arange(12.0, 29.0, 0.25):
    edges = []
    for s in ('E', 'W', 'W2'):
        for p in strands[s]:
            if abs(p['dec'] - d0) < 1e-6:
                cosd = np.cos(np.radians(d0))
                edges += [p['ra'] + p['hwhm_deg'] / cosd, p['ra'] - p['hwhm_deg'] / cosd]
    if edges:
        env.append((d0, max(edges), min(edges)))
env = np.array(env)
envelope_polygon = [[round(float(r), 3), round(float(d), 3)] for d, r, _ in env] + \
                   [[round(float(r), 3), round(float(d), 3)] for d, _, r in env[::-1]]

# ---------- Halpha patch polygon from NSNS narrow Halpha (+ SWONS red as check)
z = np.load(f'{D}/geometry_maps.npz')
ra_g, de_g = z['ra_grid'], z['dec_grid']
RA, DE = np.meshgrid(ra_g, de_g)


def patch_from(m, sig_px=3, frac=0.5, box_lim=(114.9, 116.4, 19.2, 22.4)):
    """Half-maximum region of the smoothed Halpha map around the bright Halpha/[SII] complex at the
    north end of the W strand: threshold = bg + frac*(peak - bg), bg = median over RA 114.5-118, Dec 17-24."""
    box = (RA > box_lim[0]) & (RA < box_lim[1]) & (DE > box_lim[2]) & (DE < box_lim[3])
    ref = (RA > 114.5) & (RA < 118.0) & (DE > 17.0) & (DE < 24.0)
    mm = np.where(np.isfinite(m), m, np.nanmedian(m[ref]))
    sm = ndimage.gaussian_filter(mm, sig_px)
    bg = float(np.median(sm[ref])); peak = float(np.max(sm[box]))
    thr = bg + frac * (peak - bg)
    lab, n = ndimage.label((sm > thr) & box)
    j = np.unravel_index(np.argmax(np.where(box, sm, -np.inf)), sm.shape)
    comp = ndimage.binary_fill_holes(lab == lab[j])
    return comp, dict(bg=bg, peak=peak, thr=thr, frac=frac, smoothing_sigma_deg=sig_px * 0.02, box=list(box_lim))


def polygon_from_mask(comp):
    fig, ax = plt.subplots()
    cs = ax.contour(ra_g, de_g, comp.astype(float), levels=[0.5])
    segs = cs.allsegs[0]
    plt.close(fig)
    seg = max(segs, key=len)
    # thin to ~40 vertices
    k = max(1, len(seg) // 40)
    return [[round(float(a), 3), round(float(b), 3)] for a, b in seg[::k]]


def moments(comp):
    w = comp.astype(float)
    cosd = np.cos(np.radians(np.average(DE, weights=w)))
    x = (RA - np.average(RA, weights=w)) * cosd; y = DE - np.average(DE, weights=w)
    cxx, cyy, cxy = np.average(x * x, weights=w), np.average(y * y, weights=w), np.average(x * y, weights=w)
    ev, evec = np.linalg.eigh(np.array([[cxx, cxy], [cxy, cyy]]))
    major = 4 * np.sqrt(ev[1]); minor = 4 * np.sqrt(ev[0])          # full extent ~ 4 sigma (uniform ellipse)
    pa = np.degrees(np.arctan2(evec[0, 1], evec[1, 1])) % 180       # east of north
    area = w.sum() * 0.02 * 0.02 * cosd
    return dict(centroid_ra=round(float(np.average(RA, weights=w)), 3), centroid_dec=round(float(np.average(DE, weights=w)), 3),
                major_deg=round(float(major), 2), minor_deg=round(float(minor), 2), pa_deg=round(float(pa), 1),
                area_deg2=round(float(area), 3))


ha_n, ha_n_info = patch_from(z['nsnsNs_Ha'])
ha_n30, ha_n30_info = patch_from(z['nsnsNs_Ha'], frac=0.3)
ha_sii, _ = patch_from(z['nsnsNs_SII'])
ha_s, ha_s_info = patch_from(z['swons_R'] - 0.3 * z['swons_B'])
ha_poly = polygon_from_mask(ha_n)
ha_poly30 = polygon_from_mask(ha_n30)
ha_mom = moments(ha_n); ha_mom30 = moments(ha_n30); ha_mom_s = moments(ha_s); ha_mom_sii = moments(ha_sii)
overlap = float((ha_n & ha_s).sum() / max(1, (ha_n | ha_s).sum()))
overlap_sii = float((ha_n & ha_sii).sum() / max(1, (ha_n | ha_sii).sum()))


# ---------- LAMOST [OIII] fibres (digitized Fig 4) vs strands
lam = np.genfromtxt(f'{D}/geometry_lamost_oiii_points_digitized.csv', delimiter=',', skip_header=3, names=True)


def strand_ra(s, d):
    pts = strands[s]
    dd = np.array([p['dec'] for p in pts]); rr = np.array([p['ra'] for p in pts])
    return np.where((d >= dd.min()) & (d <= dd.max()), np.interp(d, dd, rr), np.nan)


off_E = (lam['ra_deg'] - strand_ra('E', lam['dec_deg'])) * np.cos(np.radians(lam['dec_deg']))
off_W = (lam['ra_deg'] - strand_ra('W', lam['dec_deg'])) * np.cos(np.radians(lam['dec_deg']))
near = np.where(np.isnan(off_W) | (np.abs(off_E) <= np.abs(np.nan_to_num(off_W, nan=9))), 'E', 'W')
dist_near = np.where(near == 'E', off_E, off_W)
lam_stats = dict(n=int(len(lam)), dec_range=[float(lam['dec_deg'].min()), float(lam['dec_deg'].max())],
                 frac_within_0p15deg_of_a_strand=round(float(np.mean(np.abs(dist_near) < 0.15)), 3),
                 median_abs_offset_deg=round(float(np.median(np.abs(dist_near))), 3),
                 n_nearest_E=int((near == 'E').sum()), n_nearest_W=int((near == 'W').sum()))
# intensity-weighted ridge in 0.5 deg Dec bins
lam_ridge = []
for d0 in np.arange(13.5, 18.6, 0.5):
    k = (lam['dec_deg'] >= d0 - 0.25) & (lam['dec_deg'] < d0 + 0.25)
    if k.sum() >= 4:
        w = lam['oiii_rel_intensity'][k]
        lam_ridge.append(dict(dec=float(d0), ra_wmean=round(float(np.average(lam['ra_deg'][k], weights=w)), 3),
                              ra_rms=round(float(np.sqrt(np.average((lam['ra_deg'][k] - np.average(lam['ra_deg'][k], weights=w)) ** 2, weights=w))), 3),
                              n=int(k.sum()), E_ra=float(np.round(strand_ra('E', d0), 3)),
                              W_ra=None if np.isnan(strand_ra('W', d0)) else float(np.round(strand_ra('W', d0), 3))))

# ---------- regions: Table 1 plus positions of the regions as drawn on Fig 1b (registered)
drawn = json.load(open(f'{D}/geometry_fig1b_drawn_regions.json'))
regions = {}
for k, (r, d) in REGIONS.items():
    l, b = gal(r, d)
    regions[k] = dict(ra=r, dec=d, l=round(float(l[0]), 3), b=round(float(b[0]), 3), source='Wu+2026 Table 1',
                      drawn_in_fig1b=drawn.get(k))
    for s in ('E', 'W'):
        x = strand_ra(s, d)
        regions[k][f'offset_from_{s}_deg'] = None if np.isnan(x) else round(float((r - x) * np.cos(np.radians(d))), 3)

# ---------- assemble
spine = []
for p in strands['E']:
    l, b = gal(p['ra'], p['dec'])
    spine.append(dict(ra=p['ra'], dec=p['dec'], sigma_deg=p['sigma_deg'], source=p['sources'], hwhm_deg=p['hwhm_deg'],
                      n_sources=p['n_sources'], l=round(float(l[0]), 3), b=round(float(b[0]), 3), flag=p['flag']))
E_len = float(np.sum(np.degrees(np.arccos(np.clip(np.sum(unit([p['ra'] for p in spine[:-1]], [p['dec'] for p in spine[:-1]]) *
                                                         unit([p['ra'] for p in spine[1:]], [p['dec'] for p in spine[1:]]), 1), -1, 1)))))
W_len = float(np.sum(np.degrees(np.arccos(np.clip(np.sum(unit([p['ra'] for p in strands['W'][:-1]], [p['dec'] for p in strands['W'][:-1]]) *
                                                         unit([p['ra'] for p in strands['W'][1:]], [p['dec'] for p in strands['W'][1:]]), 1), -1, 1)))))
hw_E = np.array([p['hwhm_deg'] for p in strands['E']])
sep_EW = []
for p in strands['W']:
    x = strand_ra('E', p['dec'])
    if np.isfinite(x):
        sep_EW.append(((x - p['ra']) * np.cos(np.radians(p['dec'])), p['dec']))
sep_EW = np.array(sep_EW)
out = {
    'version': '1.0',
    'preliminary': False,
    'date': '2026-11-16',
    'author_tag': 'geometry',
    'spine': spine,
    'half_width_deg': round(float(np.median(hw_E)), 3),
    'half_width_note': 'HWHM of the Gaussian cross-section of the main (E) strand, median over the spine; '
                       f'16-84% range {np.percentile(hw_E, 16):.3f}-{np.percentile(hw_E, 84):.3f} deg. The full [OIII] '
                       'envelope (E+W strands) is wider: see envelope_half_width_deg and envelope_polygon.',
    'envelope_half_width_deg': {'median': round(float(np.median((env[:, 1] - env[:, 2]) * np.cos(np.radians(env[:, 0])) / 2)), 3),
                                'max': round(float(np.max((env[:, 1] - env[:, 2]) * np.cos(np.radians(env[:, 0])) / 2)), 3)},
    'envelope_polygon': envelope_polygon,
    'strands': {
        'E': {'description': 'main, longest [OIII] strand (= spine); merged single ridge south of Dec ~15',
              'points': strands['E'], 'length_deg': round(E_len, 2)},
        'W': {'description': 'bright western [OIII] strand, separates from E near Dec 15, bends west and ends at the Halpha patch',
              'points': strands['W'], 'length_deg': round(W_len, 2),
              'separation_from_E_deg': {'min': round(float(sep_EW[:, 0].min()), 3), 'max': round(float(sep_EW[:, 0].max()), 3)}},
        'W2': {'description': 'faint far-western [OIII] strand (NSNS wide view), Dec ~21.5-24', 'points': strands['W2']},
    },
    'halpha_patch': ha_poly,
    'halpha_patch_info': {'method': 'half-maximum contour (bg + 0.5*(peak-bg)) of the NSNS narrow Halpha channel (stars opened, '
                                    '0.06 deg Gaussian smoothing) in box RA 114.9-116.4, Dec 19.2-22.4; connected component containing the peak',
                          'threshold': ha_n_info, 'moments_halfmax': ha_mom,
                          'outer_polygon_30pct': ha_poly30, 'moments_30pct': ha_mom30,
                          'check_nsns_SII_halfmax_moments': ha_mom_sii, 'iou_Ha_vs_SII': round(overlap_sii, 2),
                          'check_swons_red_halfmax_moments': ha_mom_s, 'iou_nsns_vs_swons': round(overlap, 2),
                          'note': 'major/minor = 4 x sqrt(eigenvalues) of the pixel distribution (equals full axes for a uniform ellipse); '
                                  'SWONS red is truncated by the image edge at RA ~114.7 and includes broadband red continuum/halo'},
    'lamost_field': {'ra': round(HIP37797[0], 5), 'dec': round(HIP37797[1], 5), 'radius_deg': 2.5,
                     'center_star': 'HIP 37797 = HD 62523 (SIMBAD via simbad.cfa.harvard.edu mirror, 2026-11-16)',
                     'fig1b_drawn_circle': drawn.get('circle')},
    'regions': regions,
    'curvature': {'center_ra': circles['E_full']['center_ra'], 'center_dec': circles['E_full']['center_dec'],
                  'radius_deg': circles['E_full']['radius_deg'],
                  'center_l': circles['E_full']['center_l'], 'center_b': circles['E_full']['center_b'],
                  'center_ra_16_84': circles['E_full']['boot_16_84']['center_ra'],
                  'center_dec_16_84': circles['E_full']['boot_16_84']['center_dec'],
                  'radius_16_84': circles['E_full']['boot_16_84']['radius_deg'],
                  'formal_err_ra_dec_radius': circles['E_full']['err_formal'],
                  'method': 'weighted least-squares small-circle fit on the sphere to the E-strand spine points '
                            '(weights 1/sigma^2); uncertainties from 400 bootstrap resamples (16-84%)',
                  'fits': circles, 'candidate_centres_vs_E_spine': cand_out,
                  'caveat': 'a single circle is a poor description (reduced chi2 >> 1): the arc is nearly straight south of '
                            'Dec ~17.5 and curves strongly to the north-west further north (S-shaped with a junction/kink near Dec 17.5-18)'},
    'lamost_oiii_fibres_vs_strands': {'stats': lam_stats, 'intensity_weighted_ridge': lam_ridge,
                                      'file': 'data/derived/geometry_lamost_oiii_points_digitized.csv'},
    'literature_geometry': {
        'NSNS_O200.0': {'ra': 115.4516, 'dec': 19.9427, 'size_arcmin': [1015.267, 181.582], 'pa_deg': 173,
                        'desc': 'Large ~17 x 2.5 deg arc with dominant [OIII] (Ziegenbalg, simg.de/nebulae3/gem-01.html)'},
        'NSNS_O200.0a': {'ra': 115.6653, 'dec': 21.0673, 'size_arcmin': [183.859, 120.001], 'pa_deg': 173,
                         'desc': 'part of O200.0 only visible in Halpha and [SII]'},
        'NSNS_O200.0b': {'ra': 117.7208, 'dec': 12.5230, 'size_arcmin': [181.103, 43.902], 'pa_deg': 175,
                         'desc': 'very faint [OIII] filament that appears to belong to O200.0'},
        'NSNS_C0_sphere_center_if_on_limb': {'ra': 73.1139, 'dec': 6.5749},
        'paper_center': {'ra': 116.68, 'dec': 16.42},
    },
    'notes': ('Spine = ridge of the main (E) [OIII] strand measured independently in 4 registered images: '
              'SWONS colour composite (Wu+2026 Fig 1, 227 Tycho-2 stars, median resid 2.8"), Wu+2026 Fig 4 [OIII] '
              'background (144 stars, rms 19"), NSNS narrow 20"/px and NSNS wide 40"/px (published STG WCS, '
              'checked with 330/152 Tycho-2 stars, median resid 3.6"/6.5"). Per-point sigma_deg combines source-to-source '
              'scatter, Gaussian-fit error and a 0.005 deg registration floor (single-source points assume 0.02 deg scatter). '
              'RA/Dec J2000 (ICRS). See notes/geometry.md.'),
}
json.dump(out, open(f'{D}/silk_trace.json', 'w'), indent=1)
with open(f'{D}/geometry_spine_points.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['strand', 'ra_deg', 'dec_deg', 'l_deg', 'b_deg', 'sigma_deg', 'hwhm_deg', 'n_sources', 'sources', 'source_scatter_deg', 'median_snr', 'flag'])
    for s, pts in strands.items():
        for p in pts:
            l, b = gal(p['ra'], p['dec'])
            w.writerow([s, p['ra'], p['dec'], round(float(l[0]), 4), round(float(b[0]), 4), p['sigma_deg'], p['hwhm_deg'],
                        p['n_sources'], p['sources'], p['source_scatter_deg'], p['median_snr'], p['flag']])
np.save(f'{D}/geometry_halpha_mask_nsns.npy', ha_n)
print(json.dumps({k: out[k] for k in ('half_width_deg', 'envelope_half_width_deg')}, indent=1))
print('E len', E_len, 'W len', W_len, 'n E', len(strands['E']), 'n W', len(strands['W']), 'n W2', len(strands['W2']))
print('E dec range', strands['E'][0]['dec'], strands['E'][-1]['dec'], 'W', strands['W'][0]['dec'], strands['W'][-1]['dec'])
print('sigma E: median %.4f max %.4f' % (np.median([p['sigma_deg'] for p in strands['E']]), max(p['sigma_deg'] for p in strands['E'])))
print(json.dumps(circles, indent=1))
print('Ha', ha_mom, ha_mom30, 'SII', ha_mom_sii, overlap_sii, 'SWONS', ha_mom_s, 'IoU', overlap, ha_n_info)
print(json.dumps(cand_out, indent=1))
print('LAMOST', lam_stats); print(lam_ridge)
print(json.dumps(regions, indent=1))
