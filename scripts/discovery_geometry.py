"""Discovery-history geometry: quoted extents, positions and neighbours of the Silk Nebula.

Collects every published/posted geometric description of the filament (SWONS/Wu et al.,
QHYCCD press note, Ziegenbalg NSNS 2024 and 2025 catalogue entries, Fesen et al. 2026,
AstroBin plate solutions), converts them to common quantities (ends of the major axis in
RA/Dec and l/b, length, physical length at several distances), measures separations of the
neighbouring SNRs/candidates named by Fesen et al. (2026) from the NSNS 2025 major axis,
and finds where the major axis crosses the eROSITA Monogem circle(s).

Also converts the quoted [O III] surface brightnesses (mag/arcsec^2) to Rayleigh under stated
assumptions (AB zero point, all band flux = 5007 line, effective bandwidth 5-6 nm).

Inputs are literature/web numbers typed in below with their source URL (see notes/discovery.md).
Output: data/derived/discovery_extents.json
Run:    ~/.venvs/astro313/bin/python -I scripts/discovery_geometry.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import CENTER, REGIONS, MONOGEM_PULSAR, GEMINGA, load_trace  # noqa: E402

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
OUT = ROOT / 'data/derived/discovery_extents.json'


def hms(h, m, s):
    return 15.0 * (h + m / 60 + s / 3600)


def dms(d, m, s, sign=1):
    return sign * (abs(d) + m / 60 + s / 3600)


# ---- Ellipse descriptions published by Ziegenbalg (centre RA/Dec deg, major/minor arcmin, PA deg E of N)
# Source of 2024 values: GitHub stzie/simg.de-discoveries commit 2e1a2ac0 (2024-06-28), new/k.md;
# 2025 values: commit f8273e1d (2025-06-11), new2/o.md; identical numbers are encoded in the
# viewer links on https://simg.de/nebulae3/gem-01.html (accessed 2026-11-16).
NSNS = {
    'K1a_2024': dict(ra=115.9862, dec=18.6391, maj=857.354, mnr=146.401, pa=175,
                     label='NSNS K1a (2024-06-28): 14 x 2.5 deg'),
    'O200.0_2025': dict(ra=115.4516, dec=19.9427, maj=1015.267, mnr=181.582, pa=173,
                        label='NSNS O200.0 (2025-06-11): 17 x 2.5(3) deg'),
    'O200.0a_Ha_SII': dict(ra=115.6653, dec=21.0673, maj=183.859, mnr=120.001, pa=173,
                           label='NSNS O200.0a / K1b: Halpha+[SII]-only part'),
    'O200.0b_faint_south': dict(ra=117.7208, dec=12.5230, maj=181.103, mnr=43.902, pa=175,
                                label='NSNS O200.0b: very faint [OIII] filament (2025-05-22)'),
    'O205.0_east': dict(ra=118.8253, dec=16.2484, maj=309.932, mnr=99.866, pa=167,
                        label='NSNS O205.0 (ex K1c): faint filaments, maybe Monogem'),
}

# ---- Other quoted extents (no PA published)
QUOTED = [
    dict(source='SWONS AstroBin zh47r9 (2024-01-14)', length_deg=10.0, width_deg=1.0,
         url='https://www.astrobin.com/zh47r9/'),
    dict(source='QHYCCD press note (2024-07-12)', length_deg=15.0, width_deg=1.5,
         url='https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/'),
    dict(source='Ziegenbalg CN post / NSNS page v1 (2024-06-11/28)', length_deg=14.0, width_deg=2.5,
         url='https://www.cloudynights.com/forums/topic/925409-possible-discovery-of-a-large-emission-nebula-in-gemini-and-a-tiny-one/'),
    dict(source='Ziegenbalg 2025 RNAAS 9, 227 (2025-08-28)', length_deg=17.0, width_deg=3.0,
         url='https://doi.org/10.3847/2515-5172/adfec7'),
    dict(source='NSNS gem-01 page (current)', length_deg=17.0, width_deg=2.5,
         url='https://simg.de/nebulae3/gem-01.html'),
    dict(source='Sean Jackson AstroBin zwfplw (2025-05-02)', length_deg=15.0, width_deg=None,
         url='https://www.astrobin.com/zwfplw/'),
    dict(source='Fesen et al. 2026 arXiv:2608.26666 (Sect. 6.2.1)', length_deg=15.0, width_deg=None,
         url='https://arxiv.org/abs/2608.26666'),
    dict(source='Wu et al. 2026 AJ 172, 298 / arXiv:2610.13329', length_deg=10.0, width_deg=1.0,
         url='https://doi.org/10.3847/1538-3881/aeb11c'),
]

# ---- AstroBin plate solutions (AstroBin API v2, accessed 2026-11-16)
ASTROBIN_SOLUTIONS = {
    'zwfplw': dict(ra=116.019, dec=16.550, radius_deg=9.181, pixscale=3.939, orientation=266.792,
                   field_deg=(18.93, 3.27), note='Sean Jackson / Hunting Horn, RASA 11 + QHY600, OIII 576x300 s'),
}

# ---- Neighbouring SNRs / candidates (Fesen et al. 2026 Table 1 centres, J2000)
NEIGHBOURS = {
    'G206.6+6.1': dict(ra=hms(7, 2, 0), dec=dms(8, 23, 0), size=(3.6, 3.6), kind='known SNR (radio G206.7+5.9)'),
    'G190.4+12.5': dict(ra=hms(6, 57, 0), dec=dms(25, 34, 0), size=(6.1, 6.1), kind='X-ray SNR candidate (eROSITA)'),
    'G191.4+11.1': dict(ra=hms(6, 53, 0), dec=dms(24, 12, 0), size=(3.2, 2.2), kind='optical SNR candidate'),
    'G195.9+2.2': dict(ra=hms(6, 27, 0), dec=dms(16, 10, 0), size=(1.0, 1.0), kind='new optical SNR'),
    'G205.7-1.7': dict(ra=hms(6, 32, 0), dec=dms(5, 42, 0), size=(1.4, 1.4), kind='optical SNR candidate'),
}
# G203.1+6.6 and G206.9+2.3 only quoted by Galactic name in Fesen et al.; use (l,b) from the name.
NEIGHBOURS_LB = {
    'G203.1+6.6': dict(l=203.1, b=6.6, size=(1.6, 1.6), kind='known SNR (radio, Gao+2022)'),
    'G206.9+2.3': dict(l=206.9, b=2.3, size=(1.0, 1.0), kind='known SNR'),
}

# ---- Monogem circles: Knies et al. 2024 as digitised in Ziegenbalg's viewer links (diameters arcmin)
MONOGEM = {
    'single (Knies+2024)': dict(ra=107.08, dec=13.50, diam_arcmin=1436.0),
    'two-SNR: Monogem (Knies+2024)': dict(ra=104.15, dec=15.30, diam_arcmin=1320.0),
    'two-SNR: G205.6+12.4 (Knies+2024)': dict(ra=110.33, dec=12.94, diam_arcmin=1058.0),
    'Fesen+2026 Fig.31 (24 deg)': dict(ra=107.08, dec=13.50, diam_arcmin=24 * 60.0),
}

HIP37797 = dict(ra=116.25248502, dec=15.89503528,
                src='Hipparcos (I/311/hip2) via vizier.cfa.harvard.edu, accessed 2026-11-16')


def ellipse_ends(e):
    c = SkyCoord(e['ra'] * u.deg, e['dec'] * u.deg)
    half = e['maj'] / 120.0 * u.deg
    a = c.directional_offset_by(e['pa'] * u.deg, half)
    b = c.directional_offset_by((e['pa'] + 180) * u.deg, half)
    return c, a, b


def gc_segment_distance(p, a, b, n=2001):
    """Minimum angular distance (deg) from point p to great-circle segment a-b (sampled)."""
    sep = a.separation(b)
    pa = a.position_angle(b)
    t = np.linspace(0, 1, n)
    pts = a.directional_offset_by(pa, sep * t)
    d = p.separation(pts).deg
    i = int(np.argmin(d))
    return float(d[i]), pts[i], float(t[i])


def circle_crossings(cen, r_deg, a, b, n=4001):
    sep = a.separation(b)
    pa = a.position_angle(b)
    t = np.linspace(0, 1, n)
    pts = a.directional_offset_by(pa, sep * t)
    d = cen.separation(pts).deg - r_deg
    idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    return [pts[i] for i in idx], pts, d


def sb_mag_to_rayleigh(m, dlam_nm, lam_A=5007.0):
    """[O III] band surface brightness (AB mag/arcsec^2) -> Rayleigh, assuming all band flux is the
    5007 line and a top-hat effective bandwidth dlam_nm."""
    fnu = 3631e-23 * 10 ** (-0.4 * m)                      # erg/s/cm2/Hz/arcsec2
    dnu = 2.99792458e18 * (dlam_nm * 10) / lam_A ** 2       # Hz
    F = fnu * dnu                                           # erg/s/cm2/arcsec2
    E = 6.62607e-27 * 2.99792458e18 / lam_A                 # erg per photon
    R_per_arcsec2 = 1e6 / (4 * np.pi) / 4.25452e10 * E      # erg/s/cm2/arcsec2 per Rayleigh
    return F / R_per_arcsec2


def main():
    out = {'generated': '2026-11-16', 'script': 'scripts/discovery_geometry.py', 'ellipses': {},
           'quoted_extents': [], 'neighbours': {}, 'monogem_crossings': {}, 'surface_brightness': {},
           'astrobin_solutions': ASTROBIN_SOLUTIONS, 'hip37797': HIP37797}

    for k, e in NSNS.items():
        c, a, b = ellipse_ends(e)
        north, south = (a, b) if a.dec > b.dec else (b, a)
        rec = dict(e)
        rec.update(center_lb=[round(c.galactic.l.deg, 3), round(c.galactic.b.deg, 3)],
                   length_deg=round(e['maj'] / 60, 2), width_deg=round(e['mnr'] / 60, 2),
                   north_end_radec=[round(north.ra.deg, 3), round(north.dec.deg, 3)],
                   south_end_radec=[round(south.ra.deg, 3), round(south.dec.deg, 3)],
                   north_end_lb=[round(north.galactic.l.deg, 3), round(north.galactic.b.deg, 3)],
                   south_end_lb=[round(south.galactic.l.deg, 3), round(south.galactic.b.deg, 3)])
        out['ellipses'][k] = rec

    for q in QUOTED:
        L = q['length_deg']
        q = dict(q)
        q['phys_length_pc'] = {f'{d}pc': round(2 * d * np.tan(np.radians(L / 2)), 1) for d in (160, 260, 300)}
        out['quoted_extents'].append(q)

    # neighbours relative to the NSNS 2025 major axis (most complete published geometry)
    c, a, b = ellipse_ends(NSNS['O200.0_2025'])
    allnb = {}
    for k, v in NEIGHBOURS.items():
        p = SkyCoord(v['ra'] * u.deg, v['dec'] * u.deg)
        allnb[k] = (p, v)
    for k, v in NEIGHBOURS_LB.items():
        p = SkyCoord(l=v['l'] * u.deg, b=v['b'] * u.deg, frame='galactic').icrs
        allnb[k] = (p, v)
    for k, (p, v) in allnb.items():
        dmin, pt, t = gc_segment_distance(p, a, b)
        out['neighbours'][k] = dict(kind=v['kind'], ra=round(p.ra.deg, 3), dec=round(p.dec.deg, 3),
                                    l=round(p.galactic.l.deg, 2), b=round(p.galactic.b.deg, 2),
                                    size_deg=v['size'],
                                    sep_from_silk_axis_deg=round(dmin, 2),
                                    edge_gap_deg=round(dmin - max(v['size']) / 2, 2),
                                    sep_from_paper_center_deg=round(p.separation(SkyCoord(*CENTER, unit='deg')).deg, 2))
    for name, src in {'PSR B0656+14': MONOGEM_PULSAR, 'Geminga': GEMINGA}.items():
        p = SkyCoord(*src['radec'], unit='deg')
        dmin, pt, t = gc_segment_distance(p, a, b)
        out['neighbours'][name] = dict(kind='pulsar', ra=src['radec'][0], dec=src['radec'][1],
                                       sep_from_silk_axis_deg=round(dmin, 2),
                                       sep_from_paper_center_deg=round(p.separation(SkyCoord(*CENTER, unit='deg')).deg, 2))

    north, south = (a, b) if a.dec > b.dec else (b, a)
    for k, m in MONOGEM.items():
        cen = SkyCoord(m['ra'] * u.deg, m['dec'] * u.deg)
        r = m['diam_arcmin'] / 120.0
        xs, pts, d = circle_crossings(cen, r, south, north)
        inside_frac = float(np.mean(d < 0))
        out['monogem_crossings'][k] = dict(center=[m['ra'], m['dec']], radius_deg=round(r, 2),
                                           crossings_radec=[[round(x.ra.deg, 2), round(x.dec.deg, 2)] for x in xs],
                                           fraction_of_axis_inside=round(inside_frac, 3),
                                           sep_paper_center_deg=round(cen.separation(SkyCoord(*CENTER, unit='deg')).deg, 2))

    sb = {}
    for lab, m in [('Wu+2026 P1', 23.60), ('Wu+2026 P2', 23.87), ('Wu+2026 P3', 23.85), ('QHYCCD PR brightest', 24.3)]:
        sb[lab] = dict(mag_arcsec2=m, R_dlam6nm=round(sb_mag_to_rayleigh(m, 6.0), 1),
                       R_dlam5nm=round(sb_mag_to_rayleigh(m, 5.0), 1), R_dlam3nm=round(sb_mag_to_rayleigh(m, 3.0), 1))
    sb['NSNS K1a [OIII] (Ziegenbalg, estimated)'] = dict(R=2.04)
    sb["NSNS K1a' [OIII]"] = dict(R=1.53)
    sb['Fesen+2024 detection limit'] = dict(R='0.5-5')
    sb['assumptions'] = ('AB zero point 3631 Jy; entire in-band flux attributed to [OIII]5007; '
                         'top-hat bandwidth 6 nm (SWONS filter FWHM), 5 nm (QHY600 follow-up filter) or 3 nm')
    out['surface_brightness'] = sb
    tr = load_trace()
    out['trace_used'] = None if tr is None else tr.get('version')
    if tr and not tr.get('preliminary', True):
        sp = SkyCoord([p['ra'] for p in tr['spine']] * u.deg, [p['dec'] for p in tr['spine']] * u.deg)
        seg = sp[:-1].separation(sp[1:]).deg
        chord = sp[0].separation(sp[-1]).deg
        s_end, n_end = (sp[0], sp[-1]) if sp[0].dec < sp[-1].dec else (sp[-1], sp[0])
        tm = dict(version=tr['version'], n_points=len(sp), arc_length_deg=round(float(seg.sum()), 2),
                  chord_deg=round(float(chord), 2),
                  south_end_radec=[round(s_end.ra.deg, 3), round(s_end.dec.deg, 3)],
                  north_end_radec=[round(n_end.ra.deg, 3), round(n_end.dec.deg, 3)],
                  south_end_lb=[round(s_end.galactic.l.deg, 2), round(s_end.galactic.b.deg, 2)],
                  north_end_lb=[round(n_end.galactic.l.deg, 2), round(n_end.galactic.b.deg, 2)],
                  phys_arc_length_pc={f'{d}pc': round(float(np.radians(seg.sum()) * d), 1) for d in (160, 260, 300)},
                  main_strand_hwhm_deg=tr.get('half_width_deg'),
                  envelope_half_width_deg=tr.get('envelope_half_width_deg'),
                  curvature_full=dict((k, tr['curvature']['fits']['E_full'][k])
                                      for k in ('center_ra', 'center_dec', 'radius_deg', 'center_l', 'center_b')),
                  halpha_patch_moments=tr.get('halpha_patch_info', {}).get('moments'),
                  note='Spine of the main (E) [OIII] strand measured by the geometry teammate (data/derived/silk_trace.json); '
                       'its Dec range 12.0-28.25 is limited by detectability/coverage in the registered images, '
                       'NSNS O200.0b suggests faint emission continuing to Dec ~11.')
        # Monogem crossings along the measured spine
        cross = {}
        for k, m in MONOGEM.items():
            cen = SkyCoord(m['ra'] * u.deg, m['dec'] * u.deg)
            r = m['diam_arcmin'] / 120.0
            d = cen.separation(sp).deg - r
            inside = d < 0
            idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
            xs = []
            for i in idx:  # linear interpolation in Dec
                f = d[i] / (d[i] - d[i + 1])
                xs.append([round(float(sp[i].ra.deg + f * (sp[i + 1].ra.deg - sp[i].ra.deg)), 2),
                           round(float(sp[i].dec.deg + f * (sp[i + 1].dec.deg - sp[i].dec.deg)), 2)])
            L_in = float(np.sum(seg[inside[:-1] & inside[1:]]))
            cross[k] = dict(crossings_radec=xs, arc_inside_deg=round(L_in, 2),
                            arc_outside_deg=round(float(seg.sum()) - L_in, 2),
                            min_dist_to_rim_deg=round(float(np.min(np.abs(d))), 2))
        tm['monogem_crossings'] = cross
        # neighbours: min separation from the spine
        nb = {}
        for k, v in out['neighbours'].items():
            p = SkyCoord(v['ra'] * u.deg, v['dec'] * u.deg)
            dmin = float(np.min(p.separation(sp).deg))
            r = max(v.get('size_deg', (0, 0))) / 2 if 'size_deg' in v else 0.0
            nb[k] = dict(sep_from_spine_deg=round(dmin, 2), edge_gap_deg=round(dmin - r, 2))
        tm['neighbours'] = nb
        hinfo = tr.get('halpha_patch_info', {})
        hm = hinfo.get('moments') or hinfo.get('moments_halfmax') or {}
        tm['halpha_patch_moments'] = hm
        tm['halpha_patch_moments_30pct'] = hinfo.get('moments_30pct')
        hp = SkyCoord(hm['centroid_ra'] * u.deg, hm['centroid_dec'] * u.deg)
        o2a = SkyCoord(NSNS['O200.0a_Ha_SII']['ra'] * u.deg, NSNS['O200.0a_Ha_SII']['dec'] * u.deg)
        tm['halpha_patch_vs_NSNS_O200.0a_sep_deg'] = round(float(hp.separation(o2a).deg), 2)
        tm['paper_center_to_spine_mid_sep_deg'] = round(float(SkyCoord(*CENTER, unit='deg').separation(
            sp[len(sp) // 2]).deg), 2)
        cc = SkyCoord(tm['curvature_full']['center_ra'] * u.deg, tm['curvature_full']['center_dec'] * u.deg)
        refs = {'Monogem centre (Knies+2024 single)': (107.08, 13.50), 'G205.6+12.4 centre': (110.33, 12.94),
                'PSR B0656+14': MONOGEM_PULSAR['radec'], 'Geminga (now)': GEMINGA['radec'],
                'Ziegenbalg Geminga-SNR centre': (90.1412, 11.8861),
                'Ziegenbalg Geminga birthplace': (84.7088, 7.6303)}
        tm['curvature_centre_separations_deg'] = {k: round(float(cc.separation(SkyCoord(*v, unit='deg')).deg), 2)
                                                  for k, v in refs.items()}
        rim = {}
        for k, m in MONOGEM.items():
            cen = SkyCoord(m['ra'] * u.deg, m['dec'] * u.deg)
            d = cen.separation(sp).deg - m['diam_arcmin'] / 120.0
            rim[k] = {f'dec{x}': round(float(np.interp(x, sp.dec.deg, d)), 2) for x in (12, 16, 20, 24, 28)}
        tm['spine_distance_to_rim_deg_(neg=inside)'] = rim
        # arc length per source coverage: how much of the strand each image actually shows
        dec = sp.dec.deg

        def arc(lo, hi):
            m = (dec[:-1] >= lo - 1e-6) & (dec[1:] <= hi + 1e-6)
            return round(float(seg[m].sum()), 2)
        cov = {}
        for tag in ('swons', 'nsnsN', 'nsnsW', 'fig4'):
            ds = [p['dec'] for p in tr['spine'] if tag in p.get('source', '')]
            if ds:
                cov[tag] = dict(dec_range=[min(ds), max(ds)], arc_length_deg=arc(min(ds), max(ds)))
        tm['arc_length_by_detecting_image'] = cov
        tm['arc_length_split_deg'] = {'dec12.0-21.5 (bright S part, up to Halpha patch)': arc(12.0, 21.5),
                                      'dec21.5-28.25 (faint N arc)': arc(21.5, 28.25)}
        f1 = ROOT / 'data/derived/geometry_fig1_wcs.json'
        if f1.exists():
            w = json.loads(f1.read_text())
            decs = [v[1] for v in w['corners_radec'].values()]
            tm['wu_fig1a_footprint_dec'] = [round(min(decs), 2), round(max(decs), 2)]
        out['measured_trace'] = tm
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
