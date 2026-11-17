"""Assemble data/derived/context_objects.csv: one row per object in the Silk Nebula's astrophysical neighbourhood.

Combines the derived tables of the context topic (pulsars, SNRs/outlines, hot stars, PNe, catalogued nebulae)
with a few named regions (SIMBAD CfA-mirror positions, 2026-11-16). Columns:
name, class, ra, dec, l, b, dist_pc, size_deg, age, key_params, sep_spine_deg, relevance, source
"""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *

D = ROOT / 'data/derived'
tr = load_trace(); sp = tr['spine']
SRA = np.array([q['ra'] for q in sp]); SDE = np.array([q['dec'] for q in sp])


def sep_spine(ra, dec):
    return float(angsep(np.array([ra])[:, None], np.array([dec])[:, None], SRA[None], SDE[None]).min())


rows = []


def add(**k):
    l, b = gal(k['ra'], k['dec'])
    k.setdefault('l', float(l[0])); k.setdefault('b', float(b[0]))
    k['sep_spine_deg'] = round(sep_spine(k['ra'], k['dec']), 3)
    rows.append(k)


# Silk reference geometry
for nm, (ra, dec) in REGIONS.items():
    add(name=f'Silk {nm}', **{'class': 'silk_reference'}, ra=ra, dec=dec, dist_pc='', size_deg='', age='',
        key_params='Wu+2026 Table 1 region', relevance='LAMOST MRS-N region', source='Wu+2026 arXiv:2610.13329')
cv = tr['curvature']
add(name='Silk centre of curvature (E strand, single circle)', **{'class': 'silk_reference'}, ra=cv['center_ra'],
    dec=cv['center_dec'], dist_pc='', size_deg=cv['radius_deg'], age='',
    key_params=f"radius {cv['radius_deg']} deg (boot 16-84 {cv['radius_16_84']}); reduced chi2>>1, S-shaped arc",
    relevance='test locus for explosion centres', source='data/derived/silk_trace.json v1.0 (geometry)')
hp_ = tr['halpha_patch_info']['moments_halfmax']
add(name='Silk northern Halpha/[SII] patch', **{'class': 'silk_reference'}, ra=hp_['centroid_ra'], dec=hp_['centroid_dec'],
    dist_pc='', size_deg=f"{hp_['major_deg']}x{hp_['minor_deg']}", age='', key_params=f"PA {hp_['pa_deg']}",
    relevance='where [OIII] arc meets Halpha+[SII] (NSNS O200.0a)', source='silk_trace.json v1.0')

# SNRs, outlines
snr = pd.read_csv(D / 'context_snrs.csv')
for _, s in snr.iterrows():
    if not np.isfinite(s.ra):
        continue
    rel = ''
    if 'Monogem Ring (G203+12)' in s['name']:
        rel = f"{100*s.silk_inside_fraction:.0f}% of the Silk spine inside the 11.97-deg outline; spine-centre distance {s.spine_dist_from_centre_mean_deg:.1f}+-{s.spine_dist_from_centre_rms_deg:.1f} deg (not concentric)"
    elif 'G205.6' in s['name']:
        rel = f"{100*s.silk_inside_fraction:.0f}% of spine inside; Silk near its E rim"
    add(name=s['name'], **{'class': 'SNR/outline: ' + s.kind}, ra=s.ra, dec=s.dec, l=s.l, b=s.b,
        dist_pc=s.dist_pc if isinstance(s.dist_pc, str) else '', size_deg=s.diam,
        age=s.age if isinstance(s.age, str) else '', key_params='', relevance=rel, source=s.ref)

# pulsars
ps = pd.read_csv(D / 'context_pulsars_30deg.csv')
pc = json.loads((D / 'context_pulsar_coincidence.json').read_text())['pulsars']
for _, p in ps[ps.dist_best_kpc < 1.0].iterrows():
    nm = p.psrb if isinstance(p.psrb, str) else p.psrj
    key = 'Geminga' if p.psrj == 'J0633+1746' else nm
    rel = ''
    if key in pc:
        c = pc[key]
        rel = (f"min sep of 0-1 Myr LSR track from curvature centre {c['min_sep_track_from_nominal_centre_deg_p16_50_84'][1]} deg "
               f"at t={c['t_closest_myr_p16_50_84'][1]} Myr")
        if 'birth_at_tau' in c:
            rel += f"; birth locus at tau_c is {c['birth_at_tau']['sep_from_nominal_centre_deg_p16_50_84'][1]} deg from it"
    add(name=f"PSR {p.psrj}" + (f" ({p.psrb})" if isinstance(p.psrb, str) else ''), **{'class': 'pulsar'}, ra=p.ra, dec=p.dec,
        dist_pc=round(1000 * p.dist_best_kpc), size_deg='', age=f"tau_c={p.age_yr:.3g} yr" if np.isfinite(p.age_yr) else '',
        key_params=f"DM={p.dm} pc/cm3; PM=({p.pmra},{p.pmdec}) mas/yr; Edot={p.edot:.2g} erg/s; dist {p.dist_src}",
        relevance=rel or 'no PM: not back-traced', source='ATNF psrcat v2.9.0')

# hot stars
hs = pd.read_csv(D / 'context_hotstars.csv')
keep = hs[((hs.dist_pc < 1000) & (hs.teff >= 40000)) |
          hs.name.isin(['HD  60848', 'GALEX J074320.3+160752', 'NGC  2392'])]
for _, h in keep.iterrows():
    add(name=h['name'], **{'class': 'hot star (' + (str(h.sptype) if isinstance(h.sptype, str) else 'WD GF21') + ')'},
        ra=h.ra, dec=h.dec, dist_pc=round(h.dist_pc), size_deg='', age='',
        key_params=(f"Teff={h.teff:.0f} K ({h.teff_ref}); Q(H)={h.Q_H:.2g}/s; Q(>35eV)={h.Q_35eV:.2g}/s; "
                    f"R_S(n=0.1/1)={h['Rs_pc_n0.1']:.1f}/{h['Rs_pc_n1.0']:.1f} pc; L={h.L_bol_lsun:.3g} Lsun"),
        relevance=(f"d_proj to spine {h.d_proj_pc:.1f} pc; [OIII] photon budget max/required = {h.budget_ratio_max_over_req:.2g}; "
                   f"wake frac {h.wake_frac_spine_within_0p5deg_of_track:.2f}"),
        source='SIMBAD/GF21/literature (see context_hotstars.py)')

# PNe within 15 deg (SIMBAD)
s = pd.read_csv(ROOT / 'data/raw/context_catalogs/simbad_hotstars_cone25.csv')
pn = s[(s.otype == 'PN')]
for _, r in pn.iterrows():
    if angsep(r.ra, r.dec, *CENTER) > 15:
        continue
    if r.main_id in keep.name.values:
        continue
    d = round(1000 / r.plx_value) if np.isfinite(r.plx_value) and r.plx_value > 0 else ''
    add(name=r.main_id, **{'class': 'planetary nebula'}, ra=r.ra, dec=r.dec, dist_pc=d, size_deg='', age='',
        key_params=f"CSPN sp={r.sp_type if isinstance(r.sp_type, str) else ''}; Gaia plx={r.plx_value}",
        relevance='', source='SIMBAD TAP (CfA mirror)')

# catalogued Galactic nebulae near the spine / Halpha patch
nb = pd.read_csv(D / 'context_nebulae_near_trace.csv')
nb = nb[nb.otype.isin(['sh', 'Cld', 'MoC', 'HVC', 'HII', 'RNe', 'DNe', 'ISM', 'bub', 'SNR'])]
for _, r in nb.iterrows():
    add(name=r.main_id, **{'class': 'catalogued nebula/cloud (' + r.otype + ')'}, ra=r.ra, dec=r.dec, dist_pc='',
        size_deg=(r.galdim_majaxis / 60 if np.isfinite(r.galdim_majaxis) else ''), age='',
        key_params=(f"v={r.rvz_radvel:.0f} km/s" if np.isfinite(r.rvz_radvel) else ''),
        relevance=f"{r.sep_halpha_patch_deg:.2f} deg from Halpha patch", source='SIMBAD TAP (CfA mirror)')

# named regions
for nm, ra, dec, size, dist, src in [
        ('Origem Loop (G194.7+0.4; non-thermal Origem Arc)', 94.1342, 16.0724, 6.67, '~1700 (Origem Arc, Gao & Han 2013)', 'SIMBAD; Gao & Han 2013 A&A 551 A16 (arXiv:1212.1523)'),
        ('lambda Ori complex', 83.775, 9.933, '', '~400', 'SIMBAD'),
        ("Barnard's Loop", 81.875, -3.967, '', '', 'SIMBAD'),
        ('Gem OB1 molecular cloud', 94.75, 23.0, 5.0, '', 'SIMBAD (galdim 300 arcmin)'),
        ('Geminga birthplace (Pellizza+ via NSNS)', 84.7088, 7.6303, '', '90-240', 'simg.de/nebulae3/gem-01.html annotation'),
        ('NSNS sphere centre if arc on limb (C0)', 73.1139, 6.5749, '', '', 'simg.de/nebulae3/gem-01.html annotation')]:
    add(name=nm, **{'class': 'named region'}, ra=ra, dec=dec, dist_pc=dist, size_deg=size, age='', key_params='',
        relevance=f"{angsep(ra, dec, cv['center_ra'], cv['center_dec']):.2f} deg from Silk curvature centre", source=src)

df = pd.DataFrame(rows)[['name', 'class', 'ra', 'dec', 'l', 'b', 'dist_pc', 'size_deg', 'age', 'key_params',
                         'sep_spine_deg', 'relevance', 'source']]
df.to_csv(D / 'context_objects.csv', index=False, float_format='%.5g')
print(len(df), 'rows'); print(df['class'].value_counts().to_string())
