"""Known SNRs, SNR candidates, loops/superbubbles and catalogued nebulae around the Silk Nebula.

Inputs : data/raw/context_catalogs/green2025_VII297_all.tsv (Green 2025, VizieR VII/297 via vizier.cfa.harvard.edu)
         data/raw/context_catalogs/simbad_nebulae_box.csv    (SIMBAD TAP, CfA mirror; RA 100-135, Dec -2..36)
         literature values typed in below (with references)
Outputs: data/derived/context_snrs.csv          SNRs/candidates/loops with l 175-235, b -15..+40 + literature additions
         data/derived/context_nebulae_near_trace.csv  Galactic SIMBAD nebulae within 3 deg of the spine / Halpha patch
"""
import sys, re
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *

RAW = ROOT / 'data/raw/context_catalogs'

# literature objects (centre RA, Dec deg; diameter deg; distance; reference)
LIT = [
    dict(name='Monogem Ring (G203+12)', ra=107.08, dec=13.50, diam=23.94, kind='SNR (X-ray)', dist_pc='~300 (PSR B0656+14 parallax 288+33-27)',
         age='86 kyr (Plucinsky+96 via Thorsett+03); ~1.1-1.2e5 yr (Knies+24)', ref='Knies+2024 A&A (arXiv:2401.17289) Table 4 single-SNR; Fesen+2026 Table 1'),
    dict(name='Monogem Ring (two-SNR scenario, W)', ra=104.15, dec=15.30, diam=22.0, kind='SNR (X-ray)', dist_pc='300',
         age='~1.1e5 yr', ref='Knies+2024 Table 4'),
    dict(name='G205.6+12.4 (eROSITA SE hot arc)', ra=110.33, dec=12.94, diam=17.64, kind='SNR candidate (X-ray)', dist_pc='300-500',
         age='40-55 kyr', ref='Knies+2024 Table 4 / Sect. 7.5'),
    dict(name='G190.4+12.5', ra=104.24, dec=25.56, diam=6.6, kind='SNR candidate (X-ray, [OIII] filaments)', dist_pc='>1500 (Knies+24)',
         age='40-60 kyr', ref='Knies+2024; Fesen+2026 Table 1 (06:57 +25:34, 6.1 deg)'),
    dict(name='G191.4+11.1', ra=103.25, dec=24.20, diam=3.2, kind='SNR candidate (optical)', dist_pc='', age='',
         ref='Fesen+2026 Table 1 (06:53 +24:12, 3.2x2.2 deg)'),
    dict(name='G195.9+2.2', ra=96.75, dec=16.17, diam=1.0, kind='SNR (optical, new)', dist_pc='', age='',
         ref='Fesen+2026 Table 1 (06:27 +16:10)'),
    dict(name='G205.7-1.7', ra=98.0, dec=5.70, diam=1.4, kind='SNR candidate (optical)', dist_pc='', age='',
         ref='Fesen+2026 Table 1 (06:32 +05:42)'),
    dict(name='G206.6+6.1 (=G206.7+5.9)', ra=105.5, dec=8.38, diam=3.6, kind='SNR (radio+optical)', dist_pc='~440 (Gao+22 kinematic)',
         age='', ref='Fesen+2026 Table 1 (07:02 +08:23); Gao+2022 FAST'),
    dict(name='G203.1+6.6 (=G203.1+6.7)', ra=np.nan, dec=np.nan, diam=1.6, kind='SNR (radio)', dist_pc='', age='',
         ref='Gao+2022 FAST (via Fesen+2026); position via Galactic name'),
    dict(name='Gemini Halpha ring', ra=np.nan, dec=np.nan, diam=5.7, kind='wind-blown bubble (Halpha, CIV)', dist_pc='200-350',
         age='', ref='Kim+2007; Knies+2018; Fesen+2026 Sect. 6.2.1 ((l,b)~(191.5,+5.0))'),
    dict(name='Ziegenbalg possible Geminga SNR', ra=90.1412, dec=11.8861, diam=55.3, kind='hypothesis (NSNS)', dist_pc='250',
         age='342 kyr', ref='simg.de/nebulae3/gem-01.html annotation (diameter 3317 arcmin)'),
    # optical filaments previously linked to Monogem (Fesen+2026 Fig. 31 labels R20, W06)
    dict(name='R20 Halpha filament G203.7+11.5', ra=np.nan, dec=np.nan, diam=4.5, kind='optical filament (Halpha, Faraday screen)',
         dist_pc='300 (assumed)', age='', ref='Reich, Reich & Sun 2020 A&A (arXiv:2007.08221): 4.5 deg long; B_par 26+-5 uG, ne 1.6 cm^-3, RM -86 rad/m2'),
    dict(name='W06 Halpha filament', ra=112.158, dec=5.165, diam=np.nan, kind='optical filament (Halpha, no [OIII])',
         dist_pc='', age='', ref='Weinberger, Temporin & Stecklum 2006 A&A (astro-ph/0511546): spectrum at 07:28:38 +05:09:55, ne<100'),
    # HI supershell candidates whose catalogue positions fall on the Silk (VizieR, CfA mirror, 2026-11-16)
    dict(name='GS202+17+019 HI supershell (Suad+14)', ra=np.nan, dec=np.nan, diam=2 * np.sqrt(5.7 * 4.5), kind='HI supershell candidate',
         dist_pc='2500 (rotation-curve kinematic; unreliable toward the anticentre)', age='',
         ref='Suad+2014 A&A 564 A116 (VizieR J/A+A/564/A116): ellipse centre (202.4,+16.9), a,b=5.7,4.5 deg, v0=+18.6, dv=15.5 km/s, HI minimum at (203.3,+19.1)'),
    dict(name='GS202+19-019 HI shell (Ehlerova+05)', ra=np.nan, dec=np.nan, diam=6.0, kind='HI shell candidate',
         dist_pc='', age='', ref='Ehlerova & Palous 2005 A&A 437 101 (VizieR J/A+A/437/101): (202.3,+18.8), 6.5x5.5 deg, v=-18.8, dv=21.6 km/s, quality 3'),
    dict(name='GS208+15+003 HI shell (Ehlerova+05)', ra=np.nan, dec=np.nan, diam=3.75, kind='HI shell candidate',
         dist_pc='', age='', ref='Ehlerova & Palous 2005 (VizieR J/A+A/437/101): (207.6,+14.8), 3.5x4.0 deg, v=+2.5 km/s, quality 3'),
]
LIT_GAL = {'G203.1+6.6 (=G203.1+6.7)': (203.1, 6.6), 'Gemini Halpha ring': (191.5, 5.0),
           'R20 Halpha filament G203.7+11.5': (203.7, 11.5), 'GS202+17+019 HI supershell (Suad+14)': (202.4, 16.9),
           'GS202+19-019 HI shell (Ehlerova+05)': (202.3, 18.8), 'GS208+15+003 HI shell (Ehlerova+05)': (207.6, 14.8)}


def main():
    g = pd.read_csv(RAW / 'green2025_VII297_all.tsv', sep='\t', comment='#').iloc[2:].copy()
    g.columns = [c.strip() for c in g.columns]
    g['l'] = g.SNR.str.slice(1, 6).astype(float)
    g['b'] = g.SNR.str.slice(6).astype(float)
    c = SkyCoord(g.RAJ2000.str.strip(), g.DEJ2000.str.strip(), unit=(u.hourangle, u.deg))
    g['ra'], g['dec'] = c.ra.deg, c.dec.deg
    g['diam'] = pd.to_numeric(g.MajDiam, errors='coerce') / 60
    sel = g[(g.l >= 175) & (g.l <= 235) & (g.b >= -15) & (g.b <= 40)].copy()
    rows = [dict(name=r.SNR + (' ' + r.Names.strip() if isinstance(r.Names, str) and r.Names.strip() else ''),
                 ra=r.ra, dec=r.dec, l=r.l, b=r.b, diam=r.diam, kind='SNR (Green 2025, type ' + str(r.type).strip() + ')',
                 dist_pc='', age='', ref='Green 2025 (VizieR VII/297)') for _, r in sel.iterrows()]
    for d in LIT:
        d = dict(d)
        if d['name'] in LIT_GAL:
            l0, b0 = LIT_GAL[d['name']]
            cc = SkyCoord(l0 * u.deg, b0 * u.deg, frame='galactic').icrs
            d['ra'], d['dec'] = cc.ra.deg, cc.dec.deg
        ll, bb = gal(d['ra'], d['dec'])
        d['l'], d['b'] = ll[0], bb[0]
        rows.append(d)
    df = pd.DataFrame(rows)
    tr = load_trace(); sp = tr['spine']
    tra = np.array([q['ra'] for q in sp]); tde = np.array([q['dec'] for q in sp])
    df['sep_centre_to_silk_spine_min_deg'] = angsep(df.ra.values[:, None], df.dec.values[:, None], tra[None], tde[None]).min(1)
    df['spine_dist_from_centre_mean_deg'] = angsep(df.ra.values[:, None], df.dec.values[:, None], tra[None], tde[None]).mean(1)
    df['spine_dist_from_centre_rms_deg'] = angsep(df.ra.values[:, None], df.dec.values[:, None], tra[None], tde[None]).std(1)
    df['radius_deg'] = df.diam / 2
    df['silk_inside_fraction'] = (angsep(df.ra.values[:, None], df.dec.values[:, None], tra[None], tde[None]) <
                                  df.radius_deg.values[:, None]).mean(1)
    df.to_csv(ROOT / 'data/derived/context_snrs.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print(df[['name', 'l', 'b', 'radius_deg', 'sep_centre_to_silk_spine_min_deg', 'spine_dist_from_centre_mean_deg',
              'spine_dist_from_centre_rms_deg', 'silk_inside_fraction']].to_string())

    # Galactic catalogued nebulae near the spine and the Halpha patch
    s = pd.read_csv(RAW / 'simbad_nebulae_box.csv')
    extragal = s.main_id.str.contains(r'DDO|NGC 2403|\[SHK91\]|\[PSB2020\]|\[LBZ2013\]|\[ZBF2015\]|\[SRM2012\]|HIPASS|AGC|'
                                      r'\[KM98\]|\[ABB2014\]|\[FB2011\]|\[KC97c\]|\[VH2006\]|\[C51\]|\[WB2001\]|\[BRW2009\]|\[KBE97\]', regex=True)
    s = s[~extragal & (s.otype != 'LSB')].copy()
    hp = tr['halpha_patch_info']['moments_halfmax']
    s['sep_spine_deg'] = angsep(s.ra.values[:, None], s.dec.values[:, None], tra[None], tde[None]).min(1)
    s['sep_halpha_patch_deg'] = angsep(s.ra.values, s.dec.values, hp['centroid_ra'], hp['centroid_dec'])
    near = s[(s.sep_spine_deg < 3) | (s.sep_halpha_patch_deg < 3)].sort_values('sep_spine_deg')
    near.to_csv(ROOT / 'data/derived/context_nebulae_near_trace.csv', index=False, float_format='%.4f')
    print(near[['main_id', 'otype', 'ra', 'dec', 'galdim_majaxis', 'rvz_radvel', 'sep_spine_deg', 'sep_halpha_patch_deg']].to_string())


if __name__ == '__main__':
    main()
