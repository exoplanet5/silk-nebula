"""w4_distance: target list for a high-ion (C IV / Si IV / N V) absorption distance bracket of the Silk.

Steps
 1. Gaia DR3 blue stars (bp_rp < 0.10, parallax > 1.1 mas; data/raw/w4_distance_cat/gaia_hot_box.csv, Data Lab TAP)
    -> keep finite parallax, parallax_over_error > 5, 1/plx in 120-800 pc (plus 'near-miss' 100-120 / 800-1000 pc).
    NOTE: PostgreSQL treats NaN > x as TRUE, so NaN parallaxes leaked into the box query; they are removed here.
 2. Angular distance to the E strand (spine), W and W2 strands of data/derived/silk_trace.json (polylines,
    tangent-plane offsets); signed offset from the E spine (+ = east, the convex side away from Monogem).
    Keep <= 0.5 deg; flag <= 0.3 deg (primary) and <= 0.14 deg (inside the [O III] HWHM).
 3. Class from M_G = G + 5 log10(plx/100) and BP-RP: WD (M_G > 8), hot subdwarf (2.5 < M_G <= 8, BP-RP < -0.1),
    early-type MS/giant (M_G <= 2.5), other.
 4. GALEX GR6/7 (MAST Catalogs 'Galex'; the VizieR mirror gave no rows) cone 6\" per target -> FUV, NUV (AB).
 5. Known-object cross-ids from catalogues already on disk (Gentile Fusillo+21 WDs, Culpan+22 / Geier+20 hot
    subdwarfs, SIMBAD hot-star cone) within 3".
 6. MAST: the CAOM inventory already on disk (data/raw/spectra_ism_mast/mast_uv_{hst,fuse_iue_hut}.csv): every
    HST/IUE/FUSE/HUT spectrum in the field; matched to targets (<= 5") and listed for the corridor (<= 0.5 deg of any strand).
 7. HST/COS G160M time for S/N = 15 per resolution element at C IV 1550 from GALEX FUV:
       f_lambda(1550) from FUV AB (lambda_eff 1539 A, flat f_nu assumed across 1539-1550);
       A_eff(G160M, 1550 A) = 1000 cm^2 (+-30 %; approximate current-epoch value, CHECK WITH THE COS ETC);
       resel = 0.0735 A (6 px x 12.25 mA); N = 225 counts; overhead: 1st orbit 2000 s usable, then 2700 s/orbit.
    Bright-object flag: global G160M count rate > 15000 c/s (approx; FUV <~ 12) -> 'BOP: use STIS E140H/M'.
Outputs: data/derived/w4_distance_targets.csv, data/derived/w4_distance_targets_mast.csv,
         data/derived/w4_distance_targets_summary.json; raw: data/raw/w4_distance_cat/galex_*.tsv, mast_*.ecsv
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_targets.py [--nofetch]
"""
import sys, json, time, datetime, io
from pathlib import Path
import numpy as np
import pandas as pd
import requests

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/w4_distance_cat'
LOG = RAW / 'fetch_log.tsv'
VIZ = 'https://vizier.cfa.harvard.edu/viz-bin/asu-tsv'
NOFETCH = '--nofetch' in sys.argv


def log(url, fname, status, nbytes, note=''):
    with open(LOG, 'a') as f:
        f.write(f'{datetime.datetime.now(datetime.timezone.utc).isoformat()}\t{url}\t{fname}\t{status}\t{nbytes}\t{note}\n')


def read_viz_tsv(text):
    L = [l for l in text.splitlines() if l and not l.startswith('#')]
    if len(L) < 3:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO('\n'.join([L[0]] + L[3:])), sep='\t')


# ---------------- strands ----------------
tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
STR = {k: np.array([[p['ra'], p['dec']] for p in tr['strands'][k]['points']]) for k in ('E', 'W', 'W2')}
STR['E'] = np.array([[p['ra'], p['dec']] for p in tr['spine']])


def poly_dist(ra, dec, pts, dense=0.01):
    """min angular distance (deg, small-angle) from (ra,dec) to polyline pts (sorted by dec); signed (+ east)."""
    dd = np.arange(pts[:, 1].min(), pts[:, 1].max() + 1e-9, dense)
    rr = np.interp(dd, pts[:, 1], pts[:, 0])
    dx = (ra[:, None] - rr[None, :]) * np.cos(np.radians(dec[:, None]))
    dy = dec[:, None] - dd[None, :]
    d = np.hypot(dx, dy)
    j = d.argmin(1)
    sgn = np.sign(dx[np.arange(len(ra)), j])
    dmin = d[np.arange(len(ra)), j]
    # beyond the ends -> NaN-like large
    return dmin * sgn, dd[j]


g = pd.read_csv(RAW / 'gaia_hot_box.csv')
g = g[np.isfinite(g.parallax) & np.isfinite(g.parallax_over_error)].copy()
g['dist_pc'] = 1000 / g.parallax
g['dist_lo'] = 1000 / (g.parallax + g.parallax_error)
g['dist_hi'] = 1000 / np.clip(g.parallax - g.parallax_error, 1e-3, None)
g['M_G'] = g.phot_g_mean_mag + 5 * np.log10(g.parallax / 100)
for k, pts in STR.items():
    s, dd = poly_dist(g.ra.values, g.dec.values, pts)
    inside = (g.dec.values >= pts[:, 1].min() - 0.3) & (g.dec.values <= pts[:, 1].max() + 0.3)
    g[f'off_{k}'] = np.where(inside, s, np.nan)
g['sep_min'] = g[['off_E', 'off_W', 'off_W2']].abs().min(1)
g = g[np.isfinite(g.sep_min)].copy()
g['nearest'] = g[['off_E', 'off_W', 'off_W2']].abs().fillna(99).idxmin(1).str.replace('off_', '')
sel = (g.parallax_over_error > 5) & (g.dist_pc >= 100) & (g.dist_pc <= 1000) & (g.sep_min <= 0.5)
t = g[sel].copy()


def cls(r):
    if r.M_G > 8:
        return 'WD'
    if r.M_G > 2.5 and r.bp_rp < -0.1:
        return 'hot subdwarf'
    if r.M_G <= 2.5:
        return 'early-type (B/A)'
    return 'other blue'


t['class'] = t.apply(cls, axis=1)
t = t.sort_values('dist_pc').reset_index(drop=True)
print(len(t), 'candidates within 0.5 deg, plx/err>5, 100-1000 pc')

# ---------------- GALEX per target ----------------
# The CfA VizieR mirror returns headers but no rows for II/312 and II/335 (tested 2026-11-17 on a known GALEX
# source), so GALEX GR6/7 photometry comes from MAST (astroquery.mast Catalogs, catalog='Galex'), 6" cones.
from astroquery.mast import Catalogs
fu, nu, efu, enu, gsep = [], [], [], [], []
for i, r in t.iterrows():
    fn = RAW / f'galex_mast_{r.source_id}.csv'
    if not fn.exists() and not NOFETCH:
        for k in range(3):
            try:
                q = Catalogs.query_region(f'{r.ra} {r.dec}', radius=6 / 3600, catalog='Galex')
                q.to_pandas().to_csv(fn, index=False)
                log(f"astroquery.mast Catalogs.query_region('{r.ra:.6f} {r.dec:.6f}', radius=6 arcsec, catalog='Galex')", fn.name, 200, fn.stat().st_size, f'{len(q)} rows')
                break
            except Exception as e:
                log('MAST Catalogs Galex', fn.name, 'ERR', 0, repr(e)[:150]); time.sleep(5)
    df = pd.read_csv(fn) if fn.exists() and fn.stat().st_size > 5 else pd.DataFrame()
    if len(df):
        df['sep'] = np.hypot((df.ra - r.ra) * np.cos(np.radians(r.dec)), df.dec - r.dec) * 3600
        df = df.sort_values('sep')
        fu.append(df.fuv_mag.iloc[0]); efu.append(df.fuv_magerr.iloc[0] if 'fuv_magerr' in df else np.nan)
        nu.append(df.nuv_mag.iloc[0]); enu.append(df.nuv_magerr.iloc[0] if 'nuv_magerr' in df else np.nan)
        gsep.append(df.sep.iloc[0])
    else:
        fu.append(np.nan); efu.append(np.nan); nu.append(np.nan); enu.append(np.nan); gsep.append(np.nan)
t['FUV'] = pd.to_numeric(pd.Series(fu), errors='coerce').where(lambda x: x > 0).values
t['e_FUV'] = efu; t['NUV'] = pd.to_numeric(pd.Series(nu), errors='coerce').where(lambda x: x > 0).values
t['e_NUV'] = enu; t['galex_sep_arcsec'] = gsep

# ---------------- known-object cross-ids from disk catalogues ----------------
def disk_cat(fname, ra_c, de_c, name_c, extra):
    p = ROOT / 'data/raw/context_catalogs' / fname
    if not p.exists():
        return None
    df = read_viz_tsv(p.read_text())
    if not len(df):
        return None
    df = df.rename(columns={ra_c: 'ra', de_c: 'dec', name_c: 'name'})
    return df[['ra', 'dec', 'name'] + [c for c in extra if c in df.columns]]


cats = {'GF21_WD': disk_cat('gf21_wd_cone20_hot.tsv', 'RA_ICRS', 'DE_ICRS', 'WDJname', ['TeffH', 'loggH']),
        'Culpan22cand': disk_cat('culpan22_hotsd_cone20.tsv', 'RA_ICRS', 'DE_ICRS', 'GaiaEDR3', ['GMAG']),
        'Culpan22known': disk_cat('culpan22_knownhsd_cone20.tsv', 'RA_ICRS', 'DE_ICRS', 'Name', ['SpClass', 'Teff']),
        'Geier20': disk_cat('geier20_sdcat_cone20.tsv', 'RAJ2000', 'DEJ2000', 'Name', ['SpClass', 'Teff'])}
xid = [''] * len(t)
for cname, df in cats.items():
    if df is None:
        continue
    for i, r in t.iterrows():
        d = np.hypot((df.ra.astype(float) - r.ra) * np.cos(np.radians(r.dec)), df.dec.astype(float) - r.dec) * 3600
        j = np.argmin(d.values)
        if d.values[j] < 3:
            extra = ' '.join(f'{c}={df.iloc[j][c]}' for c in df.columns[3:])
            xid[i] += f'{cname}:{str(df.iloc[j]["name"]).strip()} {extra}; '
sim = ROOT / 'data/raw/context_catalogs/simbad_hotstars_cone25.csv'
if sim.exists():
    s = pd.read_csv(sim)
    rc = [c for c in s.columns if c.lower() in ('ra', 'ra_d', 'ra_icrs')][0]
    dc = [c for c in s.columns if c.lower() in ('dec', 'dec_d', 'de_icrs')][0]
    nc = [c for c in s.columns if 'main_id' in c.lower() or c.lower() == 'name'][0]
    sc = [c for c in s.columns if 'sp' in c.lower()][:1]
    for i, r in t.iterrows():
        d = np.hypot((s[rc] - r.ra) * np.cos(np.radians(r.dec)), s[dc] - r.dec) * 3600
        j = int(np.argmin(d.values))
        if d.values[j] < 3:
            xid[i] += f'SIMBAD:{s.iloc[j][nc]} {s.iloc[j][sc[0]] if sc else ""}; '
t['cross_ids'] = xid

# ---------------- COS exposure estimate ----------------
H, C = 6.626e-27, 2.998e18
LAM, AEFF, RESEL = 1550.0, 1000.0, 0.0735


def cos_time(fuv, snr=15.0, aeff=AEFF):
    if not np.isfinite(fuv):
        return np.nan, np.nan
    fnu = 10 ** (-0.4 * (fuv + 48.6))
    flam = fnu * C / LAM ** 2
    nph = flam / (H * C / LAM)                       # ph/s/cm2/A
    rate_resel = nph * aeff * RESEL
    glob = nph * aeff * 290.0                        # ~ one G160M segment
    return snr ** 2 / rate_resel, glob


tt, gl = zip(*[cos_time(f) for f in t.FUV])
t['cos_t_s_SN15'] = np.round(tt, 0)
t['cos_global_cps'] = np.round(gl, 0)


def orbits(ts):
    if not np.isfinite(ts):
        return np.nan
    return 1 if ts <= 2000 else 1 + int(np.ceil((ts - 2000) / 2700))


t['cos_orbits_SN15'] = [orbits(x) for x in t.cos_t_s_SN15]
t['cos_flag'] = np.where(t.cos_global_cps > 15000, 'BOP: too bright for COS; STIS E140H/M',
                 np.where(t.cos_orbits_SN15 <= 5, 'feasible <=5 orbits',
                 np.where(t.cos_orbits_SN15 <= 12, 'expensive 6-12 orbits',
                 np.where(np.isfinite(t.cos_orbits_SN15), 'impractical >12 orbits', 'no GALEX FUV'))))

# ---------------- MAST ----------------
# MAST: use the spectra_ism inventory already on disk (CAOM ObsCore TAP, box RA 104-128, Dec 4-34, fetched
# 2026-11-16 by scripts/spectra_ism_mast_uv.py). A fresh astroquery.mast portal query on 2026-11-17 failed
# (proxy: remote end closed connection), so it is not repeated.
MD = ROOT / 'data/raw/spectra_ism_mast'
m = pd.concat([pd.read_csv(MD / 'mast_uv_hst.csv'), pd.read_csv(MD / 'mast_uv_fuse_iue_hut.csv')], ignore_index=True)
m = m[m.dataproduct_type.astype(str).str.contains('spectrum')]
m['filters'] = (m.em_min * 1e10).round(0).astype('Int64').astype(str) + '-' + (m.em_max * 1e10).round(0).astype('Int64').astype(str) + 'A'
for k, pts in STR.items():
    s_, _ = poly_dist(m.s_ra.values, m.s_dec.values, pts)
    inside = (m.s_dec.values >= pts[:, 1].min() - 0.5) & (m.s_dec.values <= pts[:, 1].max() + 0.5)
    m[f'off_{k}'] = np.where(inside, s_, np.nan)
m['sep_min'] = m[['off_E', 'off_W', 'off_W2']].abs().min(1)
corr = m[m.sep_min <= 0.5].copy()
keep = ['obs_collection', 'instrument_name', 'filters', 'target_name', 's_ra', 's_dec', 't_exptime', 't_min',
        'obs_id', 'sep_min', 'off_E', 'off_W', 'off_W2', 'access_url']
corr = corr[[c for c in keep if c in corr.columns]].sort_values(['s_dec'])
corr.to_csv(ROOT / 'data/derived/w4_distance_targets_mast.csv', index=False)
mt = []
for i, r in t.iterrows():
    d = np.hypot((m.s_ra - r.ra) * np.cos(np.radians(r.dec)), m.s_dec - r.dec) * 3600
    hit = m[d < 5]
    mt.append('; '.join(sorted({f'{a}/{b}/{c}' for a, b, c in zip(hit.obs_collection, hit.instrument_name, hit.filters)})))
t['mast_uv_spectra'] = mt

cols = ['source_id', 'ra', 'dec', 'l', 'b', 'phot_g_mean_mag', 'bp_rp', 'M_G', 'class', 'parallax', 'parallax_error',
        'dist_pc', 'dist_lo', 'dist_hi', 'ruwe', 'teff_gspphot', 'nearest', 'sep_min', 'off_E', 'off_W', 'off_W2',
        'FUV', 'e_FUV', 'NUV', 'e_NUV', 'galex_sep_arcsec', 'cos_t_s_SN15', 'cos_orbits_SN15', 'cos_global_cps',
        'cos_flag', 'cross_ids', 'mast_uv_spectra']
t = t[cols]
t.to_csv(ROOT / 'data/derived/w4_distance_targets.csv', index=False, float_format='%.5g')
prim = t[(t.sep_min <= 0.3) & (t.dist_pc >= 120) & (t.dist_pc <= 800)]
summ = {'n_candidates_0.5deg_100_1000pc': int(len(t)), 'n_primary_0.3deg_120_800pc': int(len(prim)),
        'primary_by_class': prim['class'].value_counts().to_dict(),
        'primary_with_FUV': int(np.isfinite(prim.FUV).sum()),
        'primary_feasible_le5orbits': int((prim.cos_orbits_SN15 <= 5).sum()),
        'mast_uv_spectra_in_corridor_0.5deg': int(len(corr)),
        'mast_corridor_targets': sorted(set(corr.target_name.astype(str))),
        'cos_model': {'A_eff_cm2': AEFF, 'resel_A': RESEL, 'lambda_A': LAM, 'S/N': 15,
                      'note': 'order-of-magnitude; verify with the HST/COS ETC (sensitivity declines with time; LP-dependent)'},
        'date': '2026-11-17'}
json.dump(summ, open(ROOT / 'data/derived/w4_distance_targets_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
pd.set_option('display.width', 250)
print(prim[['source_id', 'phot_g_mean_mag', 'bp_rp', 'M_G', 'class', 'dist_pc', 'nearest', 'sep_min', 'off_E', 'FUV',
            'cos_orbits_SN15', 'cos_flag', 'cross_ids', 'mast_uv_spectra']].to_string())
