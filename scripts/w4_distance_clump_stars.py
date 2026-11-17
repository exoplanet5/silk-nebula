"""w4_distance: stellar extinction vs distance toward PGCC G199.09+19.95 (on-clump) and its surroundings (off).

Data (data/raw/w4_distance_cat/, fetched by w4_distance_fetch.py and logged):
  gaia_pgcc_r1.0.csv              Gaia DR3 (Data Lab) box +-1 deg, G < 17.8, GSP-Phot A0 / E(BP-RP) / distance
  starhorse_pgcc_r0.4.tsv          StarHorse 2021 (Anders+2022, VizieR I/354) cone 0.4 deg on the clump
  starhorse_pgccctl_{N,S,E,W}_r0.4.tsv  four 0.4-deg control cones centred 1 deg N/S/E/W of the clump
On-clump: r < R_ON from the PGCC position; off: 0.6-1.4 deg (Gaia) / 0.25-0.4 deg annulus of the same StarHorse cone
(the four StarHorse control cones came back empty from the VizieR mirror). Also an SFD-tomography regression. The SFD excess at the clump is dE(B-V) = 0.034 (A_V ~ 0.1), measured in w4_distance_dust.py.
Method: Gaia parallax distances (plx/err > 5, RUWE < 1.4); per distance bin median extinction on and off, bootstrap
errors; a step model  A_on(d) - A_off(d) = dA * H(d - d_c) fitted by a grid over d_c (100-800 pc) with free dA >= 0.
Outputs: data/derived/w4_distance_clump_stars.json, figures produced by w4_distance_figs.py
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_clump_stars.py
"""
import json, io
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/w4_distance_cat'
PG = SkyCoord(115.2847 * u.deg, 20.8334 * u.deg)
EDGES = np.array([50, 100, 150, 200, 250, 300, 375, 450, 600, 800, 1100, 1600])
rng = np.random.default_rng(1)


def read_viz(fn):
    L = [l for l in open(fn).read().splitlines() if l and not l.startswith('#')]
    if len(L) < 4:
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO('\n'.join([L[0]] + L[3:])), sep='\t')
    for c in df.columns:
        if c not in ('FlagPhoto',):
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def binned(d, a, edges, nboot=400):
    med, err, n = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (d >= lo) & (d < hi) & np.isfinite(a)
        x = a[m]
        n.append(int(m.sum()))
        if m.sum() >= 4:
            med.append(float(np.median(x)))
            bs = [np.median(rng.choice(x, len(x))) for _ in range(nboot)]
            err.append(float(np.std(bs)))
        else:
            med.append(np.nan); err.append(np.nan)
    return np.array(med), np.array(err), np.array(n)


def step_fit(dc_grid, d_on, a_on, d_off, a_off):
    """chi2 of on-star extinctions vs off-relation + step dA at dc; off relation = running median of off stars."""
    od = np.sort(d_off); oa = a_off[np.argsort(d_off)]
    # smooth off relation: median in log-distance windows
    grid = np.logspace(np.log10(60), np.log10(1500), 60)
    ref = np.array([np.nanmedian(oa[(od > g / 1.15) & (od < g * 1.15)]) if ((od > g / 1.15) & (od < g * 1.15)).sum() >= 5 else np.nan for g in grid])
    ok = np.isfinite(ref)
    base = np.interp(d_on, grid[ok], ref[ok])
    resid = a_on - base
    sig = 1.4826 * np.median(np.abs(resid - np.median(resid)))      # robust per-star scatter
    out = []
    for dc in dc_grid:
        far = d_on >= dc
        if far.sum() < 3:
            out.append((dc, np.nan, np.nan, np.nan)); continue
        dA = max(0.0, float(np.median(resid[far]) - np.median(resid[~far])) if (~far).sum() >= 3 else float(np.median(resid[far])))
        model = np.where(far, dA, 0.0) + (np.median(resid[~far]) if (~far).sum() >= 3 else 0.0)
        chi2 = float(np.sum(((resid - model) / sig) ** 2))
        out.append((dc, dA, chi2, int(far.sum())))
    return np.array(out, dtype=float), float(sig), resid


res = {'date': '2026-11-17'}
# ---------------- Gaia GSP-Phot ----------------
g = pd.read_csv(RAW / 'gaia_pgcc_r1.0.csv')
g = g[np.isfinite(g.parallax) & (g.parallax_over_error > 5) & (g.ruwe < 1.4)].copy()
g['d'] = 1000 / g.parallax
g['sep'] = SkyCoord(g.ra.values * u.deg, g.dec.values * u.deg).separation(PG).deg
gres = {}
for R_ON in (0.15, 0.25):
    on = g[(g.sep < R_ON)]; off = g[(g.sep > 0.6) & (g.sep < 1.4)]
    out = {'n_on': int(len(on)), 'n_off': int(len(off))}
    for col in ('azero_gspphot', 'ebpminrp_gspphot'):
        mo, eo, no = binned(on.d.values, on[col].values, EDGES)
        mf, ef, nf = binned(off.d.values, off[col].values, EDGES)
        out[col] = {'edges_pc': EDGES.tolist(), 'on_median': mo.round(4).tolist(), 'on_err': eo.round(4).tolist(), 'on_n': no.tolist(),
                    'off_median': mf.round(4).tolist(), 'off_err': ef.round(4).tolist(), 'off_n': nf.tolist(),
                    'diff': (mo - mf).round(4).tolist(), 'diff_err': np.hypot(eo, ef).round(4).tolist()}
        m_on = np.isfinite(on[col].values); m_off = np.isfinite(off[col].values)
        sf, sig, _ = step_fit(np.arange(100, 801, 20), on.d.values[m_on], on[col].values[m_on], off.d.values[m_off], off[col].values[m_off])
        if np.isfinite(sf[:, 2]).any():
            k = np.nanargmin(sf[:, 2])
            dchi = sf[:, 2] - sf[k, 2]
            okk = dchi <= 1.0
            out[col]['step_fit'] = {'best_dc_pc': float(sf[k, 0]), 'dA': float(sf[k, 1]), 'n_far': int(sf[k, 3]),
                                    'dc_range_dchi2<=1': [float(sf[okk, 0].min()), float(sf[okk, 0].max())],
                                    'per_star_sigma': sig,
                                    'chi2_at': {str(int(x)): round(float(c - sf[k, 2]), 2) for x, c in zip(sf[:, 0], sf[:, 2]) if int(x) in (160, 200, 250, 282, 320, 400, 500, 600)}}
    gres[f'R_on_{R_ON}'] = out
res['gaia_gspphot'] = gres

# ---------------- StarHorse (on-clump cone only; the four control cones came back empty from the mirror) ----
fon = RAW / 'starhorse_pgcc_r0.4.tsv'
s_on = read_viz(fon)
s_on = s_on[(s_on.fidelity > 0.5) & (s_on.FlagOut == 0)].copy()
s_on['d'] = s_on.dist50 * 1000
s_on['dfrac'] = (s_on.dist84 - s_on.dist16) / 2 / s_on.dist50
s_on['sep'] = SkyCoord(s_on.RA_ICRS.values * u.deg, s_on.DE_ICRS.values * u.deg).separation(PG).deg
s_on = s_on[s_on.dfrac < 0.25]
sres = {'n_stars_cone0.4_clean': int(len(s_on)), 'cuts': 'fidelity > 0.5, FlagOut == 0, (d84-d16)/2/d50 < 0.25'}
for R_ON in (0.12, 0.2):
    on = s_on[s_on.sep < R_ON]; off = s_on[(s_on.sep > 0.25) & (s_on.sep < 0.4)]
    mo, eo, no = binned(on.d.values, on.AV50.values, EDGES)
    mf, ef, nf = binned(off.d.values, off.AV50.values, EDGES)
    sres[f'R_on_{R_ON}_off_0.25-0.4'] = {'n_on': int(len(on)), 'n_off': int(len(off)), 'edges_pc': EDGES.tolist(),
                            'on_median_AV': mo.round(4).tolist(), 'on_err': eo.round(4).tolist(), 'on_n': no.tolist(),
                            'off_median_AV': mf.round(4).tolist(), 'off_err': ef.round(4).tolist(), 'off_n': nf.tolist(),
                            'diff': (mo - mf).round(4).tolist(), 'diff_err': np.hypot(eo, ef).round(4).tolist()}
res['starhorse'] = sres

# ---------------- SFD tomography: regress stellar extinction on SFD E(B-V) at the star, per distance bin ------
from astropy.io import fits
from astropy.wcs import WCS
h = fits.open(ROOT / 'data/raw/skyview_imaging_mw/eq_SFD_Dust_Map.fits')[0]
W_ = WCS(h.header)


def sfd_at(ra, dec):
    x, y = W_.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    xi, yi = np.round(x).astype(int), np.round(y).astype(int)
    ok = (xi >= 0) & (yi >= 0) & (xi < h.data.shape[1]) & (yi < h.data.shape[0])
    out = np.full(len(ra), np.nan); out[ok] = h.data[yi[ok], xi[ok]]
    return out


TB = np.array([50, 150, 250, 350, 500, 800, 1500])


def tomo(d, A, E, rmask, nboot=500):
    rows_ = []
    for lo, hi in zip(TB[:-1], TB[1:]):
        m = (d >= lo) & (d < hi) & np.isfinite(A) & np.isfinite(E) & rmask
        if m.sum() < 8:
            rows_.append({'d_lo': int(lo), 'd_hi': int(hi), 'n': int(m.sum())}); continue
        x, y = E[m], A[m]
        k = np.polyfit(x, y, 1)[0]
        bs = []
        for _ in range(nboot):
            j = rng.integers(0, len(x), len(x))
            if np.ptp(x[j]) > 0:
                bs.append(np.polyfit(x[j], y[j], 1)[0])
        rows_.append({'d_lo': int(lo), 'd_hi': int(hi), 'n': int(m.sum()), 'slope': round(float(k), 3),
                      'slope_err': round(float(np.std(bs)), 3), 'sfd_range': [round(float(x.min()), 3), round(float(x.max()), 3)],
                      'n_sfd_gt_0.065': int((x > 0.065).sum())})
    return rows_


g['sfd'] = sfd_at(g.ra.values, g.dec.values)
s_on['sfd'] = sfd_at(s_on.RA_ICRS.values, s_on.DE_ICRS.values)
res['sfd_tomography'] = {
    'note': 'slope k = dA/dE(B-V)_SFD of stars in each distance bin; ~0 in front of the dust, ~2.8-3.1 (A_V) or ~1.3 (E(BP-RP)) behind it. Only stars within the stated radius of the PGCC, where the clump dominates the SFD contrast (0.097 peak vs 0.048 background).',
    'gaia_A0_r<0.5': tomo(g.d.values, g.azero_gspphot.values, g.sfd.values, (g.sep < 0.5).values),
    'gaia_EBPRP_r<0.5': tomo(g.d.values, g.ebpminrp_gspphot.values, g.sfd.values, (g.sep < 0.5).values),
    'starhorse_AV_r<0.4': tomo(s_on.d.values, s_on.AV50.values, s_on.sfd.values, np.ones(len(s_on), bool))}
json.dump(res, open(ROOT / 'data/derived/w4_distance_clump_stars.json', 'w'), indent=1)


def show(d):
    print(json.dumps(d, indent=0)[:3500])


for k, v in res['gaia_gspphot'].items():
    print('GAIA', k, v['n_on'], v['n_off'])
    for col in ('azero_gspphot',):
        print(' diff', v[col]['diff']); print(' err ', v[col]['diff_err']); print(' n_on', v[col]['on_n'])
        print(' step', v[col].get('step_fit'))
for k, v in res['starhorse'].items():
    if k.startswith('R_on'):
        print('SH', k, v['n_on'], v['n_off']); print(' diff', v['diff']); print(' err ', v['diff_err']); print(' n_on', v['on_n'])
for k, v in res['sfd_tomography'].items():
    if k != 'note':
        print('TOMO', k)
        for r in v:
            print('   ', r)
