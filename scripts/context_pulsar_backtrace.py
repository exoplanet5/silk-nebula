"""Back-trace nearby pulsars (d < 1 kpc, within 30 deg of the Silk) over 0-1 Myr in the LSR frame and test whether
any birth locus coincides with the centre of curvature of the Silk [O III] filament.

Inputs : data/derived/context_pulsars_30deg.csv (scripts/context_psrcat.py, ATNF v2.9.0)
         data/derived/silk_trace.json (geometry teammate, v1.0: E-strand spine + curvature fits)
Outputs: data/derived/context_pulsar_tracks.csv       nominal tracks (v_r = 0, nominal distance), 10-kyr steps
         data/derived/context_pulsar_birth_mc.csv     MC loci at t = tau_c (every 10th sample)
         data/derived/context_curvature.json          centre-of-curvature regions used here
         data/derived/context_pulsar_coincidence.json summary statistics per pulsar
Centre-of-curvature regions (on the sky):
  A  'formal'     Delta-chi2 <= 6.18 region of a single small circle fitted to the E-strand spine with an
                  intrinsic scatter added so that chi2_min/dof = 1 (context_common.centre_region_grid);
  B  'systematic' union of the 16-84% bootstrap boxes of the geometry teammate's E_full, E_north and E_full+W
                  fits (silk_trace.json curvature.fits) -- RA 84.7-96.3, Dec 14.6-17.1 (piecewise boxes);
  the E_south (Dec<=17.5) part alone is almost straight (R ~ 46 deg, centre unconstrained) and is not used.
Assumptions: distance from parallax (Gaussian in parallax) else DM-YMW16 with 40% log-normal scatter;
v_r ~ N(0, 265 km/s) (1-D sigma of the Hobbs+2005 Maxwellian); solar motion (U,V,W) = (11.1, 12.24, 7.25) km/s
(Schoenrich+2010); the gas at the birth site co-moves with the LSR; Galactic differential rotation neglected
(Oort A*d < 5 km/s for d < 0.35 kpc). Geminga: parallax 4.0+-1.3 mas, PM (142.2,107.4)+-1.2 mas/yr
(Faherty, Walter & Anderson 2007) as baseline; ATNF v2.9.0 values (PX 6.4+-1.8, PM 138, 97 +-4; Caraveo+1996)
run as 'Geminga_ATNF'.
"""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *

N_MC = 3000
rng = np.random.default_rng(42)
SIG_VR = 265.0
T_GRID = np.round(np.arange(0, 1.0001, 0.005), 4)   # Myr


def region_B(tr):
    f = tr['curvature']['fits']
    boxes = []
    for k in ('E_full', 'E_north_ge17p75', 'E_full_plus_W'):
        b = f[k]['boot_16_84']
        boxes.append((b['center_ra'][0], b['center_ra'][1], b['center_dec'][0], b['center_dec'][1], k))
    return boxes


def in_boxes(boxes, ra, dec):
    m = np.zeros(np.shape(ra), bool)
    for r0, r1, d0, d1, _ in boxes:
        m |= (ra >= r0) & (ra <= r1) & (dec >= d0) & (dec <= d1)
    return m


def load_psr():
    df = pd.read_csv(ROOT / 'data/derived/context_pulsars_30deg.csv')
    df = df[(df.dist_best_kpc < 1.0) & np.isfinite(df.pmra) & np.isfinite(df.pmdec)].copy()
    rows = []
    for _, r in df.iterrows():
        name = r.psrb if isinstance(r.psrb, str) and r.psrb.strip() else r.psrj
        d = dict(name=name, psrj=r.psrj, ra=r.ra, dec=r.dec, pmra=r.pmra, pmdec=r.pmdec,
                 e_pmra=r.e_pmra if np.isfinite(r.e_pmra) else 5.0, e_pmdec=r.e_pmdec if np.isfinite(r.e_pmdec) else 5.0,
                 age=r.age_yr, px=r.px, e_px=r.e_px, d_kpc=r.dist_best_kpc, dsrc=r.dist_src, dm=r.dm, edot=r.edot)
        if r.psrj == 'J0633+1746':
            alt = dict(d); alt['name'] = 'Geminga_ATNF'; alt['dsrc'] = 'PX (ATNF v2.9.0 / Caraveo+96)'
            alt.update(px=6.4, e_px=1.8, pmra=138.0, pmdec=97.0, e_pmra=4.0, e_pmdec=4.0)
            d.update(name='Geminga', px=4.0, e_px=1.3, pmra=142.2, pmdec=107.4, e_pmra=1.2, e_pmdec=1.2,
                     dsrc='PX (Faherty+07)')
            rows += [d, alt]
            continue
        if r.psrj == 'J0659+1414':
            d.update(px=3.47, e_px=0.36, dsrc='PX (Brisken+03)')
        rows.append(d)
    return rows


def draw_dist(p, n):
    if np.isfinite(p['px']) and np.isfinite(p['e_px']) and p['px'] / p['e_px'] > 2:
        pxs = rng.normal(p['px'], p['e_px'], n)
        return 1000.0 / np.clip(pxs, 0.4 * p['px'], None)
    return 1000.0 * p['d_kpc'] * np.exp(rng.normal(0, 0.4, n))


def main():
    tr = load_trace()
    assert tr and not tr.get('preliminary', True), 'refined silk_trace.json required'
    cv = tr['curvature']
    nom = (cv['center_ra'], cv['center_dec'], cv['radius_deg'])
    regA = centre_region_grid(strand='E')
    boxes = region_B(tr)
    GR, GD = np.meshgrid(regA['ra'], regA['dec'])
    mA = regA['dchi2'] <= 6.18
    curv = dict(trace_version=tr.get('version'), nominal_geometry_fit=dict(ra=nom[0], dec=nom[1], radius_deg=nom[2],
                l=cv['center_l'], b=cv['center_b'], boot_16_84=dict(ra=cv['center_ra_16_84'], dec=cv['center_dec_16_84'],
                radius=cv['radius_16_84'])),
                regionA_formal=dict(best=[float(x) for x in regA['best']], sigma_int_deg=regA['sigma_int'],
                                    ra_range=[float(GR[mA].min()), float(GR[mA].max())],
                                    dec_range=[float(GD[mA].min()), float(GD[mA].max())],
                                    radius_range=[float(regA['radius'][mA].min()), float(regA['radius'][mA].max())]),
                regionB_systematic_boxes=[dict(ra=[b[0], b[1]], dec=[b[2], b[3]], fit=b[4]) for b in boxes],
                caveat=cv.get('caveat', ''))
    (ROOT / 'data/derived/context_curvature.json').write_text(json.dumps(curv, indent=1, default=float))
    np.savez(ROOT / 'data/derived/context_curvature_regionA.npz', ra=regA['ra'], dec=regA['dec'], dchi2=regA['dchi2'],
             radius=regA['radius'])

    psrs = load_psr()
    track_rows, birth_rows, summary = [], [], {}
    for p in psrs:
        d0 = 1000.0 / p['px'] if np.isfinite(p['px']) and p['px'] / p['e_px'] > 2 else 1000.0 * p['d_kpc']
        l, b, d = backtrack_lsr(p['ra'], p['dec'], d0, p['pmra'], p['pmdec'], 0.0, T_GRID)
        ct = SkyCoord(l[0] * u.deg, b[0] * u.deg, frame='galactic').icrs
        for k, t in enumerate(T_GRID):
            if k % 2 == 0:
                track_rows.append(dict(name=p['name'], t_myr=t, l=l[0, k], b=b[0, k], d_pc=d[0, k],
                                       ra=ct.ra.deg[k], dec=ct.dec.deg[k]))
        n = N_MC
        dd = draw_dist(p, n)
        pmra = rng.normal(p['pmra'], p['e_pmra'], n)
        pmdec = rng.normal(p['pmdec'], p['e_pmdec'], n)
        vr = rng.normal(0, SIG_VR, n)
        tau = p['age'] / 1e6 if np.isfinite(p['age']) else np.nan
        lt, bt, dt = backtrack_lsr(np.full(n, p['ra']), np.full(n, p['dec']), dd, pmra, pmdec, vr, T_GRID)
        allc = SkyCoord(lt.ravel() * u.deg, bt.ravel() * u.deg, frame='galactic').icrs
        rr, de = allc.ra.deg.reshape(n, -1), allc.dec.deg.reshape(n, -1)
        sep_nom = angsep(rr, de, nom[0], nom[1])
        inA = in_centre_region(regA, rr.ravel(), de.ravel()).reshape(n, -1)
        inB = in_boxes(boxes, rr, de)
        itau = int(np.argmin(np.abs(T_GRID - tau))) if np.isfinite(tau) and tau <= 1 else None
        for i in range(0, n, 10):
            if itau is not None:
                birth_rows.append(dict(name=p['name'], t_myr=T_GRID[itau], l=lt[i, itau], b=bt[i, itau], d_pc=dt[i, itau],
                                       ra=rr[i, itau], dec=de[i, itau], d_now=dd[i], vr=vr[i]))
        tmin = T_GRID[np.argmin(sep_nom, 1)]
        tau_win = (T_GRID >= 0.5 * tau) & (T_GRID <= 1.5 * tau) if np.isfinite(tau) else np.zeros_like(T_GRID, bool)
        pct = lambda a: np.percentile(a, [16, 50, 84]).round(2).tolist()
        s = dict(psrj=p['psrj'], tau_c_kyr=None if not np.isfinite(tau) else round(tau * 1e3, 1), edot=p['edot'],
                 dist_source=p['dsrc'], d_now_pc_p16_50_84=np.percentile(dd, [16, 50, 84]).round(0).tolist(),
                 pm_masyr=[p['pmra'], p['pmdec']],
                 sep_now_from_silk_centre_deg=round(float(angsep(p['ra'], p['dec'], *CENTER)), 2),
                 sep_now_from_curv_centre_deg=round(float(angsep(p['ra'], p['dec'], nom[0], nom[1])), 2),
                 min_sep_track_from_nominal_centre_deg_p16_50_84=pct(sep_nom.min(1)),
                 t_closest_myr_p16_50_84=pct(tmin),
                 frac_tracks_entering_regionA_0_1Myr=float(inA.any(1).mean()),
                 frac_tracks_entering_regionB_0_1Myr=float(inB.any(1).mean()),
                 frac_in_regionB_within_0p5_1p5_tau=float((inB[:, tau_win]).any(1).mean()) if tau_win.any() else None,
                 frac_in_regionA_within_0p5_1p5_tau=float((inA[:, tau_win]).any(1).mean()) if tau_win.any() else None)
        if itau is not None:
            s['birth_at_tau'] = dict(t_myr=float(T_GRID[itau]), ra_p16_50_84=pct(rr[:, itau]), dec_p16_50_84=pct(de[:, itau]),
                                     l_p16_50_84=pct(lt[:, itau]), b_p16_50_84=pct(bt[:, itau]),
                                     d_pc_p16_50_84=np.percentile(dt[:, itau], [16, 50, 84]).round(0).tolist(),
                                     sep_from_nominal_centre_deg_p16_50_84=pct(sep_nom[:, itau]))
        # 3-D: if the Silk were the limb of a sphere centred on the birth site (radius rho_nom on the sky)
        if itau is not None:
            Db = dt[:, itau]
            sp = tr['spine']
            rho_b = np.array([np.mean(angsep(np.array([q['ra'] for q in sp]), np.array([q['dec'] for q in sp]), a, b_))
                              for a, b_ in zip(rr[::30, itau], de[::30, itau])])
            rms_b = np.array([np.std(angsep(np.array([q['ra'] for q in sp]), np.array([q['dec'] for q in sp]), a, b_))
                              for a, b_ in zip(rr[::30, itau], de[::30, itau])])
            Dbs = Db[::30]
            s['limb_shell_if_centred_at_birth'] = dict(
                rho_mean_sep_spine_deg_p16_50_84=pct(rho_b), rms_sep_spine_deg_p16_50_84=pct(rms_b),
                R_pc_p16_50_84=np.percentile(Dbs * np.sin(np.radians(rho_b)), [16, 50, 84]).round(0).tolist(),
                silk_dist_pc_p16_50_84=np.percentile(Dbs * np.cos(np.radians(rho_b)), [16, 50, 84]).round(0).tolist())
        summary[p['name']] = s
        print(p['name'], json.dumps(s))
    pd.DataFrame(track_rows).to_csv(ROOT / 'data/derived/context_pulsar_tracks.csv', index=False, float_format='%.4f')
    pd.DataFrame(birth_rows).to_csv(ROOT / 'data/derived/context_pulsar_birth_mc.csv', index=False, float_format='%.4f')
    out = dict(curvature=curv, pulsars=summary,
               assumptions=dict(sigma_vr_kms=SIG_VR, sun_uvw=SUN_UVW.tolist(), n_mc=N_MC, t_grid_myr=[0, 1, 0.005],
                                dm_dist_lognormal_sigma=0.4, frame='LSR-comoving, present Sun at origin'))
    (ROOT / 'data/derived/context_pulsar_coincidence.json').write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
