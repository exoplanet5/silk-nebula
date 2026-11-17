"""Hot ionising sources around the Silk Nebula: hot WDs, hot subdwarfs, PN nuclei, O/B stars.

Inputs (data/raw/context_catalogs/, fetched by context_fetch_catalogs.py, accessed 2026-11-16):
  simbad_hotstars_cone25.csv  SIMBAD (CfA mirror) objects with O/B0-2/sdO/sdB/DO/DAO/PG1159/WR/DA types or WD/HS/PN otypes
  gf21_wd_cone20_hot.tsv      Gentile Fusillo+21 WDs, Pwd>0.75, plx>1 mas, Teff_H>20 kK (Gaia photometric Teff)
  culpan22_knownhsd_cone20.tsv Culpan+22 known hot subdwarfs with spectroscopic Teff/logg (3" match; used for sd stars)
  (geier20_sdcat_cone20.tsv and culpan22_hotsd_cone20.tsv were fetched but not used: Geier 2020 is superseded by
   the Culpan+22 known-hsd table; the Culpan+22 Gaia candidate list has no Teff)
Literature parameters added by hand (spectroscopic): WD J0736+1444 = UCAC4 524-042847 (Reindl+2023 A&A 677 A29,
  Table A.4: Teff 81319+-1599 K, log g 7.92, R 0.0148 Rsun, L 8.7 Lsun); SIMBAD mesFe_h values for the others.
Outputs: data/derived/context_hotstars.csv   (candidates within 15 deg; Q(H), Q(>35 eV), Stromgren radii, energy budget)
         data/derived/context_hotstars_top.csv (ranked shortlist)
Method: blackbody photon rates (context_common.q_bb); R from log g + mass (WD M=0.6, sd M=0.47 Msun unless given)
        or from luminosity; Stromgren radius for n = 0.1, 0.3, 1 cm^-3 (alpha_B=2.59e-13).
        Energy-budget test: maximum [O III] 5007 photon rate the star could power on the filament,
        N_max = Q(H) * f_cov * eta, with covering factor f_cov = solid angle of a 1-deg-wide x L strip at the
        filament distance seen from the star (computed from the projected separation) and eta = 2 [O III] photons
        per ionisation (generous upper bound: ~10 eV heat per photoionisation at T*~80-100 kK, ~50% radiated in
        5007), compared with the required N_req = 1e6 [ph cm^-2 s^-1 per R, column emission rate into 4 pi] * I[R] * A_proj(d)
        (1 R = 1e6/(4 pi) ph cm^-2 s^-1 sr^-1; an earlier version used 1e7 by mistake, fixed 2026-11-16).
        Also reported: budget_ratio_allsky = 2 Q(H) / N_req (f_cov = 1, absolute upper bound).
        Wake test: LSR-frame proper-motion direction vs local filament position angle.
"""
import sys, json, re
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *

RAW = ROOT / 'data/raw/context_catalogs'
G = 6.674e-8; MSUN = 1.989e33
I_OIII_R = (10.0, 13.0)          # Wu+2026 Table 3: 23.6-23.9 mag/arcsec^2 ~ 10-13 R
AREA_DEG2 = (3.0, 10.0)          # bright-filament area range: ~0.3-1 deg x 10 deg


def trace_line():
    """Dense polyline of the filament (refined trace if available, else provisional NSNS arc + P points)."""
    ra, dec, sig, lab, src = trace_points()
    o = np.argsort(dec)
    ra, dec = ra[o], dec[o]
    dd = np.linspace(dec.min(), dec.max(), 400)
    return np.interp(dd, dec, ra), dd, src


def min_sep_to_trace(ra, dec, tra, tdec):
    s = angsep(np.asarray(ra)[:, None], np.asarray(dec)[:, None], tra[None, :], tdec[None, :])
    k = np.argmin(s, 1)
    return s[np.arange(len(ra)), k], k


def local_pa(tra, tdec, k):
    k0, k1 = max(k - 5, 0), min(k + 5, len(tra) - 1)
    dra = (tra[k1] - tra[k0]) * np.cos(np.radians(tdec[k]))
    ddec = tdec[k1] - tdec[k0]
    return np.degrees(np.arctan2(dra, ddec)) % 180


def teff_from_sptype(sp):
    sp = str(sp)
    m = re.match(r'^O(\d(\.\d)?)', sp)
    if m:   # Martins+2005 observational scale, dwarfs (approx.)
        sub = float(m.group(1));
        return float(np.interp(sub, [3, 4, 5, 6, 7, 8, 9, 9.5], [44600, 42900, 41500, 38900, 36900, 34900, 31900, 31400]))
    m = re.match(r'^B(\d(\.\d)?)', sp)
    if m:
        sub = float(m.group(1)); return float(np.interp(sub, [0, 1, 2, 3], [30000, 25400, 22000, 18700]))
    if sp.startswith(('DO', 'PG1159', 'O(He)')): return 70000.0
    if sp.startswith('DAO'): return 60000.0
    if sp.startswith('sdO'): return 45000.0
    if sp.startswith('sdB'): return 28000.0
    return np.nan


def main():
    tra, tdec, tsrc = trace_line()
    cen = CENTER
    s = pd.read_csv(RAW / 'simbad_hotstars_cone25.csv')
    s['sep_c'] = angsep(s.ra.values, s.dec.values, *cen)
    s['sp'] = s.sp_type.fillna('')
    hot_sp = s.sp.str.match(r'^(O|B0|B1|B2|sdO|DO|DAO|PG ?1159|WR|\[WC|DOZ|O\(He\))') | (s.otype == 'PN')
    s = s[hot_sp & (s.sep_c < 15)].copy()
    s['src'] = 'SIMBAD'
    # GF21 hot WDs (photometric Teff); keep those not already in SIMBAD list (match 3")
    w = pd.read_csv(RAW / 'gf21_wd_cone20_hot.tsv', sep='\t', comment='#')
    w = w.iloc[2:].copy()                      # drop VizieR unit and dash rows
    w.columns = [c.strip() for c in w.columns]
    w['WDJname'] = w['WDJname'].astype(str).str.strip()
    w = w[pd.to_numeric(w._r, errors='coerce') < 15].copy()
    for c in ['RA_ICRS', 'DE_ICRS', 'Plx', 'e_Plx', 'pmRA', 'pmDE', 'TeffH', 'loggH', 'MassH', 'TeffHe', 'GMAG', 'Gmag']:
        w[c] = pd.to_numeric(w[c], errors='coerce')
    rows = []
    for _, r in s.iterrows():
        rows.append(dict(name=r.main_id, ra=r.ra, dec=r.dec, otype=r.otype, sptype=r.sp, plx=r.plx_value,
                         e_plx=r.plx_err, pmra=r.pmra, pmdec=r.pmdec, rv=r.rvz_radvel, vmag=r.vmag, cat='SIMBAD'))
    sc = SkyCoord(s.ra.values * u.deg, s.dec.values * u.deg)
    wc = SkyCoord(w.RA_ICRS.values * u.deg, w.DE_ICRS.values * u.deg)
    if len(sc):
        idx, d2d, _ = wc.match_to_catalog_sky(sc)
        new = d2d.arcsec > 3
    else:
        new = np.ones(len(w), bool)
    for (_, r), isnew in zip(w.iterrows(), new):
        if not isnew:
            continue
        rows.append(dict(name=r.WDJname, ra=r.RA_ICRS, dec=r.DE_ICRS, otype='WD(GF21)', sptype='', plx=r.Plx,
                         e_plx=r.e_Plx, pmra=r.pmRA, pmdec=r.pmDE, rv=np.nan, vmag=np.nan, cat='GF21',
                         teff_phot=r.TeffH, logg_phot=r.loggH, mass_phot=r.MassH))
    df = pd.DataFrame(rows)
    df['dist_pc'] = 1000.0 / df.plx.where(df.plx > 0)
    # spectroscopic Teff/logg of known hot subdwarfs (Culpan+22 known-hsd table, VizieR J/A+A/662/A40/knownhsd), 3" match
    kh = pd.read_csv(RAW / 'culpan22_knownhsd_cone20.tsv', sep='\t', comment='#', low_memory=False).iloc[2:].copy()
    kh.columns = [c.strip() for c in kh.columns]
    for c in ['RA_ICRS', 'DE_ICRS', 'Teff', 'logg']:
        kh[c] = pd.to_numeric(kh[c], errors='coerce')
    kc = SkyCoord(kh.RA_ICRS.values * u.deg, kh.DE_ICRS.values * u.deg)
    dc = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg)
    ik, d2k, _ = dc.match_to_catalog_sky(kc)
    okk = d2k.arcsec < 3
    df['c_teff'] = np.where(okk, kh.Teff.values[ik], np.nan)
    df['c_logg'] = np.where(okk, kh.logg.values[ik], np.nan)
    # spectroscopic parameters from SIMBAD mesFe_h / literature (hand-entered, see docstring)
    lit = {
        'UCAC4 524-042847': dict(teff=81319, e_teff=1599, logg=7.92, r_rsun=0.0148, mass=0.67, ref='Reindl+2023 A&A 677 A29 (J0736+1444)'),
        'GALEX J074507.5+144939': dict(teff=80025, logg=8.006, mass=0.6, ref='Vincent+2024 A&A 682 A5 (Gaia XP)'),
        'GALEX J074320.3+160752': dict(teff=83230, logg=7.15, mass=0.6, ref='Bedard+2020 ApJ 901 93 (also 100 kK, Kepler+15)'),
        'GALEX J075633.8+163358': dict(teff=56820, logg=7.93, mass=0.6, ref='Bedard+2020 ApJ 901 93'),
        'GALEX J075234.2+161604': dict(teff=50760, logg=5.72, mass=0.47, ref='Lei+2024 ApJS 271 21'),
        'HD  60848': dict(teff=37034, logg=4.03, r_rsun=7.0, ref='SIMBAD mesFe_h (2024A&A...681A.107R); R~7 Rsun O8V assumed'),
    }
    out = []
    for _, r in df.iterrows():
        rec = r.to_dict()
        L = lit.get(r['name'], {})
        if L:
            rec.update(teff=L['teff'], logg=L.get('logg'), mass=L.get('mass', 0.6), teff_ref=L['ref'], r_rsun=L.get('r_rsun'))
        elif np.isfinite(r.get('c_teff', np.nan)) and np.isfinite(r.get('c_logg', np.nan)):
            rec.update(teff=float(r.c_teff), logg=float(r.c_logg), mass=0.47, teff_ref='Culpan+22 known hsd (spectroscopic)')
        elif r['cat'] == 'GF21':
            rec.update(teff=r.teff_phot, logg=r.logg_phot, mass=r.mass_phot if np.isfinite(r.mass_phot) else 0.6,
                       teff_ref='GF21 Gaia photometric Teff_H (unreliable >~30 kK)')
        else:
            rec.update(teff=teff_from_sptype(r.sptype), logg=np.nan, mass=np.nan, teff_ref='spectral-type calibration')
        sp = str(r.sptype)
        if not rec.get('r_rsun') or not np.isfinite(rec.get('r_rsun') or np.nan):
            if np.isfinite(rec.get('logg') or np.nan):
                m = rec.get('mass') if np.isfinite(rec.get('mass') or np.nan) else (0.47 if sp.startswith('sd') else 0.6)
                rec['r_rsun'] = np.sqrt(G * m * MSUN / 10 ** rec['logg']) / RSUN
            elif sp.startswith(('O', 'B')):
                rec['r_rsun'] = 7.0 if sp.startswith('O') else 4.5
            elif sp.startswith(('sd',)):
                rec['r_rsun'] = 0.15
            elif sp.startswith(('DO', 'DA', 'PG')) or 'WD' in str(r.otype):
                rec['r_rsun'] = 0.015
            else:
                rec['r_rsun'] = np.nan
        T, R = rec['teff'], rec['r_rsun']
        if np.isfinite(T) and np.isfinite(R):
            rec['Q_H'] = q_bb(T, R)
            rec['Q_35eV'] = q_bb(T, R, 35.12)
            rec['L_bol_lsun'] = 4 * np.pi * (R * RSUN) ** 2 * 5.6704e-5 * T ** 4 / 3.828e33
            for n in (0.1, 0.3, 1.0):
                rec[f'Rs_pc_n{n}'] = stromgren_pc(rec['Q_H'], n)
        out.append(rec)
    df = pd.DataFrame(out)
    sep, k = min_sep_to_trace(df.ra.values, df.dec.values, tra, tdec)
    df['sep_trace_deg'] = sep
    df['sep_silk_centre_deg'] = angsep(df.ra.values, df.dec.values, *cen)
    df['l'], df['b'] = gal(df.ra.values, df.dec.values)
    df['d_proj_pc'] = df.dist_pc * np.radians(df.sep_trace_deg)
    df['Rs01_over_dproj'] = df['Rs_pc_n0.1'] / df.d_proj_pc
    # energy budget at the star's distance
    d = df.dist_pc
    A_lo = AREA_DEG2[0] * (np.radians(1) * d * 3.0857e18) ** 2
    A_hi = AREA_DEG2[1] * (np.radians(1) * d * 3.0857e18) ** 2
    df["N_OIII_req_lo"] = 1e6 * I_OIII_R[0] * A_lo
    df["N_OIII_req_hi"] = 1e6 * I_OIII_R[1] * A_hi
    # covering factor: strip of length L=10 deg and width 1 deg at distance d, seen from projected sep s (same distance)
    Lpc = np.radians(10) * d; Wpc = np.radians(1) * d; s_pc = df.d_proj_pc.clip(lower=0.5)
    df['f_cov'] = np.clip((Wpc * 2 * np.arctan(Lpc / (2 * s_pc))) / (4 * np.pi * s_pc), 0, 0.5)
    df['N_OIII_max'] = df.Q_H * df.f_cov * 2.0
    df['budget_ratio_max_over_req'] = df.N_OIII_max / df.N_OIII_req_lo
    df['budget_ratio_allsky'] = 2.0 * df.Q_H / df.N_OIII_req_lo          # f_cov = 1 upper bound
    # wake test (LSR-frame proper motion vs local filament PA)
    has_pm = np.isfinite(df.pmra) & np.isfinite(df.pmdec) & np.isfinite(df.dist_pc)
    pa_pm = np.full(len(df), np.nan); mu_lsr = np.full(len(df), np.nan)
    for i in np.where(has_pm)[0]:
        r = df.iloc[i]
        a, b_ = lsr_pm(r.ra, r.dec, r.dist_pc, r.pmra, r.pmdec, r.rv if np.isfinite(r.rv) else 0.0)
        pa_pm[i] = np.degrees(np.arctan2(a[0], b_[0])) % 360
        mu_lsr[i] = np.hypot(a[0], b_[0])
    df['pm_lsr_masyr'] = mu_lsr
    df['pa_pm_lsr'] = pa_pm
    df['pa_trace_local'] = [local_pa(tra, tdec, kk) for kk in k]
    dpa = np.abs(((df.pa_pm_lsr - df.pa_trace_local) + 90) % 180 - 90)
    df['dpa_pm_vs_trace'] = dpa
    df['v_t_lsr_kms'] = K_PM * df.pm_lsr_masyr * df.dist_pc / 1000
    # wake test 2: LSR-frame back-track (v_r = RV or 0) over 0-1 Myr; fraction of spine points within 0.5 deg of it
    tgrid = np.linspace(0, 1.0, 401)
    sp = load_trace()['spine'] if load_trace() and not load_trace().get('preliminary', True) else None
    spra = np.array([q['ra'] for q in sp]) if sp else tra
    spde = np.array([q['dec'] for q in sp]) if sp else tdec
    frac_w, t_w = np.full(len(df), np.nan), np.full(len(df), np.nan)
    for i in np.where(has_pm)[0]:
        r = df.iloc[i]
        lt, bt, _ = backtrack_lsr(r.ra, r.dec, r.dist_pc, r.pmra, r.pmdec, r.rv if np.isfinite(r.rv) else 0.0, tgrid)
        c = SkyCoord(lt[0] * u.deg, bt[0] * u.deg, frame='galactic').icrs
        dmin = angsep(spra[:, None], spde[:, None], c.ra.deg[None, :], c.dec.deg[None, :])
        near = dmin.min(1) < 0.5
        frac_w[i] = near.mean()
        if near.any():
            t_w[i] = float(np.median(tgrid[np.argmin(dmin[near], 1)]))
    df['wake_frac_spine_within_0p5deg_of_track'] = frac_w
    df['wake_t_cross_myr'] = t_w
    df['trace_source'] = tsrc
    cols = ['name', 'cat', 'otype', 'sptype', 'ra', 'dec', 'l', 'b', 'plx', 'e_plx', 'dist_pc', 'pmra', 'pmdec',
            'teff', 'logg', 'r_rsun', 'teff_ref', 'L_bol_lsun', 'Q_H', 'Q_35eV', 'Rs_pc_n0.1', 'Rs_pc_n0.3', 'Rs_pc_n1.0',
            'sep_silk_centre_deg', 'sep_trace_deg', 'd_proj_pc', 'Rs01_over_dproj', 'f_cov', 'N_OIII_max',
            'N_OIII_req_lo', 'N_OIII_req_hi', 'budget_ratio_max_over_req', 'budget_ratio_allsky', 'pm_lsr_masyr', 'pa_pm_lsr',
            'pa_trace_local', 'dpa_pm_vs_trace', 'v_t_lsr_kms', 'wake_frac_spine_within_0p5deg_of_track', 'wake_t_cross_myr', 'trace_source']
    df = df[cols].sort_values('sep_trace_deg')
    df.to_csv(ROOT / 'data/derived/context_hotstars.csv', index=False, float_format='%.5g')
    # shortlist: d < 1 kpc (plx > 1 mas) and T >= 40 kK
    top = df[(df.dist_pc < 1000) & (df.teff >= 40000)].copy()
    top['rank_metric'] = top.budget_ratio_max_over_req
    top = top.sort_values('rank_metric', ascending=False)
    top.to_csv(ROOT / 'data/derived/context_hotstars_top.csv', index=False, float_format='%.4g')
    pd.set_option('display.width', 250)
    print('trace:', tsrc, ' candidates within 15 deg:', len(df), ' hot (>=40kK) & d<1kpc:', len(top))
    print(top[['name', 'sptype', 'dist_pc', 'teff', 'r_rsun', 'Q_H', 'Q_35eV', 'Rs_pc_n0.1', 'Rs_pc_n1.0',
               'sep_trace_deg', 'd_proj_pc', 'f_cov', 'budget_ratio_max_over_req', 'budget_ratio_allsky', 'pm_lsr_masyr',
               'dpa_pm_vs_trace']].head(25).to_string())


if __name__ == '__main__':
    main()
