"""w3_northsouth_geom.py -- wave-3 (tag northsouth): geometry and kinematics of the north/south difference.

1. Spine distance theta(Dec) from candidate centres (Knies+24 Monogem 2-SNR / 1-SNR, G205.6+12.4, PSR B0656+14
   now and at birth, Geminga now, Silk curvature centre); zone statistics; crossing of the X-ray radii; sky-plane
   bend angle beta = atan(d theta / d arc) of the spine relative to the azimuth about each centre.
2. Small-circle fit to the northern spine alone (Dec >= 17.6 and >= 19.5), point bootstrap, compared with the
   candidate centres.
3. Positions of the Halpha patch and PGCC G199.09+19.95 relative to the spine and the centres.
4. LSR corrections along the spine (standard solar motion used by the paper / wave 2, and Schoenrich+2010).
5. HI (GALFA, verify_hiwham on-spine spectra) per [O III] velocity bin vs [O III] V_LSR (DESI, SDSS).
6. Expanding-shell projection models for the north/south velocity difference:
   S  one sphere about a Monogem centre that encloses the whole spine (far/near side), with its limb gain;
   T  locally tangent sheets whose normal is tilted out of the sky plane by i (v = v_pre + u sin i);
   CH a tangent curved sheet whose chord emission is one-sided (near half missing).
Inputs (read only): silk_trace.json, verify_hiwham_hi_extra.{json,npz}, w3_northsouth_desi.json,
w3_northsouth_desi_velcheck.json, w3_northsouth_sdss.json.
Output: data/derived/w3_northsouth_geom.json
"""
import json

import numpy as np
from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u
from scipy.optimize import least_squares

ROOT = '/Users/mickey/astronomy/silk-nebula'
DER = f'{ROOT}/data/derived'
RNG = np.random.default_rng(1957)
D_C = 288.0           # pc, PSR B0656+14 parallax distance (3.47 +- 0.36 mas)
D_C_RANGE = (262.0, 320.0)
SHEET_FWHM_DEG = 0.28


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)


def sep(ra1, de1, ra2, de2):
    a, b = unit(ra1, de1), unit(ra2, de2)
    return np.degrees(np.arccos(np.clip(np.sum(a * b, -1), -1, 1)))


def lsr_corrections(ra, dec):
    c = SkyCoord(ra * u.deg, dec * u.deg)
    apex = SkyCoord(270.9595 * u.deg, 30.0047 * u.deg)          # 18h03m50.29s +30d00m16.8s, 20 km/s
    std = 20.0 * np.cos(c.separation(apex).rad)
    g = c.galactic
    l, b = g.l.rad, g.b.rad
    U, V, W = 11.1, 12.24, 7.25                                 # Schoenrich, Binney & Dehnen 2010
    sch = U * np.cos(l) * np.cos(b) + V * np.sin(l) * np.cos(b) + W * np.sin(b)
    return std, sch                                            # V_LSR = V_helio + corr


def circle_fit(ra, dec, sig, x0):
    def res(p):
        return (sep(ra, dec, p[0], p[1]) - p[2]) / sig
    r = least_squares(res, x0, x_scale=[1, 1, 1])
    return r.x, float(np.sum(r.fun ** 2))


def main():
    tr = json.load(open(f'{DER}/silk_trace.json'))
    sp = tr['spine']
    ra = np.array([p['ra'] for p in sp])
    de = np.array([p['dec'] for p in sp])
    sg = np.array([p['sigma_deg'] for p in sp])
    sg = np.maximum(sg, 0.02)
    lb = SkyCoord(ra * u.deg, de * u.deg).galactic
    birth = SkyCoord(l=200.6 * u.deg, b=6.9 * u.deg, frame=Galactic).icrs
    centres = {
        'Monogem_2SNR_Knies24': (104.15, 15.30, 11.00),
        'Monogem_1SNR_Knies24': (107.08, 13.50, 11.97),
        'G205.6+12.4_Knies24': (110.33, 12.94, 8.82),
        'PSR_B0656+14_now': (104.9506, 14.2393, None),
        'PSR_B0656+14_birth_context': (float(birth.ra.deg), float(birth.dec.deg), None),
        'Geminga_now': (98.4756, 17.7703, None),
        'Silk_curvature_centre': (94.5, 16.76, 21.3),
    }
    zones = {'south_<17.6': de < 17.6, 'mid_17.6-21': (de >= 17.6) & (de < 21), 'north_>=21': de >= 21,
             'north_>=19.5': de >= 19.5}
    # arc length along the spine
    v = unit(ra, de)
    ds = np.r_[0, np.degrees(np.arccos(np.clip(np.sum(v[1:] * v[:-1], 1), -1, 1)))]
    arc = np.cumsum(ds)
    out = {'spine_dec': de.tolist(), 'spine_ra': ra.tolist(), 'spine_arc_deg': arc.tolist(), 'centres': {}}
    for nm, (cr, cd, R) in centres.items():
        th = sep(ra, de, cr, cd)
        dth = np.gradient(th, arc)
        k = 5
        dth_s = np.convolve(dth, np.ones(k) / k, mode='same')
        beta = np.degrees(np.arctan(dth_s))
        z = {}
        for zn, m in zones.items():
            z[zn] = dict(mean=float(th[m].mean()), rms=float(th[m].std()), min=float(th[m].min()),
                         max=float(th[m].max()), beta_mean_deg=float(beta[m].mean()))
        cross = {}
        msouth = de < 17.6
        thr = [R] if R else []
        if nm.startswith('Monogem_2'):
            thr += [12.0, float(th[msouth].mean() + 3 * th[msouth].std())]
        for Rr in thr:
            below = np.where(th <= Rr)[0]          # last point at or below the radius, then permanently above
            if below.size and below[-1] < th.size - 1:
                j = below[-1] + 1
                f = (Rr - th[j - 1]) / (th[j] - th[j - 1])
                cross[f'{Rr:.2f}'] = float(de[j - 1] + f * (de[j] - de[j - 1]))
        out['centres'][nm] = dict(ra=cr, dec=cd, R_deg=R, theta=th.round(3).tolist(), beta_deg=beta.round(1).tolist(),
                                  zones=z, dec_where_theta_exceeds=cross)
        print(f"{nm:28s}", {zn: (round(v_['mean'], 2), round(v_['rms'], 2)) for zn, v_ in z.items()}, cross)

    # ------------------------------------------------------------ northern small-circle fits
    fits = {}
    for nm, m in (('dec>=17.6', de >= 17.6), ('dec>=19.5', de >= 19.5), ('dec<17.6', de < 17.6)):
        best = None
        for x0 in ([104.0, 15.0, 13.0], [110.0, 13.0, 9.0], [95.0, 17.0, 21.0], [85.0, 15.0, 30.0]):
            p, c2 = circle_fit(ra[m], de[m], sg[m], x0)
            if best is None or c2 < best[1]:
                best = (p, c2)
        p0 = best[0]
        boots = []
        idx = np.where(m)[0]
        for _ in range(400):
            ii = RNG.choice(idx, idx.size)
            try:
                pb, _ = circle_fit(ra[ii], de[ii], sg[ii], p0)
                boots.append(pb)
            except Exception:  # noqa: BLE001
                pass
        boots = np.array(boots)
        seps = {c: float(sep(p0[0], p0[1], cr, cd)) for c, (cr, cd, _) in centres.items()}
        bsep = {c: np.percentile(sep(boots[:, 0], boots[:, 1], cr, cd), [2.5, 50, 97.5]).tolist()
                for c, (cr, cd, _) in centres.items()}
        # rms about each candidate centre (concentricity) for these points
        conc = {c: dict(mean=float(sep(ra[m], de[m], cr, cd).mean()), rms=float(sep(ra[m], de[m], cr, cd).std()))
                for c, (cr, cd, _) in centres.items()}
        extra = dict(frac_boot_R_lt_15=float(np.mean(boots[:, 2] < 15.0)),
                     frac_boot_R_lt_13=float(np.mean(boots[:, 2] < 13.0)))
        # fixed Knies 2-SNR centre, free radius: chi2 vs free circle
        th_k = sep(ra[m], de[m], 104.15, 15.30)
        Rk = float(np.sum(th_k / sg[m] ** 2) / np.sum(1 / sg[m] ** 2))
        extra['knies2_fixed_centre'] = dict(R=Rk, chi2=float(np.sum(((th_k - Rk) / sg[m]) ** 2)),
                                            rms=float(np.std(th_k)))
        th_c = sep(ra[m], de[m], 94.5, 16.76)
        Rc = float(np.sum(th_c / sg[m] ** 2) / np.sum(1 / sg[m] ** 2))
        extra['curvature_centre_fixed'] = dict(R=Rc, chi2=float(np.sum(((th_c - Rc) / sg[m]) ** 2)),
                                               rms=float(np.std(th_c)))
        fits[nm] = dict(**extra, centre_ra=float(p0[0]), centre_dec=float(p0[1]), R_deg=float(p0[2]), chi2=best[1],
                        n=int(m.sum()), boot_ra_p16_p84=np.percentile(boots[:, 0], [16, 84]).tolist(),
                        boot_dec_p16_p84=np.percentile(boots[:, 1], [16, 84]).tolist(),
                        boot_R_p16_p84=np.percentile(boots[:, 2], [16, 84]).tolist(),
                        sep_from_candidates=seps, sep_boot_2p5_50_97p5=bsep, concentricity=conc)
        print('circle', nm, np.round(p0, 2), 'R16-84', np.round(fits[nm]['boot_R_p16_p84'], 1),
              {c: round(s, 1) for c, s in seps.items()})
    out['circle_fits'] = fits

    # ------------------------------------------------------------ Halpha patch and PGCC
    hx = json.load(open(f'{DER}/verify_hiwham_hi_extra.json'))
    pg = hx['C_pgcc_patch']['pgcc_radec']
    pm = tr['halpha_patch_info']['moments_halfmax']
    objs = {'halpha_patch_centroid': (pm['centroid_ra'], pm['centroid_dec']), 'PGCC_G199.09+19.95': tuple(pg)}
    P = unit(ra, de)
    pos = {}
    for nm, (r0, d0) in objs.items():
        x = unit(r0, d0)
        j = int(np.argmax(P @ x))
        dmin = float(sep(r0, d0, ra[j], de[j]))
        east = np.sign(r0 - ra[j])
        pos[nm] = dict(ra=r0, dec=d0, nearest_spine_dec=float(de[j]), offset_from_spine_deg=float(east * dmin),
                       theta_from_centres={c: float(sep(r0, d0, cr, cd)) for c, (cr, cd, _) in centres.items()},
                       spine_theta_same_dec={c: float(out['centres'][c]['theta'][j]) for c in centres})
    pos['pgcc_to_patch_deg'] = float(sep(objs['halpha_patch_centroid'][0], objs['halpha_patch_centroid'][1],
                                         pg[0], pg[1]))
    out['patch_pgcc'] = pos
    print('patch/pgcc', json.dumps(pos)[:600])

    # ------------------------------------------------------------ LSR corrections
    std, sch = lsr_corrections(ra, de)
    out['lsr_corr_std'] = std.round(3).tolist()
    out['lsr_corr_schoenrich'] = sch.round(3).tolist()

    def lsr_at(dlo, dhi, which=std):
        m = (de >= dlo) & (de < dhi)
        return float(which[m].mean()) if m.any() else float(np.interp(0.5 * (dlo + dhi), de, which))

    # ------------------------------------------------------------ [O III] velocities and HI per bin
    dj = json.load(open(f'{DER}/w3_northsouth_desi.json'))
    vc = json.load(open(f'{DER}/w3_northsouth_desi_velcheck.json'))
    sj = json.load(open(f'{DER}/w3_northsouth_sdss.json'))
    hn = np.load(f'{DER}/verify_hiwham_hi_extra.npz')
    hv, hon, hdec = hn['A_v'], hn['A_on'], hn['A_dec']
    segs = hx['A_along_spine']['gauss_segments']

    def hi_bin(dlo, dhi):
        m = (hdec >= dlo) & (hdec < dhi)
        if not m.any():
            return None
        s = hon[m].mean(0)
        w = (hv > -15) & (hv < 30)
        j = np.argmax(np.where(w, s, -np.inf))
        m1 = float(np.sum(s[w] * hv[w]) / np.sum(s[w]))
        comps = [c for sgm in segs if dlo <= sgm['dec_mid'] < dhi for c in sgm['comps'] if c['NHI'] > 5e19]
        return dict(vpeak=float(hv[j]), m1=m1, comps=[(round(c['v'], 1), round(c['fwhm'], 1), float(c['NHI']))
                                                       for c in comps])

    vb = []
    edges = dj['v_edges']
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        rA = dj['vel_A_bins'][k]
        corr = lsr_at(lo, hi)
        vb.append(dict(bin=[lo, hi], lsr_corr=corr, lsr_corr_sch=lsr_at(lo, hi, sch),
                       desi_v_helio=rA['v_helio'], desi_err_formal_infl=rA['v_err_formal_infl'],
                       desi_err_boot=rA['boot']['sig'] if rA['boot'] else None, desi_reliable=rA['reliable'],
                       desi_snr=rA['amp_snr_formal_infl'], desi_v_lsr=rA['v_helio'] + corr, HI=hi_bin(lo, hi)))
    out['velocity_bins'] = vb
    zones_v = {}
    for zn, (lo, hi), key, kk in (('south_<19.5', (12, 19.5), 'halves', 'zone0'),
                                  ('north_>=19.5', (19.5, 28.6), 'halves', 'zone1'),
                                  ('south_<17.6', (12, 17.6), 'three', 'zone0'),
                                  ('mid_17.6-21', (17.6, 21), 'three', 'zone1'),
                                  ('north_>=21', (21, 28.6), 'three', 'zone2')):
        e = vc[key][kk]
        ests = {nm: (v_['v'], v_['boot']['sig']) for nm, v_ in e.items() if isinstance(v_, dict) and 'v' in v_}
        vals = np.array([x[0] for nm, x in ests.items() if nm != 'L4_4960_only'])
        corr = lsr_at(lo, hi)
        sd_key = {'south_<19.5': 'zone_dec<19.5', 'north_>=19.5': 'zone_dec>=19.5', 'south_<17.6': 'zone_dec<17.6',
                  'mid_17.6-21': 'zone_17.6-21', 'north_>=21': 'zone_>=21'}[zn]
        sv = sj['regions'][sd_key]['velocity_OIII']
        zones_v[zn] = dict(lsr_corr=corr, desi_estimators_helio=ests, desi_helio_range=[float(vals.min()),
                           float(vals.max())], desi_helio_mid=float(0.5 * (vals.min() + vals.max())),
                           desi_lsr_range=[float(vals.min() + corr), float(vals.max() + corr)],
                           sdss_helio=sv['v_corr'], sdss_err=sv['ev_boot'], sdss_lsr=sv['v_corr'] + corr,
                           HI=hi_bin(lo, hi))
        print('zone', zn, 'DESI helio', np.round(zones_v[zn]['desi_helio_range'], 1), 'LSR',
              np.round(zones_v[zn]['desi_lsr_range'], 1), 'SDSS', round(sv['v_corr'], 1), '+-', round(sv['ev_boot'], 1),
              'HI', zones_v[zn]['HI'])
    out['velocity_zones'] = zones_v
    dv = vc['halves']['diff_last_minus_first']
    out['dv_north_minus_south_helio'] = {k: [v_['p50'], v_['sig']] for k, v_ in dv.items()}
    corrS, corrN = lsr_at(12, 19.5), lsr_at(19.5, 28.6)
    out['lsr_corr_south_north'] = [corrS, corrN]
    out['dv_lsr_extra_from_solar_reflex'] = corrN - corrS

    # ------------------------------------------------------------ expanding-shell projection models
    models = {}
    dvh = np.array([v_[0] for k, v_ in out['dv_north_minus_south_helio'].items() if k != 'L4_4960_only'])
    dv_lsr = (dvh.min() + corrN - corrS, dvh.max() + corrN - corrS)
    models['observed_dv_lsr_range'] = list(dv_lsr)
    # Model S: one sphere enclosing the whole spine
    S = {}
    for c in ('Monogem_2SNR_Knies24', 'Monogem_1SNR_Knies24', 'PSR_B0656+14_now', 'PSR_B0656+14_birth_context'):
        th = np.array(out['centres'][c]['theta'])
        thR = th.max() + 0.05
        Rpc = D_C * np.sin(np.radians(thR))
        dl = D_C * np.radians(SHEET_FWHM_DEG)
        p = D_C * np.sin(np.radians(th))
        f = np.sqrt(np.clip(1 - (np.sin(np.radians(th)) / np.sin(np.radians(thR))) ** 2, 0, 1))
        Lside = np.sqrt(np.clip(Rpc ** 2 - p ** 2, 0, None)) - np.sqrt(np.clip((Rpc - dl) ** 2 - p ** 2, 0, None))
        G = Lside / dl
        Gt = np.sqrt(2 * Rpc * dl - dl ** 2) / dl
        mS, mN = de < 19.5, de >= 19.5
        fS, fN = f[mS].mean(), f[mN].mean()
        # fit v = v_sys + s * vexp * f with one common side; dv(N-S) = s*vexp*(fN - fS)
        vexp_needed = [float(x / (fN - fS)) for x in dv_lsr]
        S[c] = dict(theta_R=float(thR), R_pc=float(Rpc), f_south=float(fS), f_north=float(fN),
                    G_side_south_median=float(np.median(G[mS])), G_side_north_median=float(np.median(G[mN])),
                    G_side_tangent=float(Gt), vexp_needed_signed=vexp_needed,
                    note='negative v_exp = near (approaching) side; south would sit at -|v_exp| f_S from systemic, '
                         'double-peaked if both sides emit, and only G_south/G_tangent as limb-brightened')
        print('Model S', c, {k: (np.round(v_, 2) if isinstance(v_, (float, list)) else v_) for k, v_ in S[c].items()
                             if k != 'note'})
    models['S_single_sphere'] = S
    # Model T: tilted normal, v = v_pre + u sin(i), u = 0.85 v_s
    T = []
    for vs in (75.0, 100.0, 125.0, 150.0):
        uu = 0.85 * vs
        ii = np.degrees(np.arcsin(np.clip(np.array(dv_lsr) / uu, -1, 1)))
        th2 = np.array(out['centres']['Monogem_2SNR_Knies24']['theta'])
        rmin = D_C * np.sin(np.radians(th2[de >= 19.5].mean()))
        T.append(dict(v_s=vs, u=uu, tilt_deg=ii.tolist(), depth_offset_pc=(rmin * np.tan(np.radians(ii))).tolist(),
                      r_min_pc=float(rmin)))
    models['T_tilted_normal'] = T
    beta_n = out['centres']['Monogem_2SNR_Knies24']['zones']['north_>=19.5']['beta_mean_deg']
    models['T_sky_plane_bend_beta_north_deg'] = beta_n
    models['T_note'] = ('If the surface bent as much along the line of sight as it does in the sky plane '
                        f'(beta ~ {beta_n:.0f} deg), v_pre + u sin(beta) would be ~0.3-0.5 u = 20-60 km/s.')
    # Model CH: one-sided chord on a tangent curved sheet
    CH = []
    for Rc in (58.0, 79.0):
        for h in (0.75, 1.3, 1.7):
            for vs in (75.0, 100.0, 125.0, 150.0):
                uu = 0.85 * vs
                vmax = uu * np.sqrt(2 * h / Rc)
                CH.append(dict(R_c_pc=Rc, h_pc=h, v_s=vs, v_max=float(vmax), centroid_one_sided=float(2 / 3 * vmax)))
    models['CH_one_sided_chord'] = CH
    models['CH_note'] = ('uniform emissivity over depth 0..h behind the front, projected at offset h: '
                         'z = sqrt(2 R_c (h - s)), v_los = u z / R_c; full chord centroid 0, one-sided 2/3 v_max; '
                         'a far/near asymmetry a = (F_far - F_near)/(F_far + F_near) gives a * 2/3 * v_max')
    out['models'] = models
    json.dump(out, open(f'{DER}/w3_northsouth_geom.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
