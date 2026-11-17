"""w4_geostats task 1: chance coincidence of the southern Silk spine with the Monogem radius, one null for all centres,
plus the eRASS1 0.2-0.6 keV shell radius measured in the Silk position-angle sector itself.

(A) eRASS1 shell radius. Pixels: verify_xrayfuv_erass_pixels.npz, soft = b0204 + b0406 (cts/s/deg2), point-source
    catalogue mask applied. For each centre (K1 Knies 1-SNR, K2 Knies 2-SNR Monogem, PSR B0656+14), radial profile
    I(r) (0.2-deg bins, r 3-20 deg) in (i) the PA sector spanned by the southern spine (Dec 12-17.6) and (ii) every
    15-deg PA sector. Smoothed (Gaussian sigma 0.4 deg) estimators: r_peak (max at 6-16 deg); r_half (outer half
    maximum between the peak and the median level 3-6 deg beyond it); r_grad (steepest decline within 4 deg
    outside the peak). Errors: bootstrap over 5-deg PA sub-sectors. PA scatter over the ring = empirical "radius
    uncertainty" for matching to a circle.

(B) One null. H0: the Silk is unrelated to Monogem; its position and orientation are random within the same Galactic
    latitude band. The full 66-point trace (rigid shape) is placed with its centroid uniform in area within
    |b - b_obs| <= 3 deg (variant 6 deg), |l - 201.5| <= 45 deg, random spin (variant: random mirror too).
    Event E ("Monogem-radius tangency"): for ANY centre c in {K1, K2, PSR} (variant: + G205.6+12.4), some contiguous
    23-point window (5.5 deg of arc, the length of the observed southern section) has rms(distance) <= r0 and
    |mean distance - R_c| <= z0 * sigma_c. Primary r0 = 0.15 deg, z0 = 1. Radii: Knies K2 11.00 +- 1.10;
    K1 11.97 +- 1.20 (10 %, assumed like K2); PSR = eRASS1 ring-median outer half-max radius measured here +- 10 %
    (variant psr_wide: +- its PA scatter).
    Variant radii: all three from the eRASS1 measurement (r_half +- PA scatter).
    Reported unconditionally (per band area) and conditional on the trace reaching the rim zone of any centre
    (some point within 1.5 deg of R_c; variant 3 deg), which is the selection that made Monogem the comparison.
(C) Bayes factor BF = P(E | H1) / P(E | H0, rim zone), with P(E | H1) = P(|dbar - R| <= z0 sigma) x P(rms <= r0 |
    limb) = 0.68 x (0.5-1.0). (D) Local-edge rotation null: rotate the southern section about each centre and ask how
    often it sits as close to the local eRASS1 outer half-max edge as observed. Also reproduces the judge_phys and judge_skeptic numbers with their own definitions.
Output: data/derived/w4_geostats_registration.json, w4_geostats_registration_erass.npz
Run: ~/.venvs/astro313/bin/python -I scripts/w4_geostats_registration.py
"""
import sys
import numpy as np
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from astropy.coordinates import SkyCoord
import astropy.units as u
from w4_geostats_common import (ra_sp, dec_sp, CENTRES, sep_pa, load_erass, DER, json_dump)

RNG = np.random.default_rng(4117)
OUT = DER / 'w4_geostats_registration.json'
CEN = dict(CENTRES)
CEN['G205'] = (110.33, 12.94)            # Knies+24 second remnant G205.6+12.4, R = 8.82
SOUTH = dec_sp <= 17.6
W = 23


# ------------------------------------------------------------------------------------------------ (A) eRASS1
def radial_estimators(rc, prof, dr):
    ok = np.isfinite(prof)
    if ok.sum() < 0.7 * prof.size:
        return dict(r_peak=np.nan, r_half=np.nan, r_grad=np.nan)
    p = np.interp(rc, rc[ok], prof[ok])
    ps = gaussian_filter1d(p, 0.4 / dr, mode='nearest')
    win = (rc >= 6) & (rc <= 16)
    k = np.nonzero(win)[0][np.argmax(ps[win])]
    r_pk = rc[k]
    outz = (rc >= r_pk + 3) & (rc <= r_pk + 6)
    lev = np.median(ps[outz]) if outz.any() else np.nan
    half = 0.5 * (ps[k] + lev)
    beyond = np.nonzero((rc > r_pk) & (ps <= half))[0]
    r_half = np.nan
    if beyond.size:
        j = beyond[0]
        r_half = float(np.interp(half, [ps[j], ps[j - 1]], [rc[j], rc[j - 1]]))
    g = np.gradient(ps, dr)
    gz = (rc > r_pk) & (rc <= r_pk + 4)
    r_grad = float(rc[gz][np.argmin(g[gz])]) if gz.any() else np.nan
    return dict(r_peak=float(r_pk), r_half=r_half, r_grad=r_grad, peak=float(ps[k]), level_out=float(lev),
                contrast=float(ps[k] / lev) if lev > 0 else np.nan)


def erass_radius():
    E = load_erass()
    dr = 0.2
    edges = np.arange(3.0, 20.0001, dr); rc = 0.5 * (edges[1:] + edges[:-1])
    out = {}; store = {}
    for key in ('K1', 'K2', 'PSR'):
        r, pa = sep_pa(*CEN[key], E['ra'], E['dec'])
        _, pa_s = sep_pa(*CEN[key], ra_sp[SOUTH], dec_sp[SOUTH])
        _, pa_all = sep_pa(*CEN[key], ra_sp, dec_sp)
        d_s, _ = sep_pa(*CEN[key], ra_sp[SOUTH], dec_sp[SOUTH])
        sec = (float(pa_s.min()), float(pa_s.max()))
        ri = np.digitize(r, edges) - 1
        okr = (ri >= 0) & (ri < rc.size)
        # 5-deg sub-sector profiles over the full ring
        pa5 = np.arange(0, 360, 5.0)
        P5 = np.full((pa5.size, rc.size), np.nan)
        for i, a in enumerate(pa5):
            m = okr & (pa >= a) & (pa < a + 5)
            sw = np.bincount(ri[m], weights=E['v'][m], minlength=rc.size)
            n = np.bincount(ri[m], minlength=rc.size)
            P5[i] = np.where(n >= 3, sw / np.maximum(n, 1), np.nan)

        def sector_prof(a0, a1, idx=None):
            sel = np.nonzero((pa5 + 2.5 >= a0) & (pa5 + 2.5 <= a1))[0] if idx is None else idx
            with np.errstate(all='ignore'):
                return np.nanmean(P5[sel], axis=0), sel
        prof_s, sel_s = sector_prof(*sec)
        est_s = radial_estimators(rc, prof_s, dr)
        boots = []
        for _ in range(300):
            bi = RNG.choice(sel_s, sel_s.size, replace=True)
            boots.append(radial_estimators(rc, np.nanmean(P5[bi], axis=0), dr))
        be = {k: [b[k] for b in boots] for k in ('r_peak', 'r_half', 'r_grad')}
        est_s_err = {k: dict(sd=float(np.nanstd(v)), p16=float(np.nanpercentile(v, 16)), p84=float(np.nanpercentile(v, 84)))
                     for k, v in be.items()}
        # 15-deg sectors over the ring
        ring = []
        for a in np.arange(0, 360, 15.0):
            p, _ = sector_prof(a, a + 15)
            e = radial_estimators(rc, p, dr)
            e['pa0'] = float(a)
            ring.append(e)
        rv = {k: np.array([q[k] for q in ring]) for k in ('r_peak', 'r_half', 'r_grad')}
        # exclude the Silk sector itself from the ring statistics
        insec = np.array([(q['pa0'] + 7.5 >= sec[0]) & (q['pa0'] + 7.5 <= sec[1]) for q in ring])
        ringstat = {}
        for k, v in rv.items():
            vv = v[~insec & np.isfinite(v)]
            ringstat[k] = dict(n=int(vv.size), median=float(np.median(vv)), sd=float(np.std(vv)),
                               p16=float(np.percentile(vv, 16)), p84=float(np.percentile(vv, 84)))
        out[key] = dict(silk_sector_pa=sec, silk_south_mean_dist=float(d_s.mean()), silk_south_rms=float(d_s.std()),
                        silk_full_pa=(float(pa_all.min()), float(pa_all.max())),
                        sector=est_s, sector_boot=est_s_err,
                        silk_minus_rhalf=float(d_s.mean() - est_s['r_half']),
                        silk_minus_rpeak=float(d_s.mean() - est_s['r_peak']),
                        ring_excl_silk=ringstat, ring_sectors=ring)
        store[key] = dict(rc=rc, prof_s=prof_s, P5=P5, pa5=pa5)
        print(key, 'sector', sec, est_s, 'silk', d_s.mean(), flush=True)
    np.savez_compressed(DER / 'w4_geostats_registration_erass.npz',
                        **{f'{k}_{q}': v for k, d in store.items() for q, v in d.items()})
    return out


# ------------------------------------------------------------------------------------------------ (B) null
def gal_vec(ra, dec):
    c = SkyCoord(np.atleast_1d(ra) * u.deg, np.atleast_1d(dec) * u.deg).galactic
    l, b = np.radians(c.l.deg), np.radians(c.b.deg)
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], -1)


def frame(lb_vec):
    """Orthonormal frame (e_r, e_l, e_b) at unit vectors (N,3)."""
    er = lb_vec / np.linalg.norm(lb_vec, axis=-1, keepdims=True)
    z = np.array([0, 0, 1.0])
    el = np.cross(z, er); el /= np.linalg.norm(el, axis=-1, keepdims=True)
    eb = np.cross(er, el)
    return er, el, eb


TV = gal_vec(ra_sp, dec_sp)                      # (66, 3)
cen0 = TV.mean(0); cen0 /= np.linalg.norm(cen0)
er0, el0, eb0 = frame(cen0[None])
LOC = np.stack([TV @ er0[0], TV @ el0[0], TV @ eb0[0]], -1)     # local components (66, 3)
B_OBS = float(np.degrees(np.arcsin(cen0[2])))
L_OBS = float(np.degrees(np.arctan2(cen0[1], cen0[0])) % 360)
CV = {k: gal_vec(*v)[0] for k, v in CEN.items()}


def place(n, band=3.0, dl=45.0, l0=201.5, mirror=False):
    sb = RNG.uniform(np.sin(np.radians(B_OBS - band)), np.sin(np.radians(B_OBS + band)), n)
    b = np.arcsin(sb); l = np.radians(l0 + RNG.uniform(-dl, dl, n))
    tv = np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], -1)
    er, el, eb = frame(tv)
    phi = RNG.uniform(0, 2 * np.pi, n)
    a_l = LOC[None, :, 1]; a_b = LOC[None, :, 2]
    if mirror:
        flip = RNG.random(n) < 0.5
        a_l = np.where(flip[:, None], -a_l, a_l)
    c, s = np.cos(phi)[:, None], np.sin(phi)[:, None]
    al2 = a_l * c - a_b * s; ab2 = a_l * s + a_b * c
    P = LOC[None, :, 0, None] * er[:, None, :] + al2[..., None] * el[:, None, :] + ab2[..., None] * eb[:, None, :]
    return P, np.degrees(b), np.degrees(l)


def place_iso(n, centre, theta_max):
    """judge_phys-style: centroid uniform in area within theta_max of the centre, random spin."""
    cv = CV[centre]
    er_c, el_c, eb_c = frame(cv[None])
    rho = np.arccos(1 - RNG.random(n) * (1 - np.cos(np.radians(theta_max))))
    pa = RNG.uniform(0, 2 * np.pi, n)
    tv = (np.cos(rho)[:, None] * er_c + np.sin(rho)[:, None] * (np.cos(pa)[:, None] * eb_c + np.sin(pa)[:, None] * el_c))
    er, el, eb = frame(tv)
    phi = RNG.uniform(0, 2 * np.pi, n)
    c, s = np.cos(phi)[:, None], np.sin(phi)[:, None]
    al2 = LOC[None, :, 1] * c - LOC[None, :, 2] * s; ab2 = LOC[None, :, 1] * s + LOC[None, :, 2] * c
    P = LOC[None, :, 0, None] * er[:, None, :] + al2[..., None] * el[:, None, :] + ab2[..., None] * eb[:, None, :]
    return P


def win_stats(d):
    c = np.cumsum(np.concatenate([np.zeros((d.shape[0], 1)), d], 1), 1)
    c2 = np.cumsum(np.concatenate([np.zeros((d.shape[0], 1)), d * d], 1), 1)
    m = (c[:, W:] - c[:, :-W]) / W
    v = (c2[:, W:] - c2[:, :-W]) / W - m * m
    return m, np.sqrt(np.maximum(v, 0))


def dist(P, key):
    return np.degrees(np.arccos(np.clip(P @ CV[key], -1, 1)))


def evaluate(P, radii, r0s=(0.12, 0.15, 0.25), z0s=(0.5, 1.0, 1.5), zones=(1.5, 3.0)):
    """Returns per-placement boolean arrays: event[(key, r0, z0)], touch[(key, zone)]."""
    ev = {}; touch = {}; extra = {}
    for key, (R, sR) in radii.items():
        d = dist(P, key)
        m, r = win_stats(d)
        for r0 in r0s:
            for z0 in z0s:
                ev[(key, r0, z0)] = np.any((r <= r0) & (np.abs(m - R) <= z0 * sR), axis=1)
        for zn in zones:
            touch[(key, zn)] = np.any(np.abs(d - R) <= zn, axis=1)
        extra[key] = (m, r)
    return ev, touch, extra


def null_run(radii, n_total=600000, chunk=30000, band=3.0, mirror=False, keys_any=('K1', 'K2', 'PSR')):
    acc = None
    nb = 0
    for _ in range(n_total // chunk):
        P, b, l = place(chunk, band=band, mirror=mirror)
        ev, touch, _ = evaluate(P, radii)
        res = {}
        for (key, r0, z0), e in ev.items():
            res[f'ev|{key}|{r0}|{z0}'] = e
        for (key, zn), t in touch.items():
            res[f'touch|{key}|{zn}'] = t
        for r0 in (0.12, 0.15, 0.25):
            for z0 in (0.5, 1.0, 1.5):
                res[f'ev|any|{r0}|{z0}'] = np.any([ev[(k, r0, z0)] for k in keys_any], axis=0)
        for zn in (1.5, 3.0):
            res[f'touch|any|{zn}'] = np.any([touch[(k, zn)] for k in keys_any], axis=0)
        if acc is None:
            acc = {k: [] for k in res}
        for k, v in res.items():
            acc[k].append(v)
        nb += chunk
    acc = {k: np.concatenate(v) for k, v in acc.items()}
    out = dict(n=nb, band_deg=band, mirror=mirror, keys_any=list(keys_any))
    for r0 in (0.12, 0.15, 0.25):
        for z0 in (0.5, 1.0, 1.5):
            e = acc[f'ev|any|{r0}|{z0}']
            row = dict(p_uncond=float(e.mean()))
            for zn in (1.5, 3.0):
                t = acc[f'touch|any|{zn}']
                row[f'p_touch{zn}'] = float(t.mean())
                row[f'p_given_touch{zn}'] = float(e[t].mean())
                row[f'n_given_touch{zn}'] = int(t.sum())
            for k in keys_any:
                ek = acc[f'ev|{k}|{r0}|{z0}']; tk = acc[f'touch|{k}|1.5']
                row[f'{k}_p_given_touch1.5'] = float(ek[tk].mean()) if tk.any() else np.nan
            out[f'r0={r0}|z0={z0}'] = row
    # band area (sr) and per-deg^2 rate for rescaling to other windows
    area = np.radians(90.0) * (np.sin(np.radians(B_OBS + band)) - np.sin(np.radians(B_OBS - band))) * 2
    out['band_area_deg2'] = float(area * (180 / np.pi) ** 2)
    return out


def observed(radii, r0s=(0.12, 0.15, 0.25), z0s=(0.5, 1.0, 1.5)):
    P = TV[None]
    ev, touch, extra = evaluate(P, radii, r0s, z0s)
    out = {}
    for key, (R, sR) in radii.items():
        m, r = extra[key]
        d = dist(P, key)[0]
        k = int(np.argmin(r[0]))
        out[key] = dict(R=R, sigma_R=sR, south_mean=float(d[SOUTH].mean()), south_rms=float(d[SOUTH].std()),
                        best_window_rms=float(r[0, k]), best_window_mean=float(m[0, k]),
                        best_window_dec=(float(dec_sp[k]), float(dec_sp[k + W - 1])),
                        z_best=float((m[0, k] - R) / sR),
                        events={f'r0={r0}|z0={z0}': bool(ev[(key, r0, z0)][0]) for r0 in r0s for z0 in z0s},
                        touch15=bool(touch[(key, 1.5)][0]))
    return out


def reproduce_judges():
    out = {}
    # judge_phys: 23-pt south segment, random spin, centroid within theta_max of K2; event |mean - 12.10| <= 0.5 & rms <= 0.25
    seg = LOC[SOUTH]
    for th in (20.0, 30.0):
        n = 200000
        P = place_iso(n, 'K2', th)
        d = dist(P, 'K2')
        # use only the southern 23 points (judge_phys placed the segment, not the whole trace)
        ds = d[:, SOUTH]
        e = (np.abs(ds.mean(1) - 12.10) <= 0.5) & (ds.std(1) <= 0.25)
        out[f'phys_theta{int(th)}'] = dict(p=float(e.mean()), n=n, bf_with_limb_031=float(0.31 / max(e.mean(), 1e-9)))
    # judge_skeptic: whole trace, K2 only, rms <= 0.15, |mean - 11.0| <= 1.1, conditional on touching +-1.5
    for key, R in (('K2', 11.0), ('K1', 11.97)):
        n = 200000
        P = place_iso(n, key, 26.0)
        d = dist(P, key)
        m, r = win_stats(d)
        e = np.any((r <= 0.15) & (np.abs(m - R) <= 1.1), 1)
        t = np.any(np.abs(d - R) <= 1.5, 1)
        out[f'skeptic_{key}'] = dict(p_all=float(e.mean()), p_given_touch=float(e[t].mean()))
    del seg
    return out


def local_edge_null(A):
    """Rotate the southern section about each centre (distances to that centre are invariant) and ask how often
    the section sits as close to the LOCAL outer half-max edge of the eRASS1 shell as observed."""
    Z = np.load(DER / 'w4_geostats_registration_erass.npz')
    out = {}
    for key in ('K1', 'K2', 'PSR'):
        rc = Z[f'{key}_rc']; P5 = Z[f'{key}_P5']; pa5 = Z[f'{key}_pa5']
        dr = rc[1] - rc[0]
        a0, a1 = A[key]['silk_sector_pa']
        dbar = A[key]['silk_south_mean_dist']
        d_obs = dbar - A[key]['sector']['r_half']
        dn = []; thetas = []
        for th in np.arange(-180, 180, 5.0):
            if abs(th) < 30:
                continue
            b0 = (a0 + th) % 360; b1 = (a1 + th) % 360
            cen = (pa5 + 2.5)
            sel = np.nonzero(((cen >= b0) & (cen <= b1)) if b0 < b1 else ((cen >= b0) | (cen <= b1)))[0]
            with np.errstate(all='ignore'):
                prof = np.nanmean(P5[sel], axis=0)
            e = radial_estimators(rc, prof, dr)
            dn.append(dbar - e['r_half']); thetas.append(th)
        dn = np.array(dn)
        fin = np.isfinite(dn)
        out[key] = dict(delta_obs=float(d_obs), n=int(fin.sum()),
                        p_abs_le_obs=float((np.sum(np.abs(dn[fin]) <= abs(d_obs)) + 1) / (fin.sum() + 1)),
                        p_abs_le_0p3=float(np.mean(np.abs(dn[fin]) <= 0.3)),
                        null_sd=float(np.std(dn[fin])), null_median=float(np.median(dn[fin])),
                        theta=thetas, delta_null=dn.tolist())
        print('local edge', key, out[key]['delta_obs'], out[key]['p_abs_le_obs'], out[key]['null_sd'], flush=True)
    return out


def main():
    A = erass_radius()
    LE = local_edge_null(A)
    rx = {k: (A[k]['ring_excl_silk']['r_half']['median'], A[k]['ring_excl_silk']['r_half']['sd']) for k in A}
    radii_primary = {'K1': (11.97, 1.20), 'K2': (11.00, 1.10), 'PSR': (rx['PSR'][0], 0.1 * rx['PSR'][0]), 'G205': (8.82, 0.88)}
    radii_psrwide = dict(radii_primary, PSR=rx['PSR'])
    radii_erass = {'K1': rx['K1'], 'K2': rx['K2'], 'PSR': rx['PSR'], 'G205': (8.82, 0.88)}
    obs_p = observed(radii_primary); obs_e = observed(radii_erass)
    print('observed primary', {k: (v['south_mean'], v['events']['r0=0.15|z0=1.0']) for k, v in obs_p.items()}, flush=True)
    nulls = {}
    nulls['primary'] = null_run(radii_primary)
    print('primary', nulls['primary']['r0=0.15|z0=1.0'], flush=True)
    nulls['band6'] = null_run(radii_primary, band=6.0)
    nulls['mirror'] = null_run(radii_primary, mirror=True)
    nulls['erass_radii'] = null_run(radii_erass)
    nulls['psr_wide'] = null_run(radii_psrwide)
    nulls['with_G205'] = null_run(radii_primary, keys_any=('K1', 'K2', 'PSR', 'G205'))
    nulls['K2_only'] = null_run(radii_primary, keys_any=('K2',))
    for k, v in nulls.items():
        print(k, v['r0=0.15|z0=1.0'], flush=True)
    rep = reproduce_judges()
    print('reproduce', rep, flush=True)
    # Bayes factors (primary event r0 = 0.15, z0 = 1, rim zone 1.5 deg)
    p1 = (0.683 * 0.5, 0.683 * 1.0)
    bf = {}
    for k, v in nulls.items():
        row = v['r0=0.15|z0=1.0']
        p0 = row['p_given_touch1.5']
        bf[k] = dict(p0_given_touch=p0, bf_range=(p1[0] / p0, p1[1] / p0),
                     p0_given_touch3=row['p_given_touch3.0'], bf_touch3=(p1[0] / row['p_given_touch3.0'], p1[1] / row['p_given_touch3.0']))
    res = dict(erass_radius=A, local_edge_null=LE, radii_primary=radii_primary, radii_psrwide=radii_psrwide, radii_erass=radii_erass,
               observed_primary=obs_p, observed_erass=obs_e, nulls=nulls, reproduce_judges=rep,
               bayes=dict(p_E_given_H1=p1, note='P(E|H1) = P(|dbar-R|<=1 sigma)=0.683 x P(rms<=0.15 | limb) in 0.5-1.0', bf=bf),
               trace_centroid_lb=(L_OBS, B_OBS))
    json_dump(res, OUT)
    print('saved', OUT)


if __name__ == '__main__':
    main()
