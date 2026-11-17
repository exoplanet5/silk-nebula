"""w4_distance: expected proper motion of the Silk filament and a two-epoch imaging plan.

Inputs (all already in the project; no downloads):
  shock speed v_s 100-150 km/s (best 110-140; verify_shockmodel, w3_canon section 10), centroid = pre-shock velocity
  ([O III] V_LSR +1.2 km/s, w3_canon), spine geometry (data/derived/silk_trace.json), Monogem two-SNR centre
  (104.15, +15.30) and curvature centre (94.50, +16.76) (w3_canon).
Computes
  1. mu = v_t / (4.74 d) for d = 160 / 282 / 1000 pc (+252, 318 for the H1 band) and v_s = 100, 120, 140, 150.
  2. Reflex term: proper motion of gas at rest in the LSR (Schoenrich+10 solar motion U,V,W = 11.1, 12.24, 7.25)
     along the spine, split into components normal and parallel to the filament, plus differential Galactic
     rotation (Oort A = 15.3, B = -11.9 km/s/kpc; negligible).
  3. Direction of the outward shock normal along the spine, and the radial direction from the Monogem and
     curvature centres (two driver hypotheses -> different motion PAs).
  4. Position precision per epoch for (a) a sharp edge and (b) the smooth Gaussian south profile, by Fisher
     information, for the SWONS 2024 epoch (Takahashi Epsilon-180ED, 500 mm f/2.8, QHY600M 3.76 um -> 1.551"/px;
     36 h [O III]) and for a 0.6-m f/4 class telescope at 0.6"/px. Signal and sky models are stated in the JSON.
Outputs: data/derived/w4_distance_pm.json, data/derived/w4_distance_pm_table.csv
Usage: ~/.venvs/astro313/bin/python -I scripts/w4_distance_pm.py
"""
import json
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS, Galactocentric, LSR, CartesianDifferential

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
K = 4.74047
DISTS = [160, 252, 282, 318, 1000]
VS = [100, 120, 140, 150]
BASE = [10, 12, 20]

rows = []
for d in DISTS:
    for v in VS:
        mu = v / (K * d / 1000.0)
        rows.append({'d_pc': d, 'v_s_kms': v, 'mu_mas_yr': round(mu, 1),
                     **{f'shift_{T}yr_arcsec': round(mu * T / 1000, 3) for T in BASE},
                     'shift_10yr_SWONS_px': round(mu * 10 / 1000 / 1.551, 3),
                     'shift_10yr_0.6arcsec_px': round(mu * 10 / 1000 / 0.6, 3)})

# ---- spine geometry and directions ----
tr = json.load(open(ROOT / 'data/derived/silk_trace.json'))
sp = tr['spine']
ra = np.array([q['ra'] for q in sp]); de = np.array([q['dec'] for q in sp])
mon = SkyCoord(104.15 * u.deg, 15.30 * u.deg)
cur = SkyCoord(94.50 * u.deg, 16.76 * u.deg)
geo = []
for i in range(0, len(ra), 4):
    j0, j1 = max(i - 1, 0), min(i + 1, len(ra) - 1)
    p = SkyCoord(ra[i] * u.deg, de[i] * u.deg)
    pa_t = p.position_angle(SkyCoord(ra[j1] * u.deg, de[j1] * u.deg)).deg  # along increasing Dec
    pa_n = (pa_t + 90) % 360 if ((pa_t + 90) % 360) < 180 else (pa_t - 90) % 360   # normal pointing east (0<PA<180)
    pa_mon = mon.position_angle(p).deg % 360        # radial outward from Monogem centre
    pa_cur = cur.position_angle(p).deg % 360
    # reflex: gas at rest in LSR, seen from the Sun, at d = 282 pc (scales as 1/d)
    out = {}
    for d in (160, 282, 1000):
        g = SkyCoord(ra=ra[i] * u.deg, dec=de[i] * u.deg, distance=d * u.pc, pm_ra_cosdec=0 * u.mas / u.yr,
                     pm_dec=0 * u.mas / u.yr, radial_velocity=0 * u.km / u.s, frame=LSR)
        h = g.transform_to(ICRS())
        mra, mde = h.pm_ra_cosdec.to_value(u.mas / u.yr), h.pm_dec.to_value(u.mas / u.yr)
        pan = np.radians(pa_n)
        nrm = mra * np.sin(pan) + mde * np.cos(pan)          # component along the east normal
        par = np.hypot(mra, mde) ** 2 - nrm ** 2
        out[d] = {'pmra': round(mra, 2), 'pmdec': round(mde, 2), 'normal_east': round(nrm, 2),
                  'parallel_abs': round(float(np.sqrt(max(par, 0))), 2), 'rv_helio': round(h.radial_velocity.to_value(u.km / u.s), 2)}
    geo.append({'ra': ra[i], 'dec': de[i], 'pa_tangent': round(pa_t, 1), 'pa_normal_east': round(pa_n, 1),
                'pa_radial_from_Monogem': round(pa_mon, 1), 'pa_radial_from_curv_centre': round(pa_cur, 1),
                'reflex_LSR_gas_mas_yr': out})

# ---- precision model ----
# Signal: I([O III]) = 2 R south peak (1.4-2.5 R canon); 0.6 R north. 1 R = 1e6/(4 pi) ph s-1 cm-2 sr-1.
# SWONS epoch: 180 mm aperture, ~30 % central obstruction (area 0.7), throughput*QE 0.55, 6 nm filter, 36 h.
# Sky continuum V = 21.3 mag/arcsec^2 (rural/dark; assumption) -> f_lambda = 3.63e-9*10^(-0.4V) erg s-1 cm-2 A-1 arcsec-2.
def snr_per_px(I_R, ap_cm, obst, eff, bw_A, t_s, pix_as, Vsky=21.3, rn=3.5):
    area = np.pi * (ap_cm / 2) ** 2 * (1 - obst)
    om = (pix_as / 206265.0) ** 2
    sig = I_R * 1e6 / (4 * np.pi) * om * area * eff * t_s
    flam = 3.63e-9 * 10 ** (-0.4 * Vsky)
    sky = flam * pix_as ** 2 * bw_A / (6.626e-27 * 2.998e18 / 5007.0) * area * eff * t_s
    nsub = t_s / 300.0
    noise = np.sqrt(sig + sky + nsub * rn ** 2)
    return sig / noise, sig, sky


def sigma_edge(snr_px, w_edge_as, L_deg, pix_as):
    """Fisher: linear ramp of width w (px) over a flat zero side; per-row sigma = sqrt(w)/snr; average over L rows."""
    w = w_edge_as / pix_as
    rows_ = L_deg * 3600 / pix_as
    return np.sqrt(w) / snr_px / np.sqrt(rows_) * pix_as      # arcsec


def sigma_gauss(snr_px, fwhm_as, L_deg, pix_as):
    s = fwhm_as / 2.3548 / pix_as
    fisher_row = snr_px ** 2 * np.sqrt(np.pi) / (2 * s)          # sum (f')^2 for unit-amplitude Gaussian
    rows_ = L_deg * 3600 / pix_as
    return 1 / np.sqrt(fisher_row * rows_) * pix_as


setups = {
    'SWONS_eps180ED_2024 (epoch 1)': dict(ap_cm=18.0, obst=0.30, eff=0.55, bw_A=60.0, t_s=36 * 3600, pix_as=1.551),
    'SWONS-like repeat, 3 nm filter, 60 h': dict(ap_cm=18.0, obst=0.30, eff=0.55, bw_A=30.0, t_s=60 * 3600, pix_as=1.551),
    '0.6 m f/4, 3 nm, 30 h, 0.6"/px': dict(ap_cm=60.0, obst=0.35, eff=0.5, bw_A=30.0, t_s=30 * 3600, pix_as=0.6),
    'NSNS public product 20"/px (for reference)': dict(ap_cm=13.5 / 2.0 * 2, obst=0.0, eff=0.5, bw_A=60.0, t_s=20 * 3600, pix_as=20.0),
}
prec = {}
for name, s in setups.items():
    out = {}
    for lab, I in (('south_2R', 2.0), ('north_0.6R', 0.6)):
        snr, sig, sky = snr_per_px(I, s['ap_cm'], s['obst'], s['eff'], s['bw_A'], s['t_s'], s['pix_as'])
        out[lab] = {'snr_per_px': round(float(snr), 2), 'signal_e': round(float(sig), 1), 'sky_e': round(float(sky), 1),
                    'sigma_edge_arcsec': {f'w{w}"_L{L}deg': round(float(sigma_edge(snr, w, L, s['pix_as'])), 3)
                                          for w in (10, 60, 300) for L in (0.5, 2.0)},
                    'sigma_gauss_FWHM0.28deg_arcsec': {f'L{L}deg': round(float(sigma_gauss(snr, 0.28 * 3600, L, s['pix_as'])), 3)
                                                       for L in (1.0, 5.0)}}
    prec[name] = {'setup': s, **out}

res = {'date': '2026-11-17', 'K': K, 'table': rows, 'geometry_every_4th_spine_point': geo, 'precision_model': prec,
       'assumptions': {
           'v_s': '100-150 km/s, best 110-140 (verify_shockmodel); the [O III] zone sits ~1e15-1e16 cm behind the front, so it moves with the front (pattern speed = v_s relative to the pre-shock gas)',
           'transverse_fraction': 'limb (edge-on) geometry: motion ~ fully transverse; a tilt i from the plane of the sky reduces mu by cos i',
           'pre_shock_gas': 'assumed at rest in the LSR (V_LSR +1.2 radial measured); its transverse motion is unknown (+-10 km/s typical ISM -> +-7 mas/yr at 282 pc)',
           'SWONS_optics': 'Takahashi Epsilon-180ED (180 mm, 500 mm f/2.8; nominal spec, not from the paper) + QHY600 (3.76 um) -> 1.551 arcsec/px, 4.1 x 2.8 deg field; 36 h [O III] (+24 h Ha), 2024-01-10/11 onward per AstroBin/discovery.md',
           'sky': 'V = 21.3 mag/arcsec^2 continuum in the filter band; read noise 3.5 e-, 300 s subs; no flat-field / star-residual systematics (these will dominate below ~0.1-0.2 arcsec)',
           'edge_width': 'unknown below the 20"/px NSNS sampling: 10, 60, 300 arcsec tried',
       }}
json.dump(res, open(ROOT / 'data/derived/w4_distance_pm.json', 'w'), indent=1)
import csv
with open(ROOT / 'data/derived/w4_distance_pm_table.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
for r in rows:
    if r['v_s_kms'] in (100, 120, 150):
        print(r)
for gq in geo[::3]:
    print(gq['dec'], gq['pa_normal_east'], gq['pa_radial_from_Monogem'], gq['pa_radial_from_curv_centre'], gq['reflex_LSR_gas_mas_yr'][282])
print(json.dumps(prec, indent=1))
