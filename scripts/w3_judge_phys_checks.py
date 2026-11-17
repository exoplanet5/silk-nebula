"""w3_judge_phys: order-of-magnitude physical-consistency checks for the origin hypotheses.

Reads only project files (run with ~/.venvs/astro313/bin/python -I):
  data/derived/silk_trace.json                     spine (E strand)
  data/derived/verify_shockmodel_analysis.json     238 spectrally acceptable truncated-shock runs (closure W at 280 pc)
  data/derived/verify_shockmodel_summary.json      25 fully consistent runs
Writes data/derived/w3_judge_phys_checks.json.

Every number here is DERIVED (simple physics on quoted inputs); inputs and their sources are listed in the JSON.
"""
import json
import numpy as np

ROOT = '/Users/mickey/astronomy/silk-nebula/'
PC = 3.0857e18
YR = 3.156e7
KB = 1.380649e-16
MH = 1.6726e-24
OUT = {}

# ---------------------------------------------------------------- geometry
tr = json.load(open(ROOT + 'data/derived/silk_trace.json'))
ra = np.array([p['ra'] for p in tr['spine']])
dec = np.array([p['dec'] for p in tr['spine']])


def sep(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def pa(ra0, de0, ra1, de1):
    r0, d0, r1, d1 = map(np.radians, (ra0, de0, ra1, de1))
    y = np.sin(r1 - r0) * np.cos(d1)
    x = np.cos(d0) * np.sin(d1) - np.sin(d0) * np.cos(d1) * np.cos(r1 - r0)
    return np.degrees(np.arctan2(y, x)) % 360


K2 = (104.15, 15.30, 11.00)      # Knies+24 two-SNR Monogem centre and X-ray radius (context.md)
K1 = (107.08, 13.50, 11.97)      # Knies+24 single-SNR
south = dec < 17.6
th = sep(ra, dec, K2[0], K2[1])
pas = pa(K2[0], K2[1], ra, dec)
arc = np.sum(sep(ra[:-1], dec[:-1], ra[1:], dec[1:]))
geo = {
    'spine_arc_deg': float(arc),
    'south_mean_theta_K2_deg': float(th[south].mean()), 'south_rms_deg': float(th[south].std()),
    'PA_span_about_K2_deg': [float(pas.min()), float(pas.max())],
    'fraction_of_Monogem_circumference_covered': float((pas.max() - pas.min()) / 360.0),
    'length_pc': {str(d): float(np.radians(arc) * d) for d in (160, 252, 282, 318, 1000, 1600)},
    'width_FWHM_pc_0.28deg': {str(d): float(np.radians(0.28) * d) for d in (160, 282, 1000, 1600)},
}
OUT['geometry'] = geo

# ---------------------------------------------------------------- chance alignment (H1 vs unrelated object)
# Place the 23-point south segment at random (random spin about its centroid, centroid uniform on the
# sphere within theta_max of the K2 centre) and ask how often it lands at the Monogem radius AND tangential.
rng = np.random.default_rng(20261117)
xs = np.radians(ra[south]); ys = np.radians(dec[south])
# unit vectors of the south segment
v = np.stack([np.cos(ys) * np.cos(xs), np.cos(ys) * np.sin(xs), np.sin(ys)], 1)
c0 = v.mean(0); c0 /= np.linalg.norm(c0)
kc = np.array([np.cos(np.radians(K2[1])) * np.cos(np.radians(K2[0])),
               np.cos(np.radians(K2[1])) * np.sin(np.radians(K2[0])), np.sin(np.radians(K2[1]))])


def rotmat(axis, ang):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def rot_a_to_b(a, b):
    ax = np.cross(a, b); s = np.linalg.norm(ax); c = np.dot(a, b)
    if s < 1e-12:
        return np.eye(3)
    return rotmat(ax, np.arctan2(s, c))


obs_rbar = float(th[south].mean()); obs_rms = float(th[south].std())
chance = {}
N = 60000
for thmax in (20.0, 25.0, 30.0):
    # uniform on spherical cap: cos(theta) uniform
    cth = rng.uniform(np.cos(np.radians(thmax)), 1.0, N)
    phi = rng.uniform(0, 2 * np.pi, N)
    spin = rng.uniform(0, 2 * np.pi, N)
    # basis around kc
    e1 = np.cross(kc, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(kc, e1)
    sth = np.sqrt(1 - cth ** 2)
    tgt = (cth[:, None] * kc + sth[:, None] * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2))
    rb = np.empty(N); rr = np.empty(N)
    R0 = rot_a_to_b(c0, np.array([0, 0, 1.0]))
    vz = v @ R0.T                       # segment centred on the pole
    for i in range(N):
        Rs = rotmat(np.array([0, 0, 1.0]), spin[i])
        Rt = rot_a_to_b(np.array([0, 0, 1.0]), tgt[i])
        w = vz @ Rs.T @ Rt.T
        d = np.degrees(np.arccos(np.clip(w @ kc, -1, 1)))
        rb[i] = d.mean(); rr[i] = d.std()
    win = np.abs(rb - obs_rbar) <= 0.5
    tang = rr <= 0.25
    p_joint = float(np.mean(win & tang))       # per 1-deg-wide radius window
    # H1 likelihood for rbar: shock at or up to ~1.5 deg outside the X-ray radius 11.0 +- 1.1
    sig = np.hypot(1.1, 0.3)
    dens = np.mean([np.exp(-0.5 * ((obs_rbar - (K2[2] + off)) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
                    for off in np.linspace(0, 1.5, 31)])
    p_h1 = float(dens * 1.0)                   # same 1-deg window; P(tangential | limb) ~ 1
    chance[str(int(thmax))] = {
        'P_random_rbar_within_0.5deg_of_obs': float(np.mean(win)),
        'P_random_tangential_rms_le_0.25': float(np.mean(tang)),
        'P_random_joint': p_joint,
        'P_limb_same_window': p_h1,
        'Bayes_factor_limb_vs_random': float(p_h1 / max(p_joint, 1.0 / N)),
    }
OUT['chance_alignment_south_segment'] = {
    'obs_rbar_deg': obs_rbar, 'obs_rms_deg': obs_rms, 'N_trials': N, 'by_theta_max_deg': chance,
    'note': 'Random rigid placement of the Dec<17.6 segment within theta_max of the Knies two-SNR centre. '
            'Limb model: rbar ~ N(11.0+off, 1.14) with off uniform 0-1.5 deg. Ignores the Fesen population of '
            'other Monogem-rim [O III] filaments (which would raise the limb prior) and any selection effects.'}

# ---------------------------------------------------------------- Monogem dynamics (Sedov) and cloud shocks
mon = {}
for Rdeg_d in ((11.0, 282.0), (11.0, 300.0), (11.97, 300.0)):
    R = np.radians(Rdeg_d[0]) * Rdeg_d[1] * PC
    for t in (1.1e5, 1.2e5):
        vbw = 0.4 * R / (t * YR)                      # Sedov
        for namb in (3.1e-3, 5.5e-3, 8.6e-3):
            rho = 1.4 * MH * namb
            E = (R / 1.15) ** 5 * rho / (t * YR) ** 2   # Sedov xi=1.15 (gamma 5/3)
            Ps = 0.75 * rho * vbw ** 2 / KB              # immediate post-shock P/k (strong shock)
            tpds = 3.61e4 * (E / 1e51) ** (3 / 14) * namb ** (-4 / 7)   # Cioffi+88 shell-formation time
            key = f'R{Rdeg_d[0]}deg_d{int(Rdeg_d[1])}_t{t:.1e}_n{namb}'
            mon[key] = {'R_pc': R / PC, 'v_bw_kms': vbw / 1e5, 'E_erg': E, 'P_postshock_k': Ps,
                        't_PDS_yr': tpds, 'radiative_now': bool(tpds < t)}
vb = np.array([m['v_bw_kms'] for m in mon.values()])
cloud = {}
for chi in (2, 4, 8):
    for beta in (1.0, 2.5):     # beta: reflected-shock enhancement of the driving pressure
        cloud[f'chi{chi}_beta{beta}'] = [float(vb.min() * np.sqrt(beta / chi)), float(vb.max() * np.sqrt(beta / chi))]
OUT['monogem_sedov'] = {'cases': mon, 'v_bw_range_kms': [float(vb.min()), float(vb.max())],
                        'cloud_shock_v_kms_vbw_sqrt(beta/chi)': cloud,
                        'inputs': 'Knies+24 (context.md): R 11.0 deg (2-SNR) / 11.97 (1-SNR), age 1.1-1.2e5 yr, '
                                  'n0 f^-1/2 3.1-8.6e-3; PSR distance 288 pc; Silk limb 282 pc'}

# ---------------------------------------------------------------- closure-W distance scan of the shock models
an = json.load(open(ROOT + 'data/derived/verify_shockmodel_analysis.json'))
acc = an['acceptable']
I_lo, I_hi = 1.2, 2.5          # south ridge range (DESI 1.62; SDSS 1.4-2.5)
scan = {}
for d in (160.0, 282.0, 500.0, 1000.0, 1600.0):
    f = 280.0 / d
    rows = []
    for a in acc:
        W = a['W']
        I = W['I_pk'] * f
        ok = (I_lo <= I <= I_hi) and (W['line_fwhm'] <= 38.0) and (40 <= a['FUV_RT']['total_CU_per_R_ebv_sfd'] <= 160)
        if ok:
            rows.append((a['v'], a['b'], W['n0'] * f, W['P_ram_k'] * f, W['t_star_yr'] / f, W['t_cool_yr'] / f))
    rows = np.array(rows) if rows else np.zeros((0, 6))
    scan[str(int(d))] = {
        'n_match': int(len(rows)), 'n_acceptable': len(acc),
        'b_values': sorted(set(rows[:, 1].tolist())) if len(rows) else [],
        'v_range': [float(rows[:, 0].min()), float(rows[:, 0].max())] if len(rows) else None,
        'n0_range': [float(rows[:, 2].min()), float(rows[:, 2].max())] if len(rows) else None,
        'Pram_k_range': [float(rows[:, 3].min()), float(rows[:, 3].max())] if len(rows) else None,
        'Pram_k_median': float(np.median(rows[:, 3])) if len(rows) else None,
        't_star_yr_range': [float(rows[:, 4].min()), float(rows[:, 4].max())] if len(rows) else None,
        't_cool_yr_range': [float(rows[:, 5].min()), float(rows[:, 5].max())] if len(rows) else None,
    }
OUT['closureW_distance_scan'] = {
    'scaling': 'At fixed b and observed angular width, n0 ∝ 1/d, I_pk ∝ 1/d, P_ram ∝ 1/d, t* ∝ d (flow is scale-free in n0; '
               'G fixed because R_c and width scale together). Criteria: I_pk in [1.2, 2.5] R, line FWHM <= 38 km/s, '
               'FUV(RT, E(B-V)=0.03) 40-160 CU/R. Grid b in {1,3,5,7,10,14,20,30}: matches at the grid edge (b=30) '
               'mean the true range extends beyond the grid.',
    'by_distance_pc': scan}

# ---------------------------------------------------------------- H3: paper model (160 pc, R 30 pc, E 2e51, n0 1, PDS)
R = 30 * PC; E = 2e51
h3 = {}
for fth in (0.3, 0.5):
    P = (2 / 3) * fth * E / (4 / 3 * np.pi * R ** 3)
    h3[f'interior_P_k_fth{fth}'] = P / KB
for n0 in (1.0, 0.03):
    w = np.radians(0.28) * 160 * PC              # observed FWHM at 160 pc
    nsh = n0 * R / (3 * w)                        # swept-up shell density for thickness w
    Llimb = 2 * np.sqrt(2 * R * w)
    h3[f'limb_NH_n0_{n0}'] = nsh * Llimb
    h3[f'limb_EM_if_ionised_pc_cm-6_n0_{n0}'] = nsh ** 2 * Llimb / PC
    h3[f'limb_Halpha_R_if_ionised_1e4K_n0_{n0}'] = nsh ** 2 * Llimb / PC / 2.75
h3['limits'] = {'HI_2sigma_cm-2': '1.0-1.75e19 per 2 km/s slice; ~3e19 for 6 km/s window (verify_hiwham)',
                'dust_NH_cm-2': '< 2.3e19 (E(B-V) < 0.004)', 'Halpha_2sigma_R': '< 0.77-1.0 (south)',
                'Monogem_P_k': '1.5e4-1.6e5 (definition dependent; canon C10)'}
OUT['H3_paper_model'] = h3

# ---------------------------------------------------------------- analogue homology: Cygnus Loop R88 SE vs Silk
# R88 SE: v 140 km/s, n0 4.4, face-on F(Hb) 0.092e-5 erg/cm2/s/sr, [O III]/Hb 12.9 (energy) (verify_shockmodel §3.2)
hnu_hb = 6.626e-27 * 2.998e10 / 4861.3e-8
ph_hb = 0.092e-5 / hnu_hb * 4 * np.pi
ph_o3 = ph_hb * 12.9 * 5006.8 / 4861.3
eta_r88 = ph_o3 / (4.4 * 1.4e7)
R_unit = 1e6 / (4 * np.pi)
cons = json.load(open(ROOT + 'data/derived/verify_shockmodel_summary.json'))['consistent_ranges']
OUT['analogue_homology'] = {
    'R88_SE_face_on_O3_R': ph_o3 / (4 * np.pi) / R_unit, 'R88_SE_eta5007_photons_per_H': eta_r88,
    'Silk_model_eta5007': cons['eta5007'], 'Silk_model_I_face_R': None,
    'n0_ratio_R88_over_Silk': [4.4 / cons['n0'][2], 4.4 / cons['n0'][0]],
    'note': 'Line ratios depend on (v_s, N*, b, X) only; face-on flux ∝ eta n0 v_s; layer thickness ∝ N*/n0. '
            'The Silk is a ~100-400x lower-density version of the R88 Cygnus Loop truncated-shock filaments.'}

# ---------------------------------------------------------------- north: X-ray surface brightness scaling with pressure
OUT['north_xray_scaling'] = {'P_N_over_P_S': [0.45, 0.7], 'Xray_SB_ratio_at_fixed_T': [0.45 ** 2, 0.7 ** 2],
                             'note': 'EM ∝ n^2 L ∝ P^2 at fixed T and path; northsouth P_N/P_S needed by the slower-shock fit.'}

# ---------------------------------------------------------------- penetration vs X-ray offsets
d = 282.0
pen = {'v_s_t_star_pc': [100e5 * 3.9e4 * YR / PC, 150e5 * 1.04e5 * YR / PC],
       'sstar_pc': cons['sstar_pc'], 'sstar_deg_at_282': [float(np.degrees(s / d)) for s in cons['sstar_pc']],
       'xray_step_location_deg_south': -0.2, 'ring_peak_offset_deg': 1.1,
       'ring_peak_offset_pc': float(np.radians(1.1) * d)}
OUT['penetration'] = pen

# ---------------------------------------------------------------- kinematics expected at H4 distances
R0, TH0 = 8.18, 236.0
lsr = {}
for dk in (0.16, 0.282, 1.0, 1.6, 2.0):
    vals = []
    for l, b in ((208.1, 17.5), (202.0, 19.3), (195.0, 21.0), (191.4, 21.85)):
        lr, br = np.radians(l), np.radians(b)
        dp = dk * np.cos(br)
        Rg = np.sqrt(R0 ** 2 + dp ** 2 - 2 * R0 * dp * np.cos(lr))
        vals.append(float(TH0 * (R0 / Rg - 1) * np.sin(lr) * np.cos(br)))
    lsr[str(dk)] = vals
OUT['flat_rotation_VLSR_along_spine'] = {'l_b_points': [[208.1, 17.5], [202.0, 19.3], [195.0, 21.0], [191.4, 21.85]],
                                         'V_LSR_kms_by_d_kpc': lsr,
                                         'observed_V_LSR': 'south -2.8..-0.2 (DESI), adopted +1.2+-1.0; north +4.5..+10.4'}

json.dump(OUT, open(ROOT + 'data/derived/w3_judge_phys_checks.json', 'w'), indent=1, default=float)
print(json.dumps({k: OUT[k] for k in ('geometry', 'chance_alignment_south_segment', 'closureW_distance_scan',
                                      'H3_paper_model', 'analogue_homology', 'penetration',
                                      'flat_rotation_VLSR_along_spine')}, indent=1, default=float)[:9000])
mm = OUT['monogem_sedov']
print('v_bw range', mm['v_bw_range_kms']); print(mm['cloud_shock_v_kms_vbw_sqrt(beta/chi)'])
for k, m in list(mm['cases'].items())[:12]:
    print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in m.items()})
