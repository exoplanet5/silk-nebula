"""verify_distance (c): combine all distance evidence for the Silk Nebula into one table + figure.

Independent recomputation (own code) of:
  * Monogem-limb distance: Monte Carlo over the PSR B0656+14 parallax (3.47 +- 0.36 mas, Brisken+2003 via
    Thorsett+2003 / ATNF 'btgg03'), a centre-vs-pulsar depth offset, and the spine-centre angle measured here from
    data/derived/silk_trace.json against the Knies+2024 centres (values as given in notes/context.md).
  * Local Bubble wall from O'Neill+2024 (file downloaded by wave-1: data/raw/context_dust/...MeanDraws.fits) at the
    spine points (nearest pixel in l,b).
  * Kinematic V_LSR(d) for flat rotation (R0 8.18 kpc, Theta0 236 km/s) along the spine; LSR conversions.
  * Pulsar facts from the ATNF v2.9.0 text database on disk (data/raw/context_psr/psrcat_tar/psrcat.db).
  * HD 61717 (verify_distance_hd61717.json) and IUE (verify_distance_iue.json) results from this tag.
Outputs: data/derived/verify_distance_synthesis.json, figures/verify_distance_constraints.png
Run: ~/.venvs/astro313/bin/python -I scripts/verify_distance_synthesis.py
"""
import json
import sys
import numpy as np
from astropy.coordinates import SkyCoord, Galactic, LSR, ICRS, CartesianDifferential
from astropy.io import fits
import astropy.units as u
from scipy.spatial import cKDTree

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from silk_geometry import load_trace, ROOT  # noqa: E402

OUT = {}
rng = np.random.default_rng(288)
tr = load_trace()
sp_ra = np.array([p['ra'] for p in tr['spine']])
sp_dec = np.array([p['dec'] for p in tr['spine']])
spine = SkyCoord(sp_ra * u.deg, sp_dec * u.deg)

# ------------------------------------------------------------ 1. Monogem geometry
CEN2 = SkyCoord(104.15 * u.deg, 15.30 * u.deg)    # Knies+2024 two-SNR Monogem centre, R = 11.00 deg
CEN1 = SkyCoord(107.08 * u.deg, 13.50 * u.deg)    # Knies+2024 single-SNR centre, R = 11.97 deg
PSR = SkyCoord(104.9506 * u.deg, 14.2393 * u.deg)
th2 = spine.separation(CEN2).deg
th1 = spine.separation(CEN1).deg
thp = spine.separation(PSR).deg
south = sp_dec < 17.6
OUT['monogem_geometry'] = dict(
    spine_sep_two_snr_centre_south=[round(float(th2[south].mean()), 2), round(float(th2[south].std()), 2), int(south.sum())],
    spine_sep_two_snr_centre_north_range=[round(float(th2[~south].min()), 2), round(float(th2[~south].max()), 2)],
    spine_sep_single_centre_range=[round(float(th1.min()), 2), round(float(th1.max()), 2)],
    spine_sep_psr_south=[round(float(thp[south].mean()), 2), round(float(thp[south].std()), 2)],
    psr_sep_from_two_snr_centre=round(float(PSR.separation(CEN2).deg), 2),
    psr_sep_from_single_centre=round(float(PSR.separation(CEN1).deg), 2))

# Monte Carlo: limb distance d = D_c * cos(theta), D_c = D_psr + offset (pulsar depth offset from SNR centre)
N = 200000
plx = rng.normal(3.47, 0.36, N)
plx = plx[plx > 1.5]
Dp = 1000 / plx
theta = np.deg2rad(rng.normal(th2[south].mean(), max(th2[south].std(), 0.12), len(Dp)))
# centre depth offset: pulsar moved ~44 mas/yr * 1.1e5 yr = 1.35 deg on the sky; projected offset from centre
# 1-3 deg -> 5-15 pc; allow a radial offset of the same order (sigma 15 pc)
Dc = Dp + rng.normal(0, 15.0, len(Dp))
d_limb = Dc * np.cos(theta)
R_shell = Dc * np.sin(theta)
q = np.percentile(d_limb, [2.5, 16, 50, 84, 97.5])
OUT['monogem_limb_mc'] = dict(d_limb_pc_percentiles_2p5_16_50_84_97p5=[round(float(x)) for x in q],
                              R_shell_pc_16_50_84=[round(float(x)) for x in np.percentile(R_shell, [16, 50, 84])],
                              psr_plx_mas='3.47 +- 0.36 (Brisken+2003 via Thorsett+2003; ATNF btgg03 PX 3.5(4))',
                              centre_depth_offset_sigma_pc=15.0)
# near/far side band if the south spine is on the shell surface but projected inside a slightly larger shell
alpha = np.deg2rad(np.degrees(theta) + rng.uniform(0, 1.0, len(theta)))   # shell up to 1 deg larger than theta
Rr = Dc * np.sin(alpha)
disc = Rr ** 2 - (Dc * np.sin(theta)) ** 2
near = Dc * np.cos(theta) - np.sqrt(np.clip(disc, 0, None))
far = Dc * np.cos(theta) + np.sqrt(np.clip(disc, 0, None))
OUT['monogem_shell_near_far_mc'] = dict(near_pc_16_50_84=[round(float(x)) for x in np.percentile(near, [16, 50, 84])],
                                        far_pc_16_50_84=[round(float(x)) for x in np.percentile(far, [16, 50, 84])],
                                        note='shell angular radius drawn uniform in [theta, theta+1 deg]')

# ------------------------------------------------------------ 2. Local Bubble wall (O'Neill+2024) along the spine
with fits.open(ROOT / 'data/raw/context_dust/oneill2024_LB_ShellProperties_MeanDraws.fits', memmap=True) as h:
    t = h[1].data
    L, B = np.array(t['l'], float), np.array(t['b'], float)
    dpk, sdk, din, dout = (np.array(t[k], float) for k in ('d', 'sig_d', 'd_inner', 'd_outer'))


def unitvec(l, b):
    l, b = np.deg2rad(l), np.deg2rad(b)
    return np.c_[np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)]


tree = cKDTree(unitvec(L, B))
sg = spine.galactic
_, idx = tree.query(unitvec(sg.l.deg, sg.b.deg))
OUT['local_bubble_oneill24_spine'] = dict(peak_pc_range=[round(float(dpk[idx].min())), round(float(dpk[idx].max()))],
                                          peak_pc_median=round(float(np.median(dpk[idx]))),
                                          inner_pc_median=round(float(np.median(din[idx]))),
                                          outer_pc_median=round(float(np.median(dout[idx]))),
                                          sig_pc_median=round(float(np.median(sdk[idx]))),
                                          south_peak_median=round(float(np.median(dpk[idx][south]))),
                                          north_peak_median=round(float(np.median(dpk[idx][~south]))))

# ------------------------------------------------------------ 3. kinematics: LSR conversion and V_LSR(d)
gal_cen = spine[len(spine) // 2]


def helio_to_lsr(ra, dec, vh, std=True):
    """V_LSR - V_helio for the standard solar motion (20 km/s to RA 18h Dec +30 (B1900)) or Schoenrich+2010."""
    c = SkyCoord(ra * u.deg, dec * u.deg)
    if std:
        apex = SkyCoord(ra='18h', dec='30d', frame='fk4', equinox='B1900').icrs
        return vh + 20.0 * np.cos(c.separation(apex).rad)
    vsun = CartesianDifferential([11.1, 12.24, 7.25] * u.km / u.s)
    i = ICRS(ra=ra * u.deg, dec=dec * u.deg, distance=1 * u.pc, pm_ra_cosdec=0 * u.mas / u.yr,
             pm_dec=0 * u.mas / u.yr, radial_velocity=vh * u.km / u.s)
    return i.transform_to(LSR(v_bary=vsun)).radial_velocity.to(u.km / u.s).value


corr_std = np.array([helio_to_lsr(r, d, 0.0) for r, d in zip(sp_ra, sp_dec)])
OUT['lsr_correction_std_solar_motion_spine'] = [round(float(corr_std.min()), 2), round(float(corr_std.max()), 2)]
vel = {'[O III] DESI (+14.5 helio)': 14.5, '[O III] SDSS sky fibres (+12.5)': 12.5,
       '[O III] LAMOST Wu+2026 (+23.8..+26.8)': 25.3, 'HD 61717 Ca II (-16.88)': -16.88,
       'WHAM H-alpha main (helio ~+12)': 12.4}
cmid = float(np.median(corr_std))
OUT['vlsr'] = {k: round(v + cmid, 1) for k, v in vel.items()}
# HD 61717 specifically
c61 = SkyCoord(115.39431 * u.deg, 22.11742 * u.deg)
OUT['hd61717_vlsr_std'] = round(float(helio_to_lsr(115.39431, 22.11742, -16.88)), 2)
OUT['hd61717_vlsr_schoenrich'] = round(float(helio_to_lsr(115.39431, 22.11742, -16.88, std=False)), 2)
OUT['hd61717_lb'] = [round(float(c61.galactic.l.deg), 3), round(float(c61.galactic.b.deg), 3)]

R0, TH0 = 8.18, 236.0
dgrid = np.linspace(0.0, 3.0, 301)


def vlsr_rot(l, b, d):
    lr, br = np.deg2rad(l), np.deg2rad(b)
    R = np.sqrt(R0 ** 2 + (d * np.cos(br)) ** 2 - 2 * R0 * d * np.cos(br) * np.cos(lr))
    return (TH0 / R - TH0 / R0) * R0 * np.sin(lr) * np.cos(br)


lc, bc = sg.l.deg.mean(), sg.b.deg.mean()
vrot = np.array([vlsr_rot(sg.l.deg, sg.b.deg, d).mean() for d in dgrid])
SIG_PEC = 7.0   # km/s, typical cloud-to-cloud dispersion (assumption, stated)


def kin_range(v):
    ok = np.abs(vrot - v) < SIG_PEC
    return [round(float(dgrid[ok].min()), 2), round(float(dgrid[ok].max()), 2)] if ok.any() else None


OUT['kinematic'] = dict(dVdd_kms_per_kpc_at_1kpc=round(float(np.interp(1.0, dgrid, vrot)), 2),
                        sigma_pec_assumed=SIG_PEC,
                        ranges_kpc={k: kin_range(v) for k, v in OUT['vlsr'].items() if 'HD 61717' not in k},
                        note='flat rotation R0=8.18, Theta0=236; range = distances where |V_rot - V_obs| < sigma_pec')

# ------------------------------------------------------------ 4. pulsars (ATNF v2.9.0 on disk)
psr = {}
txt = (ROOT / 'data/raw/context_psr/psrcat_tar/psrcat.db').read_text(errors='replace')
for name in ('J0751+1807', 'J0659+1414'):
    blk = txt[txt.index('PSRJ     ' + name):]
    blk = blk[:blk.index('@-')]
    rec = {}
    for ln in blk.splitlines():
        parts = ln.split()
        if len(parts) >= 2 and parts[0] in ('RAJ', 'DECJ', 'DM', 'PX', 'RM', 'DIST_DM', 'DIST_A'):
            rec[parts[0]] = parts[1:]
    psr[name] = rec
pj = psr['J0751+1807']
ra_j = SkyCoord(pj['RAJ'][0], pj['DECJ'][0], unit=(u.hourangle, u.deg))
dsp = ra_j.separation(spine).deg
k = int(np.argmin(dsp))
pxv = float(pj['PX'][0])
pxe = 0.04   # ATNF error '4' in last quoted digit of 0.85 (abb+23)
OUT['pulsars'] = dict(
    J0751_1807=dict(px_mas=pxv, e_px=pxe, ref=pj['PX'][-1], d_kpc=round(1 / pxv, 3),
                    d_kpc_68=[round(1 / (pxv + pxe), 3), round(1 / (pxv - pxe), 3)], DM=float(pj['DM'][0]),
                    RM=float(pj['RM'][0]), RM_ref=pj['RM'][-1], min_sep_spine_deg=round(float(dsp[k]), 2),
                    dec_at_min=round(float(sp_dec[k]), 2), side='E (convex)' if ra_j.ra.deg > sp_ra[k] else 'W',
                    mean_ne=round(float(pj['DM'][0]) / (1000 / pxv), 4),
                    mean_Bpar_uG=round(1.232 * float(pj['RM'][0]) / float(pj['DM'][0]), 2),
                    note_vs_wave1='wave-1 used the 2017 VizieR B/psr snapshot: PX 1.6+-0.8 mas, no RM'),
    B0656_14=dict(px_mas=float(psr['J0659+1414']['PX'][0]), DM=float(psr['J0659+1414']['DM'][0]),
                  RM=float(psr['J0659+1414']['RM'][0]), d_pc_parallax='288 +33/-27 (3.47+-0.36 mas)'))
# expected Silk DM: from [O III] surface brightness, n_e N(O++) = 4 pi I / q5007
I_R = 2.4
T = 3.0e4
q_tot = 8.629e-6 * 2.29 / (9 * np.sqrt(T)) * np.exp(-2.514 / (8.617e-5 * T))   # cm^3/s, 1D2 excitation
q5007 = q_tot * 0.745    # branching to 5007 (5007/(5007+4959+6300...) ~ 0.745)
neNO2 = I_R * 1e6 / q5007
dm_rows = []
for ne in (0.3, 1, 3, 10, 100):
    NO2 = neNO2 / ne
    NH = NO2 / (4.9e-4 * 0.5)       # solar O/H, O++/O = 0.5 (assumption)
    dm_rows.append(dict(n_e=ne, N_O2=float(f'{NO2:.2e}'), N_H=float(f'{NH:.2e}'), L_pc=float(f'{NH / ne / 3.086e18:.3g}'),
                        DM_pc_cm3=float(f'{NH / 3.086e18:.3g}')))
OUT['silk_dm_estimate'] = dict(I_5007_R=I_R, T_e=T, ne_times_N_O2_cm5=float(f'{neNO2:.3e}'), rows=dm_rows,
                               note='DM = N_e/3.086e18; assumes O/H=4.9e-4, O++/O=0.5, edge-on path; any n_e gives DM << 30')

# ------------------------------------------------------------ 5. absorption + IUE summaries from this tag
h61 = json.loads((ROOT / 'data/derived/verify_distance_hd61717.json').read_text())
iue = json.loads((ROOT / 'data/derived/verify_distance_iue.json').read_text())
sep61 = c61.separation(spine).deg
OUT['hd61717_geometry'] = dict(min_sep_spine=round(float(sep61.min()), 2),
                               sep_two_snr_centre=round(float(c61.separation(CEN2).deg), 2),
                               sep_single_centre=round(float(c61.separation(CEN1).deg), 2),
                               spine_sep_two_snr_centre_at_same_dec=round(float(np.interp(22.12, sp_dec, th2)), 2),
                               sep_curv_centre=round(float(c61.separation(SkyCoord(94.5 * u.deg, 16.76 * u.deg)).deg), 2))
# projection factor if the -17 cloud were the approaching cap of a shell, HD 61717 delta inside the limb
for nm, cen, Rdeg in (('curvature', SkyCoord(94.5 * u.deg, 16.76 * u.deg), 21.33),):
    th = c61.separation(cen).deg
    fac = np.sqrt(max(1 - (np.sin(np.deg2rad(th)) / np.sin(np.deg2rad(Rdeg))) ** 2, 0))
    OUT['hd61717_geometry'][f'projection_factor_{nm}_shell'] = round(float(fac), 3)
    OUT['hd61717_geometry'][f'implied_vexp_{nm}_shell_kms'] = round(float(abs(OUT['hd61717_vlsr_std'] - 0.5) / fac), 0) if fac > 0 else None
# Gaia distance of HD 61717 with RUWE-inflated error
plx61, eplx61, ruwe61 = 3.6711, 0.0670, 2.143
e_infl = eplx61 * ruwe61
OUT['hd61717_distance'] = dict(plx=plx61, e_plx=eplx61, ruwe=ruwe61, e_plx_inflated=round(e_infl, 3),
                               d_pc=round(1000 / plx61, 1),
                               d_pc_68_inflated=[round(1000 / (plx61 + e_infl), 1), round(1000 / (plx61 - e_infl), 1)])

# ------------------------------------------------------------ 6. constraint table
rows = []


def row(name, kind, lo, hi, strength, basis, verdict):
    rows.append(dict(constraint=name, kind=kind, lo_pc=lo, hi_pc=hi, strength=strength, basis=basis, verdict=verdict))


ml = OUT['monogem_limb_mc']['d_limb_pc_percentiles_2p5_16_50_84_97p5']
nf = OUT['monogem_shell_near_far_mc']
row('Monogem limb (south spine at 12.1 deg from Knies centre)', 'range', ml[1], ml[3], 'strongest (morphological)',
    'PSR B0656+14 parallax x cos(theta); X-ray edge registration (eRASS1 W/E 1.32, 2.4 sigma)', 'favoured')
row('Monogem shell surface, near/far side', 'range', nf['near_pc_16_50_84'][0], nf['far_pc_16_50_84'][2], 'moderate',
    'shell up to 1 deg larger than spine radius', 'allowed')
lb = OUT['local_bubble_oneill24_spine']
row('Local Bubble wall (O\'Neill+2024)', 'reference', lb['inner_pc_median'], lb['outer_pc_median'], 'none (coincidence only)',
    'paper 160 pc sits on the wall; no dust/HI/X-ray signature ties the Silk to it', 'not a constraint')
row('HD 61717 -17 km/s Ca-rich cloud', 'upper (conditional)', None, OUT['hd61717_distance']['d_pc_68_inflated'][1], 'weak',
    'only if the cloud belongs to the Silk shell; velocity differs from the [O III] gas', 'conditional')
row('WD J0747+2106 COS: no Si IV/N V (0.16 deg)', 'lower (conditional)', 99.5, None, 'weak',
    'only if the Silk carries log N(Si IV) > 12.2 / N V > 13.0 in front of the WD', 'conditional')
row('Ca II / Na I at the [O III] velocity', 'none', None, None, 'none',
    'revised [O III] V_helio +13 equals the ubiquitous local gas; ionized gas need not carry Ca II/Na I', 'no constraint')
kr = OUT['kinematic']['ranges_kpc']
row('Kinematic, [O III] DESI/SDSS V_LSR ~ +0..+2', 'range', 0, 1000 * kr['[O III] DESI (+14.5 helio)'][1], 'none',
    f'dV/dd ~ {OUT["kinematic"]["dVdd_kms_per_kpc_at_1kpc"]} km/s/kpc << sigma_pec', 'meaningless')
row('Kinematic, LAMOST V_LSR +11..+16 (Wu+2026)', 'range', 1000 * kr['[O III] LAMOST Wu+2026 (+23.8..+26.8)'][0],
    1000 * kr['[O III] LAMOST Wu+2026 (+23.8..+26.8)'][1], 'none', 'velocity itself refuted by DESI/SDSS', 'meaningless')
row('PSR J0751+1807 DM (0.95 deg E, 1.18 kpc)', 'none', None, None, 'none',
    'not behind the filament; Silk DM << 1 pc cm^-3 for any n_e', 'no constraint')
row('Extinction / 3D dust', 'none', None, None, 'none', 'Silk has no dust; E(B-V)<0.006 to 300 pc; no 250-320 pc feature',
    'no constraint')
row('IUE HD 60848 (1.9 kpc, 2.3 deg W): high ions', 'none', None, None, 'none',
    'total Galactic C IV/Si IV only; no separable component, no |v|>70 km/s gas', 'no constraint')
row('Size / height plausibility', 'soft upper', None, 1000, 'soft',
    'd=1 kpc -> 295 pc long, 4.9 pc wide, z=330 pc; 1.6 kpc -> 470 pc long', 'disfavours >~1 kpc')
OUT['constraints'] = rows
OUT['combined'] = dict(favoured_pc=[ml[1], ml[3]], favoured_95_pc=[ml[0], ml[4]], permitted_pc=[100, 1000],
                       statement='favoured 255-315 pc (Monogem limb, south); permitted ~100 pc to ~1 kpc; 160 pc '
                                 '(Local Bubble wall) not excluded by any measurement, only by the X-ray-edge argument')
(ROOT / 'data/derived/verify_distance_synthesis.json').write_text(json.dumps(OUT, indent=1, default=float))
np.savez_compressed(ROOT / 'data/derived/verify_distance_synthesis_arrays.npz', dgrid=dgrid, vrot=vrot, d_limb=d_limb[:20000],
                    near=near[:20000], far=far[:20000], lb_peak=dpk[idx], lb_in=din[idx], lb_out=dout[idx], sp_dec=sp_dec,
                    th2=th2, th1=th1)
print(json.dumps({k: OUT[k] for k in OUT if k != 'constraints'}, indent=1, default=float))
