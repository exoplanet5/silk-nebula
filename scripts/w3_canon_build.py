"""w3_canon_build.py -- canonical values table for the Silk Nebula final report (wave 3, tag `canon`).

Builds data/derived/w3_canon_values.json from the wave-1 and wave-2 notes/JSON (wave 2 supersedes wave 1).
Every numeric entry quotes the source as written; where a machine-readable source exists, the script
re-reads it and records a `json_check` (found value, pass/fail within tolerance) so the table cannot drift
silently from the products. It also records whether the concurrent wave-3 notes (oii, northsouth) exist.

Run:  ~/.venvs/astro313/bin/python -I scripts/w3_canon_build.py
Reads only project files written by earlier agents (no network, no downloaded content).
"""
import csv
import datetime as dt
import json
import os

ROOT = '/Users/mickey/astronomy/silk-nebula'
DER = os.path.join(ROOT, 'data/derived')
OUT = os.path.join(DER, 'w3_canon_values.json')

_cache = {}


def jload(fname):
    if fname not in _cache:
        with open(os.path.join(DER, fname)) as f:
            _cache[fname] = json.load(f)
    return _cache[fname]


def jget(fname, path):
    x = jload(fname)
    for k in path:
        x = x[k]
    return x


def check(fname, path, expected, tol):
    """Return a provenance check record; never raises."""
    rec = {'file': 'data/derived/' + fname, 'path': '.'.join(str(p) for p in path), 'expected': expected}
    try:
        v = jget(fname, path)
        rec['found'] = v
        if isinstance(expected, str) or isinstance(v, str):
            rec['pass'] = (str(v) == str(expected))
        else:
            rec['pass'] = abs(float(v) - float(expected)) <= tol
    except Exception as e:  # noqa: BLE001
        rec['found'] = None
        rec['pass'] = False
        rec['error'] = repr(e)
    return rec


def csv_check(fname, key_col, key_val, col, expected):
    rec = {'file': 'data/derived/' + fname, 'row': f'{key_col}={key_val}', 'col': col, 'expected': expected}
    try:
        with open(os.path.join(DER, fname)) as f:
            for row in csv.DictReader(f):
                if row[key_col] == str(key_val):
                    rec['found'] = row[col]
                    rec['pass'] = row[col] == str(expected)
                    return rec
        rec['found'] = None
        rec['pass'] = False
    except Exception as e:  # noqa: BLE001
        rec['found'] = None
        rec['pass'] = False
        rec['error'] = repr(e)
    return rec


ENTRIES = []


def E(id_, section, quantity, adopted, unc='', units='', status='measurement', method='', prov=(),
      superseded=(), conf='medium', notes='', checks=()):
    ENTRIES.append({
        'id': id_, 'section': section, 'quantity': quantity, 'adopted': adopted, 'uncertainty': unc,
        'units': units, 'status': status, 'method': method, 'provenance': list(prov),
        'superseded': [dict(value=s[0], source=s[1], reason=s[2]) for s in superseded],
        'confidence': conf, 'notes': notes, 'json_check': list(checks),
    })


T = 'silk_trace.json'
# ============================================================== A. GEOMETRY
E('geom.spine', 'geometry', 'E-strand spine (ordered S->N)',
  '(116.5405, +12.00) -> (114.3622, +28.25); (l,b) (208.15, +17.50) -> (191.37, +21.85); 66 points at 0.25 deg Dec steps',
  'per-point sigma median 0.022 deg, max 0.20 deg (junction Dec 18.5-18.75)', 'deg (ICRS)', 'measurement',
  'ridge Gaussian fits on 4 registered images (SWONS Fig 1, Fig 4 background, NSNS narrow, NSNS wide)',
  ['data/derived/silk_trace.json (v1.0)', 'notes/geometry.md s4.2'], conf='high',
  notes='North of Dec 26.25 only NSNS-wide covers (single source). Fainter emission continues to Dec ~11.3 and ~28.9 (NSNS wide, S/N~7), not in the spine.',
  checks=[check(T, ['spine', 0, 'ra'], 116.5405, 1e-4), check(T, ['spine', -1, 'ra'], 114.3622, 1e-4),
          check(T, ['spine', 0, 'l'], 208.147, 1e-3), check(T, ['spine', -1, 'b'], 21.849, 1e-3)])
E('geom.length', 'geometry', 'Spine length', '16.94 deg arc (16.38 deg chord); bright south (Dec 12.0-21.5) 9.8 deg, faint north arc 7.2 deg',
  '~+-5% (curvature, end definition)', 'deg', 'measurement', 'path length of spine v1.0',
  ['notes/geometry.md s4.2', 'notes/discovery.md s5.3'],
  superseded=[('1 deg x 10 deg', 'Wu+2026 abstract/Fig 1 caption', 'describes only the bright southern part; the paper\'s own Fig 1 covers Dec 12.82-26.29 and shows the spine over 13.8 deg'),
              ('1.5 x 15 deg', 'QHYCCD press note 2024-07-12', 'coverage difference'),
              ('14 deg', 'NSNS K1a 2024-06', 'earlier, shallower NSNS data')],
  conf='high', notes='Physical arc length (derived here): 47.3 pc at 160 pc; 74.5/83.4/94.0 pc at 252/282/318 pc; 473 pc at 1.6 kpc.')
E('geom.width', 'geometry', 'Single-strand (E) width',
  'HWHM 0.138 deg median (16-84%: 0.076-0.191) -> FWHM ~0.28 deg; [O III] envelope (E+W+W2 +- HWHM) half-width median 0.223 deg, max 0.83 deg',
  '+-30% (JPEG stretch/filter dependent)', 'deg', 'measurement', 'Gaussian cross-sections on NSNS/SWONS renderings',
  ['data/derived/silk_trace.json: half_width_deg, envelope_half_width_deg', 'notes/geometry.md s4.2'],
  superseded=[('~1 deg width', 'Wu+2026', 'that is the multi-strand envelope, not a filament')], conf='medium',
  notes='FWHM 0.276 deg = 1.21/1.36/1.53 pc at 252/282/318 pc (derived). Shift-scan side lobes in SDSS imaging suggest real emission is broader or more structured than a 0.28-deg Gaussian (verify_precovery s2.1).',
  checks=[check(T, ['half_width_deg'], 0.138, 1e-4), check(T, ['envelope_half_width_deg', 'median'], 0.223, 1e-4),
          check(T, ['envelope_half_width_deg', 'max'], 0.83, 1e-4)])
E('geom.strands', 'geometry', 'Strands',
  'W strand: splits from E near Dec 15, bright to Dec ~19.75 (RA 116.49 -> 115.84), 4.84 deg long, 0.08-1.08 deg from E, ends at the Halpha patch, P2 lies on it. W2: RA 115.1-115.35, Dec 21.5-24.0. Junction: Dec 17.5-17.75 the E ridge steps 0.22 deg east (116.61 -> 116.83), faint (S/N 6-12) to Dec ~19.',
  'strand identification interpretive at Dec 15-19', '', 'measurement', 'ridge fits', ['notes/geometry.md s4.2', 'data/derived/silk_trace.json: strands'],
  conf='medium', notes='At Dec 15-20 the W strand carries most of the [O III]: SDSS ON-OFF W 1.93+-0.11 R vs E 0.73+-0.08 R (verify_sdss_results_base_r0.15.json); DESI W@18 node 3.88+-0.20 R.',
  checks=[check(T, ['strands', 'W', 'points', 0, 'ra'], 116.4866, 1e-4), check(T, ['strands', 'W', 'points', -1, 'ra'], 115.8371, 1e-4)])
E('geom.halpha_patch', 'geometry', 'Northern Halpha/[S II] patch',
  'half-max centroid (115.474, +20.536) = (l,b) (199.46, +20.00); 1.45 x 0.31 deg at PA 163; 0.31 deg^2; 30%-of-peak contour 2.24 x 1.28 deg; centroid 1.28 deg W of E spine',
  'threshold dependent', 'deg', 'measurement', 'NSNS narrow Halpha half-maximum contour',
  ['data/derived/silk_trace.json: halpha_patch_info', 'notes/geometry.md s4.4'],
  superseded=[('0.5 x 1 deg', 'Wu+2026', 'size depends on threshold'), ('3.06 x 2.0 deg at (115.67,+21.07)', 'NSNS O200.0a', 'catalogue ellipse')],
  conf='high', notes='PGCC G199.09+19.95 is 0.35 deg away with two CNM HI lines (V_LSR -0.74, +6.46 km/s); see mw.hi.pgcc.',
  checks=[check(T, ['halpha_patch_info', 'moments_halfmax', 'centroid_ra'], 115.474, 1e-3),
          check(T, ['halpha_patch_info', 'moments_halfmax', 'centroid_dec'], 20.536, 1e-3),
          check(T, ['halpha_patch_info', 'moments_halfmax', 'major_deg'], 1.45, 1e-3)])
E('geom.curvature', 'geometry', 'Curvature (single small circle, full E spine)',
  'centre (94.50, +16.76) = (l,b) (194.26, +0.44); R = 21.33 deg (bootstrap 20.48-22.35); chi2/dof 3744/63, rms 0.150 deg -> poor model (S-shaped)',
  'bootstrap ranges understate model error', 'deg', 'measurement (model-limited)', 'weighted small-circle fit, 400 bootstraps',
  ['data/derived/silk_trace.json: curvature', 'notes/geometry.md s4.3', 'notes/context.md s2'], conf='low',
  notes='Model-free sagitta 1.34 deg at Dec 19.75 bowing east (equiv. R ~25.8 deg). North (Dec>=17.75) R 21.8 (19.6-25.5); south (Dec<=17.5) R 45.9 (22.5-71.4), nearly straight. The curvature centre is set by the northern arc.',
  checks=[check(T, ['curvature', 'center_ra'], 94.5, 1e-3), check(T, ['curvature', 'center_dec'], 16.76, 1e-3),
          check(T, ['curvature', 'radius_deg'], 21.33, 1e-3)])
E('geom.monogem', 'geometry', 'Relation to Monogem centres',
  'South spine (Dec<17.6, 23 pts) at 12.10+-0.12 deg (rms) from the Knies+24 two-SNR centre (104.15,+15.30; R=11.00+-1.1 deg) and 11.36+-0.15 deg from PSR B0656+14; north bends to 12.39-16.03 deg. Single-SNR centre (107.08,+13.50; R=11.97): 9.08-16.23 deg, rim crossing at (116.71,+21.18), 9.25 deg of arc inside, 7.68 outside. Ring brightness peak (PA 60-100) at r~10.25 deg from the PSR: the Silk is ~1 deg outside it, on the outer flank.',
  '', 'deg', 'measurement', 'great-circle distances; eRASS1 radial profiles',
  ['data/derived/verify_distance_synthesis.json: monogem_geometry', 'notes/verify_xrayfuv.md s1', 'notes/context.md s1', 'notes/discovery.md s6'],
  conf='high', notes='Concentricity is with the Monogem RADIUS only in the south; the south section is short (5.5 deg) and nearly straight, so low rms alone does not single out a centre (context caveat). context.md gives +-0.16 deg for the PSR distance; verify_* give +-0.15.',
  checks=[check('verify_distance_synthesis.json', ['monogem_geometry', 'spine_sep_two_snr_centre_south', 0], 12.1, 1e-3),
          check('verify_distance_synthesis.json', ['monogem_geometry', 'spine_sep_two_snr_centre_south', 1], 0.12, 1e-3),
          check('verify_distance_synthesis.json', ['monogem_geometry', 'spine_sep_psr_south', 0], 11.36, 1e-3)])
E('geom.geminga', 'geometry', 'Relation to Geminga / pulsar birth loci',
  'Spine is 17.56+-0.22 deg (rms) from present-day Geminga (nearly concentric) but this is not physical (Geminga moves that rms in ~4.5 kyr; its spin-down can inflate only R~17 pc). Curvature centre is 3.9 deg from Geminga today, 12.6 deg from the Monogem single-SNR centre, 10.4 deg from PSR B0656+14, 0.8 deg from the Origem Loop centre (projection coincidence). No back-traced pulsar birth locus coincides with the curvature centre (B0656+14 birth 9.05 deg away; Geminga 13.0 deg).',
  '', 'deg', 'derived', 'LSR back-traces, 3000 MC draws', ['notes/context.md s2', 'data/derived/context_pulsar_coincidence.json'], conf='medium')
E('geom.paper_regions', 'geometry', 'Paper regions P1-P4 relative to strands',
  'P1 (116.65,+16.60) on E (-0.005 deg); P2 (116.21,+17.74) on W (-0.04 deg from W; -0.46 from E); P3 (116.45,+14.72) 0.115 deg W of the merged ridge; P4 (115.20,+16.18) 1.05-1.35 deg W of the strands, a valid background (NSNS [O III] -0.11 R). Fig 1b drawn outlines differ from Table 1 by 0.48 deg (P1) and 0.42 deg (P4).',
  '', 'deg', 'measurement', 'registered Fig 1 (227 Tycho-2 stars, rms 5.2")', ['notes/geometry.md s4.1', 'notes/imaging_mw.md (offset convention)'], conf='high',
  notes='The paper\'s "centre" (116.68,+16.42) is 3.8 deg S of the spine midpoint (discovery s5.3).')

# ============================================================== B. DISCOVERY / CREDIT
E('disc.credit', 'discovery', 'Independent discoveries and credit',
  'At least three independent detections: SWONS (2023-11 per paper; 2023-12 per G.-Y. Sun blog), S. Ziegenbalg / NSNS (~2024-01, self-reported), Schaeffer-Fesen-Drechsler / NHZ (undated; published as Monogem NE limb in Fesen+2026, arXiv 2026-08-27). Koposov/DESI claim is hearsay (S. Jackson AstroBin 2025-05-02 quoting Belokurov).',
  '', '', 'literature', 'timestamped web/API records', ['notes/discovery.md s0, s4', 'data/derived/discovery_timeline.json'],
  conf='high', notes='AJ version (2026-11-04) adds the co-discoverer acknowledgment absent from arXiv v1 but keeps "discovered and named by the SWONS team" and "first time to its full extent". DESI DR1 public spectra do show the filament at ~30 sigma (verify_desi), so the Koposov claim is technically plausible but unverified.')
E('disc.first_public', 'discovery', 'First public records',
  'First public image: AstroBin zh47r9 (Song-Yu Jia, SWONS) 2024-01-14. First public post with position + multi-line data: Ziegenbalg, Cloudy Nights, 2024-06-11. First English use of the name "Silk Nebula": QHYCCD note 2024-07-12 (Chinese name in China National Astronomy 2024 No. 4, not seen). First journal paper: Ziegenbalg 2025, RNAAS 9, 227 (2025-08-28), 14 months before Wu+ on arXiv (2026-10-09).',
  'self-reported detection dates are M-reliability', '', 'literature', '', ['notes/discovery.md s0, s3'], conf='high')

# ============================================================== D. [O III] BRIGHTNESS
E('sb.desi', 'brightness', '[O III]5007 ridge peak, DESI DR1 (fibre scale)',
  'south (Dec<19.5) median 1.62 R (variants 1.60-1.71; box global amplitude 1.62); north 0.56 R (node median) / 0.71 R (box, Dec>=19.5); E nodes 1.4-2.3 R at Dec 14-18; W 1.15-1.62 R at Dec 15-17; W@18 3.88+-0.20 R; best model-free bin 1.81+-0.05 R -> >=1.7 R (2 sigma)',
  '+-8% systematic (GFA acceptance, chromatic model)', 'R', 'measurement', 'petal forward model + box integrals; GFA point-source acceptance f_ps (median 0.566 at 5008 A)',
  ['notes/verify_desi.md "Bottom line" 2', 'data/derived/verify_desi_summary.json: brightness', 'data/derived/verify_desi_audit_extra.json'],
  superseded=[('south 2.4 R; north 0.86 R; strict lower limit 3.05+-0.14 R (>=2.8 R)', 'notes/spectra_desi.md', 'wave-1 used MEAN_PSF_TO_FIBER_SPECFLUX (Legacy 1"-seeing normalisation) as fibre acceptance; surface brightness too high by 1/0.736 (1.13-1.82 per fibre)')],
  conf='high', notes='Differential measurement: emission uniform on petal scales (>~1 deg) is invisible; sky-fibre positions not public (petal-mean proxy).',
  checks=[check('verify_desi_summary.json', ['brightness', 'south_median', 'gfa_chromatic'], 1.62, 0.01),
          check('verify_desi_summary.json', ['brightness', 'north_median', 'gfa_chromatic'], 0.56, 0.01),
          check('verify_desi_audit_extra.json', ['south_dec<19.5', 'amp_R_gfa'], 1.62, 0.01),
          check('verify_desi_audit_extra.json', ['north_dec>=19.5', 'amp_R_gfa'], 0.71, 0.01)])
E('sb.sdss_fibres', 'brightness', '[O III]5007, SDSS/BOSS/SDSS-V sky fibres (total sky)',
  'south ON-OFF 1.39+-0.07 (stat)+-0.07 (placebo) R (r_on 0.15 deg; 1.61 at r_on 0.10, 1.06 at 0.25); W strand Dec 15-20 1.93+-0.11 R; per-plate strand-profile peaks 1.6-2.5 R south, 0.4-0.8 R north; absolute P1/P2/P3 1.49+-0.06 / 1.78+-0.05 / 1.87+-0.09 R (+-0.2 R zero point); brightest single fibres 3.0-3.65 R',
  'f_eff 0.84+-0.08 (3") / 0.60+-0.06 (2")', 'R', 'measurement', 'ON-OFF of flux+sky per plate-spectrograph group; imaging-sky-anchored diffuse correction',
  ['notes/verify_sdss.md s1, s3', 'data/derived/verify_sdss_results_base_r0.15.json'],
  superseded=[('south ON-OFF 0.79+-0.02 R; W strand 1.19 R', 'notes/spectra_sdss_lamost.md', 'wave-1 used the pipeline residual `flux` only; the per-fibre sky model absorbs 36-53% of the strand (residual keeps 0.47 BOSS / 0.64 SDSS of the signal)'),
              ('f_psf 0.93 (3") / 0.66 (2")', 'notes/spectra_sdss_lamost.md', 'imaging-sky anchor gives 0.84/0.60')],
  conf='high', checks=[check('verify_sdss_results_base_r0.15.json', ['regions', 'south_dec_lt20', 'OIII:OIII5007', 'I_R'], 1.39, 0.01),
                       check('verify_sdss_results_base_r0.15.json', ['regions', 'W_strand_15_20', 'OIII:OIII5007', 'I_R'], 1.93, 0.01)])
E('sb.sdss_imaging', 'brightness', '[O III]5007 ridge peak from SDSS imaging sky levels (2002-2007, 17-arcmin field scale)',
  'g PSP 2.31+-0.38 R (0.28-deg ridge template, 19 core runs Dec 12-20.5); Frames/FramesSub 2.4-2.6 R (preferred, less model-dependent); lower bound with all confounder regressors 1.65+-0.31 R; 10-13 R excluded at >=15 sigma on 17-arcmin scales',
  '~25% template systematic; PSP temporal smoothing', 'R', 'measurement', 'field-level model per run, NSNS-free template',
  ['notes/verify_precovery.md s0, s2.1, s2.4'],
  superseded=[('A = 1.35+-0.19 x NSNS; ridge medians 2.8-2.9 R; 1-arcmin peaks 4.0-4.9 R', 'notes/imaging_optical.md s3.5', '1-arcmin peaks were NSNS morphology scaled, not SDSS; 17-arcmin peak 2.3-2.6 R')],
  conf='medium')
E('sb.nsns', 'brightness', '[O III]5007, NSNS DR0.2 (calibrated HiPS; not re-verified in wave 2)',
  'on-ridge excess P1 2.04, P2 2.06, P3 2.18 R; 1-arcmin peaks 2.95 / 3.64 / 3.37 R; south_bright 1.60 R; mid (Dec 17-19) 1.21 R; northern arc 0.58-0.62 R; Ziegenbalg own estimates K1a 2.04 R, K1a\' 1.53 R',
  'stat 0.01-0.09 R; E/W background systematic 0.05-0.15 R', 'R', 'measurement (wave 1)', 'NSNS DR0.2 HiPS order 3 (51.5")',
  ['notes/imaging_optical.md s3.1', 'notes/discovery.md s5.4'], conf='medium')
E('sb.other', 'brightness', '[O III]5007, other estimates',
  'LAMOST LRS stellar-residual ON-OFF 0.88+-0.13 R (wave 1, not re-verified); LAMOST MRS-N relative x sky-line scale 2.8-3.2 R (APO [N I]/OH ratio) or 14-16 R (UVES ratio), indeterminate; QHYCCD 24.3 mag/arcsec^2 -> 3.3-6.7 R; DSS2 J < 14 x NSNS (3 sigma, ~<40 R)',
  '', 'R', 'measurement/derived (wave 1)', '', ['notes/spectra_sdss_lamost.md s4.1', 'notes/imaging_optical.md s1', 'notes/discovery.md s5.4'], conf='low')
E('sb.paper', 'brightness', '[O III]5007, Wu+2026 Table 3 (SWONS imaging)',
  '23.60 / 23.87 / 23.85 mag/arcsec^2 -> 12.75 / 9.94 / 10.12 R (P1/P2/P3; AB, 60 A top-hat)',
  'conversion 8.3-14.9 R (50-70 A), 9.1-14.0 (AB/Vega +-0.1 mag), 7.4-9.5 (full 4959 leakage)', 'R', 'literature (arithmetic reproduced)',
  '', ['notes/verify_paper.md s6', 'data/derived/verify_paper_energetics.json'], conf='high',
  notes='Arithmetic right; value ~5-8x above every independent measurement. Stacking all multiplicative systematics gives only ~5.5-7.5 R; the f/0.95 passband shift would make SWONS UNDER-estimate. Needs an additive background/flat error of 7-12% of a 21.5 mag sky (11-20% of 22.0).')
E('sb.adopted', 'brightness', 'ADOPTED [O III]5007 surface brightness',
  'Southern ridge peak ~2 R (range 1.4-2.5 R on 0.1-0.3 deg scales; DESI 1.62 R is the reference flux-calibrated value); north 0.4-0.9 R. Paper 10-13 R is 4-8x too high.',
  'DESI vs SDSS-imaging scale differs by 1.49+-0.25 (2 sigma), unresolved', 'R', 'adopted',
  'DESI (1.62 S / 0.56-0.71 N), SDSS fibres (1.39 S ON-OFF; peaks 1.6-2.5 S / 0.4-0.8 N), SDSS imaging (2.3-2.6 on 17\'), NSNS (~2.0-2.2 on-ridge)',
  ['notes/verify_desi.md', 'notes/verify_sdss.md', 'notes/verify_precovery.md s2.4', 'notes/imaging_optical.md'],
  superseded=[('~1-4 R (provisional headline)', 'notes/RESUME.md', 'narrowed after calibration fixes'), ('10-13 R', 'Wu+2026', 'see sb.paper')],
  conf='high', notes='No variability: fixed-position SDSS imaging 2003.97-2007.14 Delta = +0.8+-1.3 R; DESI-normalised trend -0.16+-0.13 /yr over 2002.9-2007.1; <~x2 change over 2002-2022.')

# ============================================================== E. VELOCITIES
E('v.lamost_paper', 'velocity', '[O III] V_helio, Wu+2026 Table 4 (LAMOST MRS-N)',
  'P1 +26.8+-4.3, P2 +26.8+-4.1, P3 +23.8+-3.3 km/s (V_LSR 15.5/14.5/11.2)', '', 'km/s', 'literature', '',
  ['notes/paper_review.md s3'], conf='high',
  notes='REFUTED as heliocentric: topocentric, air, airglow-referenced scale. Corrected by -10.1+-1.0 km/s -> P1 +16.7, P2 +16.7, P3 +13.7 (+-3-4) km/s helio (V_LSR +4.6/+4.6/+1.1 with local corrections -12.12/-12.12/-12.56; the velocity note quotes +4.5/+4.5/+1.5 using a single -12.19).',
  checks=[check('verify_velocity_summary.json', ['table4_offset'], 10.09, 0.01), check('verify_velocity_summary.json', ['table4_corrected', 'P1'], 16.7, 0.01)])
E('v.lamost_fig3a', 'velocity', '[O III] V_helio from the paper\'s own Fig. 3a vector spectrum vs its [N I] airglow',
  '+13.6+-1.2 km/s (V_topo +21.9; v_corr -8.32)', '', 'km/s', 'measurement', 'PyMuPDF polylines (3590 px, 0.106 A); calibration-free against [N I] 5198/5200',
  ['notes/verify_velocity.md s1, s4', 'data/derived/verify_velocity_summary.json: fig3a'], conf='high',
  notes='One spectrum of unknown location/co-addition; correction applied in aggregate to Table 4.',
  checks=[check('verify_velocity_summary.json', ['fig3a', 'V_helio_corr'], 13.6, 0.05), check('verify_velocity_summary.json', ['fig3a', 'err'], 1.2, 0.05)])
E('v.desi', 'velocity', '[O III] V_helio, DESI DR1 (barycentric = helio to <0.02 km/s)',
  'verify_velocity on-ridge stack (390 fibres, predicted residual >0.5 R) +14.4+-1.2 (+-1.9 incl. +-1.5 zero point); south (Dec<17.6) +13.1+-1.2; negative-residual fibres +12.6+-1.1. verify_desi: IV regression +11.0 (10.3-14.3); positive fibres S/N>4 +12.2+-0.8 (n=346); rev-2 joint regression on both [O III] lines +11.3 (10.2-13.9); Dec<19.5 +10.6 (10.2-12.0); Dec>=19.5 +20.9 (18.5-23.1)',
  'zero point +-1.5 (stars vs Gaia +1..+3; two backup tiles -19, unused)', 'km/s', 'measurement', '',
  ['notes/verify_velocity.md s3', 'notes/verify_desi.md "Bottom line" 4, audit', 'data/derived/verify_velocity_desi_sdss.json', 'data/derived/verify_desi_audit_vel.json', 'data/derived/verify_desi_velocity.json'],
  superseded=[('+14.5 (13-16)', 'notes/spectra_desi.md', 'confirmed within errors; the estimators differ by ~3 km/s')], conf='high',
  checks=[check('verify_velocity_desi_sdss.json', ['desi', 'pred_pos_gt0.5R', 'V_helio'], 14.4, 0.05),
          check('verify_desi_audit_vel.json', ['total_all', 'Z', 'v'], 11.3, 0.05),
          check('verify_desi_audit_vel.json', ['total_dec<19.5', 'Z', 'v'], 10.6, 0.05),
          check('verify_desi_audit_vel.json', ['total_dec>=19.5', 'Z', 'v'], 20.9, 0.05),
          check('verify_desi_velocity.json', ['regspec_centroids', 'all', 'per_fibre_pos_SN4', 'wmean'], 12.2, 0.05)])
E('v.sdss', 'velocity', '[O III] V_helio, SDSS sky fibres',
  'verify_sdss: south +10.0+-1.8 (stat)+-0.8 (zp) (SDSS-I/II +8.7+-2.6; BOSS +11.2+-2.5); all +12.3+-2.5; north +16.4+-5.9; W +9.7+-2.2; E +11.3+-4.1. verify_velocity own: all +9.6+-2.8; south +7.5+-3.4; north +15.5+-3.7 (variants +6.8..+10.5)',
  '', 'km/s', 'measurement', 'free-velocity fit to group-bootstrapped ON-OFF total-sky stacks, [O I]5577 zero point',
  ['notes/verify_sdss.md s5', 'notes/verify_velocity.md s3', 'data/derived/verify_sdss_results_base_r0.15.json', 'data/derived/verify_velocity_desi_sdss.json'],
  superseded=[('+12.5+-2.5+-0.7', 'notes/spectra_sdss_lamost.md', 'confirmed within errors')], conf='high',
  checks=[check('verify_sdss_results_base_r0.15.json', ['regions', 'south_dec_lt20', 'velocity_OIII', 'v_corr'], 10.0, 0.05),
          check('verify_sdss_results_base_r0.15.json', ['regions', 'north_dec_ge20', 'velocity_OIII', 'v_corr'], 16.4, 0.05),
          check('verify_velocity_desi_sdss.json', ['sdss', 'all', 'V_helio_5007'], 9.6, 0.05)])
E('v.adopted', 'velocity', 'ADOPTED [O III] velocity (whole filament, south-dominated)',
  'V_helio = +13.3+-1.0 km/s; V_LSR = +1.2+-1.0 (standard solar motion), +1.5 (Schoenrich+10)',
  'add ~+-1 km/s estimator systematic: using verify_desi rev-2 (+11.3) and verify_sdss (+10.0) instead of +14.4/+9.6 gives +12.4+-0.9 (derived here)', 'km/s', 'adopted',
  'weighted mean of DESI +14.4+-1.9, LAMOST Fig 3a +13.6+-1.2, SDSS +9.6+-2.8 (chi2 2.2/2 dof)',
  ['notes/verify_velocity.md s1', 'data/derived/verify_velocity_summary.json: combined'],
  superseded=[('+23.8..+26.8 helio / V_LSR +11..+16', 'Wu+2026 Table 4', 'calibration artefact (-10.1+-1.0)'),
              ('+12.5 to +14.5 helio', 'notes/RESUME.md (wave-1 provisional)', 'refined')],
  conf='high', notes='Equal to the WIM (WHAM V_LSR +1.3 on-filament south), the dominant local HI (V_LSR +1.2..+2.5, Dec 12-22) and the local Ca II/Na I group (V_helio +6..+15). No kinematic distance is possible. V_LSR - V_helio runs from -13.2 (Dec 12) to -9.3 (Dec 28.25) along the spine (derived; -12.12 at P1).',
  checks=[check('verify_velocity_summary.json', ['combined', 'V_helio'], 13.3, 0.05), check('verify_velocity_summary.json', ['combined', 'err'], 1.0, 0.05),
          check('verify_velocity_summary.json', ['combined', 'V_LSR_standard'], 1.2, 0.06)])
E('v.north_south', 'velocity', '[O III] north vs south',
  'DESI Dec>=19.5 +20.9 (18.5-23.1) vs Dec<19.5 +10.6 (10.2-12.0) km/s helio (Delta 10.3); SDSS north +16.4+-5.9 (verify_sdss) / +15.5+-3.7 (verify_velocity) vs south +10.0+-1.8. In V_LSR (spine-mean corrections -10.3 north, -12.3 south; derived): north ~+10.6 (8.2-12.8), south ~-1.7 (-2.1..-0.3) -> Delta V_LSR ~12 km/s.',
  '', 'km/s', 'measurement (gradient PENDING northsouth)', '', ['notes/verify_desi.md audit velocity', 'notes/verify_sdss.md s5'],
  conf='medium', notes='Canon cross-match (not tested statistically): northern [O III] V_LSR ~+10.6 coincides with the northern HI component that becomes dominant north of Dec ~23 (V_LSR +8.7..+11.3, FWHM 10-16 km/s; verify_hiwham s2.4), as the southern value coincides with the southern HI (+1.2..+2.5). verify_hiwham wrote "no [O III] velocity data exist there" before verify_desi rev 2 measured the north.',
  checks=[check('verify_desi_audit_vel.json', ['total_dec>=19.5', 'Z', 'v_boot_p16_p84', 0], 18.54, 0.01)])
E('v.east_west', 'velocity', '[O III] E vs W strand',
  'No significant difference: DESI rev-2 joint regression E +12.0 (11.3-14.3), W +11.4 (8.6-15.6); DESI positive fibres E +15.0+-1.0, W +13.7+-1.5; SDSS W +9.7+-2.2, E +11.3+-4.1. Residual hint: Dec 15.5-17.5 E < W by 3-6 km/s (<=2 sigma; positive fibres E +10.6+-2.2, W +16.7+-2.1).',
  '', 'km/s', 'measurement', '', ['notes/verify_desi.md (d2)', 'notes/verify_sdss.md s5'],
  superseded=[('E +13.0, W +18.5', 'notes/spectra_desi.md', 'regression centroids of local subsets are not local (petal-mean subtraction)'),
              ('SDSS W vs DESI W tension 3.5 sigma', 'notes/verify_sdss.md', 'resolved by the DESI revision (W +11.4)')], conf='medium',
  checks=[check('verify_desi_audit_vel.json', ['joint_all', 'Z_E', 'v'], 12.0, 0.05), check('verify_desi_audit_vel.json', ['joint_all', 'Z_W', 'v'], 11.4, 0.05)])
E('v.local_P', 'velocity', 'DESI [O III] near the LAMOST regions (positive fibres)',
  'within 0.5 deg: P1 +11.9+-1.9 (n=61), P2 +15.2+-2.9 (n=9), P3 +9.2+-1.0 (n=128) km/s; within 0.3 deg of P1 (verify_velocity, wave-1 fibres): +7.5+-2.2 (n=66)',
  '', 'km/s', 'measurement', '', ['notes/verify_desi.md (d)', 'notes/verify_velocity.md s4.3', 'data/derived/verify_desi_velocity.json', 'data/derived/verify_velocity_desi_sdss.json'],
  conf='medium', notes='LAMOST-corrected minus DESI(0.5 deg): P1 +4.8+-4.8 (1.0 sigma), P2 +1.5+-5.1 (0.3), P3 +4.5+-3.6 (1.3) -- consistent. With the 0.3-deg P1 value: +9.2+-4.9 (1.9 sigma). Regression bootstraps collapse for these 1-2 block subsets.',
  checks=[check('verify_desi_velocity.json', ['regspec_centroids', 'within0.5deg_P1', 'per_fibre_pos_SN4', 'wmean'], 11.9, 0.05),
          check('verify_velocity_desi_sdss.json', ['desi', 'within0.3_P1_pos', 'V_helio'], 7.5, 0.06)])
E('v.wim_hi', 'velocity', 'Comparison gas',
  'WHAM Halpha main component: south on/off V_LSR +1.28/+0.54, north -0.14/+0.77, FWHM_int 30-32 km/s, same on and off (paired on-off Delta v +0.7+-1.8 S, -0.9+-1.7 N). Paper red arm [N II]/[S II] V_helio +11.7..+13.3 (V_LSR -1..+1.7), agrees with WHAM. Paper Halpha +16.3..+19.7 helio is 4-7 km/s redder: reproduced by ~33% geocoronal over-subtraction (PLAUSIBLE, model only). On-spine HI peaks V_LSR +0.4..+5.5 (median +1.8) at Dec 12-22.',
  '', 'km/s', 'measurement', '', ['notes/verify_hiwham.md s1, s3.2, s2.4', 'notes/verify_velocity.md s3'], conf='high')
E('v.width', 'velocity', '[O III] intrinsic FWHM',
  'DESI <=33 km/s (32.5, stat 28-36; 26-38 for +-2% LSF; on-ridge positive fibres 18 (0-27)); paper Fig 3a example ~50 km/s; paper Table 2 43-49 km/s (not a measurement: fits passing the S/N-6 cut return 36-47 for any true width); shock model predicts 23-38 km/s',
  'DESI LSF 98 km/s FWHM must be right to 2-3%', 'km/s', 'measurement (conflicting)', '',
  ['notes/spectra_desi.md s6 (not re-verified in wave 2)', 'notes/verify_velocity.md s4.1', 'notes/paper_review.md s4.5', 'notes/verify_shockmodel.md s7'], conf='low')

# ============================================================== F. LINE RATIOS (energy unless stated)
E('r.o3hb', 'ratios', '[O III]5007/Hbeta (south / filament)',
  'DESI: lower limit only, > 11 (box; JSON 10.53) / > 18 (Gaussian; 17.66) at 2 sigma (energy); Hb/[O III] = 0.012+-0.041 (box) / 0.021+-0.018 (Gaussian), null-inflated. SDSS south: 23.5 photon (16-84%: 12.7-47.6; median 21.6; ~22.8 energy); 2 sigma > 10 (bootstrap), > 7.5 (+placebo MAD), > 5.3 (+placebo rms). SDSS W strand 22.6 (13.6-42.2). Central value UNDETERMINED.',
  '', '', 'lower limit', 'IV regression on geometry template (DESI); ON-OFF (SDSS)',
  ['notes/verify_desi.md audit', 'data/derived/verify_desi_audit_extra.json: limits', 'notes/verify_sdss.md row 3', 'data/derived/verify_sdss_results_base_r0.15.json'],
  superseded=[('47 (>24); south 46 (>28)', 'verify_desi first pass (still in verify_desi_summary.json: ratios)', 'WITHDRAWN: uninflated bootstrap; Hb only 1-2 sigma'),
              ('20 (15-26), >11.1 at 97.5%', 'notes/spectra_desi.md', 'lower limit confirmed; central value not supported'),
              ('>13.8 (2 sigma)', 'notes/spectra_sdss_lamost.md', 'placebo systematics -> >5-7.5')],
  conf='medium', notes='altmodel caveat: in verify_desi_ratios Hgamma/Hb = 1.5 (case B 0.47) -> Hb carries ~+-0.01-0.02 systematic; if Hgamma is right, [O III]/Hb ~15. Complete shocks give <=5.',
  checks=[check('verify_desi_audit_extra.json', ['limits', 'box', 'OIII_over_Hb_lower_2sig'], 10.53, 0.01),
          check('verify_desi_audit_extra.json', ['limits', 'gauss', 'OIII_over_Hb_lower_2sig'], 17.66, 0.01),
          check('verify_sdss_results_base_r0.15.json', ['regions', 'south_dec_lt20', 'ratios_photon', 'OIII5007/Hb', 'value'], 23.5, 0.01)])
E('r.r23', 'ratios', 'R23 = ([O II]3727 + [O III]4959,5007)/Hbeta',
  'DESI > 29 (box; JSON 28.97) / > 34 (Gaussian; 34.03) at 2 sigma; SDSS south ~76 (energy; > 20 at 2 sigma)', '', '', 'lower limit', '',
  ['data/derived/verify_desi_audit_extra.json: limits', 'notes/verify_sdss.md row 3h'],
  superseded=[('~128 (>84)', 'verify_desi first pass (still in verify_desi_summary.json)', 'WITHDRAWN'), ('58 (>35)', 'notes/spectra_desi.md', 'lower limit slightly weaker after null inflation')],
  conf='medium', checks=[check('verify_desi_audit_extra.json', ['limits', 'box', 'R23_lower_2sig'], 28.97, 0.01)])
E('r.o2o3', 'ratios', '[O II]3726+3729 / [O III]5007 (energy)',
  'DESI all 1.37 (Gaussian; 68% 1.21-2.07) / 1.41 (box; 1.27-2.03); DESI south (Dec<19.5) 1.11 (1.07-1.18) / 1.16 (1.09-1.24), +-0.08; DESI north 4.4 (3.5-5.2) / 4.2. SDSS (BOSS only) south 1.84 (95%: 1.46-2.37; +-0.35 with placebo), W 2.04, E 2.01, north 7.3. SDSS imaging u/g implies ~0.8-1.0 (qualitative only).',
  '', '', 'measurement (DESI vs SDSS south in ~2-3 sigma conflict)', '',
  ['notes/verify_desi.md c2/audit', 'data/derived/verify_desi_summary.json: ratios', 'notes/verify_sdss.md s6', 'notes/verify_precovery.md s2.2'],
  superseded=[('1.49 (1.33-1.75); 1.12 (Dec<19); north arc 6.2', 'notes/spectra_desi.md', 'confirmed; -6% chromatic acceptance correction')],
  conf='medium', notes='Rises steeply to the north (x4) in both instruments. Reddening (E(B-V) 0.03-0.05) would raise by ~5%; not applied.',
  checks=[check('verify_desi_summary.json', ['ratios', 'all', 'OII3727total/OIII5007', 'value'], 1.37, 0.01),
          check('verify_desi_summary.json', ['ratios', 'south_dec<19.5', 'OII3727total/OIII5007', 'value'], 1.11, 0.01)])
E('r.heii', 'ratios', 'He II 4686 / [O III]5007', 'DESI < 0.026 (97.5%; box < 0.013); SDSS < 0.04 (bootstrap), < 0.08 with placebo', '', '', 'upper limit', '',
  ['notes/verify_desi.md c4', 'notes/verify_sdss.md row 3j'], superseded=[('< 0.034', 'notes/spectra_desi.md', 'tightened')], conf='high')
E('r.4363', 'ratios', '[O III]4363 / 5007',
  'DESI 0.083 (Gaussian IV, +-0.029 inflated), 0.091 (box, +-0.051), 0.057 (pixel regression; 16-84% 0.015-0.069); south 0.092; ~2-3 sigma (1.8 box, 2.9 Gaussian); centre 4364.55-4364.73 A vs 4364.62 [O III], Hg I 4359.55 excluded (box Hg coefficient 0.004+-0.019). SDSS 0.038+-0.035 / 0.044+-0.037 (energy; non-constraining).',
  'null box at 4381.6 A gives 0.048+-0.014 -> spurious ~0.05 possible', '', 'hint (2-3 sigma)', '',
  ['notes/verify_desi.md c5/audit', 'data/derived/verify_desi_audit_extra.json: limits', 'notes/verify_sdss.md row 3i'],
  superseded=[('0.056 (0.022-0.077), ~2 sigma', 'notes/spectra_desi.md', 'value revised up'), ('3-5 sigma', 'verify_desi first pass', 'downgraded by null calibration')],
  conf='low', checks=[check('verify_desi_audit_extra.json', ['limits', 'gauss', 'OIII4363'], 0.083, 0.001), check('verify_desi_audit_extra.json', ['limits', 'box', 'OIII4363'], 0.091, 0.001)])
E('r.7325', 'ratios', '[O II]7320+7330 / [O III]5007',
  'DESI 0.036 (0.033-0.042) masked Gaussian fit (preferred); box 0.081 is OH-contaminated. -> [O II]7325/3727 = 0.026 (0.017-0.038)', '', '',
  'measurement, but below DESI ~0.05 template-systematic floor -> treat as upper limit',
  '', ['notes/verify_desi.md c11/audit', 'notes/verify_altmodel.md s2'], superseded=[('0.037 (0.033-0.044)', 'notes/spectra_desi.md', 'confirmed')], conf='low',
  notes='verify_desi new finding 3a: "any DESI ratio below ~0.05 of [O III] is effectively an upper limit". The shockmodel x1.8 tension (model 0.059-0.068) is weaker if 0.036 is only an upper limit.')
E('r.nelines', 'ratios', '[Ne III]3869, [S III]9069, [Ar III]7136, [N II]5755 relative to [O III]5007',
  '[Ne III] 0.099 (0.057-0.109) Gaussian / 0.069 (0.033-0.087) box; [S III]9069 0.028 (0.024-0.048) Gaussian (chromatic +12%) / 0.037 box; [Ar III] 0.013 (0.011-0.020); [N II]5755/6583 0.14 (0.06-0.26, ~1 sigma)',
  '', '', '[Ne III] measurement; others below ~0.05 floor (upper-limit-like)', '', ['notes/verify_desi.md c6, c7, c11, audit'],
  superseded=[('[Ne III] 0.10; [S III] 0.031; [Ar III] 0.018', 'notes/spectra_desi.md', 'confirmed within errors')], conf='low',
  notes='[N II]5755/6583 gives no temperature (at/above the infinite-T limit; altmodel s4.1).')
E('r.lowion_filament', 'ratios', 'Filament [N II]6583/Halpha and [S II]/Halpha',
  'DESI direct: [N II]/Ha 1.05 (box) - 1.21 (Gaussian; 1.12-1.43), south 1.0-1.08; [S II]/Ha 0.70-0.75, south 0.51-0.52. SDSS south [N II]/Ha 0.70 (0.55-0.92), [S II]/Ha 0.62 (0.48-0.82); north 1.75 (1.3-2.6; 1.1-2.4 across variants; 2 sigma >~0.7), [S II]/Ha 1.86 (1.4-2.7). SYSTEMATICS-LIMITED (DESI Ha Z>0/Z<0 asymmetry; SDSS [S II]6716/6731 ON-OFF 2.5/1.8 unphysical).',
  'factor ~2', '', 'systematics-limited', '', ['notes/verify_desi.md c9', 'notes/verify_sdss.md s3, s6'],
  superseded=[('DESI [N II]/Ha 0.5-0.7, [S II]/Ha 0.3-0.55', 'notes/spectra_desi.md', 'REFUTED: came from Hb rescaling, Hb only 1-2 sigma'),
              ('SDSS north [N II]/Ha 1.85+-0.36, [S II]/Ha 1.77+-0.37', 'notes/spectra_sdss_lamost.md', 'errors omit the dominant Halpha systematic')], conf='low')
E('r.o3ha', 'ratios', '[O III]5007/Halpha (south strand)',
  'Best lower limit >~3 (2 sigma, photon; >~4 energy) co-spatial: SDSS Dec<17.6 with placebo > 3.0; DESI direct > 3.4 (raw p97.5) to > 8.8 (south); SDSS Dec<20 2.7 (2.2-3.4; +placebo 1.9-4.5). Non-co-spatial (Ha displaced up to 0.8 deg): SHASSA > 1.9 / 2.3 / 4.7 for [O III] 1.6 / 2 / 4 R. DESI nominal ~10 photon (13.7 energy).',
  '', '', 'lower limit', '', ['notes/verify_hiwham.md s4', 'data/derived/verify_hiwham_wham_summary.json: oiii_ha_lower_limits'],
  superseded=[('Halpha/[O III] < 0.10 (2 sigma) at the [O III] velocity', 'notes/spectra_ism.md', 'REFUTED: tested at the wrong velocity (V_LSR +13.6); at V_LSR ~+1 the WIM confuses'),
              ('[O III]/Ha > 13 at P3 (T_e >~1.6e4 K)', 'notes/paper_review.md s4.11', 'used the paper\'s 10 R'),
              ('4.7+-1.0', 'notes/spectra_sdss_lamost.md', 'residual-flux bias + placebo')],
  conf='medium', notes='hiwham\'s strongest single route "DESI Hb + case B > 4.8" used the first-pass uninflated [O III]/Hb > 18; with the rev-2 box limit (>10.5) it is ~2.8 photon (derived here). ">~3" survives via SDSS and DESI direct.',
  checks=[check('verify_hiwham_wham_summary.json', ['oiii_ha_lower_limits', 'sdss_dec_lt17p6_direct_placebo'], 3.05, 0.01),
          check('verify_hiwham_wham_summary.json', ['oiii_ha_lower_limits', 'desi_hb_caseB_min_subset'], 4.80, 0.01)])
E('r.o3_4959', 'ratios', '[O III]5007/4959', 'DESI 3.04 (2.94-3.10), box ~2.9; SDSS free 3.15', '', '', 'measurement (consistency check vs 2.98)', '',
  ['notes/verify_desi.md c10', 'notes/verify_sdss.md row 2'], conf='high')
E('r.oii_doublet', 'ratios', '[O II]3729/3726',
  'DESI 1.30 (1.24-1.33) all and south; subsets 1.17 (QSO/z>0.5) - 1.37 (neighbour instrument). SDSS south 1.29 (68%: 1.07-1.59); SDSS all 1.13 (0.99-1.29); SDSS Halpha patch 1.12 (1.00-1.27). Low-density limit 1.46 (CHIANTI 7.1.4) / 1.467 (shock model).',
  'deblending at R~2700 is LSF/velocity sensitive; not re-done with the DESI box method', '', 'measurement (PENDING oii)', '',
  ['notes/verify_desi.md c8', 'notes/verify_sdss.md row 3g', 'data/derived/verify_desi_summary.json: ratios.*.OII3729/3726'],
  superseded=[('1.30 (1.24-1.39) -> n_e ~110 (40-160)', 'notes/spectra_desi.md', 'confirmed; wider systematic range')], conf='low',
  notes='The often-quoted "1.29-1.30 in both DESI and SDSS" holds only for SDSS south (weak, 1.07-1.59); SDSS all-strand and patch give 1.12-1.13 (higher density). wave-1 DESI noted on-ridge 1.45 vs off-ridge 1.12.',
  checks=[check('verify_desi_summary.json', ['ratios', 'all', 'OII3729/3726', 'value'], 1.30, 0.005)])
E('r.limits_shock', 'ratios', 'Other upper limits (DESI, 97.5%, relative to [O III]5007)',
  '[Ca II]7291 < 0.004; [Fe II]7155 < 0.0014; [N I]5200 < 0.043; He I 5876 < 0.005; no [Ca II]7324 under [O II]7325', '', '', 'upper limit', '',
  ['notes/verify_desi.md "New shock diagnostics"'], conf='medium',
  notes='Argue against a complete radiative shock with a cool recombination/neutral zone. Template-correlated systematics of +-0.05 could in principle bias these.')
E('r.patch', 'ratios', 'Halpha patch (115.47,+20.54) spectrum',
  'SDSS: [O II] 3.55 R, Hb 0.37 R, Ha 1.41 R; [O II]/Hb ~13 (energy); Ha/Hb 3.8 photon = 2.84 energy (unreddened case B); [S II]/Ha 1.06 (0.95-1.19); [N II]/Ha 0.78 (0.70-0.89); [O III]/Hb 1.8; plane fit 0.99/1.29. DESI patch: [N II]/Ha 0.65 (0.50-0.71), [S II]/Ha 0.83 (0.56-0.91), [O III]/Ha 0.05. NSNS: Ha +3.76 R, [S II] +2.88 R ([S II]/Ha ~0.77), [O III] +0.40 R.',
  'SDSS systematically higher than DESI by 1-2 sigma', '', 'measurement', '',
  ['notes/verify_sdss.md s6', 'notes/spectra_desi.md s7', 'notes/imaging_optical.md s3.1'], conf='medium',
  notes='Shock-like low-ionisation spectrum ([O II]/Hb ~13 with [S II]/Ha ~1 far above photoionised WIM [O II]/Hb <~4); qualitative, no models run on the patch.')

# ============================================================== G. TEMPERATURES / DENSITIES / PRESSURES
E('phys.te_o3', 'physics', 'T_e([O III]) from 4363/5007',
  '(3-6) x 10^4 K: 4.6e4 K (3.1-5.3e4, O&F approx; 3.1-5.8e4 CHIANTI) for 0.083; 5.4e4 K south', '', 'K', 'derived from a 2-3 sigma hint',
  '', ['notes/verify_desi.md c5', 'notes/verify_altmodel.md s1'], superseded=[('~2.9-3e4 K', 'notes/spectra_desi.md; notes/verify_paper.md s5', 'input 0.056 superseded')], conf='low')
E('phys.t_o2', 'physics', 'T([O II]) from 7325/3727 = 0.026',
  '1.25e4 K (1.0-1.7e4); 1.18e4 at n_e = 100', '', 'K', 'derived (input below DESI systematic floor -> effectively an upper bound)', '',
  ['notes/verify_altmodel.md s1'], conf='low', notes='Two-temperature structure T(O III)/T(O II) ~4 is the altmodel fingerprint of a cooling flow.')
E('phys.tmin', 'physics', 'Minimum T for [O III]/Hb >= 15 (fully ionised O+/O++ at [O II]/[O III]=1.37, solar O/H)',
  '>= 1.3e4 K (13,050; 12,600-13,850); >= 2.0e4 K for 47 (withdrawn)', '', 'K', 'derived', '', ['notes/verify_altmodel.md s4.1'], conf='medium',
  notes='Photoionisation equilibrium reaches only 7-12e3 K at the required ionisation; max [O III]/Hb at any U <= 14 for any tested spectrum.')
E('phys.ne_doublet', 'physics', 'n_e from [O II]3729/3726 = 1.30',
  '88 / 100 / 125 cm^-3 at T = 1.0 / 1.25 / 2.0e4 K; 47-208 over DESI subsets 1.17-1.37 (altmodel); ~110 (40-300) (verify_desi)', '+-30% from atomic data (LDL 1.42-1.50)', 'cm^-3',
  'derived (PENDING oii)', 'CHIANTI 7.1.4', ['notes/verify_altmodel.md s1', 'notes/verify_desi.md c8'], conf='low',
  notes='Implies P/k ~2.8e6 (1.4-5.7e6) K cm^-3 and an emitting path of only 40-230 AU (filling factor ~1e-5) -- in tension with the brightness-based n_e and with every pressure-balanced model (predict 1.46-1.467).')
E('phys.ne_em', 'physics', 'n_e from emission measure (2-4 R at 280 pc)',
  '0.1-0.7 cm^-3 (edge-on shell/sheet, L 25-50 pc); 0.5-3 (filament, L = 1.37 pc); EM 0.27-10.3 pc cm^-6', '', 'cm^-3', 'derived', '5-level O++ atom; T 1-5e4 K; x(O++) 0.3-1',
  ['notes/verify_paper.md s5'], superseded=[('n_e ~1-3 (for 10 R)', 'notes/paper_review.md s4.6', 'brightness scale revised')], conf='medium',
  notes='Shock model: pre-shock n0 0.012-0.042; post-shock [O II]-zone n_e 0.14-0.25 (doublet 1.467). altmodel n_rms 0.10-0.37 on a 25.6-pc path.')
E('phys.ne_sii_paper', 'physics', '[S II] density from Wu+2026 Table 2 ratios (WIM, not filament)',
  'R = 1.03-1.19 -> 225-631 cm^-3 (5-level S+, TZ10, 8000-1e4 K); 169-400 with the paper\'s eq. 5 (189-400 at 1e4 K); P4 (background) highest', '', 'cm^-3', 'derived', '',
  ['notes/verify_paper.md s3', 'data/derived/verify_paper_sii.json'], superseded=[('64-150 cm^-3', 'Wu+2026', 'inconsistent with its own Table 2'), ('190-440', 'notes/paper_review.md', 'eq. 5 biased ~25-30% low')], conf='high',
  notes='Per-region 16-84% ranges are 11-1900 cm^-3: statistically meaningless; and the lines are WIM. Low-density limit R = 1.44-1.46 (paper ~1.42).')
E('phys.pressure', 'physics', 'Pressures (P/k, K cm^-3)',
  '[O II] doublet: ~2.8e6 (1.4-5.7e6). Shock ram pressure from the width+brightness closure: 3.1e4-1.1e5. Monogem hot gas: 1.5e4 (altmodel, n_e T), 0.8-2.9e4 (verify_paper, 2 n_e T), 3.7e4-1.6e5 x f^-1/2 (shockmodel, 9.2 n0 T)', 'definitions differ by up to x10', 'K cm^-3', 'derived',
  '', ['notes/verify_altmodel.md s1', 'notes/verify_shockmodel.md s5', 'notes/verify_paper.md s5'], conf='low',
  notes='"100-400x the Monogem pressure" uses altmodel\'s 1.5e4; against shockmodel\'s Monogem range the doublet pressure is ~18-150x.')
E('phys.timescales', 'physics', 'Timescales',
  'Shock exposure t* 3.9e4-1.04e5 yr < full cooling time 0.68-1.4e5 yr; Monogem age ~1.1-1.2e5 yr (Knies+24); fossil e-fold time 21-34 yr at the doublet density; precursor t_rec ~5 Myr, t_cool ~6 Myr at n0 0.025', '', 'yr', 'model/derived', '',
  ['data/derived/verify_shockmodel_summary.json: consistent_ranges', 'notes/verify_altmodel.md s4.6', 'notes/verify_shockmodel.md s8'], conf='medium',
  checks=[check('verify_shockmodel_summary.json', ['consistent_ranges', 't_star', 0], 39374, 1), check('verify_shockmodel_summary.json', ['consistent_ranges', 't_cool', 2], 140277, 1)])

# ============================================================== H. MULTIWAVELENGTH
E('mw.fuv', 'multiwavelength', 'GALEX FUV counterpart',
  'DETECTED: template regression on NSNS [O III]: south (Dec<17.6) 92+-10 CU per NSNS-R (z=8.7), whole spine 74+-8; 68+-8 CU per R if true [O III] = 1.35 x NSNS. Direct contrast south +120+-26 CU (z=4.6), all +67+-23. Stacked peak 145 CU at NSNS 1.55 R.',
  'null p <= 0.015 (none of 64 shifts)', 'CU per R(5007)', 'measurement', 'GR6/7 skybg blocks, Tycho-2/NSNS masks',
  ['notes/verify_xrayfuv.md s1.7, s3.3'], superseded=[('"weak excess at most"; +47..+154 CU; 2 sigma UL 110-250 CU', 'notes/imaging_mw.md', 'same numbers, higher significance; detection')],
  conf='high', notes='Predictions (E(B-V)=0.03): unsuppressed shocks 196-541 CU/R; complete-shock non-resonant floor 85-277; incomplete R88 floor 57. Requires C IV suppression >=5x (xrayfuv 6.5-26x; shockmodel RT S1548 0.16-0.27 i.e. ~3.7-6x; altmodel >=3-5x). Paper\'s 10-13 R would give 12-15 CU/R, below every non-resonant floor. GALEX stops at Dec 22.')
E('mw.nuv', 'multiwavelength', 'GALEX NUV counterpart',
  'DETECTED: 49+-7.5 (south), 57+-5.5 (whole) CU per NSNS-R; joint fit 38+-8; stacked-peak ratio ~80 -> range 36-57 (up to 80); direct contrast south +96+-23, all +68+-12', '', 'CU per R(5007)', 'measurement', '',
  ['notes/verify_xrayfuv.md s1.10'], superseded=[('NUV +64..+113 CU, nature unclear, colour ~background 0.71', 'notes/imaging_mw.md', 'background colour 0.50; excess 1.0-1.9 bluer')], conf='medium',
  notes='Matches incomplete-shock NUV lines (C III] 1909, C II] 2326, 2-photon, [Ne IV]) 47-63 CU/R; dust-scattered light < 16 CU; Mg II < 1 CU/R.')
E('mw.xray', 'multiwavelength', 'Soft X-rays (eRASS1 0.2-0.6 keV; ROSAT R1+R2)',
  'No X-ray filament (on-spine bump -0.04 all / -0.38 S / +0.07 N cts s^-1 deg^-2 on 5.4-7.1). W/E = 1.26 (S 1.22, N 1.25) eRASS1; 1.24 ROSAT. Significant vs random sky placement (p=0.007) but typical vs same-radius rotations about the PSR (p=0.37; south 0.57) -> generic to the Monogem shell at that radius. Local step at the spine: south +12% at x=-0.2 deg, p 0.013-0.08 (~2 sigma), 0.04-0.25 with absorption correction; north none. Ring peak r~10.25 deg (PA 60-100): Silk ~1 deg outside, on the outer flank; Silk-sector edge softer (x1.27/deg) than the S/SW limb (x1.43-1.57).',
  '', '', 'measurement', 'trace-relative profiles with rotation/translation/shift nulls', ['notes/verify_xrayfuv.md s1-2'],
  superseded=[('W/E 1.32, 2.4 sigma; edge 1.03+-0.20 deg inside; "Silk is the outer edge of Monogem X-rays"', 'notes/imaging_mw.md', 'parallel-shift nulls cannot separate "edge" from "at the Monogem radius"'),
              ('ROSAT R12 W/E 1.14-1.81 by Dec band', 'notes/context.md s1', 'smoothed, informal errors')], conf='high',
  notes='South interior has MORE HI (+5.0e19 HI4PI), so absorption works against the observed sense.')
E('mw.halpha', 'multiwavelength', 'Halpha imaging limits on a filament',
  'No Halpha filament. SHASSA (Dec 12-15.6, 4\'): 2 sigma UL 0.77 R (Sm) / 1.02 R (CC) on the intrinsic peak; displaced/broadened (|d0|<=0.8 deg, FWHM 0.28-0.5): < 0.85 R; F03 SHASSA part 0.78 R; F03 north of Dec 15.9 is WHAM-only (1 deg): 2 sigma 1.74 R; WHAM template south 1.0 R (+-25 km/s) / 1.25 R total, north 1.25-1.30 R.',
  '', 'R', 'upper limit', 'placement nulls (400-1258)', ['notes/verify_hiwham.md s1, s3', 'data/derived/verify_hiwham_wham_summary.json'],
  superseded=[('F03 2 sigma < 0.43 R (whole spine)', 'notes/imaging_mw.md', 'F03 is WHAM-only on 76% of the spine; 6\' correction invalid'),
              ('SHASSA 3 sigma < 0.8-1.0 R', 'notes/imaging_optical.md', 'placement nulls -> 3 sigma 1.1-1.3 R'),
              ('F03 +0.8-0.9 R at Dec 17-21', 'notes/imaging_optical.md', '1-deg WHAM feature (Halpha patch beam)')], conf='high',
  checks=[check('verify_hiwham_wham_summary.json', ['oiii_ha_lower_limits', 'rows_peak_UL_R', 'SHASSA Sm 4′'], 0.77, 0.01)])
E('mw.hi', 'multiwavelength', 'HI (GALFA-HI DR2)',
  'No HI filament: 2 sigma UL 1.0-1.4e19 cm^-2 per 2-km/s slice at V_LSR -3..+5 (whole spine), 1.5-1.75e19 (south); 1.1e19 at +11..+16; ~3.0e19 for the 6-km/s helio window; on-spine N(HI, all v) 4.0e20. No wall, cavity or velocity signature at any velocity (all |z| <= 1.4).',
  '', 'cm^-2', 'upper limit', 'coherent (120) and segment (1500) nulls', ['notes/verify_hiwham.md s1, s2'],
  superseded=[('HI wall 0.86+-0.15 deg E (outside), 11/11 segments, wall/floor 1.8', 'notes/imaging_mw.md', 'REFUTED: finder artefact (random sky gives +0.84+-0.31 deg; z=-0.2)')], conf='high')
E('mw.hi.pgcc', 'multiwavelength', 'PGCC G199.09+19.95 cold clump (0.35 deg from the Halpha patch)',
  'two CNM lines at V_LSR -0.74 and +6.46 km/s (FWHM 2.7, 3.0; dT_B ~7.0, 6.9 K; T_k <= ~160 K); dN(HI) +5.4e19; only 2.7% of random positions show a narrow excess this strong; patch polygon shows -0.5/+5.0 km/s',
  'formal +-0.05-0.1 km/s; real larger', 'km/s', 'measurement', '', ['notes/verify_hiwham.md s1 row 8b, s2.5'], conf='medium',
  notes='The +6.5 coincidence with a W-strand velocity relied on +18.5 (refuted); inconclusive now.')
E('mw.hishells', 'multiwavelength', 'Catalogued HI shells',
  'GS202+17+019: weak LAB cavity (rim-centre z 1.7-2.8 at +10..+40 km/s), NOT at the revised [O III] velocity; spine crosses its interior. GS202+19-019: not confirmed in LAB.', '', '', 'measurement', '',
  ['notes/verify_hiwham.md row 8a'], superseded=[('"centred on the Silk" (catalogue coincidence)', 'notes/context.md s3', 'velocity coincidence vanished with the revised [O III] velocity')], conf='medium')
E('mw.dust', 'multiwavelength', 'Dust / IR',
  'No counterpart: IRIS 100 um 2 sigma dI < 0.19 MJy/sr; dE(B-V) < 0.004 -> N_H(dust) < 2.3e19; Planck 857/545, SFD, IRIS 60 none. Foreground E(<300 pc) 0.0040-0.0079 along spine (Edenhofer+24; 0.0051-0.0061 at P1-P4).', '', '', 'upper limit (wave 1)', '',
  ['notes/imaging_mw.md s2', 'notes/context.md s5'], conf='high')
E('mw.radio', 'multiwavelength', 'Radio continuum (not re-verified in wave 2)',
  'No total-intensity counterpart: Haslam 408 MHz -0.07+-0.62 K -> Sigma_1GHz < 1.3e-22 W m^-2 Hz^-1 sr^-1 (2 sigma, beam-diluted, alpha=0.5); WMAP K < 1.0e-22; Planck 30 GHz < 1.3e-22; DRAO 22 MHz < 2.7e-22. Stockert 1420 MHz +19.5 mK on the spine (4.0 sigma vs nulls; Sigma ~3-5e-23) sits on a W>E step (+68 mK, 4.6 sigma) also in Planck 30 GHz PI (4.2 sigma): AMBIGUOUS. NVSS/TGSS resolve out a 17\'-wide filament.',
  '', '', 'upper limit / ambiguous', '', ['notes/imaging_mw.md s1.5, s2'], conf='medium')
E('mw.rm', 'multiwavelength', 'Faraday rotation (NVSS Taylor+09; not re-verified)',
  'North half (spine Dec 20-28): E(2-6 deg) - W(2-6 deg) = +36.8+-4.6 rad m^-2; step at +2.0 deg E (-0.6..+2.0); 70% negative RMs W vs 16% E. South: +5.3+-7.3 (none). On-off (|d|<0.5 vs 1-4 deg) +19.2+-8.1 (2.4 sigma).',
  '', 'rad m^-2', 'measurement (wave 1)', 'quadratic foreground surface', ['notes/spectra_ism.md s3.4', 'data/derived/spectra_ism_rm_results.json'], conf='low',
  notes='Not claimed as a filament signature: a large-scale foreground RM gradient can do the same.',
  checks=[check('spectra_ism_rm_results.json', ['halves', 'north_dec20_28', 'E2_6_minus_W2_6_median_RM'], 36.8, 0.01),
          check('spectra_ism_rm_results.json', ['halves', 'south_dec12_20', 'E2_6_minus_W2_6_median_RM'], 5.3, 0.01)])
E('mw.gamma', 'multiwavelength', 'Gamma rays (Fermi-LAT SkyView bands)', 'No filament contrast: 0.3-1, 1-3, 3-300 GeV at -1.0, -1.5, -0.2 sigma; TeV-halo geometry not modelled', '', '', 'non-detection (wave 1)', '',
  ['notes/imaging_mw.md s2'], conf='medium')
E('mw.uvabs', 'multiwavelength', 'UV absorption (hot/interface gas)',
  'COS WD J0747+2106 (99.5 pc, 0.16 deg): Si IV < 12.2, N V < 13.0 (log N, 3 sigma, +15..+40 km/s; wave 1). IUE HD 60848 (1.9 kpc, 2.26 deg W): S II 15.01, Al III 12.86, Si IV 12.80-12.86, C IV 13.32-13.49, N V < 12.88; high ions 5-11 km/s blueward of low ions; nothing at |v| 70-200 km/s. No Silk component separable.',
  '', 'log cm^-2', 'measurement', '', ['notes/spectra_ism.md s3.3', 'notes/verify_distance.md row b2'], conf='medium')
E('mw.pulsar', 'multiwavelength', 'Pulsars near the filament',
  'PSR J0751+1807: PX 0.85+-0.04 mas (abb+23) -> 1.18 (1.12-1.24) kpc, RM +41.6+-0.2 (sbm+22), DM 30.24, 0.95 deg E (convex side), not behind the filament. PSR B0656+14: 288 +33/-27 pc, DM 13.95-13.98, RM +22.7 (ATNF 2.9.0) / +23.5+-0.4 (B/psr 2017).', '', '', 'literature', '',
  ['notes/verify_distance.md c4', 'data/derived/verify_distance_synthesis.json: pulsars', 'notes/spectra_ism.md s3.4'],
  superseded=[('J0751+1807 PX 1.6+-0.8 mas, no RM', 'notes/spectra_ism.md', 'stale VizieR B/psr 2017 snapshot')], conf='high',
  checks=[check('verify_distance_synthesis.json', ['pulsars', 'J0751_1807', 'px_mas'], 0.85, 1e-3), check('verify_distance_synthesis.json', ['pulsars', 'J0751_1807', 'RM'], 41.59, 1e-2)])

# ============================================================== I. DISTANCE
E('dist.adopted', 'distance', 'ADOPTED distance',
  '252-318 pc (68%; median 282; 95% 225-359) IF the southern spine is the Monogem limb; otherwise ~100 pc - ~1 kpc permitted. 160 pc is not excluded by any measurement; ~1.6 kpc disfavoured (470 pc long, 7.7 pc wide, z ~0.55 kpc).',
  'conditional on the Monogem association', 'pc', 'derived (conditional)', 'MC over PSR B0656+14 parallax (3.47+-0.36 mas) x cos(theta) with 15-pc centre-depth offset',
  ['notes/verify_distance.md s4', 'data/derived/verify_distance_synthesis.json: monogem_limb_mc, combined'],
  superseded=[('~160 pc (fiducial)', 'Wu+2026', 'Table 4 is kpc: ~1.6 kpc; kinematic distances meaningless'), ('250-320 pc', 'notes/context.md', 'confirmed and quantified')],
  conf='medium', notes='Near/far side of a shell up to 1 deg larger: 237-337 pc; shell radius 61 (54-68) pc. The north (Dec>17.6) bends to 12.4-16 deg and is not on the Monogem sphere. The synthesis JSON cites the wave-1 X-ray registration (W/E 1.32, 2.4 sigma) as basis; verify_xrayfuv shows that contrast is generic at that radius, so the remaining argument is the radius coincidence (12.10 deg vs R=11.0+-1.1 deg) ~1 deg outside the ring peak.',
  checks=[check('verify_distance_synthesis.json', ['monogem_limb_mc', 'd_limb_pc_percentiles_2p5_16_50_84_97p5', 2], 282, 0.5),
          check('verify_distance_synthesis.json', ['monogem_limb_mc', 'd_limb_pc_percentiles_2p5_16_50_84_97p5', 1], 252, 0.5),
          check('verify_distance_synthesis.json', ['monogem_limb_mc', 'd_limb_pc_percentiles_2p5_16_50_84_97p5', 3], 318, 0.5)])
E('dist.kinematic', 'distance', 'Kinematic distances',
  'Meaningless: dV/dd 7.4-9.5 km/s/kpc (model-dependent) ~ sigma_pec; revised V_LSR 0..+2 gives 0-0.21 kpc (verify_paper) / 0-1.4 kpc with sigma_pec 7 (verify_distance); Wu+ [O III] values give 1.15-1.80 kpc (kd output in kpc)', '', 'kpc', 'derived', '',
  ['notes/verify_paper.md a2/a4', 'notes/verify_distance.md c5'], conf='high')
E('dist.lb', 'distance', 'Local Bubble wall (O\'Neill+24)',
  'peak 161-171 pc toward P1-P4 (context); along spine median 167 (158-223), south 165, north 178; inner 139, outer 200 pc', '', 'pc', 'literature/reference', '',
  ['notes/context.md s5', 'notes/verify_distance.md c3'], conf='high', notes='Wu+\'s 160 pc sits on the wall by coincidence; no dust/HI/X-ray feature ties the Silk to it.',
  checks=[check('verify_distance_synthesis.json', ['local_bubble_oneill24_spine', 'peak_pc_median'], 167, 0.5)])
E('dist.absorption', 'distance', 'Absorption-line constraints',
  'No bracket. [O III] velocity equals the ubiquitous local gas (+6..+15 helio, present from 63 pc). HD 61717 (SB2; 272 (262-284) pc, 0.97 deg W at Dec 22.3) -17 km/s Ca II cloud: v = -16.88+-0.04, log N(Ca II) 11.59+-0.04, b 1.49+-0.30 km/s, Ca II/Na I 23 (15-49), Ti II/Ca II < 0.16 (3 sigma) -> cool, Ti-poor, compact (<1.9 deg; absent toward 14 stars), V_LSR -28 vs [O III] ~+1: unrelated. Conditional bounds: < 272 (284) pc if that cloud were part of the shell; > 99.5 pc from WD J0747+2106 if the Silk carried log N(Si IV) > 12.2.',
  '', '', 'measurement/derived', '', ['notes/verify_distance.md s1, s2, s4'],
  superseded=[('HD 61717 is A2, RV +41.8; Ca II/Na I 22+-2 (lower limit); Routly-Spitzer shocked gas; mildly suggestive d~100-270 pc', 'notes/spectra_ism.md', 'star is an SB2 (-69/+124 km/s); ratio 23 (15-49) is a value; Ti-poor, cool -> not shock-processed')], conf='high')

# ============================================================== K. MODELS
MODELS = [
    ('model.complete_shock', 'Complete steady radiative shocks (80-200 km/s)', 'RULED OUT',
     '[O III]/Hb <= 5 when cooling completes (own code 2.4-5.0; R79 2.8-4.5, CR85 2.5-3.4, DS96 1.2-1.9, Allen+08 1.2-3.2) vs DESI > 11 (2 sigma); FUV non-resonant floor 85-277 CU/R > measured 68-92; would add an Halpha filament and cooled HI (not seen). Caveat: with SDSS alone the 2 sigma floor is only > 5.3 (placebo rms).',
     ['notes/verify_shockmodel.md s1-2', 'notes/verify_xrayfuv.md s1.8'], 'high'),
    ('model.truncated_shock', 'Incomplete (column-truncated) radiative shock driven by Monogem hot gas', 'FAVOURED (viable)',
     'v_s 100-150 km/s (best 110-140); log N* 17.6-17.9; T* 1.7-2.2e4 K; n0 0.012-0.042 cm^-3 (2-8x Monogem ambient); b = B0/sqrt(n0) 7-14 uG cm^1.5 (B0 0.8-2.9 uG); ram P/k 3.1e4-1.1e5; I_pk 1.3-4.0 R (median 2.3) for R_c = 58 pc at 280 pc (G ~20); FWHM 23-38 km/s; centroid = pre-shock gas velocity; t* 3.9e4-1.04e5 yr < t_cool; FUV 52-126 CU/R with edge-on C IV scattering (S1548 0.16-0.27); proper motion 75-113 mas/yr outward. Tensions: [O II]3729/3726 (model 1.467 vs 1.30); [O II]7325 x1.8. Predictions: sharp E edge with 1/sqrt tail toward Monogem; FUV ridge ~0.1 deg outside the [O III] ridge.',
     ['notes/verify_shockmodel.md', 'data/derived/verify_shockmodel_summary.json: consistent_ranges'], 'medium'),
    ('model.precursor', 'Photoionised shock precursor as the [O III] emitter', 'RULED OUT (as pre-ioniser: viable, part of truncated-shock picture)',
     'T_e <= 1.5-2e4 K -> 4363/5007 <= 0.02-0.034, [O III]/Hb <= 10.7; dark at n0 ~0.03 (t_rec ~5 Myr, t_cool ~6 Myr); 0.2-1.6 ionising photons per incoming H explain the fitted X(H0) <= 0.3.',
     ['notes/verify_shockmodel.md s8', 'notes/verify_altmodel.md s4'], 'medium'),
    ('model.pi_wd', 'Photoionisation by WD J0736+1444 (81 kK, 150.5 pc, 2.44 deg)', 'RULED OUT',
     'U short by 27x (n_rms, coplanar) to 9e3x (doublet n), >=1e4x if the Silk is at 282 pc; PIE [O III]/Hb <= 5-9, 4363/5007 0.006, [N II]/Ha 0.2.', ['notes/verify_altmodel.md s4.2-4.3', 'notes/context.md s4'], 'high'),
    ('model.pi_monogem', 'Photoionisation by Monogem / Local Bubble EUV', 'RULED OUT',
     'Boundary flux 340-6200 cm^-2 s^-1 -> log U <= -6.1 (needs ~-3.3); even at the right U: [O III]/Hb <= 10-14, He II/[O III] 0.06 (> 0.026), 4363/5007 0.011.', ['notes/verify_altmodel.md s4.2-4.3'], 'high'),
    ('model.pi_dig', 'Photoionisation by the general ISRF / DIG', 'RULED OUT', '[O III]/Hb <= 1.3 at any U.', ['notes/verify_altmodel.md s4.2'], 'high'),
    ('model.pi_heat', 'Photoionisation plus extra heating', 'RULED OUT', 'Needs 78-350x the photo-heating for T = 4e4 K; 7325/3727 0.076 (x3); NUV 94 CU/R; plus the U deficit.', ['notes/verify_altmodel.md s4.4'], 'medium'),
    ('model.conduction', 'Thermal conduction front (static or evaporating)', 'RULED OUT',
     '4363/5007 0.11-0.14; C IV/[O III] 1.6-16 photons (FUV 540-5300 CU/R); face-on I(5007) 2e-4 - 3e-3 R per front (>=30-640 fronts needed); doublet 1.46.', ['notes/verify_altmodel.md s4.5', 'notes/verify_xrayfuv.md s1.11'], 'high'),
    ('model.tml', 'Turbulent mixing layer (Monogem hot gas / cool gas)', 'RULED OUT',
     'Needs a 377 km/s hot-gas inflow (a shock); at 10 km/s gives ~0.05 R; FUV 254 CU/R; X-ray filament 2e3-2e4x Monogem; doublet 1.46.', ['notes/verify_altmodel.md s4.5'], 'high'),
    ('model.nei_cooling', 'Non-equilibrium radiative cooling of Monogem hot gas', 'RULED OUT',
     'Matches optical ratios only if truncated before H recombines, but radiates ~70% above 5e5 K -> X-ray filament >=1e3x the Monogem shell (none in eRASS1).', ['notes/verify_altmodel.md s4.5'], 'high'),
    ('model.warm_cie', 'Warm layer held at ~6e4 K (CIE, unspecified heating)', 'DISFAVOURED',
     'Fits [O III]/Hb (~100), [O II]/[O III], 4363 and FUV (86-92 CU/R) but 7325/3727 0.08 (x3), NUV 84-93 (x2), [N II]/Ha ~2; needs 1.7e-5 erg cm^-2 s^-1 heating with no source.', ['notes/verify_altmodel.md s4.5'], 'medium'),
    ('model.fossil', 'Fossil / recombining structure', 'RULED OUT',
     '[O III] e-fold 21-34 yr at the doublet density; matching window 9-21 yr only if started >=4e4 K (shock-heated); SDSS/NSNS show no fading over ~20 yr (Delta +0.8+-1.3 R).', ['notes/verify_altmodel.md s4.6', 'notes/verify_precovery.md s4.6'], 'high'),
    ('model.isolated_snr160', 'Isolated ~30-pc SNR at 160 pc (Wu+2026)', 'DISFAVOURED',
     'No interior of its own (X-ray-bright side is Monogem); 160 pc was a unit error; energetics inputs invalid; Balmer-bound complete-shock E <= 3e50 erg at 280 pc.', ['notes/imaging_mw.md s1', 'notes/verify_paper.md s4'], 'medium'),
    ('model.monogem_limb', 'Southern Silk = Monogem (E/NE) limb at ~250-320 pc', 'FAVOURED (morphological)',
     'South spine at the Monogem radius (12.10+-0.12 deg vs R 11.0+-1.1), 1 deg outside the 0.2-0.6 keV ring peak on its outer flank; Fesen+26 assigns the 15-deg [O III] group to the Monogem limb. Not a distinct X-ray edge (rotation null p=0.37). North (Dec>17.6) not on the Monogem sphere.',
     ['notes/verify_xrayfuv.md', 'notes/verify_distance.md s4', 'notes/context.md s1'], 'medium'),
    ('model.geminga', 'Fragment of a Geminga-progenitor SNR (Ziegenbalg 2025)', 'NOT SUPPORTED',
     'Concentricity with present-day Geminga is non-physical; no back-traced birth locus near the curvature centre (Geminga 13.0 deg; B0656+14 9.05 deg).', ['notes/context.md s2'], 'low'),
    ('model.lb_wall', 'Interface on the Local Bubble wall at ~160 pc', 'DISFAVOURED (not excluded)',
     'Would need an unseen energy source; X-ray-bright side faces Monogem; no dust/HI/X-ray feature at the wall ties the Silk to it.', ['notes/context.md hypotheses', 'notes/verify_distance.md s4'], 'low'),
    ('model.north', 'Northern Silk (Dec >= 20)', 'MODEL-CONSISTENT (not unique)',
     '[O III]/Hb ~6.4, [O II]/[O III] 4-7 fit the same truncated family at a later cut (~1e4 K, log N 18.0-18.2) or a slower 70-80 km/s shock (24 tracks chi2 < 2 on these two ratios).', ['notes/verify_shockmodel.md row 15', 'data/derived/verify_shockmodel_north.json'], 'low'),
]
for mid, name, verdict, why, prov, conf in MODELS:
    E(mid, 'models', name, verdict, '', '', 'verdict', '', prov, conf=conf, notes=why)

# ============================================================== J. PRECOVERY TABLE
PRECOVERY = [
    dict(date='1950-1996', dataset='DSS1/DSS2 photographic (POSS-I/II)', result='not detected; DSS2 J A = 4.3+-3.2 -> < 14 x NSNS (3 sigma, ~<40 R)', sigma=None, src='notes/imaging_optical.md s1'),
    dict(date='1997-2000', dataset='SHASSA / WHAM Halpha', result='no Halpha filament (not an [O III] test)', sigma=None, src='notes/verify_hiwham.md'),
    dict(date='2001-11-12 .. 2002-12-06', dataset='SDSS imaging sky, north runs (Dec >= 20.5)', result='2002-01-14 runs 2887/2888 g 3.9 sigma with u, r excess -> broadband confounder, NOT attributable', sigma=None, src='notes/verify_precovery.md s2.3'),
    dict(date='2002-02-09', dataset='SDSS plate 857-52314', result='0.18+-0.05 R formally 3.5 sigma; 0.5 sigma with placebo (verify_sdss); 1.8 sigma (verify_precovery) -> not a detection', sigma=0.5, src='notes/verify_sdss.md s7'),
    dict(date='2002-10-12', dataset='SDSS plate 927-52559', result='verify_precovery 2.7 sigma; verify_sdss total-sky ON-OFF 0.50+-0.23 R (2.1 sigma formal, not placebo-tested)', sigma=2.7, src='data/derived/verify_sdss_precovery.csv; verify_precovery_table.json'),
    dict(date='2002-10-30', dataset='SDSS plate 927-52577', result='1.03+-0.21 R, 5.0 sigma (verify_precovery per-plate fit on wave-1 residual fluxes; only 3 on-ridge fibres; not placebo-validated; absent from verify_sdss table) -> EARLIEST >=3 sigma, FRAGILE', sigma=5.0, src='data/derived/verify_precovery_table.json: first_3sigma.sdss_fibres'),
    dict(date='2002-12-31', dataset='SDSS imaging run 3560 (g, PSP), Dec 20.2', result='3.77+-1.16 R, 3.3 sigma; same run r -0.4 sigma, i -0.9 sigma -> first imaging >=3 sigma (single-run sigma uncertain by ~+-0.5)', sigma=3.3, src='notes/verify_precovery.md s3'),
    dict(date='2003-01-25 / 2003-02-01', dataset='SDSS imaging, cumulative core runs', result='3.9 sigma by 2003-01-25, 4.6 sigma by 2003-02-01 (core-only cumulative per note)', sigma=4.6, src='notes/verify_precovery.md s3'),
    dict(date='2003-01-30', dataset='SDSS plate 1203-52669 (04:42-05:20 UT)', result='1.11+-0.11+-0.31 R, 3.4 sigma (verify_sdss, placebo) / 3.5 sigma (verify_precovery) -> earliest confirmed by two independent analyses', sigma=3.4, src='notes/verify_sdss.md s7; notes/verify_precovery.md s3'),
    dict(date='2003-03-10', dataset='SDSS plate 1263-52708', result='verify_precovery 3.9 sigma; verify_sdss 0.32+-0.13 R (2.5 sigma formal, not placebo-tested) -> conflicting', sigma=3.9, src='data/derived/verify_sdss_precovery.csv'),
    dict(date='2003-2012', dataset='GALEX FUV/NUV (GR6/7 AIS/MIS; per-tile dates not extracted)', result='UV filament detected (FUV z=8.7 south) -> UV-counterpart precovery (not [O III])', sigma=8.7, src='notes/verify_xrayfuv.md'),
    dict(date='2004-02-17', dataset='SDSS imaging u (run 4469; [O II])', result='first core u-band >=3 sigma (3.1 sigma); combined u 1.40+-0.24 R [O II], 4.4 sigma coherent', sigma=3.1, src='notes/verify_precovery.md s3'),
    dict(date='2004-11-08', dataset='SDSS plate 1921-53317', result='profile fit 1.88+-0.35 R (formal 5.4 sigma) with no fibre within 0.15 deg; not placebo-tested', sigma=5.4, src='notes/verify_sdss.md s7'),
    dict(date='2005-03-02', dataset='SDSS plate 2054-53431', result='1.20+-0.07 R, 4.2 sigma (placebo)', sigma=4.2, src='notes/verify_sdss.md s7'),
    dict(date='2005-03-03', dataset='SDSS plate 2074 (MJDLIST 53432; coadd MJD 53437)', result='1.45+-0.05 R, 9.7 sigma -> FIRST >=5 sigma', sigma=9.7, src='notes/verify_sdss.md s7'),
    dict(date='2008-01-30 / 2008-02-01', dataset='SDSS plates 2890 / 2915', result='5.4 / 13.9 sigma', sigma=13.9, src='notes/verify_sdss.md s7'),
    dict(date='2010-12-06 .. 2011-10-31', dataset='BOSS plates 4456, 4497, 4492, 4502, 4488, 4470, 4473, 4462, 4464', result='0.38-2.51 R, 3.2-13.9 sigma (verify_sdss); verify_precovery first BOSS 2011-01-08 (4502, 4.7 sigma)', sigma=13.9, src='notes/verify_sdss.md s7; notes/verify_precovery.md s3'),
    dict(date='2013-03-10', dataset='LAMOST LRS stellar residuals (GAC113N27M1)', result='3.1 sigma; 2013-11-06 4.4 sigma; up to 9.3 sigma (2016-01-30)', sigma=3.1, src='notes/verify_precovery.md s3'),
    dict(date='2016-03-06', dataset='SDSS plate 8954', result='0.73+-0.06 R, 5.2 sigma', sigma=5.2, src='notes/verify_sdss.md s7'),
    dict(date='2018-10-10 .. 2019-02-01', dataset='eBOSS plates 11082-11104', result='0.48-2.33 R, 4.6-16.3 sigma', sigma=16.3, src='notes/verify_sdss.md s7'),
    dict(date='2021-01 (first night 2021-01-03)', dataset='DESI DR1', result='5.7-48 sigma per month (partly circular: regression on forward model)', sigma=48.1, src='data/derived/verify_precovery_table.json: first_3sigma.desi'),
    dict(date='2022-12-12 .. 2025-04-04', dataset='NSNS DR0.2 [O III] (first [O III] camera 2022-12; full array 2024-01)', result='detected; per-pixel epochs not public; "first clear night after full configuration" ~2024-01 (self-reported)', sigma=None, src='notes/imaging_optical.md; notes/discovery.md'),
    dict(date='2023-11 (paper) / 2023-12 (G.-Y. Sun)', dataset='SWONS 50 mm f/0.95 [O III] mosaic', result='discovery (self-reported); Chinese press 2026-11-16 "early 2023" contradicts', sigma=None, src='notes/discovery.md'),
    dict(date='2024-01-14', dataset='AstroBin zh47r9 (SWONS)', result='first public image', sigma=None, src='notes/discovery.md'),
    dict(date='2024-02-01', dataset='LAMOST MRS-N plate MN074209N160100N01 (HIP 37797)', result='paper spectra (3x900 s)', sigma=None, src='notes/discovery.md'),
    dict(date='2024-06-11', dataset='Ziegenbalg, Cloudy Nights', result='first public position + [O III]/Ha/[S II] data (14 deg arc)', sigma=None, src='notes/discovery.md'),
    dict(date='2024-07-12', dataset='QHYCCD press note', result='first English use of "Silk Nebula"', sigma=None, src='notes/discovery.md'),
    dict(date='2025-08-28', dataset='Ziegenbalg 2025, RNAAS 9, 227', result='first journal paper (17 x 2.5-3 deg)', sigma=None, src='notes/discovery.md'),
    dict(date='2026-08-27', dataset='Fesen+2026 arXiv:2608.26666', result='15-deg [O III] group = Monogem E/N limb', sigma=None, src='notes/discovery.md'),
    dict(date='2026-10-09 / 2026-11-04', dataset='Wu+2026 arXiv v1 / AJ 172, 298', result='SWONS/LAMOST paper; AJ adds co-discoverer acknowledgment', sigma=None, src='notes/discovery.md'),
]
PRECOVERY_SUMMARY = {
    'earliest_ge3sigma_any': '2002-10-30 (SDSS plate 927-52577, 5.0 sigma; fragile: 3 fibres, wave-1 residual fluxes, not placebo-tested)',
    'earliest_ge3sigma_imaging': '2002-12-31 (SDSS run 3560, g, 3.3 sigma)',
    'earliest_two_analyses': '2003-01-30 (SDSS plate 1203, 3.4 / 3.5 sigma)',
    'first_ge5sigma': '2005-03-03 (SDSS plate 2074, 9.7 sigma); plate 1921 (2004-11-08, formal 5.4 sigma profile fit) untested',
    'years_before_discovery': '~21 (2002 vs 2023-11)',
    'checks': [check('verify_precovery_table.json', ['first_3sigma', 'sdss_fibres', 0], '2002-10-30', 0),
               check('verify_precovery_table.json', ['first_3sigma', 'sdss_imaging_g_psp_gauss028_OIIIunits', 0], '2002-12-31', 0),
               csv_check('verify_sdss_precovery.csv', 'plate', 1203, 'date_ut', '2003-01-30'),
               csv_check('verify_sdss_precovery.csv', 'plate', 2074, 'date_first_ut', '2005-03-03')],
}

# ============================================================== C. PAPER CLAIMS -> VERDICTS
PAPER = [
    ('Distance ~160 pc from [O III] kinematics (Table 4 "units of 100 pc")', 'WRONG', 'Table 4 is kd output in kpc: 1.79/1.85/1.20 (our 1.785/1.743/1.154); chi2 0.05 (kpc) vs 28.9 (x100 pc); MC floors 0.00-0.02 = kd 0.01-kpc grid floor. Fiducial = 1.6 kpc. Kinematic distance meaningless anyway. AJ text unchanged; no erratum (Crossref 2026-11-16).', 'verify_paper a1-a4', 'high'),
    ('Size 1 deg x 10 deg -> 3 x 28 pc; shell radius ~30 pc', 'WRONG / INCOMPLETE', 'Spine 16.94 deg arc, strand FWHM 0.28 deg; R = 30 pc is asserted, not fitted (curvature fit 21.3 deg, poor). At 282 pc: 83 pc long, 1.4 pc wide.', 'geometry; paper_review #3', 'high'),
    ('V_s = 119+-5 km/s; E = (2.0+-0.7)e51 erg; R_PDS 18+-2.8 pc; age (6.8+-1.8)e4 yr', 'ARITHMETIC OK, INPUTS INVALID', 'E 2.01e51 reproduced; R_PDS 17.09 pc; age 7.95e4 (+1/3) or 8.72e4 (printed "+1", a typo); n0 = 1 at 280 pc gives 1.2-1.7e52, at 1.6 kpc 2.9-4.2e54 (unphysical). Balmer-bound complete shock at 280 pc: n0 <= 0.002-0.034, E <= 1.4e49-2.9e50.', 'verify_paper c1-c4', 'high'),
    ('[S II] n_e ~64 (median) to ~150 cm^-3', 'WRONG', 'Own Table 2 R = 1.03-1.19 -> 225-631 (5-level) / 169-400 (eq. 5); no noise model gives a median of 64; P4 background highest; and the [S II] is WIM.', 'verify_paper b1-b2', 'high'),
    ('Red-line diagnostics ([S II]/Ha > 0.4, [N II]/Ha, BPT/SMB) show shock excitation', 'WRONG ATTRIBUTION', 'Same values in the nebula-free P4 (|Delta| <= 0.6 sigma) and in WHAM WIM; the criterion classifies the background as "shocked" too.', 'paper_review #7-8; verify_hiwham 4a', 'high'),
    ('T_e ~9000-1e4 K; on-filament 1-2e4 vs P4 7.5e3 K (line widths)', 'NOT SUPPORTED', 'Table 2 widths give 3220-7670 K with sigma_T 7300-10600 K, P(T<0) 17-38%; no on/off difference.', 'paper_review #11', 'medium'),
    ('V_t = 29-33 km/s shock turbulence', 'NOT SUPPORTED', 'Eq. 3 applied to Table 2 gives 25.8-29.9; P4 identical -> WIM turbulence.', 'paper_review #12', 'medium'),
    ('Intrinsic [O III] FWHM 43-49 km/s, broader than Ha', 'NOT MEASURED', 'At S/N 6.1 passing fits return 36-47 km/s for any true width >= 20-30; DESI <= 33 km/s; paper Fig 3a example ~50 km/s (conflict).', 'paper_review #13; verify_velocity s4', 'medium'),
    ('Near-solar uniform O/H (N2 8.82, N2S2Ha 8.70)', 'UNINFORMATIVE', 'N2 outside PP04/Marino calibration ranges; eq. 6 is linear (cubic gives 9.09-9.19); P4 gives the same -> WIM; solar 8.66 misattributed (Allende Prieto 2001 = 8.69).', 'paper_review #16', 'high'),
    ('[O III] V_helio +23.8..+26.8 (V_LSR 11-16) offset 10-16 km/s from red lines; two shell layers', 'WRONG (artefact)', 'Blue arm re-referenced to air airglow, undoing v_corr (-8.32) + 1.8 km/s; correction -10.1+-1.0; Fig 3a vs own [N I] gives +13.6+-1.2; [O III]-[N II] = +1.1+-1.0. Same as WIM/HI.', 'verify_velocity', 'high'),
    ('[O III] surface brightness 23.6-23.9 mag/arcsec^2 ~ 10-13 R', 'ARITHMETIC OK, VALUE ~5-8x TOO HIGH', 'Independent: DESI 1.62 R, SDSS 1.4-2.5 R, SDSS imaging 2.3-2.6 R, NSNS ~2 R; FUV/[O III] ratio also excludes 10-13 R under any shock model. Additive background error suspected (7-20% of sky); not provable without SWONS frames.', 'verify_paper e1-e2; verify_xrayfuv', 'high'),
    ('Red lines are DIG along the sightline; "DIG or ionised bubble cannot be excluded"', 'RIGHT for the red lines; WRONG if applied to [O III]', 'Red lines = WIM (WHAM same on/off, P4 identical). The [O III] filament is not DIG ([O III]/Hb > 11; [O III]/Ha >~3; WIM [O III]/Ha ~0.06) and not a photoionised bubble (4363, R23, He II).', 'paper_review #7; verify_altmodel', 'high'),
    ('No NVSS/TGSS radio counterpart -> consistent with an old SNR', 'UNINFORMATIVE', 'Interferometers resolve out >~16-arcmin emission. Single-dish: Haslam Sigma_1GHz < 1.3e-22 (2 sigma); Stockert +19.5 mK ambiguous.', 'paper_review #19; imaging_mw', 'high'),
    ('"Strong [O III] 4959, 5007"', 'OVERSTATED', '4959 S/N 1.8-2.3, below the paper\'s own S/N > 3 threshold (expected from 5007 at S/N 6).', 'paper_review #15', 'high'),
    ('Detected by Fesen+26 and "part of the Monogem Ring", while adopting 160 pc', 'INTERNALLY INCONSISTENT', 'Monogem/B0656+14 at ~288-300 pc; corrected kd output 1.6 kpc; not reconciled. Our favoured reading: southern Silk = Monogem limb at 252-318 pc.', 'paper_review #20; verify_distance', 'high'),
    ('"Discovered and named by the SWONS team"; "first time to its full extent"', 'PARTLY WRONG', 'Independent NSNS and NHZ detections; NSNS published the 17-deg extent (RNAAS 2025-08-28) 14 months earlier; AJ adds an acknowledgment but keeps the sentences.', 'discovery', 'high'),
    ('LAMOST MRS-N detects [O III] 5007 in ~279 fibres along the imaged filament; none in P4', 'RIGHT', 'Confirmed by DESI (~30 sigma vs empirical nulls), SDSS fibres, SDSS imaging, NSNS (P4 -0.11 R).', 'verify_desi; imaging_optical', 'high'),
    ('"Possible SNR / shock-excited" overall verdict', 'DEFENSIBLE, FOR OTHER REASONS', 'A truncated ~100-150 km/s radiative shock is favoured from [O III]/Hb, R23, 4363, He II, FUV -- most plausibly the Monogem SNR limb, not an independent 160-pc remnant.', 'verify_shockmodel; verify_altmodel', 'medium'),
    ('LAMOST Ha 4-7 km/s redder than [N II]/[S II]', 'ARTEFACT (PLAUSIBLE)', 'Reproduced by ~33% geocoronal over-subtraction (chi2 0.4/2); template mis-placement gives the wrong sign. Model only (MRS-N spectra not public).', 'verify_velocity', 'medium'),
    ('V_LSR from V_helio (20 km/s toward 18h03m50s +30d)', 'MINOR ERRORS', 'Table 4 V_LSR - V_helio varies -10.0..-15.2 within regions (true -11.9..-12.6); P2 [O III] distance follows V_helio not the printed V_LSR; red-line traditional distances not reproducible (negative-velocity fibres dropped).', 'paper_review #6; verify_paper', 'medium'),
]

# ============================================================== L. CONFLICTS
CONFLICTS = [
    dict(id='C1', topic='[O III]/Hbeta central value', positions='DESI lower limit only (>11 box, >18 Gaussian, 2 sigma); SDSS south 23.5 photon (12.7-47.6); withdrawn DESI first pass 46-47 (still in verify_desi_summary.json); altmodel Hgamma/Hb = 1.5 anomaly suggests ~15 if Hgamma is right; shockmodel fitted 40.', status='UNRESOLVED; adopt >11 (2 sigma), plausibly 15-25', owner='none (synthesis)'),
    dict(id='C2', topic='[O II]3727/[O III]5007 in the south', positions='DESI 1.11-1.16 (+-0.08, Dec<19.5) vs SDSS (BOSS-only) 1.84 (1.46-2.37; +-0.35 incl. placebo) -> 2.0 sigma (placebo-incl.) to 3.0 sigma (bootstrap-only), derived here; SDSS imaging u/g implies 0.8-1.0. Candidate causes: Dec/strand weighting with a steep N gradient, BOSS-only groups, 0.3 R placebo floor, u-band response.', status='PENDING (likely oii / northsouth)', owner='oii or northsouth'),
    dict(id='C3', topic='Absolute [O III] scale: SDSS imaging vs revised DESI', positions='imaging 2.31+-0.38 R -> obs/DESI_rev 1.49+-0.25 (2 sigma); imaging/SDSS-fibre_rev 1.1-1.3; all agree within x1.5 at 1.6-2.5 R. Possible causes: PSP temporal smoothing (amplitude 2.8-6.0 R if smoothed templates), ~25% template systematic, DESI petal-mean proxy and f_ps model (+-8%).', status='UNRESOLVED (does not affect the 4-8x verdict on 10-13 R)', owner='none'),
    dict(id='C4', topic='DESI vs LAMOST-corrected velocity at P1', positions='DESI within 0.3 deg (wave-1 fibres) +7.5+-2.2; DESI within 0.5 deg positive S/N>4 +11.9+-1.9; LAMOST corrected +16.7+-4.4 -> 1.9 sigma vs 1.0 sigma depending on the DESI estimator (derived). P2: 0.3 sigma; P3: 1.3 sigma.', status='MOSTLY RESOLVED (aperture/estimator difference; LAMOST correction is aggregate from one spectrum)', owner='none'),
    dict(id='C5', topic='[O II] doublet density', positions='3729/3726 = 1.30 -> n_e ~100 (47-208), P/k ~3e6 vs EM density 0.1-3 and every pressure-balanced model (1.46-1.467, n_e ~0.2); SDSS all/patch 1.12-1.13 (higher still); SDSS south 1.29 (1.07-1.59) weak.', status='PENDING (oii)', owner='oii'),
    dict(id='C6', topic='North-south [O III] velocity', positions='DESI north +20.9 (18.5-23.1) vs south +10.6 helio; SDSS north +16.4+-5.9 / +15.5+-3.7; Delta V_LSR ~12 km/s after the along-spine LSR correction (derived); north V_LSR ~+10.6 coincides with the northern HI component (+8.7..+11.3, Dec>~23).', status='PENDING (northsouth)', owner='northsouth'),
    dict(id='C7', topic='[O III] line width', positions='DESI <=33 (32.5; on-ridge 18 (0-27)) vs paper Fig 3a example ~50 km/s vs paper 43-49 (biased) vs model 23-38.', status='UNRESOLVED (DESI width not re-verified in wave 2)', owner='none'),
    dict(id='C8', topic='[O II]7325 and the T(O II) / model x1.8 tension', positions='DESI 0.036 (masked fit) is below DESI\'s own ~0.05 systematic floor; box 0.081 OH-contaminated; shock model 0.059-0.068.', status='UNRESOLVED; treat 7325 as an upper limit', owner='none'),
    dict(id='C9', topic='C IV suppression factor needed for FUV', positions='xrayfuv 6.5-26x (>=5 at 1 sigma, R88 floor); shockmodel RT S1548 0.16-0.27 (~3.7-6x); altmodel >=3-5x.', status='CONSISTENT at ~5x; model-dependent', owner='none'),
    dict(id='C10', topic='Monogem pressure', positions='altmodel 1.5e4; verify_paper 0.8-2.9e4 (2 n_e T); shockmodel 3.7e4-1.6e5 x f^-1/2 (9.2 n0 T). Ratio of doublet pressure to Monogem therefore 18-400x.', status='DEFINITIONAL; quote with definition', owner='none'),
    dict(id='C11', topic='Early SDSS plates 927 / 1263 / 1921', positions='927-52577 5.0 sigma (verify_precovery, wave-1 residual fluxes, 3 fibres) not in verify_sdss table; 927-52559 2.7 vs 2.1 sigma; 1263 3.9 (precovery) vs 2.5 sigma formal (verify_sdss); 1921 formal 5.4 sigma untested.', status='UNRESOLVED; earliest secure = 2003-01-30', owner='none'),
    dict(id='C12', topic='Strength of the Monogem-limb distance', positions='verify_distance_synthesis basis cites the wave-1 eRASS1 W/E 1.32 (2.4 sigma); verify_xrayfuv shows the W/E contrast is generic at that radius (rotation p=0.37) and the local step is ~2 sigma south only.', status='The 252-318 pc range rests on the radius coincidence; keep as conditional', owner='none'),
    dict(id='C13', topic='Halpha patch low-ionisation ratios', positions='SDSS [N II]/Ha 0.78, [S II]/Ha 1.06 vs DESI 0.65 / 0.83-0.85 vs NSNS [S II]/Ha 0.77: SDSS 1-2 sigma higher.', status='MINOR', owner='none'),
    dict(id='C14', topic='NSNS [S II] vs LAMOST [S II] intensity', positions='NSNS on-ridge [S II] <= 0.35 R vs LAMOST-scaled WIM [S II] 2.4-3.9 R; likely NSNS narrow-band loses the smooth DIG level.', status='UNRESOLVED (wave 1, low impact)', owner='none'),
    dict(id='C15', topic='Halpha/Hbeta and Hgamma/Hbeta in DESI', positions='Ha/Hb 3.9 (2.2-6.1) undetermined; Hgamma/Hb 1.5 vs case B 0.47.', status='UNRESOLVED (Hb at 1-2 sigma)', owner='none'),
]

# ============================================================== M. FILE-LEVEL INCONSISTENCIES
INCONSISTENCIES = [
    ('I1', 'data/derived/verify_desi_summary.json: ratios', 'Still holds the WITHDRAWN first-pass [O III]/Hb 47.48 (all), 45.93 (south), 59.8 (petal FE) and R23 128.1 / 112.5; the revision-2 limits are only in verify_desi_audit_extra.json (limits) and notes/verify_desi.md. Readers of the summary JSON get withdrawn values.'),
    ('I2', 'notes/verify_shockmodel.md s0, s4', 'Uses the withdrawn DESI inputs ("[O III]/Hb now 46 (>28)", "R23 ~76-130", "observed >=24 (DESI 97.5%)") and fits [O III]/Hb = 40 (sigma 0.2 dex). Conclusions are robust (complete shocks <= 5; truncated family spans 10-63) but the fit should be re-run with [O III]/Hb > 11.'),
    ('I3', 'notes/verify_hiwham.md s4; verify_hiwham_wham_summary.json', '"Strongest single statement DESI Hb + case B > 4.8 (> 6.3 energy)" uses the first-pass uninflated [O III]/Hb > 18 (QSO subset p2.5); with the rev-2 box limit (> 10.5) it becomes ~2.8 photon (derived here). The adopted ">~3" survives via SDSS (3.0) and DESI direct (3.4).'),
    ('I4', 'data/derived/verify_distance_synthesis.json: combined', '"statement" says "favoured 255-315 pc" while favoured_pc = [252, 318] and the note says 252-318. Basis cites wave-1 eRASS1 W/E 1.32 (2.4 sigma), superseded by verify_xrayfuv (W/E 1.26; rotation null p = 0.37).'),
    ('I5', 'data/derived/verify_precovery_table.json', 'Field "cum_core_snr" accumulates north (Dec >= 20.5) runs from 2001-11-12 (run 2711), so it is 1.30 after run 3560, whereas the note\'s core-only cumulative is 3.3 sigma after 3560, 3.9 by 2003-01-25, 4.6 by 2003-02-01. "first_3sigma" for g Frames/FramesSub (2002-01-14, run 2887) and u (2002-01-09, run 2864) are north broadband-confounded runs that the note rejects.'),
    ('I6', 'notes/verify_sdss.md', '[O III]4363/5007 is 0.038+-0.035 (energy) in s1 and 0.044+-0.037 (energy) in table row 3i.'),
    ('I7', 'DESI global velocity', 'verify_velocity uses +14.4+-1.2 (wave-1 5220 fibres, predicted-residual > 0.5 R) in the adopted mean; verify_desi prefers +11.0 (IV) / +11.3 (rev-2) / +12.2+-0.8 (positive S/N>4) / +12.7+-0.8 (S/N>5, JSON). Adopted mean shifts from +13.3 to +12.4 with the verify_desi/verify_sdss choices (derived here).'),
    ('I8', 'SDSS velocity', 'verify_sdss south +10.0+-1.8, all +12.3+-2.5 vs verify_velocity\'s independent SDSS south +7.5+-3.4, all +9.6+-2.8 (2-5 km/s lower; unexplained, within errors).'),
    ('I9', 'DESI near P1', 'verify_velocity +7.5+-2.2 (within 0.3 deg, wave-1 fibres) vs verify_desi +11.9+-1.9 (within 0.5 deg, positive S/N>4).'),
    ('I10', 'DESI [O III]/Hb and R23 rounding', 'verify_desi.md quotes "> 11" and "> 29" for box limits that are 10.53 and 28.97 in verify_desi_audit_extra.json.'),
    ('I11', 'Distance-JSON V_LSR', 'verify_distance_synthesis.json vlsr uses a spine-mean LSR correction (~-11.2): DESI +14.5 helio -> +3.3, SDSS +12.5 -> +1.3; the P1-local correction is -12.12 (-> +2.4, +0.4). Along the spine the correction runs -13.2 (Dec 12) to -9.3 (Dec 28.25).'),
    ('I12', 'PSR B0656+14 RM/DM', 'spectra_ism (VizieR B/psr 2017) RM +23.5+-0.4, DM 13.98 vs verify_distance_synthesis.json (ATNF v2.9.0) RM 22.73, DM 13.95.'),
    ('I13', 'Spine-to-PSR separation error', 'context.md 11.36+-0.16 deg vs verify_xrayfuv / verify_distance 11.36+-0.15 deg.'),
    ('I14', 'South-only curvature radius', 'geometry.md: south (Dec <= 17.5) R = 45.9 deg (22.5-71.4); discovery.md: "southern half alone (Dec <= 20) r ~ 124 deg". Different subsets.'),
    ('I15', 'SDSS imaging run counts', 'imaging_optical: 18 core runs, 17 of 18 positive; verify_precovery: 19 runs, 16 of 19 positive (different masks/criteria).'),
    ('I16', 'FUV CU-to-R conversion', 'imaging_mw "110-250 CU ~ 0.4-0.9 R at 1549 A" (W_eff 300 A assumed) vs verify_xrayfuv "0.6-1.4 R" (W_eff 357 A, E(B-V) 0.03: 1 R = 178 CU).'),
    ('I17', '[O II]7325 constraint', 'shockmodel fitted 0.034; DESI value 0.036 (0.033-0.042).'),
    ('I18', 'Task brief vs sources', '"[O II]3729/3726 = 1.29-1.30 in BOTH DESI and SDSS" holds only for SDSS south (68% 1.07-1.59); SDSS all-strand 1.13 (0.99-1.29) and patch 1.12. "LAMOST corrected" velocity means Table 4 - 10.1.'),
    ('I19', 'verify_paper inputs', 'Uses wave-1 DESI 2.4 R (Hb bound 0.12 R) and 4363/5007 0.056 (T 2.9e4 K); with 1.62 R and 0.083 the Hb bound becomes ~0.08 R and T ~4.6e4 K, which strengthens the stated x(O++)/Balmer tension.'),
    ('I20', 'verify_sdss.md precovery vs CSV', 'Plate 927 appears in verify_sdss_precovery.csv only as MJD 52559 (first 52556 = 2002-10-09), not as 52577; verify_precovery quotes 927-52559 as 2002-10-12.'),
]


def main():
    n_checks = sum(len(e['json_check']) for e in ENTRIES) + len(PRECOVERY_SUMMARY['checks'])
    n_fail = sum(1 for e in ENTRIES for c in e['json_check'] if not c['pass']) + \
        sum(1 for c in PRECOVERY_SUMMARY['checks'] if not c['pass'])
    w3 = {t: os.path.exists(os.path.join(ROOT, 'notes', f'w3_{t}.md')) for t in ('oii', 'northsouth')}
    out = {
        'tag': 'canon', 'wave': 3, 'built': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
        'project': 'Silk Nebula ([O III] filament, Gemini)', 'paper': 'Wu et al. 2026, arXiv:2610.13329, AJ 172, 298',
        'precedence': 'wave-2 (notes/verify_*.md) supersedes wave-1 where they disagree; within wave 2, revision-2 audits supersede first passes',
        'conventions': {
            'velocity': 'heliocentric unless V_LSR; V_LSR uses the standard solar motion (20 km/s to 18h03m50.29s +30d00m16.8s J2000); V_LSR - V_helio = -13.2 (Dec 12) .. -12.1 (P1) .. -9.3 (Dec 28.25)',
            'ratios': 'energy ratios unless "photon"; SDSS/WHAM natively photon (R)', 'intensity': 'Rayleigh (photon)',
            'status_values': 'measurement | lower limit | upper limit | hint | systematics-limited | derived | model | verdict | literature | adopted',
        },
        'wave3_teammates': {k: ('note present: incorporated' if v else 'note absent at build time: PENDING') for k, v in w3.items()},
        'json_checks': {'n': n_checks, 'n_fail': n_fail},
        'entries': ENTRIES,
        'paper_claims': [dict(claim=c, verdict=v, evidence=ev, source=s, confidence=cf) for c, v, ev, s, cf in PAPER],
        'precovery_table': PRECOVERY, 'precovery_summary': PRECOVERY_SUMMARY,
        'unresolved_conflicts': CONFLICTS,
        'file_inconsistencies': [dict(id=i, where=w, issue=s) for i, w, s in INCONSISTENCIES],
    }
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False, default=float)
    print(f'wrote {OUT}: {len(ENTRIES)} entries, {len(PAPER)} paper claims, {len(PRECOVERY)} precovery rows, '
          f'{len(CONFLICTS)} conflicts, {len(INCONSISTENCIES)} inconsistencies; json checks {n_checks - n_fail}/{n_checks} pass; w3 notes {w3}')
    for e in ENTRIES:
        for c in e['json_check']:
            if not c['pass']:
                print('CHECK FAIL', e['id'], c)
    for c in PRECOVERY_SUMMARY['checks']:
        if not c['pass']:
            print('CHECK FAIL precovery', c)


if __name__ == '__main__':
    main()
