# w3_canon: canonical values for the Silk Nebula final report

Wave 3, tag `canon`, written 2026-11-16. This note reconciles every wave-1 note (`discovery`, `paper_review`, `geometry`, `context`, `imaging_optical`, `imaging_mw`, `spectra_desi`, `spectra_sdss_lamost`, `spectra_ism`) with the wave-2 adversarial checks (`verify_{paper,velocity,desi,sdss,hiwham,xrayfuv,precovery,distance,shockmodel,altmodel}`). No new measurements were made. The only new numbers are simple derivations from quoted inputs, marked **(derived here)**: LSR corrections along the spine, physical sizes, and tension significances.

**Precedence rules**
1. Wave 2 supersedes wave 1 where they disagree.
2. Within wave 2, revision-2 audits supersede first passes. This matters for `verify_desi`, whose first-pass [O III]/Hβ ≈ 47 and R23 ≈ 128 are withdrawn.
3. Where two wave-2 agents disagree, both values are given and the conflict is listed in §12.

**Machine-readable companion:** `data/derived/w3_canon_values.json`.
- 83 entries, 20 paper claims, 30 precovery rows, 15 conflicts and 20 file inconsistencies.
- Each entry has: adopted value, uncertainty, status, method, file provenance, superseded values and why, and confidence.
- Built by `scripts/w3_canon_build.py`. The script re-reads the source JSONs and verifies 70 quoted numbers against them: **70/70 pass**.

**Figure:** `figures/w3_canon_ladder.png` shows [O III] brightness and velocity by dataset, with the adopted bands.

**Wave-3 teammates:** `oii` and `northsouth` had written intermediate data products (`data/derived/w3_oii_*`, `w3_northsouth_*`) but **no notes** when this was written. Their conflicts (§12 C2, C5, C6) are marked **PENDING**. No numbers from their unfinished intermediates are used.

**Conventions**
- Velocities are heliocentric (V_helio) unless marked V_LSR.
- V_LSR uses the standard solar motion (20 km/s toward 18h03m50.29s +30°00′16.8″, the paper's convention).
- V_LSR − V_helio runs along the spine as follows **(derived here, astropy)**:
  - −13.20 km/s at Dec 12;
  - −12.12 at P1;
  - −11.33 at Dec 19.5;
  - −10.28 at Dec 24;
  - −9.27 at Dec 28.25.
  - The spine-mean values are −12.3 (south, Dec < 19.5) and −10.3 (north).
- Line ratios are energy ratios unless marked "photon". SDSS and WHAM values are natively photon (R). Intensities are in rayleighs.

**Status vocabulary:** measurement / lower limit / upper limit / hint (2–3σ) / systematics-limited / derived / model / adopted / literature. Confidence is high / medium / low.

---

## 0. Headline canonical values

| Quantity | Adopted | Status, confidence | Main source |
|---|---|---|---|
| Spine | (116.54, +12.00) → (114.36, +28.25); l 208.15→191.37, b +17.50→+21.85 | measurement, high | `silk_trace.json` v1.0 |
| Length / width | 16.94° arc (16.38° chord); E-strand FWHM ≈ 0.28° (HWHM 0.138°, 16–84 %: 0.076–0.191) | measurement, high / medium | `geometry.md` §4.2 |
| [O III]5007 brightness | south ridge ≈ **2 R (1.4–2.5 R)**, DESI reference 1.62 R; north 0.4–0.9 R. Paper's 10–13 R is 4–8× too high | adopted, high | §4 |
| [O III] velocity | **V_helio +13.3 ± 1.0 km/s; V_LSR +1.2 ± 1.0** (±1 km/s estimator systematic: +12.4 ± 0.9 with the alternative estimators, derived here) | adopted, high | `verify_velocity_summary.json` |
| North vs south | DESI +20.9 (18.5–23.1) vs +10.6 (10.2–12.0) km/s helio; ΔV_LSR ≈ 12 km/s (derived) | measurement, medium; **PENDING northsouth** | §5 |
| [O III]/Hβ | **> 11 (2σ, DESI box)**; > 18 (DESI Gaussian); SDSS south 23.5 photon (12.7–47.6). Central value undetermined | lower limit, medium | §6 |
| [O II]/[O III] (energy) | south 1.1–1.8 (DESI 1.11–1.16 vs SDSS 1.84: conflict); north 4–7 | measurement, medium | §6 |
| He II/[O III] | < 0.026 (97.5 %) | upper limit, high | `verify_desi.md` |
| [O III]4363/5007 | 0.06–0.09 at 2–3σ; T_e(O III) ≈ (3–6)×10⁴ K; Hg I 4358 excluded | hint, low | §6–7 |
| [O II]3729/3726 | 1.30 → n_e ≈ 100 (47–208) cm⁻³, P/k ≈ 3×10⁶ | measurement, low; **PENDING oii** | §6–7 |
| [O III]/Hα (south) | ≳ 3 (2σ, photon; ≳ 4 energy) | lower limit, medium | `verify_hiwham.md` §4 |
| GALEX FUV / NUV | **detected**: FUV 92 ± 10 (south) / 74 ± 8 (all) CU per NSNS-R; NUV 36–57 CU per R | measurement, high / medium | `verify_xrayfuv.md` |
| X-ray, Hα, HI, dust filament | none (Hα 2σ < 0.77–1.0 R south; HI < 1–1.75×10¹⁹ cm⁻²) | upper limits, high | §8 |
| Distance | **252–318 pc (68 %; median 282) if the south spine is the Monogem limb**; otherwise ~0.1–1 kpc; 160 pc not excluded; 1.6 kpc disfavoured | derived (conditional), medium | `verify_distance.md` |
| Excitation | column-truncated radiative shock, v_s 100–150 km/s, n₀ 0.012–0.042 cm⁻³ | model, medium | `verify_shockmodel.md` |
| Earliest ≥ 3σ precovery | 2002-10-30 (plate 927, fragile); 2002-12-31 (imaging); **2003-01-30 (plate 1203, two analyses)**; first ≥ 5σ 2005-03-03 | measurement, medium | §10 |

---

## 1. Geometry

| Quantity | Adopted value | Uncertainty / caveat | Provenance | Superseded |
|---|---|---|---|---|
| Spine (E strand) | 66 points, Dec 12.00–28.25 in 0.25° steps; (116.5405, +12.00) l,b (208.147, +17.500) → (114.3622, +28.25) (191.367, +21.849) | per-point σ median 0.022°, max 0.20° at the junction (Dec 18.5–18.75); north of Dec 26.25 only NSNS-wide | `silk_trace.json` (checked); `geometry.md` §4.2 | `0.1-prelim` trace (overwritten) |
| Extensions not in the spine | Fainter emission to Dec ~11.3 and ~28.9 (NSNS wide, S/N ~7). NSNS O200.0b ridge at RA 117.55–117.65, Dec 11.7–12.9 (8–12σ prominence) | — | `geometry.md` §4.2 | — |
| Length | 16.94° arc, 16.38° chord. Bright south part (Dec 12.0–21.5) 9.8°; faint north arc 7.2°. The paper's own Fig 1 spans Dec 12.82–26.29 and shows the spine over 13.8° | ±~5 % | `geometry.md`, `discovery.md` §5.3 | "1°×10°" (Wu+26) = bright south part; 1.5°×15° (QHYCCD); 14° (NSNS 2024) |
| Physical length (derived here) | 47.3 pc at 160 pc; **74.5 / 83.4 / 94.0 pc at 252 / 282 / 318 pc**; 473 pc at 1.6 kpc | scales with d | — | paper "28 pc" (10° at 160 pc) |
| Width | E-strand HWHM 0.138° (0.076–0.191) → FWHM ≈ 0.28° = 1.21 / 1.36 / 1.53 pc at 252 / 282 / 318 pc (derived). Envelope (E+W+W2 ± HWHM) half-width: median 0.223°, maximum 0.83° | ±30 % (JPEG stretch). SDSS-imaging shift-scan side lobes suggest broader or structured emission | `silk_trace.json` (checked) | "~1° width" = envelope |
| Height above midplane (derived) | z ≈ 105 / 115 / 127 pc at 252 / 282 / 318 pc; 553 pc at 1.6 kpc | z☉ = 20.8 pc | — | — |
| W strand | Splits from E near Dec 15. Bright to Dec ~19.75, curving west from RA 116.49 to 115.84. 4.84° long, 0.08–1.08° from E. Ends at the Hα patch; P2 lies on it. **Carries most [O III] at Dec 15–20**: SDSS ON−OFF W 1.93 ± 0.11 R vs E 0.73 ± 0.08 R; DESI W@18 node 3.88 ± 0.20 R | identification interpretive at Dec 15–19 | `geometry.md` §4.2; `verify_sdss_results_base_r0.15.json` | — |
| W2 strand; junction | W2 at RA 115.1–115.35, Dec 21.5–24.0. Junction at Dec 17.5–17.75: the E ridge steps 0.22° east (116.61 → 116.83) and stays faint (S/N 6–12) to Dec ~19 | — | `geometry.md` §4.2 | — |
| Hα/[S II] patch | Half-maximum centroid (115.474, +20.536) = (l, b) (199.46, +20.00). 1.45° × 0.31° at PA 163°, 0.31 deg². 30 %-of-peak contour 2.24° × 1.28°. Centroid 1.28° west of the E spine | threshold dependent; [S II] IoU 0.65 | `silk_trace.json` (checked) | 0.5°×1° (paper); NSNS O200.0a 3.06°×2.0° |
| Curvature | Small circle, full E spine: centre (94.50, +16.76) = (l, b) (194.26, +0.44); **R = 21.33°** (bootstrap 20.48–22.35). χ²/dof 3744/63, rms 0.150°, so a **poor model** (S-shaped arc). Sagitta 1.34° at Dec 19.75, bowing east (equivalent R ≈ 25.8°). North (Dec ≥ 17.75) R = 21.8° (19.6–25.5); south (Dec ≤ 17.5) R = 45.9° (22.5–71.4), nearly straight | bootstrap understates model error; the centre is set by the north | `silk_trace.json` (checked); `geometry.md` §4.3 | paper "R ≈ 30 pc" (asserted, not fitted) |
| Monogem centres | South spine (Dec < 17.6, 23 points): **12.10 ± 0.12° (rms) from the Knies+24 two-SNR centre** (104.15, +15.30; R = 11.00 ± 1.1°) and **11.36 ± 0.15° from PSR B0656+14**. North bends to 12.39–16.03°. Single-SNR centre (107.08, +13.50; R = 11.97°): 9.08–16.23°, rim crossing at (116.71, +21.18), 9.25° of arc inside and 7.68° outside. The eRASS1 0.2–0.6 keV ring peak (PA 60–100) lies at r ≈ 10.25° from the PSR, so the Silk is about 1° outside it, on the outer flank | the south section is short and straight, so its low rms does not single out a centre | `verify_distance_synthesis.json` (checked); `verify_xrayfuv.md` §1; `context.md` §1 | — |
| Geminga / pulsar birth loci | 17.56 ± 0.22° (rms) from present-day Geminga, nearly concentric but **not physical**. The curvature centre is 3.9° from Geminga, 12.6° from the Monogem single-SNR centre, 10.4° from B0656+14 and 0.8° from the Origem Loop. No back-traced birth locus coincides: B0656+14 is 9.05° away, Geminga 13.0° | — | `context.md` §2 | — |
| P1–P4 | P1 (116.65, +16.60) on E (−0.005°). P2 (116.21, +17.74) on W. P3 (116.45, +14.72) 0.115° west of the merged ridge. P4 (115.20, +16.18) 1.05–1.35° west, a valid background (NSNS [O III] −0.11 R). The drawn Fig 1b outlines are offset from Table 1 by 0.48° (P1) and 0.42° (P4). The paper's "centre" (116.68, +16.42) is 3.8° south of the spine midpoint | — | `geometry.md` §4.1 | — |

---

## 2. Discovery timeline and credit

**Credit verdict (high confidence).** At least three independent detections:
- **SWONS:** 2023-11 per the paper, 2023-12 per G.-Y. Sun.
- **S. Ziegenbalg / NSNS:** ~2024-01, self-reported.
- **Schaeffer–Fesen–Drechsler (NHZ):** undated; published as the Monogem north-east limb in Fesen+2026.
- The Koposov/DESI claim is hearsay. DESI DR1 public spectra do show the filament at ~30σ (`verify_desi`), so the claim is technically plausible.
- **First journal paper:** Ziegenbalg 2025 RNAAS (2025-08-28), 14 months before Wu+ on arXiv.
- The AJ version (2026-11-04) adds an acknowledgment of the co-discoverers. It keeps "discovered and named by the SWONS team" and "detected for the first time to its full extent".

| Date | Event | Reliability (`discovery.md`) |
|---|---|---|
| 2023-11 / 2023-12 / "late 2023" | SWONS detection (paper / G.-Y. Sun blog / QHYCCD). The Chinese press of 2026-11-16 says "early 2023", which contradicts all three | M |
| 2022-12 → 2024-01 | NSNS [O III] cameras: first 2022-12, full array 2024-01; the arc was seen on "the first clear night after full configuration" | H / M |
| **2024-01-14** | AstroBin zh47r9: **first public image** (1°×10°, no coordinates) | H |
| 2024-02-01 | LAMOST MRS-N plate MN074209N160100N01 (HIP 37797), 3×900 s | H |
| **2024-06-11** | Ziegenbalg on Cloudy Nights: 14° arc, first public position and multi-line data | H |
| 2024-06-12 / 07-13 | T. Schaeffer (NHZ): "part of the Monogem SNR" | H |
| **2024-07-12** | QHYCCD note: first English use of the name "Silk Nebula" (1.5°×15°, 24.3 mag/arcsec²) | H |
| 2025-05-22 / 06-11 | NSNS DR0.2; O200.0 at 17°, northern part outside Monogem, Geminga hypothesis | H |
| **2025-08-28** | Ziegenbalg 2025, RNAAS 9, 227 (17° × 2.5–3°): **first journal paper** | H |
| 2026-08-27 | Fesen+2026 arXiv:2608.26666: the 15° [O III] group is the Monogem E/N limb | H |
| 2026-10-09 / 2026-11-04 | Wu+2026 arXiv v1 / AJ 172, 298 | H |

The earlier records are statistical precoveries in archival data (§10). The physical [O III] emission is present in SDSS from 2002, about 21 years before discovery.

---

## 3. Paper claims → verdicts (Wu et al. 2026)

| # | Claim | Verdict | Evidence (wave-2 values) | Source |
|---|---|---|---|---|
| 1 | Distance ~160 pc from [O III] kinematics ("units of 100 pc") | **WRONG** | Table 4 is `kd` output in **kpc**: printed 1.79/1.85/1.20 vs our 1.785/1.743/1.154. χ² is 0.05 if kpc, 28.9 (p = 2×10⁻⁶) if ×100 pc. The MC floors 0.00–0.02 equal `kd`'s 0.01-kpc grid floor. Fiducial = **1.6 kpc**. The AJ text is unchanged and Crossref lists no erratum (2026-11-16) | `verify_paper` a1–a4 |
| 2 | Kinematic distance usable here | **WRONG** | dV/dd is 7.4–9.5 km/s/kpc, comparable to σ_pec. The revised V_LSR gives 0–0.21 kpc (or 0–1.4 kpc with σ_pec 7) | `verify_paper` a4; `verify_distance` c5 |
| 3 | Size 1°×10° → 3 × 28 pc; shell radius ~30 pc | **WRONG / incomplete** | 16.94° arc and 0.28° FWHM; R = 30 pc asserted. At 282 pc: 83 pc × 1.4 pc | `geometry` |
| 4 | V_s 119 ± 5 km/s; E = (2.0 ± 0.7)×10⁵¹ erg; R_PDS 18 ± 2.8 pc; age (6.8 ± 1.8)×10⁴ yr | **Arithmetic OK, inputs invalid** | E = 2.01×10⁵¹ reproduced; R_PDS = 17.09 pc; age 7.95×10⁴ yr with "+1/3", or 8.72×10⁴ with the printed "+1" (a typo). With n₀ = 1: 1.2–1.7×10⁵² erg at 280 pc and 2.9–4.2×10⁵⁴ at 1.6 kpc. Balmer-bound complete shock at 280 pc: n₀ ≤ 0.002–0.034, E ≤ 1.4×10⁴⁹–2.9×10⁵⁰ erg | `verify_paper` c1–c4 |
| 5 | [S II] n_e ≈ 64 (median) – 150 cm⁻³ | **WRONG** | The paper's own Table 2 (R = 1.03–1.19) gives **225–631 cm⁻³** (5-level S⁺, TZ10, 8000–10⁴ K) or 169–400 (its eq. 5). No noise model gives a median of 64. Background P4 is highest. The [S II] is WIM in any case | `verify_paper` b1–b2 |
| 6 | Red-line diagnostics ([S II]/Hα > 0.4, [N II]/Hα, BPT/SMB) ⇒ shock excitation | **Wrong attribution** | Identical in nebula-free P4 (\|Δ\| ≤ 0.6σ) and in the WHAM WIM. The criterion classifies the background as "shocked" too | `paper_review` #7–8; `verify_hiwham` 4a |
| 7 | T_e ≈ 9000–10⁴ K; 1–2×10⁴ K on the filament vs 7.5×10³ K in P4 | **Not supported** | 3220–7670 K, σ_T 7300–10 600 K, P(T < 0) = 17–38 % | `paper_review` #11 |
| 8 | V_t = 29–33 km/s shock turbulence | **Not supported** | Eq. 3 on Table 2 gives 25.8–29.9; P4 is the same (WIM) | `paper_review` #12 |
| 9 | Intrinsic [O III] FWHM 43–49 km/s, broader than Hα | **Not measured** | At S/N 6.1, fits passing the paper's cut return 36–47 km/s for any true width. DESI gives ≤ 33 km/s, the paper's Fig 3a example ~50 km/s (§12 C7) | `paper_review` #13; `verify_velocity` |
| 10 | Near-solar, uniform O/H (N2 8.82, N2S2Hα 8.70) | **Uninformative** | Outside the calibration ranges; eq. 6 is linear, not cubic; P4 gives the same; the solar value is misattributed | `paper_review` #16 |
| 11 | [O III] V_helio +23.8…+26.8 (V_LSR 11–16), 10–16 km/s off the red lines; "two shell layers" | **WRONG (artefact)** | The blue arm was re-referenced to air airglow, undoing v_corr (−8.32) plus a 1.8 km/s offset. Correction **−10.1 ± 1.0**. Fig 3a against the paper's own [N I] gives **+13.6 ± 1.2**; [O III] − [N II] = +1.1 ± 1.0 | `verify_velocity` |
| 12 | Brightness 23.6–23.9 mag/arcsec² ≈ 10–13 R | **Arithmetic OK; value ~5–8× too high** | 12.75 / 9.94 / 10.12 R reproduced. Independent values: DESI 1.62 R, SDSS 1.4–2.5 R, SDSS imaging 2.3–2.6 R, NSNS ~2 R. Multiplicative systematics reach only ~5.5–7.5 R. An additive background error of 7–20 % of the sky is required. The FUV/[O III] ratio also excludes 10–13 R | `verify_paper` e; `verify_xrayfuv` |
| 13 | Red lines are DIG along the sightline; "DIG or ionised bubble cannot be excluded" | **Right for the red lines; wrong for the [O III]** | WHAM is the same on and off the filament. The [O III] is not DIG ([O III]/Hβ > 11, [O III]/Hα ≳ 3, vs WIM ~0.06) and not a photoionised bubble (4363, R23, He II, FUV) | `paper_review`; `verify_altmodel` |
| 14 | No NVSS/TGSS radio ⇒ consistent with an old SNR | **Uninformative** | Interferometers resolve out >~16′. Single-dish Haslam gives Σ₁GHz < 1.3×10⁻²² (2σ); the Stockert +19.5 mK signal is ambiguous | `paper_review` #19; `imaging_mw` |
| 15 | "Strong [O III] λλ4959, 5007" | **Overstated** | 4959 has S/N 1.8–2.3, below the paper's S/N > 3 threshold | `paper_review` #15 |
| 16 | Part of the Monogem Ring (via Fesen) while at 160 pc | **Internally inconsistent** | Monogem is at ~288–300 pc; the corrected `kd` value is 1.6 kpc | `paper_review` #20 |
| 17 | "Discovered and named by SWONS"; "first time to its full extent" | **Partly wrong** | Independent NSNS and NHZ discoveries; NSNS published 17° 14 months earlier | `discovery` |
| 18 | [O III] detected in ~279 LAMOST fibres along the filament; none in P4 | **RIGHT** | DESI ~30σ against empirical nulls; SDSS fibres; SDSS imaging; NSNS | `verify_desi`; `imaging_optical` |
| 19 | Overall "possible SNR / shock-excited" | **Defensible, for other reasons** | A truncated 100–150 km/s shock is favoured, most plausibly the Monogem limb | `verify_shockmodel`; `verify_altmodel` |
| 20 | LAMOST Hα (+16.3…+19.7) vs [N II]/[S II] | **Artefact (plausible)** | ~33 % geocoronal over-subtraction reproduces it (χ² 0.4 for 2 dof); a mis-placed template gives the wrong sign. Model only | `verify_velocity` |
| — | Minor | **Errors** | V_LSR − V_helio varies −10.0…−15.2 within a region; the P2 distance follows V_helio; red-line "traditional" distances are not reproducible; "expansion speed ~140 km/s, Sect. 3.4" does not exist; the shock threshold is quoted as 0.4 / 0.45 / 0.45–0.5; Table 2 P3 [S II] ratio 1.19 vs intensities (1.10) | `paper_review` §4.12; `verify_paper` |

---

## 4. [O III]5007 surface brightness

| Dataset (epoch) | Value | Uncertainty / systematics | Provenance | Superseded (why) |
|---|---|---|---|---|
| **DESI DR1** (2021–22), fibre scale | **South (Dec < 19.5) median 1.62 R** (variants 1.60–1.71; box global amplitude 1.62). **North 0.56 R** (node median) / **0.71 R** (box, Dec ≥ 19.5). E nodes 1.4–2.3 R at Dec 14–18; W 1.15–1.62 R at Dec 15–17; W@18 3.88 ± 0.20 R. Best model-free bin 1.81 ± 0.05 R, so ≥ 1.7 R (2σ) | ±8 % (GFA acceptance, chromatic model). Differential only (petal-mean proxy) | `verify_desi.md` "Bottom line" 2; `verify_desi_summary.json` brightness; `verify_desi_audit_extra.json` (checked) | 2.4 R S / 0.86 R N / ≥ 2.8 R (wave 1). The Legacy-Surveys `MEAN_PSF_TO_FIBER_SPECFLUX` was used as fibre acceptance, ×1/0.736 too high |
| **SDSS/BOSS/SDSS-V sky fibres** (2001–2021), total sky | **South ON−OFF 1.39 ± 0.07 ± 0.07 R** (r_on 0.15°; 1.61 at 0.10°, 1.06 at 0.25°). W strand (Dec 15–20) 1.93 ± 0.11 R. Per-plate peaks **1.6–2.5 R south, 0.4–0.8 R north**. Absolute P1/P2/P3 1.49 ± 0.06 / 1.78 ± 0.05 / 1.87 ± 0.09 R (±0.2 R zero point). Single fibres up to 3.0–3.65 R | f_eff 0.84 ± 0.08 (3″) / 0.60 ± 0.06 (2″) | `verify_sdss.md` §1, §3; `verify_sdss_results_base_r0.15.json` (checked) | 0.79 R south, W 1.19 R (wave 1): the pipeline residual `flux` was used and the sky model absorbs 36–53 %. f 0.93/0.66 → 0.84/0.60 |
| **SDSS imaging sky levels** (2002-12 to 2007-02), 17′ scale | g PSP **2.31 ± 0.38 R** (0.28° ridge template, 19 core runs). Frames/FramesSub **2.4–2.6 R** (preferred). With all confounder regressors ≥ 1.65 ± 0.31 R. 10–13 R excluded at ≥ 15σ on 17′ scales | ~25 % template systematic; PSP temporal smoothing (smoothed templates give 2.8–6.0 R) | `verify_precovery.md` §2.1, §2.4–2.5 | A = 1.35 ± 0.19 × NSNS → 2.8–2.9 R medians and 4.0–4.9 R "1′ peaks" (wave 1; the 1′ peaks were NSNS morphology scaled) |
| **NSNS DR0.2** (2022-12 to 2025-04), not re-verified | On-ridge P1 2.04, P2 2.06, P3 2.18 R; 1′ peaks 2.95 / 3.64 / 3.37 R; south_bright 1.60; mid (Dec 17–19) 1.21; north arc 0.58–0.62 R. Ziegenbalg's own estimates: K1a 2.04, K1a′ 1.53 R | stat 0.01–0.09 R; E/W background 0.05–0.15 R | `imaging_optical.md` §3.1 | — |
| LAMOST LRS stellar residuals (wave 1) | ON−OFF 0.88 ± 0.13 R | residual method; not re-verified | `spectra_sdss_lamost.md` | — |
| LAMOST MRS-N relative × sky-line scale (wave 1) | 2.8–3.2 R (APO [N I]/OH ratio) or 14–16 R (UVES ratio) | indeterminate (×5 sky-ratio spread) | `spectra_sdss_lamost.md` §4.1 | — |
| QHYCCD note | 24.3 mag/arcsec² → 3.3–6.7 R | bandwidth assumption | `discovery.md` §5.4 | — |
| DSS2 J (1987–94) | < 14 × NSNS (3σ), ≈ < 40 R | plate structure | `imaging_optical.md` | — |
| **Wu+26 SWONS** | 12.75 / 9.94 / 10.12 R (P1/P2/P3) | 8.3–14.9 (50–70 Å), 9.1–14.0 (AB/Vega), 7.4–9.5 (4959 leak) | `verify_paper.md` §6 | — |

**Adopted.**
- **Southern ridge peak ≈ 2 R, range 1.4–2.5 R** on 0.1–0.3° scales. DESI's flux-calibrated 1.62 R is the reference value.
- **North 0.4–0.9 R.**
- The paper's 10–13 R is **4–8× too high**.

**No variability.**
- At fixed position in SDSS imaging, Δ = +0.8 ± 1.3 R over 2003.97–2007.14.
- The DESI-normalised trend is −0.16 ± 0.13 per year over 2002.9–2007.1.
- Any change over 2002–2022 is less than about a factor of 2.

**Unresolved.** The imaging / revised-DESI scale ratio is 1.49 ± 0.25 (§12 C3).

---

## 5. Velocities

### 5.1 Per dataset ([O III] unless stated)

| Dataset | V_helio (km/s) | Notes | Provenance |
|---|---|---|---|
| Wu+26 Table 4 (LAMOST MRS-N) | P1 +26.8 ± 4.3, P2 +26.8 ± 4.1, P3 +23.8 ± 3.3 (V_LSR 15.5/14.5/11.2) | **Refuted as heliocentric** (topocentric, air, airglow-referenced) | paper; `verify_velocity` |
| LAMOST corrected (Table 4 − 10.1 ± 1.0) | +16.7 / +16.7 / +13.7 (±3–4). V_LSR +4.6 / +4.6 / +1.1 with local corrections (derived; the velocity note quotes +4.5/+4.5/+1.5 using a single −12.19) | aggregate correction from one spectrum | `verify_velocity_summary.json` (checked) |
| LAMOST Fig 3a vs its own [N I] airglow | **+13.6 ± 1.2** (V_topo +21.9, v_corr −8.32) | calibration-free | `verify_velocity_summary.json` (checked) |
| DESI, on-ridge stack (`verify_velocity`) | **+14.4 ± 1.2** (±1.9 including the ±1.5 zero point); south (Dec < 17.6) +13.1 ± 1.2; negative-residual fibres +12.6 ± 1.1 | wave-1 5220 fibres, predicted residual > 0.5 R | `verify_velocity_desi_sdss.json` (checked) |
| DESI (`verify_desi`) | IV regression +11.0 (10.3–14.3). Positive fibres S/N > 4: +12.2 ± 0.8 (n = 346). **Rev-2 joint regression on both lines +11.3 (10.2–13.9)**. Dec < 19.5: **+10.6 (10.2–12.0)**. Dec ≥ 19.5: **+20.9 (18.5–23.1)** | barycentric frame confirmed; tiles 41485/41487 (−19 km/s) unused | `verify_desi_audit_vel.json`, `verify_desi_velocity.json` (checked) |
| SDSS (`verify_sdss`) | **South +10.0 ± 1.8 ± 0.8** (SDSS-I/II +8.7 ± 2.6; BOSS +11.2 ± 2.5); all +12.3 ± 2.5; north +16.4 ± 5.9; W +9.7 ± 2.2; E +11.3 ± 4.1 | [O I]5577 zero point | `verify_sdss_results_base_r0.15.json` (checked) |
| SDSS (`verify_velocity`, independent) | all +9.6 ± 2.8; south +7.5 ± 3.4; north +15.5 ± 3.7; variants +6.8…+10.5 | 2–5 km/s below `verify_sdss` (unexplained) | `verify_velocity_desi_sdss.json` |
| WHAM Hα (WIM) | V_LSR south on/off +1.28 / +0.54, north −0.14 / +0.77; at P1–P4 V_helio +11.9…+12.7. FWHM_int 30–32. Paired on−off Δv +0.7 ± 1.8 (S), −0.9 ± 1.7 (N) | identical on and off | `verify_hiwham.md` 4a |
| Paper red arm [N II]/[S II] | +11.7…+13.3 (V_LSR −1…+1.7) | agrees with WHAM | `paper_review`; `verify_velocity` |
| Paper Hα | +16.3…+19.7 | 4–7 km/s too red: geocoronal over-subtraction (plausible) | `verify_velocity` |
| Local HI on the spine (GALFA) | V_LSR peak +0.4…+5.5 (median +1.8) at Dec 12–22; a +8.7…+11.3 component dominates north of Dec ~23 | no HI velocity signature on the spine (on−off moment-1 −0.16 ± 0.73) | `verify_hiwham.md` §2.4 |
| Local Ca II/Na I | +6…+15 (main local), +18…+20 (LIC/Aur) | present from 63 pc | `spectra_ism.md`; `verify_distance` |

### 5.2 Adopted

- **V_helio = +13.3 ± 1.0 km/s; V_LSR = +1.2 ± 1.0** (standard solar motion), +1.5 (Schönrich+10).
- This is the weighted mean of DESI +14.4 ± 1.9, LAMOST Fig 3a +13.6 ± 1.2 and SDSS +9.6 ± 2.8 (χ² 2.2 for 2 dof). Source: `verify_velocity_summary.json`, checked.
- **Estimator systematic (derived here).** Using `verify_desi` rev-2 (+11.3) and `verify_sdss` (+10.0 ± 2.0) instead gives **+12.4 ± 0.9**. Using DESI positive fibres +12.2 gives +12.5 ± 0.9. Quote **±1 km/s systematic**.
- **Consequence.** The [O III] moves with the WIM (WHAM +1.3), the dominant local HI (+1.2…+2.5) and the local Ca II/Na I gas. There is no kinematic separation from the red lines and no kinematic distance.
- **Superseded:** paper +23.8…+26.8 (artefact); wave-1 provisional "+12.5 to +14.5" (refined).

### 5.3 North vs south (PENDING `northsouth`)

- DESI: Dec ≥ 19.5 gives **+20.9 (18.5–23.1)** against **+10.6 (10.2–12.0)** at Dec < 19.5, a difference of 10.3 km/s in V_helio.
- SDSS north: +16.4 ± 5.9 (`verify_sdss`) / +15.5 ± 3.7 (`verify_velocity`).
- **Derived here, LSR.** With the spine-mean corrections (−10.3 north, −12.3 south):
  - north V_LSR ≈ **+10.6 (8.2–12.8)**;
  - south ≈ **−1.7 (−2.1…−0.3)**;
  - ΔV_LSR ≈ 12 km/s.
- **Canon cross-match (not tested statistically).** The northern [O III] V_LSR coincides with the HI component that dominates north of Dec ~23 (V_LSR +8.7…+11.3; `verify_hiwham` §2.4). The southern value coincides with the southern HI and WIM. `verify_hiwham` was written before `verify_desi` rev 2 measured the north and states that no northern [O III] velocity existed.

### 5.4 East vs west strand

- **No significant difference.**
  - DESI rev-2: E +12.0 (11.3–14.3), W +11.4 (8.6–15.6).
  - DESI positive fibres: E +15.0 ± 1.0, W +13.7 ± 1.5.
  - SDSS: W +9.7 ± 2.2, E +11.3 ± 4.1.
- **Residual hint:** at Dec 15.5–17.5, E < W by 3–6 km/s at ≤ 2σ (positive fibres E +10.6 ± 2.2, W +16.7 ± 2.1).
- **Superseded:**
  - Wave-1 DESI E +13.0 / W +18.5. Regression centroids of local subsets are not local, because of the petal-mean subtraction.
  - `verify_sdss`'s 3.5σ W-strand tension with DESI. Resolved by the DESI revision.

### 5.5 Near the LAMOST regions (DESI positive fibres)

- Within 0.5°: **P1 +11.9 ± 1.9** (n = 61), **P2 +15.2 ± 2.9** (n = 9), **P3 +9.2 ± 1.0** (n = 128).
- Within 0.3° of P1 (`verify_velocity`): +7.5 ± 2.2 (n = 66).
- **LAMOST-corrected minus DESI-0.5° (derived here):**
  - P1 +4.8 ± 4.8 (1.0σ);
  - P2 +1.5 ± 5.1 (0.3σ);
  - P3 +4.5 ± 3.6 (1.3σ);
  - with the 0.3° P1 value, +9.2 ± 4.9 (1.9σ).

### 5.6 Line width (conflicting; §12 C7)

| Source | FWHM (km/s) | Status |
|---|---|---|
| DESI (wave 1, not re-verified) | ≤ 33: 32.5 (stat 28–36; 26–38 for ±2 % LSF); on-ridge 18 (0–27) | measurement |
| Paper Fig 3a example | ≈ 50 (σ_obs 0.44 Å vs LSF 0.271 Å) | measurement |
| Paper Table 2 | 43–49 | not a measurement |
| Shock model | 23–38 | prediction |

---

## 6. Line ratios (energy unless stated; relative to [O III]5007 unless stated)

**Critical systematic floor.** `verify_desi` new finding 3a:
- The sky line [O I]5577, which must give 0, gives an IV "ratio" of **0.09–0.11**.
- A null box at 4381.6 Å gives **0.048 ± 0.014**.
- So **any DESI ratio below ~0.05 of [O III] is effectively an upper limit.** This applies to Hβ, [S III], [Ar III], [O II]7325 and [N II]5755, and it makes the 4363 value marginal.

| Ratio | DESI (`verify_desi`, rev 2) | SDSS (`verify_sdss`) | Adopted / status | Superseded |
|---|---|---|---|---|
| **[O III]/Hβ** | **> 11 (box; JSON 10.53) / > 18 (Gaussian; 17.66)** at 2σ. Hβ/[O III] = 0.012 ± 0.041 (box), 0.021 ± 0.018 (Gaussian), null-inflated. Hβ is only a 1–2σ signal | South **23.5 photon** (16–84 %: 12.7–47.6; median 21.6; ≈ 22.8 energy). 2σ: > 10 (bootstrap), > 7.5 (+placebo MAD), > 5.3 (+placebo rms). W strand 22.6 (13.6–42.2). North 6.4 (3.2–12; > 2.7) | **≥ 11 (2σ); central undetermined (plausibly 15–25)**. LOWER LIMIT, medium | DESI first pass 47 (> 24), south 46 (> 28): **withdrawn** (uninflated bootstrap). Wave-1 DESI 20 (15–26): only the limit survives. Wave-1 SDSS > 13.8: placebo. Altmodel caveat: Hγ/Hβ = 1.5 (case B 0.47) → ~15 if Hγ is right |
| **R23** | **> 29 (box) / > 34 (Gaussian)** at 2σ | south ≈ 76 (> 20 at 2σ) | **> ~30**. LOWER LIMIT | ~128 (> 84) first pass withdrawn; 58 (> 35) wave 1 |
| **[O II]3727/[O III]** | all 1.37 (68 %: 1.21–2.07) / box 1.41 (1.27–2.03). **South 1.11 (1.07–1.18) / 1.16 (1.09–1.24), ±0.08**. North 4.4 (3.5–5.2) / 4.2 | BOSS only: **south 1.84 (95 %: 1.46–2.37; ±0.35 with placebo)**, W 2.04, E 2.01, **north 7.3** | **South 1.1–1.8 (conflict, §12 C2); north 4–7.** Rises ×4 to the north in both instruments | wave-1 DESI 1.49; 1.12 (Dec < 19); north arc 6.2. SDSS imaging u/g implies 0.8–1.0 (qualitative) |
| He II 4686 | < 0.026 (97.5 %); box < 0.013 | < 0.04 (bootstrap), < 0.08 with placebo | **< 0.026**. UPPER LIMIT, high | < 0.034 (wave 1) |
| [O III]4363/5007 | 0.083 (Gaussian, ±0.029 inflated), 0.091 (box, ±0.051), 0.057 (pixel regression; 16–84 % 0.015–0.069); south 0.092; 1.8σ (box) – 2.9σ (Gaussian). Centre 4364.55–4364.73 Å (expected 4364.62); Hg I 4359.55 excluded (box coefficient 0.004 ± 0.019) | 0.038 ± 0.035 or 0.044 ± 0.037 (inconsistent within the note); non-constraining; Hg cancels in ON−OFF | **0.06–0.09, HINT (2–3σ)**, low | 0.056 (0.022–0.077) wave 1; "3–5σ" first pass |
| [O II]7320+30 | **0.036 (0.033–0.042)** masked fit; box 0.081 is OH-contaminated → 7325/3727 = 0.026 (0.017–0.038) | — | below the ~0.05 floor → **treat as an upper limit**, low | 0.037 wave 1 |
| [Ne III]3869 | 0.099 (0.057–0.109) Gaussian; 0.069 (0.033–0.087) box | — | **0.07–0.10**, measurement, low–medium | 0.10 |
| [S III]9069 | 0.028 (0.024–0.048) Gaussian (chromatic +12 %); 0.037 box | — | ~0.03, below floor (upper-limit-like) | 0.031 |
| [Ar III]7136 | 0.013 (0.011–0.020) | — | below floor | 0.018 |
| [N II]5755/6583 | 0.14 (0.06–0.26), ~1σ | — | no T constraint (at or above the infinite-T limit) | — |
| [O III]5007/4959 | 3.04 (2.94–3.10); box ~2.9 | free 3.15 | consistent with 2.98 | 3.11 |
| **[N II]6583/Hα (filament)** | direct 1.05 (box) – 1.21 (Gaussian; 1.12–1.43); south 1.0–1.08 | south 0.70 (0.55–0.92); north 1.75 (1.3–2.6; 1.1–2.4 across variants; 2σ ≳ 0.7) | **~1 (0.7–1.2); SYSTEMATICS-LIMITED (×2)** | DESI 0.5–0.7 (wave 1, Hβ rescaling) **refuted**; SDSS north 1.85 ± 0.36 (errors omitted the Hα systematic) |
| **[S II]/Hα (filament)** | 0.70–0.75; south 0.51–0.52 | south 0.62 (0.48–0.82); north 1.86 (1.4–2.7). [S II]6716/6731 ON−OFF 2.5 (S) / 1.8 (N) is unphysical | **0.5–0.75 south; SYSTEMATICS-LIMITED** | DESI 0.3–0.55 (wave 1) refuted |
| **[O III]/Hα (south)** | direct > 3.4 (raw p97.5) to > 8.8 (south); nominal ~10 photon (13.7 energy). Hβ + case B > 2.8 photon with the rev-2 box limit (derived here) | Dec < 17.6 with placebo **> 3.0**; Dec < 20: 2.7 (2.2–3.4; +placebo 1.9–4.5) | **≳ 3 (2σ, photon; ≳ 4 energy)** co-spatial. Non-co-spatial (Hα displaced ≤ 0.8°): SHASSA > 1.9 / 2.3 / 4.7 for [O III] 1.6 / 2 / 4 R. LOWER LIMIT, medium | WHAM "Hα/[O III] < 0.10" **refuted** (wrong velocity). Paper_review "> 13 at P3" used 10 R. SDSS 4.7 ± 1.0 (biased). hiwham's "> 4.8 via DESI Hβ" used the withdrawn > 18 |
| Hα/Hβ (filament) | 3.9 (2.2–6.1), undetermined | — | undetermined | 1.9 (wave 1, biased) |
| **[O II]3729/3726** | **1.30 (1.24–1.33)** all and south; subsets 1.17 (QSO/z > 0.5) – 1.37 (neighbour instrument). Not re-done with the box method | south 1.29 (68 %: 1.07–1.59); **all-strand 1.13 (0.99–1.29); Hα patch 1.12 (1.00–1.27)** | **1.30 (DESI), low confidence; PENDING `oii`**. Low-density limit 1.46 (CHIANTI 7.1.4) / 1.467 (shock model) | wave-1 1.30 (1.24–1.39); wave-1 on-ridge 1.45 vs off-ridge 1.12 |
| Other limits (DESI 97.5 %) | [Ca II]7291 < 0.004; [Fe II]7155 < 0.0014; [N I]5200 < 0.043; He I 5876 < 0.005; no [Ca II]7324 | — | UPPER LIMITS. They argue against a cool recombination/neutral zone | — |

### Hα patch (115.47, +20.54) spectrum

| Source | Values |
|---|---|
| SDSS | [O II] 3.55 R, Hβ 0.37 R, Hα 1.41 R; **[O II]/Hβ ≈ 13** (energy); Hα/Hβ 3.8 photon = **2.84 energy** (case B, unreddened); **[S II]/Hα 1.06 (0.95–1.19), [N II]/Hα 0.78 (0.70–0.89)**; [O III]/Hβ 1.8; plane fit 0.99 / 1.29 |
| DESI (wave 1) | [N II]/Hα 0.65 (0.50–0.71), [S II]/Hα 0.83 (0.56–0.91), [O III]/Hα 0.05 |
| NSNS | Hα +3.76 R, [S II] +2.88 R ([S II]/Hα ≈ 0.77), [O III] +0.40 R |

Shock-like low-ionisation spectrum (qualitative; no model run). The cold clump PGCC G199.09+19.95 is 0.35° away (§8).

---

## 7. Temperatures, densities, pressures, timescales

| Quantity | Adopted | Basis | Status, confidence |
|---|---|---|---|
| T_e([O III]) | **(3–6)×10⁴ K**: 4.6×10⁴ (3.1–5.3×10⁴ O&F approx.; 3.1–5.8×10⁴ CHIANTI) for 0.083; 5.4×10⁴ south | 4363/5007 hint | derived, low. Supersedes ~2.9–3×10⁴ K (from 0.056) |
| T([O II]) | 1.25×10⁴ K (1.0–1.7×10⁴); 1.18×10⁴ at n_e = 100 | 7325/3727 = 0.026, which is below the floor | derived, low (effectively an upper bound) |
| Minimum T for [O III]/Hβ ≥ 15 | ≥ 1.3×10⁴ K (13 050; 12 600–13 850) at [O II]/[O III] = 1.37, solar O/H. Photoionisation equilibrium reaches only 7–12×10³ K | altmodel thermometer | derived, medium. Supersedes paper_review "≳ 1.6×10⁴ K" (based on 10 R) |
| Two-temperature structure | T(O III)/T(O II) ≈ 4, the signature of a cooling flow | altmodel | derived, low (both inputs marginal) |
| **n_e from the [O II] doublet** | 88 / **100** / 125 cm⁻³ at T = 1.0 / 1.25 / 2.0×10⁴ K; 47–208 over subsets (altmodel); ~110 (40–300) (verify_desi). ±30 % from atomic data (low-density limit 1.42–1.50 across datasets) | 3729/3726 = 1.30 | derived, low; **PENDING `oii`** |
| **n_e from the emission measure** | 0.1–0.7 cm⁻³ (edge-on shell/sheet, L 25–50 pc); 0.5–3 (filament, L = 1.37 pc); EM 0.27–10.3 pc cm⁻⁶ (2–4 R, 280 pc, T 1–5×10⁴ K) | 5-level O²⁺ atom | derived, medium. Supersedes n_e 1–3 (for 10 R) |
| Shock-model densities | pre-shock n₀ 0.012–0.042; post-shock [O II]-zone n_e 0.14–0.25 (doublet 1.467) | `verify_shockmodel_summary.json` | model, medium |
| [S II] n_e from Table 2 (WIM) | 225–631 cm⁻³ (5-level); 169–400 (eq. 5). Per-region 16–84 % span 11–1900, so statistically meaningless and not filament gas | verify_paper | derived, high (as a refutation of the paper's 64–150) |
| **Pressure P/k (K cm⁻³)** | doublet ~2.8×10⁶ (1.4–5.7×10⁶); shock ram pressure (width + brightness closure) 3.1×10⁴–1.1×10⁵; Monogem hot gas 1.5×10⁴ (altmodel), 0.8–2.9×10⁴ (verify_paper), 3.7×10⁴–1.6×10⁵ × f^−½ (shockmodel) | definitions differ (§12 C10) | derived, low. "100–400× Monogem" holds for altmodel's 1.5×10⁴; against shockmodel's range it is 18–150× |
| Emitting path at the doublet density | 40–230 AU (filling factor ~10⁻⁵ of a 25-pc limb path) | altmodel | derived, low |
| Timescales | shock exposure t* 3.9×10⁴–1.04×10⁵ yr < full cooling time 0.68–1.4×10⁵ yr. Monogem age ~1.1–1.2×10⁵ yr (Knies+24). Fossil e-fold time 21–34 yr at the doublet density. Precursor t_rec ~5 Myr, t_cool ~6 Myr | shockmodel; altmodel | model, medium |
| Shock temperature | T_s = 1.2–5.2×10⁵ K, so nothing above 0.2 keV | shockmodel | model |

---

## 8. Multiwavelength

| Band | Adopted result | Numbers | Superseded (why) | Source |
|---|---|---|---|---|
| **GALEX FUV** | **DETECTED**, filament-shaped, co-spatial with [O III] (south) | Template regression: **92 ± 10 CU per NSNS-R south (z = 8.7)**, 74 ± 8 whole spine; 68 ± 8 per true R if [O III] = 1.35 × NSNS. Direct contrast +120 ± 26 (south, z = 4.6), +67 ± 23 (all). Stacked peak 145 CU at NSNS 1.55 R. Null p ≤ 0.015 | "weak excess at most", 2σ UL 110–250 CU (wave 1: same numbers, now significant) | `verify_xrayfuv` §1.7 |
| FUV interpretation | Rules out complete shocks (non-resonant floor 85–277 CU/R) and the paper's 10–13 R (would be 12–15 CU/R, below every non-resonant floor). Matches an incomplete shock with **C IV suppressed ~5×** | xrayfuv 6.5–26× (≥ 5 at 1σ); shockmodel RT S₁₅₄₈ 0.16–0.27 (~3.7–6×) → 52–126 CU/R; altmodel ≥ 3–5×. R88 floor 57 CU/R | — | `verify_xrayfuv`; `verify_shockmodel` §6 |
| **GALEX NUV** | **DETECTED** | 49 ± 7.5 (south), 57 ± 5.5 (whole) CU per NSNS-R; joint fit 38 ± 8; peak ratio ~80 → **36–57 (up to 80)**. Direct +96 ± 23 (S), +68 ± 12 (all). Matches incomplete-shock NUV lines (47–63). Dust < 16 CU; excess colour 1.0–1.9 vs background 0.50 | "nature unclear; colour ≈ background 0.71" | `verify_xrayfuv` §1.10 |
| UV coverage | GALEX stops at Dec 22, so the northern Silk is untested in the UV | — | — | — |
| **Soft X-ray** | **No X-ray filament.** The Silk lies on the **outer flank** of the Monogem 0.2–0.6 keV shell, ~1° outside the ring peak (r ≈ 10.25°) | On-spine bump −0.04 (all) / −0.38 (S) / +0.07 (N) cts s⁻¹ deg⁻² on a level of 5.4–7.1. **W/E 1.26** (S 1.22, N 1.25) eRASS1; 1.24 ROSAT. **p = 0.007 vs random sky, p = 0.37 vs same-radius rotations** (south 0.57): generic to the shell. Local step south +12 % at x = −0.2°, p 0.013–0.08 (~2σ), 0.04–0.25 with absorption; north none. The Silk-sector edge is softer (×1.27 per degree) than the S/SW limb (×1.43–1.57) | W/E 1.32, 2.4σ (parallel-shift nulls), edge 1.03 ± 0.20° inside, "the Silk is the Monogem X-ray edge" (wave 1) | `verify_xrayfuv` §1–2 |
| X-ray absorption | The south interior has **more** HI (+5.0×10¹⁹, HI4PI), which works against the observed sense | — | ΔN(HI) −1.3×10¹⁹ (wave 1) | `verify_xrayfuv` |
| **Hα** | **No filament** | **SHASSA (Dec 12–15.6, 4′): 2σ < 0.77 R** (Sm) / 1.02 R (CC) on the intrinsic peak. Displaced or broadened (\|d₀\| ≤ 0.8°, FWHM 0.28–0.5°): < 0.85 R. F03 SHASSA part < 0.78 R. **F03 north of Dec 15.9 is WHAM-only (1°): 2σ < 1.74 R**. WHAM template south < 1.0 R (±25 km/s) / 1.25 R total; north 1.25–1.30 R | F03 < 0.43 R over the whole spine (invalid 6′ correction; F03 is WHAM-only on 76 % of the spine); SHASSA 3σ < 0.8–1.0 R (now 1.1–1.3 R); F03 "+0.8–0.9 R at Dec 17–21" (a 1° WHAM feature); WHAM −9.8 km/s "3σ" residual (now ~1σ) | `verify_hiwham` §3 |
| **HI** | **No filament, wall, cavity or velocity signature** | 2σ UL 1.0–1.4×10¹⁹ cm⁻² per 2-km/s slice at V_LSR −3…+5 (whole spine), 1.5–1.75×10¹⁹ (south), 1.1×10¹⁹ at +11…+16, ~3.0×10¹⁹ for the 6-km/s heliocentric window. On-spine N(HI, all v) 4.0×10²⁰. All slices \|z\| ≤ 1.4 | **"HI wall 0.86 ± 0.15° outside" refuted**: a finder artefact (random sky gives +0.84 ± 0.31°, z = −0.2). The RESUME headline should be dropped | `verify_hiwham` §2 |
| PGCC G199.09+19.95 | Cold clump 0.35° from the Hα patch | CNM lines at **V_LSR −0.74 and +6.46 km/s** (FWHM 2.7, 3.0; ΔT_B ~7 K; T_k ≲ 160 K). ΔN(HI) +5.4×10¹⁹. 2.7 % of random positions show an excess this strong | The +6.5 coincidence with a W-strand +18.5 km/s velocity is now inconclusive | `verify_hiwham` 8b |
| HI shells | GS202+17+019: weak LAB cavity (z 1.7–2.8 at +10…+40), **not** at the revised [O III] velocity. GS202+19−019: not confirmed | — | "centred on the Silk" (catalogue coincidence) | `verify_hiwham` 8a |
| Dust / IR | No counterpart | IRIS 100 µm 2σ ΔI < 0.19 MJy/sr; ΔE(B−V) < 0.004 → N_H(dust) < 2.3×10¹⁹. Foreground E(<300 pc) 0.0040–0.0079 along the spine (0.0051–0.0061 at P1–P4) | — | `imaging_mw` (wave 1); `context` §5 |
| Radio (not re-verified) | No total-intensity counterpart; Stockert ambiguous | Haslam 408 MHz −0.07 ± 0.62 K → **Σ₁GHz < 1.3×10⁻²² W m⁻² Hz⁻¹ sr⁻¹** (2σ, beam-diluted, α = 0.5); WMAP K < 1.0×10⁻²²; Planck 30 GHz < 1.3×10⁻²²; DRAO 22 MHz < 2.7×10⁻²². **Stockert 1420 MHz +19.5 mK (4.0σ)** on a W > E step (+68 mK, 4.6σ), also seen in Planck 30 GHz PI (4.2σ). NVSS/TGSS resolve it out | — | `imaging_mw` §1.5 |
| RM (not re-verified) | North-half step; not claimed as a filament signature | NVSS (Taylor+09). **North (Dec 20–28): E(2–6°) − W(2–6°) = +36.8 ± 4.6 rad m⁻²**, step at +2.0° E (−0.6…+2.0), 70 % of RMs negative on the W side vs 16 % on the E. South: +5.3 ± 7.3. On−off +19.2 ± 8.1 (2.4σ). A foreground gradient could do the same | — | `spectra_ism` §3.4; `spectra_ism_rm_results.json` (checked). The `northsouth` teammate has RM products in progress (PENDING) |
| Gamma | None | Fermi-LAT 0.3–1 / 1–3 / 3–300 GeV: −1.0 / −1.5 / −0.2σ; no TeV-halo modelling | — | `imaging_mw` |
| UV absorption | No Silk component separable | COS WD J0747+2106 (99.5 pc, 0.16°): log N(Si IV) < 12.2, N V < 13.0 (3σ; wave 1). IUE HD 60848 (1.9 kpc, 2.26° W): S II 15.01, Al III 12.86, Si IV 12.80–12.86, C IV 13.32–13.49, N V < 12.88; high ions 5–11 km/s blueward of low ions; nothing at \|v\| 70–200 | "122 unmeasured IUE spectra" (they are 14 SWP observations) | `verify_distance` b1–b2 |
| Pulsars | Not behind the filament | J0751+1807: PX 0.85 ± 0.04 mas → 1.18 (1.12–1.24) kpc, RM +41.6 ± 0.2, DM 30.24, 0.95° E. B0656+14: 288 +33/−27 pc, DM 13.95–13.98, RM +22.7 (ATNF 2.9.0) / +23.5 ± 0.4 (B/psr 2017). A Silk ΔDM ≲ 3 pc cm⁻³ for any n_e | J0751: PX 1.6 ± 0.8, no RM (stale VizieR) | `verify_distance` c4 |

---

## 9. Distance

| Constraint | Range | Strength | Source |
|---|---|---|---|
| **Monogem limb (south spine 12.10° from the two-SNR centre)** | **252–318 pc (68 %), 225–359 (95 %), median 282**; shell radius 61 (54–68) pc | strongest, but **conditional and morphological** | `verify_distance_synthesis.json` (checked) |
| Monogem shell near/far side (shell up to 1° larger) | 237–337 pc (near 237–300, far 266–337) | allowed variant | same |
| Local Bubble wall (O'Neill+24) | peak 161–171 pc toward P1–P4; spine median 167 (158–223), south 165, north 178; inner 139, outer 200 | reference only: Wu+'s 160 pc sits on it by coincidence | `context` §5; `verify_distance` c3 |
| HD 61717 −17 km/s cloud | < 272 (284) pc *if* the cloud belonged to the shell | weak, conditional. The cloud is cool (b 1.49 km/s), Ti-poor (Ti II/Ca II < 0.16), compact (absent toward 14 stars) and at V_LSR −28 vs [O III] ~+1, so unrelated. The star is an **SB2** (−69/+124 km/s); d = 272 (262–284) pc with RUWE 2.14 | `verify_distance` a1–a5 |
| WD J0747+2106 COS | > 99.5 pc *if* the Silk carried log N(Si IV) > 12.2 / N V > 13.0 | weak, conditional | `spectra_ism`; `verify_distance` |
| Ca II/Na I at the [O III] velocity | none: the [O III] velocity equals the ubiquitous local gas | none | `verify_distance` c7 |
| Kinematic | 0–0.21 kpc or 0–1.45 kpc (revised V); 0.97–3.0 kpc (LAMOST, refuted) | meaningless | `verify_paper`; `verify_distance` |
| Extinction / 3D dust | none (E(B−V) < 0.006 to 300 pc; no 250–320 pc feature) | none | `context`; `verify_distance` |
| Pulsar DMs | none | none | `verify_distance` |
| Size and height plausibility | at 1 kpc: 295 pc long, 4.9 pc wide, z = 330 pc; at 1.6 kpc: 473 pc long, 7.7 pc wide, z ~0.55 kpc | soft upper bound ~1 kpc | `verify_distance`; derived |

**Adopted.**
- **252–318 pc (68 %; median 282) if the southern spine is the Monogem limb.**
- Permitted otherwise: ~100 pc to ~1 kpc.
- **160 pc is not excluded by any measurement.**
- **~1.6 kpc** (the corrected `kd` output) is **disfavoured** by size, height and the loss of the Monogem registration.
- No absorption line brackets the ionised gas.

**Caveat (§12 C12).** The synthesis JSON gives the wave-1 X-ray registration (W/E 1.32, 2.4σ) as the basis of the distance. `verify_xrayfuv` shows that contrast is generic at that radius (p = 0.37 against rotations). What remains is the radius coincidence: 12.10° against R = 11.0 ± 1.1°, about 1° outside the ring peak, together with Fesen+26's assignment. The north half (Dec > 17.6) is not on the Monogem sphere.

---

## 10. Precovery table

| Date (UT) | Dataset | Result | Significance | Source |
|---|---|---|---|---|
| 1950–1996 | DSS1/DSS2 (POSS-I/II) | not detected; DSS2 J A = 4.3 ± 3.2 → < 14 × NSNS (≈ < 40 R) | 3σ limit | `imaging_optical` |
| 1997–2000 | SHASSA / WHAM Hα | no Hα filament (not an [O III] test) | — | `verify_hiwham` |
| 2001-11-12 – 2002-12-06 | SDSS imaging, north runs (Dec ≥ 20.5) | 2002-01-14 (runs 2887/2888) g 3.9σ with u and r excess: **broadband confounder, not attributable** | — | `verify_precovery` §2.3 |
| 2002-02-09 | SDSS plate 857-52314 | 0.18 ± 0.05 R formally 3.5σ; **0.5σ with placebo** (verify_sdss); 1.8σ (verify_precovery) | not detected | `verify_sdss` §7 |
| 2002-10-12 | SDSS plate 927-52559 | 2.7σ (verify_precovery); 0.50 ± 0.23 R, 2.1σ formal (verify_sdss total sky) | < 3σ | CSV / JSON |
| **2002-10-30** | **SDSS plate 927-52577** | 1.03 ± 0.21 R, **5.0σ** in the verify_precovery per-plate fit. Uses wave-1 residual fluxes, only 3 on-ridge fibres, not placebo-validated, absent from the verify_sdss table | **earliest ≥ 3σ: FRAGILE** | `verify_precovery_table.json` first_3sigma (checked) |
| **2002-12-31** | **SDSS imaging run 3560** (g, PSP, Dec 20.2) | 3.77 ± 1.16 R; same run r −0.4σ, i −0.9σ | **3.3σ** (single-run σ ±~0.5): **first imaging ≥ 3σ** | `verify_precovery` §3 (checked) |
| 2003-01-25 / 02-01 | SDSS imaging, cumulative core runs | 3.9σ by 2003-01-25, 4.6σ by 2003-02-01 | core-only cumulative (the note; see I5) | `verify_precovery` §3 |
| **2003-01-30** | **SDSS plate 1203-52669** (04:42–05:20 UT) | 1.11 ± 0.11 ± 0.31 R | **3.4σ (verify_sdss, placebo) / 3.5σ (verify_precovery): earliest confirmed by two independent analyses** | `verify_sdss` §7; CSV (checked) |
| 2003-03-10 | SDSS plate 1263-52708 | 3.9σ (verify_precovery) vs 0.32 ± 0.13 R, 2.5σ formal (verify_sdss) | conflicting | CSV |
| 2003–2012 | GALEX FUV/NUV (GR6/7; per-tile dates not extracted) | UV filament (FUV z = 8.7 south): a UV-counterpart precovery, not [O III] | — | `verify_xrayfuv` |
| 2004-02-17 | SDSS imaging u (run 4469, [O II]) | first core u-band ≥ 3σ; combined u 1.40 ± 0.24 R [O II], 4.4σ coherent | 3.1σ | `verify_precovery` §3 |
| 2004-11-08 | SDSS plate 1921-53317 | profile fit 1.88 ± 0.35 R, no fibre within 0.15°, not placebo-tested | formal 5.4σ (untested) | `verify_sdss` §7 |
| 2005-03-02 | SDSS plate 2054-53431 | 1.20 ± 0.07 R | 4.2σ | `verify_sdss` §7 |
| **2005-03-03** | **SDSS plate 2074** (MJDLIST 53432; coadd MJD 53437) | 1.45 ± 0.05 R | **9.7σ: first ≥ 5σ** | `verify_sdss` §7; CSV (checked) |
| 2008-01-30 / 02-01 | SDSS plates 2890 / 2915 | 1.36 / 1.44 R | 5.4σ / 13.9σ | `verify_sdss` §7 |
| 2010-12-06 – 2011-10-31 | BOSS plates 4456, 4497, 4492, 4502, 4488, 4470, 4473, 4462, 4464 | 0.38–2.51 R | 3.2–13.9σ; verify_precovery's first BOSS detection is 2011-01-08 (4502, 4.7σ) | `verify_sdss`; `verify_precovery` |
| 2013-03-10 | LAMOST LRS stellar residuals (GAC113N27M1) | 3.1σ; 4.4σ on 2013-11-06; up to 9.3σ (2016-01-30) | ≥ 3σ | `verify_precovery` |
| 2016-03-06 | SDSS plate 8954 | 0.73 ± 0.06 R | 5.2σ | `verify_sdss` |
| 2018-10-10 – 2019-02-01 | eBOSS plates 11082–11104 | 0.48–2.33 R | 4.6–16.3σ | `verify_sdss` |
| 2021-01 (first night 2021-01-03) | DESI DR1 | per month 5.7–48σ (partly circular: regression on the forward model) | ≥ 5σ | `verify_precovery` |
| 2022-12-12 – 2025-04-04 | NSNS DR0.2 | detected; per-pixel epochs not public; "first clear night after full configuration" ~2024-01 | — | `imaging_optical` |
| 2023-11 / 2023-12 | SWONS | discovery (self-reported) | — | `discovery` |
| 2024-01-14 | AstroBin zh47r9 | first public image | — | `discovery` |

**Summary.**
- **Earliest ≥ 3σ:**
  - **2002-10-30** (plate 927, fragile);
  - **2002-12-31** (imaging run 3560, 3.3σ);
  - **2003-01-30** (plate 1203, the earliest confirmed by two independent analyses; adopt as the **secure earliest**).
- **First ≥ 5σ: 2005-03-03** (plate 2074, 9.7σ). Plate 1921 (2004-11-08) would precede it if its formal 5.4σ survives a placebo test.
- The SDSS u-band sky also shows [O II].
- **No variability** over 2002–2022 beyond a factor of ~2.
- The filament was present **~21 years** before its 2023 discovery.

---

## 11. Model scenario verdicts

| Scenario | Verdict | Decisive reasons | Source |
|---|---|---|---|
| Complete steady radiative shocks (80–200 km/s) | **RULED OUT** | [O III]/Hβ ≤ 5 at completion (own code 2.4–5.0; R79 2.8–4.5, CR85 2.5–3.4, DS96 1.2–1.9, Allen+08 1.2–3.2) vs DESI > 11 (2σ). FUV non-resonant floor 85–277 CU/R > 68–92 measured. Would add an Hα filament and cooled HI. *Caveat:* with SDSS alone the 2σ floor is only > 5.3 | shockmodel; xrayfuv |
| **Incomplete (column-truncated) radiative shock, driven by Monogem hot gas** | **FAVOURED (viable)** | v_s **100–150 km/s** (best 110–140); log N* **17.6–17.9**; T* 1.7–2.2×10⁴ K; **n₀ 0.012–0.042** (2–8× the Monogem ambient); **b = B₀/√n₀ ≈ 7–14 µG cm^3/2** (B₀ 0.8–2.9 µG); ram P/k 3.1×10⁴–1.1×10⁵ (≈ Monogem). I_pk **1.3–4.0 R** (median 2.3) for R_c = 58 pc at 280 pc (G ≈ 20). FWHM 23–38 km/s. Centroid = pre-shock gas velocity. t* 3.9×10⁴–1.04×10⁵ yr < t_cool 0.68–1.4×10⁵ yr. FUV 52–126 CU/R with C IV scattering. Proper motion 75–113 mas/yr outward. **Tensions:** [O II]3729/3726 (model 1.467 vs 1.30) and [O II]7325 ×1.8 (weaker if 7325 is an upper limit). **Predictions:** sharp east edge with a 1/√h tail toward Monogem; FUV ridge ~0.1° outside the [O III] ridge. The fit used the now-withdrawn [O III]/Hβ = 40; the consistent family spans 10–63, so a re-run is advised but the verdict is unlikely to change | `verify_shockmodel`; `verify_shockmodel_summary.json` (checked) |
| Photoionised precursor as the [O III] emitter | **RULED OUT** (as a pre-ioniser: viable, part of the truncated-shock picture) | T_e ≤ 1.5–2×10⁴ K, so 4363/5007 ≤ 0.02–0.034 and [O III]/Hβ ≤ 10.7. Dark at n₀ ~0.03. Supplies 0.2–1.6 ionising photons per incoming H, hence X(H⁰) ≲ 0.3 | shockmodel; altmodel |
| Photoionisation by WD J0736+1444 (81 kK, 150.5 pc) | **RULED OUT** | U short by 27× (n_rms, coplanar) to 9×10³× (doublet n), and ≥ 10⁴× at 282 pc. PIE gives [O III]/Hβ ≤ 5–9 and 4363/5007 = 0.006 | altmodel; context |
| Photoionisation by Monogem / Local Bubble EUV | **RULED OUT** | log U ≤ −6.1 (needs ~−3.3). Even at the right U: [O III]/Hβ ≤ 10–14 and He II/[O III] = 0.06 > 0.026 | altmodel |
| ISRF / DIG photoionisation | **RULED OUT** | [O III]/Hβ ≤ 1.3 at any U | altmodel |
| Photoionisation + extra heating | **RULED OUT** | needs 78–350× the photo-heating; 7325/3727 ×3 | altmodel |
| Conduction front | **RULED OUT** | 4363/5007 0.11–0.14; FUV 540–5300 CU/R; face-on 2×10⁻⁴–3×10⁻³ R per front | altmodel; xrayfuv |
| Turbulent mixing layer | **RULED OUT** | needs a 377 km/s inflow; ~0.05 R at 10 km/s; X-ray filament 2×10³–2×10⁴× Monogem | altmodel |
| NEI cooling of Monogem hot gas | **RULED OUT** | X-ray filament ≥ 10³× Monogem predicted, none seen | altmodel |
| Warm CIE layer at ~6×10⁴ K | **DISFAVOURED** | fits [O III]/Hβ, [O II]/[O III], 4363 and FUV, but 7325/3727 ×3, NUV ×2, [N II]/Hα ×2, and no heating source | altmodel |
| Fossil / recombining structure | **RULED OUT** | e-fold time 21–34 yr at the doublet density; no fading over ~20 yr | altmodel; precovery |
| Isolated ~30-pc SNR at 160 pc (Wu+) | **DISFAVOURED** | no interior of its own; 160 pc was a unit error; invalid energetics | imaging_mw; verify_paper |
| **Southern Silk = Monogem E/NE limb (~250–320 pc)** | **FAVOURED (morphological)** | Monogem radius in the south; outer flank of the X-ray shell; Fesen+26. Not a distinct X-ray edge; the north is off the sphere | xrayfuv; distance; context |
| Geminga-progenitor SNR fragment (Ziegenbalg) | **NOT SUPPORTED** | the concentricity is non-physical; no birth locus | context |
| Local Bubble wall interface at 160 pc | **DISFAVOURED, not excluded** | needs an unseen energy source; no tie to the wall | context; distance |
| Northern Silk | **Model-consistent, not unique** | [O III]/Hβ ~6.4 and [O II]/[O III] 4–7 fit a later cut (~10⁴ K, log N 18.0–18.2) or a slower 70–80 km/s shock | shockmodel row 15 |

---

## 12. Unresolved conflicts between agents

| ID | Topic | Positions | Status |
|---|---|---|---|
| **C1** | **[O III]/Hβ central value** | DESI gives only a lower limit (> 11 box / > 18 Gaussian, 2σ). SDSS south 23.5 photon (12.7–47.6). Withdrawn DESI first pass 46–47 (still in `verify_desi_summary.json`). Altmodel's Hγ/Hβ = 1.5 anomaly implies ~15 if Hγ is right. Shockmodel fitted 40 | **UNRESOLVED.** Adopt > 11 (2σ); plausibly 15–25 |
| **C2** | **[O II]/[O III] south: DESI vs SDSS** | DESI 1.11–1.16 (±0.08, Dec < 19.5) vs SDSS (BOSS-only) 1.84 (1.46–2.37; ±0.35 with placebo). Tension **2.0σ (placebo included) to 3.0σ (bootstrap only)** (derived here). SDSS imaging u/g implies 0.8–1.0. Candidate causes: Dec/strand weighting with a steep north gradient, BOSS-only groups, the 0.3 R placebo floor, u-band response | **PENDING** (`oii` or `northsouth`) |
| **C3** | **Absolute scale: SDSS imaging vs revised DESI** | imaging 2.31 ± 0.38 R → obs/DESI_rev = 1.49 ± 0.25 (2σ); imaging / SDSS-fibre_rev 1.1–1.3. Causes: PSP smoothing (2.8–6.0 R with smoothed templates), ~25 % template systematic, DESI petal-mean proxy and the f_ps model | **UNRESOLVED.** Does not affect the 4–8× verdict on 10–13 R |
| **C4** | **DESI P1 velocity vs LAMOST-corrected** | DESI within 0.3° (wave-1 fibres) +7.5 ± 2.2; within 0.5° (positive S/N > 4) +11.9 ± 1.9; LAMOST corrected +16.7 ± 4.4. Tension 1.9σ or 1.0σ depending on the estimator; P2 0.3σ; P3 1.3σ (derived) | **Mostly resolved.** Aperture/estimator difference; the LAMOST correction is aggregate and rests on one spectrum |
| **C5** | **[O II] doublet density** | 1.30 → n_e ~100 (47–208), P/k ~3×10⁶, versus EM densities 0.1–3 and every pressure-balanced model (1.46–1.467; n_e ~0.2). SDSS all-strand and patch give 1.12–1.13 (denser still); SDSS south 1.29 (1.07–1.59) is weak. DESI subsets span 1.17–1.37 | **PENDING (`oii`)** |
| **C6** | **North–south [O III] velocity** | DESI north +20.9 (18.5–23.1) vs south +10.6; SDSS north +16.4 ± 5.9 / +15.5 ± 3.7. ΔV_LSR ≈ 12 km/s after the along-spine LSR correction. North V_LSR ~+10.6 matches the northern HI component (+8.7…+11.3) | **PENDING (`northsouth`)** |
| C7 | [O III] line width | DESI ≤ 33 (on-ridge 18, 0–27) vs paper Fig 3a ~50 vs paper 43–49 (biased) vs model 23–38 | UNRESOLVED (DESI width not re-verified in wave 2) |
| C8 | [O II]7325 / T(O II) / model ×1.8 | DESI 0.036 lies below its own ~0.05 floor; box 0.081 is OH-contaminated; model 0.059–0.068 | UNRESOLVED; treat 7325 as an upper limit |
| C9 | C IV suppression factor | xrayfuv 6.5–26× (≥ 5 at 1σ); shockmodel ~3.7–6×; altmodel ≥ 3–5× | Consistent at ~5×; model-dependent |
| C10 | Monogem pressure | 1.5×10⁴ (altmodel) / 0.8–2.9×10⁴ (verify_paper) / 3.7×10⁴–1.6×10⁵ × f^−½ (shockmodel) | Definitional. Quote with the definition; the doublet/Monogem ratio is 18–400× |
| C11 | Early SDSS plates 927 / 1263 / 1921 | 927-52577 5.0σ (wave-1 fluxes, 3 fibres), not in the verify_sdss table; 927-52559 2.7σ vs 2.1σ; 1263 3.9σ vs 2.5σ formal; 1921 formal 5.4σ, untested | UNRESOLVED. Secure earliest is 2003-01-30 |
| C12 | Strength of the Monogem-limb distance | The synthesis basis cites wave-1 W/E 1.32 (2.4σ); `verify_xrayfuv` shows it is generic (p = 0.37); the local step is ~2σ, south only | 252–318 pc rests on the radius coincidence; keep it conditional |
| C13 | Hα patch low-ionisation ratios | SDSS [N II]/Hα 0.78, [S II]/Hα 1.06 vs DESI 0.65 / 0.83–0.85 vs NSNS [S II]/Hα 0.77 | Minor (1–2σ) |
| C14 | NSNS vs LAMOST [S II] intensity | NSNS on-ridge ≤ 0.35 R vs LAMOST-scaled WIM 2.4–3.9 R (NSNS narrow-band probably loses the smooth DIG) | UNRESOLVED (wave 1, low impact) |
| C15 | DESI Hα/Hβ and Hγ/Hβ | Hα/Hβ 3.9 (2.2–6.1); Hγ/Hβ 1.5 vs case B 0.47 | UNRESOLVED (Hβ at 1–2σ) |

---

## 13. Inconsistencies found between files (numbers)

| ID | Where | Issue |
|---|---|---|
| I1 | `data/derived/verify_desi_summary.json` → `ratios` | **Still holds the withdrawn first-pass [O III]/Hβ** (47.48 all, 45.93 south, 59.8 petal FE) **and R23 128.1 / 112.5.** The rev-2 limits are only in `verify_desi_audit_extra.json` (`limits`) and the note. Use the audit file |
| I2 | `notes/verify_shockmodel.md` §0, §4 | Uses withdrawn DESI inputs ("[O III]/Hβ now 46 (> 28)", "R23 ~76–130", "observed ≥ 24") and fits [O III]/Hβ = 40 ± 0.2 dex. The conclusions are robust, but the fit should be re-run with > 11 |
| I3 | `notes/verify_hiwham.md` §4; `verify_hiwham_wham_summary.json` | "Strongest single statement: DESI Hβ + case B > 4.8" uses the uninflated first-pass > 18. With the rev-2 box limit (> 10.5) it becomes ~2.8 photon (derived). The adopted "≳ 3" survives via SDSS (3.0) and DESI direct (3.4) |
| I4 | `verify_distance_synthesis.json` → `combined` | "statement" says "favoured 255–315 pc" while `favoured_pc` = [252, 318]. The basis cites wave-1 W/E 1.32 (superseded) |
| I5 | `verify_precovery_table.json` | `cum_core_snr` accumulates north runs from 2001-11-12 (1.30σ after run 3560), whereas the note's core-only cumulative is 3.3σ after 3560, 3.9σ by 2003-01-25 and 4.6σ by 2003-02-01. `first_3sigma` for g Frames/FramesSub (2002-01-14, run 2887) and u (2002-01-09, run 2864) are north, broadband-confounded runs that the note rejects |
| I6 | `notes/verify_sdss.md` | 4363/5007 is 0.038 ± 0.035 (§1) vs 0.044 ± 0.037 (row 3i) |
| I7 | DESI global velocity | `verify_velocity` combines +14.4 ± 1.2; `verify_desi` prefers +11.0 / +11.3 / +12.2 ± 0.8 (S/N > 4) / +12.7 ± 0.8 (S/N > 5, JSON only). The adopted mean moves +13.3 → +12.4 (derived) |
| I8 | SDSS velocity | `verify_sdss` south +10.0 ± 1.8, all +12.3 ± 2.5 vs `verify_velocity` south +7.5 ± 3.4, all +9.6 ± 2.8 |
| I9 | DESI near P1 | +7.5 ± 2.2 (0.3°, wave-1 fibres) vs +11.9 ± 1.9 (0.5°, positive S/N > 4) |
| I10 | DESI limit rounding | The note says "> 11" and "> 29"; the JSON has 10.53 and 28.97 |
| I11 | `verify_distance_synthesis.json` → `vlsr` | Uses a spine-mean LSR correction (~−11.2): DESI +14.5 → +3.3. The P1-local correction is −12.12 (→ +2.4) |
| I12 | PSR B0656+14 | RM +23.5 ± 0.4, DM 13.98 (B/psr 2017; spectra_ism) vs RM 22.73, DM 13.95 (ATNF 2.9.0; verify_distance) |
| I13 | Spine–PSR separation | 11.36 ± 0.16° (context) vs ± 0.15° (verify_xrayfuv, verify_distance) |
| I14 | South-only curvature radius | 45.9° (22.5–71.4) for Dec ≤ 17.5 (geometry) vs "~124°" for Dec ≤ 20 (discovery): different subsets |
| I15 | SDSS imaging runs | 18 runs, 17 positive (imaging_optical) vs 19 runs, 16 positive (verify_precovery) |
| I16 | FUV CU → R | "110–250 CU ≈ 0.4–0.9 R" (imaging_mw, W_eff 300 Å) vs "0.6–1.4 R" (verify_xrayfuv, W_eff 357 Å, E(B−V) 0.03) |
| I17 | [O II]7325 | shockmodel constraint 0.034 vs DESI 0.036 |
| I18 | Task brief vs sources | "[O II]3729/3726 = 1.29–1.30 in both DESI and SDSS" holds only for SDSS south (68 % 1.07–1.59); SDSS all-strand 1.13, patch 1.12 |
| I19 | `verify_paper` inputs | Uses DESI 2.4 R (Hβ bound 0.12 R) and 4363/5007 0.056 (T 2.9×10⁴ K). With 1.62 R and 0.083: ~0.08 R and ~4.6×10⁴ K, which strengthens its x(O²⁺)/Balmer tension |
| I20 | Plate 927 epochs | `verify_sdss_precovery.csv` has plate 927 only at MJD 52559 (first exposure 52556 = 2002-10-09); `verify_precovery` quotes 927-52559 as 2002-10-12 and adds 927-52577 |

---

## 14. What the paper got right / wrong / could not test

**Right**
- The [O III] filament is real line emission. It is detected spectroscopically, and there is no [O III] in background P4. This is confirmed by DESI (~30σ), SDSS fibres, SDSS imaging and NSNS.
- The red lines (Hα, [N II], [S II]) are foreground/background DIG/WIM, and "DIG cannot be excluded" is correct for them.
- Shock-related excitation ("possible SNR") is defensible, though via a truncated shock and most plausibly as the Monogem limb.
- The Monogem association via Fesen is plausible.
- The arithmetic of E = 2.0×10⁵¹ erg and of the mag → R conversion is right. The eq. 2/3 coefficients are right.
- The red-arm [N II]/[S II] velocities agree with WHAM.

**Wrong**
- Distance 160 pc: Table 4 is in kpc, so the value is 1.6 kpc, and kinematic distances are meaningless here.
- Size (10°, 28 pc; R = 30 pc asserted).
- [O III] velocity (+23.8…+26.8 is −10.1 km/s too high) and the [O III]/red-line kinematic separation built on it.
- Brightness 10–13 R (~5–8× too high).
- [S II] density 64–150 cm⁻³ (its own Table 2 gives 225–631).
- Red-line SNR diagnostics, T_e, V_t and metallicity all describe the WIM.
- Energy, age and R_PDS inputs ("+1" typo; R_PDS 17.1).
- "Strong λ4959".
- Radio non-detection used as support.
- 160 pc vs Monogem unreconciled.
- Credit wording "discovered … for the first time to its full extent".
- LAMOST Hα velocity (artefact, plausible).

**Could not test with the paper's data**
- [O III]/Hβ, 4363/5007 (T_e), [O II], He II, R23 and any absolute line ratio: MRS-N has no Hβ and no flux calibration.
- The intrinsic [O III] width (S/N ~6).
- Any distance.
- UV/X-ray counterparts.

**Paper claims we could not verify**
- The SWONS photometry inputs (filter curve, camera, apertures, zero-point units), so the 10–13 R error cannot be pinned to a cause.
- Per-fibre MRS-N spectra (not public), so the −10.1 km/s correction is aggregate.
- The red-arm processing (inferred correct only from WHAM agreement).
- The `kd` aggregation rule behind the Table 4 red-line distances and the MC widths.

---

## 15. Register of superseded headline values

| Wave-1 / provisional value | Replaced by | Why |
|---|---|---|
| DESI ridge 2.4 R (S), 0.86 R (N), ≥ 2.8 R | 1.62 R, 0.56–0.71 R, ≥ 1.7 R | wrong fibre-acceptance factor |
| SDSS ON−OFF 0.79 R; W 1.19 R | 1.39 R; W 1.93 R | residual `flux` instead of total sky |
| SDSS imaging 1′ peaks 4.0–4.9 R | 17′ peak 2.3–2.6 R | NSNS morphology scaling |
| DESI [O III]/Hβ 20 (first pass 47); R23 58 (first pass 128) | > 11–18; > ~30 | Hβ is 1–2σ; null-inflated bootstrap |
| DESI filament [N II]/Hα 0.5–0.7, [S II]/Hα 0.3–0.55 | ~1.0–1.2, 0.5–0.75 (systematics-limited) | came from Hβ rescaling |
| DESI 4363 "3–5σ" | 2–3σ | null boxes |
| DESI E/W +13.0 / +18.5 | no significant difference | local regression centroids are not local |
| SDSS [O III]/Hβ > 13.8 | > 5–7.5 (2σ, placebo) | placebo systematics |
| WHAM Hα/[O III] < 0.10 | [O III]/Hα ≳ 3 | wrong velocity tested |
| F03 Hα < 0.43 R | SHASSA < 0.77 R south; F03 north < 1.74 R | F03 is WHAM-only on 76 % of the spine |
| HI wall 0.86 ± 0.15° outside | no HI structure | finder artefact |
| X-ray "outer edge" W/E 1.32, 2.4σ | W/E 1.26, generic to the shell (rotation p = 0.37); ~2σ south step | parallel-shift nulls |
| GALEX FUV "weak excess at most" | detection, 92 ± 10 CU per R (south) | same numbers, proper nulls |
| HD 61717 A2 single star; Ca II/Na I 22 ± 2 (shocked) | SB2; 23 (15–49); cool, Ti-poor, unrelated | multi-line fits, APOGEE, Gaia RUWE |
| J0751+1807 PX 1.6 ± 0.8 mas, no RM | 0.85 ± 0.04 mas (1.18 kpc), RM +41.6 | stale catalogue |
| Paper [S II] n_e 190–440 (wave-1 eq. 5) | 225–631 (5-level) | eq. 5 biased low |
| [O III] V_helio +12.5…+14.5 (provisional) | +13.3 ± 1.0 | combined with paper Fig 3a |
| LAMOST "2\|v_corr\| = 16.2 km/s slip" hypothesis | missing v_corr (8.3) + 1.8 km/s sky offset = 10.1 | Fig 3a vectors |

---

## 16. Files (tag `canon`)

- `scripts/w3_canon_build.py`: builds the JSON, with 70 provenance checks against the source JSON/CSV. Reads only project files.
- `scripts/w3_canon_fig.py`: the figure.
- `data/derived/w3_canon_values.json`: entries, paper claims, precovery table, conflicts and inconsistencies.
- `figures/w3_canon_ladder.png` (200 dpi): (left) [O III] brightness by dataset with the adopted south (1.4–2.5 R) and north (0.4–0.9 R) bands and the superseded wave-1 values marked ×; (right) [O III] V_helio by dataset with the adopted +13.3 ± 1.0 km/s band and a V_LSR (P1) axis.
- **Inputs read:**
  - all `notes/*.md` (wave 1 and 2) and `data/derived/wave{1,2}_agent_results.json`;
  - for checks: `silk_trace.json`, `verify_velocity_summary.json`, `verify_velocity_desi_sdss.json`, `verify_desi_summary.json`, `verify_desi_audit_extra.json`, `verify_desi_audit_vel.json`, `verify_desi_velocity.json`, `verify_sdss_results_base_r0.15.json`, `verify_sdss_precovery.csv`, `verify_distance_synthesis.json`, `verify_precovery_table.json`, `verify_shockmodel_summary.json`, `verify_hiwham_wham_summary.json`, `spectra_ism_rm_results.json`.

**Limitations**
- This is a reconciliation, not a re-measurement. Values not re-verified in wave 2 are marked: DESI line width, NSNS brightness, radio, RM, gamma, dust, COS limits.
- The "(derived here)" items are simple arithmetic on quoted inputs: LSR corrections, sizes, tension significances and alternative weighted means.
- The `oii` and `northsouth` results were pending when this was written. C2, C5 and C6 should be updated from `notes/w3_oii.md` and `notes/w3_northsouth.md` when those exist.
