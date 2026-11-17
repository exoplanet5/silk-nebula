# Wave-4 `family`: Monogem-family test with other periphery [O III] filaments

Date 2026-11-17. Tag `family`. Velocities heliocentric (DESI barycentric, < 0.02 km/s different) unless marked V_LSR
(standard solar motion; V_LSR − V_helio = −12.1 km/s at E1, −12.4 at the Silk south, derived here with astropy).
Intensities in rayleigh (R, photon). Ratios are **energy** ratios to [O III] 5007 unless stated.
Labels: **measured** (my code on public data), **derived**, **literature** (file named).

## 1. Verdict

**The family test is mostly blocked by coverage. Its one measurable case is related to the Silk spectrally but not
kinematically.**

1. **The true Monogem south limb (S1) cannot be measured.**
   * S1 is the bright, 11°-long [O III] arc of Fesen+26 Fig 32 (v2 numbering). I traced it in NSNS at
     RA 102.1–113.4, Dec 3.4–5.0. It lies 11.85° (rms 1.0) from the K24 two-SNR centre and 10.2° from the
     single-SNR centre, so it sits on the Monogem circle like the Silk's south spine (12.10°).
   * Its NSNS brightness is like the Silk's: stacked ridge peak 2.0 R (Gaussian fit 1.66 R), FWHM 0.10°.
   * The only spectra are DESI DR1 backup-program coadds of bright stars (1,524 coadds, 76 on the ridge). They have a
     per-fibre σ of 2.5 R. The DESI [O III] peak is 0.40 ± 0.27 R formal, but shifted-template nulls scatter by
     1.46 R rms and fake lines by ±0.9 R. The data therefore exclude nothing.
   * No SDSS sky fibre reaches S1. Plate 2846 stops at Dec 3.7; its fibres on disk start at Dec 5.0, and the ones
     near the ridge were never downloaded.
   * **S1: not testable with DESI DR1 or SDSS.**
2. **The NW limb group and the control G206.6+6.1 have no DESI DR1 or SDSS spectra.**
   * NW group: RA 97–103, Dec 13–21. The nearest DESI coadds are at Dec > 24. The nearest SDSS plates, 3334/3335 and
     2678/2696, miss the filaments.
   * G206.6+6.1: the nearest DESI coadd is 3.1° away.
   * For G206.6+6.1 the only spectrum is Fesen+26's own slit spectra (literature; peak heights read by eye from
     their Fig 28, v2):
     - P2, an Hα-selected filament: [O III] 5007 ≈ 9 and Hβ ≈ 1.5 (10⁻¹⁷ cgs), so [O III]/Hβ ≈ 6, with Hβ detected;
     - P1: no [O III] at all.
   * The control is therefore Balmer-visible where measured, but those slits were placed on Hα filaments.
3. **E1 has a firm measurement.** E1 is Fesen's "shorter and much fainter [O III] filaments a bit farther east"
   of the 15° group. The DESI main survey (1,511 coadds) gives a firm measurement:
   * [O III] peak **0.41 ± 0.06 R**; 0.25 R if a free intercept is allowed.
   * **[O II]3727/[O III] = 4.5 (68 %: 3.9–4.8; 95 %: 3.5–5.3)**.
   * [O III]/Hβ **> 6.9 (2σ)**.
   * He II/[O III] < 0.04 (97.5 %).
   * [N II]6583/[O III] 0.37 (0.28–0.45) and [S II]/[O III] 0.46.
   * Hα/[O III] 0.28 (0.12–0.40), a 1.7σ Hα detection.
   * **V_helio −12.1 ± 2.0 (stat) ± 1.5 (zp) km/s, so V_LSR −24 ± 2.5 km/s.**
4. **Comparison with the Silk, same pipeline.** My own run on the Silk south reproduces the canon:
   * [O III] 1.64 ± 0.13 R (canon 1.62);
   * V_helio +10.2 ± 1.5 (canon +10.6, 10.2–12.0);
   * [O II]/[O III] 1.09 (1.00–1.11) (canon 1.11–1.16).

   Against that check, E1:
   * is Balmer-faint with weak He II, consistent with the truncated/slow-shock class but only at > 6.9 (the Silk
     south gives > 11 to > 22);
   * has **the [O II]/[O III] of the Silk *north* (4–8), not the Silk south (1.1)**;
   * is **22 ± 2.5 km/s bluer than the Silk south** (ΔV_helio −22.3 ± 2.5; the zero point cancels because the tiles
     and pipeline are the same) and ~35 km/s bluer than the Silk north. E1 is not at the local-gas velocity
     (HI/WIM ≈ 0 to +5 V_LSR).
5. **H1 vs H2: no decisive change.** For the record (judges own the numbers), I would hold P(H1) and P(H2)
   unchanged within ±0.03. The reasons:
   * The measurable family member (E1) shows the same broad excitation class (Balmer-faint, [O II]-strong, He II
     weak). Every shock hypothesis predicts that class, so it does not discriminate.
   * Its velocity shows that the [O III] filaments east of Monogem do not share one tangential-limb velocity. The
     simplest "all limbs at V_LSR ≈ 0" picture is therefore not borne out for E1.
   * E1 lies 14.8° from the two-SNR centre, i.e. 3.8° outside the K2 radius, but on the single-SNR circle (12.15°,
     rms 0.24, against R = 11.97°).
   * The key H1 test, the S1 south-limb spectrum on the same circle as the Silk south, remains undone.

## 2. Inventory of Monogem-periphery [O III] filaments (step 1)

Source: Fesen+2026 §6.2.1 text (`data/raw/context_arxiv/fesen2026_2608.26666v1_flow.txt`) and Figs 31–32 of the v2
PDF (`data/raw/discovery_arxiv/fesen2026_2608.26666v2.pdf`, pages 36–37). The positions are my NSNS ridge traces
(`scripts/w4_family_trace.py`; NSNS DR0.2 order-3 HiPS already on disk; 1′/px; ridge S/N ≥ 4).

| ID | Fesen+26 description | NSNS trace (measured) | Distance to centres (deg): K2 / single / B0656 | NSNS peak, FWHM | Spectroscopy |
|---|---|---|---|---|---|
| Silk | "nearly vertical 15° group", E/N limb | `silk_trace.json` | 12.10 (south) / 9–16 / 11.36 | ~2 R S; 0.4–0.9 N | DESI main, SDSS (waves 1–3) |
| **S1** | S limb, "line of filaments just outside (south of) the X-ray southernmost emission" (Fig 32) | 152 points, RA 102.1–113.4, Dec 3.4–5.0; bright arc 102.3–108.5, separate piece 112–113.3 at Dec 4.0–4.2 | 11.85 ± 1.02 / 10.17 ± 0.36 / 10.72 | stack 2.0 R (fit 1.66), 0.099° | DESI **backup only** (stars); no SDSS |
| **E1** | "a bit farther east … shorter and much fainter [O III] filaments" | 65 points, RA 119.2–119.9, Dec 13.3–16.9, two close strands | 14.84 ± 0.15 / **12.15 ± 0.24** / 14.12 | stack 1.1 R (fit 0.93), 0.105° | DESI main dark/bright; SDSS/BOSS sky fibres |
| NW group | "bright and highly curved [O III] filaments" at the NW X-ray rim / Gemini Hα ring | visible in NSNS at RA ~100–103, Dec 13–20 (not traced) | ~11–12 | bright | **none** (DESI only Dec > 24; no SDSS plate) |
| G206.6+6.1 (control) | 3.6° [O III]-dominated SNR toward the Monogem centre, ~440 pc (Gao+22 via Fesen) | rim traced radially about (105.6, +8.3): 84 points, r_med 1.5° | interior | smoothed rim 2.3 R | **none**; Fesen slits only (literature) |
| G191.4+11.1, G190.4+12.5 NE rim | faint [O III], N limb | not traced | — | faint | DESI backup partial at Dec > 24 (not attempted) |

`figures/w4_family_overview.png` shows the full NSNS mosaic, the DESI DR1 tile and SDSS plate footprints, and the
Monogem circles. `figures/w4_family_traces.png` shows the traces and the aligned cross-profiles.

## 3. Measurements (step 2)

### 3.1 DESI DR1 (`scripts/w4_family_desi.py` → `data/derived/w4_family_desi.json`)

**Data.**
* 3,103 new SPARCL coadds in 2 passes: 1,603 south and 1,500 E1; 570 MB; ~5 min, because SPARCL ran at ~10 spectra/s
  this time.
* Log: `data/derived/w4_family_windows/fam_log.jsonl`.
* Queues: `w4_family_sel_queue{,_ext}.csv`.
* Also used: 45 E1 coadds and 2,127 Silk-south coadds already on disk.
* TAP queries for zpix/fiberassign/exposure of the 7 south tiles are in `data/raw/w4_family_tap/`, with the query
  text and date.

**Method.** The corrected wave-2 recipe, with code imported unchanged from `verify_desi_common.py` and
`verify_desi_linefit.py`:
* residual = flux − Redrock model;
* simultaneous quadratic + Gaussian fits at a fixed velocity: +13 km/s for S1 and the Silk, −12 km/s for E1, set from
  its first-pass stack;
* SB = F × f_ps(λ)/Ω, where f_ps is the GFA fibre fraction (Moffat-chromatic, EFFTIME-harmonic). This is the 0.736
  correction;
* unit-peak Gaussian template on the NSNS ridge (σ 0.042° for S1, 0.045° for E1; the Silk uses the wave-2 E/W/W2
  template), minus the tile-petal mean;
* weighted slope with an |χ| < 8 clip;
* bootstrap over the dominant tile-petal (46 groups for E1).

Stars get no object-line mask: Redrock z ≈ 0 would otherwise mask [O III] itself.

| Quantity | E1 | Silk south (my check) | S1 |
|---|---|---|---|
| n coadds (on ridge, T > 0.5) | 1,511 (1,139) | 2,063 (640) | 1,524 (76) |
| per-fibre σ([O III]) | 0.27 R | — | 2.5 R |
| [O III] ridge peak | **0.41 ± 0.06 R** (formal 0.010; bootstrap/formal ≈ 6). With free intercept 0.25 R (+0.11 R offset) | **1.64 ± 0.13 R** | 0.40 ± 0.27 formal; null rms 1.46 R → not constraining |
| [O III] by segment | Dec 13.2–14.5: 0.39; 14.5–15.5: 0.58; 15.5–16.9: 0.29 R (formal ±0.02) | — | RA 108–110.5: −0.06 ± 0.44; 110.5–113.5: 0.68 ± 0.35 |
| V_helio ([O III] stack) | **−12.1 ± 2.0** (boot; 16–84 %: −13.9…−10.2) | **+10.2 ± 1.5** | unmeasurable (±37) |
| V_LSR | **−24.2 ± 2.5** (incl. ±1.5 zp) | −2.1 ± 2.1 | — |
| [O II]3727/[O III] | **4.47 (3.89–4.84; 95 %: 3.52–5.30)**; two-component fit 4.31; intercept variant ~3.9 | **1.09 (1.00–1.11)** | — |
| [O II]3729/3726 (photon) | 1.58 ± ~0.5 (0.89 ± 0.18 / 0.57 ± 0.13 R; consistent with the low-density limit, uninformative) | — | — |
| Hβ/[O III] | 0.021 (−0.045…0.068) → **[O III]/Hβ > 6.9 (2σ)** | 0.027 (0.017–0.034) → > 22.7 | — |
| He II/[O III] | −0.046 (−0.086…−0.008); 97.5 % < 0.038 | — | — |
| [O III]4363/5007 | 0.004 ± 0.066 (non-constraining) | — | — |
| Hα/[O III] | 0.28 (0.12–0.40; 95 % −0.06…0.52), 1.7σ | — | — |
| [N II]6583/[O III] | 0.37 (0.28–0.45) | — | — |
| [S II]6716+6731/[O III] | 0.46 (6716 0.32, 6731 0.14) | — | — |
| [Ne III]/[O III] | 0.03 (−0.04…0.11) | — | — |
| Nulls | Fake lines 4885/4990/5030: −0.034, −0.012, −0.007 R (±0.01). [O I]5577 −0.19 ± 0.13 R | fake lines −0.014, −0.003 R | fake lines ±0.3 R formal, ±0.9 R robust ON−OFF |

**Notes on the table.**
* **[O II] normalisation trap (new, important for anyone re-using the code).** The `o2sum` group of
  `verify_desi_linefit.py` is one profile, 1.0·G(3726) + 1.3·G(3729). Its amplitude is therefore the 3726 flux only,
  and the doublet total is 2.3 × amp. This is the same trap as the [O II]7325 bug that `w3_oii` found.
  * I apply the factor 2.3 in `w4_family_desi.py`. With it, the Silk check gives 1.09, in agreement with the canon
    DESI 1.11–1.16.
  * So `verify_desi` evidently handled it, but wave-4 users of the raw `o2sum` amplitude must not forget it.
  * The independent two-component fit agrees: 4.31 against 4.47 for E1.
* **E1 brightness.** DESI E1/NSNS-E1 ≈ 0.4–0.45 (0.41 R against a 0.93–1.1 R NSNS peak). The Silk gives ≈ 0.8.
  * The intercept variant (0.25 R) shows that the E1 profile in DESI is broader than the NSNS-based template, or sits
    on a pedestal.
  * The E1 peak is therefore **0.25–0.41 R**, about 4–6× fainter than the Silk south. The ratios are less sensitive:
    [O II]/[O III] is ~3.9 with the intercept.
  * All E1 fibres lie within 0.12° of the ridge (the 3,000-spectrum budget), so off-filament leverage comes only
    from the petal-mean model. A wider E1 control sample would pin the zero point.
* **E1 velocity robustness.**
  * Same tiles, survey and pipeline as the Silk check (+10.2), so the −22 km/s difference does not depend on the
    DESI zero point.
  * Tiles 41485/41487 are not involved.
  * The fit is a free-velocity doublet fit to the Z-weighted regression stack. The [O III] doublet in the stack is
    clean (`figures/w4_family_compare.png`, panel b).

### 3.2 SDSS/BOSS sky fibres, E1 (`scripts/w4_family_sdss_e1.py` → `data/derived/w4_family_sdss_e1.json`)

**Method.** Wave-2 per-fibre TOTAL-sky differential fits (`verify_sdss_diff_base.npz`, which is linear, so a new OFF
set is exact) and the `verify_sdss_results.py` combination:
* ON: d_E1 < r_on;
* OFF: 0.45–2° from E1 and > 0.45° from every Silk strand;
* f = 0.84 / 0.60.

**Results.**

| | r_on 0.10° | r_on 0.15° |
|---|---|---|
| groups / ON fibres | 7 / 21 | 9 / 38 |
| [O III] ON−OFF | 0.37 ± 0.15 R (placebo rms 0.12) | 0.17 ± 0.08 (placebo rms 0.03) |
| Hβ | −0.11 ± 0.11 | −0.16 ± 0.07 |
| [O II] (BOSS only, 4 groups) | 1.64 R → ratio ~7 (95 %: 0.3–88) | ~9 (unconstrained) |

* Velocity: unmeasurable (±33–55 km/s).
* So SDSS shows E1 at ~2.4σ, at a level consistent with DESI (a 0.41 R peak averaged over the 0.10° ON radius gives
  ≈ 0.25 R). It adds nothing to the ratios.

## 4. Comparison and interpretation (step 3)

| | Silk south (canon) | Silk north (canon) | E1 (this note) | S1 | G206.6+6.1 (control) |
|---|---|---|---|---|---|
| Position relative to Monogem | on K2 circle (12.10°) | 1–4° outside | 3.8° outside K2; on the single-SNR circle | on K2 circle (11.85°) | interior (projected) |
| [O III] peak | 1.4–2.5 R | 0.4–0.9 R | 0.25–0.41 R | NSNS ~2 R; spectra uninformative | NSNS rim ~2.3 R |
| [O III]/Hβ | > 11 (DESI 2σ), ~20 SDSS | ~6 (SDSS; > 2.7) | > 6.9 | — | ≈ 6 (Fesen P2 Hα-selected slit, by eye) |
| [O II]/[O III] | 1.1–1.2 | 4–8 | **4.5 (3.9–4.8)** | — | — |
| He II/[O III] | < 0.026 | — | < 0.04 | — | — |
| V_LSR | +1.2 ± 1.0 (DESI S −1.7) | +10.6 (8–13) | **−24 ± 2.5** | — | — |

**Interpretation (derived; model-free).**
* **Same class?** E1 is Balmer-faint ([O III]/Hβ > 6.9), He II-weak and [O II]-strong. Its low-ionisation lines are
  about 0.3–0.5 of [O III]. That is a slow-shock spectrum, comparable to the Silk north (lower-velocity or
  more-truncated shock in the shock-model language of `verify_shockmodel`), not to the Silk south.
  * In the truncated-shock family [O II]/[O III] rises as the shock slows. E1 would then be a slower shock,
    consistent with its 4–6× lower surface brightness.
  * This is a family resemblance in excitation class, which every shock driver predicts.
* **Velocity.** Under H1 the Silk south (tangent to the K2 circle) sits at the local-gas velocity, as a tangential
  limb should.
  * E1 at V_LSR −24 is not at the local-gas velocity. If E1 is part of the Monogem shock, its shock normal must have
    a ≥ 24 km/s line-of-sight component toward us. Pre-shock gas moving at −24 km/s is the other possibility, but
    `verify_hiwham` found the local HI here at V_LSR +0…+11.
  * A ~100 km/s shock needs an inclination of ≳ 14° from the plane of the sky. That is easy for a filament projected
    3.8° outside the K2 rim (not a tangent point of a K2 sphere). It is natural on the single-SNR circle only if
    that shell is not seen exactly edge-on there.
  * The velocity therefore does not falsify H1. It does show that E1 and the Silk are not one coherent tangential
    sheet, and it argues against the stronger H1 prediction that every rim filament is at V_LSR ≈ 0.
* **What would have discriminated.**
  * S1 is the decisive analogue: same circle as the Silk south, comparable NSNS brightness, and a genuine southern
    tangent point. H1 predicts there: [O III]/Hβ ≳ 10, [O II]/[O III] ~1–2, V_LSR ≈ 0 ± 5 (plus ~−16 km/s
    V_LSR − V_helio, so V_helio ≈ +16), and FUV/[O III] ~70–90 CU/R.
  * A different driver (H2) predicts no such requirement.
  * The public DESI DR1 backup data are 6–10× too noisy, and no SDSS sky fibres exist there. **Pre-registered
    prediction for DESI DR2 (main survey does not reach b ≈ 2–11°, so likely also backup-only) or a dedicated
    long-slit/IFU spectrum on S1 at (106.0, +3.47) or (112.5, +4.05).**

## 5. Caveats

* **S1 null.** The DESI S1 "0.40 ± 0.27 R" must not be read as a detection or as a conflict with NSNS. The null
  rms is 1.46 R, chi std 2.3, and stellar-template residuals dominate.
* **NSNS calibration** is not re-verified (`imaging_optical.md`). The NSNS peaks are quoted only to compare
  filaments with each other.
* **E1 template width and intercept.** The amplitude depends on the model at the ×1.6 level (0.25–0.41 R). The
  ratios depend on it at ~±15 %.
* **No placebo for E1 in DESI.** Shifted templates find no fibres, because all E1 fibres lie within 0.12°. The
  fake-line nulls (≤ 0.03 R, 3× below the bootstrap error) and the Silk pipeline check are the controls.
* **The G206.6+6.1 numbers are read by eye** from a published low-dispersion plot, at Hα-selected slit positions.
  They are not comparable in selection to fibre averages along an [O III] ridge.
* **E1's identification** as Fesen's "farther east" filaments is my reading of the text and Fig 31. Fesen give no
  coordinates.

## 6. Files

**Scripts.**
* `scripts/w4_family_map.py`: overview mosaic and coverage.
* `scripts/w4_family_trace.py`: NSNS traces.
* `scripts/w4_family_tap.py`: Data Lab TAP.
* `scripts/w4_family_select.py`: SPARCL queue.
* `scripts/w4_family_fetch.py`: SPARCL fetch, a wrapper of the wave-1 fetcher.
* `scripts/w4_family_desi.py`: DESI regression.
* `scripts/w4_family_sdss_e1.py`: SDSS E1.
* `scripts/w4_family_figs.py`: comparison figure and summary.

**Data.**
* `data/derived/`:
  - `w4_family_nsns_wide.npz`, `w4_family_traces.json`;
  - `w4_family_sel_queue{,_ext}.csv`, `w4_family_sel_summary{,_ext}.json`;
  - `w4_family_windows/` (fam_*.npz, logs);
  - `w4_family_desi.json`, `w4_family_desi_fibres.npz`, `w4_family_desi_stacks.npz`;
  - `w4_family_sdss_e1.json`, `w4_family_summary.json`.
* Raw:
  - `data/raw/w4_family_tap/`: DESI tiles, zpix, fiberassign and exposure, with `.query.txt` files;
  - `data/raw/w4_family_sparcl/find_south.csv`;
  - `data/raw/w4_family_sdss/platex_monogem.csv`: SDSS DR17 plates, with query text.

**Figures.**
* `figures/w4_family_overview.png`
* `figures/w4_family_traces.png`
* `figures/w4_family_compare.png`
