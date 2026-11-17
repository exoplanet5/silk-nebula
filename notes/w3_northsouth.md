# Wave 3 · tag `northsouth`: the north–south difference along the Silk

Date: 2026-11-16. Velocities are heliocentric unless marked LSR. The LSR correction along the spine runs from
−12.9 (Dec 12–15) to −9.6 km/s (Dec 25.5–28.6), for the standard solar motion used by the paper and wave 2. North
minus south, the solar reflex alone therefore adds **+2.0 km/s** to any V_LSR difference. Line ratios are energy
ratios unless marked photon.

Labels:
- **measured**: new in this work, from data on disk with the wave-2-corrected methods;
- **derived**: computed from measured values plus a stated model;
- **wave-2 / literature**: taken from a named file.

No new downloads. All code is in `scripts/w3_northsouth_*.py`, and every script that reads data was run with
`python -I`.

## 1. Bottom line

1. **The north is fainter in [O III], richer in [O II], and spectrally "cooler"** (measured, DESI and SDSS agree).
   - [O III], south (trace Dec < 19.5) vs north:
     - DESI: 1.61 ± 0.03 vs 0.73 ± 0.07 R, ratio 0.45;
     - SDSS: 1.41 ± 0.08 vs 0.51 ± 0.06 R, ratio 0.36.
   - [O II]/[O III] rises smoothly with Dec, not as a step:
     - SDSS: 1.0 (Dec 12–16.5), 1.9–2.2 (16.5–19.5), 3.6 (19.5–21), 6.9–8.3 (Dec > 21);
     - DESI: 1.1–1.25 (Dec < 18), 1.9 ± 0.7 (19.5–21), 4.2–6.4 (Dec > 21);
     - halves: DESI 1.17 ± 0.08 vs 4.0 ± 0.6; SDSS 1.68 vs 7.1.
   - The absolute [O II] *increases* to the north: SDSS 1.9 → 3.0 R (photon).
   - [O III]/Hβ falls:
     - south: > 18–21 (DESI, 2σ);
     - north: SDSS central value ≈ 7–8, > 2.6–2.9 at 2σ; DESI > 4.8–8.4 at 2σ.
   - **New:** the wave-2 DESI/SDSS [O II]/[O III] tension (DESI 1.11–1.16 vs SDSS 1.84) is a Dec-weighting effect. At
     Dec < 17.6 SDSS gives **1.17** (95 %: 0.91–1.47) and DESI 1.12 ± 0.08.
2. **The northern [O III] redshift is real but smaller and more estimator-dependent than "+20.9 vs +10.6".**
   - DESI north − south is **+4 to +9 km/s heliocentric**, each estimator at 2–3.5σ (table 3b). The model-free moment
     gives +5.7 ± 2.2.
   - The wave-2 +10.3 sits at the top of that range; it uses bootstrap-σ pixel weights.
   - In V_LSR: south −2.8…−0.2, north +5.7…+10.4, so **ΔV_LSR ≈ +7…+11 (± ~2.5) km/s**.
   - SDSS north − south: +6.7 ± 6.5 km/s.
   - DESI and SDSS agree for Dec ≥ 19.5 (DESI +16.0…+20.7 vs SDSS +16.5 ± 6.3) and for Dec ≥ 21.
   - **They disagree at Dec 17.6–21** (DESI +17…+23 vs SDSS +10.0 ± 2.2, about 2–3σ). That zone holds the bright W
     strand and the Hα patch.
3. **The NVSS RM step is not a filament signature.**
   - Reproduced: +35.7 ± 4.6 rad m⁻² (bootstrap; wave 1 +36.8).
   - Rigid placements of the whole trace give |step| ≥ 35.7 in **18 %** of trials (two-sided p = 0.18, z ≈ 1.2). One-sided
     p is 0.025–0.07; on foreground-subtracted RM, p = 0.07–0.09.
   - The best step location is +2.1° east of the spine, not at it.
   - Bootstrap errors underestimate the spatially correlated RM scatter about 6×: null sd is 27–29 against the
     bootstrap 4.6.
   - **NOT CONFIRMED** as related to the Silk.
4. **Geometry: one smooth arc, not two structures.**
   - A single circle of radius 21.3° about (94.5, +16.8) fits both halves: rms 0.19° south, 0.12° north. The north
     alone gives R = 22.8° (20.5–26.1), centred 1.5° from that point.
   - The south section is significantly *straighter* than a Monogem-centred circle:
     - free fit R = 67° (16–84 %: 27–106°); 0 % of bootstrap draws have R < 15°;
     - χ² is 375 with the Knies centre fixed, against 114 for a free circle;
     - a 12° circle would have a sagitta of 0.32° over this section.
   - So the arc is *tangent* to the 12.1° Monogem circle near Dec 14–16:
     - θ (from the Knies 2-SNR centre) is 12.48° at Dec 12, 11.98° at Dec 14 and 12.1° at Dec 17;
     - north of Dec 18.5 it lies permanently > 3σ outside (12.46°);
     - it crosses the 1-SNR 11.97° rim at Dec 21.2 and ends 3.7° outside, at Dec 28.
   - Neither G205.6+12.4 nor Geminga gives a physical alternative:
     - G205.6+12.4: the spine crosses its 8.82° rim at Dec 19.0 and is not concentric (north rms 1.7–1.9°);
     - Geminga: formally concentric now (17.5 ± 0.14°), but it moves by that amount in ~4 kyr (wave 1).
5. **Physics: the north is best described as a slower shock (≈ 70–100 km/s) of the same family, at equal or
   moderately lower driving pressure. It is not a later truncation of the southern 100–150 km/s shock.**
   - Later truncation fails in verify_shockmodel's code: in **0 of the 25** consistent southern models does it reach
     [O II]/[O III] ≥ 3.5 with [O III]/Hβ ≥ 4. Truncated tracks end at 1.46–1.61, and complete runs at only 2.1–3.2.
   - Slower 70–100 km/s tracks with the same b and X do reproduce the north (table 3c):
     - [O II]/[O III] 3.7–7.2 and [O III]/Hβ 4–17;
     - an [O III] ratio N/S of 0.63–0.80 at equal pressure, so the observed 0.36–0.45 needs **P_N/P_S ≈ 0.45–0.7**;
       some 80–90 km/s solutions need P_N/P_S ≈ 1–5 instead;
     - a FUV/[O III] ratio N/S of 0.25–0.63. GALEX measures 56/92 = **0.61 ± 0.15** (verify_xrayfuv data).
   - The rising NUV/[O III] (49 → 81 CU per R) fits more Balmer/2γ and low-ion emission. That is qualitative: NUV was
     not modelled.
   - Lower pressure is plausible because the north lies outside the eROSITA rim, where verify_xrayfuv finds no X-ray
     step.
6. **The Hα patch and PGCC G199.09+19.95 sit exactly on the Monogem-radius circle.**
   - θ = 11.97° (patch) and 11.94° (PGCC), against 12.10 ± 0.12° for the southern spine.
   - They lie 1.25–1.37° inside the bent [O III] spine.
   - Their Balmer-bright, [O II]/[S II]-strong spectrum and the cold clump make a **local shock–cloud interaction**
     the natural reading: the Monogem shock meets denser gas at its nominal radius.
   - It is local (about 1.5 × 0.3°). The northern [O III] properties continue 7° beyond it.
   - The clump's +6.5 km/s CNM line and the northern HI component (+8.5…+11.3 km/s, dominant north of Dec ≈ 23)
     match the northern [O III] V_LSR.
7. **The redshift (task d): a far/near-side geometry on one sphere is excluded.** Small effects suffice.
   - Model S (one sphere enclosing the spine) needs a *near-side* |v_exp| of 26–55 km/s.
     - It would put the southern spine at a limb gain of 1.2–1.6, against 10.6–11 for a tangent sheet: 7–9× fainter
       than observed.
     - It also predicts the south, not the north, to be the more shifted part.
   - Three small effects can each produce the shift:
     - a 3–10° tilt of the shock normal out of the sky plane (4–12 pc behind the tangent plane; v_s = 150 → 75 km/s);
     - a one-sided chord (5.9–21 km/s for a fully one-sided chord);
     - the pre-shock gas itself being redder in the north, as the HI is (ΔV_HI ≈ +6…+9).
   - The sky-plane bend is β ≈ 18°. A bend of that size along the line of sight would give 20–60 km/s, so **the bend
     is mainly in the plane of the sky**.

## 2. Results table: claim → verdict → value

| # | Claim to test | Verdict | Value (this work) |
|---|---|---|---|
| 1 | North (Dec ≳ 19.5) [O III] 2–3× fainter: 0.4–0.8 vs 1.4–2.5 R | **CONFIRMED** | DESI halves 0.73 ± 0.07 / 1.61 ± 0.03 R (×0.45); SDSS 0.51 ± 0.06 / 1.41 ± 0.08 (×0.36). On the E strand the drop is at Dec 18 (DESI E-bins 2.18 → 0.52 R); the W strand stays bright to Dec 20.5 (2.47 ± 0.38 R) |
| 2 | [O II]/[O III] 4–7 north vs 1.1–1.8 south | **CONFIRMED, and it is a gradient** | DESI 1.17 ± 0.08 / 4.0 ± 0.6 (halves); 1.12 / 3.1 / 5.3 (± 0.08 / 0.9 / 0.9) in < 17.6 / 17.6–21 / ≥ 21. SDSS 1.17 (0.91–1.47) / 2.46 (2.09–3.33) / 7.68 (6.4–9.3). Absolute [O II] rises ×1.6 to the north (SDSS) |
| 2b | Wave-2 tension: SDSS [O II]/[O III] 1.84 vs DESI 1.11 (south) | **RESOLVED** (Dec weighting) | Same Dec range < 17.6: SDSS 1.17, DESI 1.12. At 16.5–18 a residual difference remains: SDSS 1.94, DESI 1.04 ± 0.19 |
| 3 | [O III]/Hβ ~6–7 north vs ≳ 15–20 south | **CONFIRMED (weakly)** | South: DESI > 18.5–21 (2σ), SDSS > 8.5–10.8 photon. North: SDSS 7–8 (16–84 %: 2.9–15), > 2.6–2.9; DESI > 8.4 (halves), > 4.8 (≥ 21). Hβ intensity is not significantly different: SDSS 0.045 ± 0.039 vs 0.064 ± 0.058 R |
| 4 | DESI [O III] V north +20.9 vs south +10.6 | **REVISED**: real, but +4…+9 | 6 estimators (table 3b): north +16.0…+20.7, south +9.9…+12.2. Δ(N−S) = +4.9…+8.8 (σ 2.2–2.7, halves) and +3.7…+8.7 (≥ 21 vs < 17.6). The split scan peaks for a split at 17.6–19.5 |
| 5 | SDSS north +16.4 ± 5.9 agrees with DESI | **CONFIRMED overall; tension at 17.6–21** | SDSS ≥ 19.5: +16.5 ± 6.3 (DESI +16.0…+20.7); ≥ 21: +18.1 ± 6.5 (DESI +14.6…+20.3); 17.6–21: SDSS +10.0 ± 2.2 vs DESI +17.3…+22.9 (~2–3σ); bin 19.5–21: +4.2 ± 5.0 vs +20.0 ± 4.8 |
| 6 | Spine bends away from the Monogem radius north of Dec 17.6, up to 16° | **CONFIRMED, REFRAMED** | θ_K2: 12.10 ± 0.12 (< 17.6), 12.8 ± 0.3 (17.6–21), 14.4 ± 0.8 (≥ 21), 15.8 at Dec 28. But the south is itself straighter than the 12° circle (R_curv 67°, 27–106°). The arc (R ≈ 21.3°) is tangent to the Monogem circle near Dec 14–16. Persistent departure starts at Dec 18.5 |
| 7 | Crosses the 12° eROSITA rim near Dec 21.25 | **CONFIRMED** (1-SNR rim) | The 1-SNR circle (R = 11.97°) is crossed at Dec 21.18. The 2-SNR X-ray radius (11.0°) is exceeded everywhere |
| 8 | Hα + [S II] patch next to PGCC (CNM −0.7, +6.5) | **CONFIRMED; new geometric link** | Patch centroid θ_K2 = 11.97°, PGCC 11.94° (spine at the same Dec 13.19°), 1.25 / 1.37° W of the spine, 0.35° apart. SDSS patch [O III]/Hβ 1.8, [O II]/[O III] 11.8 (9.6–14.8), [S II]/Hα 1.06, [N II]/Hα 0.78 |
| 9 | NVSS RM step +37 rad m⁻² in the north only | **NOT CONFIRMED** as Silk-related | +35.7 ± 4.6 (boot); placement null sd 27–29, p(\|S\| ≥ obs) = 0.18 (translations; spin + translation 0.185), one-sided 0.025–0.07, residual RM p = 0.07–0.09. Step centre +2.1° E. South +5.4 (p 0.71–0.82) |
| 10 | GALEX covers only Dec < 22 | CONFIRMED (verify_xrayfuv) | GALEX north (17.6–22): FUV 56 ± 12, NUV 81 ± 8 CU per NSNS-R, against south 92 ± 10 and 49 ± 7.5 (wave-2 file, not re-measured). Dec 22–28 is untested |
| 11 | X-ray outer flank seen in the south only | CONFIRMED (verify_xrayfuv) | Local step +12 % south (~2σ), none north (+1.5 %, p ≈ 0.3). The north lies outside the 11.97° rim |
| 12 | North = same shock, later cooling stage | **DISFAVOURED** (in this model code) | 0 / 25 consistent southern models reach the northern ratios. Truncated-track ends reach [O II]/[O III] 1.46–1.61; complete runs 2.1–3.2 (Allen+08 100 km/s: 3.3). Later truncation also keeps FUV/[O III] constant, against 0.61 observed |
| 13 | North = slower shock (70–80 km/s, shockmodel) | **FAVOURED** | 70–100 km/s tracks (same b, X) reproduce the ratios. Brightness N/S 0.63–0.80 at equal P, so P_N/P_S ≈ 0.45–0.7 (some solutions ≈ 1–5); FUV/[O III] N/S 0.25–0.63 vs 0.61 ± 0.15 |
| 14 | Shock–cloud interaction at the PGCC clump | **PLAUSIBLE, local** | Patch on the Monogem radius; a Balmer-bright, complete-like shock (code [S II] is unreliable, see §4). It cannot explain the 7° of northern [O III] |
| 15 | Separate structure (G205.6+12.4 / Geminga) | **NOT SUPPORTED** | G205.6: rim crossed at Dec 19.0, north 12.8 ± 1.7° from its centre (not concentric). Geminga: unphysical concentricity. One 21° arc fits both halves with no break |
| 16 | North redshift = far-side expansion of a shell bending away from the tangent point | **Single-sphere version EXCLUDED**; small tilt or one-sided chord allowed; pre-shock gas velocity equally good | See §6 |

## 3. Measurements (task a)

### 3a. Methods

**DESI** (`w3_northsouth_desi.py`, `w3_northsouth_desi_velcheck.py`):
- Inputs are verify_desi revision-2 products: box-integral line fluxes, the GFA acceptance and the chromatic
  correction (`verify_desi_audit_boxfluxes.npz`, `verify_desi_coadd_meta.npz`).
- **New: Dec-binned templates.** Each strand sample (E, W, W2 polylines with the trace HWHM) is assigned to the bin of
  its own trace Dec. Then Z_k = T_k − Σ w_t ⟨T_k⟩ over (tile, petal), i.e. DESI's per-petal sky mean.
  - Check: Σ_k Z_k reproduces the wave-2 Z to 6×10⁻¹⁴.
  - All bins are fitted **jointly** to all 6,979 fibres, so the petal-mean leakage between bins is modelled. This
    removes the wave-2 caveat that "local-subset regression is not local".
- **Brightness and ratios:**
  - brightness is a WLS on [1, Z_1…Z_K];
  - ratios use an IV-type slope ratio with common weights;
  - as a check, the matched-filter amplitudes reproduce the audit's 1.462 / 1.617 / 0.707 R.
- **Errors:**
  - Formal WLS errors (scaled by the reduced χ²) are multiplied by placebo inflation factors from 60 rigid,
    non-overlapping placements (|r| < 0.05) fitted jointly with the true bins: [O III] 1.51, [O II] 1.49, Hβ 1.49,
    Hα 2.06, [N II] 1.88, [S II] 1.82.
  - HEALPix nside-32 and nside-64 block bootstraps are stored in the JSON. They are **unstable for 1.5° bins**, which
    span 1–3 nside-32 blocks (robust σ up to 0.6 R), so the placebo-inflated formal error is the primary one.
- **Velocity:**
  - per-pixel regression of the continuum-subtracted residual around [O III] 4960 and 5008 (verify_desi_audit_vel
    `prep`) on the zone templates;
  - doublet χ² grid at 0.25 km/s, nside-32 bootstrap;
  - six estimators: bootstrap-σ weights, formal-σ weights, uniform weights, 5008 only, 4960 only, and a model-free
    first moment.

**SDSS** (`w3_northsouth_sdss.py`):
- The verify_sdss pipeline functions are imported read-only (`region_analysis`, `oiii_velocity`):
  - total sky, per-fibre differential fits against same plate–spectrograph OFF fibres (0.45–2.0°, Hα patch excluded);
  - f = 0.84 / 0.60;
  - group bootstrap (2,000);
  - [O III] velocity from group ON−OFF stacks minus the per-group [O I] 5577 zero point.
- ON = within 0.15° of a strand, outside the patch 30 % contour.
- Ratio intervals are bootstrap-only. The verify_sdss placebo floors ([O II] 0.32 R south, 0.14 R north) are *not*
  included.

### 3b. Numbers

DESI joint Dec bins, all strands per bin (error = formal × placebo):

| Trace Dec | [O III] (R) | [O II]/[O III] | Hβ/[O III] | V_helio (km/s) | V_LSR |
|---|---|---|---|---|---|
| 12–15 | 1.88 ± 0.06 | 1.25 ± 0.13 | −0.009 ± 0.030 | 11.2 ± 1.0 | −1.7 |
| 15–16.5 | 1.50 ± 0.04 | 1.08 ± 0.11 | 0.006 ± 0.024 | 10.1 ± 0.9 | −2.3 |
| 16.5–18 | 1.81 ± 0.07 | 1.04 ± 0.19 | 0.065 ± 0.036 | 11.2 ± 1.2 | −0.8 |
| 18–19.5 | 0.62 ± 0.17 | undetermined | (0.78 ± 0.44) | 19.3 ± 6.2 | +7.8 |
| 19.5–21 | 1.00 ± 0.14 | 1.90 ± 0.69 | −0.10 ± 0.14 | 20.0 ± 4.8 | +8.8 |
| 21–22.5 | 0.54 ± 0.12 | 4.69 ± 1.54 | −0.10 ± 0.23 | 16.8 ± 6.8 | +6.0 |
| 22.5–24 / 22.5–25.5 (v) | 0.55 ± 0.19 | 4.17 ± 2.36 | −0.09 ± 0.35 | 8.9 ± 10.2 | −1.4 |
| 24–25.5 | 0.81 ± 0.22 | 6.40 ± 2.58 | 0.14 ± 0.28 | (in 22.5–25.5) | |
| 25.5–27 / 25.5–28.6 (v) | 0.86 ± 0.17 | 5.71 ± 1.59 | −0.13 ± 0.19 | 16.1 ± 6.1 | +6.5 |
| 27–28.6 | 0.55 ± 0.21 | 5.85 ± 3.44 | 0.19 ± 0.37 | (in 25.5–28.6) | |

- Velocities use coarser bins (`v_edges`).
- E-strand-only bins: 1.89, 1.69, 2.18, 0.52, 0.99, 1.07, 1.18, 0.81, 0.88, 0.56 R.
  W strand: 1.44 ± 0.05 (12–17.5) and 2.47 ± 0.38 (17.5–20.5). W2: 0.20 ± 0.12 R.
- Excluding the Hα-patch fibres (variant C) changes brightness by ≤ 0.02 R and ratios by ≲ 1σ (largest: [O II]/[O III]
  at 19.5–21, 1.90 → 2.51 ± 0.7).

SDSS (ON−OFF, f-corrected):

| Dec | [O III] (R) | [O II]/[O III] (16–84 %) | [O II] (R, photon) | V_helio |
|---|---|---|---|---|
| 12–16.5 | 1.65 ± 0.15 | 1.00 (0.86–1.09) | 1.22 | +7.6 ± 2.4 |
| 16.5–18 | 1.43 ± 0.09 | 1.94 (1.83–2.02) | 1.91 | +6.5 ± 2.6 |
| 18–19.5 | 1.32 ± 0.13 | 2.18 (2.01–2.46) | 2.48 | +13.3 ± 2.4 |
| 19.5–21 | 0.64 ± 0.08 | 3.60 (3.46–3.74) | 1.85 | +4.2 ± 5.0 |
| 21–22.5 | 0.59 ± 0.14 | 6.86 (5.60–8.45) | 3.37 | +12.9 ± 7.2 |
| 22.5–25 | 0.56 ± 0.10 | 8.29 (6.84–9.43) | 3.99 | +24.9 ± 10.2 |
| 25–28.6 | 0.35 ± 0.07 | 7.71 (7.00–9.21) | 2.27 | +17.3 ± 14.5 |

Zones (DESI estimator range / SDSS):

| Zone | DESI V_helio | DESI V_LSR | SDSS V_helio | HI on spine (V_LSR) |
|---|---|---|---|---|
| < 17.6 | +9.7…+12.2 | −2.8…−0.4 | +7.4 ± 2.5 | peak +1.8; comps +1.8…+2.7 (narrow +1…+2) |
| 17.6–21 | +17.3…+22.9 | +5.9…+11.5 | +10.0 ± 2.2 | peak +1.8; +1.4…+5.3; PGCC CNM −0.7, +6.5 |
| ≥ 21 | +14.6…+20.3 | +4.5…+10.1 | +18.1 ± 6.5 | peak +8.5; red comp. +8.5…+11.3 (FWHM 10–16) |

North − south ΔV_helio (DESI halves; bootstrap over common resamples):

| Estimator | ΔV_helio (km/s) |
|---|---|
| bootstrap-σ weights | +8.8 ± 2.4 |
| formal-σ weights | +5.4 ± 2.6 |
| uniform weights | +4.9 ± 2.7 |
| 5008 only | +7.9 ± 2.5 |
| 4960 only | +8.8 ± 4.4 |
| model-free moment | +5.7 ± 2.2 |

- The audit-style reproduction (fibre-Dec subsets, global template) gives +10.7 / +16.7 with formal weights, against
  the wave-2 +10.6 / +20.9. **The weighting, not the sample, drives the difference.**
- Split scan (formal weights): ΔV is largest for a split at 17.6 (+6.9) and drops to +3.1…+3.3 for splits at
  21.5–23. The DESI redshift is therefore concentrated at Dec ≈ 17.6–21, which is exactly where SDSS disagrees.
- A south 5008/4960 offset of 4 km/s appears in both halves (12.2 vs 8.1). It is a common systematic and cancels in ΔV.

## 4. Scenario tests (task c)

Shock tracks are read from `verify_shockmodel_tracks.npz` and nothing is re-run. Scenario L = later truncation on the
same track; scenario V = a slower shock with the same b and X at equal driving pressure (n₀ ∝ P/v²). The FUV proxy is
0.21 × C IV + O III] + He II, with the same scattering suppression north and south.

**L: later truncation.**
- 0 of 25 consistent southern models (v_s 100–150 km/s) reach [O II]/[O III] ≥ 3.5 with [O III]/Hβ ≥ 4.
- At the 10⁴ K end of the tracks [O II]/[O III] is only 1.46–1.61 ([O III]/Hβ 9–39).
- Complete code runs give 2.1–3.2; Allen+08 at 100 km/s, as transcribed by verify_shockmodel, gives 3.3.
- [O III] is 93–99 % complete at N*, so a later cut cannot make the north 2–3× fainter without a lower n₀G.
- C IV forms near the front, so FUV/[O III] would stay unchanged.

**V: slower shock** (from `w3_northsouth_shock.json`):

| v_N | Matches | [O II]/[O III] | [O III]/Hβ | [O III] N/S at equal P | P_N/P_S needed | FUV/[O III] N/S |
|---|---|---|---|---|---|---|
| 70 | 6 | 4.2–4.5 | 12.8–15.8 | 0.63–0.79 | 0.45–0.72 | 0.35–0.63 |
| 80 | 8 | 3.6–7.2 | 5.9–17 | 0.09–0.80 | 0.45–5.0 | 0.28–0.50 |
| 90 | 3 | 5.9–6.2 | 4.0–7.0 | 0.14–0.37 | 0.98–3.2 | 0.25–0.54 |
| 100 | 3 | 3.7 | 7.1 | 0.70–0.74 | 0.48–0.64 | 0.26–0.41 |

The observed values are [O III] N/S 0.36–0.45, FUV/[O III] N/S 0.61 ± 0.15, and [O II] N/S 1.6 (SDSS). The
near-complete slow solutions predict:
- [O II] N/S of 2.6–3.1 × P_N/P_S = 1.2–2.2 (observed 1.6);
- Hβ N/S of 2.0–3.1 × P = 0.9–2.2 (observed ~1.4, poorly measured).

**Reading.** The north is the same incomplete-shock family, but slower: about 70–100 km/s against 100–150 km/s south.
It is driven at equal or somewhat lower pressure, with pre-shock density of order equal or up to ~1.5× higher.
- This fits the north lying outside the X-ray rim, beyond Dec 21.2, with no X-ray step there.
- It fits the larger O⁺ fraction, more Balmer emission, and the drop in FUV/[O III].
- The NUV/[O III] rise (×1.66) is not modelled. It is consistent in direction with more H two-photon and C II]/C III]
  emission.
- Within the code's systematics (§7), the FUV prediction sits at the low end of the observed 0.61 ± 0.15.

**Hα patch / PGCC** (task c, shock–cloud):
- The patch spectrum (SDSS; DESI gives a similar [N II]/Hα 0.65 and [S II]/Hα 0.85) is Balmer-bright:
  - [O III]/Hβ 1.8;
  - Hα/Hβ 2.84, unreddened case B;
  - [O II]/Hβ ≈ 13;
  - [S II]/Hα 1.06 and [N II]/Hα 0.78.
- The closest code spectra are 10⁴ K ends of partly pre-ionised 80–140 km/s runs:
  - [O III]/Hβ 1.8–3.5, [N II]/Hα 0.74–1.0, Hα/Hβ 2.71–2.75;
  - but [O II]/Hβ 20–31 against 13, and [S II]/Hα 0.5–0.75 against 1.06.
- The code's [S II] is not trustworthy: it has no S⁺/S²⁺ dielectronic recombination and an 8,000 K floor, so its
  complete runs give 0.33–0.35 where published complete shocks give 0.5–1.2.
- Qualitatively, this is a radiative shock in denser gas that has completed or nearly completed its cooling.
- **Geometry:** the patch (θ_K2 = 11.97°) and the PGCC (11.94°) lie exactly on the circle traced by the southern
  spine (12.10 ± 0.12°). The [O III] spine at the same Dec is at 13.19°.
- **Reading:** where the Monogem shock reaches the PGCC complex at its nominal radius, it has become a slow,
  Balmer-bright cloud shock. North of Dec 18.5 the faint [O III] sheet lies 1–4° further out.
- The patch is local (1.45° × 0.31°). The northern [O III] change persists 7° beyond it, so the clump is a marker of
  denser and redder ambient gas, not the cause of the whole northern change.

**Separate structure:**
- **G205.6+12.4** (Knies+24 second remnant, R = 8.82°):
  - the south spine lies 6.5° from its centre, inside its rim;
  - the spine crosses the rim at Dec 19.0 and is not concentric in the north (12.8 ± 1.7°);
  - a G205.6 limb would need θ ≈ 8.8°, which holds only near Dec 19.
- **Geminga (present position):** north θ = 17.5 ± 0.14°, but Geminga moves by more than 0.14° in a few kyr (wave 1),
  so this is chance.
- **One arc:** a single circle of R = 21.3° about (94.5, +16.8) fits both halves (rms 0.19° / 0.12°). The north alone
  gives a centre 1.5° from that point, with R = 22.8° (20.5–26.1). There is no geometric break, so the "north" is the
  continuation of one filament.

## 5. RM placement null (task b)

`w3_northsouth_rm.py`; input `spectra_ism_rm_sources.csv` (Taylor+09). My spine coordinates agree with wave 1: median
|Δ offset| 0.0003°, signs agree for 96 %.

| Half | S_true (rad m⁻²) | Null sd (translate / spin) | p two-sided | p one-sided | p (residual RM) |
|---|---|---|---|---|---|
| North (trace Dec 20–28) | +35.7 ± 4.6 | 29.2 / 27.1 | 0.18 / 0.185 | 0.025 / 0.07 | 0.07 / 0.09 |
| South (12–20) | +5.4 ± 7.1 | 21.1 / 26.9 | 0.71 / 0.82 | 0.26 / 0.35 | 0.30 / 0.46 |

- 600–634 accepted placements per case, of which 367–527 do not overlap the true side bands. The non-overlapping
  subset gives the same p (0.18).
- The L1 step location in the north is +2.1° east of the spine. Placements give L1 steps at least as large 24–27 % of
  the time.
- A smooth planar or quadratic gradient fitted at ≥ 7° predicts −7.4 for the step, so the step is a local RM feature.
  Such features are common at this angular scale.
- **Conclusion:** the "north-only RM step" is not specific to the filament.

## 6. Shell-projection models for the north redshift (task d)

`w3_northsouth_geom.py`; D_c = 288 pc; sheet FWHM 0.28°. The observed ΔV_LSR(N−S) is +7.0…+10.8 km/s: the
estimator range, plus the +2.0 km/s solar reflex, with σ ≈ 2.5.

- **S. One sphere about a Monogem centre, enclosing the whole spine.**
  - θ_R = 16.1–17.5°, so R = 80–86 pc.
  - The far/near factor f = √(1 − (sin θ / sin θ_R)²) is 0.65–0.80 in the south and 0.44–0.53 in the north.
  - So the *south* is the more displaced part, and the observed sign requires the near side with |v_exp| = 26–55 km/s.
    The systemic velocity would then be ≈ +19…+35 km/s LSR (south V_LSR + |v_exp| f_S), far from all local gas.
  - Per-side limb gain: 1.2–1.6 in the south against 10.6–11.0 for a tangent sheet. The sharp southern ridge would be
    7–9× fainter than observed, and double-peaked if both sides emit. **Excluded.**
- **T. Locally tangent sheet whose normal is tilted by i out of the sky plane** (v = v_pre + u sin i, u = 0.85 v_s):

  | v_s (km/s) | i needed | Depth offset behind the tangent plane (r_min = 70 pc) |
  |---|---|---|
  | 75 | 6.3–9.8° | 7.7–12.1 pc |
  | 100 | 4.7–7.3° | 5.8–9.0 pc |
  | 125 | 3.8–5.8° | 4.6–7.2 pc |
  | 150 | 3.1–4.9° | 3.9–6.0 pc |

  The sky-plane bend relative to the azimuth about the 2-SNR centre is β ≈ 18° in the north. A comparable
  line-of-sight bend would give 0.31u = 20–40 km/s, far above the observed shift. So any far-side component is small,
  and the bend lies mostly in the plane of the sky.
- **CH. Tangent curved sheet with a one-sided chord** (the near half missing):
  - centroid 2/3 v_max = 5.9–20.6 km/s for R_c = 58–79 pc, emitting depth 0.75–1.7 pc and v_s = 75–150;
  - a far/near flux asymmetry of ≈ 0.4–1 reproduces the observed shift.
- **Pre-shock gas.** verify_shockmodel's tangent geometry puts the centroid at the pre-shock gas velocity.
  - The on-spine HI shifts redward in the north: a +8.5…+11.3 km/s component is dominant north of Dec ≈ 23, and the
    PGCC CNM sits at +6.5. ΔV_HI ≈ +6…+9 against ΔV_[O III] ≈ +7…+11.
  - [O III] − HI (relevant component) is about −2…−4 km/s in *both* halves. That is a common zero-point offset: the
    DESI south here is 1–3 km/s below the wave-2 best +13.3.
  - **This explains the redshift without any far-side geometry.** It is weakest at Dec 19.5–22.5, where the red HI
    component is not yet dominant and only the PGCC CNM line is at +6.5.

**Verdict (d).** A shell that bends away from the tangent point *on one sphere* cannot give a northern redshift
together with a sharp southern limb. A few-degree normal tilt or a one-sided chord can. Redder pre-shock gas, seen
directly in HI, is the most economical explanation.

**Discriminating test.** Spectra at R ≳ 20,000 (σ ≲ 6 km/s) of the north:
- pre-shock velocity: a symmetric profile centred on the HI component;
- tilt or one-sided chord: an asymmetric profile with a sharp edge at the pre-shock velocity and a red tail to
  v_max ≈ 10–25 km/s.

## 7. Limitations

- **DESI fine bins:**
  - The 18–19.5 bin is degenerate for ratios: the [O III] slope under [O II] weights is near zero, and the
    null-[O II] box gives 6.6 there.
  - Bins with < 1 R in the north carry 15–40 % errors.
  - Placebo inflation is one factor per line, applied to every bin.
  - Template shape is fixed: the trace HWHM, with Gaussian cross-sections.
- **DESI velocities** depend on the pixel weighting at the ±2–4 km/s level in the north. I quote the estimator range,
  not a single value. The absolute zero point carries ±1.5 km/s.
- **SDSS:**
  - Bins hold 6–16 groups and only 7 BOSS groups at Dec < 17.6, for [O II].
  - Ratio intervals omit the placebo floors.
  - Velocities in the north have ±5–15 km/s errors.
  - [O II] 3729/3726 is 1.05–1.08 (16–84 %: 0.91–1.25) north against 1.29–1.31 (1.06–1.67) south. That is suggestive
    of higher n_e north but not significant, and not claimed.
- **Shock-track inferences** inherit every verify_shockmodel limitation:
  - no radiative transfer;
  - Hβ uncertain ×2;
  - [S II] unreliable;
  - truncated runs end at 10⁴ K;
  - FUV and NUV not modelled beyond C IV / O III] / He II;
  - b gridded.
  The "slower shock" range (70–100 km/s) and the pressure ratio are indicative.
- **GALEX** north/south ratios come from verify_xrayfuv: NSNS-scaled, north only to Dec 22, joint-fit significance in
  the north only z = 2.7. I did not re-measure them.
- **HI comparison:** the HI used is the on-spine spectrum, which includes foreground and background gas along the
  line of sight. Which HI component is "pre-shock" is not determined.
- **Geometry:** spine σ per point is 0.005–0.02°, but the arc wiggles, so the circle χ² values are dominated by
  structure (χ²/dof ≈ 5). Curvature comparisons use bootstrap ranges and χ² differences, not formal probabilities.
- **Shell models:** a single D_c = 288 pc. Model T and model CH are kinematic sketches, not hydrodynamics.

## 8. New findings

1. The DESI–SDSS [O II]/[O III] tension is a Dec-weighting effect: at Dec < 17.6, SDSS 1.17 against DESI 1.12.
2. [O II]/[O III] rises smoothly along the arc, while absolute [O II] rises and [O III] falls. The north is a slower
   shock (70–100 km/s) rather than a later truncation of the southern one.
3. The DESI northern redshift is +4…+9 km/s depending on the estimator. It is concentrated at Dec 17.6–21, where SDSS
   disagrees (2–3σ).
4. The NVSS RM step is not filament-specific: placement p = 0.18.
5. The southern spine is tangent to, not concentric with, the Monogem 12° circle.
   - The Silk is one R ≈ 21° arc.
   - This weakens slightly the "on the Monogem limb" registration used for the 252–318 pc distance: only the mean
     radius matches, near Dec 14–16. It does not refute it.
6. The Hα patch and the PGCC clump lie exactly on the Monogem-radius circle (11.94–11.97°), 1.25–1.37° inside the
   [O III] spine. That points to a local shock–cloud interaction at the nominal shock radius.

## 9. Files

**Scripts** (run in this order with `~/.venvs/astro313/bin/python -I`):
- `scripts/w3_northsouth_desi.py`
- `scripts/w3_northsouth_desi_velcheck.py`
- `scripts/w3_northsouth_sdss.py`
- `scripts/w3_northsouth_rm.py`
- `scripts/w3_northsouth_geom.py`
- `scripts/w3_northsouth_shock.py`
- `scripts/w3_northsouth_figs.py`
- `scripts/w3_northsouth_summary.py`

**Data:**
- `data/derived/w3_northsouth_desi.json`, `w3_northsouth_desi_spec.npz`
- `w3_northsouth_desi_velcheck.json`
- `w3_northsouth_sdss.json`
- `w3_northsouth_rm.json`, `w3_northsouth_rm_nulls.npz`
- `w3_northsouth_geom.json`
- `w3_northsouth_shock.json`
- `w3_northsouth_summary.json`, which collates all key numbers

**Figures** (dpi 200):
- `figures/w3_northsouth_profiles.png`: SB, [O II]/[O III], [O III]/Hβ, V_LSR with HI, θ − R vs Dec
- `figures/w3_northsouth_velocity.png`: regression spectra, estimator spread, shell models
- `figures/w3_northsouth_rm.png`
- `figures/w3_northsouth_shock.png`: sky geometry, ratio tracks, slower-shock predictions
