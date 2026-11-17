# w4_geostats: unified geometry statistics for the Silk Nebula

Wave 4, tag `geostats`, written 2026-11-17. This note re-derives, with one defined null per question, the geometry statistics that the three judges computed inconsistently (`w3_judge_phys` §2.4, `w3_judge_skeptic` §3A–B, `w3_judge_obs` §4).

**Inputs.** Project files only, all on disk. No downloads were made, and every script ran with `python -I`.
- `silk_trace.json` v1.0
- `w3_figs_f1_maps.npz` (NSNS DR0.2 [O III], 1′ CAR)
- GALEX 1′ skybg blocks via `verify_xrayfuv_galex_measure.prep_band`
- `verify_desi_linefluxes.npz` and `verify_desi_coadd_meta.npz` (o3 and o2sum × kR, `good_fit` × `ok_contrib`: 7,158 of 7,620 fibres usable)
- `verify_xrayfuv_erass_pixels.npz` (0.2–0.6 keV = b0204 + b0406, catalogue mask applied)

**Conventions**
- x is the perpendicular offset from the E spine; + is east (outside, away from Monogem).
- Segments:
  - S1, Dec 12–15: one strand, the only clean section;
  - S12, Dec 12–17.6: includes the W strand;
  - N, Dec 20–28: W2 strand and the Hα patch; the patch is masked within 0.9°.
- "p" is an empirical placement/injection probability, (k+1)/(n+1). Its floor is about 0.02–0.03 for maps and about 0.003 for DESI.
- Gaussian z values are also given, but the nulls have heavy tails, so **the empirical p is the significance to quote**.

---

## 0. Bottom line

| Question | Judges said | Unified result here |
|---|---|---|
| **Chance that the south spine sits at the Monogem radius** | phys: Bayes factor (BF) 60–120 (10–20 after look-elsewhere). skeptic: p = 0.03–0.07, K2 only | One null: rigid trace placed at random in the same b band, with a match to **any** of K1/K2/PSR counted. **P(match \| trace reaches a rim) = 0.092** (variants 0.072–0.122). The unconditional value in the band is 0.030. **BF 3.7–7.4** (2.8–9.4 over variants); 8–16 if K2 is the only a-priori centre. The judges' numbers are reproduced with their own definitions (§1.4) |
| **Does the coincidence depend on the Knies decomposition?** | skeptic: yes, K1 gives 9.41° vs R = 11.97° | **No.** In the Silk's own PA sector the eRASS1 shell's outer half-maximum lies at the Silk in **every** centre frame. Silk − r½ = **+0.19 ± 0.09° (K2), +0.08 ± 0.09° (PSR), −0.06 ± 0.14° (K1)**. The Silk is 0.9–1.2° outside the sector's peak. Against rotations about the same centre, p = 0.03–0.05 |
| **FUV/NUV east of [O III]** | skeptic 1.4–2σ (jackknife); obs z 0.9–1.5 | **FUV δ = +0.067° (S1), p = 0.09** (z 2.4 with heavy tails). NUV +0.039°, p = 0.14. With the W strand (S12): FUV +0.085°, p = 0.12; NUV +0.081°, p = 0.11. **The FUV peak is on the [O III] ridge** (mode −0.014 ± 0.045°). The "offset" is an **east wing** (FUV σ_W/σ_E = 0.49; vs null p = 0.13) |
| **DESI [O II] vs [O III]** | skeptic +0.04°, wrong sign, "not robust" | **+0.043° (S1), p = 0.27; +0.051° (S12), p = 0.12; north +0.014°, p = 0.37** relative to co-spatial. The shock prediction (−0.045°, [O II] downstream) is disfavoured: **p = 0.067 (S1), < 0.003 (S12), 0.02 (N)**. DESI and NSNS [O III] ridges agree to −0.006 ± 0.027° |
| **X-ray ridge or edge at the spine** | — | No ridge: −0.30 ± 0.60 cts s⁻¹ deg⁻² (−4 ± 9 % of the level). The step at x = −0.2° is not Silk-specific: amplitude p = 0.37 vs rotations, and 21 % of rotations put a step within 0.3° of the same place |
| **South profile symmetry** | skeptic σ_W/σ_E 0.98 ± 0.05 vs limb 2.8–3.0; obs 1.16, p = 0.21 | **1.07 (block bootstrap 1.00–1.16).** Estimator variants 0.98–1.13 span both judges. The real-sky scatter for a symmetric truth is **±0.17**, not ±0.05. Symmetric null p = 0.18. A full smooth limb (q ≫ 1) gives r ≤ observed in 2/44 placements (**p ≈ 0.04–0.09**, about 2σ). A flat ribbon, or a curved sheet with LOS sagitta ≲ its depth (q ≲ 1, limb gain G ≲ 15), fits |
| **North one-sided** | obs: W/E 1.4–3.2, p 0.02–0.05 | **N 2.28 (2.03–2.62) vs symmetric null 1.01 ± 0.08: p = 0.023 (floor; Gaussian z ≈ 16).** N3 1.93 (p = 0.025), N2 2.07 (p = 0.025). N1 (Dec 20–21.5, patch masked) gives 1.13, p = 0.21, which is **not** significant here. Flat ribbon (q ≤ 1) excluded. A curved sheet with q ≈ 2–4 (G ≈ 17–24) fits. **A one-sided chord has the identical shape** |

---

## 1. Chance coincidence with the Monogem radius (`w4_geostats_registration.py`)

### 1.1 eRASS1 0.2–0.6 keV shell radius in the Silk sector (new measurement)

**Method**
- Radial profiles in 0.2° bins in the PA sector spanned by the southern spine (Dec 12–17.6) about each centre.
- Smoothed with σ = 0.4°.
- Three radius estimators:
  - r_peak: the maximum at 6–16°;
  - r½: the outer half-maximum between the peak and the median level 3–6° beyond it;
  - r_grad: the steepest decline.
- Errors: bootstrap over 5° PA sub-sectors.

| Centre | Sector PA (°) | r_peak | r½ (16–84 %) | r_grad | Peak / outside | Silk south mean | Silk − r½ | Silk − r_peak | Ring r½ in other 15° sectors, median ± PA sd |
|---|---|---|---|---|---|---|---|---|---|
| K2 (104.15, +15.30) | 78–104 | 10.9 | **11.91 (11.87–12.04)** | 11.7 | 1.30 | 12.10 | **+0.19** | +1.2 | 10.6 ± 2.2 (n = 22) |
| PSR B0656+14 | 72–100 | 10.3 | **11.27 (11.19–11.37)** | 11.1 | 1.33 | 11.36 | **+0.08** | +1.06 | 10.46 ± 2.32 (22) |
| K1 (107.08, +13.50) | 65–98 | 8.5 | **9.47 (9.32–9.59)** | 9.1 | 1.31 | 9.41 | **−0.06** | +0.9 | 9.97 ± 2.0 (21) |

**Readings**
- The southern spine runs along the **outer half-maximum contour of the X-ray shell in its own sector**, 0.9–1.2° beyond the brightness peak. This confirms `verify_xrayfuv`'s "outer flank", now located precisely.
- This is a property of the sky contour, so every centre gives it.
- The skeptic's "registration fails for K1" comes from comparing with Knies' **global** single-SNR radius (11.97°), which does not describe the shell edge in this sector (9.47° from K1).
- The Knies two-SNR global R = 11.0° is 0.9° inside the sector's r½. The Silk sector is one of the larger-radius parts of a clearly non-circular shell: the PA scatter of r½ is ±2.0–2.3°, partly estimator noise.

**Local-edge rotation null.**
- Rotating the southern section about each centre keeps its distance to that centre fixed and changes only the PA. The sector edge is then re-measured.
- P(|Silk − r½| ≤ observed) = **0.032 (K1), 0.048 (K2), 0.048 (PSR)**, from 61 rotations at 5° steps with |θ| ≥ 30°.
- P(|Δ| ≤ 0.3°) = 0.08–0.13.
- The rotated sectors overlap; roughly 12 are independent. Read this as **p ≈ 0.03–0.05, coarse, about 2σ**. The three centres are not independent tests of it.

### 1.2 The unified null

**H0.** The Silk is unrelated to Monogem. The full 66-point trace (rigid shape) is placed with:
- centroid uniform in area within |b − 20.46°| ≤ 3° (the observed centroid is (l, b) = (200.13, +20.46));
- |l − 201.5| ≤ 45°;
- uniform spin.

**Event E.** For any centre c ∈ {K1, K2, PSR}, some contiguous 23-point window (5.5°, the length of the southern section) has:
- rms(distance) ≤ r₀ = 0.15°, and
- |mean − R_c| ≤ z₀σ_c with z₀ = 1.

**Radii**
- K2: 11.00 ± 1.10 (Knies).
- K1: 11.97 ± 1.20. The 10 % uncertainty is **assumed** like K2's; Knies quote none in our notes.
- PSR: 10.46 ± 1.05. This is the eRASS1 ring-median r½ measured here, with ±10 %.

**Observed**
- K2: event. Best window Dec 12.25–17.75, mean 12.09°, rms 0.110°, z = +0.99.
- PSR: no event at r₀ = 0.15. Mean 11.36°, rms 0.154°. It is an event at r₀ = 0.25.
- K1: no event (9.41°, z = −2.1).

**Results** (600,000 placements per variant; Monte Carlo error ≲ 0.001)

| Variant | P(E), unconditional in band | P(rim zone ±1.5°) | **P(E \| rim zone ±1.5°)** | P(E \| ±3°) | Per centre (±1.5°) K1 / K2 / PSR |
|---|---|---|---|---|---|
| **primary** (K1/K2/PSR) | 0.030 | 0.33 | **0.092** | 0.080 | 0.054 / 0.042 / 0.037 |
| mirror images allowed | 0.030 | 0.33 | 0.091 | 0.079 | — |
| b band ±6° | 0.022 | 0.31 | 0.072 | 0.062 | — |
| PSR σ = PA scatter (2.3°) | 0.036 | 0.33 | 0.110 | 0.096 | PSR 0.059 |
| all radii from eRASS1 r½ ± PA sd | 0.035 | 0.27 | 0.122 | 0.108 | — |
| + G205.6+12.4 (8.82 ± 0.88) | 0.036 | 0.33 | 0.109 | 0.095 | G205 0.035 |
| K2 only (a-priori Monogem centre) | 0.010 | 0.24 | **0.043** | 0.034 | — |

**Sensitivity** (primary, conditional on ±1.5°):

| | z₀ = 0.5 | z₀ = 1 | z₀ = 1.5 |
|---|---|---|---|
| r₀ = 0.12 | 0.032 | 0.059 | 0.081 |
| r₀ = 0.15 | 0.051 | 0.092 | 0.125 |
| r₀ = 0.25 | 0.105 | 0.181 | 0.238 |

### 1.3 Bayes factor

**Definition.** BF = P(E | H1) / P(E | H0, rim zone).
- P(E | H1) = P(|d̄ − R| ≤ 1σ) × P(rms ≤ 0.15 | limb filament) = 0.683 × (0.5–1.0) = 0.34–0.68.
- The second factor is a judgement, because no limb-filament population was simulated.

**Results**
- **BF = 3.7–7.4 (primary)**, and 2.8–9.4 across the variants.
- It is 8–16 if K2 is the only a-priori centre.
- Using the ±3° zone gives 4.2–8.5.

**Why condition on the rim zone.** Every judge's prior was "given a shock next to Monogem's rim". The rim-zone conditioning matches that prior, so the fact that the Silk lies near Monogem is not counted twice.

**Combined with the local-edge test.** The two are correlated: a filament on the local edge will usually also match a nearby circle. A combined likelihood ratio (H1 vs unrelated) of **~4–10** is the honest summary, about a 2σ coincidence.

### 1.4 Why the judges disagreed (reproduced with their own definitions)

**judge_phys**
- Definition: the 23-point south segment only, centroid uniform within θ_max of K2, event |mean − 12.10| ≤ 0.5° and rms ≤ 0.25°.
- Reproduced here: **P = 0.0054 (θ_max 20°) / 0.0026 (30°) → BF 58 / 122**, matching their 60–120.
- The value is large because:
  1. the null is **not conditioned on reaching the rim**: most placements are far from Monogem, which adds a factor 1/P(rim) ≈ 3–10 that belongs in the prior;
  2. it uses a 1°-wide window around the **observed** value against a limb density of 0.31 per degree;
  3. it applies look-elsewhere afterwards as a crude ÷6.

**judge_skeptic**
- Definition: K2 (or K1) alone, whole trace, rms ≤ 0.15, |mean − R| ≤ 1.1, conditional on touching the rim at ±1.5°.
- Reproduced here: **P = 0.031 (K2), 0.030 (K1)**, matching their 0.031–0.069.
- Their "p ≈ 0.06–0.15 for 2–3 centres" is close to the 0.09 (0.07–0.12) found here by counting the centres inside the null.
- Their K1 "no event" is a global-radius artefact (§1.1).

---

## 2. Ionisation order across the filament (`w4_geostats_order.py`)

**Estimator** (all bands)
- Block-stacked profile: 0.25° blocks, 0.025° bins (maps); a single block and 0.1° bins (DESI).
- Template shift δ against the NSNS [O III] split-normal fit of the same segment (T_N centred on the NSNS ridge x₀_N), with a linear baseline over |x| ≤ 0.9° (DESI 1.2°).
- Variants:
  - free Gaussian centre;
  - fit windows of 0.6 and 1.2°;
  - the split-normal mode.

**Nulls**
- **GALEX: placement + injection.**
  - The trace is shifted to off-Silk sky: dRA ±0.8…3.0 in 0.2° steps × dDec −1.5/0/+1.5; pixels within 0.6° of the spine are excluded; ≥ 85 % of bins must be finite. This leaves 31 (FUV S1), 36 (NUV S1), 15/17 (S12) usable placements.
  - Three pixel-level injections at the observed amplitude:
    - **coshape**: the NSNS shape on the [O III] ridge (the judge_obs definition; tests "distributed like [O III]");
    - **ownshape**: the band's own fitted shape with its mode on the ridge (tests mode displacement);
    - **pred**: the NSNS shape displaced by the shock prediction (FUV +0.065°, the 50 % C IV mix of `w3_judge_skeptic_modelprof`; NUV +0.035°).
- **NSNS ridge itself.** Off-Silk placement with self-injection gives a scatter of x₀_N of 0.007–0.009° (robust), 0.018–0.041° (sd). This is negligible.
- **DESI.** Fibres lie at |x| ≲ 2°, so off-Silk placements are impossible. Instead:
  1. petal-group bootstrap (300);
  2. **petal-orientation injection null**: fit a co-spatial model, keep each fibre's residual, mirror x for a random half of the petal groups, and re-inject the co-spatial model (300). This includes petal-coherent baseline systematics;
  3. along-spine jackknife.
- **eRASS1.** 34 rotations of the trace about PSR B0656+14 (same radius, |θ| 20–180°).

| Segment / band | Observed δ (°) | Co-spatial null median ± sd | **p (one-sided, toward the shock prediction)** | Other variants | Prediction injection: median δ; P(δ ≥ obs) |
|---|---|---|---|---|---|
| S1 FUV (coshape δ) | **+0.067** (windows 0.042–0.067; Gaussian centre +0.076) | −0.004 ± 0.029 (n = 31) | **0.09** (Gaussian-centre variant 0.06); z = 2.4 | jackknife ± 0.030 | +0.057; 0.39 |
| S1 FUV mode | −0.014 | −0.001 ± 0.045 | 0.69 | pred mode +0.069 ± 0.033 → LR(mode displaced / co-spatial) = 0.23 | — |
| S1 FUV σ_W/σ_E (east wing) | 0.49 | 1.15 ± 0.42 | 0.13 | NSNS 1.07 | — |
| S1 NUV | +0.039 | +0.004 ± 0.030 (n = 36) | 0.14 (z 1.2) | mode −0.001, p = 0.68 | +0.032; 0.47 |
| S12 FUV | +0.085 | −0.006 ± 0.032 (n = 15) | 0.12 (floor 0.06) (z 2.8) | mode +0.015, p = 0.44 | +0.055; 0.33 |
| S12 NUV | +0.081 | +0.018 ± 0.035 (n = 17) | 0.11 (z 1.8) | mode +0.054, p = 0.17 | +0.047; 0.12 |
| **S1 DESI [O II] − [O III]** (DESI-[O III] template) | **+0.043** | orientation −0.027 ± 0.055; bootstrap sd 0.025; jackknife ± 0.020 | two-sided 0.27 vs co-spatial | — | pred −0.045 → median −0.067; **P(≥ obs) = 0.067** |
| **S12 DESI [O II] − [O III]** | **+0.051** | −0.001 ± 0.029; bootstrap 0.019 | two-sided 0.12 (bootstrap z 2.7) | — | median −0.022; **P(≥ obs) < 0.003 (0/300)** |
| N DESI [O II] − [O III] | +0.014 | −0.004 ± 0.019 | two-sided 0.37 | — | median −0.026; **P = 0.02** |
| DESI [O III] − NSNS | −0.006 (S1); −0.031 (S12); +0.006 (N) | ± 0.027 / 0.024 / 0.013 | 0.81 / 0.15 / 0.78 | — | registration check |
| eRASS1 ridge at the [O III] ridge (S1 / S12) | −0.30 / −0.42 cts s⁻¹ deg⁻² | −0.15 ± 0.60 / −0.19 ± 0.50 (rotations) | none (z −0.2 / −0.5) | −4 % / −6 % of the level (null sd 12 % / 10 %) | — |
| eRASS1 step (S1 / S12) | x_step −0.2°; step 0.95 / 1.10 | step 0.77 ± 0.50 / 0.58 ± 0.42 | 0.37 / 0.20 | 21 % / 12 % of rotations have x_step within 0.3° | — |

**Readings**
1. **FUV/NUV.** The eastward FUV excess is a **~1.3–1.7σ effect (p = 0.06–0.13)**, not 2–4σ.
   - The Gaussian z (2.4–2.8) overstates it because the off-Silk GALEX nulls have heavy tails.
   - Using the judge_obs definition, the significance rises slightly from their z ≈ 1.4 because 31 placements now replace 24.
   - The **FUV mode is co-spatial with [O III]**. What moves the template shift is an east wing.
   - That is the shape expected if part of the FUV (C IV) comes from a zone just behind an east-facing front while the rest (O III], 2γ) is co-spatial. The skeptic's 50 %-C IV mix gives +0.035…+0.065° with this estimator.
   - A wholesale displacement of the whole FUV ridge by +0.065° is disfavoured by ×4 (mode test).
   - Likelihood ratio for the predicted mixed shift vs co-spatial: Gaussian 14 (S1), tail-based about 4 (P(≥ obs) 0.39 vs 0.09). **Mildly favours the shock sign; not a detection.**
2. **[O II].** DESI puts [O II] **east** of [O III] by +0.04–0.05°. This is significant only at p = 0.12–0.27 against co-spatial.
   - The simple 1-D shock stratification ([O II] 0.045° **west**, downstream) is disfavoured at p = 0.067 (S1, clean), < 0.003 (S12) and 0.02 (N).
   - The W strand has the same [O II]/[O III] as E in SDSS (2.04 vs 2.01), so S12 is not obviously biased.
   - DESI [O III] is registered to NSNS within 0.006 ± 0.027°, so this is not a DESI astrometry problem.
   - It remains vulnerable to residual WIM-[O II] structure. The WIM [O II] is 2–3 R against about 1.8 R from the Silk, and it is removed only as a petal mean plus a linear baseline.
   - **This is the one ordering statistic that disagrees with the shock picture at ≥ 2σ.** It is consistent with FUV and [O II] both sitting slightly east, i.e. with no monotonic stratification. The skeptic's "two components or a systematic" stands.
3. **X-ray.** No X-ray ridge or Silk-specific edge exists at the spine. The X-ray structure there is the generic outer flank of the shell (§1.1).

---

## 3. Cross-profile shape (`w4_geostats_shape.py`)

**Estimator.** Primary: block medians (0.25° blocks, 0.025° bins), split normal plus linear baseline, unweighted, |x| ≤ 0.9° → r = σ_W/σ_E. R_asym = ⟨I−b⟩(−0.45…−0.15) / ⟨I−b⟩(+0.15…+0.45). Variants reproduce each judge's estimator.

**Geometric model** (small-angle). The emitting layer has depth w (top hat) and sits on a sheet with curvature R_c, extending l along the line of sight from the tangent point:
- h = u + z²/2R_c, so N(h) ∝ √min(h, qw) − √max(h − w, 0).
- q = l²/(2R_c w) is the LOS sagitta over the layer depth.
- **q → 0 is a flat ribbon; q ≫ 1 is a curved limb** with a √h rise and a 1/√h tail toward the centre of curvature.
- **A one-sided chord (z ∈ [0, l]) gives exactly the same N(h)** at half the column. Profile shape therefore cannot separate it from a full limb; only the gain (×2) and the velocity signature (`w3_northsouth`) can.
- The profile is convolved with a blur σ_b (NSNS PSF + trace error ≈ 0.03–0.06°; 0.10° tested). w is solved so the model matches the observed FWHM.
- Limb gain G = 2√(2qR_c/w) for q ≤ 1, with R_c = 58 pc (the Monogem limb used in closure W) and w in pc at 282 pc.

**Nulls.** Each model, plus a symmetric Gaussian, is injected with the observed amplitude into 39–44 off-Silk NSNS placements (dRA ±1.5…5, dDec −2/0/+2; pixels within 1° of the spine excluded).

| Segment | **σ_W/σ_E** (block bootstrap 16–84 %) | Variants (single-block median / 0.02° weighted / \|x\| ≤ 0.6 / block mean) | Symmetric-truth null mean ± sd; p(≥ obs) | Full limb (q = 10; σ_b 0.03 / 0.06) injected mean; P(≤ obs) | Allowed q (rel. likelihood > 0.1; σ_b 0.03–0.06) | G at allowed q |
|---|---|---|---|---|---|---|
| **S1 Dec 12–15** | **1.07 (1.00–1.16)**; R_asym 1.12 | 0.98 / 1.13 / 1.00 / 1.11 | 0.98 ± 0.17; **p = 0.18** | 3.24 / 3.46; **P = 0.044 / 0.044** (σ_b 0.10: 0.044; q = 2 at σ_b 0.10: 0.089) | **0 – ~1.5** (q ≥ 2: 0.02–0.06) | **≲ 15** (q = 1: 14.9–15.3) |
| **N Dec 20–28** | **2.28 (2.03–2.62)**; R_asym 2.61 | 2.32 / 2.86 / 1.65 / 1.90 | 1.01 ± 0.08; **p = 0.023 (floor)** | 4.05 / 4.21; P(≤ obs) = 0.023 (all 42 injections above the observed value) | **≈ 2–4** (q ≤ 1: < 10⁻³) | ≈ 20–24 |
| N1 Dec 20–21.5 (patch masked) | 1.13 (1.02–1.39) | 1.16 / 1.25 / 1.00 / 1.14 | 1.01 ± 0.10; p = 0.21 | — | 0–1 | — |
| N3 Dec 21.5–24.25 | 1.93 (1.75–2.18); R_asym 5.3 | 1.89 / 2.18 / 1.69 / 1.95 | 1.02 ± 0.10; p = 0.025 (floor) | — | ≈ 2 | ≈ 23–25 |
| N2 Dec 24.25–28 | 2.07 (1.77–2.56) | 2.34 / 2.25 / 1.96 / 1.79 | 0.98 ± 0.18; p = 0.025 (floor) | 2.55 / 2.76; P(≤ obs) 0.25 / 0.13 | ≥ 2 | ≈ 16–17 |

**Readings**
1. **South (S1) is symmetric**, and both judges' numbers are estimator variants of the same profile (0.98–1.16).
   - The real-sky scatter of the estimator for a symmetric truth is ±0.17. For a limb truth it is ±1.0–1.7, because split-normal fits of limb profiles on structured sky are unstable.
   - The skeptic's ±0.05 jackknife therefore overstated the precision about 3.5×.
   - **A full smooth limb (the skeptic's 2.8–3.0) is disfavoured at p ≈ 0.04–0.09 (about 2σ), not excluded.**
   - The data allow a **flat ribbon, or a sheet with Monogem curvature whose LOS half-length is ≲ √(2R_c w) ≈ 15 pc** (q ≲ 1–1.5, with w ≈ 0.40–0.42°, i.e. 2.0–2.1 pc at 282 pc, at q = 1).
   - Such a sheet still has **G ≲ 15**, against the closure-W G ≈ 20. Closure W's pressure match is therefore weakened by ≲ ×1.3–2, not invalidated. For a true flat ribbon G is a free parameter, which is the skeptic's point.
2. **North is one-sided**, with the sharp edge on the east, in N, N3 and N2 at the placement floor (p ≤ 0.025; Gaussian z 6–16).
   - A flat ribbon is excluded there.
   - A curved sheet with q ≈ 2–4 and G ≈ 17–24 fits.
   - Note that q ≈ 2 is preferred over q ≫ 1 when the blur is small: a full limb with σ_b ≤ 0.06 would be more one-sided (about 4) than observed.
   - **Limb and one-sided chord are degenerate in shape.** The judge_obs N1 result (1.40, p at floor) is **not reproduced** with a 0.9° patch mask (1.13, p = 0.21). The one-sidedness rests on Dec ≥ 21.5.
3. **Geometry.** In the south the Silk behaves like a nearly flat (R_curv 67°, `w3_northsouth`) ribbon or a gently curved sheet seen with limited LOS extent. In the north it is a curved limb or chord, concave to the west.
   - This fits judge_phys's H1b: a flatter sheet touched tangentially by the blast wave in the south, with the curvature appearing where the arc bends away.
   - It equally fits a sheet with its own curvature (H2).

---

## 4. Impact on the hypotheses

**Distance registration (H1)**
- Real but modest: a likelihood ratio of about **4–10** vs an unrelated filament (BF 3.7–7.4 any-centre, 8–16 K2 a priori; local-edge p ≈ 0.03–0.05).
- That is about **3–10× weaker than judge_phys's 60–120** and similar to the skeptic's 3–5.
- **The skeptic's "K2-only" objection is removed.** The Silk lies on the local X-ray half-max edge in every centre frame, so the coincidence is with the sky contour, not with one decomposition.

**FUV order.** Supports the shock sign only mildly (p ≈ 0.09, LR ≈ 4–14), as an east wing rather than a displaced ridge.

**[O II] order.** Opposite to simple 1-D stratification at p = 0.003–0.07. This is the main open tension for any single-front model, and it applies equally to H1 and H2.

**South profile**
- Excludes a smooth full limb only at about 2σ.
- G ≲ 15 is allowed, so closure W is under-determined by ≲ ×2.
- The skeptic's "no limb ⇒ removes the pressure coincidence" is too strong. "Weakens it" is right.

**North profile.** A one-sided curved limb or chord with G ≈ 20 is confirmed robustly. It cannot separate H1 from H2.

---

## 5. Limitations

- The placement nulls use rigid RA/Dec shifts (GALEX, NSNS) or rotations (eRASS1).
- GALEX coverage limits S12 to 15–17 placements, so its p floor is 0.06.
- DESI has no off-Silk placements. Its orientation null assumes petal systematics are not east–west aligned across petals.
- The eRASS1 sector r½ is a smoothed-profile estimator. Its ±2° PA scatter mixes real asphericity with estimator noise (foreground structure, coverage).
- No exposure or absorption correction was applied here: the `verify_xrayfuv` absorption tests found the south interior has more HI, which works against the observed sense.
- P(E | H1) uses a judgemental factor of 0.5–1 for the tangency of a limb filament. The K1 radius uncertainty (10 %) is assumed.
- The sheet model uses a top-hat layer and a Gaussian blur. Shock emissivity profiles, ripples and superposed sheets were not modelled. The blur (0.03–0.06°) is estimated, not measured.
- The shock-predicted offsets (+0.065° FUV, −0.045° [O II]) are taken from `w3_judge_skeptic_modelprof.json` (curved-limb projection). A flat-ribbon projection gives similar [O II] values (−0.04…−0.055°).

---

## 6. Files (tag `geostats`)

**Scripts** (run with `~/.venvs/astro313/bin/python -I`):
- `scripts/w4_geostats_common.py`: loaders, stack estimator, fits, placements.
- `scripts/w4_geostats_registration.py`: §1, about 1.5 min.
- `scripts/w4_geostats_order.py`: §2, about 4 min.
- `scripts/w4_geostats_shape.py`: §3, about 4 min.
- `scripts/w4_geostats_fig.py`

**Data:**
- `data/derived/w4_geostats_registration.json` and `w4_geostats_registration_erass.npz` (sector and 5° sub-sector radial profiles)
- `data/derived/w4_geostats_order.json` and `w4_geostats_order_profiles.npz`
- `data/derived/w4_geostats_shape.json` and `w4_geostats_shape_profiles.npz`

**Figure:** `figures/w4_geostats_summary.png` (dpi 200, `silk_style`, inspected):
- (a) spine distance from K1/K2/PSR, with the sector X-ray half-max edges and peaks;
- (b) eRASS1 radial profiles in the Silk sector;
- (c) P(E | rim) per null variant, with the judges' definitions;
- (d) S1 cross-profiles in NSNS, FUV, NUV and DESI [O III]/[O II];
- (e) offsets against co-spatial injection nulls (±1σ, ±2σ), with the injected shock predictions;
- (f) σ_W/σ_E against the sheet parameter q for S1, N and N2, with the observed bands and the S1 limb gain G.
