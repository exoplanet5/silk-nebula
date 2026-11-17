# w3_judge_skeptic: devil's-advocate ranking of origin hypotheses for the Silk Nebula

Wave 3, judge tag `judge_skeptic`, written 2026-11-17. The brief was to assume the team consensus is wrong and build the strongest case for the alternatives. The consensus is an incompletely cooled (column-truncated) radiative shock on the Monogem limb at about 280 pc.

**Inputs.** All `notes/*.md`, with wave 2 superseding wave 1 and wave 3 (`w3_oii`, `w3_northsouth`, `w3_canon`, `w3_figs`) superseding wave 2 where they disagree.

**New work.** Three small tests use data already on disk (§3). No downloads were made. Every script that reads data was run with `python -I`.

**Labels.**
- **measured**: computed here from data on disk.
- **derived**: arithmetic on quoted inputs.
- **teammate**: quoted from a named note or file.
- **literature**: from memory or secondhand; approximate.

---

## 1. Bottom line

1. **The mechanism is probably right; the location is less secure than the team says.**
   - A shock-heated, Balmer-dark, column-truncated cooling layer still explains the most constraints. My probability that the [O III] is shock-heated (any driver) is about 0.91.
   - My probability that it is driven by **Monogem** and lies at **220–360 pc** is about 0.6–0.7. The modellers' implicit value is about 0.9.
   - The distance rests on one morphological coincidence. I quantify it here as a roughly 2σ effect (chance probability 0.03–0.07), and it holds only for one of the two Knies+24 decompositions.
2. **New against the limb geometry used in the closure** (measured). On the single-strand southern Silk (Dec 12–15), the stacked NSNS [O III] cross-profile is **symmetric**:
   - observed σ_W/σ_E = 0.98 ± 0.05 (range 0.89–1.09 over 12 estimator variants);
   - predicted by the consensus shock models projected on a 58-pc limb, with the same estimator and matched blur: 2.8–3.0.
   - A locally flat sheet predicts 1.0, which fits.

   The Silk is therefore not seen as the limb of a smooth sphere. That removes the "brightness + width ⇒ Monogem pressure" coincidence (`verify_shockmodel` closure W) as independent support for H1. It does not exclude a shock.
3. **New, mildly for an outward-moving shock** (measured). GALEX FUV sits **east of (outside) the [O III] ridge**:
   - Dec 12–15: +0.04…+0.065° ± 0.03 (1.4–2.1σ);
   - Dec 12–17.6: +0.07…+0.09° ± 0.03 (2.5–3.9σ, but the W strand biases this one);
   - predicted by the shock model: +0.035…+0.065° (50 % C IV) to +0.09° (pure C IV);
   - a hot-gas interface with Monogem's hot gas on the inside predicts the opposite sign.

   The DESI [O II]-vs-[O III] offset has the "wrong" sign (+0.04 ± 0.03°; the shock predicts −0.04°). It rests on structured, petal-relative DESI baselines and is not robust.
4. **Wave 3 reopened one non-shock option.** `w3_oii` showed two things:
   - The [O II] doublet does not establish n_e ≈ 100.
   - The [O II]7325 flux had been low by a factor of 1.8 (7325/3727 = 0.060, so T(O II) ≈ 3.5 × 10⁴ K).

   This removes `verify_altmodel`'s "two-temperature" and "high-pressure" fingerprints. A **single-temperature heated warm layer (5–6 × 10⁴ K)** now fails only:
   - NUV, by ×1.5–2.5;
   - [N II]/Hα, by ×2 (a systematics-limited measurement);
   - the lack of a heat source. It needs about 1–2 × 10⁻⁵ erg cm⁻² s⁻¹ face-on, which is a shock-sized energy flux.

   It reproduces the FUV **without** invoking C IV resonance-scattering suppression. I raise "other non-shock" (H6) from about 1 % to about 8 %. The altmodel verdicts have not been recomputed with the corrected inputs.
5. **My posteriors:**

   | Hypothesis | Posterior |
   |---|---|
   | H1 Monogem-limb truncated shock | 0.60 |
   | H2 other shock driver, d ≲ 1 kpc | 0.27 |
   | H3 paper's 160-pc SNR | 0.005 |
   | H4 superbubble wall at 1–2 kpc | 0.035 |
   | H5 photoionised / DIG | 0.01 |
   | H6 other non-shock | 0.08 |

   The machine-readable version is `data/derived/w3_judge_skeptic_posterior.json`.

---

## 2. What changed since wave 2 that matters for the ranking

| Change (source) | Effect on hypotheses |
|---|---|
| [O II] 3729/3726 = 1.31 ± 0.07 (DESI) is within 2σ of the low-density limit. The SDSS "confirmation" was a blue-end calibration artefact (`w3_oii`) | The "P/k ~ 3 × 10⁶, 100–400× Monogem" tension with H1 disappears. **But** the doublet argument that excluded every pressure-balanced alternative (interfaces, warm layers, fossils, precursors) also disappears |
| [O II]7320+30 was a normalisation bug (×1.8). 7325/3727 = 0.060 (0.055–0.073); T(O II) ≈ 2.4–5.5 × 10⁴ K (`w3_oii`) | Shock tension 10 is resolved. The altmodel two-temperature fingerprint is gone. The warm-layer, conduction-front and fossil strikes based on 7325 also go |
| North–south: one smooth R ≈ 21° arc; the south is straight (R_curv 67°, 27–106°) and only **tangent** to the 12° Monogem circle. The north lies 1–4° outside the X-ray rim. The north is slower (70–100 km/s) at lower pressure. The DESI N–S ΔV is +4…+9 km/s. The NVSS RM step is not Silk-specific (p = 0.18) (`w3_northsouth`) | Against a spherical-limb H1. Neutral to mildly for a separate coherent structure (H2). H1 survives only with a non-spherical blast wave and extra parameters for the north |
| [O II]/[O III] DESI–SDSS conflict resolved as a Dec-weighting effect (`w3_northsouth`) | Removes a systematics worry; neutral |
| Hα patch and PGCC clump lie exactly on the Monogem-radius circle (θ_K2 11.94–11.97°) (`w3_northsouth`) | A second, partly independent coincidence in favour of H1. It is 1.3° **inside** the [O III] spine at the same Dec |

---

## 3. New tests (this judge)

Scripts:
- `scripts/w3_judge_skeptic_tests.py`
- `scripts/w3_judge_skeptic_modelprof.py`
- `scripts/w3_judge_skeptic_robust.py`
- `scripts/w3_judge_skeptic_fig.py`

Outputs:
- `data/derived/w3_judge_skeptic_{tests,modelprof,robust,posterior}.json`
- `figures/w3_judge_skeptic_summary.png` (main)
- `figures/w3_judge_skeptic_tests.png` (diagnostic; one annotation runs slightly past the right axis edge)

### 3A. How special is the Monogem-limb registration? (measured)

**Observed.**
- The 23 southern spine points lie at 12.10 ± 0.12° (rms) from K2 (R = 11.00 ± 1.1°).
- The spine is tangent to the K2 circle within a median 7°; it is 54–69° off tangent at the very southern tip.
- Measured from **K1** (Knies single-SNR, R = 11.97°), the same points are at 9.41 ± 0.29°, i.e. 2.6° inside the rim. The registration exists **only** in the two-SNR decomposition.

**Null test.**
- I kept the real trace and drew random circle centres around its midpoint. This is equivalent to a random rigid translation and spin of the trace near a fixed Monogem.
- Event: some contiguous 23-point window has rms ≤ 0.15° and |mean − R| ≤ 1.1°.
- 40,000 draws per centre.

| Centre | P(event), all draws | P(event \| trace reaches rim ±1.5°) | P(event \| trace crosses rim) |
|---|---|---|---|
| K2 | 0.013 (centre ≤ 26°) – 0.057 (10–16°) | **0.031–0.069** | 0.024–0.050 |
| K1 | 0.013–0.057 | 0.028–0.064 | 0.021–0.042 (observed: no event) |

**Reading.** Suppose a long filament reaches Monogem's rim. A tangent, rim-radius section appears by chance 3–7 % of the time. Allowing for the choice among two or three centres, p ≈ 0.06–0.15. This is real but modest evidence, worth a likelihood ratio of about 3–5. It is the **only** positive distance evidence for 252–318 pc.

### 3B. Cross-filament profile shape and ordering (measured vs model)

**Method.**
- Data, using the `w3_figs` F6 machinery:
  - NSNS [O III], 1′ map;
  - GALEX FUV and NUV 1′ skybg blocks, with per-tile demeaning;
  - DESI fibre line fluxes.
- Fits: split-normal plus linear baseline.
- Shifts: a template shift against the NSNS [O III] fit.
- Errors: along-spine jackknife over 10 chunks.
- Model: the same estimators applied to the 25 "consistent" `verify_shockmodel` tracks.
  - Each track is projected exactly as `verify_shockmodel_analysis.profile()` does (curved limb, R_c = 58 pc, 280 pc).
  - It is also projected as a flat edge-on sheet.
  - Both are convolved to the observed stacked FWHM (0.41°).

| Statistic (Dec 12–15, one strand) | Observed | Shock, curved limb | Shock, flat sheet | Hot-side-west interface |
|---|---|---|---|---|
| [O III] σ_W/σ_E | **0.98 ± 0.05** (variants 0.89–1.09) | 2.8–3.0 | 1.0 | symmetric or west-extended |
| Non-parametric R_asym = ⟨I⟩(−0.45…−0.15) / ⟨I⟩(+0.15…+0.45) | 1.04–1.13 ± 0.08 | 3.4–4.6 | 0.85–1.66 | — |
| Δ(FUV − [O III]) | **+0.040…+0.065° ± 0.03** | C IV +0.065…+0.135; 50 % C IV mix +0.035…+0.065 | C IV +0.05…+0.145 | < 0 (illustrative) |
| Δ(NUV − [O III]) | +0.025…+0.035 ± 0.016 | not modelled (C III]/C II]/2γ lie between C IV and [O III]) | — | < 0 |
| Δ(DESI [O II] − [O III]) | +0.035…+0.055 ± 0.03 (moments +0.05…+0.11) | −0.035…−0.055 | −0.04…−0.055 | > 0 |

For Dec 12–17.6, which includes the W strand:
- NSNS σ_W/σ_E is 1.3–2.0 (the W strand adds west emission);
- FUV shift +0.07…+0.09 ± 0.03;
- NUV shift +0.065…+0.08 ± 0.02;
- DESI [O II] shift +0.05 ± 0.01–0.03.

**Readings.**
1. **No limb signature.** The emitting structure is a locally flat sheet, ribbon or tube, or the sum of many small ripples. It is not the edge of a smooth ~60 pc sphere.
   - This invalidates the G ≈ 20 used in closure W, which gave "1.3–4 R at the Monogem pressure for b = 7–14". It also affects G_CIV, which sets the C IV-scattering suppression.
   - For a flat ribbon G = L/Δ is a free parameter, so pressure, n₀ and FUV suppression are each uncertain by ×2–5.
   - The coincidence that "supports H1" becomes a fit with one more free parameter.
2. **FUV outside [O III]** matches the sign and size of an outward shock (front → C IV → [O III]) driven from the concave (west) side. Both Monogem (H1) and any west-side driver (H2) qualify. It argues against interfaces heated from Monogem's side and against a static heated layer (H6), at about 2σ in the clean segment.
3. **DESI [O II] east of [O III]** is opposite to the shock prediction.
   - Taken together with (2), it would mean no monotonic stratification. That is impossible for any single 1-D structure, which suggests either two components or a systematic.
   - The DESI east-side baseline is structured: about −0.7 R at +0.35…+0.65° in the [O III] profile.
   - [O II] carries about 2–3 R of WIM foreground. The DESI Hα profile shows no matching structure, but it is noisy.
   - **Not robust; flagged for a dedicated test.**

**Caveats.**
- The NSNS DR0.2 processing (star removal, background subtraction) is not documented here. A large-scale background filter could suppress a long tail, but not the near-peak asymmetry within ±0.3°, which is what σ_W/σ_E measures.
- The GALEX skybg maps are smoothed by the pipeline. A symmetric smoothing kernel cannot create a shift.
- The jackknife ignores structure that is coherent along the spine, and no placement or injection null was run for the shift statistics.
- The "interface" bands in the figure are qualitative signs, not models.

### 3C. Order-of-magnitude checks used in §4 (derived)

| Check | Calculation | Result |
|---|---|---|
| Energy flux any mechanism must supply | 2 R at G ≈ 19 gives about 1.1 × 10⁵ photons cm⁻² s⁻¹ face-on at 5007 (altmodel). With f(5007) = 2.5 % (warm layer) | 1.7 × 10⁻⁵ erg cm⁻² s⁻¹ |
| Same, for a shock | ½ρv³ at n₀ = 0.025, v = 120 km/s | 5 × 10⁻⁵ erg cm⁻² s⁻¹ |
| Implication | | Any alternative needs a shock-sized energy flux |
| Brightness closure vs distance | `verify_shockmodel` gives 1.5–4.5 R at 250 pc and 1.1–3.5 R at 320 pc, so I_pk ∝ 1/d at fixed angular width (curved-limb G) | 160 pc: 2.3–7 R; 1 kpc: 0.4–1.1 R; 1.6 kpc: 0.23–0.7 R. Observed: 1.4–2.5 R south |
| Implied distance window | | ≈ 150–800 pc. Weakens further now that the curved-limb G is in doubt |
| Driving pressure at 160 pc | n₀ ∝ 1/d at fixed width and N* | P/k ≈ 5.5 × 10⁴ – 1.9 × 10⁵, i.e. 5–20× the Local Bubble's P/k ≈ 1 × 10⁴ (literature, approximate). The LB cannot drive the Silk without a recent energy injection |
| Charge exchange | O-ion flux available from Monogem gas, n_O ≈ 2 × 10⁻⁶ cm⁻³, against the ≈ 1 × 10⁵ cm⁻² s⁻¹ needed (at most 1 photon per ion) | Would need v > c. **Excluded** |
| Cosmic-ray streaming heating | v_A ΔP_cr with v_A ≈ 32 km/s (B 3 µG, n 0.03) gives 5 × 10⁻⁶ (ΔP_cr / 1 eV cm⁻³) erg cm⁻² s⁻¹ | Needs a CR pressure drop of ≈ 3 eV cm⁻³ across ≲ 1 pc. Energetically conceivable at an SNR's CR-escape boundary; sharpness unexplained |
| κ-distributed electrons | e–e equilibration takes ~10⁴ s at n_e = 1, T = 10⁴ K | Needs continuous suprathermal injection at the same power as above. Not an energy source |

---

## 4. Hypothesis scorecard

| # | Hypothesis | Posterior | Verdict | Strongest evidence for | Strongest evidence against |
|---|---|---|---|---|---|
| **H1** | Truncated radiative shock on the outer limb of the Monogem SNR, 250–320 pc (here 220–360, including non-spherical variants) | **0.60** | **favoured, not established** | Fesen+26 find [O III] filaments around most of Monogem's X-ray periphery, with the southern ones also just outside the X-ray edge. The south spine sits at the K2 radius (chance 3–7 %, §3A). The Hα patch and PGCC also sit on the K2 radius. The X-ray-bright side faces Monogem. FUV lies outside [O III] (outward shock, §3B). Truncation is natural (t* < t_cool, comparable to the 1.1 × 10⁵ yr age). Dust at 300–600 pc rises north and 4° east (E 0.010–0.014 vs 0.003–0.005 south; density peaks 349–367 pc), consistent with denser gas beyond the limb (context CSV, my reading) | One-sided morphology only. The registration fails for K1. The south is tangent but not concentric (R_curv 67°). The north 7° of the same smooth arc lies 1–4° outside the X-ray rim with no X-ray step. No limb profile (§3B). The X-ray edge is generic (rotation p = 0.37). The pressure "coincidence" used a tuned b and the now-doubtful G. No kinematic or proper-motion evidence |
| **H2** | Truncated shock from a different driver, distance unconstrained. Split: unknown older SNR or bubble at 0.1–0.8 kpc 0.13; G205.6+12.4 (Knies' second, hotter remnant) 0.06; **Local Bubble wall 160–170 pc** 0.05; Geminga-progenitor fragment 0.02; other 0.01 | **0.27** | **viable** | One coherent R ≈ 21° arc across both halves with no break at the rim. That is natural for a structure with its own geometry, unnatural for a Monogem limb. 160 pc is not excluded by any measurement. The [O III] velocity equals local gas at any d < 1 kpc. No absorption bracket. The brightness closure allows 150–800 pc | No identified driver: no second X-ray interior, no pulsar birth locus, no OB group. G205.6 is not concentric (rim crossed at Dec 19). Geminga concentricity is unphysical. The LB cannot supply the ≥ 5 × 10⁴ K cm⁻³ ram pressure. No dust/HI feature at 160–170 pc along the Silk. The FUV order requires a west-side driver, and Monogem is the only known one |
| **H3** | The paper's model: isolated SNR at 160 pc, R ≈ 30 pc, E ≈ 2 × 10⁵¹ erg, 6.8 × 10⁴ yr, PDS phase | **0.005** | **ruled out** | None independent. 160 pc happens to sit on the LB wall | 160 pc is a unit error (Table 4 is in kpc). The inputs (10–13 R, T_e, V_t, n₀ = 1) are WIM-contaminated or 5–8× wrong. A PDS-phase radiative shell is a **complete** shock, which needs [O III]/Hβ ≤ 5 and an Hα filament; neither is seen. R = 30 pc is asserted (at 160 pc the 21° arc implies R ≈ 60 pc). No interior X-ray of its own |
| **H4** | Superbubble or chimney wall at 1–2 kpc (paper's kd output at face value) | **0.035** | **disfavoured** | Not excluded by absorption or extinction (E(B−V) to 1.25 kpc is only 0.02–0.05). The curvature centre (l, b) = (194.3, +0.4) lies in the plane, 0.8° from the Origem Loop centre, toward Gem OB1. Kpc-high walls and chimneys exist | The kd distance is an artefact (revised V_LSR ≈ +1 gives 0–0.2 kpc). 470 pc long and 7.8 pc wide at z ≈ 0.55 kpc, Balmer-dark and 1-pc-sharp only at small d. The closure predicts 0.23–0.7 R at 1.6 kpc (×2–11 too faint). No HI supershell or H II region at the [O III] velocity. The Origem Loop (6.7° across, ~200 pc diameter at ~1.7 kpc) is far smaller than the ~620 pc radius the 21° arc would imply there |
| **H5** | Diffuse ionised gas or photoionised bubble wall (the paper's alternative) | **0.010** | **ruled out** (U budget); spectrally about 2–3σ | [O III]/Hβ is only a lower limit (> 11 at 2σ). PIE with the hardest spectra reaches 10–14. The 4363 detection is only 2–3σ | No ionising source: U is short by ≥ 27× even for a coplanar WD J0736+1444, and by ~10⁴× at 282 pc. He II/[O III] < 0.026 (hot-gas PIE 0.06). [O III]/Hα ≳ 3 against WIM ~0.06. FUV per R 1.3–1.8× short. A transient-flash fossil needs a ≥ 4 × 10⁴ K start, a pre-existing filamentary density structure and an unrecorded event ≲ 4 kyr ago within ~0.3 kpc |
| **H6** | Other non-shock. Split: heated warm layer (5–6 × 10⁴ K; CR-streaming, turbulent or reconnection heating) 0.05; conduction front 0.01; mixing layer 0.005; NEI cooling of hot gas 0.003; κ electrons 0.005; photoionised fossil 0.003; charge exchange 0.001; unlisted 0.003 | **0.08** | **disfavoured** (warm-layer variant: viable at a low level) | Warm CIE layer: [O III]/Hβ ≈ 100 (SDSS 23.5, 16–84 % 12.7–47.6, so 1–2σ), [O II]/[O III] 1.2–1.6, 4363/5007 0.097, **FUV 86–92 CU/R with no C IV suppression needed**, 7325/3727 0.08 against 0.060 ± 0.015. Static gas at the ambient velocity matches the centroid. The symmetric profile fits a static layer | Needs an unidentified heat flux of ~2 × 10⁻⁵ erg cm⁻² s⁻¹ (a shock-sized budget) at the thermally unstable part of the cooling curve. NUV ×1.5–2.5 and [N II]/Hα ×2 too high. FUV east of [O III] (§3B) has no natural explanation. Conduction fronts are 30–1000× too faint per front with FUV ≥ 540 CU/R. A mixing layer needs a 377 km/s inflow. NEI cooling predicts an X-ray filament ≥ 10³× Monogem, not seen. CX is excluded by > 10⁴ (§3C) |

**Marginals** (derived from the table):

| Quantity | Probability |
|---|---|
| Shock-heated, any driver (H1 + H2 + H3 + H4) | 0.91 |
| Monogem-related, any mechanism | ≈ 0.69 |
| d = 220–360 pc | ≈ 0.70 |
| d ≈ 150–175 pc (LB wall) | ≈ 0.055 |
| d > 1 kpc | ≈ 0.04 |

**Calibration of my H1 number.** As a rough Bayes ledger, I start from a prior of 0.6 for Monogem-driven given a shock next to Monogem's rim. The likelihood ratios are:
- registration: ≈ 3–4;
- Hα patch on the radius: ≈ 1.5;
- only known driver: ≈ 1.5;
- north outside the rim, curvature mismatch, one smooth arc: ≈ 0.4;
- no limb profile: ≈ 0.8.

That gives P(Monogem | shock) ≈ 0.75–0.8. I then multiply by P(shock) ≈ 0.88 (excluding H3 and H4) and shade down for the systematics track record (§5 W3). The result is 0.60.

---

## 5. Weakest links in the team's analysis

1. **W1: Distance = one geometric coincidence, now quantified.**
   - Chance probability 0.03–0.07 (§3A), and it exists only for the K2 decomposition.
   - The south is straight and tangent, not concentric. The north 7° of the same smooth arc is off-sphere and outside the X-ray rim.
   - The X-ray "edge" is generic (p = 0.37), and the ~2σ local step peaks 0.5° inside the spine.
   - No absorption, kinematic or proper-motion constraint brackets the gas.
   - Yet 252–318 pc is propagated into every physical number: size, n₀, pressure, energetics.
2. **W2: The geometry behind the "Monogem pressure coincidence" is contradicted by the profile shape (new, §3B).**
   - The closure assumes a smooth 58-pc limb, G ≈ 20, and tunes b = 7–14 so that ρv² matches the Monogem pressure.
   - The Monogem pressure itself spans 0.8 × 10⁴ to 1.6 × 10⁵ across the team's definitions (canon C10).
   - The observed symmetric profile excludes the smooth-limb projection. The same G sets the C IV scattering suppression (G_CIV 19–35) used to rescue the FUV.
   - The shock model is not wrong, but two of its "successes" (pressure match, FUV after scattering) are now under-determined.
3. **W3: Key diagnostics are marginal, and they have been revised repeatedly by large factors.**
   - Hβ is a 1–2σ signal. [O III]/Hβ has only a lower limit (> 11). The withdrawn value of 46–47 was used to fit the shock model (fit target 40).
   - 4363 is a 2–3σ hint.
   - The doublet density has gone from 100 cm⁻³ to "uninformative".
   - 7325 changed by ×1.8 because of a code bug.
   - The DESI brightness went 2.4 → 1.62 R; the N–S ΔV went 10 → 4–9 km/s.
   - [N II]/Hα and [S II]/Hα are systematics-limited by ×2.
   - **Neither `verify_shockmodel` nor `verify_altmodel` has been re-run with the current inputs.** The altmodel exclusions of the warm layer, fossil and interfaces partly relied on the now-retracted n_e ≈ 100 and T(O II) = 1.25 × 10⁴ K.
4. **W4: No direct kinematic evidence of a shock.**
   - The intrinsic width is unresolved: the DESI LSF is ~98 km/s FWHM at 5007 against a predicted 23–38 km/s; LAMOST's 43–49 km/s is not a measurement.
   - The centroid equals the ambient gas, which any static layer would also give.
   - No proper motion has been measured.
   - No spectrum of the filament exists at R > 7500.
5. **W5: The UV argument is model- and calibration-dependent and only half-covered.**
   - FUV per R uses the NSNS scale, and the absolute [O III] scale is uncertain by 1.35–1.5× (canon C3).
   - C IV suppression is invoked post hoc (factors 3.7–26 across agents).
   - GALEX stops at Dec 22, so the northern third has no UV test.
   - The NUV per R (36–80) is not modelled by the shock code.

---

## 6. Decisive future observations

Predictions are given per hypothesis. "Shock" means H1 or H2, which differ only in distance and driver.

| # | Observation (instrument, line, depth or resolution) | H1 (Monogem shock, ~280 pc) | H2-LB (shock at 160 pc) / H2 other | H4 (1–2 kpc) | H6 warm heated layer | H5 photoionised |
|---|---|---|---|---|---|---|
| D1 | **[O III] 5007 line profile** on the sharp south ridge (Dec 12–15). Fabry–Pérot, spatial-heterodyne or echelle at R ≥ 50,000 (σ ≤ 3 km/s), S/N ≥ 20 per 2 km/s on a 2 R line, ~0.1–0.25° aperture | Flat-topped or double-horned, FWHM 23–38 km/s, centred on the pre-shock gas (V_LSR ≈ +1) | same as H1 | same, possibly broader | Gaussian, FWHM ≈ 13 km/s thermal (6 × 10⁴ K) plus turbulence ≲ 20 | FWHM ≈ 6 km/s thermal plus turbulence |
| D2 | **Proper motion** of arcsecond-scale [O III] substructure. Deep narrow-band imaging at ≤ 2″ now (the 2022–25 NSNS and NHZ frames are epoch 0) and again in ≥ 10 yr; needs ≤ 0.3″ registration | 75–112 mas/yr outward (E/NE): 0.75–1.1″ per decade | 132–198 mas/yr (1.3–2″ per decade) | 13–20 mas/yr (0.13–0.2″ per decade) | ≈ 0 | ≈ 0 |
| D2+D1 | v_s from D1 plus μ from D2 gives an **expansion-parallax distance** d = v_s / 4.74μ: the single most decisive distance test | ≈ 280 pc | ≈ 160 pc / other | ≥ 1 kpc | — | — |
| D3 | **FUV spectroscopy** of the south ridge: archival FIMS/SPEAR (STSAT-1) C IV 1548/1551, Si IV 1400, O III] 1664, He II 1640, binned over ~0.5 × 5°. First check Kim+2007, whose C IV was "bright only in a small high-latitude part of the ring" (secondhand). The Silk is Monogem's highest-latitude limb | C IV after scattering ~1.6–7 × 10⁴ LU at 2 R [O III]. 1548/1551 → 1.0–1.3 (thick). C IV halo broader than [O III] and centred ~0.05–0.1° outside | same as H1 | same | C IV ~6 × 10³ LU (C IV/[O III] 0.03–0.04 photons); Si IV ≈ 2× O III] | FUV dominated by O III] and 2γ, no C IV |
| D4 | **Hβ to ≥ 4σ.** Long-slit or IFU mosaic along the ridge on a 2–4 m telescope (2σ ≈ 0.03 R), or more DESI/4MOST sky fibres | [O III]/Hβ 10–63 (median ~20–49) | same | same | ≈ 100 | ≤ 10–14 |
| D5 | **[O II] 3726/3729 at R ≥ 10,000** with ±0.03 precision, plus Hα/[O II] co-spatial | ≈ 1.46 (low density) | ≈ 1.46 | ≈ 1.46 | ≈ 1.46 | ≈ 1.46 |
| D5 (falsifier) | | A value ≤ 1.35 at ≥ 4σ would mean n_e ≳ 60 and P/k ≳ 10⁶. That falsifies **every** pressure-balanced reading and demands dense clumps | | | | |
| D6 | **Absorption bracket.** HST/COS C IV, Si IV, N V at ~15 km/s resolution toward 3–6 Gaia hot stars (WD, sdB, BHB) within 0.3° of the south ridge at 150–450 pc, plus a foreground control. Predicted face-on N(C IV) ~ 1–3 × 10¹³ at V_LSR ≈ +1 (shock model) | Absorption only in stars beyond ~250–260 pc | Absorption already at ≥ 175 pc / at the driver's distance | Only in stars beyond ~1 kpc | Weaker C IV, strong Si IV | N(C IV) small |
| D7 | **Ordering test with nulls.** GALEX 5″ FUV ridge (re-extract at native resolution, not skybg) plus a matched-resolution [O III] map and DESI/SDSS [O II] positions, with injection and placement nulls | FUV +0.05–0.1° outside; [O II] 0.04° inside | FUV outside if the driver is west | — | No offset | No offset |
| D8 | **Faraday rotation / depolarisation** along the ridge: Effelsberg λ11/21 cm, Urumqi λ6 cm, GMIMS, CHIME or DRAO polarisation; compare with the R20 Monogem filament | Compressed B (x ≈ 5–12) gives RM ~ 50–100 rad m⁻² through an edge-on chord (derived, rough: 0.81 × 0.24 cm⁻³ × 12 µG × 25 pc); a depolarisation filament at λ11–21 cm | same | smaller, diluted | Uncompressed B, RM ≲ 15 | ≲ 15 |
| D9 | **eRASS:4–5 spectral mapping**: hot-gas pressure immediately west of the south ridge versus beyond Dec 21 | P/k 3–11 × 10⁴ behind the south ridge; a fainter, lower-P component reaching the north arc | A separate hot interior not following Monogem's radial profile | none | — | — |
| D10 | **3D dust (Edenhofer+24) slices at 200–450 pc** and the distance of PGCC G199.09+19.95. Test the dust rise at 349–367 pc north and east of the ridge | Pre-shock or denser gas just outside the ridge at ~280–340 pc; the cloud the patch shock hits at the same distance | Nothing specific at 280 pc; the LB wall is featureless along the Silk | — | — | — |

**Falsifying the consensus.** The most direct routes are:
- D1 shows a Gaussian FWHM ≤ 15 km/s.
- D2 shows μ ≈ 0, or μ ≥ 130 mas/yr with v_s ≈ 100–150 km/s (that is, d ≈ 160 pc).
- D6 shows Silk-velocity high ions toward a star at < 230 pc.
- D3 shows C IV/O III] ≲ 0.2 together with Si IV ≥ O III].

Any one of these overturns H1. A ≥ 4σ Hβ with [O III]/Hβ ≥ 80 (D4) would move most of H1's weight to H6.

---

## 7. Important things not yet done

- **Re-run `verify_shockmodel`** with an [O III]/Hβ lower-limit likelihood (not 40), the corrected 7325, the low-density doublet and a **flat-ribbon or rippled geometry**. Then redo closure W, G and G_CIV, and check whether slower (70–90 km/s) or more complete shocks re-enter.
- **Re-run `verify_altmodel`** with the `w3_oii` inputs, especially the warm heated layer, the shock-heated fossil (now allowed at n ≈ 0.5–1, e-fold ~kyr) and conduction fronts. Add a CR-streaming-heated layer model with an explicit thermal-stability analysis.
- **Placement and injection nulls** for the §3B statistics, and a native-resolution GALEX re-extraction (not skybg). A DESI [O II] ordering test that models the WIM [O II] foreground explicitly.
- **Literature check of Kim+2007 (FIMS/SPEAR)**: does the "small high-latitude part" with bright C IV coincide with the Silk? If so, it is a direct C IV measurement for D3.
- **Gaia DR3 search for COS-able hot stars** within 0.3° of the south ridge at 150–450 pc (for D6).
- **3D-dust distance of the PGCC/Hα-patch cloud** and of the 300–600 pc dust rising toward the north and east (`context_dust_los.csv`: E(300–600) 0.0028 at Dec 12 → 0.0135 at Dec 28; 4° E of P1 0.0124, peak 367 pc; 4° W 0.0033). Done only at the sightline level here.
- **No R > 7500 spectrum of the filament exists.** Every kinematic claim (width, shock velocity, N–S shift) rests on R ≲ 5000 data or on LAMOST at S/N ~6.
- **The north (Dec 22–28.6) has no UV data and no X-ray step**, yet it is 40 % of the arc. Its excitation (a slower shock) is inferred only from optical ratios with ×2 systematics.
- **No robustness test of the Knies two-SNR versus one-SNR choice**, which decides whether the registration exists at all. The eRASS1 shell radius should be measured along the Silk sector itself.

---

## 8. Limitations of this judgement

- The posteriors are subjective and use coarse likelihood ratios (§4). They are rounded to 0.005.
- The §3A null randomises the circle position around the real trace. It does not model the prior that [O III] filaments occur preferentially near SNR rims (Fesen's population), which is folded in separately as a prior.
- The §3B model comparison uses the shock tracks as stored: no radiative transfer, and truncated runs end at 10⁴ K. The convolution kernel is a single Gaussian matched to the stacked FWHM. The "flat sheet" is the emissivity profile across the layer with a uniform line-of-sight length.
- The "interface" predictions in §3B and Fig. (c) are qualitative signs, not computed models.
- The NSNS map and GALEX skybg products were used as processed by earlier waves (`w3_figs` F6 inputs).
- Order-of-magnitude numbers in §3C use approximate literature values: the LB pressure ~10⁴ K cm⁻³ and the CR pressure ~1 eV cm⁻³.

## 9. Files (tag `judge_skeptic`)

**Scripts** (run with `~/.venvs/astro313/bin/python -I`, in this order):
- `scripts/w3_judge_skeptic_tests.py`: registration null (§3A) and profile ordering (§3B); ~3 min.
- `scripts/w3_judge_skeptic_modelprof.py`: shock-model projected profiles with the same estimators.
- `scripts/w3_judge_skeptic_robust.py`: estimator variants plus the flat-sheet model.
- `scripts/w3_judge_skeptic_fig.py`: summary figure.

**Data:**
- `data/derived/w3_judge_skeptic_tests.json`
- `data/derived/w3_judge_skeptic_modelprof.json`
- `data/derived/w3_judge_skeptic_robust.json`
- `data/derived/w3_judge_skeptic_posterior.json`

**Figures** (dpi 200, `silk_style`, inspected):
- `figures/w3_judge_skeptic_summary.png`:
  - (a) spine distance from the K2 and K1 centres, with the chance probability;
  - (b) Dec 12–15 NSNS [O III] and GALEX FUV profiles against the shock model as a curved limb (dashed; plotted without baseline fitting) and as a flat sheet (dotted);
  - (c) observed σ_W/σ_E, ΔFUV and Δ[O II] against the model bands;
  - (d) posteriors.
- `figures/w3_judge_skeptic_tests.png`: diagnostic profiles for both segments, including NUV and DESI.

**Inputs read** (read-only):
- `silk_trace.json`
- `w3_figs_f1_maps.npz`
- the GALEX blocks via `verify_xrayfuv_galex_measure.prep_band`
- `verify_desi_linefluxes.npz`, `verify_desi_coadd_meta.npz`, `verify_desi_detect.json`
- `verify_shockmodel_summary.json`, `verify_shockmodel_tracks.npz`
- `context_dust_los.csv`
- all notes listed in the header
