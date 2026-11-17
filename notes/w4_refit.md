# w4_refit: shock and alternative models refitted with the corrected inputs

Wave 4 (close-out), tag `refit`. Written 2026-11-17.

**Precedence.** w3 > verify > wave 1. Inputs come from `w3_canon`, `w3_oii`, `w3_northsouth` and the three w3 judge notes.

**Code reuse.** All model physics is the wave-2 code, imported read-only:
- `verify_shockmodel_flow` and `_fit`: the 1-D NEI truncated-shock flow, with no radiative transfer.
- `verify_shockmodel_analysis` and `_rt`: projection, line profile and the C IV Monte Carlo.
- `verify_altmodel_lines`: CHIANTI 7.1.4 UV and optical emissivities.
- The altmodel JSON products, for the non-shock scenarios.

**What is new here**
- The likelihood (corrected inputs).
- An extended grid that adds UV lines.
- The flat and rippled geometries.
- The north pairing.

**Labels**
- **model**: my run of the wave-2 code.
- **derived**: arithmetic on model output.
- **input**: a canonical measured value.

---

## 0. Bottom line

1. **Truncated shocks still fit, but the optical spectrum now constrains them much less (model).** I used the corrected inputs:
   - [O III]/Hβ: a one-sided limit from the DESI Hβ box, 0.012 ± 0.041 (> 11 at 2σ).
   - [O II]/[O III] = 1.14 ± 0.08.
   - 4363/5007 = 0.075 ± 0.03, with a "floor" variant of 0.050 ± 0.035.
   - 7325/[O III] = 0.076 ± 0.019.
   - He II and He I as limits.

   Results:
   - The best χ² is **0.3 for 10 terms**. **279 of 480** truncated runs lie within Δχ² = 4.
   - The ratios bound v_s only from below: **v_s ≥ 80 km/s**. Best χ² by v_s is 5.0 at 60, 2.6 at 70, 0.7 at 80, and 0.3–0.5 for 90–200 km/s.
   - b (1–45 µG cm^3/2) is not constrained.
   - Neutral pre-shock gas is mildly disfavoured: best χ² is 1.8 at X(H⁰) = 1 against 0.3 at X(H⁰) = 0.
   - Truncation columns span log N* 17.4–18.7 and T* 1.0–3.2 × 10⁴ K.
   - Wave 2 found T* 1.7–2.2 × 10⁴ K, but that was forced by the bugged 7325 value. The corrected 7325 allows later cuts.
   - Adding SDSS [O III]/Hβ ≈ 23 as a weak central value (variant B), or the 4363 floor, changes almost nothing: χ²_min 0.8–1.0 and 280–283 acceptable runs.
2. **Complete shocks remain excluded (model).** Their χ² is ≥ 37.4 for 10 terms (≥ 41.7 in variant B, ≥ 38.3 with the 4363 floor). Several independent terms contribute:

   | Term | Contribution to χ² |
   |---|---|
   | Hα/[O III] | 9.4 |
   | He I 5876 | 8.5 |
   | Hβ | 8.4 |
   | [O II]/[O III] (2.1 vs 1.14) | 3.3 |
   | [S II] | 3.0 |
   | [N II] | 2.4 |

   Without the approximate He I term the χ² is still ≥ 29.
3. **Geometry (task 3).** The smooth 58-pc curved limb is excluded for every acceptable model.
   - It predicts σ_W/σ_E of 1.9–6.0 with my estimator (16–84 %; median 4.8). The observed value is 0.98 ± 0.05.
   - Even before applying the profile test, only 9 of 279 models pass brightness, width, FUV and age.

   **An edge-on flat sheet or ribbon works.** 35 models pass FUV, age ≤ 1.2 × 10⁵ yr, σ_W/σ_E 0.83–1.13 and a line-of-sight depth L ≤ 83 pc (the sky length at 282 pc). Their properties:

   | Quantity | Range (median) |
   |---|---|
   | v_s | 80–160 (110) km/s |
   | X(H⁰) | 0–0.3 (0) |
   | b | 5–45 (14) µG cm^3/2 |
   | log N* | 17.5–18.1 |
   | T* | 1.7–3.2 × 10⁴ K |
   | n₀ | 0.012–0.13 (0.024) cm⁻³ |
   | P/k | 1.2 × 10⁴ – 5.5 × 10⁵ (5.5 × 10⁴) K cm⁻³ |
   | L | 3.5–59 (21) pc |
   | Brightness gain G = I_obs/I_face | 3–47 (16) |
   | t* | 2–11 × 10⁴ yr |
   | σ_W/σ_E | 0.95–1.12 |

   How the closure works in this geometry:
   - The width plus the spectrum fix n₀, and so P (P ∝ 1/(W d)).
   - The brightness fixes only the depth L.
   - So the pressure match with Monogem survives as a range, but it is neither distance-specific nor a test (it is 6.1 → 4.9 × 10⁴ from 252 to 318 pc).
4. **C IV scattering in a flat sheet.** The suppression is stronger than in the curved limb, but the GALEX FUV fit is now set by the *non-resonant* lines.
   - G_CIV = L/Δ_CIV is 11–129 (median 45) and τ1548,face is 0.09–0.85.
   - S_eff is 0.12–0.57 (median 0.19, about ×5 suppression). That gives **FUV 68–143 CU/R (median 109)**, against 68–92 measured per NSNS-R (about 85–115 on the DESI R scale).
   - **New: with the full FUV line list the FUV caps v_s from above.** The list is C IV, Si IV, N IV], O IV], O III], He II and H 2γ; wave 2 used only C IV, O III] and He II.

   | v_s (km/s) | 80 | 90 | 100 | 110 | 120 | 130 | ≥ 140 |
   |---|---|---|---|---|---|---|---|
   | Median edge-on FUV (CU/R) | 85 | 98 | 114 | 132 | 144 | 163 | ≥ 206 |

   - The non-resonant floor alone reaches about 110 CU/R at v_s ≥ 140 km/s.
   - So the FUV favours **v_s ≈ 80–120 km/s**.
   - Caveat: this code's non-resonant floor is ~1.5–1.9× Raymond+88's (57 CU/R per R), so treat the cap as a soft bound of ≲ 150 km/s.
5. **Predicted FUV ridge offset.** In the flat sheet the FUV ridge leads the [O III] ridge by **+0.043…+0.103° (median +0.069°)**. The observed offset is +0.04…+0.09° (judges; 1.4–2σ). C IV alone leads by +0.09…+0.15°. [O II] trails by −0.01…−0.08°. The DESI measurement has the opposite sign, but it is flagged as not robust.
6. **Line width: what it tests has changed.**
   - An edge-on flat sheet has no line-of-sight flow, so the [O III] FWHM is thermal only: **11–13 km/s** (plus turbulence).
   - A warm 6 × 10⁴ K static layer gives 13 km/s.
   - A 5° or 10° tilt shifts the centroid by 5–10 or 10–19 km/s without broadening it. That is the size of the N–S shift.
   - Rippled sheets (10–20 pc ripples) give 20–50 km/s. The excluded 58-pc limb gave 20–26 km/s.
   - **A line width at R ≥ 30k now tests the geometry (flat vs rippled), not shock vs static layer.**
7. **NUV no longer separates the warm layer from the shock (new).**
   - The shocks' own C III] 1909, C II] 2326, Si III] and 2γ give NUV 47–99 CU/R (flat survivors 79–99).
   - The warm CIE layer gives 94 CU/R.
   - Both exceed the measured 36–57 (80 at the profile peak) by a similar ×1.5–2.
8. **Alternatives (task 2).** On the corrected optical likelihood:
   - Spectrally acceptable (p 0.27–1.0): warm CIE, NEI cooling (truncated), photoionisation with extra heating, the fossil, and conduction fronts (marginal).
   - Excluded on the spectrum alone: only photoionisation *without* extra heating (χ² 21–51 for 8 terms; p ≤ 0.02 even with the 4363 floor). The drivers are He II, Hβ/Hα and 4363.
   - The other exclusions now rest on physics: the X-ray filament (NEI cooling, TML), the ionising-photon budget (PIE), the heating source and power (warm layer, PIE + heating), the FUV floor plus width/pressure plus stratification order (conduction), and a recent hot start (fossil). See §3.
9. **North.**
   - Later truncation of the southern sheet is disfavoured. Only 3 of 35 survivors reach the northern ratios, and at a fixed n₀·G they would be 1.23–1.35× *brighter* than the south; the observed ratio is 0.36–0.45.
   - A slower shock fits. 52 north-acceptable runs span 60–180 km/s, mostly ≤ 110, at T* ≈ 10⁴ K, with [O III]/Hβ 8–26 and [O II]/[O III] 2.3–7.8.
   - It needs (P·G)_N/(P·G)_S = 0.10–3.3 (median 0.51; 0.37–0.44 for v_N = 60–90).
   - Its FUV/[O III] N/S is 0.34–1.72 (median 0.82); 43 % of pairs fall within 0.61 ± 0.15.
   - Because the north profile is a one-sided curved limb and the south is flat, G_N/G_S is unknown and P_N/P_S is not separately determined.
10. **The [O II] doublet is a common tension.** Every low-density scenario predicts 1.455–1.50, against 1.31 ± 0.07. That adds χ² = 4.3–7.4 (about 2.1–2.7σ) to all of them equally, so it does not rank them.

---

## 1. Inputs and method

### Likelihood (`scripts/w4_refit_common.py`)

Linear terms use σ² = σ_obs² + (f_model × prediction)²; log terms add in dex.

| Term | Value | f_model |
|---|---|---|
| Hβ/[O III] | 0.0122 ± 0.0414 (DESI box, null-inflated, `verify_desi_audit_extra.json`) | 0.25 |
| [O III]/Hβ, variant B only | SDSS 22.8 × / ÷ 0.30 dex | 0.10 dex |
| [O II]/[O III] | 1.14 ± 0.08 | 0.25 |
| 4363/5007 | 0.075 ± 0.030; floor variant 0.050 ± 0.035 | 0.25 |
| 7320+30/[O III] | 0.076 ± 0.019 | 0.30 |
| He II/[O III] | 0 ± 0.013 | — |
| He I/[O III] | 0 ± 0.0025 | — |
| [Ne III]/[O III] | 0.085 ± 0.02 | 0.35 |
| Hα/[O III] (energy) | 0.07 ± 0.09 | — |
| [N II]/[O III] | 0.12 ± 0.30 dex | — |
| [S II]/[O III] | 0.05 ± 0.30 dex | — |
| Doublet (separate common term) | 1.31 ± 0.07 | — |
| UV | FUV 80 ± 15 and NUV 47 ± 15 CU/R | 0.2–0.3 |
| North | [O III]/Hβ 7.5 × / ÷ 0.30 dex, plus DESI Hβ/[O III] 0 ± 0.06; [O II]/[O III] 5.5 × / ÷ 0.15 dex | — |

Model Hβ is scored both natively and ×2.3 ("R88-calibrated", the wave-2 correction), as in wave 2. "Acceptable" means within Δχ² = 4 of either variant's minimum.

### Grid (`w4_refit_grid.py`, 708 s)

- **Runs:** 532 in total.
  - 480 truncated: v_s = 60–200 km/s (13 values) × X(H⁰) = 0 / 0.05 / 0.3 / 1 × b = 1–45 (10 values).
  - 52 complete: b = 3 and 10.
- **New outputs:** cumulative UV lines from the flow's own NEI ion fractions, using the altmodel CHIANTI tables.
- **Cross-check:** for the same gas, the synthesized lines agree with the flow's own emissivities:
  - [O III]5007: 1.04–1.05×;
  - [O II]3727: 0.98×;
  - [O II]7325: 1.09×;
  - Hβ: 0.89×.

### Geometry (`w4_refit_geom.py`)

| Geometry | How n₀ is set | What brightness determines | Other details |
|---|---|---|---|
| Flat sheet | from the FWHM of the emissivity profile j(s) | L = I_obs/j_max | profile blurred to the stacked 0.41° and fitted with a split normal; C IV Monte Carlo with G_CIV = L/(2 s_CIV), LOS velocity span 0 or 1 Doppler width |
| Curved limb | wave-2 closure W (R_c = 11.9° at d) | predicted I_pk | wave-2 projection |
| Rippled sheet | from the brightness, for N_t = 1, 2, 4 tangencies of R_r = 10 or 20 pc | — | the width is set by the ripple amplitude; tangency width must be ≤ 0.28° |

All geometries use d = 282 pc (variants 252 and 318), W = 0.28° (variants 0.15 and 0.38°) and I_obs = 1.62 R.

### Other scripts

- **Alternatives (`w4_refit_alt.py`):** each altmodel family is scanned over its own parameter, and its best member is scored on the same likelihood. The parameters are T (warm CIE), log U (PIE), T_fix (PIE + heating), epoch (fossil), stop T (NEI) and P or R_cl (fronts).
- **Summary (`w4_refit_summary.py`):** applies the filters, runs the north pairing and collects the numbers. `w4_refit_widths.py` gives per-line widths and FUV lines for the flat survivors.

---

## 2. Results tables

### 2a. Which truncated shocks survive (south)

| Stage | n | v_s (km/s) | X(H⁰) | b | P/k (K cm⁻³) | Key reason for losses |
|---|---|---|---|---|---|---|
| Optical, likelihood A | 279 / 480 | 70–200 | 0–1 | 1–45 | — | v ≤ 60; X = 1 disfavoured |
| + curved 58-pc limb (I 1.2–2.5 R, FWHM ≤ 38, FUV, t* ≤ 1.2 × 10⁵ yr) | 9 | 90–110 | 0–0.3 | 7–20 | 2.6–9.6 × 10⁴ | I (215 fail), FUV (225), age (158), width (119) |
| + symmetric profile (0.98 ± 0.05) | **0** | — | — | — | — | σ_W/σ_E 1.9–6.0 (16–84 %) |
| **Flat sheet, L ≤ 83 pc** | **35** | **80–160 (110)** | 0–0.3 (0) | 5–45 (14) | **1.2 × 10⁴ – 5.5 × 10⁵ (5.5 × 10⁴)** | FUV (194 fail), age (206), L (108), asymmetry of back-weighted j(s) (136) |
| Flat sheet, L ≤ 166 pc | 36 | 70–160 | 0–0.3 | 5–45 | 0.9 × 10⁴ – 5.5 × 10⁵ | — |
| Ripple, R_r = 20 pc, N_t = 2 (tangency ≤ 0.28°, FUV) | 16 | 70–100 | — | — | 1.5–3.5 × 10⁴ | FWHM 20–40 km/s |
| Ripple, R_r = 10 pc, N_t = 2 | 16 | 70–100 | — | — | 1.8–4.4 × 10⁴ | FWHM 26–50 km/s |

**Flat-survivor model ratios**

| Ratio | Model range (median) | Observed (input) |
|---|---|---|
| [O III]/Hβ | 17–61 (40) | > 11 |
| [O II]/[O III] | 1.12–1.77 (1.18) | 1.14 |
| 4363/5007 | 0.082–0.087 | 0.06–0.09 |
| 7325/[O III] | 0.063–0.119 (0.069) | 0.071–0.082 |
| 7325/3727 | 0.056–0.067 | 0.060 |
| He II/[O III] | ≤ 0.003 | < 0.026 |
| [O III]/Hα (energy) | 6.5–21 | — |
| [N II]/Hα | 0.7–2.3 | — |

**Other flat-survivor properties**
- The emission-weighted T is 4.7–5.5 × 10⁴ K for [O III] and 3.6–4.6 × 10⁴ K for [O II]. That matches T(O II) ≈ 3.5 × 10⁴ K from the corrected 7325 (`w3_oii`).
- Proper motion at 282 pc: 60–120 mas/yr.
- **Distance:** P/k medians are 6.1 × 10⁴ at 252 pc, 5.5 × 10⁴ at 282 pc and 4.9 × 10⁴ at 318 pc; L medians are 17 / 21 / 27 pc.
- **Width:** W = 0.15° gives P ×1.9 and L ×0.29; W = 0.38° gives P ×0.74 and L ×1.85.

### 2b. UV per R (south, CU per R of [O III]5007; flat survivors unless stated)

| Quantity | Model | Observed |
|---|---|---|
| FUV, thin (all lines) | 102–354 (173) | 68–92 |
| FUV, non-resonant floor (no C IV / Si IV / N V) | 40–80 (52); all 279 acceptable: up to 129 | — |
| C IV thin → scattered | 9–39 CU/R after S_eff 0.12–0.57 | — |
| Si IV (unscattered) | 17–34 | — |
| O III] 1663 | ~21 | — |
| FUV, flat edge-on (span 0) | **68–143 (109)** | 68–92 (85–115 on the DESI scale) |
| FUV, curved limb (9 survivors) | 84–136 | — |
| FUV, ripple | 97–157 | — |
| NUV | 79–99 (85); all shocks 47–99 | 36–57 (80 at the peak) |

### 2c. Complete shocks

χ² is 37–50 for v = 130–200 and 50 to > 10³ for slower runs (Fig. a). [O III]/Hβ ≤ 5.0 and [O II]/[O III] ≥ 2.1.

**Verdict: still excluded**, with no reliance on the withdrawn [O III]/Hβ ≈ 40. The Hα and [O II]/[O III] terms alone give χ² ≈ 13 for 2 terms (derived from the per-term values in §0). The exclusion therefore does not rest on Hβ or He I.

### 2d. North

| Test | Result |
|---|---|
| North fit | χ²_min 2.0 (4 terms); 52 runs within +4: v 60–180 (median 80), X ≤ 0.3, T* ≈ 1.0 × 10⁴ K (to 2.7 × 10⁴) |
| Predicted north ratios | [O III]/Hβ 7.7–26 (14), [O II]/[O III] 2.3–7.8 (4.0), **7325/[O III] 0.11–0.50 (0.17)** (unmeasured; a prediction) |
| Later truncation of the south survivors | 3 / 35 within +4. Brightness at fixed n₀G is **×1.23–1.35** (needs ×0.36–0.45). Disfavoured |
| Slower shock (all north-acc. × south-flat pairs) | equal-P brightness ratio 0.12–3.9 (0.80); required (PG)_N/(PG)_S 0.10–3.3 (0.51); FUV/[O III] N/S 0.34–1.72 (0.82), 43 % within 0.61 ± 0.15; NUV N/S 0.59–3.7 (1.5), cf. observed 81/49 = 1.65 |

---

## 3. Alternatives re-scored (task 2)

Rows give the best member of each family on likelihood A. "Floor" means the 4363 floor variant. UV is per R; altmodel FUV values are optically thin.

| Scenario (best member) | Optical χ² / terms (p) | Floor χ² | Where it misses | What still excludes it |
|---|---|---|---|---|
| **Warm CIE layer**, T = 6.1 × 10⁴ K | 6.2 / 9 (0.72); B 10.2 | 7.1 | [O III]/Hβ 97 (2σ above SDSS 23 in B); 7325/3727 0.082 against 0.060 ± 0.015 (1.5σ); [N II]/Hα 1.6 | Viable as a spectrum. **Heating:** a face-on flux of 2.6 / 1.3 / 0.5 × 10⁻⁵ erg cm⁻² s⁻¹ for G = 10 / 20 / 50. That equals the ½ρv³ of a 110 km/s shock into n₀ = 0.024, ≈ 4 × 10⁻⁵ (derived). No identified source, and the layer sits on the thermally unstable rise of Λ(T) (altmodel). A single-T layer gives no FUV–[O III] offset (observed +0.04…+0.09°, 1.4–2σ). **NUV no longer counts against it** (94 against shocks' 79–99). FWHM 13 km/s, the same as a flat-sheet shock |
| **PIE + extra heating** to 4 × 10⁴ K | 1.0–2.4 / 9 (≥ 0.98) | 1.4–2.8 | none | Needs 78–351× the photo-heating *and* the ionisation-parameter (U) deficit below. It is a heated layer plus a missing ionising source |
| **NEI cooling of hot gas / TML** (isobaric, stopped at 2 × 10⁴ K) | 5.2–5.5 / 9 (0.8) | 6.1–6.4 | [Ne III] (GS07 Ne) | **X-ray:** cooling from 1.6 MK predicts an X-ray filament 1.9 × 10³ – 2.3 × 10⁴× the Monogem surface brightness. This ratio is independent of G; none is seen. The TML needs an inflow of 377 × (18.7/G) km/s: 141 km/s even at G = 50, near transonic (derived). Starting at 3 × 10⁵ K is the cooling zone of a shock. Truncation at 2 × 10⁴ K needs a cause. FUV 227–245 thin |
| **Conduction front**, static P/k 2 × 10⁶ / evaporating | 11.1 / 9 (0.27); 10.4 / 9 (0.32) | 12.2; 11.6 | [O II]/[O III] 1.96–2.54 (2.2–3σ); 7325 0.14–0.18; [Ne III] | **FUV floor:** the non-resonant part alone is 179–213 CU/R, 2–2.7× the measurement even with C IV and Si IV fully scattered. **Width/pressure:** the front is 86 pc thick at P/k 1.5 × 10⁴ (∝ 1/P), so the 1.38-pc width needs P/k ≈ 9 × 10⁵, about 60× Monogem. Then one front gives 0.17 R face-on (needs G ≈ 9), but the driving ~10⁶ K gas would be X-ray bright. **Order:** hot side west puts the FUV *west* of [O III], opposite to the observed sign |
| **PIE**, best spectrum (0.14 keV hot gas / 150 kK / WD 81 kK / precursor) | 21–29 / 8 (0.007–3 × 10⁻⁴); ISRF 51 | 18–25 (p ≤ 0.02) | He II 0.03–0.06 (χ² 3–8); Hβ, Hα (2.4–7); 4363 ≤ 0.012 (4.4–5.3; ~1–1.5 with the floor); FUV 49–53 | Spectrally excluded at p ≤ 0.02, plus the U budget. The WD J0736 shortfall is 27× at L = 25.6 pc; scaling as L^−½ it is ~18–30× for the flat-sheet L of 21–59 pc (derived). Monogem EUV is short by 630–3400× |
| **Fossil**, hot start T₀ = 4 × 10⁴ K | 0.3 / 5 (1.0) | 0.2 | none (5 terms only) | Since `w3_oii` (n_e not ~100), the e-fold time at n_H = 0.5 is 6.9 × 10³ yr, so the 20-yr constancy **no longer excludes it**. It still needs gas heated to ≥ 4 × 10⁴ K within the last 1.5–4.2 × 10³ yr along 17°: a recent shock under another name. It predicts no FUV offset; its FUV was not computed |
| Fossil, T₀ = 2 × 10⁴ K (photoionised start) | 3.4–3.8 / 5 (0.6) | 1.4–1.8 | 4363 0.028 | Fails only on 4363 (2.4 nominal, ~1 with the floor). A photoionised start also needs the missing ionising source |

**Verdict per scenario**

| Scenario | Verdict | Basis |
|---|---|---|
| Complete shocks | ruled out | spectrum |
| PIE (no heating) | ruled out | spectrum and U budget |
| NEI cooling / TML | ruled out | X-ray |
| Conduction front | ruled out | FUV floor, width/pressure, order |
| Warm heated layer and PIE + heating | disfavoured | a shock-sized heating power with no source; no stratification. Not excluded by the spectrum or NUV |
| Fossil | disfavoured | needs a recent shock; equivalent to the shock in practice |

---

## 4. Best discriminating measurements, with predictions per surviving scenario

The surviving scenarios are:
- **FS:** flat-sheet truncated shock (35 models).
- **RS:** rippled-sheet shock.
- **WL:** warm heated layer at 6 × 10⁴ K, including PIE + heating at 4 × 10⁴ K.
- **FO:** hot-start fossil at 4 × 10⁴ K.

The curved limb and complete shocks are out.

| # | Measurement (precision needed) | FS | RS | WL | FO |
|---|---|---|---|---|---|
| D1 | **Hβ depth**: Hβ/[O III] to ±0.005 (≈ 0.008 R at 1.62 R, 1σ) | [O III]/Hβ 17–61 (40); Hβ 0.027–0.095 R | same family as FS (70–100 km/s: 17–60) | **97** (Hβ 0.017 R); PIE + heat 45 | 39–40 |
| D2 | **FUV spectrum** (FIMS/SPEAR archive or UVEX-class; C IV 1548/1551, Si IV, O III] 1663) | C IV 0.5–2.3 × 10⁴ LU at 1.62 R (9–39 CU/R); **1548/1551 → 1.0–1.3** (τ_edge ≫ 1); C IV/Si IV 0.4–1.6 (CU) | C IV larger (S 0.3–0.6): FUV 97–157 CU/R | C IV/Si IV 0.6 (8 / 14 CU/R); **1548/1551 ≈ 2** (thin); O III] dominant | no C IV, Si IV (T ≤ 4 × 10⁴ K) |
| D3 | **FUV–[O III] ridge offset** (GALEX 5″, injection/placement nulls; ±0.015°) | **+0.043…+0.103°** east; C IV alone +0.09…+0.15; [O II] −0.01…−0.08 | ~+0.1° per tangency, diluted | 0 | 0 |
| D4 | **[O III]5007 width** at R ≥ 30 k (σ_inst ≤ 4 km/s), plus the centroid | FWHM 11–13 (thermal); centroid = pre-shock gas ± 5–19 km/s for a 5–10° tilt; [O II]/[O III] width 0.81–0.98 | **20–50 km/s**, flat-topped | 13 km/s; [O II]/[O III] width 1.0 | 11–12 km/s |
| D5 | **[O II]7325/3727** OH-clean to ±0.007 | 0.056–0.067 (0.058) | ~0.06 | **0.082** | ~0.066 |
| D6 | **4363/5007** with Hg resolved, ±0.01 | 0.082–0.087 | 0.082–0.091 | 0.098 (PIE + heat 0.075) | 0.062 |
| D7 | **Proper motion** of the east edge (≥ 10 yr baseline, ≤ 0.3″ registration) | v/4.74d = **60–120 mas/yr** at 282 pc (0.6–1.2″ per decade) outward | 50–80 mas/yr | ≈ 0 | ≈ 0 |
| D8 | **North 7325/[O III]** (new prediction of the slower-shock north) | 0.11–0.50 (0.17) | — | single-T: ~0.1 (like the south) | — |
| D9 | **[O II] doublet** to ±0.03 | 1.45–1.47 | 1.45–1.47 | 1.46 | 1.46 (n ~0.5) |

D9 does not discriminate between scenarios. A firm ≤ 1.35 would falsify all of them.

**Single best discriminators**
- **D3 and D7 against WL/FO:** stratification and motion. D3 is doable now with existing GALEX data.
- **D1 against WL:** a factor of 2.5 in Hβ.
- **D2's C IV doublet ratio against WL:** optically thick against thin.
- **D4 for FS against RS,** and for the tilt behind the N–S velocity shift.
- **D7 plus a v_s from D2/D1:** an expansion-parallax distance d = v_s/(4.74 μ).

---

## 5. Changes relative to wave 2 and wave 3

| Item | Before | Now |
|---|---|---|
| Truncated-shock v_s | 100–150 (best 110–140), T* 1.7–2.2 × 10⁴ K | Optical ≥ 80 only. FUV with the full line list favours **80–120 (≲ 150)**; flat survivors 80–160. T* 1.7–3.2 × 10⁴ K |
| X(H⁰) | ≲ 0.3 | ≲ 0.3 (X = 1 is still Δχ² +1.5, and fails the FUV) |
| b | 7–14 (set by a 58-pc closure) | **unconstrained** (5–45): the flat closure has no curvature term |
| n₀ / P/k | 0.012–0.042; 3–11 × 10⁴ | flat: 0.012–0.13 (median 0.024); 1.2 × 10⁴ – 5.5 × 10⁵ (5.5 × 10⁴). ∝ 1/d; P ∝ 1/W |
| Limb gain G | 20–21 (curved) | flat G = I_obs/I_face = 3–47 (16), via L = 3.5–59 pc (21) |
| FUV per R | 52–126 (C IV + O III] + He II only) | 68–143 flat, full line list; the non-resonant floor bounds v_s |
| [O III] FWHM | 23–38 km/s (curved limb) | **11–13 km/s flat** (thermal); 20–50 rippled |
| "Two failures" (doublet, 7325 ×1.8) | tensions | 7325 is consistent (model 0.063–0.119 against 0.071–0.082). The doublet is a common 2.1–2.7σ term |
| Altmodel two-temperature, n_e ~ 100 | excluded single-T models | gone. The warm layer, PIE + heating and fossil are spectrally viable; exclusion is by budget and source |
| NUV against the warm layer | ×2 | not discriminating (shock 79–99 vs warm 94) |

---

## 6. Limitations

- **Wave-2 code limitations carry over:**
  - no radiative transfer, so Hβ at intermediate columns is ~2× below R88 (handled by the ×2.3 variant);
  - truncated runs end at 10⁴ K, so no recombination zone;
  - Fe cooling is in CIE;
  - b is gridded;
  - S⁺/S²⁺ have no dielectronic recombination ([S II] is unreliable).
- **UV lines are from CHIANTI 7.1.4 on the flow's ion fractions.** Validation: the optical lines agree within 4–11 %. The non-resonant FUV floor is ~1.5–1.9× R88's; the NUV is ~1.4× R88's 47–63. **The FUV cap on v_s and the NUV comparison carry that systematic.** Si IV (a resonance doublet) was not scattered, so FUV is an upper side.
- **Flat-sheet model:** exactly edge-on, uniform along L, emissivity profile = j(s). Scattering escape through the faces is assumed lost. No turbulence is in the line widths.
- **Rippled sheet:** simplified as N_t independent tangencies. The overlap factor is folded into N_t.
- **Asymmetry statistic:** my split-normal estimator differs from `judge_skeptic`'s; curved-limb values are 1.9–6.0 here against 2.8–3.0 there. Both are ≫ 0.98.
- **Alternatives:** read from the altmodel JSONs; not recomputed. Warm-CIE UV comes from the altmodel "own" CIE list (FUV 62); its window table gives 86–92. Fossil and PIE FUV were not computed. The altmodel geometry-dependent budgets were rescaled by G (derived).
- **Observational inputs are canonical, not re-measured.** UV per R carries the NSNS-vs-DESI scale factor (×1.0–1.35).
- **The ±0.015 OH systematic on 7325** dominates the D5 error budget.
- **The north pairing** compares mixed geometries (north curved limb, south flat). Only (P·G) ratios are constrained.

## 7. Files (tag `refit`)

**Scripts** (all run with `~/.venvs/astro313/bin/python -I`, in this order)
- `scripts/w4_refit_grid.py`
- `w4_refit_fit.py`
- `w4_refit_alt.py`
- `w4_refit_geom.py`
- `w4_refit_summary.py`
- `w4_refit_widths.py`
- `w4_refit_fig.py`
- `w4_refit_common.py` (likelihood)

**Data**
- `data/derived/w4_refit_tracks.npz` (532 runs) and `w4_refit_grid.json`, `w4_refit_grid_log.txt`
- `w4_refit_fit.json` (per-run χ² for 8 variants, UV)
- `w4_refit_geom.json` (279 models × 3 geometries, d/W variants), `w4_refit_geom_log.txt`
- `w4_refit_alt.json`
- `w4_refit_summary.json`
- `w4_refit_widths.json`

**Figure:** `figures/w4_refit_summary.png` (200 dpi, inspected).
- (a) χ² against v_s, truncated and complete.
- (b) Hβ/[O III] against [O II]/[O III] for shocks, alternatives and the S/N boxes.
- (c) Best-model cross-profiles, flat against curved.
- (d) FUV/NUV per R against v_s.
- (e) Optical χ² of the alternatives.
- (f) Predicted FWHM by geometry.

No downloads were made. No earlier files were modified.
