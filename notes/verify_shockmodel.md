# verify_shockmodel: radiative-shock models tested against the Silk's spectrum

Wave-2 modelling agent, tag `shockmodel`, 2026-11-16.

All model numbers come from my own code (`scripts/verify_shockmodel_*.py`). I read no wave-1 scripts and reused none of them. Observational inputs are the **wave-2 revised values**, taken from `notes/verify_desi.md`, `verify_sdss.md`, `verify_xrayfuv.md` and `verify_hiwham.md`, and I mark them as such. Wave 2 superseded several wave-1 numbers quoted in the task:
- [O III]/Hβ is now 46 (>28) in the DESI south and 23 in the SDSS south.
- [O III]4363/5007 is now 0.083–0.092, and it is not Hg.
- R23 is now ~76–130.
- The ridge peak is now ~1.6 R (1.2–2.3).
- GALEX FUV is now a 92 ± 10 CU per R detection, not an upper limit.
- The HI wall was refuted.

Labels used below: **measured** (by another agent; cited), **model** (my code), **derived** (model plus a stated assumption), **literature**.

---

## 1. Bottom line

| Scenario | Verdict | Why |
|---|---|---|
| (1) Complete steady radiative shocks, 80–200 km/s | **Ruled out** | When cooling completes, [O III]/Hβ is 2.4–5.0 in my code and 1.2–4.5 in every published grid (R79, CR85, DS96, Allen+08). The observed value is ≥24 (DESI 97.5 %) to 46. A complete shock would also add an Hα filament and cooled HI, and neither is seen. |
| (2) Incomplete, column-truncated shock (Cox & Raymond 1985; Raymond+88) | **Favoured, viable** | Truncated just before H recombination (T* ≈ 1.7–2.2×10⁴ K, log N* ≈ 17.6–17.9), the model reproduces [O III]/Hβ, [O II]/[O III], 4363/5007, [Ne III], [N II], [S II] and He II. With the filament's 0.28° width it also gives 1.3–4.0 R at the Monogem pressure, the observed line width and centroid, and the GALEX FUV once edge-on C IV resonance scattering is included. **Two failures:** the [O II] 3729/3726 density, and [O II]7325 (×1.8). |
| (3) Shock plus photoionized precursor | **Precursor as the [O III] emitter: ruled out.** Precursor as a pre-ioniser: viable, and part of (2) | Photoionized gas cannot give T_e([O III]) ≈ 4–6×10⁴ K or R23 ≳ 80. At n₀ ≈ 0.02–0.04 cm⁻³ the precursor neither recombines nor cools within its flow time (t_rec ≈ 5 Myr, t_cool ≈ 6 Myr), so it is dark. Its EUV (0.2–1.6 ionizing photons per incoming H at 100–140 km/s) explains why the fits prefer pre-ionized gas (X(H⁰) ≲ 0.3). |

**Best physical picture.** A v_s ≈ 100–150 km/s (best 110–140) shock driven by the Monogem hot gas into modestly dense gas:
- pre-shock n₀ ≈ 0.012–0.042 cm⁻³, i.e. 2–8× the Monogem ambient 5×10⁻³;
- B₀ ≈ 0.8–2.9 µG (b = B₀/√n₀ ≈ 7–14 µG cm^{3/2});
- ram pressure P/k ≈ 3–11×10⁴ K cm⁻³;
- seen edge-on on a ~58 pc radius sheet (limb-brightening G ≈ 20);
- the shock entered this gas t* ≈ 4–10×10⁴ yr ago, which is less than its full cooling time of 0.7–1.4×10⁵ yr and comparable to the remnant age of 1.1×10⁵ yr.

So truncation is the natural state, not a fine-tuning.

## 2. Results table: claim → verdict → my value

| # | Claim tested | Verdict | My value / reason |
|---|---|---|---|
| 1 | Complete radiative shocks 80–200 km/s can produce the spectrum | **REFUTED** | Complete-flow [O III]/Hβ: my code 2.4–5.0 (Tfloor 8000 K, b = 3). Published: R79 2.8–4.5, CR85 2.5–3.4, DS96 1.2–1.9, Allen+08 1.2–3.2. The observed value is ≥24. χ² at completion is 41–4200, against 4.9 for truncated models. |
| 2 | Truncated shocks reproduce [O III]/Hβ ≈ 20 (now 46) and R23 ≈ 58 (now ~80–130) with weak He II | **CONFIRMED** (with revised inputs) | Best χ² is 4.9–5.2 for 11 terms, about 9 dof (native Hβ). Consistent models have [O III]/Hβ 10–63 (median 49) and R23 23–155 (median 117). He II/[O III] is 0.001–0.005 (limit 0.026). |
| 2a | Which v_s? | REVISED → **v_s ≳ 100 km/s** from the optical ratios; **100–150** with all tests | χ²_min by v: 70: 14.9; 80: 9.6; 90: 7.4; 100: 6.2; 110–200: 4.9–5.5. The optical ratios do not bound v from above. Line width (≤ 38 km/s) and brightness do. FUV gives v ≳ 100–110 (the C IV-weak slow shocks fall below the measured FUV). |
| 2b | Which truncation? | **log N* = 17.6–17.9 cm⁻², T* = 1.7–2.2×10⁴ K** | This cut keeps 93–99 % of the [O III] the flow makes down to 10⁴ K, but stops before the H recombination zone. The [O III] is emitted at a 5007-weighted ⟨T⟩ = 5.5–6.4×10⁴ K, compared with the measured T_e([O III]) 4.6 (3.1–5.3)×10⁴ K (verify_desi). 4363/5007 in the model is 0.086–0.094; measured 0.083–0.092. |
| 2c | Pre-shock ionization | **X(H⁰) ≲ 0.3** | χ²_min by X (native / R88-calibrated Hβ): 0: 5.7/4.6; 0.05: 4.9/4.9; 0.3: 5.0/6.7; 1.0: 6.1/9.1. Neutral pre-shock gas makes Balmer emission from collisional excitation at the front, which gives [O III]/Hβ ≲ 10–20. |
| 3 | (b) 2–4 R reachable at d ≈ 280 pc under Monogem pressure balance | **CONFIRMED** (b ≈ 7–14 only) | I set n₀ from the observed width 0.28° (R_c = 58 pc). Then I_pk = 1.3–4.0 R (median 2.3) and the implied ρ₀v²/k = 3.1×10⁴–1.1×10⁵. The Monogem estimate is P_s/k = 9.2 n₀T = 3.7×10⁴–1.6×10⁵ (Knies+24 n₀f^−½ 3.1–8.6×10⁻³, kT 0.11–0.17 keV). Ranges: d = 250 pc gives 1.5–4.5 R; 320 pc gives 1.1–3.5 R; a 15 pc ripple radius gives 0.7–2.0 R. b ≤ 5 is too faint or under-pressured; b ≥ 20 is too bright or over-pressured (≥ 2×10⁵). |
| 3a | n₀ from pressure balance with Monogem | **CONFIRMED** | n₀ = 0.012–0.042 cm⁻³. The alternative closure, n₀ fixed by P/k = 10⁵, gives 2.8–8.2 R and widths of 0.09–0.30°, still within a factor ~2. P/k = 3×10⁴ gives 0.5–1.4 R. P/k = 4×10⁵ gives 22–66 R and widths of 0.02–0.08°, which is excluded. |
| 3b | n₀ from the "HI wall" (N_HI 5–9×10¹⁹ over 0.9°) | **REFUTED** as the pre-shock medium | verify_hiwham showed the wall is a finder artefact. Even if it were real, n ≈ 3.7–6.6 cm⁻³ would allow only v_s ≈ 10–25 km/s at P/k = 10⁵–4×10⁵, so no [O III]. |
| 4 | (c) Is truncation natural, given cooling time vs Monogem age (1.1×10⁵ yr)? | **CONFIRMED** | t* = N*/(n₀v_s) = 3.9×10⁴–1.04×10⁵ yr, against t_cool(→10⁴ K) = 0.68–1.4×10⁵ yr. The shock has been in this gas for a large fraction of the remnant age but less than its cooling time. Penetration v_s·t* ≈ 5–14 pc; the emitting layer is s* ≈ 1.0–1.7 pc. |
| 5 | (d) FUV/C IV versus GALEX | **CONFIRMED** (resonance scattering required) | Optically thin, the model gives 168–403 CU per R(5007): C IV 145–377, O III] 22–25 (cf. xrayfuv's independent 22.6 for O III]+He II), He II 0.3–1.8. That is 2–4× the measured 92 ± 10 (south). The edge-on Monte Carlo with face-on τ₁₅₄₈ = 0.42–1.13 (N(C IV) = 1.1–3.5×10¹³ cm⁻²), G_CIV ≈ 19–35 and b = 11–14 km/s gives S₁₅₄₈ = 0.16–0.27, so 59–143 CU/R, consistent. Omitted lines (N IV], Si IV, C II, N III], 2γ) would add ≲ 35 CU/R (xrayfuv's R88 floor 57 minus 22.6). The old "UL 0.4–0.9 R at C IV" is superseded by the detection. |
| 6 | (e) V_LSR ≈ 0 and intrinsic FWHM ≲ 33–40 km/s | **CONFIRMED** | At the tangent the flow (u = 0.79–0.89 v_s) projects symmetrically, so the centroid sits at the pre-shock gas velocity. The local HI and WIM are at V_LSR +1…+2 (verify_hiwham), against [O III] V_LSR ≈ −2…+2. The profile is flat-topped or double-horned, |v_los| ≤ 12–22 km/s, plus thermal σ_O ≈ 5–6 km/s. FWHM 23–38 km/s for consistent models; ≥ 39–48 km/s for v ≥ 150–160, which is disfavoured. |
| 7 | (f) Proper motion, and is it testable? | **Derived** | μ = v_s/(4.74 d) = 75–113 mas/yr at 280 pc, i.e. 1.5–2.3″ in 20 yr, outward (E/NE, away from Monogem). The projected rise behind the front is ~0.15° wide, so the smooth ridge cannot show it. It needs arcsecond-sharp substructure in deep images ≥ 10–20 yr apart (e.g. the 2024–25 Fesen/NSNS frames re-observed ~2035–2045). The ratio evolution (timescale ≥ 10⁴ yr) is not testable. |
| 8 | [O II]3727/[O III] 1.1–1.5 | **CONFIRMED (south)** | Model 0.98–1.16, against DESI-south 1.11 (1.07–1.18). It is low against DESI-all 1.37 and SDSS-south 1.84 (1.46–2.37). |
| 9 | [O II]3729/3726 = 1.30 → n_e ≈ 110 (40–300) | **INCONSISTENT with every pressure-balanced shock** (new tension) | Model n_e in the [O II] zone is 0.14–0.25 cm⁻³, so the ratio is 1.467 (low-density limit). n_e ≈ 150 at ~1.5×10⁴ K needs P/k ≈ 5×10⁶, 50–150× the Monogem pressure. Such a shock would be ≳ 10² R bright, or the [O II] gas would need a volume filling factor ~10⁻³. Either the doublet ratio carries a ~10 % systematic (SDSS-all gives 1.13, and the SDSS [S II] ON−OFF ratio of 2.5 is unphysical, so such systematics exist at this S/N), or the [O II] is not from this shock. |
| 10 | [O II]7320+30/[O III] = 0.034 (DESI) | **Model ×1.8 too high** | Model 0.059–0.068. The model's [O II] forms at T ≈ 2×10⁴ K, where 7325/3727 is ~0.06; the measurement (0.026) implies T([O II]) ≈ 1.4×10⁴ K. |
| 11 | [N II]5755/6583 = 0.14 (0.06–0.26) | CONFIRMED (within errors) | Model 0.078–0.085, so T([N II]) ≳ 2–3×10⁴ K in both. |
| 12 | South [O III]/Hα ≈ 5; [N II]/Hα 0.5–1; [S II]/Hα 0.3–0.5 | CONFIRMED (broad) | Model [O III]/Hα 3.8–22 (median 17; DESI-south 13.7, SDSS-south 2.7 (1.9–4.5)). [N II]/Hα 0.39–2.3, [S II]/Hα 0.17–1.05, Hα/Hβ 2.6–2.9 (DESI 3.9, 2.2–6.1). |
| 13 | No X-ray, Hα or HI filament | **CONFIRMED** (consistent) | T_s = 1.2–5.2×10⁵ K, so nothing above 0.2 keV. Predicted I(Hα) is 0.09–1.1 R (median 0.2) against the 2σ limits of 1.0–1.25 R (WHAM template) and 0.43 R (F03, south). The two X = 1 models at 0.6–1.1 R are marginal. No cooled HI. Pre-shock X ≤ 0.3 at n₀ ≈ 0.03 gives N(HI) ≲ 6×10¹⁷ cm⁻² along a 20 pc chord, far below the 1–1.75×10¹⁹ limit. |
| 14 | Soft X-rays brighter just inside (Monogem, kT ≈ 0.14 keV) | CONFIRMED (consistent) | The required driving pressure (3–11×10⁴) lies inside the Monogem post-shock range. The density contrast with the Monogem ambient is n₀/n_amb ≈ 2–8, so v_s ≈ v_bw/√(n₀/n_amb) ≈ 100–200 km/s for v_bw ≈ 290 km/s (kT_s ≈ 0.1 keV). |
| 15 | North (Dec ≥ 20): [O III]/Hβ ≈ 6.4, [O II]/[O III] 4.4–7.3 | **Model-consistent** (not unique) | Along the same tracks this needs either a slower shock (70–80 km/s, earlier cut) or a later cut at ~10⁴ K with recombination starting (log N ≈ 18.0–18.2). 24 tracks reach χ² < 2 on these two ratios. |

## 3. The model and its validation

### 3.1 Code (`verify_shockmodel_atomic.py`, `verify_shockmodel_flow.py`)

**Flow and energy.**
- Plane-parallel steady flow behind a perpendicular MHD shock: exact mass, momentum and energy fluxes, with B ∝ ρ (transverse field frozen in).
- Rankine–Hugoniot jump from the same relations; T_e = T_i immediately.
- Cooling: n_e Σ n_ion Λ_ion(T) with the Gnat & Ferland (2012) ion-by-ion efficiencies. Fe is in CIE, cooling only.

**Ionization (non-equilibrium, all stages of H, He, C, N, O, Ne, Mg, Si, S).**
- Voronov (1997) collisional ionization, parsed from Verner's `cfit.f`.
- Badnell RR and DR fits, with Verner `rrfit.f` as fallback. S⁺ and S²⁺ have no DR.
- H charge transfer from Kingdon & Ferland, Stancil+99 (O⁺) and Barragán+06 (O²⁺/N²⁺), parsed from the Cloudy c17 source listing.

**Line emission.**
- Forbidden lines from 5- or 6-level atoms with PyNeb atomic-data files: O III AK99/FFT04-SZ00; O II Kal09/FFT04; N II T11; S II TZ10; Ne III McL11; C IV AK04/WFD96.
- H I and He II recombination from the Péquignot+91 case-B fits.
- H I collisional Balmer excitation from CHIANTI h_1 (Anderson+ Υ) with case-A cascades.
- He I 5876: approximate, ε = 1.40 (T/10⁴)^−0.15 ε_Hβ.

**Approximations.**
- No radiative transfer: no self-photoionization of the cooling gas, and no self-consistent precursor (the pre-shock state is an input).
- Optically thin lines, except the separate C IV Monte Carlo.
- Truncated runs stop at 10⁴ K. Complete runs hold T at 8000 K until x(H⁺) < 0.01.

**Grid.**
- 404 truncated runs: v = 70–200 km/s × X(H⁰) = 0/0.05/0.3/1 × b = 1–30 µG cm^{3/2}.
- 8 runs with refractory-depleted abundances; 12 complete runs.
- Every swept column along a run is a candidate truncation.
- The flow is scale-free in n₀ at fixed b (s, t ∝ 1/n₀, flux ∝ n₀), apart from the [O II] doublet. Tracks are computed at n₀ = 0.1 and rescaled.

### 3.2 Validation (`verify_shockmodel_validate.py`, `verify_shockmodel_r88.py`; `figures/verify_shockmodel_validation.png`)

**CIE.** Ion fractions agree with Gnat & Sternberg 2007 (VizieR, Asplund05) to ≲ 0.02 absolute for He, C, N and O at 2×10⁴–10⁶ K. Example, O at 1.01×10⁵ K: 0.034/0.646/0.318 vs 0.034/0.64/0.32.

**Isobaric non-equilibrium cooling from 5×10⁶ K.**
- O ion fractions match GS07 within ≈ 0.03.
- Λ is 5–16 % low at 2×10⁴–2×10⁵ K and 30–40 % low at 3–5×10⁵ K (Fe in CIE).

**Complete shocks.**
- **Allen+08 Table 6** (100 km/s, n = 1, B = 3.23 µG):
  - Hβ flux log −5.300 vs −5.317;
  - [O III]/Hβ 5.2 vs 3.2; [O II] 11.6 vs 10.7; C IV 13.1 vs 15.0; Hα 2.96 vs 3.21; 4363 0.43 vs 0.29.
  - At 200 km/s my complete model fails: Hβ is 3.7× low and [O III]/Hβ is 6.4 vs 1.2. The photoionized recombination zone of fast shocks is not modelled. This is irrelevant for truncated flows, but my "complete" values at v ≳ 150 should not be trusted.
- **R79 (n₀ = 10)**, [O III]5007/Hβ (mine vs R79):

  | v_s (km/s) | 100 | 120 | 141 | 160 | 200 |
  |---|---|---|---|---|---|
  | mine | 3.5 | 4.3 | 5.4 | 5.4 | 5.0 |
  | R79 | 4.1 | 4.3 | 4.5 | 3.7 | 2.8 |

**Truncated shocks: Raymond+88 (R88).** I read R88 Tables 5 and 7 myself from the ADS scan. They agree with the xrayfuv transcription.
- **SE** (140 km/s, n₀ 4.4, B 11 µG, log N 18.20) and **M** (120 km/s, n₀ 6.0, B 15 µG, log N 18.25):
  - F(Hβ) face-on: 0.073/0.096 vs 0.092/0.107 ×10⁻⁵ erg cm⁻² s⁻¹ sr⁻¹;
  - absolute [O III] 1.22×/1.11× R88; C IV/Hβ 69/60 vs 68/43; Hα/Hβ 2.89 vs 3.03;
  - [O III]/Hβ 19.8/13.0 vs 12.9/10.5; [O II]/Hβ 37/32 vs 17/17;
  - 4363, O III]1664 and [O II]7325 are 1.4–2.8× high.
  - My code reaches 10⁴ K before log N = 18.2 and then sits on a 10⁴ K plateau, which inflates [O II].
- **R88 Fig. 8a** (80 km/s, n₀ 8.6, B 18 µG):
  - the peak *positions* in log N of [O III]/Hβ, [O II]/Hβ and O III]/Hβ match within 0.03–0.08 dex;
  - cumulative [O III] reaches 50 % at log N 17.13 vs ≈ 17.05;
  - the complete end values match within 0.05–0.13 dex;
  - but the **peak line/Hβ ratios are 0.28–0.40 dex higher**.

**Interpretation.** The column scale and the [O III], C IV and Hβ totals are right. At intermediate columns my Hβ is ~2–2.5× lower than R88's. R88's truncated spectra are slices of complete models, so they include photoionization of H in the cooling zone by the flow's own EUV. I therefore carry a 0.32 dex model term on every Hβ or Hα ratio, and I repeat every fit with model Hβ ×2.3 ("R88-calibrated"). The conclusions survive: global χ²_min is 4.6, the preferred X moves to 0–0.05, and log N* moves by ≲ 0.1 dex.

## 4. Fitting (task a): `verify_shockmodel_fit.py`

**Observed constraints** (energy ratios, south-weighted, wave-2 values), each as value (σ_obs; σ_model) in dex:

| Ratio | Value | σ_obs | σ_model |
|---|---|---|---|
| [O III]/Hβ | 40 | 0.20 | 0.32 |
| [O II]/[O III] | 1.35 | 0.12 | 0.10 |
| [O III]/Hα | 7 | 0.30 | 0.32 |
| [N II]6583/[O III] | 0.12 | 0.25 | 0.14 |
| [S II]/[O III] | 0.05 | 0.30 | 0.20 |
| [Ne III]/[O III] | 0.10 | 0.12 | 0.14 |
| 4363/5007 | 0.087 | 0.10 | 0.10 |
| [O II]7325/[O III] | 0.034 | 0.10 | 0.10 |
| [N II]5755/6583 | 0.13 | 0.25 | 0.10 |

Upper limits: He II/[O III] < 0.026 and He I 5876/[O III] < 0.005, each with a one-sided 0.1 dex penalty.

**Results.**
- 238 truncated runs are acceptable (χ² ≤ χ²_min + 4).
- The best fits all cut the flow at T ≈ 1.7–2.5×10⁴ K, where O²⁺ has started to recombine (so [O II] ≈ [O III]) but H has not.
- The dominant χ² residual is [O II]7325, at 2.5–4.
- Figure `figures/verify_shockmodel_tracks.png` shows the cumulative ratios against N for 100–160 km/s, a neutral-pre-shock case, the published complete values and R88.

## 5. Brightness, pressure and geometry (b): `verify_shockmodel_analysis.py`

**Projection.** I project the model emissivity j(s) on a curved sheet (thin-shell chord z² = 2R_c(h − s)) and integrate along the line of sight. The ridge sits ≈ 0.15° behind the front. The profile is asymmetric: a sharp outer (east) edge and a long 1/√h tail toward Monogem. G = I_pk/I_face = 20.2–20.8 for R_c = 58 pc (11.9° at 280 pc). Face-on intensities are 0.06–0.20 R, and the 5007 yield is η = 0.35–0.53 photons per H nucleus swept.

**Two closures.**
- **(W) n₀ from the observed 0.28° FWHM.** Width and brightness then fix n₀ without any pressure assumption. For b = 7–14 the implied ρ₀v² equals the Monogem pressure, which is a non-trivial coincidence and the main support for scenario 2. Map in `figures/verify_shockmodel_closure.png`.
- **(P) n₀ from P/k = 10⁵.** This gives brightness within ×2 and widths of 0.09–0.30°.

**Monogem pressure** (derived): P_s/k = n_tot T = 9.2 n₀ T_X, using Knies+24 Table 4 (n₀f^−½ = 3.1–7.2×10⁻³ single-SNR or 5.5–8.6×10⁻³ two-SNR; kT 0.11–0.17 keV). This gives 3.7×10⁴–1.6×10⁵ K cm⁻³, ×f^−½. A cloud shock has ρ_c v_c² ≈ (1–2.5)(4/3)P_s. The fits sit at the low end, so the Silk does not need reflected-shock enhancement.

## 6. FUV and resonance scattering (d): `verify_shockmodel_rt.py`

**Method.**
- A 2-D Monte Carlo: a ribbon of thickness D and line-of-sight length G·D, uniform emission, complete redistribution over a Doppler profile, isotropic scattering, no destruction (E(B−V) to 280 pc ≈ 0.006).
- It includes the line-of-sight velocity gradient of the curved flow (span 1.2–2.6 Doppler widths).
- Peel-off estimator toward the observer.
- Sanity limits: S → 1 as τ → 0; at τ_face = 1 and G = 30, S = 0.14 (static) or 0.20 (with gradient).

**C IV layer.** The C IV layer is thinner and nearer the front (⟨T⟩ ≈ 0.9–1.4×10⁵ K), so G_CIV ≈ 19–35. That gives τ_edge ≈ 10–40.

**Result.** The predicted FUV per R(5007) is 59–143 CU (E(B−V) 0) and 52–126 CU (E(B−V) 0.03). Measured: 92 ± 10.

**Prediction.** The C IV/FUV ridge should lie ≈ 0.1° *outside* (east of) the [O III] ridge (`figures/verify_shockmodel_fuv.png`, middle). That is testable with GALEX at 5″ resolution.

## 7. Kinematics (e) and proper motion (f)

**Line profile at the projected peak.** v_los = u(s)·z/R_c along the chord, weighted by j(s), plus thermal O broadening. The result is flat-topped, |v_los,max| = 12–22 km/s, FWHM 23–38 km/s. Convolved with the DESI LSF (98 km/s) it is unresolved, consistent with DESI's ≤ 33 km/s (26–38).

**Centroid.** The centroid equals the pre-shock velocity. A face-on or oblique part of the sheet would be shifted by up to ±(0.8–0.9)v_s ≈ ±80–130 km/s. The single observed velocity (DESI and SDSS within ~3 km/s along the south) therefore means we only see the tangent limb, as limb-brightening requires.

**Proper motion.** 75–113 mas/yr at 280 pc (66–127 for 250–320 pc).

## 8. Precursor (scenario 3): `verify_shockmodel_precursor.py`

**Ionizing photon budget.**
- Energy radiated above 5×10⁴ K: 27 / 51 / 80 eV per H at 100 / 120 / 140 km/s (total ½μv² = 71 / 102 / 139 eV).
- Half goes upstream, of which 50–80 % is ionizing at ⟨hν⟩ 40–20 eV.
- That gives φ = 0.17–0.54 / 0.32–1.0 / 0.50–1.6 ionizing photons per incoming H (derived, order of magnitude).

**Timescales at n₀ = 0.025** (all derived):

| Quantity | Value |
|---|---|
| Precursor length, X = 0.3 | 7–74 pc |
| Precursor length, X = 0.05 | 41–440 pc |
| Flow time | 0.05–4 Myr |
| Recombination time | 4.9 Myr |
| Cooling time (Λ ~ 10⁻²⁴) | 6 Myr |

**Consequences.**
- The precursor makes no Hβ and no [O III] worth counting.
- It cannot be the [O III] source in any case: T_e of photoionized gas is ≤ 1.5–2×10⁴ K, which gives 4363/5007 ≤ 0.02–0.034 and R23 ≤ 60–115 (my atom, all O as O²⁺, solar O/H). The measured values are 0.083–0.092 and 76–130.

## 9. New findings and tensions

1. **The magnetic field is set by the closure.** Matching the 0.28° width, 1–4 R and the Monogem pressure together needs b ≈ 7–14 µG cm^{3/2}, i.e. B₀ ≈ 0.8–2.9 µG. Compression is magnetically limited at x ≈ 5.6–12 (median 7.8) in the [O III] zone. Weaker fields give thin, over-dense layers; stronger fields give thick, faint ones.
2. **[O II] doublet tension** (row 9): n_e(doublet) ≈ 100 vs model 0.2. This is the single strongest observational argument against a pressure-balanced shock, and it should be re-measured at higher resolution.
3. **Resonance scattering explains the C IV deficit quantitatively.** It needs no extinction (E(B−V) foreground ≈ 0.006) and no ad hoc C depletion.
4. **Predicted morphology.** Sharp east edge and tail toward Monogem; FUV ridge ~0.1° outside the [O III] ridge.
5. **The north fits a more advanced (or slower) stage** of the same family.
6. **The paper's 10–13 R** would need P/k ≈ 2–4×10⁵ with b ≈ 20–30 at the observed width. Pressure alone does not forbid that, but the measured 1.2–2.3 R (DESI/SDSS) does.

## 10. Limitations

- **No radiative transfer.** Hβ at intermediate N is ~2× below R88-type models, which I treat as a 0.32 dex systematic. The precursor ionization is an input, not computed.
- **Fe cooling** is in CIE, so Λ is 30–40 % low at 3–5×10⁵ K. This affects v ≳ 140 km/s cooling columns at the ~0.1–0.2 dex level.
- **Geometry.** A single spherical sheet of R_c = 58 pc is assumed. Ripples (Hester 1987) give G ∝ √R_c; R_c = 15 pc halves the brightness.
- **C IV Monte Carlo.** Slab geometry, uniform emissivity; scattering out of the sheet is assumed to escape (no re-absorption elsewhere).
- **FUV prediction** omits N IV], Si IV, C II, N III] and 2γ (≲ +35 CU/R). NUV is not modelled.
- **Abundances.** Asplund09 throughout, with one depleted variant. N/O, S/O and Ne/O uncertainties enter via the 0.1 dex term. S⁺/S²⁺ have no DR.
- **Input revisions.** The observed ratios combine DESI and SDSS values that disagree at the 2σ level ([O II]/[O III] 1.1 vs 1.8; [O III]/Hα 2.7 vs 13.7).
- **Grid resolution.** b is gridded 1/3/5/7/10/14/20/30, so the "consistent" ranges are grid-limited.
- **Literature values transcribed by eye.** R88 Fig. 8a values were read by eye from the scan (±0.05 dex). The R79, CR85 and DS96 numbers are the xrayfuv agent's transcriptions; I re-read only R88 Tables 5 and 7.

## 11. Files

**Scripts.**
- `scripts/verify_shockmodel_atomic.py`, `verify_shockmodel_flow.py`: atomic data, the 1-D model.
- `verify_shockmodel_validate.py`, `verify_shockmodel_r88.py`: validation.
- `verify_shockmodel_fit.py`: grid and χ².
- `verify_shockmodel_rt.py`: C IV Monte Carlo.
- `verify_shockmodel_analysis.py`: closures b–f.
- `verify_shockmodel_summary.py`, `verify_shockmodel_north.py`, `verify_shockmodel_precursor.py`: summaries, north scan, precursor budget.
- `verify_shockmodel_figures.py`, `verify_shockmodel_valfig.py`: figures.
- `verify_shockmodel_grid.py`: an earlier draft from the interrupted run. It is superseded by `verify_shockmodel_fit.py`, was never run, and produced no output.

**Data.**
- `data/derived/verify_shockmodel_validation.json`, `_r88.json`: validation results.
- `_fit.json`, `_tracks.npz`, `_fit_log.txt`: grid output.
- `_analysis.json`, `_analysis_log.txt`: closures.
- `_summary.json`, `_north.json`, `_precursor.json`: summaries.

**Raw atomic data and literature** (fetched 2026-11-16 by the earlier, interrupted run of this agent):
- `data/raw/verify_shockmodel_atomic/`. URL manifests exist only for `pyneb/` (raw.githubusercontent PyNeb_devel) and `chianti/` (sohoftp CHIANTI dbase). The download URLs of the root-level files were **not logged**, so their provenance rests on their own headers:
  - Verner `cfit.f` (Voronov fits, version 2, 1997);
  - `rrfit.f` (version 4, 1999);
  - Badnell RR/DR clists ((C) 2023);
  - Gnat & Ferland 2012 "Ion-by-Ion Cooling Efficiencies" data tables, plus arXiv:1111.6980;
  - Cloudy c17 `atmdat_char_tran` doxygen listing.
- `data/raw/verify_shockmodel_lit/`: GS07 VizieR tables (URLs in `MANIFEST.txt`); Allen+08 arXiv:0805.0204 PDF and text.

**Figures.**
- `figures/verify_shockmodel_validation.png`
- `figures/verify_shockmodel_tracks.png`
- `figures/verify_shockmodel_closure.png`
- `figures/verify_shockmodel_fuv.png`
