# verify_xrayfuv: Monogem X-ray edge null test, and GALEX FUV/NUV against UV-line predictions

Wave-2 adversarial check, tag `xrayfuv`, written 2026-11-16. Every number below is my own measurement or derivation unless marked *literature*. I wrote all the code (`scripts/verify_xrayfuv_*.py`). Wave-1 scripts were read only to look for bugs, never rerun. Wave-1 data products were read only as inputs: the trace, the NSNS map, and two wave-1 EQ grids (IRIS 100 µm and GALFA N(HI)).

**Sign convention.** The spine trace is `silk_trace.json` v1.0 (E spine). The offset x is positive to the EAST. I checked this: P4 gives −1.35°, as in wave-1. So x < 0 is the west side, toward the Monogem centres.

---

## 1. Bottom line

### X-ray: the Silk lies on the outer flank of the Monogem shell, but no specific edge at the spine is established

1. **Where the Silk sits.** The southern Silk (Dec < 17.6) lies on the declining outer flank of the Monogem 0.2–0.6 keV shell:
   - It is 11.36 ± 0.15° from PSR B0656+14 and 12.10 ± 0.12° from the Knies+24 two-SNR centre. This reproduces the context note.
   - The ring's brightness peak in PA 60–100° is at r ≈ 10.25°, so the Silk is ~1° outside it.
2. **The W/E brightness ratio is real but generic to the shell.** It is W/E = 1.26 in eRASS1 (1.22 south, 1.25 north) and 1.24 in ROSAT R1+R2.
   - Against random sky translations of the trace it is significant: p = 0.007 (z ≈ 3.2).
   - Against rigid rotations of the trace about the PSR, i.e. the *same* Monogem radius at other position angles, it is entirely typical: p = 0.37 (z = 0.3); south p = 0.57.
   - The wave-1 "2.4σ" used parallel-shift nulls. Those cannot tell "Silk = edge" apart from "Silk lies at the Monogem radius".
3. **A local step at the spine is marginal (~2σ), south only.** I fitted a step beyond a linear gradient at x = 0.
   - **South:** +12 % (+0.78 cts s⁻¹ deg⁻²), best location x = −0.2°.
     - p = 0.013–0.033 against sky translations and 0.027–0.08 against rotations, depending on the statistic.
     - It is still ~2σ against a smooth-shell model (p = 0.014–0.04). With high-latitude-only nulls it is p = 0.011–0.08.
     - Not robust: it drops to p ≈ 0.04–0.25 with a realistic absorption correction, and it is not significant after the ~3–9 statistic × segment trials.
   - **North:** none (+1.5 %, p ≈ 0.3).
   - **ROSAT:** +6 % south, p ≈ 0.15–0.23.
   - **Radial scan:** the step peaks 0.5° *inside* the spine (+19 %), not at it. ROSAT R1+R2 south shows a larger step 4° outside.
4. **The Silk edge is softer than other parts of the shell.** In PA 60–80° the steepest 1° drop (×1.27) is centred at r = 11.25°, i.e. at the southern Silk. Elsewhere on the shell the edges are sharper, for example ×1.43–1.57 per degree at PA 160–240°.
5. **No X-ray filament on the spine.** The on-spine bump relative to the ±0.3–0.6° interpolation is −0.04 (all), −0.38 (south) and +0.07 (north) cts s⁻¹ deg⁻², on a level of 5.4–7.1. ROSAT gives −3 %.
6. **Systematics checked.**
   - Exposure (eRASS1 catalogue MLExp1 map, median 78 s): no step (+0.3 %).
   - N(HI) has no step that could make the X-ray step. In the south the interior has *more* HI (+5×10¹⁹ cm⁻²), so absorption works against the observed sense.
   - The 1.0–2.3 keV control band shows nothing.
   - Results are stable to the point-source masking variant (none, catalogue, catalogue+clip) and to the choice of centre (PSR vs K2).

### FUV/NUV: a significant UV filament, consistent with an incomplete shock with suppressed C IV, and not with complete shocks or interfaces

7. **The GALEX FUV and NUV excess is a clear detection.** It is filament-shaped and co-spatial with the NSNS [O III] in the south (Fig. `verify_xrayfuv_uv.png` a–c).
   - Template regression on NSNS [O III], with Tycho-2 and NSNS star masks:

     | Band | South (Dec < 17.6), CU per NSNS-R | Whole spine, CU per NSNS-R |
     |---|---|---|
     | FUV | **92 ± 10** (z = 8.7) | 74 ± 8 |
     | NUV | **49 ± 7.5** | 57 ± 5.5 |

   - The direct contrast on the E spine in the south is FUV +120 ± 26 CU (z = 4.6) and NUV +96 ± 23 CU (z = 4.2).
   - The stacked peak at |x| < 0.1 is 145 CU (FUV) and 124 CU (NUV) at NSNS 1.55 R.
   - Wave-1 called the FUV signal "a weak excess at most". Its numbers (+47…+154 CU) are close to mine, but the significance is higher here (≥ 4.5σ in the south).
8. **The FUV/[O III] ratio rules out complete shocks and requires C IV suppression ≥ 5× for incomplete ones.**
   - Measured FUV per R(5007): 92 ± 10 CU (NSNS scale), or 68 ± 8 CU if the true [O III] is 1.35× NSNS (the SDSS precovery amplitude).
   - Predictions at E(B−V) = 0.03:
     - Unsuppressed shocks, complete or incomplete: 196–541 CU per R.
     - Non-resonant floor: 57 CU per R for the incomplete R88 SE/M models, which are also the best optical analogues ([O III]/Hβ 10.5–12.9, [O II]/[O III] 1.34–1.64). It is 85–277 CU per R for complete shocks.
   - The measurement sits 1.2–1.6× above the R88 floor. It needs C IV suppressed by 6.5–26× (≥ 5 at 1σ). For comparison, Kim+2010 required ≥ 20× in G65.3+5.7 R2 (*literature*).
9. **The paper's 10–13 R does not survive this test.** If the NSNS scale were 6× low, the FUV per true R would be only 12–15 CU. That is 4–8× *below the non-resonant floor of every shock model*, and non-resonant lines cannot be scattered away. Closing the gap would need E(B−V) ≈ 0.18, far above SFD's 0.03–0.04.
10. **The NUV excess is consistent with incomplete-shock NUV lines.** Measured 36–49 CU per R (joint fit 38 ± 8; the stacked-profile peak ratio gives ~80 CU per R, so that is the systematic range). The R88 prediction is 47–63 CU per R.
    - Composition for R88 SE: C III] 1909 ≈ 35 %, C II] 2326 ≈ 24 %, H two-photon ≈ 22 %, [Ne IV] 2423 ≈ 13 %, Si III] 1892 ≈ 4 %, [O II] 2470 ≈ 3 %.
    - Mg II 2800 contributes < 1 CU per R in R79 because of the low NUV response at 2800 Å.
    - So the NUV excess is *not* mainly [O II] 2470 or Mg II.
    - Dust-scattered light is bounded at < 16 CU (NUV) and < 34 CU (FUV). This uses my measured local dust response (160 CU and 69 CU per MJy sr⁻¹ of high-passed IRIS 100 µm) and the wave-1 2σ limit ΔI₁₀₀ < 0.19 MJy sr⁻¹. Both bounds are well below the measured excess.
    - Colour: FUV/NUV of the excess is 1.0–1.9, depending on the estimator. That is 2–4× bluer than the diffuse background (0.50), but similar to small-scale dust-scattered structure (2.3), so colour alone is not diagnostic.
11. **Conductive interfaces fail on brightness.** I bounded [O III] 5007 from the BBF90 O III] 1663 fluxes using the smallest shock-table N(1663)/N(5007) = 0.135 photons.
    - That gives face-on I(5007) ≤ 3×10⁻⁴ to 1.6×10⁻³ R per front.
    - Reaching 1 R would need ≥ 640–2900× limb-brightening or that many fronts along the line of sight.
    - At fixed 5007 they also give FUV ≥ 119–715 CU per R.
    - Slavin 1989's C IV is only 0.0035 R (*literature*).

---

## 2. Claim → verdict table

| # | Wave-1 claim | Verdict | My value |
|---|---|---|---|
| 7a | Spine on the outer (E/N) edge of the Monogem soft X-rays; eRASS1 0.2–0.6 keV W/E 1.32, 2.4σ vs parallel-shift nulls | **REVISED** | W/E 1.26 (S 1.22, N 1.25). Significant vs random sky placement (p = 0.007), typical vs same-radius rotations (p = 0.37). The Silk is on the Monogem shell's outer flank, not a distinct edge. |
| 7b | ROSAT R12 W/E 1.31, 1.4σ | REVISED | R1+R2 W/E 1.24 (S 1.27, N 1.21). Sky p = 0.026, rotation p = 0.51. |
| 7c | X-ray half-drop edge 0.3–1.0° inside the spine | CONFIRMED (roughly) | Best step location −0.2° (S), −0.5° (all), −1.0° (N) in eRASS1; −0.7…−1.3° in ROSAT. The radial scan peaks 0.5° inside. |
| 7d | (implicit) a genuine edge at the Silk beyond the smooth shell decline | **not established, ~2σ, south only** | +12 % local step (p 0.01–0.08 raw, 0.04–0.25 absorption-corrected); against a smooth-shell model p 0.014–0.04; none in the north. |
| 7e | No X-ray filament on the spine | CONFIRMED | On-spine bump −0.04 (all) to −0.38 (S) cts s⁻¹ deg⁻² on a level of 5.4–7.1. |
| 7f | Absorption cannot explain the interior excess (ΔN(HI) = −1.3×10¹⁹) | CONFIRMED, with a twist | HI4PI ΔN(W−E): all −0.7×10¹⁹, **south +5.0×10¹⁹** (interior has more HI, wrong sign), north −2.4×10¹⁹. Implied W/E change at f = 0.15–0.3: S −4…−14 %, N +2…+7 %. |
| 8a | South spine at the Monogem radius (12.10 ± 0.12° K2; 11.36° PSR) | CONFIRMED | 12.10 ± 0.12° and 11.36 ± 0.15°. |
| 7g | GALEX FUV on-filament +47…+154 CU (2σ UL 110–250 CU ≈ 0.4–0.9 R at 1549 Å) | **REVISED (detection, not just a limit)** | Contrast +67 ± 23 (all), **+120 ± 26 (south, z = 4.6)**. Template β = 74 ± 8 (all) and 92 ± 10 CU per NSNS-R (south). The conversion is 1 R at 1549 Å = 223 CU (W_eff = 357 Å), or 178 CU after E(B−V) = 0.03, so 110–250 CU ≈ 0.6–1.4 R. |
| 7h | Predicted 900–4400 CU incl. C IV and 210–300 CU floor (O III]+He II) for 10–13 R | REVISED | O III]+He II alone: 22.6 CU per R (226–294 at 10–13 R; reproduced). Adding N IV] 1486, N III] 1750 and two-photon gives a floor of 57 CU per R (R88), i.e. 570–740 at 10–13 R. Totals 2900–4700 (R88), 1960–5400 (R79). |
| 7i | NUV narrow excess +64…+113 CU, nature unclear; FUV/NUV colour ≈ background (0.71) | **REVISED** | NUV +68 ± 12 (all), +96 ± 23 (S), β 49–57 CU per NSNS-R. This matches incomplete-shock NUV lines (C III], C II], 2γ, [Ne IV]): 47–63 CU per R. Excess colour 1.0–1.9 vs background 0.50 (not 0.71). Dust < 16 CU. |
| 1 | "Weak excess at most, well below a fast radiative shock" | REVISED | Strong UV filament. Below *unsuppressed* shocks by 3–6×, but consistent with an incomplete 120–140 km/s shock if C IV is suppressed ≥ 5×. |

---

## 3. Methods

### 3.1 X-ray (`verify_xrayfuv_erass_prep.py`, `_xray_edge.py`, `_xray_absorb.py`, `_xray_highb.py`)

**Data.**
- Raw eRASS1 DR1 half-sky AIT rate maps (cts s⁻¹ per 3′ pixel, TM1–4+6), cut to l 180.05–226, b −12…+38 and converted to cts s⁻¹ deg⁻² (×400). The bands are summed to 0.2–0.6 keV; 0.6–1.0 and 1.0–2.3 keV serve as controls.
- The product is a smoothed image: adjacent-pixel correlation is 0.73 at lag 1, 0.33 at lag 2 and 0.16 at lag 3.

**Point sources.**
- The eRASS1 Main catalogue (Merloni+24, J/A+A/682/A34/erass1-m; 84,548 rows) was fetched in RA strips from the CfA VizieR mirror. Five strips first failed with a "database not reachable" error and were refetched. The failures are logged as `*.failed_*.txt`.
- Median-stacking the map at catalogue positions gave the map's response to a source. Mask radii follow from it: 2/3/4/6/15 px for MLRate1 0.1–0.3/0.3–1/1–5/5–50/>50 cts s⁻¹, plus 3×EXT for extended sources. This masks 21 % of pixels.
- I then clipped pixels above 5 robust σ relative to a 1° running median (an extra 1.4 %). Three variants were compared: none, catalogue, catalogue+clip.

**Ancillary maps.**
- Exposure: the median MLExp1 of catalogue sources in 1° cells.
- N(HI): HI4PI from a wide SkyView CAR map (l 203, b 13; 50° × 50°, 0.2°). SFD E(B−V) was fetched the same way.
- ROSAT: RASS diffuse R1+R2 and R4+R5, wide SkyView maps on the same grid.

**Statistics.** On the cross-spine profile (0.1° bins for eRASS, 0.2° for ROSAT, x ∈ ±2.5°):
- W/E = mean(−2.5…−0.5)/mean(+0.5…+2.5);
- c = step at x = 0 in I = a + bx + c[x < 0]; cfrac = c/a;
- the best step location, scanned over −2…+2°.
- Segments: all, south (Dec < 17.6), north.

**Nulls.**
- (R) Rigid rotations of the spine about PSR B0656+14 (or K2) in 4° steps; n = 73 (68 for K2).
- (T) 150 random rigid translations within ±22° RA × ±18° Dec.
- (S) A radial scan with the trace moved by Δr = −6…+6°.
- Placements were rejected if they lie outside eROSITA-DE (l < 180.3), have < 50 % valid pixels, or overlap the real spine by > 10 %.
- I quote z = (obs − mean)/σ_null and the one-sided empirical p = (k+1)/(n+1).

**Smooth-shell model.** I(r, PA) = B(PA) + A(PA)·g(r), with g(r) the median normalised radial profile over all 15° sectors except the Silk's, fitted per sector. The residual map was then put through the same tests.

**Absorption.**
- (i) The N(HI) profile across the spine.
- (ii) Corrections exp(σ f N(HI)) with σf = 1.5×10⁻²¹ and 3×10⁻²¹ cm². This assumes σ_eff ≈ 0.5–1×10⁻²⁰ cm², which is an approximate band-averaged literature-order value, not computed from a response matrix. The foreground fraction f is 0.15–0.3 from the context note's 3D dust (E(B−V) < 300 pc ≈ 0.005 against SFD 0.03–0.04).
- The f = 1 runs (σ = 5×10⁻²¹ and 10⁻²⁰) are kept in the JSON only as an unphysical extreme: they amplify the map by e²–e⁵.
- (iii) High-latitude-only nulls (min b > 8°).

### 3.2 UV predictions (`verify_xrayfuv_uv_predict.py`)

**Model line tables.** I transcribed them **by hand from the ADS scans** in `data/raw/verify_xrayfuv_lit/`. All are energy ratios to Hβ.
- Raymond 1979 (ApJS 39, 1), Tables 1/3/4, models E–I: v_s = 100/120/141/160/200 km/s, n₀ = 10, complete.
- Cox & Raymond 1985 (ApJ 298, 651), Table 3, runs 1–5: 200/186/140/130/119 km/s, complete, with preionization families. This includes the two-photon continuum, 18–33 × recombination Hβ.
- Dopita & Sutherland 1996 (ApJS 102, 161), Tables 8/9, shock only: s15010–s20014, v = 150 and 200 km/s, B/√n = 0–4.
- Raymond et al. 1988 (ApJ 324, 869), Tables 5 and 7, the **incomplete/truncated** Cygnus Loop spur models:
  - SE: 140 km/s, n₀ = 4.4, log N_H(trunc) = 18.20;
  - M: 120 km/s, n₀ = 6.0, log N_H = 18.25.
- Borkowski, Balbus & Fristrom 1990 (ApJ 355, 501), Table 1: conductive fronts (photon fluxes).
- Slavin 1989 (ApJ 346, 718), Table 2: the C IV of the Local Cloud boundary.
- Hartigan, Raymond & Hartmann 1987 was downloaded but not transcribed (its planar grid is for n₀ = 100–1000 HH conditions).
- Allen et al. 2008 / MAPPINGS III library: not obtained. The ITERA host is parked and the GitHub API returned nothing.

**Conversion to CU.**
- GALEX effective areas come from the SVO FPS (GALEX team curves; Morrissey+2007).
- A line of N photons s⁻¹ cm⁻² sr⁻¹ gives N·A(λ)·λ_eff/∫Aλ dλ CU. The resulting W_eff values (Å) are:

  | Band | Line | W_eff (Å) |
  |---|---|---|
  | FUV | C IV 1549 | 357 |
  | FUV | He II 1640 | 516 |
  | FUV | O III] 1664 | 636 |
  | FUV | N IV] 1486 | 261 |
  | FUV | 1400 blend | 779 |
  | NUV | C III] 1909 | 1645 |
  | NUV | C II] 2326 | 913 |
  | NUV | [O II] 2470 | 865 |
  | NUV | Mg II 2800 | 4238 |

- The curve-based cps → F_λ factors match the published 1.40e-15 / 2.06e-16 to 2 % (FUV) and 8 % (NUV).
- Photons are converted from energy ratios with λ_line/5007.

**Extinction.** CCM89 with R_V = 3.1, implemented by me. Check: A_FUV/E(B−V) = 8.15 at 1528 Å; A(1549)/E(B−V) = 8.1. E(B−V) = 0.03 is the default and 0.01 is also tabulated.

**He II and two-photon continuum.**
- He II 1640 = (5.6–12.5) × 4686, the CR85 case A/B ratios, where 1640 is not tabulated (R79, DS96).
- H two-photon spectrum: the Nussbaumer & Schmutz (1984) fit A(y), which I checked integrates to 2 photons and hν_Lyα per decay. The total comes from the tables (R88, CR85) or is taken as (14–21) × Hβ (R79, DS96).

**C IV handling.** C IV is treated as a separate resonant term with suppression ×1/3.5/20. The "non-resonant" floor excludes C IV and the Si IV/O IV] 1400 blend.

**Predicted CU per 1 R of [O III] 5007** (E(B−V) = 0.03; multiply by 1, 2, 4 for I(5007) = 1, 2, 4 R):

| model | v_s | [O III]/Hβ | FUV (C IV ×1) | ÷3.5 | ÷20 | FUV non-res. | NUV |
|---|---|---|---|---|---|---|---|
| R79 E–I (complete) | 100–200 | 2.8–4.5 | 196–414 | 123–187 | 99–148 | 85–125 | 67–93 |
| CR85 runs 1/3/5 (complete) | 200/140/119 | 1.9–3.2 | 311–363 | 187–270 | 146–240 | 122–206 | 116–159 |
| DS96 150, 200 (B 0–4) | 150–200 | 1.2–1.9 | 446–541 | 298–384 | 249–332 | 196–277 | 115–187 |
| **R88 SE (incomplete)** | 140 | 12.9 | **360** | 152 | 83 | **57** | **47** |
| **R88 M (incomplete)** | 120 | 10.5 | **293** | 130 | 77 | **57** | **63** |
| BBF90 interface, 2.2/7.4 Myr | — | — | ≥ 171 / ≥ 119 | — | — | ≥ 59 / ≥ 44 | ≥ 78 / ≥ 51 |

So for I(5007) = 1/2/4 R with the R88 models: FUV is 293–360 / 586–720 / 1170–1440 CU unsuppressed, with a floor of 57 / 114 / 228 CU; NUV is 47–63 / 94–126 / 188–252 CU.

With E(B−V) = 0.01 the totals rise by ~15 %.

### 3.3 GALEX measurement (`verify_xrayfuv_galex_blocks.py`, `_galex_measure.py`, `_galex_dust.py`)

**Blocks.**
- All 267 GR6/7 skybg maps (127 FUV, 140 NUV) were block-averaged to 1′ (40 × 40 px). Only fully valid blocks within 0.55° of the tile centre were kept: 916k blocks.
- Conversion: 2.036×10⁶ (FUV) and 4.453×10⁵ (NUV) CU per cps/pixel. This matches wave-1.

**Cleaning.**
- Per-tile median offset.
- Mean fractional radial (flat-field/vignetting) pattern removed.
- Tycho-2 VT < 11 star masks, r = max(2′, 12′ × 10^(−0.2(VT−5))), capped at 15′.
- NSNS bright-star residual mask (variance > 98th percentile).
- 4-MAD clipping within each tile. About 10 % of blocks are star-masked and 18 % masked in total.

**Estimators.**
- (a) Regression on the NSNS DR0.2 [O III] map. The template is high-passed with a 1° median filter and smoothed to 3′. Per-tile demeaning removes airglow and zodiacal offsets.
- (a2) A joint fit with the NSNS Hα and the star-subtracted continuum ("vc") as nuisance templates.
- (b) Direct contrast within tiles: |x| ≤ 0.15° minus 0.30–0.55°.

**Nulls.**
- 48 RA shifts (±0.7…±3°) plus 16 RA+Dec shifts.
- For (a), none of the 64 nulls reached the observed β. That gives p ≤ 0.015, the minimum possible with 64 nulls.

**Dust check.** Regression of the same blocks on high-passed IRIS 100 µm (wave-1 EQ grid) gives:
- FUV: 160 ± 8 CU per MJy sr⁻¹;
- NUV: 69 ± 6.5 CU per MJy sr⁻¹.

---

## 4. New findings

1. **A significant GALEX FUV + NUV filament on the southern Silk.** It is co-spatial with the [O III] in its cross-section, with matching width and peak position (Fig. c). This is the first quantitative UV counterpart of the Silk. FUV/[O III] ≈ 70–90 CU per R and NUV/[O III] ≈ 36–57 CU per R.
2. **The FUV/[O III] ratio is a calibration-independent test.** It rules out complete shocks (whose non-resonant floor alone is 85–277 CU per R, against the measured 68–92). It also rules out the paper's 10–13 R scale under any shock model (12–15 CU per R would sit below every non-resonant floor). The ratio favours an *incomplete* 120–140 km/s shock (R88-like), with C IV suppressed ≥ 5× by resonance scattering in an edge-on sheet or by a different geometry. Note that R88 SE independently matches the Silk's optical diagnostics:
   - [O III]/Hβ 12.9 vs > 14–20 observed;
   - [O II]/[O III] 1.34 vs 1.1–1.5;
   - 4363/5007 0.071 vs 0.056 (0.022–0.077).
3. **The NSNS Hα map also correlates with GALEX** (z ≈ 7–11, including the south where there is no Hα filament). Likely causes: diffuse Galactic light or WIM structure in the NSNS continuum-subtracted maps, which also appears in dust-scattered UV. The joint fit keeps the [O III] term significant in the south, but this is a caveat for template-based UV/optical ratios.
4. **The Monogem S/SW limb has the sharpest X-ray edge** (×1.4–1.6 per degree at PA 160–240°). The Silk sector is softer (×1.27), so "the Silk marks the outer edge" is not supported beyond the generic shell geometry.

## 5. Limitations

- **eRASS1 depth.** Exposure is ~78 s here, so Poisson noise is negligible next to sky-structure systematics. The nulls capture those, but rotated placements sample other parts of the shell and of the Galactic plane. The high-b-only rotation null has only n = 11.
- **Absorption.** It is uncertain (unknown foreground fraction f; approximate σ_eff). Absorption-corrected nulls are broader, and with them the W/E significance against the sky drops to ~0. I read this as the correction being uncertain, not as a demonstration that absorption makes the edge.
- **Trials.** About 3 segments × 3 statistics × 2 surveys were examined; the 2σ south step should be read with that in mind.
- **NSNS calibration.** The per-R UV ratios assume the NSNS Rayleigh scale (DESI and SDSS suggest the true values are ×1.0–1.5). The FUV conclusion holds across that range.
- **NUV estimator spread.** The NUV per-R value depends on the estimator: 38–57 from regressions, ~80 from the stacked-profile peak ratio. The FUV/NUV excess colour ranges 1.0–1.9.
- **Coverage.** GALEX skybg coverage stops at Dec 22, so the northern Silk (Dec 22–28) is untested.
- **Model tables.**
  - They were hand-transcribed from scans. I double-read the key rows, but transcription errors of a few per cent are possible.
  - R88's "1400" is the O IV]+Si IV blend.
  - He II 1640 and the two-photon totals for R79/DS96 are scaled assumptions (ranges in the JSON).
  - No magnetic or precursor-photoionization spectra (DS96 precursor tables) were included.
  - The BBF90 [O III] value is a bound derived from O III], not a model output.
  - Allen+08 and Sutherland & Dopita 2017 were not obtained.
- **Interface geometry.** Face-on interface intensities assume F/4π. A one-sided emission convention would double them, which does not change the ≥ 300× shortfall.

## 6. Files

**Scripts.** All `scripts/verify_xrayfuv_*.py`:
- `common`
- `fetch_skyview`, `fetch_erass1cat`
- `erass_prep`, `xray_edge`, `xray_absorb`, `xray_highb`, `fig_xray`
- `uv_predict`
- `galex_blocks`, `galex_measure`, `galex_dust`, `fig_uv`

**Derived products** (`data/derived/`):
- `verify_xrayfuv_erass_pixels.npz`
- `verify_xrayfuv_xray_edge.json`, `verify_xrayfuv_xray_nulls.npz`
- `verify_xrayfuv_shellresid_{PSR,K2}.npy`
- `verify_xrayfuv_xray_absorb.json`, `verify_xrayfuv_xray_highb.json`
- `verify_xrayfuv_uv_predictions.json`
- `verify_xrayfuv_galex_blocks.npz`, `verify_xrayfuv_galex_measure.json`, `verify_xrayfuv_galex_dust.json`

**Raw downloads.** Each folder has a `MANIFEST.jsonl` with URL, date (2026-11-16), bytes and sha256.
- `data/raw/verify_xrayfuv_lit/`: 7 ADS scans and the SVO GALEX curves.
- `data/raw/verify_xrayfuv_skyview/`: wide HI4PI, SFD, RASS R1/R2/R4–R7.
- `data/raw/verify_xrayfuv_erass1cat/`: eRASS1 Main catalogue strips.

**Figures.**
- `figures/verify_xrayfuv_xray_edge.png`: the eRASS1 map with null placements; cross-spine profiles; null distributions of the step and of W/E; radial scan; Monogem sector profiles.
- `figures/verify_xrayfuv_uv.png`: FUV and NUV residual maps with [O III] contours; stacked cross-spine profiles; predicted vs measured FUV and NUV per R; the FUV–NUV plane.
