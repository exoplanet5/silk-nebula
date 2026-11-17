# verify_distance: interstellar absorption and the distance of the Silk

Tag `distance` (wave 2, adversarial check). Date 2026-11-16. This work used my own code only. I read the wave-1 notes and scripts to understand their method and to look for bugs, but I reran none of them.

**Scripts**
- `scripts/verify_distance_uves.py`: shared reader.
- `scripts/verify_distance_hd61717.py`: part (a).
- `scripts/verify_distance_sophie_check.py`: part (a).
- `scripts/verify_distance_iuelib.py`: IUE reader.
- `scripts/verify_distance_iue.py`: part (b).
- `scripts/verify_distance_synthesis.py`: part (c).
- `scripts/verify_distance_figures.py`: figures.

**Derived outputs:** `data/derived/verify_distance_*`.

**Figures**
- `figures/verify_distance_constraints.png`: the main figure with all constraints.
- `figures/verify_distance_hd61717.png`
- `figures/verify_distance_iue.png`

**New downloads.** Each is logged with its URL in a `fetch_log.tsv` in the same folder.
- `data/raw/verify_distance_iue/`: 34 IUE NEWSIPS SWP high-dispersion `.mxhi.gz` files, 22 MB, from `https://archive.stsci.edu/missions/iue/data/swp/...`.
- `data/raw/verify_distance_meta/`, small tables:
  - Gaia DR3 row for HD 61717 (VizieR CfA);
  - APOGEE DR17 star and visit rows (SDSS SkyServer DR17 SQL);
  - SIMBAD rows for the UV stars (CfA TAP).

Velocities are heliocentric unless labelled LSR. "Measured" means computed in this work from spectra. "Literature" means a catalogue value.

## 1. Results table

| # | Wave-1 claim | Verdict | My value / reason |
|---|---|---|---|
| a1 | HD 61717 has a narrow Ca II component at V_helio −17 km/s | **CONFIRMED** | v = −16.88 ± 0.04 (stat) ± ~0.5 (wavelength scale) km/s; log N(Ca II) = 11.59 ± 0.04; b = 1.49 ± 0.30 km/s. K and H fitted jointly with optical depth, LSF and continuum varied in a Monte Carlo. The same line appears in K and H with W_K/W_H = 1.70, so it is real absorption. V_LSR = −27.7 km/s. |
| a2 | N(Ca II)/N(Na I) = 22 ± 2, a lower limit because of saturation | **REVISED** | **23, with 68 % range 15–49.** Ca II is corrected for saturation (+0.04–0.07 dex over the AOD value), so the ratio is now a value, not a lower limit. Na I is weak: W(D2) = 3.0 ± 0.9 mÅ, W(D1) = 1.9 ± 1.4 mÅ, log N = 10.23 ± 0.38 (joint fit with a HARPS telluric template). D1 sits on a weak telluric line (0.8 mÅ removed). The ±2 error was underestimated 5–10×. The local gas on the same sightline gives 2.5. |
| a3 | The −17 component has a Routly–Spitzer (shocked-gas) signature | **REVISED (not supported)** | See below the table. |
| a4 | HD 61717 is an A2 star with RV +41.8 (used for stellar-blend vetting) | **REVISED** | HD 61717 is a **double-lined binary (SB2)**. See below the table. |
| a5 | The component is absent toward HD 64705 (280 pc) and HD 61038 (149 pc) | **CONFIRMED and extended** | A narrow-line matched filter (15 km/s high-pass) at −16.88 km/s finds nothing toward 14 other stars: log N(Ca II) < 9.9–10.9 (3σ), against 11.59 toward HD 61717. The full list is below the table. The cloud is compact, smaller than ~1.9° (< 9 pc at 272 pc). |
| a6 | No LISM cloud is at this velocity (Redfield & Linsky calculator) | **CONFIRMED** | All 15 model clouds project to positive velocities: +4.8 (Hyades) to +52 km/s. The sightline traverses the LIC (+19.3) and Gem (+33.3) and passes near Aur (+19.5). The wave-1 query was at (l, b) = (198.42, 20.81), 0.6° from HD 61717's true (197.87, 20.53); that shifts projections by < 1 km/s. The −17 gas is 22–50 km/s from every LISM cloud, so it lies beyond the CLIC. |
| a7 | Multi-epoch stability of the −17 line | **UNVERIFIABLE** | Only one high-resolution epoch exists in ESO, ELODIE and SOPHIE (UVES, 2023-11-12, prog. 112.25VN.001). APOGEE has 22 epochs, but the H band has no Ca II. |
| b1 | HD 60848 has 122 unmeasured IUE spectra | **REVISED** | The 122 are MAST product rows. There are **14 SWP observations (13 high-dispersion)** plus LWR, all now measured: 12 large-aperture and 1 small-aperture (swp02897) spectra. |
| b2 | C IV / Si IV / N V toward HD 60848 (no claim in wave 1) | **new measurement** | Measured log N, 13 co-added spectra: S II 15.01 ± 0.05, Al III 12.86 ± 0.02, Si IV 12.80–12.86 (± 0.06), C IV 13.32 ± 0.07 (1550 AOD) to 13.49 ± 0.05 (joint fit), **N V < 12.88 (3σ)**. Kinematics and limits are below the table. |
| c1 | The south spine is at 12.10 ± 0.12° from the Knies two-SNR centre and 11.36° from B0656+14 | **CONFIRMED** | 12.10 ± 0.12° (23 points with Dec < 17.6) and 11.36 ± 0.15°. The north spine is at 12.4–16.0°. |
| c2 | If on the Monogem limb, the Silk is at 248–320 pc | **CONFIRMED** | Monte Carlo over the PSR parallax (3.47 ± 0.36 mas) and a 15 pc centre-depth offset: **limb 282 pc, 68 % 252–318, 95 % 225–359**. Near/far side of a shell up to 1° larger: 237–337 pc. Shell radius 61 (54–68) pc. |
| c3 | The Local Bubble wall toward P1–P4 is at 161–171 pc | **CONFIRMED** | O'Neill+2024 file, nearest pixel at every spine point: peak median 167 pc (range 158–223), south 165, north 178; inner and outer edges 139 and 200 pc. |
| c4 | PSR J0751+1807: parallax 1.6 ± 0.8 mas, no RM | **REVISED** | The ATNF v2.9.0 database already on disk gives **PX = 0.85 ± 0.04 mas (ref. abb+23), so d = 1.18 (1.12–1.24) kpc, and RM = +41.6 ± 0.2 rad m⁻² (sbm+22)**. DM 30.24, <n_e> = 0.026 cm⁻³, <B∥> = 1.7 μG. Wave 1 used the 2017 VizieR B/psr snapshot. The pulsar is 0.95° E of the spine, on the convex side, so it is not behind the filament. |
| c5 | Kinematic distances are meaningless here | **CONFIRMED** | Flat rotation gives dV_LSR/dd ≈ 7.4 km/s per kpc at 1 kpc. With σ_pec = 7 km/s: DESI/SDSS [O III] (V_LSR ≈ +1…+3) allows 0–1.4 kpc; LAMOST (+14) allows ≥ 1 kpc. Neither is physical. |
| c6 | Extinction gives no constraint | **CONFIRMED (by reasoning)** | The filament is dust-free, and Hβ is undetected, so there is no Balmer decrement. Any distance below ~600 pc has E(B−V) < 0.01. I did not re-read the Edenhofer cube; the values come from `notes/context.md`. |
| c7 | Ca II / Na I absorption cannot bracket the [O III] gas | **CONFIRMED and strengthened** | The revised [O III] V_helio (+12.5…+14.5) coincides with the local gas (+6…+15), which is present toward every star from 63 pc outward. The ionized, Balmer-dark gas also need not carry Ca II or Na I. |
| c8 | COS WD J0747+2106: no Si IV / N V within 100 pc | **not re-checked** | This rests on the wave-1 measurement. I use it only as a conditional lower limit. |

**a3: why the shock signature is not supported.**
- **High Ca II/Na I is not specific to shocks.** It is also the normal state of warm, low-density gas inside the Local Bubble.
- **Ti II is low:** **Ti II/Ca II < 0.16 (3σ)** at −17 km/s, against ≈0.48 (±~0.1) in the local gas on the same sightline. The measurements are log N(Ti II, −17) < 10.79 and log N(Ti II, local) = 11.29.
  - My expectation, from the literature and not re-verified here, is that shock-processed high-velocity gas usually shows Ti II enhanced, not depleted.
- **The line is narrow:** b = 1.49 km/s means T < 5300 K if purely thermal, or < 1.05×10⁴ K at +2σ, with little turbulence.
- **Reading:** a cool, quiescent, Ti-depleted, Ca-rich cloudlet at V_LSR −28 km/s. Its link to a shock is unproven.

**a4: evidence that HD 61717 is a binary.**
- **Two components.** Fits to Mg II 4481, Fe II 5018/5316, Si II 6347 and Ti II 4468 at the UVES epoch give components at **−69 and +124 km/s**, v sin i ≈ 113/115 km/s each. Against one rotating star (v sin i 227), ΔBIC = 3053.
- **APOGEE DR17:** 22 visits from 2011 to 2019 give V_helio 35–60 km/s (mean 43.6, σ 6.3). The visits are flagged RV_SUSPECT + SUSPECT_BROAD_LINES.
- **Gaia DR3:** RV 35.7 ± 1.3, Vbroad 160, **RUWE 2.14**.
- **Systemic velocity:** γ ≈ +28…+44 km/s. The −17 line is 45–61 km/s from γ and is at neither star's velocity. A circumstellar or shell origin is therefore unlikely: the lines below show no narrow shell absorption at −17 km/s:

  | Line | EW at −17 km/s (mÅ) |
  |---|---|
  | Ti II 3759 | −1.0 ± 0.4 |
  | Ti II 3761 | −0.6 ± 0.4 |
  | Fe II 4233 | −0.3 ± 0.3 |
  | Fe II 5169 | −0.9 ± 0.4 |

  Multi-epoch data are still needed to exclude it.
- **Distance:** with the RUWE-inflated parallax error, d = 272 (262–284) pc.

**a5: stars with no narrow −17 km/s line.**
- SOPHIE: HD 64705, HD 61038, HD 55869 (423 pc, 4.8° W), HD 53588, HD 66552, HD 65900, HD 67959, 8 Cnc, 64 Gem.
- UVES: HD 52124, HD 53588, HD 63432, HD 63774, HD 68332.
- The plain polynomial-continuum filter does find "signals" in HD 52124 and HD 53588. These are broad (~20 km/s) stellar or broad shoulders and vanish under the narrow-line filter.

**b2: HD 60848 kinematics and limits.**
- **Velocities (NEWSIPS zero point):** low ions +1.5 to +3.4 km/s. **The high ions are 5–11 km/s blueward of the low ions**: Si IV −3.4 ± 2.1, C IV −7.0 ± 2.7.
- **Si IV at +15…+40 km/s:** 1393 Å mean depth 0.007 ± 0.022. Nothing.
- **High velocity (|v| = 70–200 km/s):** none in C IV or N V. C IV EW < ~33 mÅ, so log N < ~12.9 (3σ).
- **C IV 1548 wings** are inconsistent with the 1550 line in both wings: at −60…−25 the depths are 0.21 against 0.04. They are therefore stellar-wind contamination.
- **C IV/Si IV ≈ 3–4.**

## 2. Part (a): HD 61717 in detail

### Frame and resolution
- **Frame.** UVES phase-3 products are labelled SPECSYS = TOPOCENT.
  - I cross-correlated the telluric O₂ γ-band and H₂O lines of HD 61717 against 11 CMi (BERV −29.24, from astropy).
  - The result is −0.5 km/s. If both products were barycentric, the shift would be +56.2. Both are therefore topocentric; I added BERV = +26.92.
  - The local-gas components (+7.4, +14.1 km/s in Ca II; +7.0, +13.6 in Na I) agree with the barycentric HARPS and SOPHIE stars, which confirms the frame.
- **Resolution.** 36 telluric lines in the red arm have FWHM 3.99 km/s (16–84 %: 3.74–4.17), consistent with the header R = 74,450. The adopted blue-arm LSF is 4.36 km/s FWHM, against 4.61 km/s for the header R = 65,030.
- **b depends on the LSF:** b = 1.01 / 1.28 / 1.55 / 1.78 / 2.09 km/s for R = 55k / 65k / 70k / 75k / 85k. log N changes by only 0.17 dex across that range.
- **Noise.** UVES phase-3 ERR overestimates the empirical scatter by 2.5–3× near Ca II (0.0046 against 0.0135 per pixel). I used the empirical noise.

### Ca II components
Fit χ²/dof = 706/319; block-bootstrap Monte Carlo, n = 300.

| v (km/s) | log N | b (km/s) |
|---|---|---|
| −16.88 ± 0.04 | 11.59 ± 0.04 | 1.49 ± 0.30 |
| +7.4 ± 2.2 | 11.22 ± 0.15 | 9.4 ± 2.0 |
| +14.08 ± 0.23 | 11.38 ± 0.08 | 4.1 ± 0.5 |

- AOD over −24…−10 km/s: log N_a(K) = 11.518, log N_a(H) = 11.555. The higher H value shows mild saturation, which the fit corrects.

### Other species at −17 km/s
3σ upper limits on log N from a matched filter with the Ca II line shape:

| Species | 3σ upper limit on log N |
|---|---|
| Ti II 3384 | 10.79 |
| Ca I 4227 | 9.55 |
| CH⁺ 4232 | 11.97 |
| CH 4300 | 11.93 |
| CN 3874 | 11.20 |

- Fe I 3860 has a +1.6 ± 0.25 mÅ residual, but stellar-blend structure is present there, so I do not claim it.
- Ti II is clearly detected in the local gas at +12 km/s: 7.1 ± 0.8 mÅ, log N = 11.29.

### What the −17 cloud implies
- **If it is interstellar** (most likely), there is a compact (< 2°), cool, Ca-rich and Ti-poor cloudlet at V_LSR −28 km/s in front of 272 pc. It lies 0.97° inside the concave side of the arc at Dec +22.3.
  - Relative to the Silk curvature centre, HD 61717 sits 0.93° inside the R = 21.3° circle. The approaching cap of an expanding shell would then show only 0.29 of its expansion speed, so V_exp would be about 100 km/s.
  - HD 61717 is 12.6° from the Monogem two-SNR centre (X-ray radius 11.0°), outside that disc. It is 11.69° from the single-SNR centre, just inside that model's 11.97° radius, so Monogem membership is ambiguous.
- **Its velocity does not match the [O III] gas:** V_LSR −28 against +1…+3 (DESI/SDSS) or +11…+16 (LAMOST). It is not the emitting gas.
- **Distance consequence:** "Silk < 272 pc" holds only if the cloudlet is a cooled fragment of the same shell. One sightline cannot show that. I treat it as a conditional upper limit, not a measurement.

## 3. Part (b): IUE high-dispersion spectra

**Stars measured.** Distances are from SIMBAD/Gaia DR3 parallaxes; separations are from the spine.

| Star | Type | Distance | Separation | Spectra |
|---|---|---|---|---|
| HD 60848 | O8:V:pe | 1.9 kpc | 2.26° W | 13 |
| NGC 2392 | O6f CSPN | 1.8 kpc | 4.0° W | 7 |
| NGC 2371 | CSPN | 1.7 kpc | 2.9° W | 1 |
| HD 58050 | B2Ve | 452 pc | 5.3° W | 4 |
| HD 52559 | B2IV-V | 811 pc | 12.7° W | 2 |
| HD 51354 | B3ne | 559 pc | 11.6° W | 2 |
| HD 53929 | B9III | 225 pc | 12.0° W | 2 |
| HD 66665 | B1V | 1.8 kpc | 7.4° E, outside the arc (control) | 3 |

Two header names map to these stars: HR 2633 = HD 52559 and HR 2817 = HD 58050.

**Reader and normalisation (my own).**
- NEWSIPS MXHI: STARTPIX is 1-based (I checked it against NET ≠ 0). The wavelengths are vacuum and heliocentric (HISTORY "NET CORRECTION VECTOR TO HELIOCENTRIC VELOCITY").
- The heliocentric correction is applied: for HD 60848 the per-spectrum low-ion centroids scatter by 2.1 km/s although LRADVELO spans −30 to +32 km/s. Mg II 1239.9 appears at the expected place in the N V window.
- Each line is normalised with a quadratic continuum, the orders are co-added, and the spectra are aligned on the star's low-ion centroid and co-added.
- I measured AOD columns and joint doublet fits (Gaussian τ convolved with a 25 km/s LSF), with block-bootstrap errors.

**Which sightlines are usable.**
- **Usable for interstellar high ions: HD 60848 only.** Its values are in table row b2.
- **NGC 2392:** C IV, N V and Si IV are saturated and extend 100–200 km/s blueward (W ≈ 240–300 mÅ). This is the CSPN wind and nebular absorption, not the ISM.
- **B2–B9 stars:** photospheric Si IV, C IV and Al III dominate. HD 58050's doublet fits go to the b = 60 km/s bound, so they are unusable.
- **HD 66665 (control):** too noisy at C IV. Its "N V 1238" dip has no 1242 counterpart, so it is not N V.
- **NGC 2371:** one noisy spectrum; unusable.

**Interpretation.**
- **HD 60848's high-ion columns are plausible for a ~1.9 kpc disc path.** I did not check this against a specific survey. No component is separable at the Silk's velocity.
  - The revised [O III] velocity is within ~2 km/s of the low ions, so IUE's ~25 km/s resolution cannot separate the two anyway.
  - The old LAMOST velocity (low ions +13 km/s) is not seen in Si IV: the high ions are blueward, not redward.
- **No high-velocity C IV or N V.**
- **The Silk's high-ion column is unconstrained.** Any contribution below ~10¹³ cm⁻² is hidden in the Galactic column, and the IUE data cannot isolate a Silk component.
- **Scale estimate (derived).** For I(5007) = 2.4 R at T = 3×10⁴ K, n_e·N(O²⁺) ≈ 6.7×10¹⁴ cm⁻⁵. This gives N(O²⁺) ≈ 6.7×10¹⁴ / n_e along the filament; a face-on crossing 2° inside the limb would be ~10× less. I leave this for the modelling agents.
- **Conclusion: no evidence for hot or interface gas tied to the Silk on these sightlines, and no distance constraint.**

## 4. Part (c): distance synthesis
See `figures/verify_distance_constraints.png` and `data/derived/verify_distance_synthesis.json`.

| Constraint | Range | Strength / status |
|---|---|---|
| Monogem limb (south spine at 12.10° from the Knies centre; PSR B0656+14 parallax) | **252–318 pc (68 %), 225–359 (95 %)**, median 282 | Strongest. It is morphological: the X-ray edge registration in the eRASS1 W/E ratio of 1.32 is only 2.4σ. Valid only in the south: the north (Dec > 17.6) bends to 12.4–16° and is not on the Monogem sphere. |
| Monogem shell near/far side (shell up to 1° larger) | 237–337 pc | Allowed variant. |
| Local Bubble wall (O'Neill+2024) | peak 167 pc (inner 139, outer 200) | Reference only. Wu+2026's "160 pc" sits on the wall by coincidence; no dust, HI or X-ray feature ties the Silk to it. |
| HD 61717 −17 km/s cloud, if part of the Silk shell | < 272 (284) pc | Conditional and weak; the velocity does not match the [O III] gas. |
| COS WD J0747+2106 (0.16° from the spine): no Si IV / N V | > 99.5 pc | Conditional: needs Silk log N(Si IV) > 12.2 or N V > 13.0. |
| Ca II / Na I at the [O III] velocity | none | The [O III] velocity equals the local gas, and the [O III] gas is ionized. Ca II/Na I trace low ions in neutral or warm gas and cannot see a Balmer-dark, T ≳ 2×10⁴ K emitter. |
| Kinematic (DESI/SDSS V_LSR +1…+3; LAMOST +11…+16) | 0–1.4 kpc / ≥ 1 kpc | Meaningless: the rotation gradient (7 km/s per kpc) is about σ_pec. |
| Pulsar DMs | none | J0751+1807 (1.18 kpc) is 0.95° E, not behind the filament. B0656+14 (288 pc, DM 13.95) lies inside Monogem. A Silk ΔDM ≈ 0.09–3 pc cm⁻³ for n_e = 10–0.3 is undetectable against DM 30. |
| Extinction | none | The Silk is dust-free. E(B−V) < 0.006 to 300 pc, with no dust feature at 250–320 pc. |
| IUE HD 60848 | none | Only total Galactic high ions. |
| Size / height plausibility | soft upper limit, ~1 kpc | At 1 kpc the spine would be 295 pc long, 4.9 pc wide (FWHM) and at z = 330 pc. At 1.6 kpc it would be 470 pc long. |

**Combined distance.**
- **Favoured: 250–320 pc**, with the south spine on the Monogem limb (68 %: 252–318; 95 %: 225–359).
- **Permitted: ~100 pc to ~1 kpc.**
- 160 pc is *not excluded by any measurement*. It is disfavoured only by the coincidence of the spine with the Monogem X-ray edge.
- ~1.6 kpc (Wu+2026's kd output) is disfavoured by its implausible size and z-height and by the loss of the Monogem registration, not by a direct measurement.
- **No absorption measurement brackets the ionized gas itself.** The only absorber with a distinct velocity, HD 61717's −17 km/s cloud, is not kinematically related to the [O III] gas.

## 5. New findings
1. **HD 61717 is a double-lined spectroscopic binary.** Two A-type components at −69 / +124 km/s on 2023-11-12. APOGEE RVs vary from 35 to 60 km/s over 22 visits. Gaia RUWE is 2.14, so its parallax error is probably underestimated: 272 (262–284) pc.
2. **The −17 km/s cloud is Ti-poor** (Ti II/Ca II < 0.16, against 0.48 locally), cool (b ≈ 1.5 km/s) and compact (absent toward 14 stars ≥ 1.9° away). That argues against the "shocked, grain-destroyed gas" reading.
3. **The pulsar J0751+1807 has an accurate parallax** (0.85 ± 0.04 mas, 1.18 kpc) and an RM (+41.6 rad m⁻²) in the ATNF v2.9.0 database already on disk. Wave 1 used stale values.
4. **IUE toward HD 60848: no high-velocity or N V absorption**, and the high ions are 5–11 km/s blueward of the low ions. The interstellar Si IV, C IV and Al III columns of a 1.9 kpc O-star sightline 2.3° inside the arc are given in row b2.

## 6. Limitations
- **HD 61717:** one epoch only, so stability and circumstellar origin cannot be tested directly. The binary's systemic velocity is uncertain (+28…+44 km/s). b depends on the blue-arm LSF, which I scaled from the red arm.
- **Na I at −17 km/s:** only 2–3σ per line. The telluric template comes from La Silla HARPS stars observed 2006–2014, applied to Paranal 2023.
- **IUE:** resolution ~25 km/s and a NEWSIPS absolute velocity zero point good to ±5–8 km/s. Only the relative high-ion minus low-ion velocity is robust. Only one usable high-ion sightline exists; the B stars and PN central stars are contaminated by photospheric, wind or nebular lines. No FUSE O VI analysis was done (FUSE data exist for NGC 2392 and NGC 2371, but PN contamination is expected).
- **Monogem-limb distance:** it assumes a spherical shell centred near PSR B0656+14's distance. The 15 pc centre-depth offset is my assumption, and so is the Knies+2024 centre and radius as given in `notes/context.md`. The north half of the spine is not on that sphere.
- **COS white dwarf** and **3D dust** values were taken from the wave-1 notes and not re-derived.
- **σ_pec = 7 km/s** in the kinematic test is an assumed typical value.
