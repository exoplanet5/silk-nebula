# verify_hiwham: HI and WHAM Hα at the revised [O III] velocity (wave-2 verification)

Tag `hiwham`. Written 2026-11-16; completed on resume the same day. Everything below was measured with my own code (`scripts/verify_hiwham_*.py`) from the raw files on disk.
- **Wave-1 inputs.** The only wave-1 product I read is the trace geometry, `data/derived/silk_trace.json` v1.0. I read the wave-1 scripts only to look for bugs, and I re-implemented the wave-1 wall finder solely to test it against nulls. Numbers I quote from wave-1 notes, or derive from them, are marked *(wave-1)*.
- **Wave-2 inputs.** Numbers taken from other wave-2 verifiers are marked *(verify_sdss)*, *(verify_desi)* or *(verify_velocity)*. I propagate them and have not re-derived them.
- No new downloads were made.

**Velocity convention.** V_LSR − V_helio = 20 km/s × cos(angle to 18h, +30°, B1900). This is the standard solar motion used by both WHAM and GALFA. Along the spine the correction runs from −13.2 km/s (Dec 12) to −9.3 km/s (Dec 28.25).

The revised [O III] velocity, V_helio +13.5 ± 1 (DESI E strand +13.0, DESI ridge +14.5, SDSS +12.5), corresponds to:
- V_LSR +0.3 at the south end;
- V_LSR +1.5 at P1;
- V_LSR +4.2 at the north end.

This agrees with the wave-2 best value V_helio +13.3 ± 1.0, i.e. V_LSR +1.2 ± 1.0 *(verify_velocity)*.

The W strand was taken at V_helio +18.5 *(wave-1)*, which is V_LSR ≈ +6.6 to +7.1. That velocity is now contested: see §2.4.

**Null philosophy.** Every significance comes from placing the trace at many unrelated positions and measuring the same statistic there:
- **Coherent nulls:** rigid copies of the whole E spine (translated, rotated by ±20° or 180°±20°, and E–W mirrored).
- **Segment nulls:** spine pieces placed independently at random positions and orientations (0–360°).
- **WHAM nulls:** the trace moved in Galactic longitude by 8–66°, at the same latitude band.

Widths are robust: σ_r = (p84 − p16)/2. Pixel or beam bootstraps underestimate the errors by a factor of 2–3 (shown below).

---

## 1. Results table

| # | Claim (wave-1) | Verdict | My value |
|---|---|---|---|
| 7a | GALFA-HI at V_LSR +11…+16 shows an **HI wall 0.86 ± 0.15° east (outside)** of the spine, found in 11/11 segments, wall/floor 1.8 | **REFUTED** (not a detection) | My re-implementation of the wave-1 finder reproduces the segment offsets (mean +0.77 ± 0.17°, 10/11 found; wave-1 +0.88 ± 0.15). The same finder on random sky gives +0.84 ± 0.31° for the 11-segment mean, and a single-segment median of +0.73° (16–84 %: −0.08…+1.84). The real offset sits on the null expectation (z = −0.2). Wall/floor 1.87 vs null 1.58 (P = 0.09). 6/11 segments in 0.4–1.4° (P = 0.06 independent, 0.08 coherent). Independent step fit: +0.87 × 10¹⁹ cm⁻² at d₀ = +0.40°, z = −0.1. Outside−inside step S_loc = +0.81 ± 0.98 × 10¹⁹ (z = +1.0). |
| 7b | At −1…+6 km/s a wall appears in only 6/11 segments, so no comparable structure | **REVISED**: the finder's count of failures carries no information; there is no structure at any velocity | Revised [O III] window (V_helio +10.5…+16.5 ≈ V_LSR −2…+6): S_loc = +0.28 ± 3.05 × 10¹⁹ (z = +0.4), S_wide = −0.49 ± 5.7 × 10¹⁹ (z = +0.4). 2-km/s slices from −15 to +25: all \|z\| ≤ 1.4 against coherent nulls. |
| 7c | No HI filament on the Silk; 2σ ΔN(HI) < 1.3 × 10¹⁹ at +11…+16 | **CONFIRMED**, and extended to the true velocity | +11…+16: A_fil = −0.63 ± 0.55 × 10¹⁹, so UL₂ = 1.1 × 10¹⁹. Revised velocity, whole spine: 2σ UL 1.0–1.4 × 10¹⁹ per 2-km/s slice at V_LSR −3…+5 (south: 1.5–1.75 × 10¹⁹), and ≈ 3.0 × 10¹⁹ for the 6-km/s heliocentric window. On-spine N(HI)(all v) is 4.0 × 10²⁰ for comparison. |
| 7d | "Cavity inside, wall outside" morphology (Monogem boundary in HI) | **REFUTED / not supported** at any velocity | Outside−inside steps are \|z\| ≤ 1.4 in every slice. The largest (+9…+19 km/s, north) is z = +1.0…+1.4. |
| 4a | WHAM Hα main component at V_LSR −0.4…+1.0, FWHM ≈ 30, 2.9–4.4 R, the same on and off filament | **CONFIRMED** | South on (6 beams) / off (22): I_g = 2.59 / 2.75 R; v = +1.28 / +0.54 km/s (beam-to-beam std 0.42 / 1.28); FWHM_int = 30.4 / 31.9 km/s. North on (16) / off (41): 2.90 / 2.80 R; −0.14 / +0.77 km/s; 31.4 / 31.1 km/s. Paired on−off vs 1258 nulls: ΔI_tot = +0.04 ± 0.35 R (S), +0.14 ± 0.27 R (N); Δv = +0.7 ± 1.8 (S), −0.9 ± 1.7 (N) km/s; ΔFWHM = −1.6 ± 3.0 (S), +0.4 ± 2.6 (N) km/s. All \|z\| < 0.7. |
| 4b (implied) | WHAM "Hα/[O III] < 0.10 (2σ) at the [O III] velocity" (spectra_ism; it used V_LSR +13.6) | **REFUTED at the revised velocity**: the limit no longer applies | The [O III] now sits on top of the WIM line: on-filament WIM V_LSR +1.3, [O III] +1.2 ± 1.0 *(verify_velocity)*. Template-regression 2σ UL on the *intrinsic* filament Hα peak (FWHM 0.276°): **1.0 R** south in a ±25 km/s window (corrected for 86 % capture of a 40 km/s-FWHM line), **1.25 R** south for total Hα, and 1.25–1.30 R north. Beam dilution D = 0.11–0.36 (median 0.24). |
| — | spectra_ism: on−off WHAM residual of +0.089 ± 0.028 R at V_LSR −9.8 (3σ, bootstrap) | **REVISED** to ≈ 1σ | Window −30…−2 km/s, on−off: north +0.13 ± 0.13 R (z = +1.2); south −0.02 ± 0.16 R. Beam bootstraps underestimate the WIM-dominated error by 2.3–2.6× (south: 0.15 vs 0.35 R; north: 0.10 vs 0.27 R). |
| 7e | Finkbeiner (2003) Hα: no filament, 2σ < 0.43 R (beam-corrected for 6′) | **REVISED** | No filament, confirmed, but the limit is weaker. The F03 bitmask (bit definitions read from the FITS header) puts SHASSA under only **24 %** of the spine (Dec ≤ 15.75 at the trace nodes; 15.94 on the densified trace). VTSS covers 0 % of the spine; it reaches only 2–12 % at d > +1.25° at the north end. Everything else is **WHAM-only, at 1° resolution**, so a 6′ correction is invalid there. SHASSA part (Dec 12–15.6, 6′): A = +0.18 R, z_r = +0.9 vs 400 placement nulls, **2σ UL 0.78 R** on the intrinsic peak (1.19 R using the plain std, which is inflated by two outlier placements). WHAM-only part (Dec 15.9–28.25): z_r = +0.4…+0.7. Its 2σ UL is 0.30 R if treated as 6′ (the wave-1 error) but **1.74 R** with the true ≈ 1° beam. |
| 7f | SHASSA: 3σ < 0.8–1.0 R in the south (imaging_optical) | **REVISED** (no detection; limit about 30 % weaker) | SHASSA Sm (4′, Dec 12–15.6): A = +0.07 R, z_r = +0.39 vs 400 placement nulls, **2σ UL 0.77 R** on the intrinsic Hα+[N II] peak (3σ 1.1 R). SHASSA CC (unsmoothed, binned medians): A = +0.29 R, z_r = +1.1, 2σ UL 1.02 R (3σ 1.34 R). Why weaker: the placement-null σ of the amplitude is 0.31 R. Sky structure is correlated, so σ does not fall with piece length; the true 3σ is therefore 1.1–1.3 R, not 0.8–1.0. **New:** displaced or broader Hα is also excluded. With the template centred anywhere in d₀ = −0.8…+0.8° (FWHM 0.28° or 0.5°), the 2σ UL is ≤ 0.85 R at every d₀. Largest z = +0.55 (global p = 0.82); deepest dip z = −2.3 at d₀ = −0.7 (global p = 0.15). |
| 8a | GS202+17+019 and GS202+19−019 HI shells are centred on the Silk | **REVISED** | GS202+17+019 is a modest cavity in LAB at +10…+40 km/s (rim−centre z = +1.7, +2.7, +2.8 vs ellipses shifted in l). GALFA N(HI) at v₀ ± 4 km/s rises from 3.1 × 10¹⁹ (ρ < 0.4) to 5.0–5.5 × 10¹⁹ (ρ > 0.8). It is **not** at the revised [O III] velocity: there is no central deficit at V_LSR −2…+4. The spine cuts through its interior: ρ = 0.55–0.6 at Dec 15.75–17, rising to 2.2–2.4 at Dec 28, for either orientation convention. Only 23–33 % of the spine lies at ρ 0.8–1.2. GS202+19−019 (−18.8 km/s) is not significant in LAB (rim−centre z ≤ 0.7 at −20…−10, ≤ 1.4 elsewhere). The paper's [O III] velocity (+11…+16) coincided with the GS202+17+019 range; at the revised velocity that coincidence disappears. |
| 8b | PGCC G199.09+19.95, 0.35° from the Hα patch | **CONFIRMED, plus new HI detail** | GALFA shows two narrow (CNM) HI lines at the clump (r < 0.1°, minus a 0.3–0.8° annulus): V_LSR −0.74 and +6.46 km/s, FWHM 2.7 and 3.0 km/s, ΔT_B ≈ 7.0 and 6.9 K. Formal fit errors are ±0.05–0.1 km/s and ±0.2 K; the real uncertainties are larger. ΔN(HI)(all v) = +5.4 × 10¹⁹ cm⁻². Only 2.7 % of 2000 random positions show a narrow (≤ 5 km/s) excess this strong. The Hα-patch polygon shows the same pair, at −0.5 and +5.0 km/s (ΔT 2.5 and 4.4 K). |
| 9 (ISM link) | Local Ca II/Na I gas at V_helio +6…+15, the same as WHAM | **CONFIRMED** in HI | On-spine HI peaks at V_LSR +0.4…+5.5 (median +1.8) for Dec 12–22. That is V_helio ≈ +13…+14, equal to the revised [O III] velocity. |

---

## 2. HI (GALFA-HI DR2, 2′ mosaic)

### 2.1 Data and method
- **Mosaic** (`verify_hiwham_common.galfa_cube`):
  - The 9 raw DR2 *wide* cutouts are placed on the 1′ lattice by index arithmetic, with a lattice-alignment assertion, and then 2×2 block-averaged.
  - Channels are 0.736 km/s, V_LSR −29.8…+40.1.
  - The RA 119.33–125, Dec 10.5–30 tiles are missing (HTTP 429 in wave 1). East-side profiles beyond about +2.3° in the south therefore have fewer pixels. Nulls require at least 75–80 % coverage.
- **Maps** (`verify_hiwham_hi.py`):
  - 20 LSR slices of 2 km/s each, from −15 to +25.
  - A per-pixel heliocentric window V_helio +10.5…+16.5 (the revised [O III]) and a W-strand window, +16…+21 helio.
  - The wave-1 windows (+11…+16 and −1…+6 LSR), and all velocities.
  - Channels are weighted fractionally at the window edges.
- **Profiles:** 0.05° bins in signed offset d (d > 0 = east/convex side), for |d| < 3°. They are computed for the whole spine, the south (spine Dec ≤ 17.6), the north, and 11 segments of 1.5°.
- **Statistics:**
  - Gaussian-filament amplitude A_fil (FWHM 0.284° = the 0.276° filament convolved with the 4′ beam) over a linear baseline fitted to |d| < 1.2°.
  - Core contrast.
  - Outside−inside steps: S_loc (0.2–1.2°) and S_wide (0.5–2.5°).
  - A smoothed-step (erf, σ = 0.1°) edge fit scanned over d₀ = −1.5…+2.0°.
  - My re-implementation of the wave-1 half-rise wall finder.
- **Nulls:**
  - 120 coherent whole-trace copies, chosen from 2573 accepted placements, all overlapping the real |d| < 1.5° band by less than 5 %.
  - 1500 independent segment placements.

### 2.2 The +11…+16 km/s "wall" is a property of the finder
- **How the finder works.** It takes the floor as the minimum in −1…+1°, the peak as the maximum in 0…+3°, and returns the first mid-level crossing between them.
- **What it returns.** By construction it gives a position between floor and peak, typically +0.7…+0.9°. On random sky it finds a "wall" in 75 % of single segments, with median offset +0.73°.
- **Reproduction.** My re-implementation matches the wave-1 per-segment offsets to about 0.1°:
  - mine: [n/a, 1.51, 1.32, 0.62, 0.03, 1.29, 1.26, 0.71, 0.40, 0.62, −0.08];
  - wave-1: [1.60, 1.49, 1.21, 0.68, −0.04, 1.33, 1.26, 0.86, 0.46, 0.60, 0.27].
- **The mean is not evidence of coherence.** The 11-segment mean of +0.77 ± 0.17 (wave-1 +0.88 ± 0.15) sits on the null distribution of the same quantity, +0.84 ± 0.31. The wave-1 "± 0.15" is a standard error of 11 values that scatter by 0.5°.
- **The only hint of anything real** is the clustering of offsets in 0.4–1.4°:
  - 6/11 segments at +11…+16 (P = 0.06–0.08);
  - 7/11 at +13…+15 (P = 0.013 independent, 1/120 coherent).
  - With a look-elsewhere factor of about 20 slices, this is not significant.
- **Mean profile.** The mean whole-spine profile at +11…+16 rises by about 1 × 10¹⁹ from inside to outside. That rise lies inside the 16–84 % coherent-null band (`figures/verify_hiwham_hi.png` panel d).

### 2.3 At the revised velocity: no wall, cavity or filament
- **True-velocity window** (V_LSR −2…+6, tracking heliocentric +13.5 per pixel):
  - A_fil = +0.09 ± 1.45 × 10¹⁹ (z = +0.1);
  - south +1.10 ± 1.62 × 10¹⁹ (z = +0.6);
  - S_loc = +0.28 ± 3.05 × 10¹⁹;
  - steps and edges are consistent with the nulls. The step-amplitude z is −1.4, i.e. *weaker* than typical random edges.
- **Per 2-km/s slice:** whole-spine A_fil z = +0.1 (−1…+1), +0.2 (+1…+3) and −0.3 (+3…+5). The largest |z| over all 20 slices is 1.4 (A_fil, +15…+19 km/s, a slight HI deficit on the spine).
- **Segment-combined statistics** (Σz/√11): A_fil −2.8 at +15…+19 km/s, driven by segment 9 (Dec ≈ 25.6, z = −3.5). This is not at the [O III] velocity, and it is not significant after the trials factor.
- **2σ upper limits on an HI filament** (peak ΔN(HI) at the filament scale): 1.0–1.4 × 10¹⁹ cm⁻² per 2-km/s slice at V_LSR −3…+5 for the whole spine, and 1.5–1.75 × 10¹⁹ for the south (panel f). That is ≲ 3 % of the local N(HI) in this velocity range (~1.0 × 10²⁰ in the 6-km/s window).

### 2.4 (c) Velocity structure along the filament, and the W strand
`verify_hiwham_hi_extra.py`; `figures/verify_hiwham_hi_extra.png`.

**On-spine HI** (|d| < 0.15°):
- Moment-1 (−15…+30 km/s) runs from about +4 km/s (south) to about +6 km/s (north). The slope is +0.14 ± 0.03 km/s deg⁻¹ along Dec.
- The peak-velocity slope is +0.61 ± 0.12 km/s deg⁻¹. It is driven by a +8…+12 km/s component that becomes dominant north of Dec ≈ 23 (Gaussian decomposition: +8.7…+11.3 km/s, FWHM 10–16 km/s, 1.2–3.2 × 10²⁰ cm⁻²).
- For Dec 12–22 the dominant component is at V_LSR +1.2…+2.5 (FWHM 10–18 km/s), with narrow CNM pieces at +1…+2 km/s (FWHM 4–6 km/s).

**Comparison with [O III].**
- At constant V_helio = +13.5 the [O III] maps to V_LSR +0.3…+4.2. That slope of about 0.24 km/s deg⁻¹ is pure solar reflex.
- In the south, the revised [O III] velocity coincides with the dominant local HI, the WHAM WIM and the Ca II/Na I "local" group.
- In the north the HI shifts redward by about 6–8 km/s. No [O III] velocity data exist there to compare with; the DESI/SDSS velocities are southern.

**No HI velocity signature on the spine.** On − off (0.4–1.2°) moment-1 = −0.16 ± 0.73 km/s (scatter of 34 bins of 0.5°). The difference spectra (panel b) show ±2–3 K channel structure, with no coherent feature along the [O III] track.

**W strand.** Its velocity is now disputed between wave-2 verifiers:
- DESI positive fibres: +13.7 ± 1.5 overall and +16.7 ± 2.1 at Dec 15.5–17.5 *(verify_desi)*;
- SDSS: +9.7 ± 2.2 *(verify_sdss)*;
- wave-1 DESI regression: +18.5. verify_desi shows that regression centroids of local subsets are not local.

So the W strand lies somewhere in V_helio ≈ +10…+19, i.e. V_LSR ≈ −1…+8. I tested both windows that cover it, along the W strand's own polyline (A_fil, 600 random placement nulls):
- V_helio +16…+21: −0.43 × 10¹⁹ (z = −0.55), 2σ UL 1.4 × 10¹⁹;
- V_helio +10.5…+16.5: z = −0.65;
- +11…+16 LSR: z = −0.76, UL 1.0 × 10¹⁹;
- all velocities: z = −0.67.

There is no HI along the W strand at any of its candidate velocities.

### 2.5 PGCC G199.09+19.95 and the Hα patch
See table row 8b.

**Velocity coincidence with the W strand (now weaker).**
- If the W strand were at V_helio +18.5, it would be at V_LSR +7.1 at the patch, within about 1 km/s of the +6.46 km/s CNM line.
- Even at that velocity, 12–13.5 % of random GALFA positions have their HI-excess peak within ±1.5 km/s of +7.1, so the coincidence is not significant (about 1.1σ).
- At the wave-2 W-strand values (V_helio +9.7…+16.7, i.e. V_LSR ≈ −1.6…+5.4) the W strand lies *between* or *on* the −0.7 and +6.5 lines. So this coincidence test cannot discriminate.

**Column and temperature of the clump.**
- ΔN(HI) ≈ 5 × 10¹⁹ cm⁻² from the two CNM lines.
- FWHM 2.7 km/s gives T_k ≤ 21.9 × 2.7² ≈ 160 K *(derived)*.
- This is a genuine cold clump at the western tip of the Hα patch. Its two velocities bracket the local gas.

### 2.6 Catalogued HI shells
`verify_hiwham_hi_extra.py` §D uses LAB 10-km/s HEALPix slices plus GALFA.

**Method.** Velocity–radius diagrams (panels d, e), and LAB rim−centre contrasts compared against 66 ellipses shifted by Δl = ±15…80° and Δb = 0, ±4°.

**Parameters.**
- GS202+17+019: centre (202.4, +16.9), semi-axes 5.7 × 4.5°, φ = 29.6°, v₀ = +18.6 km/s.
- GS202+19−019: centre (202.3, +18.8); Δl × Δb = 6.5 × 5.5°, treated both as full sizes and as semi-axes; v = −18.8 km/s.

**Results.** See table row 8a. The "HI minimum" listed for GS202+17+019, at (l, b) = (203.3, 19.1) = RA 116.10, Dec 16.85, lies only 0.54° west of the spine at Dec 17. It is inside the arc (the concave/Monogem side).

---

## 3. Hα: WHAM, SHASSA and Finkbeiner

### 3.1 WHAM method (`verify_hiwham_wham.py`, `verify_hiwham_fig_wham.py`)
- **Data:** WHAM-SS DR1 v161116.
  - 9129 beams at l 125–285, b −4…+50.
  - Spectra are interpolated to a 2-km/s grid.
  - The STAR flag excludes no beam here.
- **Per-beam measures:**
  - INTEN (−80…+80 km/s);
  - the window v_c ± 25 km/s (v_c = V_helio +13.5, converted per beam);
  - a narrow core window, v_c ± 12;
  - a blue window, −30…−2;
  - baseline-subtracted moments;
  - a single-Gaussian + linear-baseline fit of the main component.
- **On/off:** on = |d| < 0.5°; off = 1.5–3° on both sides; matched within |Δs| < 0.75°.
- **Template regression:**
  - D_b = the beam-averaged unit-peak Gaussian filament (FWHM 0.276°) in a 1° top-hat beam (437-point sunflower sampling);
  - the model is y_b = quadratic(d, s) + I_peak · D_b, fitted over beams with |d| < 3.5°.
- **Nulls:** 1258 copies (Δl = ±8…66°, Δb = −6…+9°, rotations 0/180/±20°, mirrors; centre b = 12–30°). They are correlated, so I use robust widths.
- **Injection test:** a 1-R-peak filament (FWHM 40 km/s at v_c) is recovered with I_peak = 1.000 (total) and 0.859 (±25-km/s window, equal to the analytic capture fraction). In the on−off difference it gives +0.31 R (south) and +0.23 R (north).

### 3.2 WHAM results
- **Same on and off.** The main WIM component is identical on and off the filament in intensity, centroid and width (table row 4a).
- **Main-component velocity.** In the south the on-filament V_LSR is +1.3 (beam std 0.4) km/s; off-filament it is +0.5 (std 1.3). This equals the revised [O III] V_LSR (+0.3…+1.5 at Dec 12–17; +1.2 ± 1.0 *(verify_velocity)*).
  - The Hα/[N II]/[S II] (WIM) and the [O III] are therefore **no longer kinematically distinct**.
  - The paper argued the red lines are a separate DIG because their velocity differs. That rests on the LAMOST [O III] velocity, which is about 10–12 km/s too red (wave 1; verify_velocity). Any distinction now rests only on spatial morphology and line ratios.
- **Filament Hα peak limits (2σ):**
  - south: 1.0 R (window) and 1.25 R (total);
  - north: 1.25 R and 1.30 R.
  - These are about 10× weaker than the wave-1 approach implied at the old velocity. The WIM line, not photon noise, sets the limit: its beam-to-beam fluctuations are 0.25–0.35 R at 1° (in on−off) and 0.43–0.62 R per unit template (σ of I_peak).
- **Narrower line?** A filament line narrower than the WIM would raise the on-beam peak and narrow its width. The on−off widths (ΔFWHM = −1.6 ± 3.0 and +0.4 ± 2.6 km/s) exclude only gross effects. The injection shows that a 1-R-peak filament changes the fitted FWHM by just +0.4…+0.5 km/s.

### 3.3 SHASSA and Finkbeiner (`verify_hiwham_hamaps.py`, `verify_hiwham_shassa_offset.py`, `verify_hiwham_f03bits.py`)
- **Units.** BUNIT = dR, confirmed in the SkyView header; I divide by 10.
- **SHASSA coverage.** SHASSA covers Dec < +15 nominally; the cutout has data to 15.98, so the analysed piece is spine Dec 12–15.6. That is 3.6° of the 16.9° spine.
- **Filter.** The 32 Å filter passes Hα plus partially transmitted [N II] 6548/6583. The limits therefore apply to Hα + τ[N II], which is conservative for Hα alone.
- **Statistic.** The amplitude of a Gaussian (FWHM = 0.276° ⊕ beam) over a linear baseline, on 0.05°-binned profiles. The intrinsic peak is A × FWHM_obs/0.276.
- **Nulls.** The same spine piece at random centres, orientations (0–360°) and mirrors within the footprint, avoiding the real band. 400 nulls (Sm, F03), 200 (CC), 500 (offset scan).
- **Displaced-Hα scan.** In a radiative shock the H recombination zone lies downstream of the [O III] zone. Every aperture-matched spectroscopic limit (DESI regression on an [O III] template; SDSS ON with r_on = 0.15°) is blind to an Hα ridge displaced by a few tenths of a degree.
  - I scanned the template centre over d₀ = −0.8…+0.8°, at FWHM 0.28° and 0.5°.
  - No offset shows a significant signal.
  - The 2σ UL is ≤ 0.85 R at every d₀, and 0.80 R for the 0.5° width.
  - 0.8° corresponds to 3.5 pc at 250 pc and 4.2 pc at 300 pc *(derived)*.
- **F03 provenance.** Bits are taken from the header: 0 WHAM, 1 VTSS, 2 SHASSA, 3 STAR. Along the spine: WHAM 100 %, SHASSA 24 % (Dec ≤ 15.75), VTSS 0 %, STAR 0 %. North of Dec ≈ 15.9 the composite is a 1° WHAM interpolation. Any arcminute-scale Hα statement from F03 there is invalid, including:
  - the imaging_mw "< 0.43 R" over the whole spine;
  - the imaging_optical "+0.8–0.9 R at Dec 17–21", which is a 1°-scale WHAM feature, consistent with the WHAM Hα-patch beam.

---

## 4. Best lower limit on [O III]5007/Hα for the south strand

**Units and [O III] scale.**
- All ratios are photon (Rayleigh) ratios; multiply by 6563/5007 = 1.31 for energy units.
- For the imaging and WHAM routes, the ratio depends on the assumed [O III] *peak*. I give it for the task's 2–4 R *(wave-1)* and for the wave-2 revised DESI south ridge median, **1.6 R** *(verify_desi; E strand 1.4–2.3 R at Dec 14–18)*.
- The aperture-matched spectroscopic routes do not need the [O III] peak.

**Beam dilution.**
- WHAM 1° beams: D = 0.11–0.36 (median 0.24). The template regression already includes it, so the WHAM numbers are limits on the intrinsic peak.
- SHASSA 4′: D = 0.97. F03 6′ (SHASSA part): D = 0.94.
- SDSS and DESI fibres: no dilution.

| Route | Aperture / assumption | Hα constraint (2σ) | [O III]/Hα lower limit |
|---|---|---|---|
| WHAM template, ±25 km/s at the [O III] velocity (this work) | 1°, D ≈ 0.24; filament profile; velocity-resolved | I_peak < 1.00 R | > 1.6 (1.6 R) / 2.0 (2 R) / 4.0 (4 R) |
| WHAM template, total Hα (this work) | as above | < 1.25 R | > 1.3 / 1.6 / 3.2 |
| SHASSA Sm, d₀ = 0 (this work) | 4′, D = 0.97; Hα + [N II] | < 0.77 R | > 2.1 / 2.6 / 5.2 |
| SHASSA CC, binned medians (this work) | as above | < 1.02 R | > 1.6 / 2.0 / 3.9 |
| F03, SHASSA part (this work) | 6′, D = 0.94 | < 0.78 R | > 2.1 / 2.6 / 5.1 |
| **SHASSA Sm, Hα displaced \|d₀\| ≤ 0.8°, FWHM 0.28–0.5° (this work)** | **no co-spatiality assumed** | **< 0.85 R** | **> 1.9 / 2.3 / 4.7** |
| SDSS ON−OFF, Dec < 17.6 *(verify_sdss inputs; placebo rms added here)* | same fibres; r_on 0.15° | Hα 0.08 ± 0.22 R (boot 0.12 ⊕ placebo 0.18) vs [O III] 1.57 R | **> 3.0** |
| SDSS Hβ + case B, Dec < 17.6 *(verify_sdss inputs)* | same fibres | Hβ 0.03 ± 0.11 R → Hα_rec < 0.99 R | > 1.6 |
| SDSS ON−OFF, Dec < 20 *(verify_sdss)* | includes Dec 17.6–20, near the Hα patch | Hα 0.52 ± 0.21 R vs [O III] 1.39 R | > 1.5 (central 2.7) |
| DESI Hα, IV regression *(verify_desi)* | [O III]-template-matched | Hα/[O III] (energy): all-fibre p97.5 0.226 (skewed); south 0.073, p97.5 0.082 | **> 3.4** (all, raw p97.5) to > 8.8 (south, bootstrap × 1.5); nominal south 10.5 |
| DESI Hβ + case B *(verify_desi)* | [O III]-template-matched | [O III]/Hβ (energy) > 18 (smallest p2.5 over 5 subsets) | **> 4.8** (Hα/Hβ = 2.86 energy); > 3.9 if Hα/Hβ ≤ 3.5 (collisional enhancement) |
| wave-1 SDSS Dec < 20 (superseded) | — | Hα 0.17 ± 0.03 R | > 3.4 |

**Best lower limit: [O III]5007/Hα ≳ 3 (2σ, photon units; ≳ 4 in energy units) for the south strand**, with Hα co-spatial with the [O III]. Three independent aperture-matched routes all clear 3:
- SDSS Dec < 17.6 with the placebo floor: > 3.0;
- DESI direct Hα, most conservative tail: > 3.4;
- DESI Hβ + case B: > 3.9–4.8.

The strongest single statement is DESI Hβ + case B, **> 4.8 (> 6.3 energy)**, if Hα/Hβ ≤ 2.86. The DESI nominal value is about 10 photon (13.7 energy).

**Hα not co-spatial with the [O III].** This includes a recombination zone displaced by up to 0.8°, or one twice as broad. The only valid limit is then SHASSA's, **> 1.9 (for [O III] 1.6 R) to > 2.3 (2 R) to > 4.7 (4 R)**. WHAM gives ≥ 1.6–2.0 at 1.6–2 R and is velocity-resolved, but it is limited by WIM confusion at the same velocity.

**Change from the previous draft.** Wave-1's "≳ 3.4 from SDSS" survives numerically, but only through the Dec < 17.6 subset with placebo errors. The Dec < 20 SDSS Hα is 3× larger, 0.52 R *(verify_sdss)*, and probably picks up Hα structure at Dec 17.6–20 near the patch. The wave-1 "> 13.8 [O III]/Hβ" is weakened to > 7 by the placebo floor *(verify_sdss)*.

---

## 5. New findings
1. **The wave-1 HI wall is an artefact of the edge finder.** The finder returns about 0.8° on random sky. The "0.86 ± 0.15°" is a standard error of values that scatter by 0.5°. The RESUME headline "with an HI wall just outside" should be dropped.
2. **The revised [O III] velocity coincides with the dominant local HI and the WIM.** On-spine HI peaks at V_LSR +1…+2.5, the WHAM WIM sits at +1.3, and the local Ca II/Na I group at V_helio +6…+15. Two consequences:
   - The ionized filament is at rest relative to the local ISM along the line of sight, consistent with a tangential (limb) geometry.
   - The paper's kinematic separation of [O III] from the red lines no longer holds.
3. **The WHAM limit on filament Hα degrades about 10× at the revised velocity.** The "Hα/[O III] < 0.10" result (spectra_ism) is no longer valid. The best [O III]/Hα constraints now come from aperture-matched DESI/SDSS spectra (≳ 3–5) and, without the co-spatiality assumption, from SHASSA (≳ 2).
4. **Finkbeiner (2003) is WHAM-only along 76 % of the spine.** Its header-defined bitmask shows SHASSA only at Dec ≤ 15.75–15.94, and VTSS nowhere on the spine. Any 6′ beam correction applied over the whole spine is invalid. The true northern F03 limit is 2σ < 1.7 R, not 0.43 R. This affects imaging_mw's "< 0.43 R", and any F03 value north of Dec 15.9 (including the "+0.8–0.9 R" at Dec 17–21, which is 1°-resolution WHAM).
5. **No displaced or broadened Hα ridge** near the southern [O III] filament. SHASSA excludes a > 0.85 R (2σ) Hα+[N II] ridge anywhere within ±0.8° (≈ 3.5–4.2 pc at 250–300 pc) of the spine, for widths 0.28–0.5°. This closes the loophole in the fibre-based limits, at the 2σ ≈ 0.8 R level.
6. **The cold clump at the Hα patch has two narrow HI lines** (−0.7 and +6.5 km/s, T_k ≲ 160 K). The earlier +6.5 coincidence relied on a W-strand V_helio of +18.5. That velocity is now disputed (SDSS +9.7, DESI +13.7…+16.7), so the test is inconclusive.
7. **GS202+17+019 is real-ish but at the wrong velocity.** It is a weak LAB cavity (z = 1.7–2.8 at +10…+40 km/s), and the spine crosses its interior. GS202+19−019 is not confirmed in LAB.
8. **Beam bootstrap errors underestimate the WHAM on−off uncertainty by 2.3–2.6×** relative to placement nulls. The wave-1 WHAM significances (e.g. the "3σ" −9.8 km/s residual) should be read as about 1σ. The same applies to SHASSA: the wave-1 3σ of 0.8–1.0 R becomes 1.1–1.3 R with placement nulls.

## 6. Limitations
- **GALFA coverage.** The eastern tiles (RA > 119.33 at Dec > 10.5) are missing, so the southern east-side profiles stop at about +2.3–2.8°. The coherent nulls are confined to the GALFA footprint and partly overlap each other: 120 copies from about 25 effectively independent footprints. The independent-segment nulls (1500) complement them.
- **Optical depth.** N(HI) is computed optically thin. The CNM lines at the PGCC (T_B up to 17 K) are mildly optically thick, so their N(HI) is a lower limit.
- **WHAM resolution.** The 1° beams sample the 0.28° filament with only 1–2 beams across it. The template assumes the E-spine Gaussian (FWHM 0.276°). The W strand and the Hα patch are not in the template; the patch lies about 1.3° west, near the off band in the north.
- **WHAM nulls.** The 1258 copies are correlated (overlapping footprints and shared WIM structure). Their tails include copies that cross bright H II regions, which is why I use robust widths. The near-latitude subset (|Δl| ≤ 36°, |Δb| ≤ 6°) gives similar widths (see `verify_hiwham_wham_summary.json`).
- **Line width.** The assumed filament Hα line width (40 km/s FWHM, for a 32 km/s [O III] FWHM_int that is mostly non-thermal, plus thermal H broadening) sets the window capture fraction: 0.78–0.93 for FWHM 48–32 km/s.
- **SHASSA** covers only Dec 12–15.6 (21 % of the spine). North of that the only arcminute Hα data on disk are NSNS *(wave-1)*. SHASSA continuum-subtraction residuals and field seams are absorbed in the null σ, not modelled.
- **SDSS/DESI inputs.** These are the wave-2 verifiers' products, propagated not re-derived. The SDSS placebo rms (0.18 R Hα, 0.09 R Hβ) was computed for the southern region with ≥ 20 ON fibres. The Dec < 17.6 subset has 54 ON fibres and may have a larger floor. The DESI ratio bootstraps are skewed (a few blocks dominate), hence the range > 3.4…8.8. The case-B route assumes Hα/Hβ ≤ 2.86 (energy); collisional excitation in a shock could raise this to about 3–3.5, giving > 3.9.
- **Shell parameters.** I did not verify against the papers the orientation convention of φ (Suad+2014), or whether Ehlerová & Palouš Δl, Δb are radii or diameters. I tested both choices. LAB slices are 10 km/s wide with a 36′ beam.
- **No new data** were downloaded. HI4PI and EBHIS cubes are not on disk (only N(HI) maps, which GALFA supersedes).

## 7. Files
- **Scripts** (`scripts/`):
  - `verify_hiwham_common.py` (geometry, rigid trace transforms, LSR, GALFA mosaic)
  - `verify_hiwham_hi.py`, `verify_hiwham_hi_summary.py`, `verify_hiwham_hi_extra.py`, `verify_hiwham_hi_extra2.py`
  - `verify_hiwham_wham.py`
  - `verify_hiwham_hamaps.py` (SHASSA / F03 limits, d₀ = 0)
  - `verify_hiwham_shassa_offset.py` (displaced/broad Hα scan; also merges the hamaps part files)
  - `verify_hiwham_f03bits.py` (F03 survey bitmask along the spine)
  - `verify_hiwham_fig_hi.py`, `verify_hiwham_fig_hi_extra.py`, `verify_hiwham_fig_wham.py`
- **Derived products** (`data/derived/`):
  - `verify_hiwham_hi_profiles.npz`, `verify_hiwham_hi_results.json`
  - `verify_hiwham_hi_extra.{json,npz}`, `verify_hiwham_hi_extra2.json`
  - `verify_hiwham_wham_{beams.csv,results.json,spec.npz,summary.json}`; the summary holds `oiii_ha_lower_limits`
  - `verify_hiwham_hamaps_results.json` (merged), `verify_hiwham_hamaps_part_{Sm,CC,F03}.{json,npz}`
  - `verify_hiwham_shassa_offset.{json,npz}`, `verify_hiwham_f03bits.json`
- **Figures (200 dpi):**
  - `figures/verify_hiwham_hi.png`
  - `figures/verify_hiwham_hi_extra.png`
  - `figures/verify_hiwham_wham.png`; panel f is the [O III]/Hα lower-limit summary
- **Cache:** the 2′ GALFA mosaic is kept outside the project, in the session scratchpad (`hiwham/galfa_mosaic_2arcmin.npz`). It is rebuilt automatically from the raw cutouts.
