# w3_judge_phys: origin hypotheses for the Silk Nebula, ranked for astrophysical consistency

Wave 3, judge panel, tag `judge_phys`. Written 2026-11-17. I am one of three independent judges. My lens is physical self-consistency: energetics, pressure balance, timescales, geometry, Galactic context and analogue objects.

**What I read**
- All wave-1 and wave-2 notes.
- `w3_canon.md` with `w3_canon_values.json`, plus `w3_oii.md`, `w3_northsouth.md` and `w3_figs.md`.
- From the literature on disk: Fesen+2026 (`data/raw/context_arxiv/fesen2026_2608.26666v1.txt`) and Knies+2024 (`…/2401.17289.txt`).
- Two web searches, for Antlia and Nereides only. These are marked *web, not re-verified*.

**Precedence.** Wave 2 supersedes wave 1. Where `w3_oii` and `w3_northsouth` revise wave 2, I adopt their values; three of their findings change the physics:
1. **[O II] doublet.** 3729/3726 = 1.31 ± 0.07 from DESI, within 2σ of the low-density limit. The SDSS value of 1.29 is an instrumental artefact, so **n_e ≈ 100 cm⁻³ is not established**.
2. **[O II]7325.** The wave-1/2 flux was low by a factor of 1.8, a normalisation bug. The corrected 7325/3727 is **0.060 (0.055–0.073), ± 0.015 from OH**, which gives T(O II) ≈ 3.5 × 10⁴ K, not 1.25 × 10⁴ K.
3. **North versus south.** The north is a *slower* shock (70–100 km/s) at lower driving pressure, P_N/P_S ≈ 0.45–0.7. The DESI/SDSS [O II]/[O III] conflict was a Dec-weighting effect. The NVSS RM step is not specific to the filament (placement p = 0.18). The southern spine is *tangent to*, not concentric with, the Monogem circle.

**Labels and files**
- Labels: **derived here** marks my own arithmetic, from `scripts/w3_judge_phys_checks.py` → `data/derived/w3_judge_phys_checks.json`. Other numbers are cited from the named teammate note.
- Figure: `figures/w3_judge_phys_posteriors.png`.

---

## 0. Verdict

| # | Hypothesis | Posterior | Verdict |
|---|---|---|---|
| H1 | Incomplete (column-truncated) radiative shock driven by the **Monogem Ring** hot gas, at ~250–320 pc. This includes **H1b**, the Monogem blast wave striking an external, flatter sheet; see §2.4 | **0.64** | **favoured** |
| H2 | Incomplete radiative shock driven by a **different** remnant or bubble: G205.6+12.4, a Geminga-progenitor SNR, an unknown old SNR, or the Local Bubble wall. Distance unconstrained | **0.26** | viable |
| H4 | Superbubble or chimney wall at 1–2 kpc (shock-excited), taking the paper's kinematic distance at face value | 0.04 | disfavoured |
| H6 | Other non-shock mechanisms: conduction front, turbulent mixing layer, NEI cooling of hot gas, photoionisation by a hot WD or Monogem EUV, fossil | 0.025 (summed) | disfavoured; most members ruled out individually |
| H9 | *Added:* unanticipated mechanism or a combination | 0.012 | — |
| H5 | DIG or photoionised ionised-bubble wall (the paper's stated alternative) | 0.01 | ruled out |
| H3 | The paper's specific model: independent SNR at 160 pc, R ≈ 30 pc, E ≈ 2 × 10⁵¹ erg, 6.8 × 10⁴ yr, PDS phase | 0.005 | ruled out |
| H7 | *Added:* bow shock or wake of a runaway star or pulsar | 0.005 | ruled out (in practice) |
| H8 | *Added:* intermediate- or high-velocity cloud (IVC/HVC) interaction or infall shock | 0.003 | ruled out |

**Main points**
1. **The excitation question is effectively settled.** P(shock-excited, any driver) ≈ 0.95.
2. **The driver and distance are not.** The Monogem identification (H1) is preferred by about 2.5 : 1 over all other drivers combined (H2). That preference rests on morphology and on where the driving hot gas must sit (§2.5). It does not rest on a measured distance.
3. **The paper's quantitative model (H3) is excluded by more than a factor of 20 in three independent observables** (§2.6). Its 160-pc distance came from a unit error. Its R ≈ 30 pc at 160 pc is 10.7°, which is essentially the Monogem angular radius placed at the wrong distance.

---

## 1. Evidence base I used (canonical values)

**Excitation**
- [O III]/Hβ > 11 (2σ, DESI); SDSS ~23 (13–48).
- [O III]4363/5007 = 0.06–0.09 (2–3σ), so T(O III) ≈ 3–6 × 10⁴ K.
- 7325/3727 = 0.060 ± 0.015 (`w3_oii`), so T(O II) ≈ 2.4–5.5 × 10⁴ K.
- He II/[O III] < 0.026.
- [O II]/[O III]: 1.1–1.2 south, rising smoothly to 4–8 north.
- [O II]3729/3726 = 1.31 ± 0.07.
- GALEX FUV 68–92 CU per R and NUV 36–57 CU per R (south). No Hα (< 0.8–1.0 R, 2σ), HI (< 1–3 × 10¹⁹ cm⁻²), dust, X-ray or radio filament.

**Kinematics**
- V_LSR(south) = +1.2 ± 1.0 km/s, the same as the WIM and the local HI.
- North is redder by ΔV_LSR ≈ +7…+11 km/s. This is estimator-dependent, and DESI and SDSS disagree at 2–3σ at Dec 17.6–21.
- [O III] FWHM: DESI ≤ 33 km/s (not re-verified) against ~50 km/s in the paper's Fig. 3a.

**Brightness and geometry**
- South ridge ≈ 1.6 R (1.4–2.5); north 0.4–0.9 R.
- Arc 16.94° long. Strand FWHM 0.28°. Southern spine at 12.10 ± 0.12° from the Knies two-SNR centre, against an X-ray radius of 11.0 ± 1.1°.
- About 1° outside the eRASS1 ring peak. The local X-ray step is ~2σ, south only.

**Monogem (Knies+24 Table 4)**
- kT 0.11–0.17 keV; n₀ f^−½ = 3.1–8.6 × 10⁻³ cm⁻³; R = 58–63 pc at 300 pc.
- Age 1.1–1.2 × 10⁵ yr; E₀ ≈ (0.1–0.4) × 10⁵¹ erg.

---

## 2. Physical-consistency checks (derived here unless cited)

### 2.1 Does Monogem drive 100–150 km/s shocks into 0.01–0.04 cm⁻³ gas?

**Blast-wave speed and energy (Sedov).** With R = 54–58 pc (11.0° at 282–300 pc), t = 1.1–1.2 × 10⁵ yr and n_amb = 3–9 × 10⁻³:
- v_bw = 0.4 R/t = **176–223 km/s**.
- Sedov E = 3 × 10⁴⁹ – 2.3 × 10⁵⁰ erg. That is the same order as Knies' (0.1–0.4) × 10⁵¹; definitions differ by ~2×.
- Post-shock pressure P/k = 1.2–5.4 × 10⁴ K cm⁻³.

**The ambient blast wave is non-radiative.** Using the Cioffi+88 shell-formation time, t_PDS = 3.3–5.7 × 10⁵ yr, which is much longer than the age. So Monogem's ambient blast wave cannot emit [O III] on its own. Radiative [O III] can arise only where the shock enters gas a few times denser than ambient.

**The transmitted cloud shock lands in the right range.** Its speed is v_c ≈ v_bw √(β/χ), where χ is the density contrast and β = 1–2.5 is the reflected-shock enhancement. This gives **100–150 km/s for χ ≈ 2–8 (β = 1–2.5)**. That is exactly the contrast the shock fit requires: n₀/n_amb = 2–8 (n₀ = 0.012–0.042).

**The same picture explains Fesen's morphology.** Fesen+26 find [O III] filaments around much of the Monogem X-ray boundary (E, S, NW), and Hα+[S II] filaments where the shock meets denser clouds. In this picture:
- χ ≲ 2 gives an invisible, non-radiative shock;
- χ ≈ 2–8 gives incomplete [O III] shocks (the Silk);
- χ ≫ 10 gives slow, complete, Balmer-bright shocks (the Hα/[S II] patch beside PGCC G199.09+19.95, and the NE/"DNA" Hα filaments).

**Energy budget.** At 282 pc a sheet 83 pc long with a 25 pc limb path has an area of ~2 × 10³ pc², about 5 % of Monogem's surface.
- It processes ½ρ₀v³ × area × t* = **0.5–6 × 10⁴⁸ erg**, i.e. 0.1–6 % of E₀.
- That is consistent with its share of the area.
- The same sheet at 1.6 kpc needs 1 × 10⁵⁰ – 1 × 10⁵¹ erg, which is superbubble-class (H4).

### 2.2 Pressure balance: the "Monogem-pressure coincidence" is not distance-specific

`verify_shockmodel` presented the following as "a non-trivial coincidence and the main support" for H1: width plus brightness give ρ₀v² ≈ Monogem pressure for b = 7–14.

**Why it carries no distance information.** For an edge-on truncated shock:
- I_pk = G·η·n₀v_s/4π (photon units), where G ≈ 20 is the limb-brightening gain and η = 0.35–0.53 is the number of 5007 photons per swept-up H atom.
- So ρ₀v_s² = 1.4 m_H · 4π I_pk · v_s/(G η). This depends only on the observed brightness, the velocity (from the line ratios), G and η.
- G depends on the ratio R_c/Δ, which is distance-independent when both scale with d.

**Scan of the 238 spectrally acceptable runs** (`w3_judge_phys_checks.json` → `closureW_distance_scan`). I rescaled each run to distance d (n₀ ∝ 1/d, I ∝ 1/d at fixed b) and required:
- I_pk = 1.2–2.5 R;
- FWHM ≤ 38 km/s;
- FUV 40–160 CU per R.

| d (pc) | Matching runs | b (µG cm^3/2) | P_ram/k (K cm⁻³) | t* (yr) |
|---|---|---|---|---|
| 160 | 20 | 3–7 | 2.1–5.8 × 10⁴ | 0.4–0.8 × 10⁵ |
| 282 | 25 | 5–10 | 2.8–6.2 × 10⁴ | 0.5–1.9 × 10⁵ |
| 500 | 23 | 7–14 | 2.5–6.5 × 10⁴ | 0.7–2.3 × 10⁵ |
| 1000 | 18 | 10–20 | 2.7–6.9 × 10⁴ | 1.0–3.6 × 10⁵ |
| 1600 | 12 | 14–30 (grid edge) | 2.9–9.3 × 10⁴ | 1.1–4.2 × 10⁵ |

**Consequences**
1. The closure constrains **b(d)**, not d. At 1.6 kpc it needs B₀ = b√n₀ ≈ 2–4.5 µG, which is not implausible.
2. **The useful physical statement is that the driver must supply P/k ≈ 3–9 × 10⁴ K cm⁻³** (for G ≈ 20).
   - That is typical of an SNR interior: Monogem gives 1.5 × 10⁴ – 1.6 × 10⁵ depending on definition (canon C10).
   - It is ≥3× above the Local Bubble hot-gas pressure (~10⁴; literature order, not re-measured).
   - It is ~10× above the quiescent thermal pressure expected at z ≈ 0.3–0.55 kpc (literature order).
   - So H4 needs an *active*, high-pressure superbubble interior at high z, and a Local-Bubble-driven shock (part of H2) is marginal.
3. The definitional ×10 spread in "Monogem pressure" (C10) means the phrase "pressure balance with Monogem" is currently unfalsifiable. It is consistent, not confirming.

### 2.3 Timescales

**Truncation is natural.** t* = 0.4–1.0 × 10⁵ yr, less than t_cool = 0.7–1.4 × 10⁵ yr, which is at most the Monogem age of 1.1–1.2 × 10⁵ yr (and the PSR B0656+14 τ_c of 111 kyr). Truncation therefore needs no fine-tuning, provided the shock entered the denser sheet partway through the remnant's life.

**Uniformity in the south needs nearly simultaneous contact.** Along the southern ~30 pc, DESI [O II]/[O III] = 1.04–1.25 in four bins, so N* is uniform to ≲ 0.1 dex. That requires the shock to have reached the sheet at nearly the same time along its length (ΔN*/N* ≈ Δt/t*).
- For a straight sheet tangent to the blast sphere, the contact-time spread over the southern section is ~1.6 pc / 200 km/s ≈ 8 × 10³ yr (the 0.32° sagitta at 282 pc).
- That gives ΔN*/N* ≈ 0.1–0.2, consistent with the data.

**The X-ray offsets are consistent with the layer thickness.**
- In H1 the hot gas abuts the cooled side of the layer, so the X-ray edge should sit within s* = 1.0–1.7 pc (**0.20–0.35°** at 282 pc) west of the [O III] ridge.
- The best eRASS1 step location in the south is −0.2° (2σ only).
- The ring peak 1.1° (5.4 pc) inside is ordinary thick-shell structure. It also matches the lower end of the pre-shock penetration depth v_s t* = 4–16 pc.

**A fossil is excluded only by needing a recent hot start.** No variability above ×2 is seen over 2002–2022. Once the doublet density falls (`w3_oii`), a fossil's e-fold time at n ≈ 0.5–1 is 4–7 × 10³ yr, so the 20-yr constancy no longer excludes it. It is still excluded because it needs a start temperature ≥ 4 × 10⁴ K reached within ~10⁴ yr, i.e. a recent shock, which is H1/H2 under another name.

### 2.4 Geometry: length, width, curvature, and a refinement of H1 (H1b)

**Size at 282 pc**
- Length 83 pc (75–94 pc at 252–318 pc). It spans PA 34°–104° about the Knies centre, i.e. 19 % of the Monogem circumference.
- FWHM 1.38 pc. This equals the cooling-layer thickness s* = 1.0–1.7 pc in the closure, so the width is physical, not instrumental.

**The arc is not a spherical Monogem limb.** It is one R ≈ 21° arc (103 pc at 282 pc). It is tangent to the 12.1° Monogem circle near Dec 14–16 and bends 1–4° outside it in the north (`w3_northsouth`).
- A straight sheet of half-length 41.7 pc tangent to a 59.6 pc sphere has its ends at 72.7 pc, i.e. **14.8°**. The observed ends are at 12.5° and 15.8°.
- So geometrically the Silk behaves like a **flat-ish sheet (R_curv ~100 pc) touched tangentially by the Monogem blast wave**.

**H1b, a variant of H1 I propose.** The Monogem blast wave strikes a pre-existing, flatter sheet of χ ≈ 2–8 gas, and the transmitted shock takes the sheet's shape, not the blast sphere's.
- The bright southern 9.8° (48 pc) is then the contact chord of a blast wave of radius R_bw = √(59.6² + 24.1²) ≈ **64 pc (13.1°)**. That is ~1–2° beyond the X-ray-peak radius, which is plausible because X-ray radii trace emission peaks.
- Contact is most oblique, and the driving pressure lowest, toward the north. That is consistent with the slower 70–100 km/s northern shock at P_N/P_S ≈ 0.45–0.7.
- At fixed T, X-ray surface brightness scales as P², so a 0.2–0.5× fainter northern rim explains why no northern X-ray step is seen.
- H1b keeps the 250–320 pc distance; the tangency point fixes it. Its weakness is that the faint northern [O III] reaching 15.8° (79 pc) needs the blast wave to extend that far in the north, i.e. an aspherical remnant expanding toward higher |b| (lower density).

**A Local Bubble (LB) wall shock at 160 pc (part of H2) is geometrically disfavoured.**
- From inside a bubble the wall is seen roughly face-on. The O'Neill+24 wall sits at a smooth 158–167 pc along the spine (Dec ≤ 22), so no edge-on fold is indicated.
- A face-on sheet has G ≈ 1, so it would need ~20× the emissivity, i.e. P/k ~10⁶.
- It would also appear as a diffuse patch, not a 17° × 0.28° filament.

**Chance alignment** (`chance_alignment_south_segment`). I placed the 23-point southern segment randomly (random spin; centroid uniform within θ_max of the Knies centre).
- The probability that it lands within ±0.5° of the observed mean radius *and* is as tangential (rms ≤ 0.25°) is 0.27–0.53 % for θ_max = 30–20°.
- A Monogem-limb model gives 0.31 per degree of mean radius in the same window.
- **The Bayes factor, limb against random placement, is 60–120.**
- Look-elsewhere: about six candidate shells were tested (1-SNR, 2-SNR, PSR, G205.6+12.4, Geminga, Origem). Correcting for that reduces it to an **effective ~10–20**.
- The south segment's straightness relative to the 12° circle (χ² 375 vs 114, `w3_northsouth`) costs H1a, but not H1b.

### 2.5 Where must the driver be?

**The driver is on the X-ray-bright side.**
- In every shock or interface picture, the driving high-pressure gas lies on the concave (west) side. That side is X-ray-brighter (W/E 1.26).
- Hot gas at P/k ~5 × 10⁴ and 1–3 MK on the convex (east) side would have an emission measure comparable to Monogem's. It would appear as X-ray emission east of the Silk, which is the faint side.
- So the driver is X-ray-bright gas to the west, which is Monogem or something projected onto it.
- The alternative is gas too distant to see in the 0.2–0.6 keV band: N_H ≈ 2–3.5 × 10²⁰ to 1–2 kpc gives τ(0.25 keV) ≈ 2–3, which is H4.

**Consequences for H2**
- **G205.6+12.4** (kT 0.28–0.45 keV, ~5 × 10⁴ yr) contains the southern spine *inside* its projected disk: 6.5° from its centre against R = 8.82°. The southern spine therefore cannot be its tangent limb, which needs G ≈ 20.
- **A Geminga-progenitor remnant** is energetically possible: a 10⁵¹-erg Sedov remnant in n = 0.005 reaches ~150 pc at 3.4 × 10⁵ yr with v_bw ≈ 170 km/s (derived). But the spine is not concentric with the back-traced birth locus (rms 1.1°), and no ~60°-diameter hot shell is seen.
- These are the reasons H2 sits well below H1. It is not zero, because an unknown remnant hidden inside the Monogem X-ray glow cannot be excluded.

### 2.6 H3, the paper's model, fails by large factors

Inputs: n₀ = 1, R = 30 pc, E = 2 × 10⁵¹ erg, PDS shell, FWHM 0.28° = 0.78 pc at 160 pc.

| Observable | H3 prediction | Measured / limit | Factor |
|---|---|---|---|
| Interior pressure (thermal fraction 0.3–0.5) | P/k ≈ 0.9–1.5 × 10⁶ K cm⁻³, 6–100× Monogem | — | Implies an X-ray-bright, 21°-diameter shell (EM ~10²–10⁴× Monogem) that is not seen |
| Swept-up shell column at the limb | N_H ≈ 5.4 × 10²⁰ cm⁻² | HI < 1–3 × 10¹⁹; dust N_H < 2.3 × 10¹⁹ | ≥ 20× |
| Hα if the shell is ionised | EM ≈ 2200 pc cm⁻⁶, i.e. Hα ≈ 800 R | < 1 R | ~10³× |
| [O III]/Hβ of a complete (PDS) radiative shell | ≤ 5 | > 11 | — |

- With n₀ = 0.03 the shell column (1.6 × 10¹⁹) and Hα (0.7 R) just reach the limits, but the paper's own formula then gives E ~4 × 10⁴⁹ erg. **H3 then collapses into H1 or H2.**
- The paper's "R ≈ 30 pc at 160 pc" corresponds to 10.74° (`verify_paper`). At 282 pc that is 53 pc, the Monogem radius.

### 2.7 Kinematics against distance (H4)

Expected V_LSR under flat rotation (Θ₀ = 236 km/s, R₀ = 8.18 kpc), at the south end (l 208°) and north end (l 191°) of the spine:

| d (kpc) | South (km/s) | North (km/s) |
|---|---|---|
| 0.282 | +3.0 | +1.3 |
| 1.0 | +10.0 | +4.3 |
| 1.6 | +15.2 | +6.6 |
| 2.0 | +18.4 | +7.9 |

- **Observed:** south −2.8…+2, north +4.5…+10.4.
- At 1–2 kpc the south is 10–18 km/s too low.
- Rotation predicts a north-minus-south gradient of **−6 to −10 km/s at 1–2 kpc**, against **+7 to +11 km/s observed**. The sign is opposite at every distance, so the gradient is local: pre-shock gas motion or sheet tilt.
- With σ_pec ≈ 7 km/s per position, each discrepancy is ~1.5–2σ against H4. They are not decisive, but they point the same way.

### 2.8 Analogue objects

**Cygnus Loop truncated shocks (Raymond+88 SE).** v = 140 km/s, n₀ = 4.4 cm⁻³, log N = 18.2.
- Face-on [O III] ≈ 38 R, i.e. η ≈ 0.61 photons per swept-up H atom (derived from the face-on Hβ and ratio quoted in `verify_shockmodel`). The Silk models give η = 0.35–0.53.
- Line ratios depend only on (v_s, N*, b, pre-shock ionisation X). Face-on flux scales as n₀, and layer thickness as N*/n₀.
- **The Silk is the same physics at ~100–370× lower pre-shock density:**
  - face-on brightness 0.06–0.2 R instead of tens of R;
  - emitting layer ~1–1.7 pc instead of ~0.01 pc;
  - hence a 0.28°-wide, ~2 R ridge after G ≈ 20.
- R88 SE also matches the optical diagnostics: [O III]/Hβ 12.9, [O II]/[O III] 1.34, 4363/5007 0.071.
- This homology is a strong consistency argument for the truncated-shock family. It cannot separate H1, H2 and H4.

**C IV suppression.**
- G65.3+5.7 needed ≥ 20× C IV suppression (Kim+2010, via `verify_xrayfuv`). The Silk needs ~4–6×.
- Resonance scattering in an edge-on sheet with τ₁₅₄₈(edge) ≈ 10–40 is the natural cause, with precedent.

**[O III]-dominated, radio-faint large remnants**
- **G206.6+6.1:** 3.6°, ~440 pc HI-kinematic distance, [O III]-dominated with diffuse Hα; faint FAST radio shell with Σ₁GHz ≈ (2.5–3.8) × 10⁻²³ (Fesen+26, on disk).
- **G239.9+7.0:** [O III]-dominated; 110–130 km/s inferred (Fesen+26).
- **Nereides (G107.5−5.1):** ~3° [O III]-filamentary shell; Hα mostly diffuse; no radio or ROSAT X-ray detection (Fesen+2024 ApJS 272, 36, *web, not re-verified*).
- These show that Balmer-faint, radio-quiet, X-ray-faint [O III] shells are an established class of old-SNR shock. The Silk's radio limit (Σ₁GHz < 1.3 × 10⁻²²) and its lack of an X-ray filament are typical, not anomalous.

**Antlia SNR** (Fesen+2021, arXiv:2102.12599; *web, not re-verified*).
- About 20° × 26°, ~250 pc, b = +18°. Its ~20°-long FUV filaments along the NE limb correlate with Hα. Shocks are 70 to > 100 km/s.
- It is the closest analogue in angular size, distance and latitude. The contrast is instructive: Antlia's limb is Balmer-bright (complete or slower shocks), whereas the Silk is Balmer-dark.
- The difference fits the χ-dependence in §2.1: an older remnant, or denser gas, completes its shocks.

**Within Monogem itself.** Fesen+26 put [O III] where v ≳ 100 km/s and Hα+[S II] where the shock is slower in denser clouds. The Silk and its Balmer-bright patch (shock–cloud, on the same Monogem-radius circle: 11.94–11.97° against 12.10 ± 0.12°) reproduce that dichotomy on a small scale.

---

## 3. Scorecard by hypothesis

### H1: Monogem-driven incomplete radiative shock at ~250–320 pc. P = 0.64, favoured

**For**
1. **Spectrum:**
   - [O III]/Hβ > 11;
   - T(O III) and T(O II) both ≳ 3 × 10⁴ K (4363 and the corrected 7325, two independent marginal thermometers that now agree);
   - He II < 0.026; Balmer-dark; FUV 68–92 CU per R with C IV scattering.
   All of these are reproduced by column-truncated 100–150 km/s shocks and by nothing else tested.
2. **Dynamics** (§2.1). Monogem's Sedov blast wave (~175–225 km/s, non-radiative, t_PDS ≫ age) naturally drives 100–150 km/s shocks into χ ≈ 2–8 (β = 1–2.5) gas, the contrast the fit requires.
3. **Timescales** (§2.3). t* < t_cool ≲ age, so truncation is natural. The southern uniformity matches near-simultaneous contact.
4. **Geometry and context.**
   - The southern spine sits at the Monogem radius and is tangential (effective Bayes factor ~10–20 against an unrelated placement).
   - It lies on the outer flank of the X-ray shell, with the driver on the X-ray-bright side (§2.5).
   - It belongs to Fesen's population of Monogem-rim [O III] filaments.
   - Its Balmer-bright shock–cloud patch lies on the same circle.
5. **Energetics.** The sheet processes 0.1–6 % of E₀, matching its ~5 % share of the surface area.
6. **Kinematics.** The centroid is at the local pre-shock gas velocity, as a tangent limb requires. The northern redshift tracks the redder northern HI.

**Against**
1. No distance is measured. 252–318 pc is conditional on the morphology.
2. The south is tangent to, not concentric with, the 12° circle, and the north lies 1–4° (5–20 pc) outside the Monogem sphere with no X-ray step. H1b and pressure scaling soften this, but at the cost of an aspherical blast wave.
3. The local X-ray step is only ~2σ; W/E is generic to the shell.
4. The FUV agreement needs ~4–6× C IV suppression, a geometry-dependent factor not measured directly.
5. The pressure match is not distance-specific (§2.2).

**Decisive tests:** T1, T2, T3, T5, T7 and T8 in §4.

### H2: incomplete shock driven by another remnant or bubble. P = 0.26, viable

**For**
- The same spectral and analogue physics as H1 applies.
- The distance is genuinely unmeasured.
- The arc as a whole is not a Monogem circle.
- The region holds several remnants: G205.6+12.4 at 300–500 pc, G206.6+6.1 at ~440 pc, and others.

**Against**
- Each named variant fails a specific test:
  - G205.6+12.4: the spine is inside its disk.
  - Geminga: concentricity is not physical, and no 60° hot shell is seen.
  - Local Bubble wall: face-on geometry, P/k too low.
- An unnamed driver must hide inside Monogem's X-ray glow.
- It also loses the radius and tangency coincidence (effective Bayes factor ~10–20).

**Decisive tests**
- An absorption or proper-motion distance outside 225–359 pc (T1, T2).
- Deep X-ray showing a gap of more than s* between the Monogem hot gas and the Silk, or a distinct hot component next to it (T5).
- A distance for the PGCC clump at the Hα patch inconsistent with ~280 pc (T7).

### H3: the paper's model. P = 0.005, ruled out

**For:** a shock origin and the "possible SNR" are right in kind.

**Against** (§2.6):
- The distance is a unit error.
- E needs n₀ = 1 taken from WIM [S II] (whose own Table 2 gives 225–631 cm⁻³, not 64–150).
- A PDS shell at n₀ = 1 predicts an HI or dust limb ≥ 20× above the limits, or ~800 R of Hα if ionised.
- It predicts an interior at P/k ~10⁶ that would be one of the brightest soft-X-ray objects in the sky.
- It predicts complete-shock ratios ([O III]/Hβ ≤ 5) against > 11 measured.
- R ≈ 30 pc is asserted, not fitted.

**Decisive tests:** none needed. Any version at 160 pc with physical parameters is an H2 sub-case, testable by T1 and T2.

### H4: superbubble or chimney wall at 1–2 kpc. P = 0.04, disfavoured

**For**
- The paper's kinematic distance, but with the refuted +25 km/s velocity: the revised V_LSR ≈ +1 gives 0–0.2 kpc.
- The arc's curvature centre lies in the plane at (l, b) = (194.3, +0.4), 0.8° from the Origem Loop centre and ~6° from Gem OB1, which is suggestive of a structure rising from the disk.
- Distant hot gas would be X-ray-absorbed, so it is hard to exclude via X-rays.
- The shock closure can be satisfied at 1.6 kpc with b ≈ 14–30 (§2.2).

**Against**
- Length 470 pc at z ≈ 0.55 kpc.
- It needs an active driver at P/k 3–9 × 10⁴ at high z, ~10× the quiescent value.
- It needs 10⁵⁰–10⁵¹ erg processed through the sheet.
- The velocity is too low by 10–18 km/s, and the N–S gradient has the wrong sign (~1.5–2σ each, §2.7).
- No HI supershell is seen at the [O III] velocity.
- It loses the Monogem coincidence.

**Decisive tests**
- A proper motion ≤ 0.2″ per decade against ~1″ for H1 (T2).
- C IV/Si IV absorption absent toward every star within 1 kpc (T1).

### H5: DIG or photoionised bubble wall. P = 0.01, ruled out

**For**
- It is right for the *red* lines: [S II]/Hα, [N II]/Hα and the paper's T_e, V_t and O/H all describe the WIM, identical in P4.
- At the 2σ floor ([O III]/Hβ ≈ 11), a hard (150 kK or 0.14 keV) spectrum at log U ≈ −3.3 reaches [O III]/Hβ ≈ 9–14.

**Against**
- There is no ionising source:
  - log U available ≤ −6.1 from Monogem/Local Bubble EUV, against −3.3 needed;
  - WD J0736+1444 is short by 27× (coplanar, n_rms) to 10⁴× (at 282 pc);
  - the ISRF/DIG gives [O III]/Hβ ≤ 1.3 at any U.
- Photoionisation equilibrium gives 4363/5007 ≤ 0.011–0.034 against 0.06–0.09, and 7325/3727 ≈ 0.02–0.04 against 0.060.
- Hard spectra overproduce He II (0.06 against < 0.026).
- The FUV is 1.3–1.9× low (49–53 CU per R).

**Decisive test:** T4. Deep Hβ and 4363 (photoionisation predicts [O III]/Hβ ≤ 14 and 4363/5007 ≤ 0.012), plus line width (photoionisation predicts FWHM ≤ 12 km/s).

### H6: other non-shock mechanisms. P = 0.025 summed, disfavoured

**Conduction front.** This is the least-dead member.
- The wave-2 FUV exclusion is weakened, because the same C IV scattering invoked for H1 applies here.
- The brightness exclusion is pressure-dependent. A static radiative front's flux scales as P (derived): at P/k = 1.5 × 10⁴ / 5 × 10⁴ / 1.6 × 10⁵ one needs 31 / 9 / 3 superposed edge-on fronts for 1.6 R.
- The front predicts [O II]/[O III] ≥ 2.0 against 1.1 in the south. It predicts 4363/5007 = 0.11–0.14 against 0.06–0.09 (1–2σ).
- Perpendicular conduction is suppressed in a field tangential to the boundary.
- **Robust discriminator:** the order of ionisation stages across the filament (T3).

**Other members**
- Turbulent mixing layer: needs a 377 km/s inflow, which is a shock.
- NEI cooling of hot gas: predicts an X-ray filament ≥ 10³× Monogem; none is seen. This is robust.
- WD or Monogem-EUV photoionisation: as H5.
- Fossil: needs a hot start ≤ 10⁴ yr ago, which is a shock (§2.3).

### Added hypotheses

**H7: runaway-star or pulsar bow shock or wake. P = 0.005.**
- The context wake test found no hot star whose 0–1 Myr track follows the spine; the best cases are 20–31 kK WDs with negligible Q(H).
- Pulsar bow shocks are arcsecond-to-arcminute, Balmer-dominated structures.
- Neither B0656+14 nor Geminga moves along the arc.

**H8: IVC/HVC–halo interaction. P = 0.003.**
- The GALFA compact clouds on the line of sight are at +52…+149 km/s. The [O III] sits at the local-gas velocity (V_LSR ≈ +1), not the cloud velocity.

**H9: unanticipated. P = 0.012.**
- For example, MHD-wave or cosmic-ray heating of a magnetised sheet, or a shock whose driver is something not yet catalogued.
- It is kept small because the spectrum is so specifically that of a truncated cooling flow.

---

## 4. Decisive future observations (concrete; predictions per hypothesis)

### T1. UV absorption-line distance bracket

- **Instrument:** HST/COS G130M + G160M (R ≈ 16–20 k).
- **Lines:** C IV 1548/1550, Si IV 1393/1402, N V 1238/1242, and Si III 1206 / C II 1334 for the low ions.
- **Targets:** 4–6 UV-bright stars (hot WDs, sdB/O, late-B) within 0.3° of the southern ridge, with Gaia distances bracketing 160, 250, 320 and 600 pc. The target list can be built now from Gaia DR3 × GALEX; this has not been done.
- **Required depth:** log N(C IV) ≈ 12.7 at 3σ, i.e. S/N ≳ 15 per resolution element.

**Predictions**
- **H1:** for stars beyond ~320 pc and sightlines crossing the layer (≤ 0.3° east of the [O III] ridge), N(C IV) ≈ 10¹³·⁵–10¹⁴·⁵ cm⁻² near V_LSR 0 ± 20 km/s (edge-on through the C IV zone: face-on 1–3.5 × 10¹³ × G_CIV 19–35). Face-on parts of the sheet would show components at ±50–130 km/s. Nothing for stars < 250 pc.
- **H3 / Local Bubble H2:** absorption already in stars at 170–250 pc.
- **H4:** nothing for any star < 1 kpc.
- **H5:** no C IV; low-ion absorption only.
- **Conduction front:** C IV ≈ 10¹²·⁵ per front, with O VI ≫ C IV.
- **Existing constraint:** IUE HD 60848 (1.9 kpc, 2.26° W) shows no |v| = 70–200 km/s C IV, at log N < 12.9. This argues against a complete radiative shell over Monogem's face-on cap, consistent with radiative emission only in χ ≈ 2–8 sheets.

### T2. Two-epoch proper motion of the east edge

- **Instrument:** deep [O III] imaging at ≤ 1.5″ FWHM with a 3 nm filter; 2 R is about 25.6 AB mag arcsec⁻² (60 Å). Use a 0.5–1 m telescope with ~30–50 h per field on 2–3 southern fields; epoch 1 in 2026–27, epoch 2 in ≥ 2036.
- **Measurement:** register stars to ≤ 0.05″ and measure arcsec-scale edge structure. The [O III] onset behind the front is ~0.01 pc ≈ 7″ at 282 pc if the sheet is seen exactly edge-on.

**Predictions: μ = v_s/(4.74 d)**

| Hypothesis | Distance | Proper motion |
|---|---|---|
| H1 | 282 pc | 75–112 mas/yr = **0.75–1.1″ per decade**, radially away from (104.15, +15.30) |
| H3 | 160 pc | 1.3–2.0″ per decade |
| H4 | 1.6 kpc | 0.13–0.2″ per decade |
| H5 / photoionisation | — | ≲ 0.1″ per decade (bulk ≤ 10 km/s) |

Combined with v_s from T3/T4 line profiles, this is an expansion-parallax distance good to ±25 %. The motion direction separates Monogem (PA of motion ≈ 60–100°) from the curvature centre (motion due east).

### T3. Ionisation stratification across the filament

This needs mostly existing data, and has not been done.
- **Measurement:** cross-spine centroid offsets of FUV (GALEX AIS 5″ data or stacked skybg), [O III] (NSNS or DESI), [O II] (a DESI cross-spine regression with offset templates) and the X-ray edge.
- **Shock moving east (H1, H2, H4):** from east to west, C IV/FUV → [O III] → [O II] → hot gas. FUV leads [O III] by ≈ 0.1° (`verify_shockmodel`), and [O II] trails by ≈ 0.05–0.1°.
- **Interface heated from the hot side (conduction, mixing layer, NEI):** the order reverses. From west to east, hot gas → C IV → [O III] → [O II].
- **Photoionisation by an external source:** stratified toward the source, with no FUV ridge.
- **Value:** this is the single cheapest test that can kill the whole H6 family and confirm the outward-moving shock.
- **Tension to resolve:** `verify_xrayfuv` reports the FUV peak "matching" the [O III] position; a 0.1° offset was never measured.

### T4. Deep optical spectroscopy of the southern ridge

- **Instrument:** an R ≈ 4000 large-étendue IFU or fibre array (HET/VIRUS-class, or an LVM-class instrument pointed off-plane), ≥ 10 h on-ridge, with interleaved off-ridge sky and sky subtraction at the ≤ 0.1 % level.
- **Required depth:** 1σ ≈ 0.005–0.01 R at Hβ and 4363; S/N ≳ 30 on the [O II] doublet.

**Predictions**

| Quantity | H1 / H2 / H4 (truncated shock) | Photoionised (H5 / H6) | Complete shock |
|---|---|---|---|
| [O III]/Hβ | 15–60 (model family 10–63) | ≤ 10–14 | ≤ 5 |
| 4363/5007 (Hg 4358 resolved) | 0.07–0.09 | ≤ 0.012–0.034 | — |
| He II/[O III] | ≤ 0.005 | 0.025–0.06 | — |
| [O II]3729/3726 | **1.45–1.47** (n_e ≲ 10) | 1.46 | — |

- σ(ratio) ≈ 0.03–0.04 separates 1.31 from 1.46 at ≥ 4σ.
- A firm 1.30 would revive the high-pressure problem for every pressure-balanced model.

### T5. Soft X-ray structure at the ridge

- **Instrument:** an eROSITA multi-survey stack (if released) or a pointed XMM-Newton strip across the southern spine.
- **H1:** the 0.2–0.6 keV emission ends within s* ≈ 0.2–0.35° west of the [O III] ridge. Thermal P/k at the inner flank ≈ (1/β) ρ₀v_s² ≈ (1–9) × 10⁴. Nothing east of the ridge. No X-ray filament.
- **H2 with a separate driver:** a distinct plasma component (different kT, abundances or absorption) adjacent to the ridge.
- **Mixing layer or NEI cooling:** an X-ray filament on the ridge (already excluded at ≥ 10³×).

### T6. High-resolution line profiles

- **Instrument:** R ≥ 30 000 (≤ 10 km/s) on [O III]5007, Hα and [O II]3726/29, with large étendue: a spatial-heterodyne spectrometer or a 1 m-class Fabry–Pérot integrating ≥ 1 arcmin².
- **Truncated shock:** flat-topped or double-horned profile, |v_los| ≤ 12–22 km/s, FWHM 23–38 km/s, centroid = pre-shock gas (local HI or WIM). This also resolves canon conflict C7 (DESI ≤ 33 against paper ~50).
- **Photoionisation or conduction:** Gaussian with FWHM ≈ 5–12 km/s (thermal O at 10⁴–5 × 10⁴ K) plus modest turbulence.
- **North:** a symmetric profile about the +8.5…+11 km/s HI component (pre-shock-gas explanation) against an asymmetric red tail (tilted sheet or one-sided chord; `w3_northsouth`).

### T7. Distance of the shock–cloud patch

- **Method:** Gaia DR3 XP or Gaia+2MASS extinction against distance toward stars within 0.3° of PGCC G199.09+19.95 and the Hα/[S II] patch. The clump should give a resolvable A_V step; the Edenhofer map at 0.46° pixels shows no significant feature, peak S/N ≈ 1.4 at Dec 21.
- **Predictions:** H1, a step at 250–320 pc. A step at ~160 pc would mean the clump is LB-wall foreground or that H2/LB is right.
- **Status:** cheap; not done.

### T8. Family test within Monogem

- **Method:** apply the DESI/SDSS sky-fibre regression already built to Fesen's *other* Monogem [O III] limb filaments (S and NW groups).
- **H1:** the same truncated-shock spectra ([O III]/Hβ ≳ 10, [O II]/[O III] ≈ 1–2, 4363 high), each at its local pre-shock velocity, with brightness scaling with local χ.
- **H2:** no such family requirement.
- **Status:** uses existing data; not done.

### Secondary tests

**FUV spectroscopy** (e.g. a UVEX-class long slit, 1150–2650 Å).
- H1 predicts:
  - a C IV 1548/1550 ratio near 1–1.3 (optically thick);
  - C IV/O III]1663 suppressed ~5× relative to optically thin;
  - O III]/[O III] ≈ 0.5.
- The northern FUV (Dec > 22, never covered by GALEX) should be 0.25–0.63× the southern FUV/[O III].

**Rotation-measure grid.**
- A shock-compressed field (B ≈ 5–30 µG in the layer, n_e ≈ 0.1–0.2, path ≈ 25 pc) gives an RM step of ~10–50 rad m⁻² confined to ≤ 0.3–0.5° of the ridge.
- This needs ≥ 5 polarised sources per deg² (VLASS/LoTSS/next generation).
- The NVSS step is not filament-specific (p = 0.18).

---

## 5. Weakest links in the team's analysis

1. **The distance, and so the H1/H2 split, is a single morphological coincidence.**
   - It rests on a 5.5°, nearly straight segment whose mean radius matches an X-ray radius that is itself uncertain by ±1.1° (Knies' assumed 10 %).
   - Fesen+26's assignment uses the same morphology and is not independent.
   - The synthesis JSON still cites the superseded W/E 2.4σ edge (canon C12).
   - My chance-alignment estimate gives an effective Bayes factor of only ~10–20 after look-elsewhere.
   - Quote 252–318 pc as *conditional*, never as a measurement.
2. **The "pressure balance with Monogem" support is overstated.**
   - P_ram follows from the observed brightness and v_s for any distance (§2.2).
   - "Monogem pressure" spans ×10 by definition (C10).
   - The shock grid was fitted with withdrawn inputs ([O III]/Hβ = 40, 7325/[O III] = 0.034 instead of 0.07–0.08) and was never run as a distance scan. A re-fit with the wave-3 inputs, b > 30 and d free is needed.
3. **The exclusion of non-shock models leans partly on marginal or revised lines.**
   - 4363 is a 2–3σ detection 5 Å from Hg 4358.
   - 7325 was mis-normalised in both earlier waves and carries a ±0.015 OH systematic.
   - The doublet density used to exclude every pressure-balanced alternative was an artefact.
   - The FUV argument against conduction fronts is weakened by the same C IV scattering invoked for the shock.
   - **What remains robust:** the ionising-photon (U) budget, the X-ray non-detection against NEI cooling and mixing layers, and the [O III]/Hβ > 11 floor. These carry H5 and H6 to ≤ 3 %, not the line-ratio thermometers alone.
4. **The C IV suppression factor (~4–6×) is a model-dependent free parameter.**
   - It comes from an idealised slab Monte Carlo with no re-absorption, and the FUV/[O III] ratio inherits the NSNS Rayleigh scale (×1.0–1.5).
   - GALEX also correlates with the NSNS Hα map at z ≈ 7–11, a confusion term in the template ratios.
5. **The north–south kinematics and ratios are estimator-dependent.**
   - ΔV shifts by ~4 km/s with pixel weighting.
   - DESI and SDSS disagree at Dec 17.6–21 (2–3σ), exactly where the W strand and Hα patch lie.
   - The "pre-shock HI" identification is assumed.
   - The line width, a first-order shock-versus-photoionisation discriminator, was never re-verified in wave 2 (C7).

## 6. Important things not yet done

- **Cross-spine ionisation stratification (T3)** with existing GALEX, DESI and eRASS1 data. This is the cheapest decisive test of shock against interface.
- **Re-fit of the shock grid** with wave-3 inputs: DESI south ratios, 7325/[O III] ≈ 0.07, doublet 1.31 ± 0.07, separate north/south fits, an explicit d and b scan, and an R88-calibrated Hβ.
- **Spectra of other Monogem [O III] filaments** from existing DESI/SDSS sky fibres (T8): the family test for H1.
- **Extinction distance of the PGCC G199.09+19.95 clump and the Hα patch** from Gaia (T7).
- **A UV-bright background-star list** within 0.5° of the southern spine at 150–700 pc for a COS proposal (T1).
- **GALEX-based check of the predicted 0.1° FUV lead** over [O III] (part of T3), and of FUV in the north to Dec 22.
- **Geometry for H1b:** fit a "plane or low-curvature sheet tangent to a sphere" model to the full spine, as an alternative to a single circle and to a Monogem circle.

## 7. Where I differ from the modellers

- **The pressure coincidence.** `verify_shockmodel` calls it the main support for the Monogem picture. I find it distance-independent (§2.2): it supports "an SNR-interior driver," not "Monogem at 280 pc."
- **altmodel's high-pressure conclusions.** Its "shock into n₀ ~ 1 cm⁻³" and "two-temperature cooling-flow fingerprint" were both driven by the doublet and the mis-normalised 7325. Neither survives `w3_oii`. The low-density truncated shock of `verify_shockmodel` is now internally consistent, with no remaining ×1.8 or doublet tension.
- **H1 as a spherical limb.** I would not state H1 as "the Monogem limb" (a spherical shell). The arc is tangent to the Monogem sphere and departs from it in the north. H1b (Monogem blast wave on a flatter external sheet) fits the geometry, the uniform southern truncation and the slower northern shock better.
- **Non-shock exclusions.** I rate the H5 and H6 exclusions as secure because of the energy and photon budgets and the X-ray limits, not because of 4363 or 7325.

---

## 8. Files (tag `judge_phys`)

**Scripts** (run with `~/.venvs/astro313/bin/python -I`)
- `scripts/w3_judge_phys_checks.py`: Sedov dynamics, cloud-shock speeds, the closure-W distance scan of the 238 `verify_shockmodel` runs, the chance-alignment Monte Carlo (60 000 placements per θ_max), the H3 shell and pressure numbers, Raymond+88 homology, penetration depths, and flat-rotation V_LSR.
- `scripts/w3_judge_phys_fig.py`: the figure.

**Output**
- `data/derived/w3_judge_phys_checks.json`.
- `figures/w3_judge_phys_posteriors.png` (200 dpi, inspected). Left: posteriors on a log axis, coloured by verdict. Right: driving P/k against assumed distance for shock runs matching 1.2–2.5 R, with the required b labelled; the closure is satisfied at every distance from 160 pc to 1.6 kpc.

**Inputs:** project files only, plus two web searches (Antlia: arXiv:2102.12599; Nereides: Fesen+2024 ApJS 272, 36 / AAS Nova 2024-07-15). Those two facts are not re-verified on disk.

**Limitations**
- Posteriors are judgement informed by the derived Bayes factors. They are not a formal likelihood.
- The Sedov, PDS, cloud-shock, H3 and LB/halo pressure numbers are order-of-magnitude. Thick-disk and LB pressures are literature-order values quoted from memory.
- The chance-alignment Monte Carlo ignores selection effects and Fesen's population prior.
- The distance scan assumes that `verify_shockmodel`'s scale-free rescaling holds, and that G is fixed at ~20.
