# w3_judge_obs: origin hypotheses for the Silk Nebula, judged on the observations

Wave 3, judge stage, tag `judge_obs`. Written 2026-11-17.

**Lens.** I weight each hypothesis only by measurements that survive null and placebo systematics and that are replicated across independent instruments (DESI, SDSS/BOSS, LAMOST, SDSS imaging, NSNS, GALEX, SHASSA/WHAM, GALFA, eRASS1/ROSAT). Arguments that need a model to turn into a constraint are discounted: the shock-model pressure closure, the Monogem pressure, C IV resonance scattering, the U budgets and the energetics. I note where they matter.

**Inputs read.**
- Wave 1: all nine notes.
- Wave 2: `verify_{paper,velocity,desi,sdss,hiwham,xrayfuv,precovery,distance,shockmodel,altmodel}`.
- Wave 3 stage 1: `w3_oii`, `w3_northsouth`, `w3_canon` (with `w3_canon_values.json`) and `w3_figs`.

**Precedence.** Wave 2 over wave 1. Wave-3 `oii` and `northsouth` over wave 2 where they re-measured the same thing:
- the [O II] doublet;
- [O II]7325;
- the DESI/SDSS [O II]/[O III] tension;
- the north–south velocity;
- the RM step.

**New work here.** Two scripts, `scripts/w3_judge_obs_edge.py` and `w3_judge_obs_edge2.py`. They make a first test, from existing NSNS and GALEX data, of two model-light geometric predictions of an edge-on shock (§4). Everything else is a judgment on teammates' numbers, cited by note.

---

## 0. Verdict

**Robust, multiply replicated measurements establish three things:**
1. A real, steady (2002–2025), 17°-long, 0.28°-wide [O III] sheet.
2. It is ~1.4–2.5 R in the south and 0.4–0.9 R in the north. It moves with the local gas (V_LSR ≈ 0), and it is **Balmer-dark**: [O III]/Hβ ≳ 8–11 at 2σ by three independent routes.
3. It has no Hα, HI, dust or X-ray counterpart.

**That is enough to exclude:**
- complete radiative shocks;
- DIG / soft photoionisation;
- the paper's specific SNR model;
- thermal-conduction fronts and hot-gas cooling flows (they would be FUV- or X-ray-bright).

**Two classes remain:**
- an **incomplete (column- or time-truncated) radiative shock**;
- a **statically heated warm layer or hard-spectrum photoionisation**.

Moderately robust evidence favours the first by about 8:1:
- the single-instrument GALEX FUV/NUV counterpart;
- the [O III]/Hβ lower limits sitting at or above the photoionisation maximum;
- the absence of any ionising or heating source.

The fragile 4363 hint points the same way.

**Who drives the shock is the weakest part of the case.** The Monogem limb rests on one geometric fact: the southern 5.5° of the spine lies at the Monogem radius, on the outer flank of its X-ray shell. Fesen+2026's independent assignment and the absence of any other adjacent hot-gas source support it. Against it:
- the arc is tangent to the Monogem circle rather than concentric with it;
- the north runs 1–4° outside it;
- the whole arc is concentric with a different point (R ≈ 21°, rms 0.15°);
- there is no Silk-specific X-ray edge.

I therefore give the Monogem limb (H1) **0.50**, not the ≥ 0.7 a modeller-led reading would suggest. A shock with some other driver (H2) gets **0.29**.

| # | Hypothesis | Posterior | Verdict |
|---|---|---|---|
| H1 | Incomplete radiative shock on the outer limb of Monogem, 250–320 pc | **0.50** | favoured |
| H2 | Incomplete radiative shock of another remnant/bubble (G205.6+12.4, Geminga-progenitor fragment, unknown older SNR, Local-Bubble-wall interaction), d ≈ 0.1–1 kpc | **0.29** | viable |
| H3 | Paper's model: independent SNR at 160 pc, R ≈ 30 pc, E ≈ 2×10⁵¹ erg, 6.8×10⁴ yr, PDS phase | **0.01** | ruled out (as specified) |
| H4 | Superbubble/chimney wall at 1–2 kpc | **0.04** | disfavoured |
| H5 | DIG / photoionised ionised-bubble wall (the paper's alternative) | **0.02** | disfavoured (DIG version ruled out) |
| H6 | Other non-shock: conduction front, turbulent mixing layer, NEI cooling of hot gas, photoionisation by WD J0736+1444 or Monogem EUV, fossil | **0.04** | disfavoured |
| H7 | *(added)* Non-shock-heated warm (~5–6×10⁴ K) layer: unspecified heating (cosmic rays, MHD-wave or reconnection heating) | **0.05** | disfavoured on physics, **not** excluded by data |
| H8 | *(added)* None of the above / unrecognised or composite mechanism | **0.05** | — |

Grouped:
- shock family (H1 + H2 + H3 + H4): **0.84**;
- static or non-shock (H5 + H6 + H7): **0.11**;
- unknown: **0.05**.

The Monogem fraction of the shock family is 0.60.

**New here (§4).**
- The northern [O III] cross-ridge profile is a **one-sided limb**: a sharp east edge and a broad west side, W/E = 1.4–2.2, p ≈ 0.02–0.05 in each of three disjoint sub-segments. This confirms an edge-on sheet curved about a centre to the west. It does not separate shock from static, or Monogem from another western driver.
- The clean southern section shows no significant asymmetry at 1′.
- The GALEX FUV/NUV ridges sit 0.04–0.10° **east** of the [O III] ridge, the sign and size expected for a front. Against GALEX sky structure this is only z ≈ 0.9–1.5. The ranking is unchanged.

---

## 1. Evidence ledger: how robust is each measurement?

**Tier A**: replicated by ≥ 2 independent instruments or methods, ≥ 4–5σ after nulls/placebos, or a robust upper limit.
**Tier B**: 2–5σ, single instrument, or a lower limit near a model boundary.
**Tier C**: ≤ 3σ, at a known systematic floor, recently revised by a pipeline bug, or estimator-dependent.

| Tier | Measurement | Replication / significance | Source |
|---|---|---|---|
| A1 | [O III]5007 filament exists: 16.9° arc, E-strand FWHM 0.28° | NSNS (z = 33 ridge), DESI (~30σ vs empirical nulls), SDSS/BOSS/eBOSS plates (up to 9.7–16σ), SDSS imaging sky (3.9σ coherent nulls), LAMOST, ≥ 3 independent amateur discoveries | `verify_desi`, `verify_sdss`, `verify_precovery`, `w3_figs` F6 |
| A2 | Brightness: south 1.4–2.5 R, north 0.4–0.9 R; paper's 10–13 R wrong | DESI 1.61 ± 0.03 / 0.73 ± 0.07 R; SDSS 1.41 ± 0.08 / 0.51 ± 0.06; SDSS imaging 2.3 ± 0.4 R (10–13 R excluded ≥ 15σ); NSNS ~2 R. Scale conflict C3 (×1.49 ± 0.25) unresolved | `w3_northsouth` §2, `verify_precovery` |
| A3 | V_helio ≈ +10…+14 km/s (V_LSR ≈ −2…+2 south); paper's +23.8…+26.8 is an artefact | DESI, SDSS, the paper's own Fig 3a against [N I] airglow (+13.6 ± 1.2); identical to WHAM WIM, local HI, Ca II/Na I | `verify_velocity`, `verify_desi`, `verify_sdss` |
| A4 | No Hα, HI, dust or X-ray filament | SHASSA 2σ < 0.77 R (Dec 12–15.6; 400 placement nulls); WHAM template < 1.0 R; GALFA < 1–1.75×10¹⁹ cm⁻²; IRIS ΔI < 0.19 MJy sr⁻¹; eRASS1 on-spine bump −0.04 ± ~0.4 on 7 cts s⁻¹ deg⁻² | `verify_hiwham`, `verify_xrayfuv`, `w3_figs` F6 |
| A5 | [O II]3727/[O III] ≈ 1.1–1.2 at Dec < 17.6, rising smoothly to 4–8 north | DESI 1.12 ± 0.08 and SDSS 1.17 (0.91–1.47) on the same Dec range. The wave-2 1.11 vs 1.84 conflict was Dec weighting. Both instruments show the gradient | `w3_northsouth` §1 |
| A6 | Present and non-varying since ≤ 2003 | Plate 1203 (2003-01-30) by two analyses; ≥ 5σ by 2005-03-03; Δ = +0.8 ± 1.3 R over 2004–07 at fixed position; DESI 2021–22 consistent | `verify_precovery`, `verify_sdss` |
| A7 | Southern spine (Dec < 17.6) lies at 12.10 ± 0.12° from the Knies two-SNR centre (Monogem R = 11.0 ± 1.1°), ~1° outside the 0.2–0.6 keV ring peak, on the outer flank | geometry (exact); X-ray W/E 1.26 is real (p = 0.007 vs sky) but generic at that radius (p = 0.37 vs rotations) | `verify_xrayfuv` |
| A8 | The paper's distance is a unit slip (kd kpc output); kinematic distances are meaningless | Table 4 reproduced to χ² = 0.05 as kpc | `verify_paper` |
| B1 | **Balmer-dark: [O III]/Hβ ≳ 8–11 (2σ)** | Three independent routes, each ≈ 2σ: (i) DESI Hβ box > 10.5, Gaussian > 18, south > 18–21 (`w3_northsouth`); (ii) SDSS [O III]/Hα > 3.0 photon with placebo → > 11 energy for case B; (iii) SHASSA Hα < 0.77 R against 1.6–2 R [O III] → > 8–10 energy for case B. SDSS Hβ alone: > 5.3–7.5 with placebo. Central value undetermined (plausibly 15–25) | `verify_desi`, `verify_hiwham` §4 |
| B2 | GALEX FUV and NUV ridge co-spatial with [O III] (south) | FUV direct contrast z = 4.6; template β 92 ± 10 CU per NSNS-R; none of 64 nulls reach β (p ≤ 0.015, the floor). F6 ridge z = 5.2 (FUV), 4.1 (NUV) vs 20 shifts. **Single instrument**; the NSNS Hα template also correlates with GALEX (z 7–11), so template z values are inflated; per-R values depend on the [O III] scale (×1.0–1.5) | `verify_xrayfuv` |
| B3 | He II 4686/[O III] < 0.026 | DESI (97.5 %); SDSS < 0.04–0.08 | `verify_desi`, `verify_sdss` |
| B4 | [Ne III]/[O III] 0.07–0.10 | DESI only | `verify_desi` |
| B5 | North [O III] redder by +4…+9 km/s (helio) | DESI, 2–3.5σ per estimator. SDSS agrees at Dec ≥ 21 and disagrees at 17.6–21 (2–3σ) | `w3_northsouth` |
| B6 | Arc geometry: south straighter than a 12° circle (R_curv 67°, 27–106°); one R ≈ 21.3° arc fits both halves (rms 0.19° / 0.12°); north 12.8–15.8° from the Monogem centre | geometric | `w3_northsouth` §4, `context` |
| B7 | **New (§4): the north [O III] profile is a one-sided limb**: sharp east edge, broad west side | W/E 1.4–2.2 in three disjoint sub-segments (W2 avoided, Hα patch masked), each at the empirical-null floor (p 0.02–0.05). The clean south shows none (p = 0.21) | `w3_judge_obs_edge2.json` |
| C1 | [O III]4363/5007 = 0.06–0.09 | 1.8σ (box) to 2.9σ (Gaussian). A null box at 4381.6 Å gives 0.048 ± 0.014. The pixel regression's 2.5th percentile is negative. SDSS 0.038 ± 0.035 is non-constraining | `verify_desi` c5 |
| C2 | [O II]3729/3726 = 1.31 ± 0.07 (DESI only) | 2.0σ from the low-density limit. The SDSS "confirmation" is a blue-end calibration artefact (its WIM control gives 1.02–1.10) | `w3_oii` |
| C3 | [O II]7325/3727 = 0.060 (0.055–0.073) ± 0.015 | Wave-1/2 value 0.026 came from a ×1.8 normalisation bug; OH-degenerate | `w3_oii` §1.4 |
| C4 | Hβ, Hγ, Hα/Hβ of the filament | Hβ at 1–2σ; Hγ/Hβ = 1.5 vs case B 0.47 | `verify_desi`, `verify_altmodel` |
| C5 | Filament [N II]/Hα ≈ 1, [S II]/Hα 0.5–0.75 | systematics-limited ×2; SDSS [S II] 6716/6731 unphysical | `verify_desi`, `verify_sdss` |
| C6 | X-ray local step at the spine | +12 % south, ~2σ, 0.04–0.25 with absorption; none north | `verify_xrayfuv` |
| C7 | [O III] intrinsic width | DESI ≤ 33 km/s (not re-verified, LSF-limited at ~98 km/s FWHM) vs paper Fig 3a ~50 km/s (S/N ~6) | `w3_canon` C7 |
| C8 | NVSS RM step in the north | placement p = 0.18: **not** Silk-specific | `w3_northsouth` §5 |
| C9 | Pre-2003 precoveries | plate 927 (2002-10-30) fragile; imaging run 3560 at 3.3σ | `verify_precovery` |

**Two wave-2 "model-independent" facts have since collapsed.** `verify_altmodel` built its exclusions on two-temperature structure (T(O II) 1.25×10⁴ vs T(O III) 4.6×10⁴ K) and high pressure (n_e ≈ 100, P/k ≈ 3×10⁶). `w3_oii` showed that:
- 7325/3727 is 0.060, not 0.026, giving T(O II) ≈ 3.5×10⁴ K, consistent with T(O III);
- the doublet is within 2σ of the low-density limit.

Neither can be used now. The non-shock verdicts must stand on FUV, X-ray, Balmer and source budgets alone (§3).

---

## 2. Hypothesis scorecard

### H1. Incomplete radiative shock on the outer limb of Monogem, 250–320 pc: **0.50, favoured**

**For (robust → less robust)**
- *Spectral class.* A Balmer-dark (B1) collisionally excited spectrum with FUV/NUV lines (B2), [O II]/[O III] ≈ 1.1 (A5), weak He II (B3), and no Hα/HI/X-ray (A4). Together these define an incomplete, ~100–150 km/s shock. Complete shocks give [O III]/Hβ ≤ 5 in every published grid, and their non-resonant FUV floor of 85–277 CU per R exceeds the measured 68–92. This applies equally to H2 and H4.
- *Position (A7).* The southern spine sits at the Monogem radius, oriented tangentially, on the outer flank of the X-ray shell. Monogem's hot gas is the only known source of the required driving pressure adjacent to the Silk.
- *Independent assignment.* Fesen+2026 put the 15° [O III] group at the Monogem E/N limb, part of a family of [O III] filaments on the E, S and NW limbs (`context` §1).
- *Hα patch / PGCC G199.09+19.95.* They lie exactly on the Monogem-radius circle (θ = 11.97° / 11.94°), a natural shock–cloud site. This is post hoc (geometric).
- *Model-dependent, discounted.* The width–brightness closure gives ρ₀v² equal to the Monogem pressure. It depends on the assumed d = 280 pc, R_c = 58 pc, G ≈ 20 and a Monogem pressure uncertain ×10 (C10).

**Against**
- *Geometry (B6).*
  - The south is straighter than a 12° circle (R_curv 67°, 0 % of bootstraps below 15°): it is **tangent to, not concentric with**, the Monogem circle.
  - North of Dec 18.5 the arc lies > 3σ outside, ending 3.7° beyond the single-SNR rim.
  - The whole arc is concentric with (94.5, +16.8), (l, b) = (194.3, +0.4), with rms 0.15° (`context` §1).

  H1 therefore needs an aspherical, locally flattened blast wave, which is common for old SNRs but is an extra assumption.
- *No Silk-specific X-ray signature.* The W/E contrast is generic at that radius (rotation p = 0.37). The local step is ~2σ and south only (C6). The Silk-sector edge is softer than other parts of the shell.
- *The north is outside the X-ray rim.* It has a different spectrum ([O II]/[O III] 4–8, [O III]/Hβ ~6–8) and a redder velocity. H1 accommodates this with a slower shock at lower pressure; that is model-consistent but not predicted in advance.
- No chance-coincidence probability has been computed for the radius match (§7).

### H2. Incomplete shock of another remnant or bubble, distance ≈ 0.1–1 kpc: **0.29, viable**

**For**
- The same robust spectral class as H1, which says nothing about the driver.
- The arc's own geometry (one R ≈ 21° arc, concentric about a point near the Galactic plane, tangent rather than concentric to Monogem) is what an unrelated, larger structure would look like.
- 160 pc is not excluded by any measurement (`verify_distance` §9). A Local-Bubble-wall interaction at 160–170 pc is allowed.

**Against**
- No second driver is seen: no catalogued SNR other than Monogem/G205.6+12.4 contains the spine (`context` §3), and no independent X-ray interior.
- G205.6+12.4 is not concentric (the spine crosses its 8.8° rim at Dec 19.0).
- The Geminga concentricity is non-physical (Geminga moves 0.2° in ~4 kyr), and no pulsar birth locus matches.
- The Origem Loop centre sits 0.8° from the curvature centre, but at the Origem Arc's ~1.7 kpc a 21° radius is ~640 pc. I read that as chance.
- The Monogem radius match in the south would have to be a coincidence.

### H3. The paper's model (160 pc, R ≈ 30 pc, 2×10⁵¹ erg, 6.8×10⁴ yr, PDS phase): **0.01, ruled out as specified**

- Every input is invalidated by robust measurements:
  - 160 pc is a unit slip for 1.6 kpc (A8);
  - the [O III]-vs-red-line velocity split is an artefact (A3);
  - the brightness is 5–8× too high (A2);
  - the red-line "shock" diagnostics are the WIM (identical in the nebula-free P4 and in WHAM);
  - R ≈ 30 pc is asserted, not fitted (the arc is 47 pc long at 160 pc, with curvature R ≈ 21° ≈ 59 pc).
- A **PDS-phase** shell has a cooled, dense, Balmer-bright shell and swept-up HI. The Silk is Balmer-dark (B1) with no Hα (< 0.77 R) and no HI (< 1–1.75×10¹⁹). That contradicts the defining property of the model.
- What survives is only the generic idea of an "independent old SNR near 160 pc", which is part of H2.

### H4. Superbubble/chimney wall at 1–2 kpc: **0.04, disfavoured**

**For**
- No absorption bracket exists. Surface brightness is distance-independent. Extinction to 1.25 kpc is only E(B−V) ≈ 0.02–0.05.

**Against**
- Its only basis was the artefact velocity (A3, A8).
- The corrected V_LSR ≈ 0 in the south is ~10 km/s below the Galactic-rotation expectation at 1.6 kpc toward l ≈ 200°, b ≈ +20° (my flat-rotation estimate: ≈ +10 km/s). That is only ~1.4σ for σ_pec ≈ 7 km/s.
- At 1.6 kpc the sheet is 470 pc long and 7.7 pc wide at z ≈ 0.55 kpc, with no catalogued superbubble or HI shell at the [O III] velocity. GS202+17+019 is not at the revised velocity.
- The south's Monogem radius match becomes a chance alignment.

### H5. DIG / photoionised ionised-bubble wall: **0.02, disfavoured; DIG ruled out**

**DIG / soft-field photoionisation.** [O III]/Hβ ≤ 1.3 at any U against ≳ 8–11 (B1), and [O III]/Hα ≳ 3 against the WIM's ~0.06. Excluded by a factor ≥ 6 through three independent routes, and also by the identical WHAM spectrum on and off the filament. **Ruled out.**

**Hard-source photoionised wall** (150 kK, or a 0.14 keV hot-gas spectrum):
- Its maximum [O III]/Hβ is ≈ 14 at any U and ≈ 9–10 at the observed [O II]/[O III]. That is only ~2σ below each of the three Balmer routes; combined, ~3σ.
- The 4363 hint (C1) adds another 2–3σ: photoionisation gives ≤ 0.011–0.02.
- GALEX FUV is above the photoionisation prediction (49–53 CU per R) by 1.5–4σ, depending on the [O III] scale.
- He II/[O III] < 0.026 excludes the hot-gas-EUV spectrum (0.06).
- The decisive strike is the absence of any source with enough ionisation parameter: WD J0736+1444 falls 27× short even coplanar. That is a catalogue-and-geometry argument and is not totally secure.

### H6. Conduction front, mixing layer, NEI cooling of hot gas, WD/EUV photoionisation, fossil: **0.04 combined, disfavoured**

- **Conduction front.**
  - Predicted FUV 540–5300 CU per R against 68–92 measured (B2). Even an optically thick C IV, which requires a column ~10× larger than these fronts carry, would not close a ×6–60 gap.
  - Face-on 2×10⁻⁴–3×10⁻³ R per front (model).
  - The 4363 strike (0.11–0.14) is only 1–2σ given C1, and the 7325 strike no longer exists (C3 = 0.060 against a predicted 0.07).
  - Disfavoured, mainly by FUV.
- **Mixing layer / NEI cooling of Monogem gas.** An X-ray filament 10³–10⁴× the Monogem surface brightness would be required; eRASS1 sees none (A4, robust). A TML at realistic inflow gives ~0.05 R (model). Ruled out.
- **WD J0736 / Monogem EUV photoionisation.** As H5 (hard). Disfavoured.
- **Fossil.** With the doublet now consistent with low density (C2), the 21–34 yr e-folding argument lapses: at n ≈ 0.5 the e-fold time is 4–7×10³ yr, and no fading is expected over 20 yr. What remains:
  - a matching window requires a start at ≥ 4×10⁴ K, i.e. shock heating, which makes it a shock remnant;
  - photoionised starts need the fragile 4363;
  - there is no ionising event in the last ~10⁴ yr (the Monogem SN was 1.1×10⁵ yr ago).

  Disfavoured; the robust evidence against it is weaker than wave 2 stated.

### H7 (added). Non-shock-heated warm layer at ~5–6×10⁴ K: **0.05, disfavoured on physics only**

**Why I added it.** `verify_altmodel` found that a warm CIE layer at 5.8–6.1×10⁴ K reproduces:
- [O III]/Hβ ~85–100 (consistent with B1, which is a lower limit);
- [O II]/[O III] 1.2–1.6;
- 4363/5007 0.097;
- FUV 86–92 CU per R, which matches B2 *without* needing C IV scattering.

**Its observational strikes were:**
- 7325/3727 ×3: now 0.08 predicted vs 0.060 ± 0.015 measured, ~1.3σ;
- NUV ×2: 84–93 predicted vs 36–57 measured, but up to ~80 with the stacked-peak estimator;
- [N II]/Hα ×2 (C5, systematics-limited ×2).

None of these is robust. What really disfavours H7 is physical: it needs ~1.7×10⁻⁵ erg cm⁻² s⁻¹ of face-on heating, at the thermally unstable peak of the cooling curve, from an unidentified source. Cosmic-ray heating falls short by orders of magnitude (my estimate).

**Discriminators:**
- line width: thermal ≈ 12–15 km/s vs the shock's 23–38;
- FUV spectrum: Si IV-dominated vs C IV;
- proper motion: none vs 75–113 mas/yr.

### H8 (added). None of the above / composite: **0.05**

- This reserves probability for model incompleteness, given that every temperature/density diagnostic is in tier C.
- A *composite*, with the south as Monogem's limb and the north a separate structure, is disfavoured geometrically: one smooth arc with no break, and smooth spectral gradients (`w3_northsouth` §1, item 4).
- Solar-system or heliospheric origins are excluded: the filament is fixed on the sky over 2002–2025 in images taken in all seasons, and the velocity is consistent across epochs in the heliocentric frame.

---

## 3. Which conclusions survive if the least-robust measurements are wrong?

| If this is wrong… | …what depends on it | What still stands (tier A/B only) |
|---|---|---|
| **4363/5007** (C1): true value ≤ 0.02 | T_e 3–6×10⁴ K; strongest anti-photoionisation line; the truncated-shock fit predicts 0.086–0.094 | Complete shocks and DIG are still out (B1, three routes). Conduction front and NEI cooling are still out (B2, A4). Hard photoionisation is then excluded only at ~2–3σ (B1) plus the source budget, and **both** the truncated-shock and warm-layer models would be in tension. My P(shock family) would fall from 0.84 to ~0.6 and P(H5) would rise to ~0.1. |
| **4363/5007** (C1): true value ≥ 0.06 at ≥ 5σ | — | Photoionisation of any kind is ruled out. Shock family + H7 → ~0.93. |
| **[O II] doublet** (C2): confirmed at 1.30 ± 0.02 | n_e ≈ 100, P/k ≈ 3×10⁶ | The spectral class survives. H1's pressure-balanced picture fails (×30–100 the Monogem pressure). That would favour a dense-cloud shock or a non-steady event (H2-like) or H8. |
| **[O II] doublet** (C2): low-density limit | — | No change: it was already discounted. |
| **Hβ central value** (C4) | model fits; the ratio's central value | Only lower limits are used. A ratio of 11–15 still excludes complete shocks; a ratio of 50–100 would shift weight toward H7 or very early truncation. |
| **X-ray step** (C6) | the "Silk = X-ray edge" statement | Already discounted. The Monogem registration rests on the radius alone. |
| **GALEX FUV/NUV** (B2) spurious | collisional-UV argument; the C IV suppression requirement; one route against 10–13 R | Complete shocks are still out (B1). Brightness is still 1.4–2.5 R (A2). The conduction-front exclusion actually strengthens. The case against hard photoionisation weakens: P(H5) ~0.05–0.08. |
| **Monogem radius coincidence** (A7) is chance | H1; the 252–318 pc distance | Spectral class unchanged. Weight moves H1 → H2 and H4. Distance becomes 0.1–1 kpc. |
| **Line width** (C7): FWHM ~50 km/s | v_s ≲ 150 km/s bound | v_s ≳ 160 km/s: mild model tension, no change in ranking. FWHM ≤ 15 km/s would instead favour H5–H7 strongly (§5 T1). |
| **North–south velocity** (B5), **RM** (C8), **[N II]/[S II]** (C5), **7325** (C3) | north physics; H7 strikes | No change in ranking. |

**Robust bottom line**, true even if every tier-C item is wrong:
1. The Silk is a real, steady, Balmer-dark [O III] sheet at local-gas velocity, with no X-ray, Hα, HI or dust counterpart.
2. Complete shocks, DIG, conduction fronts, hot-gas cooling flows and the paper's model are excluded.
3. Its southern spine sits at the Monogem radius.

Everything finer, including "column-truncated shock at 282 pc with n₀ ≈ 0.03 and b ≈ 10", is a model-consistent interpretation, not a measurement.

---

## 4. New test here: is the [O III] ridge one-sided, and is the FUV ridge outside it?

`scripts/w3_judge_obs_edge.py` → `data/derived/w3_judge_obs_edge.json`, `figures/w3_judge_obs_edge.png`.

**Why.** An edge-on shock sheet makes two geometric predictions that need no line-ratio modelling.
- The [O III] profile is one-sided: an abrupt outer edge at the front, and a √h rise to a 1/√h tail on the post-shock side.
- C IV/FUV forms nearer the front, so the FUV ridge sits ≈ 0.1° outside the [O III] ridge (`verify_shockmodel` §6).

If the driver is on the concave, Monogem (west) side, the front faces **east**.

**What each test can tell.**
- The one-sided profile is **generic to the limb of any thin curved sheet**, shocked or not. It tells which side is concave, and that the sheet is a thin limb rather than a tube. It cannot separate a shock from a static layer.
- Both Monogem and the arc's own curvature centre lie to the west, so the profile cannot separate H1 from H2 either.
- The FUV offset is the more physical test:
  - a shock predicts FUV **outside** (east of) the [O III] ridge;
  - a static warm layer predicts co-spatial UV;
  - a wall photoionised from the concave side predicts its high-ionisation zone on the **source (west)** side.

**Method.**
- Data: NSNS DR0.2 [O III] (1′ CAR, Rayleigh) and GALEX FUV/NUV 1′ blocks (`verify_xrayfuv` preparation).
- Coordinates: offsets from `CurveCoords` (+ = east).
- Profiles: medians in 0.02° bins within 0.25° along-spine blocks, then the mean over blocks.
- Models fitted at |x| ≤ 0.9°, each with a linear baseline:
  - a symmetric Gaussian;
  - a split Gaussian;
  - a projected uniform sheet with its front on the east (SHE) or the west (SHW), Gaussian-smoothed with free width.
- Errors: block bootstrap.
- **Calibration of the statistics on real sky.** The best symmetric Gaussian, and separately the best SHE sheet, are injected into 42 off-Silk placements of the same estimator.
- GALEX: tile-demeaned profiles, ridge centres from free-centre Gaussians, tile bootstrap.
- Segments:
  - **S1, Dec 12–15**: the only clean section; the W strand splits off at Dec 15;
  - S2, Dec 15–17.6: the W strand lies 0.22–0.34° west and carries more [O III] than E, so it is uninterpretable and reported only for completeness;
  - N, Dec 20–28: W2 at 0.6–1.3° west.

**Results (measured here).** The robustness checks are in `scripts/w3_judge_obs_edge2.py` → `data/derived/w3_judge_obs_edge2.json`.

*[O III] cross-ridge shape (NSNS).*
- W/E is the split-Gaussian σ_W/σ_E.
- The symmetric null injects the best symmetric Gaussian into 41–42 off-Silk placements; p is the empirical one-sided probability, with a floor of 0.023–0.024.
- Δχ²(W−E) is χ²(front-west sheet) − χ²(front-east sheet).

| Segment | W/E (bootstrap 16–84 %) | σ_E / σ_W (deg) | Symmetric null W/E 16–84 %, p | Δχ²(W−E); p vs symmetric null | Sheet vs Gaussian |
|---|---|---|---|---|---|
| **S1 Dec 12–15** (clean) | 1.16 (1.05–1.23) | 0.22 / 0.26 | 0.87–1.21; **p = 0.21** | +52; p = 0.44 | Gaussian better by Δχ² = 193. A front-east sheet as fitted would be preferred over a Gaussian in 95 % of sheet-truth injections (p = 0.047) |
| S2 Dec 15–17.6 | 0.64 | — | — | −366 | two strands; uninterpretable |
| N Dec 20–28 | **3.17** (2.67–3.56) | 0.08 / 0.26 | 0.88–1.13; p = 0.023 (floor) | +1092; p = 0.023 | sheet better by 440 |
| N1 Dec 20–21.5, Hα patch masked | **1.40** (1.29–1.65) | 0.084 / 0.118 | 0.88–1.21; p = 0.024 (floor) | +2617; p = 0.024 | sheet better by 216 |
| N2 Dec 24.25–28 (beyond W2) | **2.21** (1.88–2.72) | 0.145 / 0.32 | 0.84–1.28; p = 0.048 | +195; p = 0.048 | sheet better by 62 |
| N3 Dec 21.5–24.25 (W2 0.6–1.3° W) | **2.15** (1.95–2.34) | 0.078 / 0.168 | 0.89–1.15; p = 0.024 (floor) | +2596; p = 0.024 | ≈ equal |

*GALEX ridge centre minus the NSNS [O III] ridge centre (+ = east), same bins and fit.*
- Variants span fit windows of ±0.5/0.7/0.9° with and without a linear baseline.
- The null injects a **co-spatial** copy of the observed NSNS ridge shape, scaled to the GALEX amplitude, into 24 off-Silk GALEX placements, and refits.

| Segment | FUV Δx (variants) | Null median, sd; p, z | NUV Δx (variants) | Null median, sd; p, z |
|---|---|---|---|---|
| S1 Dec 12–15 | **+0.085°** (+0.046…+0.098) | +0.002, 0.058; p = 0.12, z = 1.4 | **+0.063°** (+0.045…+0.063) | +0.004, 0.069; p = 0.12, z = 0.85 |
| S1+S2 Dec 12–17.6 | **+0.078°** (+0.042…+0.078) | −0.010, 0.057; p = 0.12, z = 1.5 | **+0.064°** (+0.034…+0.073) | −0.011, 0.051; p = 0.12, z = 1.5 |

Dust and Hα controls: across 0.1–0.5° east minus west of the spine, there is no IRIS 100 µm excess (−0.05 vs null −0.09 ± 0.24 MJy sr⁻¹) and no SHASSA Hα excess (−0.16 vs −0.04 ± 0.50 R), from `w3_figs` F6 profiles.

**Reading.**
1. **The north is a one-sided limb.** At Dec 20–28 the [O III] profile has a sharp east edge (σ_E = 0.08–0.15°) and a broad west side (σ_W = 0.12–0.32°). This holds in three disjoint sub-segments, including one with the Hα patch masked and one beyond the W2 strand, each at the empirical-null floor (p ≈ 0.02–0.05).
   - This is the projection signature of a thin sheet curved about a centre on the **west**, seen at its limb. It confirms the edge-on-limb geometry (and the limb gain G ≫ 1) that the brightness argument assumes.
   - It is new. Wave 2 had only the curvature fit.
   - It does **not** distinguish shock from static, nor Monogem from the arc's own curvature centre: both lie west.
2. **The clean south (Dec 12–15) shows no significant one-sidedness.** W/E = 1.16 (p = 0.21). The profile is a narrow core on a broad, roughly symmetric pedestal. A uniform-emissivity limb of the fitted width is mildly disfavoured there (p ≈ 0.05).
   - Possible causes: the south is nearly straight (R_curv ≈ 67°), several ripples or sheets are superposed, or the 1′ NSNS map with a stacked, wiggly ridge blurs a ≤ 0.1° front.
   - This is neutral for every hypothesis.
3. **FUV and NUV ridges sit 0.04–0.10° east of the [O III] ridge** in every variant and both segments, which is the sign and size predicted for shock stratification (+0.1°). Against GALEX sky structure the offset is only **z ≈ 0.9–1.5 per band and segment**, and the two bands and segments are not independent.
   - It is suggestive, not a detection.
   - It argues weakly against a wall photoionised from the west, which predicts the opposite sign.
   - Reaching ≥ 3σ needs ~3× better centroid precision (T4).

**Effect on the ranking: none beyond the rounding of §0.** I do not update on signals below 2σ. The north asymmetry confirms a geometry that all west-driven hypotheses (H1, H2, and H5 illuminated from the west) already assume.

---

## 5. Decisive future observations

The "expected" columns refer to the hypotheses H1–H8 above; "static" means H5–H7.

| # | Observation (instrument / survey) | Requirement | H1 / H2 / H4 (incomplete shock) | Static (H5, H6-photo, H7) | Decides |
|---|---|---|---|---|---|
| T1 | **[O III]5007 line profile at R ≥ 30 000.** Options: WHAM Fabry–Pérot on/off pointings (1° beam dilutes the ridge to ~0.4 R; on−off removes the WIM [O III] of ~0.1–0.2 R at the same velocity), an arcminute Fabry–Pérot or spatial-heterodyne spectrometer, or a long-slit echelle crossing the ridge on a 4–10 m telescope | σ_instr ≤ 4–5 km/s; S/N ≥ 20 at the south ridge | Intrinsic FWHM **23–38 km/s**, flat-topped or double-horned, \|v\| ≤ 12–22 km/s about the pre-shock velocity (south V_LSR ≈ 0; north on the HI component at +6…+10). A one-sided red tail in the north would indicate tilt or chord geometry | FWHM ≈ 12–15 km/s (H7: thermal O at 5–6×10⁴ K) or ≤ 10 km/s (photoionised, fossil) | **Shock vs static, in one measurement.** It also settles C7 (DESI ≤ 33 vs paper ~50) |
| T2 | **Deep blue IFU/long-slit spectroscopy of the south ridge**, e.g. SDSS-V LVM (R ≈ 4000, 0.5° IFU) or WEAVE-LIFU, with off-ridge sky fields | σ(Hβ) ≤ 0.02 R and σ(4363) ≤ 0.01 R in ≥ 0.05 deg²; Hg I 4358 resolved (R ≥ 2000); Hα to σ ≤ 0.05 R | [O III]/Hβ 10–60; 4363/5007 0.07–0.10; Hα 0.1–1 R ([O III]/Hα ≥ 2) | Photoionised: [O III]/Hβ ≤ 10–14, 4363 ≤ 0.02, Hα ≥ 0.4 R. H7: [O III]/Hβ ~85–100, 4363 ≈ 0.097. Complete shock: [O III]/Hβ ≤ 5 | Converts the tier-B/C ratios to tier A. 4363 alone separates photoionised from collisional at > 5σ |
| T3 | **[O II]3726/3729 doublet**, R ≥ 5000 (e.g. WEAVE-LIFU, Keck/KCWI) | S/N ≥ 30 per component | 1.45–1.47 (low-density limit; n_e ≲ 10) | Pressure-balanced alternatives also give 1.46 | Tests the n_e ≈ 100 hint (C2). 1.30 ± 0.02 would break H1's Monogem-pressure picture |
| T4 | **FUV vs [O III] ridge position.** GALEX now gives +0.04…+0.10° at z ≈ 0.9–1.5 (§4); UVEX (planned ~2030) at higher depth, also covering Dec 22–28.25 | ≤ 0.02° centroid precision per 3° of spine | FUV/C IV ridge ≈ **+0.1° (east)**, outside the [O III] ridge | Co-spatial | Model-light geometry test of a front |
| T5 | **FUV spectroscopy** (the planned UVEX long-slit spectrograph, or a diffuse-UV smallsat or sounding rocket) | Resolve C IV 1548/1551, Si IV 1394/1403, O III] 1661/1666, He II 1640; σ ≈ 0.05 R | C IV suppressed ~4–6× relative to thin emission, doublet ratio → ~1 (optically thick), O III] ≈ 20–25 CU per R | H7: C IV/[O III] ≈ 0.03–0.04 photons, Si IV dominant (≈ 45 CU per R). Conduction front: C IV/[O III] 1.6–16 photons | Tests resonance scattering directly. Separates H1/H2 from H7 and H6 |
| T6 | **UV absorption toward background stars.** HST/COS G130M + G160M on hot WDs, sdB/sdO or early-B stars (FUV ≲ 16–17 mag) within ≲ 0.05–0.1° of the ridge, with Gaia distances spanning 150–600 pc | log N(C IV) sensitivity ≈ 13 | Extra C IV/Si IV (log N ~14+ along the tangent) at V_LSR ≈ 0 ± 20 km/s **only for stars beyond ~320 pc**, none within ~250 pc (H1). H3-like 160 pc: present beyond ~170 pc. H4: none below 1 kpc | Static warm layer: Si IV, weak C IV. Photoionised: weak high ions | **The only route to a measured distance.** WD J0747+2106 at 99.5 pc, 0.16° away, already shows nothing (log N(Si IV) < 12.2), consistent with all |
| T7 | **Proper motion** from deep narrow-band [O III] imaging at ≤ 2″ sampling, two epochs ≥ 10 yr apart (2024–26 NSNS, Fesen and amateur frames → ≳ 2035) | sharp sub-filaments; ≤ 0.2″ relative astrometry | 66–127 mas yr⁻¹ outward (E/NE) at 250–320 pc (H1); 130–200 at 160 pc; 13–20 at 1.6 kpc (H4) | 0 | Motion proves a shock. With T1's v_s it gives a kinematic-parallax distance |
| T8 | **Monogem family test (doable now).** DESI DR1 and SDSS sky fibres, with the verified pipeline, on Fesen+26's S and NW limb [O III] filaments and on random control arcs | same estimators and placebos as `verify_desi`/`verify_sdss` | H1: same class elsewhere on the rim ([O III]/Hβ ≳ 10, [O II]/[O III] ~1–2, 1–3 R, V ≈ local gas), FUV/[O III] similar | — | Raises or lowers H1 vs H2 without any model |
| T9 | **Chance-coincidence Monte Carlo (doable now)** | rigid rotations and translations of the Silk arc over the Monogem field | — | — | Quantifies the only H1-specific datum (A7) |
| T10 | **Deeper soft X-ray** (eRASS:4/5 when released) on the spine, absorption-corrected | step at x ≈ 0 to ±5 % | Hot gas bounded at or just outside the spine; a ≥ 10 % step in the south | No feature expected | Weak; supports H1 if positive |

**Ranked by decisive power per unit effort:**
1. T1 (line profile; shock vs static);
2. T2 (4363 and Hβ at S/N ≥ 5);
3. T8 and T9 (archive; H1 vs H2);
4. T6 (distance);
5. T7 (proper motion; slow but definitive).

---

## 6. Weakest links in the team analysis

1. **The Monogem registration carries the whole distance (252–318 pc) and the "favoured" H1.** It is one geometric coincidence over 5.5° of a 17° arc. It has no chance-coincidence statistic, a tangent rather than concentric geometry, and an X-ray contrast that is generic to the shell (p = 0.37). `w3_canon`'s headline "Excitation: column-truncated radiative shock … Monogem" and the distance row should stay explicitly conditional. The `verify_distance_synthesis.json` basis still cites the superseded wave-1 W/E 1.32 at 2.4σ (`w3_canon` I4, C12).
2. **Every temperature/density diagnostic sits at the DESI systematic floor**, and the faint-line pipeline has produced a run of late corrections:
   - the fibre-acceptance factor (×1/0.736);
   - SDSS residual vs total sky (×0.55);
   - uninflated bootstraps (the withdrawn [O III]/Hβ ≈ 47, R23 ≈ 128);
   - the [O II]7325 normalisation (×1.8);
   - the SDSS doublet blue-end calibration;
   - the DESI subset pulls (×1.8, cause unidentified).

   Both modelling notes fitted or argued from now-withdrawn inputs (`verify_shockmodel` used [O III]/Hβ = 40 and 7325 = 0.034; `verify_altmodel` used the two-temperature and n_e ≈ 100 arguments). More undetected errors of this kind are likely.
3. **The single most diagnostic kinematic quantity, the intrinsic [O III] width, was never re-verified.** DESI's ≤ 33 km/s is measured against a ~98 km/s LSF; the paper's ~50 km/s comes from one S/N ≈ 6 spectrum. The shock-vs-static distinction (T1) and the v_s ≤ 150 km/s bound both hinge on it.
4. **The GALEX FUV/NUV detection is one instrument with a template dependence.**
   - It uses NSNS-template regression, and the NSNS Hα template also correlates with GALEX (z 7–11) where there is no Hα filament.
   - The p value is floored at 0.015 by 64 nulls.
   - The per-R values move ×1.35 with the unresolved [O III] scale (C3).
   - The ~5× C IV suppression it requires is a fitted rescue, not a measurement.

   It nevertheless carries much of the case against photoionisation and H7.
5. **There is a common-mode dependence on NSNS DR0.2.** Its processing (continuum subtraction, large-scale filtering, PSF) is undocumented here. It is the template for the SDSS imaging precovery, the GALEX regression, the trace geometry and the edge test in §4. DESI and SDSS fibre results are independent of it; the morphology-based results are not.

---

## 7. Important things not yet done

- **T8 (Monogem family spectroscopy)** and **T9 (chance-coincidence Monte Carlo)**. Both use archival data and the existing pipelines, and they bear directly on H1 vs H2, the main open question.
- **Re-run the shock-model χ² with the corrected inputs.** These are: [O III]/Hβ as a one-sided > 11 limit, [O II]/[O III] 1.12 (south), 7325/[O III] 0.071–0.082, the doublet 1.31 ± 0.07, and 4363 with its null floor. `verify_shockmodel`'s acceptable set was built on the withdrawn values.
- **A candidate list of UV-bright background stars** (Gaia DR3 × GALEX GUVcat) within 0.1° of the ridge at 150–600 pc, to make T6 proposal-ready.
- **An independent, null-calibrated [O III] width measurement** from DESI (LSF from sky lines near 5007) and, if accessible, LAMOST MRS public spectra. This would resolve C7 before T1 is attempted.
- **H7 versus NUV.** The NUV per-R estimator spread (38–57 regression vs ~80 stacked-peak) is the only data-based strike against the warm-layer alternative, and it is unresolved.
- **A pre-registered prediction list for the north** (Dec 22–28.25, no UV data): H1 predicts FUV/[O III] 0.25–0.63 of the south (slower shock). This is testable with UVEX.
- **The edge test of §4 at higher resolution.** The first pass here is at NSNS 1′. The south needs it most, since it showed no one-sidedness at 1′. Deep amateur frames at 2–10″ already exist (e.g. AstroBin zh47r9 and the 2024–25 Fesen/NSNS data).
- **The GALEX offset test with a better centroid estimator.** One option is a joint fit of a shifted [O III] template per tile. Another is a fit to the Dec 12–15 ridge only after masking structured background. The present z ≈ 1–1.5 is limited by GALEX sky structure.

---

## 8. How these probabilities differ from the modellers'

- `verify_shockmodel` and `verify_altmodel` call the truncated shock "favoured" and nearly every alternative "ruled out". I agree on the excitation class but give the non-shock family 0.11, not ≈ 0, for three reasons:
  - their two "model-independent" facts (two temperatures; high pressure) were removed by `w3_oii`;
  - the warm-layer fit is observationally viable;
  - the hard-photoionisation exclusion rests on 2–3σ lines plus a source-budget argument.
- `w3_canon` and `context` call the Monogem limb "favoured (morphological)". I keep it favoured but at 0.60 of the shock family. The geometry is tangent rather than concentric, and the strongest supporting datum is uncalibrated against chance.
- I treat H3 as effectively dead because its defining PDS-phase property contradicts the robust Balmer-dark, HI-free observation. That does not depend on the unit slip.

---

## 9. Files (tag `judge_obs`)

- `scripts/w3_judge_obs_edge.py`: §4 test. Run with `~/.venvs/astro313/bin/python -I scripts/w3_judge_obs_edge.py`. It reads only project files (`w3_figs_f1_maps.npz`, `verify_xrayfuv_galex_blocks.npz` via `verify_xrayfuv_galex_measure.prep_band`, `silk_trace.json`). No downloads.
- `scripts/w3_judge_obs_edge2.py`: robustness checks. North sub-segments with W2 avoided and the Hα patch masked; GALEX centre offsets calibrated by co-spatial-template injection into 24 off-Silk placements, with window and baseline variants. Run with `~/.venvs/astro313/bin/python -I scripts/w3_judge_obs_edge2.py`.
- `data/derived/w3_judge_obs_edge.json`: fit parameters, bootstrap and null-injection distributions per segment, and GALEX ridge centres.
- `data/derived/w3_judge_obs_edge2.json`: north sub-segment statistics and GALEX offset null calibration.
- `figures/w3_judge_obs_edge.png` (200 dpi): NSNS cross-ridge profiles with the Gaussian and the two sheet orientations (S1, S2, N), and ridge-normalised NSNS/FUV/NUV profiles with fitted centres (Dec 12–15 and 12–17.6).
- `scripts/w3_judge_obs_edge_fig.py` → `figures/w3_judge_obs_edge_summary.png` (200 dpi). Observed W/E ratios and GALEX offsets against their null bands (16–84 %), with p and z values.
  The text panel of `w3_judge_obs_edge.png` quotes tile-bootstrap errors only. The null-calibrated values are those in the summary figure and in §4.
- `notes/w3_judge_obs.md`: this note.
