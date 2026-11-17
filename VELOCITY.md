# Silk Nebula velocities: frames, corrections, and who is right

*中文版：[VELOCITY_zh.md](VELOCITY_zh.md)*

*2026-11-17. This page puts every velocity for the Silk Nebula in one place: [O III], Hα, [N II] and [S II] from Wu et al. (2026, arXiv:2610.13329; "the paper"), DESI, SDSS/BOSS, WHAM and GALFA-HI. It states the reference frame of each value and says which values are right. All velocities are in km/s; positive means receding. The derivations are in [DISPUTE.md](DISPUTE.md) §1–4.*

**Bottom line.**
- On the southern filament (P1–P3), every line in every dataset sits at **V_helio ≈ +12…+14, V_LSR ≈ 0…+2**, once two LAMOST errors are removed.
- **The paper's [O III] is 10.1 ± 1.0 km/s too high.** The Gaussian fits are fine; the wavelength scale is wrong. It is in air wavelengths and in the telescope's frame, so the −8.3 km/s Earth-motion correction was never applied, and a +1.8 km/s airglow offset was added.
- **The paper's Hα is 5.5 ± 1.1 km/s too high** compared with its own [N II] and [S II] and with WHAM. The frame is fine, but the line centre is biased, most plausibly by over-subtracted geocoronal Hα.
- **The paper's [N II] and [S II] are right.** So are **DESI, SDSS/BOSS, WHAM and GALFA-HI.**
- **The Sun-motion (LSR) correction is not needed to decide who is right.** It adds the same constant to every line on a sightline, so it cannot create or remove an offset. It matters only for comparison with LSR-based surveys and for kinematic distances.

---

## 1. Definitions

| Quantity | Meaning | How it is obtained |
|---|---|---|
| λ_rest (air or vacuum) | Laboratory wavelength of the line. Air values are 0.028 % shorter | Must be in the same system as the spectrum. Mixing the two shifts velocities by **≈ 83 km/s** (§3) |
| V_topo | Velocity relative to the telescope | c (λ_obs − λ_rest)/λ_rest on an uncorrected scale. Airglow lines sit at V_topo = 0 by definition |
| v_corr | Earth's motion toward the target (orbit + spin) | astropy `SkyCoord.radial_velocity_correction`. The LAMOST and SDSS headers store `HELIO_RV` = −v_corr |
| V_helio, V_bary | Velocity relative to the Sun, or to the Solar-System barycentre | **V_helio = V_topo + v_corr.** Here the two differ by 0.01 km/s, so they are interchangeable |
| V_LSR | Velocity relative to the Local Standard of Rest | **V_LSR = V_helio + U⊙ cos l cos b + V⊙ sin l cos b + W⊙ sin b** |
| Standard solar motion | 20 km/s toward 18h, +30° (B1900), i.e. 18h03m50.29s, +30°00′16.8″ (J2000); (U, V, W)⊙ = (10.27, 15.32, 7.74) | Used by the paper, by `kd`, by WHAM and by GALFA-HI |

## 2. The two corrections

### 2.1 Earth motion: required, and the source of the [O III] error

- **Size on the paper's plate.** For 2024-02-01 at Xinglong (exposure mid-times 13:58–15:12 UT), **v_corr = −8.32 km/s**. It ranges from −8.03 to −8.59 over P1–P3 and the three exposures.
- **Components.** Almost all of it is Earth's orbit (≈ −8.2). Earth's spin adds at most 0.34 km/s at Xinglong.
- **Why it must be applied before comparing surveys.** The filament lies on the ecliptic (ecliptic latitude −4.5°), so v_corr swings between about −30 and +30 km/s over the year. LAMOST, SDSS and DESI observed on different dates, so their spectra can only be compared after the correction.
- **Who applies it.** The LAMOST, SDSS and DESI pipelines all apply it. WHAM and HI are published as V_LSR, so it is already included.

### 2.2 Sun motion (LSR): not needed to decide who is right

- **No effect on offsets.** On one sightline, the LSR step adds the same constant to every line. Offsets between lines and between surveys are therefore identical in V_helio and V_LSR. Every comparison below is made in V_helio.
- **When it is needed.** Only to compare with WHAM and HI, which are published in V_LSR, and for Galactic rotation and kinematic distances.
- **Choice of solar motion.** The standard and Schönrich et al. (2010) values differ by at most 0.4 km/s here:

| Region | RA, Dec (°) | l, b (°) | V_LSR − V_helio, standard | Schönrich et al. 2010 |
|---|---|---|---|---|
| P1 | 116.65, +16.60 | 203.8, +19.5 | −12.10 | −11.81 |
| P2 | 116.21, +17.74 | 202.5, +19.6 | −11.87 | −11.65 |
| P3 | 116.45, +14.72 | 205.5, +18.5 | −12.58 | −12.19 |
| P4 (no nebula) | 115.20, +16.18 | 203.6, +18.0 | −12.38 | −12.08 |

- **One trap: the correction varies along the spine.** It runs from −13.2 km/s at Dec +12° to −9.3 km/s at Dec +28.25°. North–south comparisons made in V_LSR therefore pick up as much as 4 km/s that is not in V_helio.

## 3. Lines and rest wavelengths

| Line | λ_air (Å) | λ_vac (Å) | LAMOST arm | Role here |
|---|---|---|---|---|
| [O III] 4959 | 4958.911 | 4960.295 | blue | Used in the joint DESI fit. In LAMOST its S/N is 1.8–2.3, below the paper's own S/N > 3 cut |
| [O III] 5007 | 5006.843 | 5008.240 | blue | The filament line |
| [N I] 5198, 5200 | 5197.902, 5200.257 | 5199.350, 5201.705 | blue | Airglow. It marks V_topo = 0 in the paper's Fig. 3a |
| Hα 6563 | 6562.80 | 6564.61 | red | WIM, blended with the geocoronal sky line |
| [N II] 6548, 6583 | 6548.05, 6583.45 | 6549.86, 6585.27 | red | WIM |
| [S II] 6716, 6731 | 6716.44, 6730.82 | 6718.29, 6732.68 | red | WIM |

The air-to-vacuum shift is 82.8–83.7 km/s for every line here.

## 4. What each dataset delivers

| Dataset | Observed | Wavelengths | Velocity frame | How we checked the frame |
|---|---|---|---|---|
| LAMOST MRS-N pipeline products | 2024-02-01 | vacuum | heliocentric | Headers: `VACUUM = T`, `HELIO = T`, `HELIO_RV` = +8.33 (= −v_corr). Stellar RVs on this plate: blue arm − red arm = +0.7 ± 0.5 (76 stars); both arms agree with Gaia |
| **Paper's blue arm** (Fig. 3; Table 4 [O III]) | same | **air** | **topocentric**, +1.8 offset | The [N I] airglow sits at its *air* rest wavelengths on the plotted scale. Every other frame is excluded at ≥ 7.8σ |
| Paper's red arm (Table 4 [N II], [S II], Hα) | same | not stated | heliocentric | [N II] − WHAM = −0.8 ± 1.0 |
| DESI DR1 | 2021–22 | vacuum | barycentric | Sky lines move with v_corr at slope 0.98 |
| SDSS-I/II and BOSS sky fibres | 2001–14 | vacuum | heliocentric | `HELIO_RV` = −v_corr to within 0.9 km/s on 168 groups; [O I] 5577 zero point −0.0 ± 0.3 |
| WHAM-SS DR1 Hα | — | — | V_LSR, standard solar motion | Published frame |
| GALFA-HI | — | — | V_LSR, standard solar motion | Published frame |

DESI and SDSS measure the filament as an excess over the surrounding sky. The smooth WIM largely cancels in sky subtraction and in ON−OFF, so these two surveys test the [O III] velocity, not the WIM lines. For Hα, [N II] and [S II], the independent check is WHAM.

## 5. All velocities side by side

### 5.1 By region

V_helio, with V_LSR in brackets. The V_LSR values are recomputed from V_helio with the per-region corrections of §2.2, except for the row of printed values.

| Measurement | P1 | P2 | P3 | P4 (no nebula) | Status |
|---|---|---|---|---|---|
| Paper [O III], Table 4 RV | +26.8 ± 4.3 | +26.8 ± 4.1 | +23.8 ± 3.3 | — | **wrong frame (+10.1)** |
| Paper [O III], Table 4 V_LSR as printed | [+15.5] | [+14.5] | [+11.2] | — | inherits the +10.1 |
| [O III] corrected (RV − 10.1) | +16.7 ± 4.4 [+4.6] | +16.7 ± 4.2 [+4.8] | +13.7 ± 3.5 [+1.1] | — | right |
| DESI [O III], fibres within 0.5° (n) | +11.9 ± 1.9 [−0.2] (61) | +15.2 ± 2.9 [+3.3] (9) | +9.2 ± 1.0 [−3.4] (128) | — | right |
| Paper [N II], Table 4 | +11.8 ± 1.6 [−0.3] | +12.0 ± 1.6 [+0.1] | +11.7 ± 1.9 [−0.9] | +13.3 ± 2.3 [+0.9] | right |
| Paper [S II], Table 4 | +11.8 ± 2.5 [−0.3] | +12.0 ± 1.6 [+0.1] | +11.7 ± 1.9 [−0.9] | +12.7 ± 2.1 [+0.3] | right |
| Paper Hα, Table 4 RV | +16.3 ± 1.2 [+4.2] | +17.0 ± 1.2 [+5.1] | +18.4 ± 1.4 [+5.8] | +19.7 ± 1.4 [+7.3] | **biased (+5.5)** |
| Paper Hα, Table 4 V_LSR as printed | [+4.7] | [+4.0] | [+6.1] | [+4.5] | biased; arithmetic off (§6.3) |
| WHAM Hα, nearest 1° beam | +13.4 ± 0.4 [+1.3] | +11.6 ± 1.1 [−0.3] | +13.4 ± 0.3 [+0.9] | +12.9 ± 0.5 [+0.5] | right |

P2 lies at the edge of the redder northern zone (§8), and its DESI value rests on only 9 fibres.

### 5.2 Whole filament

V_LSR here uses the mean southern correction, −12.2 km/s.

| Measurement | V_helio | V_LSR | Status |
|---|---|---|---|
| Paper, [O III] P1–P3 weighted (Table 4) | +25.5 ± 2.2 | +13.3 | **wrong frame** |
| Paper, Fig. 3a [O III] read off the plotted scale | +24.1 ± 0.9 | — | **wrong frame** (the scale Table 4 was measured on) |
| **Paper, Fig. 3a [O III] measured against its own [N I] airglow** | **+13.6 ± 1.2** | +1.4 | right, and needs no calibration |
| **DESI, on-ridge stack** (390 fibres; includes ±1.5 zero point) | **+14.4 ± 1.9** | +2.2 | right |
| DESI, joint 4959 + 5007 regression | +11.3 (10.2–13.9) | −0.9 | right |
| **SDSS/BOSS, ON−OFF, all strands** | **+9.6 ± 2.8** | −2.6 | right |
| SDSS/BOSS, second analysis, all strands | +12.3 ± 2.5 | +0.1 | right |
| SDSS/BOSS, second analysis, south (Dec < 20°) | +10.0 ± 1.8 ± 0.8 | −2.2 | right |
| **Adopted: weighted mean of the three bold values** | **+13.3 ± 1.0** | **+1.2 ± 1.0** (+1.5 with Schönrich) | χ² = 2.2 for 2 dof |
| WHAM Hα, on-filament, south (6 beams) | +13.8 ± 0.8 | +1.3 | local WIM |
| GALFA-HI on-spine peak, Dec +12…+22° | +14.2 (13.1–15.5) | +1.8 | local neutral gas |
| Ca II / Na I absorption, local group | +6…+15 | — | local gas |

![Velocity ladder in V_helio, with V_LSR on the top axis. The paper's Table 4 [O III] sits about 10 km/s above everything else. Once corrected, and when re-measured against the paper's own airglow, it joins DESI, SDSS, the paper's [N II] and [S II], WHAM Hα and local HI at +13.3 ± 1.0 km/s. The paper's Hα sits about 5 km/s high.](figures/w3_figs_f2_velocity.png)

*Figure 1. Every velocity in V_helio (bottom axis) and V_LSR = V_helio − 12.2 (top axis). The band marks the adopted +13.3 ± 1.0 km/s. Code: `scripts/w3_figs_f2_velocity.py`.*

## 6. What went wrong, line by line

### 6.1 [O III]: right fit, wrong frame

The paper's Table 4 values relate to the true velocity by

**V_Table4 = V_helio − v_corr + δ_sky = V_helio + 8.32 + 1.8 = V_helio + 10.1 ± 1.0**

where δ_sky is the airglow offset on the plotted scale. The evidence:
1. **The scale is air and topocentric.** On the paper's own Fig. 3a scale, the [N I] 5198/5200 airglow sits at +1.8 ± 0.8 km/s from its *air* rest wavelengths. On any other scale it would sit elsewhere:
   - air, heliocentric: −8.3 km/s (excluded at 12σ);
   - vacuum: +75 to +83 km/s.
2. **Table 4 was measured on that scale.** Re-measuring the Fig. 3b spectra on it reproduces Table 4 to −1.5 ± 2.4 km/s. The Gaussian fits are therefore fine.
3. **A measurement that needs no calibration gives the right value.** [O III] measured against the airglow in the same spectrum gives V_topo = +21.9 km/s. Adding v_corr gives **V_helio = +13.6 ± 1.2**.
4. **The LAMOST pipeline is not at fault.** Its stellar RVs are correct in both arms. We infer, but cannot prove, that the blue arm was re-referenced to air sky-line wavelengths after the pipeline, which also removed the heliocentric correction. The paper uses the 5198 Å sky line for flux scaling; the re-referencing step itself is undocumented.

**Result:** [O III] − [N II] goes from +13.7 ± 2.4 to **+1.1 ± 1.3 km/s**.

### 6.2 Hα: right frame, biased centre

- **The offset is internal to the red arm.**
  - Hα − [N II] = **+5.5 ± 1.1** (4.9σ); Hα − [S II] = +5.8 ± 1.2; Hα − WHAM = +4.8 ± 0.7.
  - By contrast, [N II] − [S II] = +0.1 ± 1.4 and [N II] − WHAM = −0.8 ± 1.0.
  - The offset is just as large in the nebula-free P4.
- **It is not a frame error.** A frame error would move [N II] and [S II] as well. It does not, so no frame correction can fix Hα.
- **The likely cause is geocoronal Hα.** On this night the geocoronal line sits at V_helio −10.7 to −11.4, only 22–25 km/s from the WIM, which is less than one resolution element. The sky template subtracted from it is 0.86–1.27× as bright as the published line. A model with about 30 % (24–33 %) over-subtraction reproduces both the shift and the line widths.
  - This is a model only: the spectra are unreleased.
- **What to use instead.** Take the [N II] and [S II] or WHAM velocity: V_helio ≈ +12…+13, V_LSR ≈ 0. The paper's +4…+6 is too high.

### 6.3 The V_LSR column: right formula, inconsistent numbers

- **The method is right.** The paper's solar motion is the standard one.
- **The arithmetic is not.** V_LSR − RV should be a single value per region (−11.87 to −12.58). Table 4 has −10.0 to −15.2, which leaves residuals of −2.8 to +2.6 km/s (rms 1.2).
  - In P1–P3, [S II] and [N II] have identical RVs but different printed V_LSR.
  - P4 Hα is printed as +4.5, but its RV gives +7.3.
- **Likely reason.** The two columns were probably built from different fibre sets or statistics. This is minor, but it adds 1–3 km/s of scatter on top of the [O III] and Hα errors the column already carries.

### 6.4 Consequences for the paper's interpretation

- **The line offset is not physical.** The "∼5–15 km/s (V_helio) / ∼10–16 km/s (V_LSR)" offset between [O III] and the red lines is the +10.1 frame error, with Hα's +5.5 on top. Nothing physical needs explaining, so the two-layer shell and the inhomogeneous bubble fall away.
- **No usable kinematic distance.** The paper's [O III] distances used V_LSR = +11…+16; corrected, V_LSR is +1…+5 (adopted +1.2).
  - Galactic rotation gives only ~1.6 km/s at 160 pc and ~9 km/s per kpc on this sightline, while peculiar motions are 7–10 km/s.
  - The printed distances are also in kpc, not units of 100 pc (DISPUTE.md §1–2).

## 7. Who is right and who is wrong

| Item | Verdict | Size |
|---|---|---|
| LAMOST pipeline: wavelength scale and heliocentric correction | Right | blue − red arm, stars: +0.7 ± 0.5 |
| Paper: [N II] and [S II] RVs | **Right** | vs WHAM: −0.8 ± 1.0 |
| Paper: [O III] Gaussian fits | Right | reproduced to −1.5 ± 2.4 |
| Paper: [O III] RVs labelled heliocentric | **Wrong**: air, topocentric scale | **+10.1 ± 1.0** |
| Paper: Hα RVs | **Wrong**: biased centre, most plausibly geocoronal | **+5.5 ± 1.1** |
| Paper: choice of solar motion | Right | standard |
| Paper: V_LSR column | Inconsistent arithmetic, and inherits both errors above | ±1–3, plus the above |
| Paper: [O III]-to-red-line offset; two-layer shell | **Wrong**: artefact | — |
| Paper: [O III] kinematic distance, "∼160 pc" | **Wrong**: kpc read as 100 pc, and the method gives no distance here | ×10 |
| DESI DR1 [O III] | Right | +11…+14 |
| SDSS/BOSS [O III] | Right, within errors | +10…+12 |
| WHAM Hα; GALFA-HI | Right | V_LSR ≈ +1…+2 |
| This project's −10.1 correction | Agrees with DESI, SDSS and the airglow test; cannot be checked per fibre | ±1.0 |

## 8. Open points

- **The −10.1 correction is applied in aggregate.** Fig. 3a, whose position is not stated, is the only published spectrum that shows the airglow, and the per-fibre MRS-N spectra are unreleased. Plate zero points (0.2–4 km/s) cannot be checked.
- **The Hα cause is a model.** The geocoronal intensity on this night is unknown. The bias itself does not depend on the cause.
- **The filament's velocity changes along its length.**
  - DESI gives +10.6 (10.2–12.0) at Dec < 19.5° and +20.9 (18.5–23.1) at Dec ≥ 19.5°. P2 sits at that edge.
  - In the north, WHAM gives +10.3 ± 1.3, so the northern [O III] is redder than the WIM there; it follows the northern HI.
  - SDSS and DESI disagree on the W strand: +9.7 ± 2.2 against +18.5, a difference of 8.8 ± 2.5.
  - None of this affects the +10.1 error, which belongs to the LAMOST wavelength scale.

## Sources

- [DISPUTE.md](DISPUTE.md): §1 units, §2 kinematic distance, §3 [O III] (Table C4, Figure C3), §4 Hα (Table C5, Figure C4).
- Notes: `notes/verify_velocity.md`, `notes/verify_desi.md`, `notes/verify_sdss.md`, `notes/verify_hiwham.md`, `notes/spectra_ism.md`, `notes/paper_review.md` (V_LSR column, item 6).
- Code: `scripts/dispute_oiii_values.py`, `scripts/dispute_halpha_analysis.py`, `scripts/verify_velocity_desi_sdss.py`, `scripts/verify_desi_velocity.py`, `scripts/w3_figs_f2_velocity.py`.
- Values: `data/derived/verify_velocity_summary.json`, `data/derived/dispute_oiii_values.json`, `data/derived/dispute_halpha_results.json`.
