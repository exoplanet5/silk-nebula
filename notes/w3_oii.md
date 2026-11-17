# Wave-3 `oii`: forensic test of the [O II] 3729/3726 ≈ 1.30 doublet ratio

Date: 2026-11-17. Tag `oii`. Vacuum wavelengths and heliocentric (= barycentric) velocities throughout. Ratios are energy
ratios; for 3729/3726 and 6716/6731 photon and energy ratios differ by < 0.2 %. LDL means low-density limit.
Labels: **measured** (from public spectra by my own code), **derived** (computed from measurements plus a stated
assumption), **simulated**, and **literature / teammate** (with the file named).

## 1. Verdict

**n_e ≈ 100 cm⁻³ is not established.** The "independent confirmation" by SDSS is an artefact. DESI by itself gives
3729/3726 = 1.31 ± 0.04 (bootstrap) ± ~0.07 (realistic). That is 2.0 σ below the low-density limit, which these
data cannot exclude. The two DESI [O II] items behind the wave-2 "tension" with the pressure-balanced shock model were
both wrong or overstated.

1. **SDSS/BOSS 1.29 is instrumental (measured).** At 3727 Å the BOSS doublet is blended: the separation is
   2.7 LSF σ. A fixed-velocity fit is therefore degenerate with the line centroid; it moves by −0.0095 in ratio per
   km/s (simulated). Two calibration errors at the blue end bias the wave-2 setup low:
   * the centroids sit 20–40 km/s blueward of the expected velocity;
   * the true LSF is 14–35 % broader than `wdisp`.

   The decisive control is the Galactic **WIM [O II]**, isolated in the absolute total sky by regression on WHAM Hα.
   The WIM must be at the LDL, and its own [S II] is (1.40–1.53). With the wave-2 setup its [O II] gives
   **1.02 (0.88–1.18) and 1.10 (1.00–1.28)**.

   The Hα patch has [S II] = 1.43–1.47, also low density, yet gives 1.03–1.14 with the same setup.

   With v and σ free, every SDSS region returns 1.33–1.67: south 1.41 (1.26–1.61), all 1.50, north 1.56, WIM 1.33 and
   1.48. **SDSS therefore constrains nothing below ±0.15 and does not confirm DESI.**
2. **DESI 1.30 is instrumentally robust but statistically marginal (measured + simulated).**
   * My own regression pipeline reproduces it: 1.31 for all fibres, south 1.31, north 1.31. The result is
     insensitive to the fit model:
     - the wave-2 replica gives 1.304;
     - LSF ±10 % gives 1.310–1.312;
     - v ± 5 km/s gives 1.299–1.319;
     - free v and σ give 1.311, with σ/σ_R = 1.04 and v = +10 (equal to the [O III] velocity);
     - the robust (clipped) variant gives 1.293.
   * Simulations at the DESI LSF (separation 3.9 σ) show that no plausible instrumental effect moves a 1.46 doublet
     below ~1.40:

     | Effect | Fitted ratio (true 1.46) |
     |---|---|
     | LSF ±20 % | 1.445–1.461 |
     | v ± 10 km/s | 1.435–1.483 |
     | E/W strand blending, Δv ≤ 20 km/s | 1.459 |
     | ON−OFF/WIM width, λ or throughput mismatch | 1.43–1.49 |
     | Fe I 3727.6 continuum dip | +0.01 |
     | Spectroperfectionism sub-pixel phase (1D specter replica) | 1.443–1.468 |

     Only a skewed LSF (Gauss–Hermite h3 = −0.1 gives 1.40) or σ(3729)/σ(3726) ≥ 1.1 comes close to 1.30.
   * Disjoint fibre subsets scatter far beyond their bootstraps. QSO/z > 0.5 gives 1.18 against 1.62 for z < 0.5
     galaxies (pull 3.2). Dark gives 1.24 against bright 1.54. Dec bands give 1.10–1.48.
   * The rms pull is 1.8. For the same subsets the DESI filament [S II] 6716/6731 moves the **opposite** way: QSO 1.69
     against galaxies 0.76. The scatter is therefore systematic, not a physical density difference.
   * Realistic result: **R = 1.31 ± 0.07 (68 %)**, so n_e = 100–160 (+80/−60) cm⁻³ at 1.25 × 10⁴ K for the four
     atomic-data sets. At 2 σ, n_e < 260–360 cm⁻³, and the LDL is within 2.0 σ. The profile likelihood (v, σ free)
     gives a 68 % interval of 1.21–1.42, with Δχ²(1.46) = 1.6.
   * The DESI Hα-correlated (WIM-like) component, a low-density control, gives [O II] 1.64 (1.39–1.91) and
     [S II] 1.40 (1.33–1.49). This shows no low bias in DESI, but it is only a ±0.25 test.
3. **[S II] cannot rescue either side (measured).**
   * The filament's own [S II] 6716/6731 is 1.34 (1.15–1.55) in DESI (all; south 1.52, 1.16–2.11). That means
     n_e ≲ 260 cm⁻³ at 1 σ: weak, and consistent with both.
   * The SDSS filament [S II] is unphysical: south 2.56, north 1.63. These values are systematics-dominated (OH),
     as wave 2 found.
   * The Hα patch (SDSS 1.43–1.47) and the WIM (1.40–1.53) are at the [S II] LDL.
4. **New: the wave-1/2 [O II]7325 flux is low by ×1.8, a normalisation bug.**
   * Both waves fitted λ7320 and λ7330 as one profile, 1.0·G(7320) + 0.8·G(7330), and reported its amplitude, which
     is the λ7320 flux alone, as "[O II]7320+7330" (`scripts/spectra_desi_regspec.py`, `scripts/verify_desi_linefit.py`
     group `o2_7325`).
   * A four-line fit with atomic internal ratios and free OH residual terms gives:
     - [O II]7325/[O III] = 0.079–0.082 (all) and 0.071 (south), against wave 2's 0.036 and the box value 0.081;
     - **7325/3727 = 0.060 (0.055–0.073; south 0.064)**, not 0.026.
   * So T(O II) ≈ 3.5 (2.4–5.5) × 10⁴ K at low density, or 3.1 × 10⁴ K at n_e = 100. That is not 1.25 × 10⁴ K.
   * This removes the truncated-shock model's second "failure" ([O II]7325 × 1.8; the model gives 0.059–0.068 of
     [O III]). It also weakens the 7325-based arguments in `verify_altmodel` (see §6).
   * OH 7316/7329 residuals correlate with the filament template at comparable amplitude, so this ratio carries a
     ±0.015 systematic.
5. **If n_e ≈ 100 were real (derived):**
   * P/k ≈ 2.2 n_e T ≈ 3 × 10⁶ K cm⁻³, 30–100 × the Monogem pressure.
   * The 1.5–2.5 R of [O II] would come from a layer only **19–31 AU** deep. That is a filling factor of ~10⁻⁴ of the
     1.37-pc filament width.
   * Recombination and cooling times would be ~1,500 and ~1,000 yr, so the dense gas would need continuous
     re-supply.
   * The driver would have to be a ram pressure of ≳ 100 km/s into n₀ ≳ 2 cm⁻³ clumps. In a pressure-balanced model
     the surface brightness scales as n₀v³, so this would require a covering factor of ~10⁻².
   * At n_e ≈ 1, by contrast, the same [O II] needs **0.9–1.5 pc** of path, which is the observed filament width, at
     P/k ≈ 3 × 10⁴ (the Monogem pressure).
   * The low-density reading is self-consistent; the high-density one is contrived.

## 2. Data and methods (own code; wave-2 products used only as inputs)

| Step | Script | Output |
|---|---|---|
| DESI pixel-level IV regression spectra (3727, 5007, [S II], 7325, Hα windows), geometry and neighbour-Hα instruments, 23 + 5 subsets, two weightings, HEALPix-nside-32 block bootstrap (400), robust variant | `scripts/w3_oii_desi_regress.py [robust]` | `data/derived/w3_oii_desi_regspec{,_robust}.npz`, `w3_oii_desi_lsf.json` |
| Doublet fits: [O III] kinematics per subset, then [O II]/[S II] models w2/kin/±σ/±v/vfree/vsfree, bootstrap refits, white-noise MC, profile likelihood | `scripts/w3_oii_desi_fit.py [robust]` | `w3_oii_desi_fits{,_robust}.json`, logs `w3_oii_desi_fit*_log.txt` |
| SDSS/BOSS ON−OFF stacks rebuilt from total sky (wave-2 ON/OFF definitions), fits incl. ε×OFF template | `scripts/w3_oii_sdss.py` | `w3_oii_sdss.json`, `w3_oii_sdss_stacks.npz`, `w3_oii_sdss_lsf.json` |
| WIM control: total-sky regression on WHAM-SS DR1 Hα (between and within groups) | `scripts/w3_oii_sdss_wim.py` | `w3_oii_sdss_wim.{json,npz}` |
| BOSS blue-end wavelength zero point and LSF from street-lamp Hg I 3650/3663/4047/4358 and [O I]5577 | `scripts/w3_oii_sdss_wavecal.py` | `w3_oii_sdss_wavecal.json` |
| Fit simulations (DESI, BOSS) and 1D specter ex2d replica | `scripts/w3_oii_sim.py` | `w3_oii_sim.json` |
| 5-level O⁺/S⁺ atoms (PyNeb data files on disk: Kal09/T07 × FFT04/Z82-WFD96; TZ10 + IFF05), inversions, OH-aware 7325 re-fit, implications | `scripts/w3_oii_atomic.py` (shared code `w3_oii_common.py`) | `w3_oii_atomic.json` |
| Key numbers, subset-pull inflation factor | `scripts/w3_oii_summary.py` | `w3_oii_summary.json` |
| Figures | `scripts/w3_oii_figures.py` | `figures/w3_oii_{spectra,systematics,summary}.png` |

Notes on method:
* **DESI LSF.** I used the median SPARCL `wave_sigma` of 430 pilot coadds: **0.692 Å at 3727** (16–84 %: 0.655–0.711),
  0.696 at 5008 and 0.644 at 6724. The per-fibre 3729/3726 σ ratio is 1.000 ± 0.006.
  The fitted [O II] widths (σ/σ_R = 1.04–1.07) match σ_R plus the [O III] intrinsic σ of 10–17 km/s.
* **DESI line representation.** In specter `ex2d`, `PSF.projection_matrix` uses delta-function spots at the grid
  wavelengths; R comes from √iCov (source text read from GitHub, `data/raw/w3_oii_src/`). An exact 1D replica shows
  that an off-grid line comes out as a point-sampled Gaussian at the true centre. The "two-bin split" picture, which
  would bias 1.46 to 1.34 at the south velocity (3726 on a pixel centre, 3729 on a pixel edge), is **wrong**: the
  replica gives 1.443–1.468 at all phases.
* **Weighting lesson.** A GLS fit with the shrunk bootstrap covariance, applied to the same replicates, gives
  artificially narrow intervals and a spurious −0.03 per km/s velocity dependence; my first pass did this. All
  numbers here use diagonal weights (bootstrap variance after removing the common line-strength mode). With these the
  DESI ratio is velocity-independent, as simulated.
* **BOSS LSF.** Median `wdisp` gives σ = 1.006 Å at 3728 (0.858 Å pixels). Fitted widths of narrow lines in the coadded
  total sky are larger than `wdisp`:
  - [O I]5577: × 1.13;
  - Hg 4358: × 1.31;
  - Hg 4047: × 1.47;
  - Hg 3663: ~ × 1.5 (low S/N);
  - [O II] doublet (v and σ free): × 1.14 (filament, all regions), × 1.26–1.34 (WIM), × 1.28–1.32 (patch).

  Street-lamp Hg lines may be pressure-broadened, so 5577 and the [O II] doublet itself are the clean measures.
* **BOSS blue-end zero point** (centroid minus expected, km/s):

  | Line | Offset (km/s) |
  |---|---|
  | [O I]5577 | −1 ± 4 |
  | Hg 4358 | +1 ± 1 |
  | Hg 4047 | +8 ± 2 |
  | Hg 3663 | −36 ± 5 |
  | Hg 3650 | −11 ± 5 |
  | [O II] doublet, v free, WIM within | ≤ −42 (grid edge) |
  | [O II] doublet, v free, WIM between | −18 |
  | [O II] doublet, v free, patch | −31 |
  | [O II] doublet, v free, north | −40 |
  | [O II] doublet, v free, all | −29 |
  | [O II] doublet, v free, south | −6 |

  Expected velocities for the [O II] rows come from [S II] (+11…+16, red arm) or from DESI. The [O II] and [S II] of
  the same patch gas cannot differ physically by 30–45 km/s. The blue-end offset varies between plate sets (the south
  groups are least affected).

## 3. Results tables

**DESI, geometry instrument, model `kin`** (v and σ_int from the subset's own [O III]). The bootstrap interval is 68 %.

| Subset | n | v_[O III] (km/s) | 3729/3726 | 68 % (boot) |
|---|---|---|---|---|
| all | 7,345 | +11.0 | **1.310** | 1.261–1.343 |
| south (Dec < 19.5) | 3,852 | +10.5 | 1.311 | 1.206–1.347 |
| north | 3,493 | +15.5 | 1.315 | 1.246–1.431 |
| Z > 0 / Z < 0 | 1,945 / 3,900 | | 1.39 / 1.24 | 1.32–1.46 / 1.14–1.39 |
| QSO + z > 0.5 / z < 0.5 galaxies | 4,216 / 3,129 | | 1.18 / 1.62 | 1.12–1.22 / 1.54–1.80 |
| dark / bright | 4,651 / 2,694 | | 1.24 / 1.54 | 1.17–1.29 / 1.45–1.87 |
| Dec 12–15.5 / 15.5–17.5 / 17.5–19.5 / 19.5–22 / 22–25 / 25–29 | | | 1.30 / 1.18 / 1.28 / 1.10 / 1.48 / 1.47 | 8–17 blocks each |

The wave-2 per-fibre IV with the same data gives 1.298 (all) and 1.301 (south), reproduced. Wave 1's "on-ridge 1.45
vs off-ridge 1.12" is not significant: Z > 0 against Z < 0 gives a pull of 1.0.

**SDSS/BOSS ON−OFF stacks** (fixed = wave-2 setup: v = +13, σ = `wdisp`):

| Region | Fixed | v free | v, σ free (v, σ/wd) |
|---|---|---|---|
| south | 1.29 (1.15–1.47) | 1.36 (v = +4) | 1.41 (1.26–1.61) (+2, 1.14) |
| all | 1.16 (1.10–1.24) | 1.44 (−17) | 1.50 (1.44–1.58) (−18, 1.14) |
| north | 1.13 (1.06–1.21) | 1.49 (−24) | 1.56 (1.47–1.63) (−26, 1.14) |
| Hα patch 30 % | 1.14 (1.05–1.24) | 1.45 (−17) | 1.67 (1.49–1.79) (−24, 1.28) |
| WIM, between groups | 1.10 (1.00–1.28) | 1.23 (−6) | 1.33 (1.20–1.53) (−12, 1.34) |
| WIM, within groups | 1.02 (0.88–1.18) | 1.44 (−30, edge) | 1.48 (1.17–1.83) (−30, 1.26) |

* The ε × OFF-sky template changes nothing (1.29, 1.16), so ON/OFF throughput or continuum mismatch is excluded.
* The fitted WIM [O II] is 0.63–0.89 R per R of WHAM Hα (as calibrated). Using the same regression's Hα slope gives
  [O II]3727/Hα ≈ 1.0 (energy), a normal WIM value.

**[O II] 7320+7330 (DESI, four-line fit; ratio to 3726+3729 of the same subset)**:

| Subset | Model | 7325/[O III] | 7325/3727 | T(O II), n_e = 1 (Kal09) |
|---|---|---|---|---|
| all | no OH / OH free | 0.079 / 0.082 | 0.058 / 0.060 (0.055–0.073) | ~3.3–3.5 × 10⁴ K |
| south | no OH / OH free | 0.072 / 0.071 | 0.065 / 0.064 (0.047–0.067) | ~3.5–4 × 10⁴ K |

Wave-1/2 values: 0.037/0.036 (λ7320 amplitude only) give 0.026. Multiplied by 1.8 that becomes 0.047.

## 4. Simulations (noiseless unless stated; `data/derived/w3_oii_sim.json`)

Fitted ratio for a true 1.46, at true v = +11 km/s and σ_int = 14 km/s:

| Effect | DESI (sep = 3.9 σ) | BOSS (sep = 2.7 σ) |
|---|---|---|
| wave-2 fit setup | 1.454 | 1.437 |
| LSF true/fit 0.90 / 1.10 / 1.20 | 1.461 / 1.454 / 1.445 | 1.490 / 1.429 / 1.398 |
| v_fit − v_true = +10 / +20 km/s | 1.435 / 1.405 | 1.366 / 1.274 |
| σ(3729)/σ(3726) = 1.05 | 1.416 | 1.391 |
| two strands Δv = 20 km/s | 1.459 | 1.458 |
| LSF skew h3 = −0.10 / +0.10 | 1.402 / 1.521 | 1.303 / 1.643 |
| broad wings (20 %, 2σ) | 1.436 | 1.416 |
| ON−OFF with WIM 3×, OFF λ offset ±0.05 Å | 1.432–1.489 | 1.352–1.583 |
| ON−OFF with WIM 3×, OFF LSF ±5 % | 1.457–1.466 | 1.413–1.505 |
| noise, peak S/N 17: median (16–84 %) | 1.461 (1.38–1.55) | 1.461 (1.38–1.54) |

The same table for a true 1.30 is in the JSON (e.g. BOSS +10 km/s gives 1.219). The ratio estimator is unbiased in
the median at the observed S/N.

## 5. What changes for the team

* **Paper/judges.** The [O II] doublet cannot be used to claim n_e ~ 100 or P/k ~ 3 × 10⁶. State it as
  "3729/3726 = 1.31 ± 0.07 (DESI), within 2σ of the low-density limit; SDSS uninformative (blue-end calibration)".
* **verify_shockmodel.** Tension 9 (doublet) falls from "inconsistent with every pressure-balanced shock" to a ≤ 2σ
  tension. Tension 10 ([O II]7325 × 1.8) is **resolved**: it was the normalisation bug. Measured 7325/[O III] is
  0.071–0.082 against the model's 0.059–0.068.
* **verify_altmodel.** The "two temperatures" argument, T(O II) = 1.25 × 10⁴ against T(O III) = 4.6 × 10⁴, is no
  longer supported: T(O II) ≈ 2.4–5.5 × 10⁴ K is consistent with T(O III) ≈ 3–6 × 10⁴ within the errors.
  - The 7325/3727-based strikes against the conduction front (predicted 0.07), the mixing layer (× 2) and the fossil
    (0.07) disappear; the measured value is ~0.06.
  - The doublet-based strikes against all pressure-balanced alternatives also weaken to ≤ 2σ.
  - Their other constraints (4363, FUV, [O III]/Hβ, X-ray, energetics) are unaffected and should carry those verdicts.
* **verify_sdss (side finding).** The BOSS coadd LSF is ~13 % broader than `wdisp` even at 5577. Single-line fluxes
  fitted with σ = `wdisp` are therefore ~6 % low at 5577, and more in the blue. This is a small correction to the SDSS
  [O III] brightness (1.39 → ~1.5 R) that I did not propagate.
  - WHAM regression gives an Hα slope of 1.06 (0.96–1.15) within groups and 1.29 (0.99–1.73) between groups,
    as-calibrated R per WHAM R. With wdisp-σ Hα fits this suggests f_BOSS ≈ 0.7–0.9 rather than 0.60.
  - That is an unverified by-product, limited by geocoronal Hα, the WHAM beam and the solar absorption model.
* **The SDSS [O II]/[O III] ratios** (south 1.84) come from fixed-v fits. The doublet *sum* is insensitive to the
  centroid error: v-free sums change by < 1 %. σ-free sums are 4–10 % higher because the LSF is broader than `wdisp`.
  The SDSS [O II] intensities are therefore at most ~5–10 % low, comparable to the [O III] effect above.

## 6. Limitations

* The DESI uncertainty inflation (× 1.8) rests on five disjoint-subset pairs. The class/program-dependent systematic
  that drives it is not identified: Redrock residuals at a fixed observed wavelength, or bright-time sky-subtraction
  structure, are candidates. Spatial n_e variations along the filament cannot be excluded either: Dec bands span
  1.10–1.48 with 8–17 blocks each.
* An instrumental LSF skew in DESI's blue camera (h3 ≈ −0.1 to −0.2) would bias 1.46 toward 1.30–1.40. I had no DESI
  arc/sky line at 3727 Å to test LSF asymmetry: SPARCL provides only `wave_sigma`, and DESI sky spectra/resolution
  matrices are at NERSC, which is blocked.
* The BOSS calibration offsets rely on weak Hg 3650/3663 blends and on the doublet centroids themselves. Sky-model or
  coadd resampling effects were not modelled from first principles.
* The [O II]7325 lines lie 1.5–3.5 Å from OH 7316/7329, whose residuals correlate with the filament template. The OH
  terms are partly degenerate with λ7330, so 7325/3727 carries ±0.015. Kal09 Υ extend to 10⁵ K; T07 stops at 2 × 10⁴ K
  and is not used above that.
* The atomic-data spread is LDL 1.455 (T07) to 1.50 (Kal09); n_e(1.31) = 101–157 cm⁻³ at 1.25 × 10⁴ K.
* Implication numbers assume solar O/H, all O as O⁺, and a Λ ~ 2 × 10⁻²⁴ cooling coefficient. They are order-of-
  magnitude.

## 7. Provenance

* **Inputs on disk:**
  - DESI windows: wave-1 `spectra_desi_windows/`, wave-2 `verify_desi_windows/`;
  - wave-2 `verify_desi_linefluxes.npz` ([O III] and Hα amplitudes) and `verify_desi_coadd_meta.npz` (Z, f(λ));
  - SDSS `verify_sdss_windows.npz`, fibre and group tables, `verify_sdss_diff_base.npz` (group weights);
  - spec-lite files `data/raw/sdss_dr1{7,9}_speclite/` (FITS tables only, `python -I`);
  - WHAM-SS DR1 `data/raw/wham_ss_dr1/`;
  - PyNeb data files `data/raw/verify_shockmodel_atomic/pyneb/`.
* **New downloads (small, logged):**
  - `data/raw/w3_oii_src/`: specter `psf.py` and `ex2d.py`, desispec `resolution.py` (raw.githubusercontent.com,
    2026-11-16, `MANIFEST.txt`);
  - `data/raw/w3_oii_hanuschik/`: VizieR UVES sky-atlas tables 1–5. The tables start at 3740 Å, so they give no
    information at 3727; not used.
* **Note on an abandoned approach:** a direct absolute fit of the total sky with the CALSPEC solar template fails
  (χ²/dof ≈ 23, template mismatch). It is kept only behind `DO_ABS = False` in `w3_oii_sdss.py`; the WHAM regression
  replaced it.
