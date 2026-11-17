# Wave-2 verification `sdss`: independent re-analysis of the SDSS sky-fibre spectra

Date 2026-11-16. Vacuum wavelengths; heliocentric velocities unless marked V_LSR (V_LSR = V_helio − 12.2 km/s here).
Intensities in Rayleigh (photon units). Line ratios are **photon** ratios unless marked "energy"
(energy = photon × λ_den/λ_num). Labels: **measured** (fitted here from the public spectra), **derived** (computed from
measured values plus a stated assumption), **literature** (a paper or a teammate's file, named).

I wrote new code for every step (`scripts/verify_sdss_*.py`). I read the wave-1 scripts only to look for bugs, and I
re-ran none of them.

## 1. Bottom line

1. **Bug in wave 1, which matters.** Wave 1 formed the ON−OFF differences and the per-fibre Δ values from the spec-lite
   `flux` column only (`spectra_sdss_lamost_sdss_stack.py` line 52: `f = d['flux'][idx]`). For a SKY fibre `flux` is the
   residual after pipeline sky subtraction. The pipeline sky model varies from fibre to fibre, and it absorbs part of
   the strand.
   * Per plate–spectrograph group, the residual keeps a median **0.47** (BOSS) or **0.64** (SDSS-I/II) of the [O III]
     ON−OFF signal that is present in the total sky `flux + sky` (26 southern groups;
     `figures/verify_sdss_skymodel_bias.png`).
   * So wave 1's strand contrasts and ON−OFF line intensities are low by about ×0.55. Its absolute per-fibre values near
     P1–P3 are roughly unaffected, because they are anchored to group means.
   * Using the total sky, the southern ON−OFF [O III] is **1.39 ± 0.07 (stat) ± 0.07 (placebo) R**, not 0.79 R. The
     W strand is **1.93 ± 0.11 R**, not 1.19 R. Per-plate strand-profile peaks are 1.6–2.5 R in the south and 0.4–0.8 R in
     the north, close to DESI (2.4 / 0.9 R) and NSNS (2.6–3.6 R). This is still 4–8× below the paper's 10–13 R.
2. **Absolute calibration (a). Units and fibre areas are confirmed. The diffuse-light correction goes the right way, but
   wave 1's factor is about 10 % too large.** Against SDSS imaging sky (DR18 `Field.sky_g/r`), the spectroscopic sky is
   brighter by 0.17–0.20 mag (SDSS-I/II) and 0.49–0.60 mag (BOSS). This gives f_eff = **0.84 ± 0.08** for 3″ fibres and
   **0.60 ± 0.06** for 2″ fibres. Wave 1 used 0.93 and 0.66.
   * Continuum, [O I]5577 and total-sky Hα each give the same BOSS/SDSS as-calibrated ratio, 1.33–1.44.
   * After my correction both instruments give [O I]5577 = 200 R (median).
   * Both f values match one effective Gaussian seeing of about 1.8″.
3. **Velocity (b). Wave 1 is confirmed within errors.** Southern [O III] V_helio = **+10.0 ± 1.8 (stat) ± 0.8 (zero
   point) km/s** (SDSS-I/II +8.7 ± 2.6, BOSS +11.2 ± 2.5), so V_LSR ≈ −2 ± 2. Wave 1 had +12.5 ± 2.6.
   * HELIO_RV = −v_corr(astropy) to within 0.9 km/s on 168 groups.
   * Zero points relative to the topocentric position: [O I]5577 −0.0 ± 0.3 km/s, Hg I 4358 +1.5 ± 0.3, [O I]6300
     −1.9 ± 0.2. [N I]5198/5200 is unusable at R ≈ 2000, because its centroid moves from −3 to +10 km/s for an assumed
     doublet ratio of 1.3–2.0.
   * The paper's +23.8…+26.8 km/s is excluded by 14–17 km/s (≈ 4σ).
   * **New:** the SDSS W strand is at +9.7 ± 2.2 km/s, which disagrees with DESI's W strand (+18.5) by 8.8 ± 2.5 km/s
     (3.5σ).
4. **Ratios (c).**
   * [O III]/Hβ (south) is about 20 (bootstrap 16–84 %: 13–48). Including the placebo Hβ systematic, the 2σ lower limit
     is **> 7**, not > 13.8.
   * **[O II]3726+3729 is strongly detected** in BOSS sky fibres: 2.0 ± 0.2 ± 0.3 R in the south.
     [O II]/[O III]5007 is **1.84 energy (95 %: 1.46–2.37)** in the south and 7.3 in the north, against DESI's
     1.12 (Dec < 19), 1.49 (all) and 6.2 (northern arc).
   * [O II]3729/3726 = 1.29 (68 %: 1.07–1.59), consistent with DESI's 1.30.
   * R23 ≈ 76 (> 20 at 2σ, derived), consistent with DESI's 58.
   * [O III]4363/5007 = 0.038 ± 0.035 energy, which neither confirms nor contradicts DESI's 0.056. Hg 4358 cancels in
     ON−OFF (0.16 ± 0.15 R).
   * He II/[O III] < 0.04 (bootstrap 97.5 %).
5. **North (d).** SDSS-north and the DESI patch measure different structures, so they cannot be compared directly.
   Measured on the Hα patch itself, SDSS gives [N II]/Hα = 0.78 (0.70–0.89) and [S II]/Hα = 1.06 (0.95–1.19), close to
   DESI's 0.65 and 0.85.
   * The high north-strand ratios (1.75 and 1.86 here; 1.85 and 1.77 in wave 1) come from a tiny Hα difference
     (0.1–0.4 R). That difference is no larger than the Hα placebo floor (rms 0.21 R), and it moves with the OFF choice:
     [N II]/Hα ranges 1.1–2.4 across variants.
   * The 2σ statement is only [N II]/Hα ≳ 0.7. This is consistent with DESI's northern arc, ≈ 1.0.
   * OH residuals (< 2 % of the OH lines, ≤ 0.015 R at 6555/6579) and geocoronal Hα (uniform within an exposure, so it
     cancels) are excluded as causes.
6. **Precovery (e).** The earliest plate with ≥ 3σ on-strand [O III], including a per-plate placebo systematic, is
   **SDSS plate 1203, MJD 52669 = 2003-01-30 (UT 04:42–05:20)** at (116.48, +21.89), with ON−OFF
   1.11 ± 0.11 ± 0.31 R (3.4σ). Plate 857 (2002-02-09) is formally 3.5σ but falls to 0.5σ once placebo scatter is
   included. Highly significant detections begin with plates 2054 (2005-03-02, 4.2σ) and 2074 (2005-03-03, 9.7σ). The
   full list is in §7.

## 2. What I did

Pipeline (all in `scripts/`, run with `python -I`; only FITS tables are read from downloaded files):

1. `verify_sdss_extract.py`
   * Reads all 9,092 spec-lite SKY-fibre files on disk: 7,958 DR17 and 1,134 DR19, of which 2,452 are SDSS-I/II 3″
     fibres and 6,640 are BOSS-spectrograph 2″ fibres.
   * Keeps the **total sky `flux + sky`**, the sky model, ivar and wdisp in eight windows (3690–6800 Å) on a common
     10⁻⁴-dex grid.
   * Computes synthetic g and r fluxes of the total sky (Doi+2010 responses).
   * DR19 FPS fields have an empty PLATEID, so the FIELD number is used. Plates with SPECTROGRAPHID = −1 are split by
     fibre number.
2. `verify_sdss_groups.py`
   * Measures fibre distances to trace v1.0 (E, W and W2 strands) with my own great-circle code.
   * Defines 302 release × plate × MJD × spectrograph groups. 250 are "dark": no more than 0.75 mag brighter than the
     instrument median at 5070 Å.
   * For each group's clipped-mean total sky, fits:
     - [O I]5577 and 6300 (centroid, intensity), continuum and g/r surface brightness;
     - a red absolute model: CALSPEC solar model broadened to the wdisp LSF, plus OH lines at topocentric positions I
       identified empirically (6524, 6535, 6546, 6555.4, 6579.1, 6599, 6606.2, 6636.5, 6670, 6702, 6709.5, 6746 Å), plus
       Hα, [N II] and [S II];
     - blue continuum templates.
   * Computes the astropy heliocentric correction wherever TAI is in the header.
3. `verify_sdss_diff.py`
   * For each fibre: total sky − clipped mean of the same group's OFF fibres.
   * Fits fixed-velocity Gaussians (+13 km/s; LSF from wdisp), a linear baseline and ε_c × (group continuum template,
     which carries the solar Hβ/Hα absorption). In the red it adds ε_OH × (group OH template). Near 4363 it adds a
     topocentric Hg I 4358 term.
   * Lines: [O II]3726, 3729; Hγ; [O III]4363; He II 4686; Hβ; [O III]4959 (tied photon 1/3.01) and 5007; Hα; [N II]6548
     (tied 1/3.06) and 6583; [S II]6716, 6731.
4. `verify_sdss_results.py`
   * ON = within r_on of the nearest strand. OFF = 0.45–2.0° away in the same group, with the Hα-patch 30 % contour
     excluded.
   * Each group's error comes from the robust scatter of its OFF fibres. Groups are combined with common weights.
   * Errors and ratio percentiles come from a 2,000× bootstrap over groups.
   * [O III] velocity: a free-velocity fit to the weighted group stack, minus the weighted [O I]5577 zero point, with a
     400× group bootstrap.
5. `verify_sdss_zeropoint.py`: Hg 4358, [N I], [O I]5577 and 6300 zero points per group.
6. `verify_sdss_abscal.py` and `verify_sdss_geocorona.py`: imaging-sky, zodiacal and BOSS/SDSS anchors; [O I]5577;
   total-sky Hα against solar depression.
7. `verify_sdss_robust.py`: placebo strands (trace shifted ±1…±3.5° in RA), plane + ON fits, and the per-plate
   precovery table. `verify_sdss_plate_placebo.py`: per-plate placebo (±0.6…±1.2°).
8. `verify_sdss_north.py`: eight selection variants (ON 0.10/0.15/0.25°; OFF 0.3–1.0 / 0.45–2.0 / 1.0–2.5°; OFF from
   the whole plate; Hα patch kept in OFF; no templates). `verify_sdss_ohtest.py`: topocentric-frame OH residual test.
   `verify_sdss_check4363.py`: stack-level 4363/Hg fit. `verify_sdss_pregions.py`: absolute values at P1–P4.
9. `verify_sdss_figures.py`: figures.

Systematic floors: placebo rms of ON−OFF from shifted traces with at least 20 ON fibres
(`data/derived/verify_sdss_placebo_summary.json`).

| Region | [O III] | Hβ | Hα | [N II] | [S II] | [O II] | 4363 | He II |
|---|---|---|---|---|---|---|---|---|
| South (R) | 0.07 | 0.09 | 0.18 | 0.05 | 0.10 | 0.32 | 0.04 | 0.05 |
| North (R) | 0.05 | 0.08 | 0.21 | 0.09 | 0.10 | 0.14 | 0.04 | 0.03 |

Bootstrap errors alone underestimate the uncertainty of every low-ionisation line.

## 3. Results: claim → verdict → my value

| # | Wave-1 claim | Verdict | My value (diffuse-corrected with f = 0.84 / 0.60 unless noted) |
|---|---|---|---|
| a1 | Spec-lite units 10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹ per fibre, vacuum, heliocentric | CONFIRMED | BUNIT, VACUUM = T. 1 R fills a fibre with 5.24e-17 (3″) or 2.33e-17 (2″) erg s⁻¹ cm⁻² at 5007 Å (derived) |
| a2 | Fibre diameters 3″ (SDSS-I/II) and 2″ (BOSS, SDSS-V plates and FPS) | CONFIRMED | As-calibrated BOSS/SDSS ratio: continuum 1.44, g 1.46, r 1.33, [O I]5577 1.40, Hα 1.33 (measured), the aperture-loss ratio expected for ~1.8″ seeing |
| a3 | Diffuse correction × f_psf, median 0.78 (SDSS 0.93, BOSS 0.66), ±20 % | REVISED | Direction correct. f_SDSS = 0.84 ± 0.08, f_BOSS = 0.60 ± 0.06 (imaging-sky anchor, g/r mean; ± covers epoch differences). Wave-1 intensities ×0.90 from this alone |
| a4 | Corrected sky continuum 21.75 AB at 5040 Å = dark APO sky | CONFIRMED, weak test | 21.75 (SDSS) and 21.72 (BOSS) AB at 5070 Å. The field is on the ecliptic: zodiacal light alone is V ≈ 22.1–22.35 at the observed elongations (Leinert+98 via STScI ETC), about 60 % of the continuum, so this test does not discriminate f at the 20 % level |
| a5 | [O I]5577 median 286 R | REVISED | 201 R (SDSS) and 200 R (BOSS) with my f (zenith-reduced ~170 R): normal |
| a6 | WHAM Hα 2.9–4.4 R test | UNVERIFIABLE | Total-sky Hα (Galactic + geocoronal) after the solar-absorption model is 3.3–3.4 R, and falls with solar depression to 1.8 R at > 60°. That is below WHAM Galactic alone, so the absolute Hα model is biased low by 1–2 R (template-dependent: 4.1 vs 5.5 R as-calibrated between two OH line lists). Not usable as an anchor |
| 1 | Southern ON−OFF [O III] 0.79 ± 0.02 R | REFUTED (bias) | **1.39 ± 0.07 ± 0.07 R** (r_on = 0.15°). 1.61 (r_on 0.10), 1.06 (0.25), 1.24–1.42 across OFF choices. Plane + ON fit: 1.05 ± 0.03 |
| 1b | W strand 1.19 (ON−OFF), amplitude 1.40 R | REVISED | 1.93 ± 0.11 R ON−OFF. Per-plate profile peaks 2.4 R (plates 4497, 4488, 11083) |
| 1c | P1/P2/P3 mean within 0.3°: 1.42 / 1.80 / 1.70 R | CONFIRMED | 1.49 ± 0.06 / 1.78 ± 0.05 / 1.87 ± 0.09 R (absolute, ±0.2 R template zero point). Brightest single fibres 3.0–3.65 R. Off-strand (1–3°) median 0.07 R (wave 1: 0.32–0.38) |
| 2 | V_helio [O III] south +12.5 ± 2.5 ± 0.7 | CONFIRMED | **+10.0 ± 1.8 ± 0.8 km/s**. W +9.7 ± 2.2, E +11.3 ± 4.1, north +16.4 ± 5.9, all +12.3 ± 2.5. 5007/4959 free ratio 3.15 |
| 2b | Zero point from [O I]5577 (blue −0.2, red −5.3 km/s) | REVISED | 5577: −0.0 ± 0.3; Hg 4358: +1.5 ± 0.3 (interpolated to 5007: +0.8 ± 0.3); [O I]6300: −1.9 ± 0.2 (not −5.3). The OLS slope of v5577 on HELIO_RV is −1.58 because of outliers; Theil–Sen gives −1.05 |
| 3 | [O III]/Hβ south > 13.8 (2σ) | REVISED | Hβ = 0.059 ± 0.038 ± 0.09 R. Ratio 23 (16–84 %: 13–48). 2σ lower limit 10 (bootstrap only), 7.5 (+ placebo MAD), 5.3 (+ placebo rms); plane + ON fit > 14.6 (formal). Balmer-dark: yes |
| 3b | [O III]/Hα south 4.7 ± 1.0 | REVISED | 2.7 (16–84 %: 2.2–3.4; + placebo: 1.9–4.5). Plane + ON fit: Hα 0.10 ± 0.06, so ≳ 6. Poorly constrained |
| 3c | [N II]/Hα south 0.96 ± 0.25; [S II]/Hα 0.62 ± 0.23 | CONFIRMED (large errors) | 0.70 (0.55–0.92) and 0.62 (0.48–0.82), plus systematics. [S II]6716/6731 ON−OFF = 2.5 (south) and 1.8 (north), above the physical maximum of about 1.45, so [S II] carries a ~0.1 R systematic |
| 3d | North: [O III]/Hβ 7.4 ± 2.5; [N II]/Hα 1.85 ± 0.36; [S II]/Hα 1.77 ± 0.37 | REVISED (overstated precision) | [O III]/Hβ 6.4 (3.2–12; > 2.7 at 2σ). [N II]/Hα 1.75 (1.3–2.6 boot), 1.1–2.4 across variants; 2σ with placebo ≳ 0.7. [S II]/Hα 1.86 (1.4–2.7), 1.4–2.2 across variants |
| 3e | All-strand stack 5007 0.58, Hβ 0.04, Hα 0.17, [N II] 0.23 R | REVISED | 0.92 ± 0.09, 0.072 ± 0.035, 0.42 ± 0.08, 0.49 ± 0.05 R |
| 3f | DESI [O II]/[O III] 1.12–1.49 (energy) | REVISED (independent test) | SDSS south 1.84 energy (1.46–2.37 at 95 %; ± 0.35 with placebo); W strand 2.04; E strand 2.01; north 7.3 (DESI northern arc 6.2). Consistent with DESI-all (1.49), about 2σ above DESI Dec < 19 (1.12) |
| 3g | DESI [O II]3729/3726 1.30 | CONFIRMED | 1.29 (68 %: 1.07–1.59) south; 1.13 (0.99–1.29) all; 1.12 (1.00–1.27) patch |
| 3h | DESI R23 ≈ 58 | CONFIRMED | R23(energy) ≈ 76 south; > 20 at 2σ with Hβ systematics (derived; [O II] from BOSS groups only) |
| 3i | DESI 4363/5007 0.056; Hg 4358 risk | UNVERIFIABLE (not constraining) | Stack fit with topocentric Hg: 4363 = 0.053 ± 0.020 (boot) ± 0.04 (placebo) R. Ratio 0.044 ± 0.037 energy. Hg ON−OFF 0.16 ± 0.15 R, so street-lamp Hg cancels |
| 3j | DESI He II/[O III] < 0.034 | CONFIRMED | Photon 0.005 (97.5 % < 0.036); < 0.08 with placebo |
| 4 | LAMOST LRS stack [O III] 0.88, Hβ −0.08 | not re-tested | — |
| 5 | Precovery: (not given in wave 1) | NEW | Earliest ≥ 3σ: 2003-01-30 (plate 1203); see §7 |

## 4. (a) Absolute calibration in detail

* **Spectrophotometry convention.** SDSS/BOSS spectra are calibrated to point-source total flux (standard stars versus
  PSF magnitudes). A uniform source is therefore over-estimated by 1/f, where f is the fraction of a standard star's
  light entering the fibre during the exposure (seeing, guiding, ADR). Multiplying by f < 1 is the correct direction.
* **Imaging anchor** (`verify_sdss_abscal.py`, T1). For 250 dark groups, the spectroscopic g/r sky (synthetic, as
  calibrated, divided by fibre area) is compared with the median DR18 imaging `Field.sky_g/sky_r` within 1° (fields
  2000–2007).
  - SDSS-I/II: −0.17 (g) and −0.20 (r) mag, so f = 0.85 and 0.83.
  - BOSS: −0.60 (g) and −0.49 (r) mag, so f = 0.57 and 0.63.
  - f is larger in r, as expected from better red seeing.
  - Epoch differences add up to ±0.1 mag: imaging is 2000–07; SDSS-I/II spectra 2001–08; BOSS 2010–14, mostly weaker
    solar maximum 24.
* **Internal consistency.**
  - A Gaussian PSF of FWHM 1.83″ gives f(3″) = 0.84, and 1.74″ gives f(2″) = 0.60. One effective seeing explains both.
  - The as-calibrated BOSS/SDSS ratios of continuum, [O I]5577 and Hα (1.33–1.46) match f_SDSS/f_BOSS = 1.40.
  - With these f, [O I]5577 is 201 R (SDSS) and 200 R (BOSS).
  - Wave 1 assumed 1.5″ for all SDSS-I/II plates (PlateX `seeing50` = 0), so f = 0.93. That is too high.
* **Zodiacal anchor** (T2). I tried to measure the solar-spectrum component of the continuum from Fraunhofer-line
  depths and regress it on the Leinert+98 zodiacal prediction for each date (λ′ = 73–177°, β = 0–14°). It **failed as
  a calibrator**:
  - correlation coefficient only 0.3–0.5;
  - SDSS slope 2.4 ± 0.5 with a negative intercept;
  - the line-depth estimator depends on the LSF.

  I report it only as context: corrected continuum 21.72–21.75 AB at 5070 Å against zodiacal V ≈ 22.1–22.35. So the
  "dark-sky 21.7–22" check of wave 1 is not a strong test.
* **Hα and geocorona** (T4). The total-sky Hα falls with solar depression:

  | Solar depression | 18–30° | 30–45° | 45–60° | 60–90° |
  |---|---|---|---|---|
  | Total-sky Hα (R) | ≈ 6.5 | 3.6 | 2.6 | 1.8 |

  This is the geocoronal trend. But the 60–90° value lies below WHAM's Galactic 2.9–4.4 R, so the absolute Hα (in the
  solar Hα absorption trough) is biased low by 1–2 R. The test cannot constrain f.

## 5. (b) Velocity

* **Sign convention.** HELIO_RV + v_corr(astropy) = −0.06 km/s (median), with |Δ| ≤ 0.9 km/s on 168 groups. Topocentric
  lines sit at λ₀(1 − HELIO_RV/c); the Theil–Sen slope of v5577 on HELIO_RV is −1.05.
* **Zero points** (median over dark groups, measured − expected, km/s):

  | Line | Zero point (km/s) |
  |---|---|
  | [O I]5577 | −0.0 ± 0.3 (MAD 3.7) |
  | Hg I 4358 | +1.5 ± 0.3 |
  | [O I]6300 | −1.9 ± 0.2 |
  | [N I] (ratio 1.6) | +4.1 ± 0.8 |

  The [N I] value changes from −3.0 to +9.8 for an assumed ratio of 1.3–2.0. Interpolating linearly between Hg 4358 and
  5577 gives +0.8 ± 0.3 km/s at 5007 Å. I use 5577 and carry 0.8 km/s as the systematic.
* **[O III] centroids.** I fitted group-bootstrapped stacks (total-sky ON−OFF; [O III]4959 tied):

  | Stack | V_helio (km/s) |
  |---|---|
  | South | +10.0 ± 1.8 |
  | South, SDSS-I/II | +8.7 ± 2.6 |
  | South, BOSS | +11.2 ± 2.5 |
  | W strand | +9.7 ± 2.2 |
  | E strand | +11.3 ± 4.1 |
  | North | +16.4 ± 5.9 |
  | All | +12.3 ± 2.5 |

* **Comparisons.**
  - WHAM Hα: V_helio ≈ +12.3 (literature, ISM teammate).
  - Paper LAMOST [O III]: +23.8…+26.8 (literature). Excluded by 14–17 km/s.
  - DESI (literature, DESI teammate): E +13.0 agrees. The W strand at +18.5 does not: SDSS minus DESI = −8.8 ± 2.5 km/s.
    Both SDSS instruments agree with each other on W (+9.1 ± 3.1 and +10.6 ± 2.4).

## 6. (c, d) Ratios, [O II], 4363, He II and the north

* **Stacks** (`figures/verify_sdss_stacks.png`). [O II]3726/3729 is a clean, partly resolved doublet in all three
  regions (south, north, patch), with a flat baseline. [O III] is strong in the south. Hβ is barely visible on the strand
  and clear on the patch.
* **[O II]** comes only from BOSS and SDSS-V (2″ fibres, from 3600 Å).
  - South: 0.87 + 1.12 = 1.99 R against [O III] 1.45 R in the same groups. Photon ratio 1.37 (95 %: 1.08–1.77), energy
    1.84 (1.46–2.37).
  - The placebo [O II] floor in the south is 0.32 R (one shift reached +0.9 R), so energy ratio 1.84 ± 0.35 overall.
  - This agrees with DESI's all-fibre 1.49 and is about 2σ above DESI's Dec < 19 value (1.12). The ratio rises steeply
    to the north: 7.3 energy, against DESI's northern arc 6.2.
* **Hα patch (new).**
  - [O II] 3.55 R, Hβ 0.37 R, Hα 1.41 R, so [O II]/Hβ (energy) ≈ 13.
  - Hα/Hβ = 3.8 photon = 2.84 energy, i.e. case B with no reddening.
  - [S II]/Hα 1.06, [N II]/Hα 0.78, [O III]/Hβ 1.8.
  - [O II]/Hβ ≈ 13 together with [S II]/Hα ≈ 1 is far above photoionised WIM ([O II]/Hβ ≲ 4). It is typical of
    radiative-shock spectra. This is my qualitative comparison; I did not run models.
* **North reconciliation.**
  - Variants (`figures/verify_sdss_north.png`, `data/derived/verify_sdss_variants.json`): north [N II]/Hα runs
    1.11–2.36 and [S II]/Hα 1.42–2.18. The Hα difference is 0.25–0.74 R, compared with a placebo rms of 0.21 R.
  - The plane + ON fit gives Hα 0.14 ± 0.06, [N II] 0.40 ± 0.03 and [S II] 0.44 ± 0.04 R.
  - OH test (topocentric-frame stack): the ON−OFF residuals of OH 6555 and 6579 are −0.009 ± 0.015 and
    −0.001 ± 0.015 R. That is 0.2 % of the lines. Adding or removing OH terms changes Hα and [N II] by < 0.005 R.
  - Geocoronal Hα is common to ON and OFF fibres in the same exposure.
  - DESI's 0.65 / 0.85 is the Hα-patch component. On the patch, SDSS gives 0.78 (0.70–0.89) and 1.06 (0.95–1.19) from
    ON−OFF, and 0.99 / 1.29 from the plane fit. That is within 1–2σ of DESI, with SDSS systematically higher.
  - The [O III]-strand north has a genuine [N II] and [S II] excess (5–7σ over placebo) and an Hα that is not
    significantly detected. Its ratios are therefore lower limits of order 0.7–1 (2σ). This agrees with DESI's northern
    arc ([N II]/Hα ≈ 1.0, [S II]/Hα ≈ 1.05, photon).
  - **Neither number is wrong. Wave 1 compared different structures, and quoted the SDSS north ratios with errors that
    omit the dominant Hα systematic.**

## 7. (e) Spectroscopic precovery

Per plate–MJD (`data/derived/verify_sdss_precovery.csv`, `verify_sdss_precovery_placebo.csv`;
`figures/verify_sdss_precovery.png`). S/N_pl includes a per-plate placebo rms from traces shifted ±0.6…±1.2°.

| Plate | MJD (first exposure) | UT date | ON−OFF [O III] (R) | Placebo rms (R) | S/N_pl | Profile fit (R) |
|---|---|---|---|---|---|---|
| 857 | 52314 | 2002-02-09 | 0.18 ± 0.05 | 0.34 | 0.5 | 0.16 ± 0.19 (not detected) |
| **1203** | **52669** | **2003-01-30** | **1.11 ± 0.11** | **0.31** | **3.4** | **1.14 ± 0.18** |
| 2054 | 53431 | 2005-03-02 | 1.20 ± 0.07 | 0.28 | 4.2 | 1.45 ± 0.29 |
| 2074 | 53432 (coadd MJD 53437) | 2005-03-03 | 1.45 ± 0.05 | 0.14 | 9.7 | 1.86 ± 0.12 |
| 2890 | 54495 | 2008-01-30 | 1.36 ± 0.06 | 0.24 | 5.4 | 1.86 ± 0.25 |
| 2915 | 54497 | 2008-02-01 | 1.44 ± 0.05 | 0.09 | 13.9 | 1.66 ± 0.17 |
| 4456 | 55536 | 2010-12-06 | 0.38 ± 0.04 | 0.11 | 3.3 | 0.45 ± 0.07 |
| 4497, 4492, 4502, 4488, 4470, 4473, 4462, 4464 | 55564–55865 | 2011-01-03 … 2011-10-31 | 0.54–2.51 | 0.05–0.18 | 3.2–13.9 | — |
| 8954 | 57453 | 2016-03-06 | 0.73 ± 0.06 | 0.13 | 5.2 | — |
| 11102, 11084, 11104, 11085, 11103, 11082, 11083 | 58401–58515 | 2018-10-10 … 2019-02-01 | 0.48–2.33 | 0.05–0.25 | 4.6–16.3 | — |

* **Below 3σ (S/N_pl):** 857, 888, 2078, 2079, 4480, 4485, 8955, 10285 and 11095.
* **Untested:** 927, 1263, 1582 and 1921, because too few placebo positions fell on the plate. 1921 (2004-11-08) has a
  profile fit of 1.88 ± 0.35 R (formal 5.4σ) with no fibre within 0.15°.
* **Earliest:** **SDSS plate 1203, MJD 52669 (2003-01-30, exposures 04:42–05:20 UT, three exposures of about 12.4 min),
  3.4σ**. The first ≥ 5σ detection is plate 2074 (MJDLIST starts 53432 = 2005-03-03).
* SDSS-V DR19 groups did not pass the dark-sky cut and do not appear in the table.
* For context, the precovery teammate's imaging-sky detection covers 2002–2007.

## 8. New findings

1. **Sky-model absorption bias in wave 1** (§1.1). It affects every wave-1 SDSS ON−OFF intensity and strand amplitude,
   by ×0.47 for BOSS and ×0.64 for SDSS. Ratios are less affected, but not immune.
2. **SDSS [O II] is an independent detection** of strong [O II] on the strand ([O II]/[O III] ≈ 1.8 energy in the south,
   ≈ 7 in the north). It supports DESI's high R23 and the "Balmer-dark, collisionally excited" picture. The
   north-to-south gradient in [O II]/[O III] is large: ×4.
3. **The Hα patch has a shock-like spectrum** ([O II]/Hβ ≈ 13 energy, [S II]/Hα ≈ 1.1, unreddened Balmer decrement).
4. **The W-strand [O III] velocity** (+9.7 ± 2.2) disagrees with DESI's +18.5 at 3.5σ. The DESI verifier should check
   this.
5. **[S II]6716/6731 ON−OFF exceeds the physical limit** (1.8–2.5) on the strand. With Hα at the placebo floor, SDSS
   cannot measure the strand's own low-ionisation ratios better than about a factor of 2.
6. **Spectroscopic precovery from 2003-01-30**, secure (≥ 5σ) from 2005-03-03.

## 9. Limitations

* f_eff rests on a comparison of imaging sky and spectroscopic sky at different epochs (±0.1 mag). No same-night
  comparison was possible.
* All differential low-ionisation intensities (Hβ, Hα, [N II], [S II], [O II]) carry placebo systematics comparable to
  or larger than the bootstrap errors. Ratios quoted with bootstrap ranges alone are optimistic. The placebo sample is
  small (10–12 shifts) and some shifts may hit real unrelated structure.
* ON/OFF classes depend on trace v1.0. The plane + ON estimator gives a lower [O III] (1.05 vs 1.39 R) because the
  0.15–0.45° wings enter the background.
* R ≈ 1800–2200: line widths are not measurable; the 3726/3729 doublet is only partly resolved; the [N I] zero point is
  ratio-dependent.
* The [O II]/[O III] ratios use BOSS groups only. R23 mixes [O II] from BOSS groups with Hβ from all groups.
* The zodiacal table was transcribed from the STScI ETC page (§10) through a web-fetch summary. Four rows were checked
  twice.

## 10. Files

* **Scripts:** `scripts/verify_sdss_{lib,extract,groups,zeropoint,diff,results,abscal,geocorona,robust,plate_placebo,north,ohtest,check4363,pregions,figures}.py`.
* **Derived data** (`data/derived/verify_sdss_*`):
  - windows: `windows.npz`
  - fibre and group tables: `fibres.csv`, `groups.csv`, `templates.npz`, `zeropoints.csv`
  - differential fits: `diff_{base,offnear,offfar,plate,inclpatch,notmpl}.npz`
  - region results: `results_*_r*.json`, `stack_*.npy`
  - calibration: `abscal.json`, `abscal_groups.npz`, `geocorona.json`
  - systematics: `robust_base.json`, `placebo_summary.json`, `variants.json`, `ohtest.json`, `ohtest_*.npy`,
    `check4363.json`, `skymodel_bias.json`
  - precovery: `precovery.csv`, `precovery_placebo.csv`
  - P regions: `pregions.json`
  - zodiacal table: `zodi_table_stsci_etc.json`, transcribed from
    https://etc.stsci.edu/etcstatic/users_guide/1_ref_9_background.html (accessed 2026-11-16; Table 2, derived from
    Leinert et al. 1998 Table 17)
* **Figures:** `figures/verify_sdss_{skymodel_bias,abscal,stacks,velocity,north,precovery}.png`.
* **Inputs used (no new downloads except the zodiacal table text):**
  - `data/raw/sdss_dr17_speclite/`, `data/raw/sdss_dr19_speclite/`
  - `data/raw/sdss_dr18/fields_box.csv`
  - `data/raw/filter_curves/`, `data/raw/calspec_sun/`
  - `data/derived/silk_trace.json`
