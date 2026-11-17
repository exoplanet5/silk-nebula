# verify_precovery: independent check of the SDSS imaging-sky precovery, and a dated precovery table

Wave-2 adversarial verification, tag `precovery`. Date 2026-11-16 (resumed run; the first pass of this agent wrote the
scripts and data products below but no note; this pass re-checked them, added steps 7-8 and wrote this note).

Labels: **measurement** = my fit to archival data; **derived** = computed from measurements plus a stated
assumption; **literature/teammate** = a published value or another agent's file (named).

## 0. Verdicts at a glance

| Claim (wave 1) | Verdict | My value |
|---|---|---|
| 10. SDSS imaging `sky_g` (2002-07) tracks the [O III] filament, ≈ 4σ against RA-shift nulls (A = 1.35 ± 0.19 × NSNS) | **CONFIRMED, with a lower and more template-dependent significance** | Without NSNS, using unit E+W Gaussian ridges (FWHM 0.28°), 19 runs at Dec 12-20.5: **A = 2.31 ± 0.38 R** peak [O III]5007 (PSP sky). Against coherent same-shift nulls: **3.9σ**; against the off-peak shift scan: 4.4σ; against random off-window nulls: 6.6σ. Across estimators, templates and masks the coherent-null significance is 1.4-4.6σ. The r and i controls are null. |
| 3 (part). SDSS imaging-sky brightness 2.8-4.9 R | **REVISED** | Peak of a 17′-FWHM ridge per strand: **2.3-2.6 R** (PSP / FramesSub / Frames). Runs at Dec 17.2-17.7 give 3-4 R. The 4.0-4.9 R "1′ peaks" come from scaling NSNS morphology, not from SDSS. Lower bound with every confounder regressor: 1.65 ± 0.31 R. Upper range 4-6 R only if the PSP sky is modelled as temporally smoothed. 10-13 R is excluded at ≥ 15σ on 17′ scales. |
| Imaging S/N "17 of 18 runs positive" | CONFIRMED (16 of 19) | Three negative runs are all at Dec 12.5-14.4, the filament's southern tip. |
| Epoch subsets agree (2002-04: 1.40 ± 0.21; 2004-07: 1.11 ± 0.43 × NSNS) | **CONFIRMED** | At fixed position (Dec 17.2-17.7): 2003-12/2004-02 gives 3.15 ± 0.79 R; 2006-11/2007-02 gives 3.95 ± 1.05 R; Δ = +0.8 ± 1.3 R. DESI-normalised trend over 2002.9-2007.1: −0.16 ± 0.13 per yr. My 2004-05 "season" is negative, but only because those runs cross Dec 12.5-13.4. |
| Wave-1 error model ("√2 for pairwise-repeated `sky_g` values"; placebo scatter) | **REFUTED as stated, conservative in effect** | DR18 `sky_g` values are not repeated: 0 of 1482 consecutive pairs are identical. PSP is, however, temporally smooth (high-frequency rms 0.36 × Frames). Wave-1 placebos were fitted alone in the same window, so they pick up the real ridge's side lobes and inflate the null scatter. That explains wave 1's χ²/dof = 0.41. |
| PSP estimator responds linearly to a 17′ filament (assumed by wave 1, untested) | **PARTLY CONFIRMED / UNVERIFIABLE at frame level** | Run 4458 (2004-02-16): PSP field sky vs my own frame-level field medians has slope 0.72 (Frames 0.87; regression-diluted lower limits). The frame-level ridge fit of a single run is inconclusive: A = 7.2 R against re-fitted shift nulls of 4.5 ± 3.2 R, and the same structure appears in r and i. |
| Earliest ≥ 3σ precovery (task d; not stated by wave 1) | **NEW** | Spectroscopic: **2002-10-30** (SDSS plate 927, MJD 52577, 5.0σ in my per-plate fit, only 3 on-ridge fibres). Imaging: **2002-12-31** (SDSS run 3560, 3.3σ in g, r/i null). Confirmed by two independent analyses: **2003-01-30** (plate 1203, MJD 52669; 3.5σ here, 3.4σ in verify_sdss). |

## 1. What I did (all code mine; no wave-1 imaging code imported)

Scripts, in run order:
- `scripts/verify_precovery_lib.py`: geometry, ridge models, filter responses, robust least squares.
- `verify_precovery_prep.py`: per-field templates, masks, observing conditions.
- `verify_precovery_download.py`: frame downloader.
- `verify_precovery_fit.py`: Field-level fits, 70+ variants.
- `verify_precovery_frames.py` and `verify_precovery_frames_sed.py`: frame-level test of run 4458.
- `verify_precovery_stats.py`: jackknives, seasons, joint g+u.
- `verify_precovery_scan.py`: RA-shift scan.
- `verify_precovery_table.py`: dated precovery table.
- `verify_precovery_extra.py`: cross-calibration, SED, null calibration, epochs, confounders, fibres.
- `verify_precovery_pspkernel.py`: PSP temporal response.
- `verify_precovery_figs.py`: figures.

Data (downloaded by me, logged):
- SkyServer DR18 `Field` table, all five bands and all three sky estimators (`sky_b` = PSP, `skyFrames_b`, `skyFramesSub_b`), RA 109-124, Dec 6-30, 20 435 fields. Query URL and access time (2026-11-16T11:36:40Z) are in `data/raw/verify_precovery_sdss/fields_box_full.json`; the column descriptions are in `field_columns.csv`.
- SDSS DR17 `frame-[ugri]-004458-*.fits.bz2` from https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/4458/: 708 files, 2.2 GB, 0 errors. The URL lists and per-file UTC log are in `data/raw/verify_precovery_sdss/frames_4458/{urls*.txt,download_log.jsonl}`.
- speclite SDSS-2010 u, i, z response curves in `data/raw/verify_precovery_sdss/filter_curves/`. g and r came from wave-1's `data/raw/filter_curves/`.

### 1.1 Template that does not use NSNS (task a)
- Unit-amplitude Gaussian ridges along the silk_trace v1.0 E spine and W strand, FWHM 0.28°, with a Gaussian roll-off (σ = 0.15°) beyond the strand ends.
- Each field's template is the **mean over a 9 × 7 grid covering its real 13.5′ × 9.8′ footprint**, oriented along the scan direction derived from neighbouring field centres.
- Variants:
  - E only; E+W+W2; FWHM 0.6° and 1.0°;
  - footprint median, or the centre value;
  - the **DESI forward-model amplitudes** (`spectra_desi_results.json`, trace HWHM widths);
  - the **SDSS sky-fibre strand amplitudes** (`spectra_sdss_lamost_sdss_onoff.json`).
- Amplitudes are in R of [O III]5007 photons, counting 4959 at 1/2.98. The g response is 1.800 × 10⁻³ nmgy arcsec⁻² per R; wave 1 also had 1.80 × 10⁻³.
- The u band is reported in R of the [O II]3726+3729 doublet (2.377 × 10⁻³ nmgy arcsec⁻² per R). r, i and z use the full filament line spectrum (DESI line ratios).

### 1.2 Field-level model (per run)
- y[f, c] = S_f + g_c + Σ_{p=1,2} h_{c,p} t^p + A·T[f, c]:
  - S_f is a free common mode per field number (the six camcols are simultaneous);
  - g_c are camcol offsets;
  - h_{c,p} are quadratic per-camcol drifts.
- Iterative 4σ MAD clipping.
- Window: |d| ≤ 2.5° from the E spine. Quality ≥ 2.
- Excluded fields:
  - near Tycho-2 VT < 7 stars, plus Pollux, Procyon, Castor, Gomeisa and δ Gem; radius 0.12-1° by magnitude;
  - within 1° of a planet.
- A run is used only if it truly crosses the ridge: ≥ 6 fields at d < −0.7°, ≥ 6 at d > +0.7° and ≥ 4 at |d| < 0.3°.
- That leaves 28 crossing runs, 2001-11-12 to 2007-02-21: 19 "core" runs (Dec 12-20.5) and 9 "north".

### 1.3 Nulls and errors (task b)
- **Off-window nulls.** For each run, the whole analysis is repeated with template and window centred on the spine shifted by |s| = 3-6° in RA. Fields within 1.2° of the real strands are excluded. This gives 10-20 nulls per run, 398 in all.
- **In-window nulls.** The template is shifted by 0.8-2.0° and fitted *jointly* with the real one.
- Per-run error e = max(formal, std of the off-window nulls). With these errors the core χ²/dof is 1.14, so they are adopted.
- Run combinations are inverse-variance weighted, and errors are scaled by √(χ²/dof) when that exceeds 1.
- Three combined null tests:
  1. *random*: each run draws an independent null; 20 000 draws;
  2. *same-shift*: the same sky offset for every run, which keeps any coherent large-scale structure. This is the most conservative;
  3. a *shift scan*: template and window moved together from −3° to +6°.
- Injection-recovery:
  - Field-level injections at the off-window positions (1 and 3 R) are recovered exactly. The model is linear and the clipping mask does not change, so this is a bookkeeping check and not an estimator test.
  - The informative injection is a **template-mismatch test**. The DESI model (varying amplitude and width) is injected at 315 off-window run positions and recovered with the fixed gauss028 template: recovered peak 1.86 ± 0.26 R (16-84 %: 1.64-2.05). Over the same Dec range the DESI model's mean peak is ≈ 2.3-2.5 R, so the fixed template under-recovers by about 10-25 %.

### 1.4 Controls and confounders (tasks b, c)
- **Bands:** u ([O II]3727 plus [Ne III]), r (Hα+[N II], weak here), i and z (no strong lines) as controls.
- **Estimators:**
  - PSP `sky_b`;
  - `skyFrames_b`, the frames-pipeline sky per frame;
  - `skyFramesSub_b`, the same after object subtraction.
- **Regressors added per run** (`conf_*` cases):
  - airmass;
  - twilight (10^(0.4(h_sun+18)));
  - zodiacal light (exp(−|β|/5°) and exp(−|β|/15°));
  - gegenschein (Gaussian, σ = 10°, about the antisolar point);
  - IRIS 100 µm;
  - Finkbeiner Hα (r band).
- Observing geometry (Moon, Sun, planets, ecliptic) was computed with astropy at APO for each field's `mjd_g`.
- Enlarged star masks (`bigmask`): 2° for V < 1.5, 1° for V < 4.5, 0.5° for V < 6, 0.25° for V < 7. These remove **81 Gem** (V = 4.9, at (116.53, +18.51), *on* the ridge) and **κ Gem** (V = 3.6, 0.25° east of the spine at Dec 24.4).

## 2. Results

### 2.1 Core detection, Dec 12-20.5, 19 runs, 2002-12-31 to 2007-02-21

Measurement; `verify_precovery_fit.json`, `verify_precovery_stats.json`, `verify_precovery_scan.json`.

| Case | A (R) | χ²/dof | random null | same-shift null | shift scan |
|---|---|---|---|---|---|
| g PSP, gauss028 ([O III] units) | **2.31 ± 0.38** | 1.14 | 6.6σ | **3.9σ** | 4.4σ |
| g FramesSub, gauss028 | 2.39 ± 0.46 | 1.08 | 5.5σ | 2.6σ | 3.1σ |
| g Frames, gauss028 | 2.55 ± 0.53 | 1.28 | 5.5σ | 2.7σ | — |
| g PSP, bigmask | 1.96 ± 0.34 | 0.91 | 6.2σ | 4.6σ | — |
| g PSP, all confounder regressors | 1.65 ± 0.31 | 1.12 | 5.9σ | 3.5σ | — |
| g PSP, window ±1.5° | 1.32 ± 0.32 | 0.72 | 3.8σ | 2.6σ | — |
| g PSP, quality = 3 only | 2.31 ± 0.37 | 1.11 | 6.6σ | 4.1σ | — |
| g PSP, FWHM 0.6° / 1.0° | 2.12 / 2.26 (peak) | 1.3 / 1.6 | 8.4 / 9.1σ | 5.2 / 5.6σ | — |
| g PSP, DESI-shaped template (× DESI) | 0.55 ± 0.16 | 0.93 | 3.4σ | 1.6σ | — |
| g Frames / FramesSub, DESI template | 0.88 ± 0.22 / 0.77 ± 0.19 | 0.7 / 0.9 | 4.0 / 4.1σ | 1.5 / 1.4σ | — |
| **u PSP, [O II] units** | **1.40 ± 0.24** | 1.28 | 6.7σ | 4.4σ | 4.7σ |
| u FramesSub, [O II] units | 1.72 ± 0.31 | 1.14 | 6.0σ | 5.4σ | 4.4σ |
| u PSP, bigmask | 1.37 ± 0.21 | 1.48 | 7.8σ | 6.4σ | — |
| joint g+u (PSP) | — | — | 7.7σ | 4.5σ | — |
| r PSP (R5007-equivalent) | +0.42 ± 2.01 | 1.49 | 0.3σ | 0.2σ | — |
| i PSP | −8.7 ± 15.0 | 0.97 | −0.6σ | −0.4σ | — |
| r / i FramesSub | −1.05 ± 2.24 / −15.9 ± 16.7 | 1.1 / 1.0 | −0.5 / −1.0σ | — | −0.3 / −0.3σ |
| z PSP | +124 ± 33 | 0.82 | (4.7σ) | (2.7σ) | — |

Notes on the table:
- **z** is unstable. Its jackknife error is 98. Dropping run 3614 alone moves it from +124 to +23. I disregard it.
- **Jackknife (g PSP):** jackknife error 0.34. Leave-one-out range +2.16 (without run 3560) to +2.45 (without run 3614). No single run drives the result.
- **Shift scan:** a peak at s = 0 and negative side lobes at ±0.8-1.4°.
  - The side lobes are deeper (−1.2…−1.7 R) than the response expected from the ridge alone (−0.35…−0.6 R).
  - The excess lobe depth is ≈ 2σ. It suggests the real emission is broader or more structured than the 0.28° template. Templates of FWHM 0.6-1.0° give similar peaks and higher significance (table).
  - Off-peak (|s| ≥ 2°) the mean is 0.12 R and the scatter 0.53 R for g PSP. Wave 1 found 0.13 ± 0.25 × NSNS.
- **Significance depends on the null.** Against coherent same-shift nulls the g detection is 2.6-4.6σ for gauss028-type templates (1.4-1.6σ for the DESI-shaped template). Wave 1's "≈ 4σ" is the upper end of this range, not a floor.

### 2.2 Is the core excess [O III] line light? Band SED (derived)

Core ridge amplitudes in 10⁻³ nmgy arcsec⁻² (PSP):
- u 3.33 ± 0.57;
- g 4.16 ± 0.68;
- r 0.29 ± 1.41;
- i −1.68 ± 2.92.

Fit of two spectral shapes to u, g, r and i:
- The DESI emission-line spectrum predicts u : g : r : i = 1.44 : 1 : 0.37 : 0.10. It fits with χ² = 6.3 (3 dof).
- A flat f_ν continuum gives χ² = 9.1. FramesSub gives 4.9 vs 13.0.
- With g fixed, a flat continuum over-predicts r by 2.5σ and i by 1.9σ; FramesSub gives 2.8σ and 2.2σ.
- Scattered starlight or cirrus, which is redder still, is therefore disfavoured at ≈ 3σ. The g excess behaves like emission lines.
- The IRIS 100 µm regressor reduces A by only 10 %.

**u/g is low.** Observed 0.80 ± 0.19 (PSP) or 0.95 ± 0.25 (FramesSub), against 1.44 predicted by the DESI ratios.
- This implies [O II]3727/[O III]5007 ≈ 0.8-1.0 in energy (derived, ±0.25 stat).
- Comparisons: DESI 1.12-1.49 (teammate); BOSS 1.84, 95 % range 1.46-2.37 (verify_sdss).
- Likely biases: [O II] sits on the steep red edge of the u filter; the u PSP sky is the most temporally smoothed estimator (§2.5); a 3727 Å line in a broadband sky estimate carries a larger calibration uncertainty.
- I treat the u result as a **qualitative [O II] precovery**, not a ratio measurement.

### 2.3 North (Dec ≥ 20.5, 9 runs): a broadband confounder, not [O III]

North ridge amplitudes in 10⁻³ nmgy arcsec⁻² (PSP):
- g 0.01 ± 1.47;
- u 6.4 ± 1.3;
- r 10.0 ± 4.4;
- i 15.6 ± 9.5.

What this shows:
- Both SED shapes fail (χ² 13-14).
- Run 2864 (2002-01-09, Dec 24.9) is +5.0σ in u, +4.3σ in r and +4.4σ in i, but 0.1σ in g.
- The night of 2002-01-14 (runs 2887+2888, Dec 23) has g = 5.5 ± 1.4 R, u +5.4σ and r +2.6σ.
- The spectroscopic [O III] in the north is only 0.56-0.9 R (verify_desi / wave 1), so a 5.5 R g excess is not [O III].
- The bright stars κ Gem (0.25° from the spine), Pollux (northern tip) and Castor are the obvious suspects. The enlarged masks leave only 6 north runs and the u/i excess persists, so the source is not identified.
- **Conclusion:** sky-level photometry cannot be used at Dec ≥ 20.5. This agrees with wave 1, which found the north noisy and unconstraining.

### 2.4 Absolute scale: imaging vs spectroscopy (task a, cross-calibration)

I ran noise-free models through exactly the same pipeline (same runs, fields, nuisance model, footprint averaging and gauss028 template); `verify_precovery_extra.json`.

| Comparison | Value | Comment |
|---|---|---|
| Imaging A_obs (g PSP) | 2.31 ± 0.38 R | measurement |
| Expected from the wave-1 DESI forward model | 2.10 R → obs/DESI = **1.10 ± 0.18** | derived |
| With verify_desi's −26 % recalibration (south median 1.6 R) | obs/DESI_rev ≈ 1.49 ± 0.25 | derived; about 2σ tension |
| Expected from the wave-1 SDSS fibre strand model | 0.97 R → obs/fibre = 2.37 ± 0.39 | derived |
| Wave-1 per-fibre `I_OIII_cal`, gauss028, Dec 12-20.5 (my fit, plate-bootstrap error) | 1.12 ± 0.11 R → imaging/fibre = 2.06 ± 0.40 | measurement |
| verify_sdss: wave-1 fibre contrast low by × ≈ 0.55 (residual vs total sky) | fibre_rev ≈ 2.0 R → imaging/fibre_rev ≈ 1.1-1.3 | derived |

- The SDSS imaging (2002-07), BOSS/SDSS total-sky fibres and DESI agree to within a factor ≈ 1.5 at **≈ 1.6-2.5 R**. They independently support verify_sdss's correction factor.
- All three are a factor 4-8 below the paper's 10-13 R on 10′-20′ scales.
- SDSS imaging cannot resolve the arcminute structure. A paper-level 1′ peak would need a filling factor of ≲ 0.25 on 17′ scales.

### 2.5 Estimator response and PSP smoothing (task c)

- **Field-table definitions** (`field_columns.csv`):
  - `sky_b` is the "PSP estimate of sky";
  - `skyFrames_b` is the frames-pipeline sky value;
  - `skyFramesSub_b` is the same after object subtraction.
- **PSP is temporally smoothed** (measurement; `verify_precovery_pspkernel.json`):
  - The high-frequency (second-difference) rms of PSP is 0.31 / 0.36 / 0.46 / 0.59 / 0.75 of Frames in u / g / r / i / z.
  - Regressing PSP on lagged Frames differences gives a kernel spread over about ±2 fields (sum 0.77 in g, 0.46 in u).
  - Part of that spread can come from noise in Frames (Wiener-filter effect), so the kernel is not a clean measurement of PSP's smoothing.
- **Time-smoothed templates:**
  - Smoothing the template along the scan (Gaussian σ = 0.75 / 1.5 / 2.5 fields) improves the PSP fit strongly (Δχ² = −9 / −55 / −107) and raises A_PSP to 2.8 / 4.0 / 6.0 R.
  - For Frames it does not help until σ = 2.5 (Δχ² = +8 / +8 / −21).
  - So the PSP amplitude depends on an unknown temporal response. **The Frames and FramesSub amplitudes, 2.4-2.6 R, are the less model-dependent values**, and they agree with PSP at σ = 0.
- **Frame-level test of run 4458:**
  - Method: I took total sky = calibrated image + the HDU2 global sky model × calibvec, formed 64 × 64 px (25″) block medians with stars masked, and back-fitted common mode, column offsets, camcol drifts and the ridge.
  - The Field-table PSP sky tracks my frame-level field medians with slope 0.72 (Frames 0.87, FramesSub 0.87; correlation 0.83-0.90). Noise in my values dilutes these slopes, so they are lower limits.
  - The single-run ridge fit is **inconclusive**: g gives A = 7.2 R against re-fitted shift nulls of 4.5 ± 3.2 R.
  - A ≈ 1°-wide enhancement centred slightly west of the spine also appears in r and i. In nmgy it has u/g = 0.07 ± 0.03, r/g = 0.72 ± 0.12 and i/g = 1.39 ± 0.71: broadband, not [O III].
  - The first pass labelled the figure "6.9 ± 0.8 (shift nulls)". Those nulls were projection-only (no re-fit) and much too small. **Corrected** in `figures/verify_precovery_frames4458.png`.
  - About 40 % (3.0 of 6.9 R) of the frame-level ridge-correlated g excess sits in the Blanton+2011 global sky model (HDU2). Sky-subtracted SDSS images have therefore partly removed it, as wave 1 found for PS1.
  - In the Field table, run 4458 is +2.6σ in r (8.7 ± 3.3), so this run carries some broadband excess at field level too.
  - The 19-run core result does not rely on run 4458: leave-one-out without 4458 gives +2.35 (Frames).

### 2.6 Confounders (task c)

- **Moon.** Below the horizon (altitude −30° to −80°) for all 19 core runs. It is excluded as a confounder. The nominal correlation of S/N with lunar phase (Spearman 0.48, p = 0.04, one of five tests) is meaningless with the Moon set.
- **Zodiacal light and gegenschein.** The ecliptic is at Dec ≈ +21° here, so β ≈ Dec − 21° for every crossing. Per-run S/N correlates with β (ρ = 0.74, p = 3 × 10⁻⁴), but that is the same thing as Dec: the southern-tip crossings (Dec 12.5-14.4) are weak. The ZL regressors change A by −15 % to +6 % (`conf_zl5` 1.87, `conf_gegen` 1.95 vs 2.20 in full-spectrum units). ZL varies on > 10° scales and cannot make a 17′ ridge that recurs at the spine in runs taken at different times and epochs.
- **Airglow.** Absorbed by the per-field common mode, because the six camcols are simultaneous. Airglow gradients across camcols are absorbed by the per-camcol drifts. The r/i nulls argue against OH, Na and continuum airglow. [O I]5577 varies in time and is removed by S_f.
- **Twilight.** The regressor has no effect on A: `conf_twi` gives exactly the base value.
- **Galactic cirrus.** IRIS 100 µm reduces A by 10 %. Wave 1 and the teammates find no dust or IR counterpart.
- **Bright stars.** 81 Gem lies on the ridge at Dec 18.5. The enlarged masks give 1.96 ± 0.34 R at 4.6σ (same-shift). FramesSub, which subtracts objects, gives the same as PSP.
- **All regressors together.** A = 1.65 ± 0.31 R in full-spectrum units (−25 %), still 5.9σ against random nulls and 3.5σ against same-shift nulls.

### 2.7 Wave-1 code review (`scripts/imaging_optical_sdss_sky.py`)

Error model:
- (i) The formal error is multiplied by √2 "for the pairwise-repeated sky_g values". DR18 `sky_g` is not pairwise-repeated: 0 of 1482 consecutive same-camcol pairs are identical in runs 3614, 4458 and 6597. It is smooth, so the residuals are correlated and some inflation is justified, but the stated reason is wrong.
- (ii) Placebo templates are fitted alone in the same ±2.5° window. At ±0.8-2.0° they absorb the real ridge's negative side lobes, which inflates the placebo scatter. Wave 1 reported a χ²/dof of 0.41.

Both effects make wave 1's errors **conservative**. Wave 1's significance is therefore not inflated by its error model; my independent ≈ 4σ (same-shift) agrees.

Other differences:
- (iii) Wave 1 samples the template at the field centre, not over the footprint. My `centre` variant gives 1.43 vs 1.55 (in-window error model), so this does not matter.
- (iv) Wave 1 uses linear drifts and I use quadratic. My `drift1` variant is consistent.

I found no sign error or unit error. The g response, 1.80 × 10⁻³ nmgy arcsec⁻² per R, agrees exactly.

## 3. Dated precovery table (task d)

Per-epoch significances come from my own fits (`verify_precovery_table.json`, `verify_precovery_extra.json`, `figures/verify_precovery_timeline.png`) unless marked.

| Dataset | Epochs covering the filament | Result per epoch | First ≥ 3σ single epoch |
|---|---|---|---|
| DSS / POSS plates (1950s-1990s) | various | not detected; A < 14 × NSNS (3σ) (wave 1, imaging_optical) | none |
| SHASSA Hα (1997-2000) | — | no Hα (wave 1); not an [O III] test | n/a |
| **SDSS imaging sky, g ([O III])**, core Dec 12-20.5 | 19 runs, 2002-12-31 … 2007-02-21 (seasons 2002-03, 2003-04, 2004-05, 2006-07) | per run −1.2…+3.3σ; 16 of 19 positive; cumulative 3.3σ after the first core run, 3.9σ by 2003-01-25, 4.6σ by 2003-02-01, 6.6σ total (random null) | **2002-12-31, run 3560, Dec 20.2, 3.77 ± 1.16 R (3.3σ); same run r −0.4σ, i −0.9σ** |
| SDSS imaging sky, g, north Dec ≥ 20.5 | 9 runs, 2001-11-12 … 2002-12-06 | 2002-01-14 (runs 2887/2888) 3.9σ in g, with u, r excess; others ≤ 1σ | none attributable (broadband confounder) |
| **SDSS imaging sky, u ([O II])**, core | same 19 runs | 1.40 ± 0.24 R combined | 2004-02-17, run 4469 (3.1σ). Earlier north u peaks (2002-01-09, 5.0σ) are confounded |
| SDSS-I/II sky fibres ([O III]) | plates from 2001-11-27 (754, no on-ridge fibres) and 2002-02-09 (857) … 2008-02-01 | 857-52314 (2002-02-09) 1.8σ; 888-52339 0.5σ; 927-52559 (2002-10-12) 2.7σ | **2002-10-30, plate 927-52577, Dec 23.1: 1.03 ± 0.21 R, 5.0σ** (only 3 on-ridge fibres; wave-1 residual-flux values; not placebo-tested by verify_sdss). Next: 1203-52669 (2003-01-30) 3.5σ, which verify_sdss independently finds at 3.4σ, and 1263-52708 (2003-03-10) 3.9σ |
| BOSS / SDSS-V-era plates | 2010-12-07 … 2020-12-29 | up to 8.1σ (4470-55587) | 2011-01-08 (4502-55569, 4.7σ) |
| LAMOST LRS stellar residuals ([O III]5007) | 77 plan-dates 2011-11-09 … 2024-01-13 (39 with ≥ 5 on/off stars) | up to 9.3σ (2016-01-30) | 2013-03-10 (GAC113N27M1, 3.1σ); 2013-11-06 4.4σ |
| DESI DR1 | 2021-01-03 … 2022 (monthly bins with ≥ 30 fibres) | 5.7-48σ per month (slope on the forward-model prediction; partly circular) | 2021-01 (first night in bin 2021-01-03) |
| NSNS [O III] | DR0.2 data 2022-12-12 … 2025-04-04 (first [O III] camera 2022-12; full array 2024-01) | detected; per-pixel epochs not public | ≤ 2024-01 (self-reported "first clear night" with the full array); undated in 2022-12 … 2023-12 |
| SWONS | 2023-11 (or 12), self-reported | discovery | 2023-11 |
| Amateur public images | AstroBin zh47r9 **2024-01-14** (first public image); Ziegenbalg Cloudy Nights post 2024-06-11 | — | 2024-01-14 |

**Earliest ≥ 3σ:**
- **2002-10-30**: SDSS spectroscopic plate 927, MJD 52577, [O III] at the right wavelength, 5.0σ in my per-plate fit. This rests on 3 on-ridge fibres and is not placebo-validated.
- **2002-12-31**: the first imaging detection (run 3560, 3.3σ).
- **2003-01-30**: the earliest epoch confirmed by two independent analyses (plate 1203).
- In every reading the filament was present about **21 years before** its 2023-11 discovery.
- Data from 2001-11 to 2002-03 (SDSS runs 2711-2888 and plates 754, 857, 888) give no attributable ≥ 3σ detection.

## 4. New findings

1. **[O II] precovery in the SDSS u-band sky.** u PSP gives 1.40 ± 0.24 R and FramesSub 1.72 ± 0.31 R ([O II] doublet) on the core ridge, 4.4-6.4σ against same-shift nulls. The joint g+u detection is 4.5σ against coherent nulls (7.7σ random). The low u/g ratio (§2.2) is a calibration-level tension, not a physical measurement.
2. **The core g excess is line-like.** The r and i controls are null, and a flat or red continuum is disfavoured at ≈ 3σ. That excludes scattered starlight, cirrus and airglow continuum as the source.
3. **The north sky-level signal is a broadband artefact.** It shows in u, r and i without g. It must not be counted as a precovery.
4. **The PSP `sky_b` estimator is temporally smoothed** relative to Frames (high-frequency rms 0.36 in g, 0.31 in u). Amplitudes from PSP depend on the assumed temporal response. Use `skyFrames`/`skyFramesSub` for quantitative work.
5. **Imaging cross-calibration** sides with the total-sky SDSS fibres (verify_sdss) and is 1.1-1.5 × DESI. On 17′ scales the [O III] ridge is ≈ 2-2.5 R from three independent instruments over 2002-2022.
6. **No variability.** At fixed position over 2003.97-2007.14, Δ = +0.8 ± 1.3 R. Taken with DESI, the [O III] has not changed by more than a factor of about 2 over 2002-2022.
7. **Sky-subtracted SDSS images already remove ≈ 40 % of the ridge-correlated sky** into the global sky model. Precovery must use sky levels or the global-sky HDU, not the sky-subtracted images.

## 5. Limitations

- I deliberately avoided NSNS. The morphology is the silk_trace E+W ridge model, so amplitudes are template-dependent:
  - FWHM 0.28-1.0° gives 2.1-2.3 R;
  - the DESI-shaped template gives 0.55-0.88 × DESI, with low coherent-null significance (1.4-1.6σ);
  - the shift-scan side lobes suggest broader or more structured emission.
- Per-run errors rest on 10-20 off-window nulls each. The per-run tail cannot be calibrated independently, because the null z-distribution is |z| < 3 by construction. Single-run significances such as 3.3σ for run 3560 are uncertain by about ±0.5σ.
- Field-level injections are linear and cannot test the PSP algorithm. Only one run (4458) was examined at frame level, and that test is inconclusive. A multi-run frame-level analysis (about 2 GB per run) would be needed to validate the estimator response properly.
- The SDSS fibre per-plate significances use wave-1 per-fibre `I_OIII_cal`. verify_sdss showed these underestimate strand contrast by about 0.55 (residual vs total sky), so per-plate S/N can be biased. Plate 927 rests on 3 on-ridge fibres.
- The DESI per-month significance regresses on wave 1's forward model, which was fitted to all epochs, so it is partly circular. The DESI absolute scale is under revision (verify_desi). LAMOST per-plan values are ON−OFF medians of wave-1 stellar residuals.
- The z band is unstable, and the u-band line response is uncertain at the 10-20 % level.

## 6. Files

- Scripts: `scripts/verify_precovery_{lib,prep,download,fit,frames,frames_sed,stats,scan,table,extra,pspkernel,figs}.py`
- Derived data:
  - `data/derived/verify_precovery_fields.npz`: per-field templates, masks and conditions;
  - `verify_precovery_fit.json` (+ `_fit_log.txt`): all cases, per-run results and nulls;
  - `verify_precovery_stats.json` (+ log): jackknife, seasons, null significance, joint g+u, Dec profile, per-run tables;
  - `verify_precovery_scan.json`;
  - `verify_precovery_table.json` (+ log): dated table;
  - `verify_precovery_extra.json`: cross-calibration, SED, null calibration, temporal, confounders, fibres, per-night;
  - `verify_precovery_pspkernel.json`;
  - `verify_precovery_frames.json`, `verify_precovery_frames_sed.json` (+ logs), `verify_precovery_frames_blocks*.npz`, `verify_precovery_frames_resid_{g,i}.npz`, `verify_precovery_frames_fieldcompare_{g,i}.npy`.
- Raw data: `data/raw/verify_precovery_sdss/` (Field table CSV and query JSON, the column list, filter curves, and `frames_4458/` with its URL lists and download log).
- Figures:
  - `figures/verify_precovery_shiftscan.png`: combined amplitude vs RA shift for g, u, r and i;
  - `figures/verify_precovery_dec.png`: per-run amplitude vs crossing Dec, compared with DESI and fibres (both at wave-1 calibration) and the paper;
  - `figures/verify_precovery_frames4458.png`: frame-level maps and profiles for run 4458, with corrected null labels;
  - `figures/verify_precovery_timeline.png`: per-epoch significance against date for all datasets.
