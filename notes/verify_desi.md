# Wave-2 verification: DESI DR1 residual measurements (tag `desi`)

Date: 2026-11-16. Verifier tag `desi`. Velocities are heliocentric (= barycentric to < 0.02 km/s);
V_LSR = V_helio - 12.1 km/s. Ratios are energy ratios to [O III]5007 unless stated.

Inputs:
* the wave-1 DESI products listed in the task (5,220 SPARCL coadds);
* **2,400 new coadds** fetched for this check (1,225 near P3 / Dec 12-15.5, plus 1,175 controls more than
  5 deg from every strand).

Everything was re-measured with my own line fitter, template, statistics and calibration.

**Revision 2 (resume, 22:40-23:30 JST).** A second-pass audit (section "Second-pass audit") re-measured every line
with a different method: box integrals of the residual with an empirical sky mask, in place of Gaussian fits.
It calibrated the bootstrap with 11 line-free null boxes and the [O I]5577 sky line.
* Three first-pass claims were overstated and are corrected below:
  - [O III]/Hb ~ 47 is now a lower limit, > 11-18 (2 sigma);
  - R23 ~ 130 is now > ~30;
  - [O III]4363 at "3-5 sigma" is now about 2-3 sigma.
* The audit confirmed:
  - the calibration revision (1.62 R south by both methods);
  - [O II]/[O III], He II and the velocities;
  - that Hg cannot produce the 4363 feature.

## Bottom line

1. **Detection: confirmed and robust.** A geometry-only matched filter gives S/N = 154 (107 after |chi| < 8
   clipping). The empirical nulls are all far below that:
   * rigid rotations/translations of the trace that do not overlap it (|r| < 0.05): std 4.7, max 27.5 (n = 227);
   * label shuffles within a tile-petal: 12 +- 4.5, max 27.7;
   * 17 fake lines: |S/N| <= 2.4.

   The ridge is found:
   * in every Dec band (11.6-132 sigma);
   * in dark and bright, main and SV;
   * in QSOs alone (55 sigma);
   * in the 2,328 new fibres alone (125 sigma);
   * in every tile-dropping jackknife (min S/N 81) and every HEALPix jackknife (min S/N 82).

   The revision-2 box-integral fluxes give S/N 106 with the same template.

   Formal errors understate the null scatter by about 3-5x, so wave 1's "> 60 sigma" is generous as a number.
   The detection is still beyond doubt: about 30 sigma against the empirical nulls.
2. **Absolute brightness: REVISED down by about 26 %, because of a calibration error in wave 1.**
   * Wave 1 converted calibrated flux to surface brightness with `MEAN_PSF_TO_FIBER_SPECFLUX` (median 0.77). That
     factor is the Legacy-Surveys "fibre-flux" normalisation; it is the constant 0.789 for point-source targets.
   * The correct factor is the true point-source fibre acceptance of each exposure, from the GFA:
     median 0.566 at 5008 A (dark 0.6, bright 0.45).
   * Per fibre, the wave-1 surface brightness is too high by 1/0.736 (median; 1.13-1.82).
   * Revised ridge peaks (forward model):
     - south median **1.6 R** (1.60-1.71 over the calibration variants; wave 1 2.4);
     - E strand 1.4-2.3 R at Dec 14-18;
     - W strand 1.2-1.6 R, plus 3.9 +- 0.2 R at Dec 18 near P2;
     - north median **0.56 R** (wave-1 definition; wave 1 0.86).
   * The model-free lower limit of 3.05 +- 0.14 R becomes **1.31 +- 0.05 R** in the same bin. The brightest bin
     anywhere is 1.81 +- 0.05 R at (116.59, +15.48), so SB >= 1.7 R (2 sigma).
   * Revision 2 check with box fluxes and one global matched-filter amplitude: **1.62 R south, 0.71 R north**
     (2.35 / 0.97 R with wave 1's conversion). The corrected SDSS re-analysis (`notes/verify_sdss.md`) gives
     1.4-2.5 R south.
   * The paper's 10-13 R is 4-8x above every DESI ridge value.
3. **Line ratios** (revision 2 values; the first-pass values are given for reference).
   The filament is Balmer-dark, but Hb is not detected well enough to give a central [O III]/Hb.
   * **Hb/[O III] = 0.012-0.021, +- 0.018 (Gaussian) to +- 0.04 (box) after null-calibrated inflation.**
     So **[O III]/Hb > 18 (Gaussian) / > 11 (box) at 2 sigma**. Wave 1 had 20 (> 11); that lower limit is
     **confirmed**. The first pass reported 47 (> 24), which used the bootstrap uninflated and is withdrawn.
   * [O II]3727/[O III] = **1.37-1.41** (68 %: 1.21-2.07); 1.11-1.16 (+- 0.08) at Dec < 19.5; 4.2-4.4 in the
     north. Confirmed by both methods. The blue end is OK after a -6 % chromatic-acceptance correction.
   * R23 > **~30** (2 sigma: 29 box, 34 Gaussian). Wave 1 had 58 (> 35). The first pass had about 130, which is
     withdrawn.
   * He II/[O III] < 0.026 (97.5 %). Confirmed.
   * [O III]4363/5007 = **0.06-0.09**: 0.083 Gaussian IV, 0.091 box IV, 0.057 pixel regression.
     * Significance is **about 2-3 sigma** after null calibration (1.8 sigma box, 2.9 sigma Gaussian). One null box
       at 4381.6 A gives 0.048 +- 0.014, so this region produces spurious ratios of about 0.05.
     * The feature is centred at 4364.55-4364.73 A, against 4364.62 A expected for [O III] and 4359.55 A for Hg I.
       The box Hg coefficient is 0.004 +- 0.019. **Hg is excluded as the cause.**
     * T_e([O III]) ~ (3-5) x 10^4 K is favoured but not secure.
   * [N II]/Ha ~ 0.5 and [S II]/Ha ~ 0.4 (wave 1): **refuted**. Direct differential values:
     * [N II]6583/Ha = 1.05-1.21 (box 1.05; Gaussian 1.21, 68 % 1.12-1.43);
     * [S II]/Ha = 0.70-0.75;
     * south: 1.0-1.08 and 0.51-0.52.

     Ha is the least reliable line. In the box method Ha is positive in Z > 0 fibres (0.115) and near zero in
     Z < 0 fibres (0.001 +- 0.09). Wave 1's sign split found the same asymmetry, but here it is only about 1 sigma.
4. **Velocity.**
   * [O III] V_helio = **+11 to +12 km/s** overall: IV regression +11.0 (10.3-14.3); positive fibres +12.2 +- 0.8.
     This confirms wave 1's +14.5 within the errors.
   * Near the LAMOST regions DESI is 12-15 km/s below LAMOST: P1 +11.9 +- 1.9, P2 +15.2 +- 2.9, P3 +9.2 +- 1.0,
     against LAMOST +26.8, +26.8, +23.8.
   * The E/W strand difference (wave 1: +13.0 vs +18.5) is **not reproduced**.
     * Positive fibres: E +15.0 +- 1.0, W +13.7 +- 1.5.
     * Revision 2 joint E/W pixel regression on both [O III] lines: E +12.0 (11.3-14.3), W +11.4 (8.6-15.6).
     * A 3-6 km/s E < W difference appears only at Dec 15.5-17.5, at <= 2 sigma.
   * **New (revision 2):** Dec >= 19.5 gives +20.9 (18.5-23.1) km/s, against +10.6 (10.2-12.0) at Dec < 19.5. The
     northern half is about 10 km/s redder, partly reproducing wave 1's Dec-band trend.
   * The spectra are barycentric: night-sky lines move with v_corr at slope 0.98 after estimator calibration.
   * No [O III] coadd uses the -19 km/s backup tiles 41485/41487 (they lie at Dec +5 to +6.5).

## Results table

| # | Wave-1 claim | Verdict | My value (68 % unless stated) |
|---|---|---|---|
| a | Ridge detection, Delta chi^2 = 8,903, "> 60 sigma", filament-specific | **CONFIRMED** (significance overstated) | S/N 154 vs null std 3.2-4.7 (non-overlapping rigid rotations) -> 30-48 sigma; shuffle z = 31.5; fake lines \|S/N\| <= 2.4; all jackknives >= 81 |
| b1 | Ridge peak south median 2.4 R | **REVISED** | 1.62 R (variants 1.60-1.71; +-8 % sys). With wave-1 conversion my model gives 2.50 R, so the model agrees and the calibration differs |
| b2 | North median 0.86 R | **REVISED** | 0.56 R (same node definition; 0.72 R with wave-1 conversion) |
| b3 | Strict lower limit 3.05 +- 0.14 R at (116.27, +17.28); >= 2.8 R (2 sigma) | **REVISED** | Same bin 1.31 +- 0.05 R (median 1.37). Best bin 1.81 +- 0.05 R at (116.59, +15.48) -> >= 1.7 R (2 sigma). Look-elsewhere: fake-line max bin S/N <= 5.2 vs 35 |
| b4 | 1 R = 7.43e-18 erg s^-1 cm^-2 arcsec^-2; 1.52" fibre | **CONFIRMED** | 7.419e-18; plate-scale fit gives fibre-weighted mean 1.520" (1.585" centre to 1.454" edge) |
| b5 | Conversion uses PSF_TO_FIBER_SPECFLUX (+-15 % sys) | **REFUTED** (wrong factor) | True acceptance / PSF_TO_FIBER = 0.736 median (0.55-0.88), 0.57 for bright-program coadds |
| c1 | [O III]/Hb = 20 (15-26), > 11 (97.5 %) | **CONFIRMED as a lower limit**; central value not determined (rev. 2) | Hb/[O III] = 0.021 +- 0.018 (Gaussian, x1.5 inflated) / 0.012 +- 0.041 (box, x1.68) -> [O III]/Hb > 18 / > 11 (2 sigma). First-pass "47 (> 24)" withdrawn (uninflated bootstrap; Hb only ~1-2 sigma) |
| c2 | [O II]3727/[O III] = 1.49; 1.12 at Dec < 19 | **CONFIRMED** | 1.37 (1.21-2.07); south 1.11 (1.07-1.18); north 4.4 (3.5-5.2); QSO/z > 0.5 1.44. Achromatic would be 6 % higher |
| c3 | R23 ~ 58 (> 35) | **CONFIRMED as a lower limit**, slightly weaker (rev. 2) | R23 > 34 (Gaussian) / > 29 (box) at 2 sigma; central value not determined. First-pass "128 (> 84)" withdrawn |
| c4 | He II/[O III] < 0.034 | **CONFIRMED** | -0.008 (-0.024...0.000), < 0.026 (97.5 %) |
| c5 | [O III]4363/5007 = 0.056 (0.022-0.077), T_e ~ 3e4 K (~2 sigma) | **CONFIRMED as a 2-3 sigma hint; not Hg** (rev. 2) | 0.083 Gaussian (+-0.029 inflated), 0.091 box (+-0.051), 0.057 pixel regression (16-84 %: 0.015-0.069). Centre 4364.55-4364.73 vs 4364.62 [O III] and 4359.55 Hg. Hg coefficient 0.004 +- 0.019 (box), 0.009 +- 0.015 (Gaussian). T_e ~ (3-5)e4 K favoured, not secure |
| c6 | [Ne III]/[O III] = 0.10 | **CONFIRMED** | 0.099 (0.057-0.109) |
| c7 | [S III]9069/[O III] = 0.031 | **CONFIRMED** | 0.028 (0.024-0.048), chromatic-corrected (+12 %) |
| c8 | [O II]3729/3726 = 1.30 (1.24-1.39), n_e ~ 110 | **CONFIRMED**, wider systematic range | 1.30 (1.24-1.33). Subsets 1.17 (QSO/z > 0.5) to 1.37 (neighbour instrument), so n_e ~ 40-300 cm^-3 |
| c9 | Filament [N II]/Ha ~ 0.5-0.7, [S II]/Ha ~ 0.3-0.55 (Ha biased, Ha/Hb = 1.9) | **REFUTED** (systematics-limited) | Direct: 1.21 (1.12-1.43) Gaussian, 1.05 box; [S II]/Ha 0.75 / 0.70. South 1.0-1.08 / 0.51. Z > 0 fibres only: 1.08 (box). Ha/Hb undetermined (Hb ~1-2 sigma). The wave-1 values came from Hb rescaling |
| c10 | [O III]5007/4959 = 3.11 | **CONFIRMED** | 3.04 (2.94-3.10) |
| c11 | [O II]7325/[O III] = 0.037; [Ar III]/[O III] = 0.018 | **CONFIRMED** | 0.036 (0.033-0.042); 0.013 (0.011-0.020) |
| d1 | [O III] V_helio +14.5 (13-16) | **CONFIRMED** (slightly lower) | +11.0 (10.3-14.3) regression; +12.2 +- 0.8 positive fibres; rev. 2: +11.3 (10.2-13.9), both [O III] lines; south +10.6 (10.2-12.0), north +20.9 (18.5-23.1) |
| d2 | E strand +13.0, W +18.5 | **REFUTED** (no significant difference) | Positive fibres: E +15.0 +- 1.0, W +13.7 +- 1.5. Rev. 2 joint regression: E +12.0 (11.3-14.3), W +11.4 (8.6-15.6). Only Dec 15.5-17.5 shows E < W by 3-6 km/s (<= 2 sigma) |
| d3 | Barycentric frame; v-vs-v_corr slope +0.20 +- 0.22 | **CONFIRMED** | Sky-line test: slope 0.869 +- 0.005, i.e. 0.98 after estimator calibration (Hg 4358: 0.98). [O III] slope +0.49 +- 0.13, or +0.30 +- 0.15 with Dec fixed effects (position-confounded) |
| d4 | No [O III] fibre from tiles 41485/41487 | **CONFIRMED** | 0 of 7,345 usable coadds. The tiles are main/backup at Dec +6.5 / +5.1 |
| d5 | LAMOST [O III] higher than DESI by 13-16 (P1), 8 (P3) km/s | **REVISED** (larger at P3) | P1 -14.9 +- 4.7, P2 -11.6 +- 5.0, P3 -14.6 +- 3.4 km/s (DESI minus LAMOST) |
| e | Coverage near P3 / Dec 12-15 sparse; no > 5 deg controls | Done | +1,225 near and +1,175 controls. Dec < 15.5 now has 2,207 fibres. DESI has no coverage at Dec 12-13.7 within 0.5 deg of the trace. Controls: -0.007 +- 0.008 R, chi std 1.09 |

## (b) How DESI calibrates a fibre, and what that means for a uniform source

Sources read as text, never executed:
* `data/raw/desispec_src/desispec_fluxcalibration.py.txt` (the wave-1 copy);
* fetched 2026-11-16 from github.com/desihub (main branch), saved in `data/raw/verify_desi_src/` with URLs in
  `fetch_date.txt`: `desispec_fiberfluxcorr.py.txt`, `desimodel_fastfiberacceptance.py.txt`,
  `desispec_coaddition.py.txt`.

1. **Flat-fielding.** The fibre flat makes a uniform surface brightness give equal counts in every fibre. It
   therefore contains the fibre solid angle, which shrinks towards the field edge (see the comment in
   `flat_to_psf_flux_correction`).
2. **Point-source correction.** `fluxcalibration.py` (lines 1063-1083) multiplies the flat-fielded flux by
   FLAT_TO_PSF = 1/(f_ps,i x platescale_i^2), median-normalised. Here f_ps,i is the point-source acceptance from
   `FastFiberAcceptance` at the exposure seeing and the fibre's offset.
3. **Calibration vector.** C is fitted on standard stars, so it absorbs their real absolute acceptance. The
   calibrated flux is `flux / C` (line 1781). The comment at line 1170 reads: "correction for point sources already
   applied to the flux vector".
4. **Consequence.** A point source receives its total flux. A uniform surface brightness receives
   **F_cal = SB x Omega_i / f_ps,i**, so **SB = F_cal x f_ps / Omega**. The docstring of
   `psf_to_fiber_flux_correction` says the same:
   "psf calibrated flux ~= current_fiber_frac * total_flux / current_fiber_frac_point_source".
   Ignoring f_ps would overestimate a diffuse source by about 1.7x.
5. **PSF_TO_FIBER_SPECFLUX is not f_ps.** It converts PSF-calibrated flux to the Legacy-Surveys "fibre flux" of a
   1.5" fibre in 1" Gaussian seeing, scaled so that a point source gets exactly 0.789. In practice:
   * it is 0.789 for PSF-morphology targets (91 % of the QSOs here have 0.789 +- 0.002);
   * it is 0.789 x f_ps/f_gal x f_gal,nom/f_ps,nom for galaxies (median 0.77);
   * it carries almost no information about the real seeing.

   **This is the wave-1 bug.**
6. **DESI's own sky-brightness code.** `skymag.py` divides the calibration by the real acceptance `ffracflux_wave`
   (default 0.6, DESI-6043) and a 1.52" fibre. That is the same formula as item 4.
7. **The real acceptance per exposure** is `desi_dr1.exposure.fiber_fracflux_gfa`
   (Data Lab TAP, `data/raw/verify_desi_tap/exposure_gfa_box.csv`, script `verify_desi_tap.py`).
   * 196 of the 238 exposures have it: median 0.595 (dark 0.637, bright 0.446).
   * By seeing: 0.745 at 0.8", 0.555 at 1.1", 0.36 at 1.6".
   * The 42 exposures with seeing = 0 were filled with their program median.
   * A Moffat(beta = 3.5) profile in a 0.76" radius reproduces the GFA values to 3 % (scale 0.968).
8. **Per-coadd value.** f_ps is the EFFTIME-weighted *harmonic* mean over the coadd's exposures (correct for
   averaging F_e = SB Omega/f_e). It is moved from the GFA band to the line wavelength with the Moffat model, with
   the whole FWHM scaling as lambda^-0.2.
   * Median f_ps(5008 A) = 0.566 (16-84 %: 0.42-0.68).
   * Ratio to wave 1's factor: 0.736 median (0.55-0.88). Dark coadds 0.77, bright 0.57.
   * Wave 1's own systematic allowance of +-15 % did not cover this, because the bias is one-sided.
9. **Fibre angular size** (`verify_desi_platescale.py`). I fitted the 716,386 fiberassign focal-plane radii against
   the targets' angular offsets from their tile centres (odd polynomial, rms 14 um).
   * Plate scale: 67.5 um/" at the centre; at the edge 76.8 um/" radial and 70.5 um/" azimuthal.
   * The 107 um core subtends 1.585" at the centre and 1.454" at the edge.
   * Fibre-weighted mean 1.520": the 1.52" mean is **confirmed**. Solid-angle rms 4.7 %, which averages out.
10. **Systematics of my conversion** (in `verify_desi_brightness.json`):
    * chromatic model (full / half / none): 1.62 / 1.66 / 1.71 R south median;
    * GFA effective wavelength 6100-6700 A: +-1 %;
    * fibre positioning offsets (zpix median 4.7 um = 0.07"): <= 2 % lower;
    * total about +-8 %.

    For line ratios the chromatic acceptance matters:
    * f(3727)/f(5008) = 0.94 (half-model 0.97);
    * f(6563)/f(5008) = 1.06;
    * f(9069)/f(5008) = 1.12.

### Brightness results (`verify_desi_brightness.py`, my own forward model)

I used unit-peak Gaussian ridges on E/W/W2 (trace HWHM) with triangular 1-deg Dec nodes, minus the EFFTIME-weighted
mean over all fiberassign positions of each tile-petal, fitted by WLS with 4-sigma clipping.
* GFA conversion: chi^2/dof = 10522/7204. Wave-1 conversion: 10424/7207.
* Fitted with wave 1's conversion, my model reproduces the wave-1 nodes to about 10 %. For example E@16 is 2.60 vs
  2.52, E@17 2.50 vs 2.53 and W@18 4.63 vs 4.98, so the downward revision is purely the calibration.

| Node | GFA conversion (R) | my model, wave-1 conversion | wave 1 |
|---|---|---|---|
| E@14 | 2.30 +- 0.05 | 3.42 +- 0.07 | 4.00 +- 0.16 |
| E@15 | 1.44 +- 0.05 | 2.01 +- 0.07 | 1.69 +- 0.30 |
| E@16 | 1.70 +- 0.06 | 2.60 +- 0.07 | 2.52 +- 0.08 |
| E@17 | 1.90 +- 0.08 | 2.50 +- 0.12 | 2.53 +- 0.12 |
| E@18 | 2.26 +- 0.17 | 2.63 +- 0.37 | 2.24 +- 0.38 |
| E@20-22 | 1.18-1.24 (+-0.11) | 1.71-2.07 | 1.72-2.14 |
| E@25-27 | 0.56-0.94 (+-0.12) | 0.72-1.28 | 0.86-1.19 |
| W@15-17 | 1.15-1.62 (+-0.06) | 1.91-2.21 | 1.59-2.42 |
| W@18 | 3.88 +- 0.20 | 4.63 +- 0.36 | 4.98 +- 0.38 |
| W2@22-23 | 0.23-0.33 (+-0.10) | 0.32-0.48 | 0.32-0.40 |

* Model sky SB at the paper's regions (GFA conversion): P1 2.1 R, P2 3.0 R, P3 1.5 R, P4 0.
* Measured fibre residuals within 0.1 deg (lower limits): P1 1.18 +- 0.05 R (n = 17), P2 1.33 +- 0.44 (n = 6),
  P3 0.74 +- 0.04 (n = 23).
* Controls more than 5 deg from the trace (n = 1,148): mean residual -0.007 +- 0.008 R, chi std 1.09, 0.35 % above
  +3 sigma. There is no zero-point offset and no false-positive structure far away.
* The binned lower limit has a winner's-curse risk. The same maximum-bin search on 17 fake lines gives a median
  maximum S/N of 2.6 and a largest of 5.2, against 35 for [O III]. Selection bias is negligible.

## Methods (own code; wave-1 scripts read for bugs, not reused)

* **Data.** The wave-1 windows plus my fetch (`verify_desi_fetch.py`).
  * The fetcher reuses the wave-1 SPARCL worker but writes to `data/derived/verify_desi_windows/`, so the wave-1
    globs are unaffected.
  * Endpoint: `https://astrosparcl.datalab.noirlab.edu/api/spectras/`, DESI-DR1, 2026-11-16. Transfer 0.47 GB,
    logged per batch in `ext_log.jsonl` and `ext_fetch_stdout.log`.
  * Queue `verify_desi_sel_extend.csv` interleaves two sets:
    - "near": GALAXY/QSO coadds with ZWARN = 0 at Dec 12-15.5 within 1 deg of a strand, or within 0.7 deg of P3,
      sorted by distance;
    - "ctrl": random coadds more than 5 deg from all strands.
  * I stopped after about 45 min with 2,400 of 3,005 fetched: 1,225 of 1,805 near and 1,175 of 1,200 controls. The
    unfetched near coadds were the farthest ones, 0.6-1 deg from the trace.
* **Line fitter** (`verify_desi_linefit.py`).
  * Residual = flux - Redrock model. Every line group is fitted *simultaneously* with a quadratic continuum over
    +-30 A, so continuum errors propagate.
  * Masks: other nebular lines +-3 A; sky lines +-2.5 A (Hg I 4047/4358/5461/5770/5791, Na D, [O I]5577, at
    topocentric wavelength x (1 + v_corr/c)); the object's own redshifted lines +-8 A. One 5-sigma clip.
  * Profile: sigma^2 = LSF^2 + (0.8 A)^2/12 + (12 km/s)^2. The LSF is my median of the SPARCL `wave_sigma` arrays of
    430 pilot spectra.
  * Nebular lines are fitted at V_NEB = +14 km/s. Air-to-vacuum uses Edlen/IAU; the values agree with wave 1 to
    <= 0.005 A.
  * The [O III] doublet is tied at 1:2.98. There is a free-velocity scan for [O III].
  * Hg I 4358 is fitted jointly with [O III]4363 (and [Fe II]4359).
  * 17 fake lines give an MAD of amp/err of 0.98-1.04 per fibre.
* **Per-coadd metadata** (`verify_desi_setup.py`):
  * tile-petal contributions and EFFTIME weights (`desi_dr1.fiberassign` / `exposure`);
  * EFFTIME-weighted barycentric correction (astropy, KPNO, every exposure);
  * f_ps(lambda) and k_R(lambda).
* **Template** (`verify_desi_common.Template`). Unit-peak Gaussians about densely resampled E/W/W2 polylines,
  using nearest-sample distance on the sphere, so any rotation works.
  * Z_i = T(x_i) - sum_t w_t <T>_(tile, petal).
  * The petal mean is over all fiberassign positions, because sky-fibre positions are not in
    `desi_dr1.fiberassign`: only FA_TYPE 1/2/3 rows exist, which confirms the wave-1 statement.
  * No data enter Z.

## (a) Detection (`verify_desi_detect.py`, `verify_desi_detect.json`, `figures/verify_desi_detection.png`)

Statistic: matched filter A = sum w Z O / sum w Z^2, with S/N = A sqrt(sum w Z^2). The per-fibre errors are scaled
by the fake-line MAD (1.01).

| Test | Result |
|---|---|
| True trace (7,345 fibres) | S/N 153.9; 106.9 with \|chi\| < 8 |
| Rigid random rotations + translations (<= 4 deg), 300 draws | Non-overlapping (\|r\| < 0.05, n = 227): mean ~0, std 4.7, max 27.5. \|r\| < 0.02 (n = 126): std 3.2. None of 284 with \|r\| < 0.2 reaches the true value |
| Shuffle O within dominant tile-petal (1,000) | 12.0 +- 4.5, max 27.7. True value is 31.5 sigma above. The non-zero mean is the petal-level part of the signal; the excess is the within-petal spatial match |
| 17 fake lines, true template | S/N -2.35 ... +2.09 |
| RA shifts of the trace | +-0.5 deg: -74 / -26 (anti-correlated troughs); +-1 deg: -38 / -13; +-2-4 deg: -3.8 ... +2.6 |
| Splits | Dec < 15.5: 132; 15.5-17.5: 80; 17.5-19.5: 21; 19.5-22: 24; 22-25: 12; >= 25: 17; dark 140; bright 65; main 50; SV 149; QSO 55; z > 0.5 or QSO 134; new fibres only 125 |
| Jackknives | HEALPix nside-4 min 82; drop-one-tile min 81 (tile 80683); drop-one-tile-petal min 133 |

The detection is specific to the mapped trace and is not driven by any tile, petal, region, target class or
program. The formal S/N overstates significance by the null-scatter factor 3-5 (correlated per-petal residuals and
heavy tails), so the honest significance is about 30-48 sigma.

## (c) Line ratios (`verify_desi_ratios.py`, `verify_desi_regspec.py`; `figures/verify_desi_ratios.png`, `verify_desi_regspec.png`)

First-pass (Gaussian-fitter) results. The [O III]/Hb, R23 and 4363-significance statements in this section are
**superseded by the second-pass audit below**: the bootstrap intervals here are not inflated, and null boxes show
they need x1.5-1.8.

**Estimator.** I used an instrumental-variable slope:

rho_L = sum u (Zf - <Zf>)(c L - <cL>) / sum u (Zf - <Zf>)(O - <O>)

* Zf is the geometry template in flux units and c = f(lambda_L)/f(5008) is the chromatic acceptance.
* u is the inverse variance of the line.
* There is an intercept, or petal fixed effects.
* This differs from wave 1's per-pixel regression on a data-derived instrument.

**Errors.** A HEALPix nside-32 block bootstrap (400 resamples), joint across lines.

**Checks.**
* The 17 fake lines scatter 1.5x their bootstrap sigma (1.8x near 4363), so bootstrap intervals should be inflated
  by about 1.5.
* Variants: petal fixed effects; a neighbour instrument (mean [O III] of fibres in *other* tile-petals within
  0.15 deg); QSO or z > 0.5 targets only; smooth continua; dark only; south/north; Ha patch removed; wave-1 fibres
  only; nside-64 blocks.

| Line / [O III]5007 | all | south (Dec < 19.5) | petal FE | QSO or z > 0.5 | neighbour instr. |
|---|---|---|---|---|---|
| Hb | 0.021 (0.004-0.028) | 0.022 (0.010-0.028) | 0.017 | 0.028 | 0.021 |
| [O III]4959 | 0.332 (0.326-0.342) | 0.334 | 0.333 | 0.340 | 0.331 |
| [O III]4363 (joint Hg, [Fe II]) | 0.083 (0.059-0.091) | 0.092 (0.081-0.096) | 0.083 | 0.078 | 0.098 |
| Hg I 4358 sky (joint) | 0.009 (-0.005-0.019) | 0.009 | 0.010 | 0.008 | 0.005 |
| He II 4686 | -0.008 (< 0.026 at 97.5 %) | -0.008 | -0.011 | -0.009 | -0.007 |
| [O II]3726 + 3729 | 1.37 (1.21-2.07) | 1.11 (1.07-1.18) | 1.27 | 1.44 | 1.21 |
| [Ne III]3869 | 0.099 (0.057-0.109) | 0.111 | 0.099 | 0.096 | 0.106 |
| Ha | 0.082 (0.069-0.101) | 0.073 | 0.058 | 0.077 | 0.096 |
| [N II]6583 | 0.100 (0.087-0.140) | 0.079 | 0.089 | 0.096 | 0.099 |
| [S II]6716 + 6731 | 0.062 | 0.037 | 0.047 | 0.067 | 0.058 |
| [S III]9069 | 0.028 (0.024-0.048) | 0.022 | 0.027 | 0.028 | 0.025 |
| [Ar III]7136 | 0.013 | 0.011 | 0.013 | 0.011 | 0.010 |
| [O II]7320 + 7330 | 0.036 | 0.032 | 0.035 | 0.035 | 0.037 |
| [N II]5755 | 0.014 (0.006-0.026) | 0.010 | 0.015 | 0.018 | 0.012 |
| [Ca II]7291 | -0.003 (< 0.004 at 97.5 %) | -0.003 | -0.004 | -0.005 | -0.003 |
| [Fe II]7155 | -0.004 (< 0.0014 at 97.5 %) | -0.004 | -0.003 | -0.006 | 0.001 |
| [N I]5200 | -0.004 (< 0.043) | -0.006 | -0.004 | -0.005 | 0.000 |
| He I 5876 | -0.005 (< 0.005) | -0.004 | -0.006 | -0.007 | -0.006 |
| [O I]5577 sky (control) | 0.09 (0.04-0.15; 95 % -0.02...0.47) | 0.09 | 0.10 | 0.06 | 0.10 |

Derived values (first pass, uninflated; see the audit):
* **[O III]/Hb** = 47 (all), 46 (south), 60 (petal FE), 35 (QSO/z > 0.5), 47 (neighbour).
  97.5 % lower limits: 24, 28, 21, 18, 26.
* **R23** = 128 (all; > 84 at 97.5 %) and 112 (south; > 69).
* **[O II]3729/3726** = 1.30 (1.24-1.33).
* **5007/4959** = 3.04.
* **Ha/Hb** = 3.9 (2.2-6.1).
* **[N II]/Ha** = 1.21 (1.12-1.43) and **[S II]/Ha** = 0.75 (0.60-1.22) all. South: 1.08 (0.93-1.17) and 0.51
  (0.42-0.65).

**Hb bias tests.** The objects are extragalactic, so stellar Hb absorption at z = 0 does not apply. Redrock mismatch
is random in observed wavelength. Restricting to QSOs and z > 0.5 galaxies (rest-frame Hb window < 3240 A) gives
Hb = 0.028 (0.014-0.040). Smooth continua give 0.012. Both are consistent. Hb lies in a sky-variance excess:
1/ivar is 17-24 % above its median at 4862-4863 A. The excess is probably diffuse Galactic Hb in the sky fibres. It
adds noise, not bias, because the residuals are differential.

**[O III]4363 vs the Kitt Peak Hg I 4358 line.** These are the CRITICAL checks.
1. The Hg line **is present** in the DESI sky. The median 1/ivar of the coadds rises by +50 % at 4359.2-4360.0 A
   (Hg I 4358.33 air = 4359.55 vac, shifted by v_corr). The stacked variance profiles move with v_corr at slope 1.09
   (see (d)).
2. In the pixel-level IV regression spectrum, the 4363 feature is centred at **4364.55 A (4364.15-4364.85)**
   (south: 4364.70, 4364.40-4364.95). [O III]4363 at +14 km/s is expected at 4364.64 A; Hg would be at 4359.52 A,
   5.1 A away. That is 7 LSF sigma, so the line profiles do not overlap.
3. Fitted jointly, the Hg coefficient is 0.011 +- 0.025 (pixel fit) and 0.009 +- 0.015 (per-fibre IV). Hg sky
   residuals do not correlate with the filament instrument.
4. Masking or modelling Hg changes 4363 by < 0.005: alone 0.079; with Hg 0.083; with Hg and [Fe II] 0.083.
5. Six fake lines at 4330-4395 A give ratios with rms 0.015, so 4363 is about 5 sigma above that empirical scatter.
   It is about 3 sigma with the 1.8x inflated bootstrap. The bootstrap 95 % interval reaches 0: a few blocks carry
   most of the weight (3 % of resamples <= 0).

Conclusion: the 4363 signal is at the [O III] wavelength, not Hg. The first pass put its significance at 3-5 sigma;
the audit (null boxes, box method, pixel regression) puts it at **about 2-3 sigma**.
* T_e([O III]): R = (5007 + 4959)/4363 = 16.0 gives about 4.6e4 K (68 %: 3.1-5.3e4 K) for all fibres and 5.4e4 K in
  the south. I used the low-density O&F (2006) eq. 5.4, which is only approximate above 3e4 K.
* So T_e([O III]) ~ (3-6) x 10^4 K: hotter than photoionised gas. This supports collisional excitation and
  confirms the wave-1 hint.

**[O II]3727 at the blue end.**
* The doublet is cleanly resolved in the regression spectrum (`verify_desi_regspec.png`). There is no sky feature in
  1/ivar at 3725-3732 A; the nearest excess is +30 % at 3700-3708 A, outside the line masks.
* The chromatic fibre acceptance lowers [O II]/[O III] by 6 % (full model) or 3 % (half model) relative to
  wave 1's achromatic conversion.
* The ratio is stable across QSO/z > 0.5, dark, petal FE and neighbour subsets (1.21-1.44). It rises steeply to the
  north (4.4), as wave 1 found (6.2).
* Reddening (E(B-V) ~ 0.03-0.05) would raise it by about 5 %; I did not correct for it.

**Red lines.** Wave 1 discarded its direct differential Ha because Ha/Hb = 1.9. Its own all-fibre [N II]/Ha was
1.08 (0.94-1.26) and [S II]/Ha 0.83. It then rescaled with Hb (0.049 of [O III]), which gave 0.5-0.7 and 0.3-0.55.
* My Hb is 2.3x weaker, so Ha/Hb = 3.9 is plausible and the direct ratios stand. Per the audit, Hb is only a
  1-2 sigma signal, so Ha/Hb is undetermined; the direct [N II]/Ha ~ 1.0-1.2 does not depend on Hb.
* The filament's [N II]/Ha ~ 1.1-1.2 and [S II]/Ha ~ 0.5-0.75 are therefore *comparable to* LAMOST's on-nebula
  values (0.88-1.11, 0.83-0.93), not lower.
* Ha suffers OH residuals (6553/6568/6577 A), so these are the least reliable ratios.

**New shock diagnostics.**
* [Ca II]7291 < 0.004, [Fe II]7155 < 0.0014 and [N I]5200 < 0.04 of [O III] (97.5 %).
* No [Ca II]7324 under [O II]7325 (-0.001).
* [N II]5755/6583 = 0.14 (0.06-0.26) hints at T_e([N II]) >~ 2e4 K (1 sigma).
* Together with [O III]/Hb > 11-18 and a 2-3 sigma hint of T_e([O III]) ~ (3-5)e4 K, the absence of [N I], [Ca II] and [Fe II] fits a
  shock with little or no cool recombination/neutral zone (an incomplete radiative shock, or an evaporative or
  mixing interface). It does not fit a complete radiative shock or photoionised gas.

## (d) Velocity (`verify_desi_velocity.py`, `verify_desi_velocity.json`, `figures/verify_desi_velocity.png`)

| Sample | IV regression centroid (km/s) | Positive fibres, S/N > 4 (km/s) |
|---|---|---|
| all | +11.0 (10.3-14.3) | +12.2 +- 0.8 (n = 346) |
| E strand / W strand | +13.0 / +8.5 | +15.0 +- 1.0 (123) / +13.7 +- 1.5 (82) |
| E / W, Dec 15.5-17.5 | +9.3 / +16.5 | +10.6 +- 2.2 (20) / +16.7 +- 2.1 (41) |
| Dec 12-15.5 / 15.5-17.5 / 17.5-19.5 | +8.5 / +13.3 / +18.5 | +11.5 +- 0.8 / +12.2 +- 1.5 / (3 fibres) |
| within 0.5 deg of P1 / P2 / P3 | (not local, see below) | **+11.9 +- 1.9 / +15.2 +- 2.9 / +9.2 +- 1.0** |
| dark / bright; main / SV | +11.8 / +9.3; +15.0 / +10.3 | +12.8 +- 1.0 / +10.0 +- 1.1; +16.9 +- 3.0 / +11.5 +- 0.7 |

* **Method caveat.** Regression centroids of *local* subsets are not local. The negative residuals of fibres in a
  petal carry the sky-fibre mean of the filament from elsewhere in that petal. A local subset may also span only
  1-2 bootstrap blocks, in which case the bootstrap collapses (P1 and P3 gave zero-width intervals). Local and strand
  velocities should therefore come from positive-residual fibres. Wave 1's strand and band velocities came from
  regression spectra and share this issue.
* **Comparison with LAMOST (Table 4).** DESI minus LAMOST: P1 -14.9 +- 4.7, P2 -11.6 +- 5.0, P3 -14.6 +- 3.4 km/s.
  The newly covered P3 now gives the strongest test. DESI V_LSR at P3 is -2.9 +- 1.0 km/s.
* **Frame test with the night sky** (single-night coadds, v_corr -23 to +30 km/s):
  * The variance-weighted centroid of the [O I]5577 variance peak follows v_corr with slope 0.869 +- 0.005
    (n = 5,342; rms 6.5 km/s).
  * A simulation of the same estimator on the 0.8 A grid, with the measured profile width (sigma ~ 1.0 A), gives a
    response of 0.883. The calibrated slope is therefore **0.98**.
  * Hg I 4358 gives 0.77 measured, 0.79 simulated, so **0.98**.
  * Gaussian fits to profiles stacked in six v_corr bins give 0.87 (5577) and 1.09 (Hg).
  * Topocentric wavelengths (slope 0) are excluded, so the DR1 coadds are barycentric.
* **[O III] velocity vs v_corr.** The slope is +0.49 +- 0.13, falling to +0.30 +- 0.15 with Dec fixed effects
  (v_corr correlates with Dec at r = 0.42). This reflects position confounding, not a frame error: a missing
  correction would give -1.
* **Tiles 41485/41487.** These are main/backup tiles at (111.56, +6.50) and (114.60, +5.08), observed 2022-03-16.
  None of the 7,345 usable coadds in the [O III] analysis include them; only 2 backup-program control coadds do,
  with no [O III]. Their 5577 residuals (+12 +- 3 km/s, n = 2) cannot test the -19 km/s stellar offset.

## (e) Coverage extension

* Fetched 2,400 coadds in about 45 min (SPARCL through the proxy at about 0.9 spectra/s, with frequent proxy
  disconnects and automatic retries).
* Near the trace at Dec < 15.5 there are now 2,207 usable fibres, against 649 in wave 1, and P3 is densely covered
  (`figures/verify_desi_extension_map.png`): the E and W strands show as positive ridges from Dec 14 to 15.5.
* DESI DR1 has essentially no coadds within 0.5 deg of the trace at Dec 12-13.7, so the southern end cannot be
  covered.
* 1,175 control coadds lie more than 5 deg from all strands. They show no structure: -0.007 +- 0.008 R,
  chi std 1.09.

## Second-pass audit (revision 2, resume; `verify_desi_audit*.py`, `figures/verify_desi_audit.png`)

Purpose: test the first-pass verdicts that overturned wave 1 with a method that shares no line-fitting code with
them. The audit reuses only the per-coadd metadata: the geometry template Z and the GFA acceptance f(lambda).

**Method.**
* Line fluxes are box integrals (+-2.4 A; +-2.0 A for 4363; doublet sums for [O II]) of flux - Redrock model, minus
  a weighted linear continuum from side bands 6-30 A away.
* The masks are:
  - an empirical sky-line mask (median 1/ivar > 1.25x the local median);
  - other nebular lines +-3 A;
  - the object's own redshifted lines +-8 A.
* Topocentric lines (Hg I 4358, [O I]5577) use boxes shifted by v_corr. LSF capture fractions are >= 0.999.
* Ratios use the IV estimator with Z as instrument and per-line inverse-variance weights. Errors come from a
  HEALPix nside-32 and nside-16 block bootstrap (400 resamples).
* Bootstrap calibration: 11 line-free null boxes near the key lines give ratio rms 0.020, which is
  **1.25x the nside-32 sigma and 1.68x the nside-16 sigma**. The worst null is 4381.6 A (0.048 +- 0.014). The sky
  line [O I]5577, which must give 0, gives **0.111 (0.079-0.147)** (first pass 0.09). So template-correlated
  systematics reach ~0.05 of [O III] in quiet regions and ~0.1 at strong sky lines.
* The first-pass Gaussian-fit ratios were calibrated with 17 fake lines (x1.5; x1.8 near 4363).

**Results (box method; 68 % bootstrap; `data/derived/verify_desi_audit.json`, `_extra.json`).**

| Quantity | all | Dec < 19.5 | Z > 0 | Z < 0 | first pass (Gaussian) |
|---|---|---|---|---|---|
| [O III] amplitude (R, GFA conversion) | 1.46 | **1.62** | | | 1.62 south median |
| same, wave-1 conversion | 2.12 | 2.35 | | | 2.4 (wave 1) |
| Hb | 0.012 (-0.015-0.031) | 0.018 | 0.027 | 0.047 | 0.021 (0.004-0.028) |
| Ha | 0.094 (0.077-0.169) | 0.082 | 0.115 | 0.001 | 0.082 |
| [N II]6583 | 0.099 (0.087-0.129) | 0.082 | 0.124 | 0.100 | 0.100 |
| [S II]6716 + 6731 | 0.066 | 0.043 | 0.095 | 0.067 | 0.062 |
| [O II]3727 sum | 1.41 (1.27-2.03) | 1.16 (1.09-1.24) | 1.34 | 1.44 | 1.37 |
| [O III]4363 | 0.091 (0.026-0.103) | 0.101 | 0.093 | 0.152 | 0.083 |
| Hg I 4358 (sky) | 0.004 (-0.025-0.013) | 0.012 | 0.047 | -0.020 | 0.009 |
| He II 4686 | -0.011 (< 0.013 at 97.5 %) | -0.010 | -0.015 | -0.009 | -0.008 |
| [Ne III]3869 | 0.069 (0.033-0.087) | 0.070 | | | 0.099 |
| [O III]4959 | 0.346 (0.336-0.381) | 0.344 | | | 0.332 |
| [S III]9069 | 0.037 | 0.031 | | | 0.028 |
| [O II]7325 sum | 0.081 (OH-contaminated box) | 0.079 | | | 0.036 (masked fit, preferred) |

Subsets disagree beyond their bootstraps. For example, Hb is 0.041 in the wave-1 fibres and -0.011 in the new
ones, and 4363 is 0.048 and 0.118. This supports inflating the errors.

**Pixel regression near 4363.** I regressed the continuum-subtracted residual at each pixel on Z, normalised by
the [O III] box slope, and bootstrapped it.
* With free centres, the 4363 component sits at **4364.73 A (16-84 %: 4363.68-4365.18)**. Expected: 4364.62
  ([O III] at +13 km/s); Hg I is at 4359.55 (5.1 A = 7 LSF sigma away); [Fe II]4359 is at 4360.6.
* With fixed centres:
  - 4363 = 0.057 (2.5/16/50/84/97.5 %: -0.047/0.015/0.055/0.069/0.076);
  - Hg = 0.024 (-0.052/0.011/0.025/0.037/0.075).
* A box at 4362.6-4366.6 A cannot receive Hg light at DESI resolution. Hg *residuals* also do not correlate with Z
  in the full sample.

**Verdict changes from this audit.**
* [O III]/Hb: from "47 (> 24)" to a lower limit of **> 11-18 (2 sigma)**, i.e. wave 1's > 11 is confirmed. Hb is at
  most a 1-2 sigma signal, so neither 20 nor 47 is a measurement.
* R23: from "~130" to **> ~30**.
* 4363: from "3-5 sigma" to **~2-3 sigma**. Hg remains excluded.
* Ha: the asymmetry between Z > 0 and Z < 0 fibres (0.115 vs 0.001, each +- ~0.09) is ~1 sigma. [N II]/Ha
  ~ 1.0-1.2 holds in both methods, but the red-line ratios are systematics-limited.
* Confirmed by the independent method:
  - [O II]/[O III] (1.41; south 1.16; north 4.2);
  - He II (< 0.013-0.026);
  - the GFA calibration revision (south 1.62 R by both methods);
  - no strong-[O III] fibre on tiles 41485/41487 (0 of 225 with box S/N > 4).

**Velocity (`verify_desi_audit_vel.py`, `data/derived/verify_desi_audit_vel.json`).**
* Method: per-pixel regression of continuum-subtracted residuals around **both** [O III]4960 and 5008 (4959 tied
  at 1/2.98) on Z_E, Z_W and Z_W2 jointly. The velocity comes from a chi^2 grid (0.25 km/s) of a pixel-integrated
  Gaussian, using 5,405 fibres with clean windows and an nside-32 block bootstrap.
* Results:
  - total +11.3 (10.2-13.9);
  - E +12.0 (11.3-14.3), W +11.4 (8.6-15.6);
  - Dec < 19.5: +10.6 (10.2-12.0); Dec >= 19.5: **+20.9 (18.5-23.1)**;
  - Dec 15.5-17.5: E +10.1 (5.9-10.1), W +13.0 (12.7-16.9), from only 9 blocks.
* Local on-ridge subsets (P1-P3) span 1-2 blocks, so their bootstraps collapse. For those I keep the first-pass
  positive-fibre values.

**Frame check (variance peaks).** Pixel-integrated Gaussian centroids of the 1/ivar peaks in coadds with a v_corr
spread < 3 km/s give slopes against v_corr of:
* 0.63 +- 0.015 ([O I]5577, n = 5,413);
* 1.16 +- 0.011 (Hg I 4358, n = 4,305).

The estimator is biased in different directions for the two lines (sampling at 0.8 A, non-Gaussian variance
profiles). The first pass's simulation-calibrated slope of 0.98 is the better number. Both raw slopes exclude a
topocentric frame (slope 0) and a doubled correction (slope 2), so the barycentric frame stands.

**Calibration bookkeeping (from the GFA table directly).**
* FIBER_FRACFLUX_GFA median 0.595 (dark 0.637, bright 0.446, backup 0.521; 196/238 exposures). The Legacy
  1"-Gaussian / 1.5"-fibre fraction is 0.790.
* Coadd f(5008) median 0.566. Ratio to MEAN_PSF_TO_FIBER_SPECFLUX: 0.736 (16-84 %: 0.549-0.885).
* desispec `psf_to_fiber_flux_correction` docstring (read as text): "flux one would collect ... in a fiber of 1.5
  arcsec diameter, for a 1 arcsec seeing ... (same definition as for the Legacy Surveys)". This confirms that it is
  not the exposure's acceptance.
* `fluxcalibration.py` divides the per-fibre calibration by FLAT_TO_PSF (`ccalibration /= point_source_correction`).
  So diffuse SB = F_cal x f_ps / Omega.

**Cross-instrument check (literature = teammate file `notes/verify_sdss.md`, read only).**
* The independent SDSS re-analysis finds a southern ON-OFF [O III] of 1.39 +- 0.10 R and per-plate strand peaks of
  1.6-2.5 R (south) and 0.4-0.8 R (north). Its fibre-acceptance correction f_eff = 0.84 (3") / 0.60 (2") was
  measured against imaging sky.
* My revised DESI values, 1.62 R south and 0.56-0.71 R north, agree. Before the two corrections the instruments
  differed by about 3x (wave-1 SDSS 0.79 R vs wave-1 DESI 2.4 R).
* Open tension: SDSS [O II]/[O III] south = 1.84 (95 %: 1.46-2.37) vs DESI 1.11-1.16 +- 0.08 (about 2.5-3 sigma).
  The ratio rises steeply to the north (DESI 4.2-4.4), so different Dec/strand weighting could explain part of it.

## New findings

1. **Calibration error in wave 1** (PSF_TO_FIBER_SPECFLUX used as the fibre acceptance). Every DESI surface
   brightness in wave 1 is about 1.36x too high (1.13-1.82 per fibre; worst for bright-program coadds). Line
   *ratios* are unaffected except through the chromatic term, which wave 1 also omitted: [O II]/[O III] -6 %,
   red lines/[O III] +5-12 %. After this correction DESI (1.6 R south) agrees with the corrected SDSS (1.4-2.5 R).
2. **[O III]/Hb > 11-18 and R23 > ~30 (2 sigma), with [O III]4363/5007 ~ 0.06-0.09 at about 2-3 sigma.** Hg is
   ruled out. This points to collisionally heated gas with weak Balmer emission. The first-pass central values
   ([O III]/Hb ~ 50, R23 ~ 130, 4363 at 3-5 sigma) are withdrawn after the null-box calibration.
3. **The filament's own [N II]/Ha (~1.0-1.2) and [S II]/Ha (~0.5-0.75) are not lower than the DIG/WIM values.**
   Wave 1's "weaker low-ionisation lines" conclusion came from its Hb-based rescaling, and Hb is not detected
   significantly. The Ha differential remains systematics-limited.
3a. **Template-correlated systematics are measurable:** the sky line [O I]5577 gives an IV "ratio" of 0.09-0.11
   (3 sigma), and a null box at 4381.6 A gives 0.048. Any DESI ratio below ~0.05 of [O III] is therefore an upper
   limit, whatever its formal significance.
3b. **[O III] velocity gradient:** Dec >= 19.5 is at +20.9 (18.5-23.1) km/s, against +10.6 (10.2-12.0) in the
   south.
4. **Upper limits** for [Ca II]7291, [Fe II]7155 and [N I]5200 (< 0.4 %, < 0.14 %, < 4 % of [O III]). These argue
   against a complete radiative shock with a cool recombination zone.
5. **The DESI-LAMOST [O III] velocity offset is about 12-15 km/s at all three regions**, including the newly
   covered P3 (4.3 sigma).
6. **Method caveat:** regression-spectrum velocities and ratios for small sky regions are not local, because the
   petal-mean subtraction carries the filament spectrum from elsewhere in the petal.
7. **The Kitt Peak Hg I 4358 line is present in DR1 sky variance** (+50 %). It is harmless for [O III]4363 at
   DESI's resolution.

## Limitations

* f_ps comes from DESI's GFA PSF model (centred source). Offsets (median 0.07"), the wavelength scaling and the
  GFA effective wavelength leave about +-8 % in absolute surface brightness. The forward model also inherits wave 1's
  assumptions: Gaussian cross-sections at the trace HWHM, and the science-fibre mean as a proxy for the sky-fibre
  mean (sky-fibre positions are not public). chi^2/dof is 1.46.
* Ratio bootstraps underestimate the fake-line scatter by about 1.5x, and the distributions are skewed (a few blocks
  dominate), so 95 % intervals are wide. The Hb, 4363 and red-line ratios depend on the instrument and subset at the
  +-30 % level. Ha is the least reliable line.
* The [O III] T_e formula is approximate at T > 3e4 K. The [O II] doublet deblending depends on the LSF (taken from
  the SPARCL `wave_sigma` of 430 pilot spectra).
* I did not re-check the [O III] line width, the Ha-structure component, the northern Ha patch ratios, or the
  stellar-RV zero points.
* The fetch was stopped at 2,400 of 3,005 coadds. The unfetched coadds are the farthest "near" ones and 25 controls.
  I did not resume it in revision 2: the remaining 580 "near" coadds lie 0.6-1 deg from the trace, mostly outside
  the ridge cross-section, and DR1 has no coadds within 0.5 deg of the trace at Dec 12-13.7.
* Revision 2: the box method is less efficient than Gaussian fits (matched-filter S/N 106 vs 154), and its [O II]7325
  box contains OH residuals. The audit's local (P1-P3) velocity bootstraps collapse (1-2 blocks), and its
  variance-peak frame estimator is biased (slopes 0.63 / 1.16), so the first-pass values are kept for those.
* The null calibration uses 11 boxes, so the inflation factor (1.25-1.68) is itself uncertain by about +-25 %.

## Files

* Scripts (run in this order, with `~/.venvs/astro313/bin/python -I`): `verify_desi_tap.py`,
  `verify_desi_fetch.py select|fetch`, `verify_desi_linefit.py 14`, `verify_desi_setup.py`,
  `verify_desi_detect.py`, `verify_desi_brightness.py`, `verify_desi_ratios.py`, `verify_desi_regspec.py`,
  `verify_desi_velocity.py`, `verify_desi_platescale.py`, `verify_desi_summary.py`, `verify_desi_figs.py`.
  Revision 2 audit: `verify_desi_audit.py` (box fluxes, IV ratios, nulls, 4363 pixel regression, calibration),
  `verify_desi_audit_vel.py` (joint E/W velocity, frame), `verify_desi_audit_fig.py`.
  Shared code is in `verify_desi_common.py`.
* Data products:
  * `data/derived/verify_desi_summary.json` (all key numbers);
  * `verify_desi_linefluxes.npz`, `verify_desi_coadd_meta.npz`;
  * `verify_desi_detect.json` / `_nulls.npz`, `verify_desi_brightness.json`, `verify_desi_ratios.json` /
    `_ratio_boots.npz`;
  * `verify_desi_regspec.npz` / `_fits.json`, `verify_desi_velocity.json` / `_sky.npz`,
    `verify_desi_platescale.json`, `verify_desi_sel_extend.csv`;
  * `verify_desi_windows/` (2,400 new ±40 A windows and fetch logs);
  * revision 2: `verify_desi_audit.json`, `verify_desi_audit_extra.json`, `verify_desi_audit_vel.json`,
    `verify_desi_audit_regspec.npz`, `verify_desi_audit_boxfluxes.npz`.
* Raw downloads:
  * `data/raw/verify_desi_tap/exposure_gfa_box.csv` (+ query);
  * `data/raw/verify_desi_src/` (desispec/desimodel source as text, URLs and date).
* Figures: `figures/verify_desi_detection.png`, `verify_desi_calibration.png`, `verify_desi_regspec.png`,
  `verify_desi_ratios.png`, `verify_desi_velocity.png`, `verify_desi_extension_map.png`, `verify_desi_audit.png`.
