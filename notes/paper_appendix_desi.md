# Appendix material: [O III] detection and line fitting in individual DESI DR1 spectra near P1 and P3

Date: 2026-11-18. Script: `scripts/paper_appendix_desi.py` (one run, about 30 s):
`SILK_THEME=light SILK_PAPER=1 ~/.venvs/astro313/bin/python -I scripts/paper_appendix_desi.py`.

Products:
- `data/derived/paper_appendix_desi_fibers.csv`: 860 rows, one per fiber per region (330 for P1, 530 for P3).
- `data/derived/paper_appendix_desi_summary.json`: every number below, with uncertainties and N.
- `data/derived/paper_appendix_desi_stacks.npz`: the stacked spectra.
- `paper/figures/paper_appendix_desi.pdf` / `.png`: 7.1 x 7.35 in, 12 panels (a)-(l), all text >= 6.5 pt.

Conventions: V_helio (DESI wavelengths are vacuum and barycentric; barycentric = heliocentric to < 0.02 km/s).
1 R = 10^6/(4 pi) photons cm^-2 s^-1 sr^-1. Ratios are energy ratios unless marked. Regions: P1 (116.65, +16.60),
P3 (116.45, +14.72), radius 0.5 deg (flat-sky separation, as in the paper).

## 0. What the spectra are (and are not)

**No DESI sky-fiber spectra were used.** Those files are hosted at NERSC and were not reachable. The inputs are
public DESI DR1 healpix **coadds of science targets** (galaxies and QSOs with ZWARN = 0), fetched through NOIRLab
SPARCL as +-40 A windows around the nebular lines (7,620 coadds in the whole study). Each spectrum is the pipeline
sky-subtracted, flux-calibrated flux minus the Redrock best-fit model of the target ("residual"). The DESI pipeline
subtracts one sky model per petal and exposure, built from that petal's sky fibers. A Galactic filament therefore
survives as a **positive** line residual in fibers on the ridge and a **negative** one in fibers beside it on the
same petal: every residual is relative to the petal's sky-fiber mean.

## 1. Bottom line

| Quantity | P1 | P3 |
|---|---|---|
| Coadds within 0.5 deg (fetched / good [O III] fit) | 330 / 326 | 530 / 505 |
| Positive-residual fibers with S/N > 4 | 61 | 128 |
| **Table 2 velocity, reproduced exactly** | **+11.92 +- 1.95 (61)** | **+9.21 +- 1.03 (128)** |
| ON-OFF stack velocity (petal mean cancels) | +12.3 +- 1.1 | +8.0 +- 0.9 |
| Other estimators (range) | +11.7 to +13.7 | +8.0 to +9.7 |
| W26 Table 4 minus Table 2 value | +14.9 +- 4.7 | +14.6 +- 3.5 |
| Model sky brightness at the region (global forward model) | 2.07 +- 0.06 R (paper: 2.1) | 1.54 +- 0.04 R (paper: 1.5) |
| E-ridge peak (global nodes / local B / local C) | 1.82 / 1.83 / 1.69 R | 1.68 / 1.65 / 1.34 R |
| ON-OFF mean brightness, on-ridge minus off-ridge | 1.36 +- 0.05 R | 1.65 +- 0.05 R |
| 5007/4959 (three estimators) | 2.87-2.97 (+-0.15-0.18) | 3.08-3.29 (+-0.17-0.19) |
| [O II]3727/[O III]5007 | 1.00-1.05 (+-0.13-0.14) | 0.89-0.99 (+-0.09-0.11) |
| [O III]5007/Hb (2 sigma lower limit) | > 8.5 (8.5-12.0) | > 18 (18.1-22.2) |
| He II 4686/[O III] (2 sigma) | < 0.08 | < 0.045 |
| Ha, [N II]6583, [S II] over [O III] (2 sigma) | < 0.16, < 0.13, < 0.12 | < 0.13, < 0.11, < 0.06 |

Three points the appendix should state plainly:
1. Table 2 is reproduced to the last digit. The quoted errors are the larger of the formal error and std/sqrt(n).
2. Per-fiber velocities scatter more than their formal errors (robust scatter 13-14 km/s against a median formal
   error of 9-11 km/s). Part of this is real structure amplified by the petal-mean subtraction (section 5.3).
3. **Calibration caveat found here.** The two SV1 tiles that supply 85 % (P1) and 100 % (P3) of the fibers have
   no guide-camera acceptance in the DR1 exposure table; program medians were used (section 6.4). The two tiles
   then disagree in surface brightness by a factor 1.6-1.8, while their calibrated fluxes agree. Near P1 and P3
   the absolute brightness is therefore uncertain by about +-20 %, more than the +-8 % quoted in the paper.

## 2. Constants

| Constant | Value |
|---|---|
| Rest wavelengths (vacuum, A; Edlen 1966 / IAU air-to-vacuum) | [O III] 5008.240, 4960.295; Hb 4862.683; He II 4687.021; [O II] 3727.092, 3729.875; Ha 6564.613; [N II] 6585.278, 6549.849; [S II] 6718.294, 6732.668 |
| Wavelength grid | 3600 + 0.8 i A |
| LSF sigma at 5008 A | 0.6951 A (median SPARCL `wave_sigma` of 430 spectra) |
| Model line sigma at 5008 A | sqrt(0.6951^2 + (lambda 12 km/s / c)^2 + 0.8^2/12) = 0.7594 A = 45.5 km/s (fixed, not fitted) |
| [O III] doublet | 5007 : 4959 = 2.98 : 1 (tied) |
| Fixed nebular velocity for fluxes | V_NEB = +14 km/s |
| Velocity scan | -120 to +150 km/s in 2 km/s steps, parabolic interpolation |
| Error scale from 17 line-free windows | 0.995 |
| Fiber | 1.52 arcsec diameter, Omega = 1.8146 arcsec^2 |
| 1 R at 5008.24 A | 7.4188e-18 erg s^-1 cm^-2 arcsec^-2 |
| Flux to rayleighs | SB [R] = F [1e-17 erg s^-1 cm^-2] x f x 0.7428; median f = 0.566 gives 0.4207 R per 1e-17 |

## 3. Method, step by step

1. **Sample.** From the SPARCL catalog of DR1 coadds, keep GALAXY and QSO with ZWARN = 0 (stars are excluded
   because their z = 0 Balmer and Mg b model mismatch falls on the lines). Take every fetched coadd with
   hypot((RA - RA0) cos Dec0, Dec - Dec0) < 0.5 deg.
2. **Residual.** r(lambda) = coadd flux - Redrock model, with the pipeline inverse variance.
3. **Per-fiber fit** (`verify_desi_linefit.fit_group`, the paper's fitter, re-run here on the region fibers; it
   returns the stored values exactly, max |difference| = 0).
   - Window: [4960.3 - 30, 5008.2 + 30] A.
   - Model: quadratic continuum + A [G(5007) + G(4959)/2.98], fitted simultaneously by weighted linear least
     squares. G is a unit-integral Gaussian of the fixed width above, so A is the integrated 5007 flux.
   - Masks: other nebular lines +-3 A; sky lines (Hg I 4047/4358/5461/5770/5791, Na D, [O I] 5577) +-2.5 A; the
     target's own redshifted lines +-8 A; one 5 sigma clip outside +-3 A of the fitted lines.
   - Flux F and error sigma_F at V_NEB = +14 km/s; S/N = F / (0.995 sigma_F).
   - Velocity: chi^2 of the same linear fit on the velocity grid; V is the parabola minimum and its error is
     sqrt(2 / curvature). Centroid lambda = 5008.240 (1 + V/c).
4. **Quality mask** (`good_fit`): finite fit, error < 4 x the sample median error, >= 60 % of the window pixels
   unmasked, |S/N| < 30.
5. **Surface brightness.** DESI calibrates each fiber for point sources, so a uniform source gives
   F_cal = SB Omega / f, hence **SB = F f / Omega**. f is the point-source fiber acceptance of the coadd: the
   EFFTIME-weighted harmonic mean of the guide-camera (GFA) values of its exposures, moved to the line wavelength
   with a Moffat (beta = 3.5) profile whose FWHM scales as lambda^-0.2.
6. **Ridge template and petal term.** T(x) = sum over the E, W and W2 strands of unit-peak Gaussian
   cross-sections with the traced half width. The predicted petal-relative residual of fiber i is
   Z_i = T(x_i) - sum_t w_t <T>_(tile t, petal), where the mean is over all fiberassign positions of the
   tile-petal (a proxy for its sky fibers) and w_t are EFFTIME weights.
7. **Table 2 velocity.** Select good fibers with S/N > 4 (positive residual), finite V and sigma_V < 30 km/s.
   Weights w = 1/(sigma_V^2 + 5^2); V = sum(w V)/sum(w); error = max(sum(w)^-1/2, std(V)/sqrt(n)).
8. **Stacks.** Each residual is converted to energy surface brightness in rayleigh-equivalents at 5008 A,
   e_i = r_i f_i(lambda) 1e-17 / (Omega R_5008), so that ratios of stack fluxes are energy ratios and [O III] 5007
   is in true R. Weights: **w_i(lambda) = ivar_i / kappa_i^2** (inverse variance in those units), zero within
   +-8 A of the target's own lines. Stack = sum(w e)/sum(w) per pixel.
   - ON: good fibers in the region with T >= 0.5.
   - OFF: good fibers on the **same tile-petals** with T < 0.05, anywhere on that petal (median distance from
     the region centre 0.48 deg for P1 and 0.70 deg for P3).
   - ON-OFF is formed per tile-petal (at least 5 ON and 5 OFF fibers) and combined over petals with per-pixel
     inverse-variance weights. The petal sky mean cancels exactly in ON-OFF.
   - Lines in the stacks are fitted with the same fitter.
9. **Ridge amplitude.** Four estimates (section 6).
10. **Line ratios.** Two estimators (section 7), each with bootstrap over fibers and null tests.

## 4. Inventory of coadds

| | P1 | P3 |
|---|---|---|
| DR1 coadds in the SPARCL catalog within 0.5 deg | 1,467 (1,063 GALAXY, 181 QSO, 223 STAR) | 841 (613, 44, 184) |
| GALAXY/QSO with ZWARN = 0 | 1,034 | 531 |
| Fetched and analyzed | 330 (32 %; the fetch queue was cut short) | 530 (99.8 %) |
| Good [O III] fit | 326 | 505 |
| Program | 197 dark, 133 bright | 259 dark, 271 bright |
| Survey | 280 SV1, 50 main | 530 SV1 |
| Target class | 295 GALAXY, 35 QSO | 488 GALAXY, 42 QSO |
| Redshift (min / median / max) | 0.0 / 0.54 / 3.27 | 0.0 / 0.42 / 4.11 |
| Coadd exposure time (s) | 648 / 4,501 / 5,406 | 600 / 4,501 / 5,361 |
| EFFTIME (s) | 169 / 1,606 / 3,494 | 1,606 / 1,606 / 3,494 |
| Offset from the E spine (deg) | -0.49 to +0.49 | -0.61 to +0.38 |
| Median distance to the nearest strand (deg) | 0.115 | 0.185 |
| Fibers with T >= 0.5 / T < 0.05 | 169 / 103 | 292 / 77 |
| Template at the region centre (E, W) | 0.994, 0.178 | 0.808, 0.155 |
| Region centre relative to the E spine; local HWHM | 0.007 deg W; 0.071 deg | 0.115 deg W; 0.207 deg |

Tiles: 80683 (SV1 dark; 8 exposures on 4 nights, 2021-01-10 to 02-06; EFFTIME 3,494 s) and 80644 (SV1 bright;
15 exposures on 6 nights, 2021-01-02 to 03-28; EFFTIME 1,606 s) share the tile centre (116.0, +15.5). They give
280 of 330 fibers near P1 (petals 6 and 5) and all 530 near P3 (petals 9, 0 and 8). Near P1 there are also 50
main-survey fibers on single-night tiles 2385, 7836, 23703 and 4703. The barycentric corrections of the SV1
exposures span 34.5 km/s (80644) and 13.9 km/s (80683), so topocentric sky residuals are smeared in these coadds
while nebular lines are not.

## 5. Velocities

### 5.1 Per-fiber fits

| | P1 | P3 |
|---|---|---|
| Median flux error (1e-17 erg s^-1 cm^-2 / R) | 0.573 / 0.214 | 0.628 / 0.198 |
| Fibers with S/N > 4 / S/N < -4 / in between | 61 / 39 / 226 | 128 / 11 / 366 |
| Largest S/N | 19.4 | 15.9 |
| Range of fiber brightness (R, petal-relative) | -1.68 to +3.31 | -1.13 to +3.40 |
| Median reduced chi^2 of the fits | 0.96 | 0.96 |

### 5.2 Table 2 values and their distribution

| | P1 (n = 61) | P3 (n = 128) |
|---|---|---|
| **Weighted mean (paper's estimator)** | **+11.92** | **+9.21** |
| Formal error / std/sqrt(n) / bootstrap over fibers | 1.35 / 1.95 / 1.87 | 0.87 / 1.03 / 1.14 |
| **Quoted error = max(formal, scatter)** | **1.95** | **1.03** |
| Median +- 1.25 sigma_MAD/sqrt(n), sigma_MAD = 1.4826 MAD | +13.7 +- 2.2 | +8.5 +- 1.4 |
| Robust scatter (1.4826 MAD) / std | 14.0 / 15.2 | 12.9 / 11.7 |
| 16-84 % of the fiber velocities | -4.9 to +24.3 | -3.1 to +19.8 |
| Median formal error | 10.6 | 9.3 |
| Reduced chi^2 (formal errors / with the 5 km/s floor) | 2.52 / 1.85 | 2.18 / 1.38 |
| Extra scatter that gives reduced chi^2 = 1 | 11.0 | 7.5 |
| Mean with that extra scatter | +11.9 +- 1.9 | +9.1 +- 1.0 |
| Median S/N of the selected fibers | 5.6 | 6.1 |

The formal errors are too small by a factor of about 1.3-1.6. The quoted errors already allow for this, because
they use the scatter. They do not include the +-1.5 km/s zero-point uncertainty.

### 5.3 Alternative estimators

| Estimator | P1 | P3 |
|---|---|---|
| S/N > 3 / > 5 / > 6 / > 8 (same weights) | +11.7 +- 1.8 / +12.8 +- 2.6 / +13.6 +- 2.4 / +12.3 +- 2.2 | +9.1 +- 1.0 / +9.7 +- 1.1 / +10.3 +- 1.3 / +9.2 +- 1.8 |
| Stack of the Table 2 fibers, free velocity (bootstrap) | +13.0 +- 2.0 | +9.2 +- 1.1 |
| **ON-OFF stack, all on-ridge fibers** (161+175; 262+157 fibers) | **+12.3 +- 1.1** | **+8.0 +- 0.9** |
| ON-OFF with the Table 2 fibers as ON | +12.6 +- 1.3 | +9.7 +- 0.8 |
| ON-only stack of all T >= 0.5 fibers (biased, see below) | +10.4 +- 1.7 | +5.9 +- 1.2 |
| OFF stack (negative line = petal sky mean) | +13.4 +- 1.7 | +10.6 +- 1.7 |
| Negative-residual fibers, S/N < -4 | +11.9 +- 2.0 (39) | +13.2 +- 4.4 (11) |

**Why ON-only velocities are not local.** A residual is a_i G(v_i) - abar G(vbar), where abar and vbar are the
amplitude and velocity of the petal sky mean. Its centroid is (a_i v_i - abar vbar)/(a_i - abar), so any
difference between the fiber and the petal mean is amplified by T/(T - <T>_petal). For the Table 2 fibers the
median T is 0.91 (P1) and 1.01 (P3), the petal mean is 0.54 and 0.60, and the median amplification is 1.9 and
2.2. ON-OFF removes the petal term and is the cleanest local estimator. It agrees with Table 2 within 1 sigma
in both regions.

**Structure inside the regions (ON-OFF, so not amplified).**
- P3: east of the E spine +11.2 +- 1.1, west +5.5 +- 1.0 (difference 5.7 +- 1.4, about 4 sigma); south half
  +9.9 +- 1.1, north half +6.7 +- 1.0; fibers on the E strand alone +10.3 +- 1.2.
- P1: E strand alone +12.4 +- 1.9, W strand alone +14.3 +- 1.5; east +10.9 +- 1.7, west +12.6 +- 1.1.
- Per-fiber (amplified) values show the same pattern more strongly, for example P3 west +4.7 +- 1.2 and east
  +14.8 +- 1.6, and a slope of +52 +- 10 km/s per degree of offset. This explains much of the excess scatter.

The unbiased estimators lie 13-16 km/s below the W26 Table 4 values (+26.8 and +23.8), and at or below the
adopted +13.3 +- 1.0 km/s.

### 5.4 Zero point and frame

- **Tile to tile** (independent nights, ON-OFF): 80644 minus 80683 = +2.9 +- 2.9 km/s (P1) and -1.5 +- 2.1 km/s
  (P3). Main-survey tiles minus SV1 tiles near P1: -1.6 +- 3.0 km/s. No tile-level offset is seen.
- **Frame** (quoted from `verify_desi_velocity.json`, not recomputed): the [O I] 5577 variance peak follows the
  barycentric correction with slope 0.869 +- 0.005 in 5,342 single-night coadds; the estimator's simulated
  response is 0.883, so the calibrated slope is 0.98. Topocentric wavelengths (slope 0) are excluded.
- **Stars** (quoted from `spectra_desi_rvcheck.json`): DESI minus Gaia DR3, through LAMOST LRS, is -0.3 +- 1.0
  km/s (SV1 bright, 82 stars), +2.2 +- 1.8 (SV1 dark, 45), +1.2 +- 0.2 (main bright) and +2.7 +- 0.6 (main dark).
  The paper carries +-1.5 km/s.
- A local [O I] 5577 check is **not** possible: per-tile centroids of the variance peak scatter by about 10 km/s
  around the expected value for the single-night tiles (for example +17.9 +- 4.8 against +6.4 for tile 2385), and
  the SV1 coadds mix nights. The numbers are in the JSON for completeness.

### 5.5 Line width

The [O III] 5007 line in the ON-OFF stacks is unresolved: observed sigma = 0.753 +- 0.016 A (P1) and
0.715 +- 0.017 A (P3), against an instrumental 0.733 A (LSF plus pixel). Formal intrinsic FWHM: 25 (11-33) km/s
at P1 and 0 (< 22 at 97.5 %) at P3. The LSF is known to 2-3 %, which alone corresponds to about +-10 km/s, so
this is an upper limit of roughly 40 km/s, consistent with the study's earlier global value of <= 33 km/s
(`notes/spectra_desi.md`; the comment paper itself quotes no DESI line width).

## 6. Line strength

### 6.1 Per-fiber brightness against ridge offset

Weighted means of the petal-relative residuals (lower limits to the sky brightness):

| | P1 | P3 |
|---|---|---|
| Within 0.1 deg of the region centre | 1.18 +- 0.05 R (17) | 0.74 +- 0.04 R (23) |
| Within 0.2 deg | 0.44 +- 0.03 R (62) | 0.73 +- 0.02 R (84) |
| On-ridge fibers (T >= 0.5) | 0.70 +- 0.02 R (166) | 0.79 +- 0.01 R (271) |
| Off-ridge fibers (T < 0.05) in the region | -0.75 +- 0.02 R (102) | -0.16 +- 0.02 R (77) |
| Brightest single fiber | 3.3 +- 0.2 R | 3.4 +- 0.6 R |

The 0.1 deg values equal those in `notes/verify_desi.md`. Panel (i) of the figure shows the binned profile
against offset from the E spine with the forward-model prediction: a single ridge at P3, and the E and W ridges
(0.27 deg apart) at P1.

### 6.2 Global forward model (the paper's), reproduced

Unit-peak ridges times triangular 1 deg Dec nodes, minus the tile-petal mean, weighted least squares with 4 sigma
clipping over 7,345 fibers: chi^2/dof = 10,522/7,204 and a southern node median of **1.616 R** (paper: 1.62 R).
Nodes near the regions: E@14 2.31 +- 0.05, E@15 1.44 +- 0.05, E@16 1.70 +- 0.06, E@17 1.90 +- 0.08; W@16
1.24 +- 0.06, W@17 1.62 +- 0.06 R.

| From the global model | P1 | P3 |
|---|---|---|
| Sky brightness at the region centre | 2.07 +- 0.06 R | 1.54 +- 0.04 R |
| E-ridge peak at the region's Dec | 1.82 +- 0.06 R | 1.68 +- 0.05 R |
| W-ridge peak at the region's Dec | 1.47 +- 0.05 R | none (W starts at Dec 15) |
| Median petal term of the region fibers | 0.85 R | 0.97 R |

Errors are bootstraps over fibers of the global fit. W26 Table 3 corresponds to 12.75 R (P1) and 10.12 R (P3),
so the ratio to the model sky brightness is 6.2 and 6.6.

### 6.3 Local fits restricted to each region (bootstrap over fibers, 1,000 resamples)

| Fit | P1 (E, W amplitudes) | P3 (one amplitude) |
|---|---|---|
| B: y_i + S_i = a T_i, petal term S_i from the global model | 1.83 +- 0.06, 1.41 +- 0.05 R | 1.65 +- 0.03 R |
| C: y_i = a T_i + c_petal, free constant per tile-petal | 1.69 +- 0.07, 1.22 +- 0.06 R | 1.34 +- 0.05 R |
| C with one amplitude for E + W | 1.42 +- 0.05 R | 1.34 +- 0.05 R |
| D: y_i = a Z_i (uniform amplitude over the petals) | 1.44 +- 0.05 R | 1.38 +- 0.04 R |
| ON-OFF stack, mean on-ridge minus off-ridge | 1.36 +- 0.05 R | 1.65 +- 0.05 R |

chi^2/dof is 2.2-2.9 for all local fits, so the Gaussian-ridge template is only approximate; the bootstrap errors
allow for this but the spread between B and C (0.15-0.3 R) is the better guide to the model uncertainty.

### 6.4 Calibration test: the two SV1 tiles

The DR1 exposure table has `fiber_fracflux_gfa` = 0 (not measured) for all 23 exposures of tiles 80644 and 80683.
The study filled them with program medians, which give f(5008) = 0.419 (bright) and 0.608 (dark). Among tiles
with measurements, tile-level acceptances span 0.26-0.70 (bright, 16-84 %, 33 tiles) and 0.49-0.69 (dark, 54
tiles). Filled values apply to 24 % of all usable fibers and 46 % of those at Dec < 19.5.

| Fixed-effects amplitude per tile | P1 | P3 |
|---|---|---|
| Tile 80644 (R, with f = 0.419) | 0.94 +- 0.06 | 0.99 +- 0.07 |
| Tile 80683 (R, with f = 0.608) | 1.70 +- 0.07 | 1.53 +- 0.07 |
| Flux ratio 80644/80683 (acceptance-free) | 0.81 +- 0.06 | 0.94 +- 0.08 |
| Expected flux ratio for the adopted acceptances | 1.45 | 1.45 |
| Implied f(80644)/f(80683) | 1.24 +- 0.09 | 1.07 +- 0.09 |
| Tiles with measured acceptance (2385, 7836; 43 fibers) | 1.62 +- 0.16 R | none |
| Both SV1 tiles with one acceptance: 0.608 / 0.514 / 0.419 | 1.62 / 1.37 / 1.12 R | 1.52 / 1.29 / 1.05 R |

The two tiles see the same sky and return the same calibrated flux within 6-20 %, so their true acceptances were
similar, not 0.42 and 0.61. The measured-acceptance tiles near P1 imply a common value of 0.61 +- 0.06, which
would raise the one-amplitude value at P1 from 1.42 to about 1.6 R. With the adopted fills the two errors partly
cancel in the regional mean. **Adopt +-20 % for the absolute brightness near P1 and P3.** None of this changes the
factor of about 6 to W26.

## 7. Line ratios

### 7.1 Estimators

1. **Regression (the paper's family).** rho_L = sum u (Zf - <Zf>)(c L - <c L>) / sum u (Zf - <Zf>)(O - <O>),
   with per-fiber line fluxes L at V_NEB, O the tied [O III] 5007 flux, Zf = Z/k_R the geometry-only template in
   flux units, c = f(lambda_L)/f(5008) the chromatic acceptance, u = 1/sigma(cL)^2, and weighted means <.>
   (intercept). A second version de-means within each tile-petal (petal fixed effects).
2. **Direct ON-OFF stack.** Ratio of the Gaussian-fit fluxes in the ON-OFF stack of section 3.

Errors: bootstrap over fibers (1,000 resamples for the regression, 500 for the stack; the paper's HEALPix-block
bootstrap collapses on a 1 deg region).

Null tests:
- the 17 line-free windows of the paper, processed identically. Their ratios have rms 0.032 / 0.040 / 0.029 (P1;
  intercept, fixed effects, stack) and 0.018 / 0.022 / 0.022 (P3), with largest values 0.07-0.13 and 0.04-0.06;
  their scatter is 0.98-1.12 times the bootstrap sigma, so the bootstrap is realistic here;
- 300 random OFF-minus-OFF splits per region (no ridge signal). Line fluxes scatter by 1.0-1.5 times the formal
  errors, and by 1.9-2.3 times for Ha.

Adopted uncertainty of each ratio: max(bootstrap sigma x null inflation, null-window rms, split-null sigma).
Limits are max(ratio, 0) + 2 sigma.

Check added by the independent re-derivation (2026-11-18): the 17 null windows are not alike. The seven blue
windows (3760-4395 A) scatter by 0.04-0.06 (P1) and 0.02-0.03 (P3), the five green ones (4700-5215 A) by
0.017-0.019 and 0.008-0.019, the five red ones (6615-7350 A) by 0.009-0.012 and 0.012-0.016. The single floor
is therefore set by the blue windows. For [N II], [S II], Hb and He II the adopted sigma is **larger than the
band-matched one** by up to a factor of 2 at P1 and 1.5 at P3 (conservative). The one exception is Ha at P3: the
red windows scatter 1.3-1.4 times their bootstrap sigma, so the matched sigma (0.026-0.030) exceeds the adopted
0.022-0.023 and the 2 sigma limit would be < 0.14 instead of < 0.13. The adopted values and limits are kept. The
JSON also carries the same recipe with the null windows of the line's own band (`sig_matched`, `snr_matched`,
`_nulls.by_band`): matched sigma for Hb 0.019-0.025 (P1) and 0.017-0.021 (P3); for [N II] 6583 0.015-0.020 and
0.015-0.016; for Ha 0.027-0.037 and 0.026-0.030.

### 7.2 Results (ratio to [O III] 5007; regression with intercept / with petal fixed effects / ON-OFF stack)

| Ratio | P1 | P3 |
|---|---|---|
| [O III] 4959 | 0.341 / 0.336 / 0.347 (+-0.016-0.018 bootstrap) | 0.325 / 0.329 / 0.310 (+-0.016-0.017) |
| 5007/4959 (single-line fits) | 2.91 / 2.97 / 2.87 (+-0.15-0.18) | 3.12 / 3.08 / 3.29 (+-0.17-0.19) |
| [O II] 3726+3729 | 1.00 / 1.03 / 1.05 (+-0.13-0.14) | 0.92 / 0.89 / 0.99 (+-0.09-0.11) |
| [O II] 3729/3726 | 1.16 / 1.07 / 1.05 (+-0.19-0.24) | 1.28 / 1.26 / 1.32 (+-0.22-0.24) |
| Hb | 0.039 / 0.038 / 0.026 (+-0.029-0.040) | 0.010 / -0.010 / 0.012 (+-0.019-0.023) |
| [O III]/Hb, 2 sigma | > 9.8 / > 8.5 / > 12.0 | > 21.3 / > 22.2 / > 18.1 |
| He II 4686 | -0.027 / -0.049 / -0.042; < 0.06-0.08 | -0.008 / 0.000 / -0.020; < 0.04-0.045 |
| Ha | 0.071 / 0.084 / 0.036 (+-0.032-0.040); < 0.16 | 0.080 / 0.013 / 0.039 (+-0.023); < 0.13 |
| [N II] 6583 | 0.037 / 0.049 / 0.038 (+-0.029-0.040); < 0.13 | 0.063 / 0.045 / 0.068 (+-0.019-0.023); < 0.11 |
| [S II] 6716+6731 | 0.037 / 0.039 / 0.061 (+-0.029-0.040); < 0.12 | 0.027 / 0.005 / -0.014 (+-0.019-0.023); < 0.06 |
| [O III]/Ha, photons, 2 sigma | > 4.6 | > 6.1 |

Reading the table:
- **Detected:** the [O III] doublet (5007/4959 consistent with 2.98 in both regions) and [O II] (7-11 sigma).
  [O II]/[O III] is 1.0-1.05 at P1 and 0.9-1.0 at P3, within about 1 sigma of the paper's 1.11-1.16 for
  Dec < 19.5. [O II] 3729/3726 is consistent with the paper's 1.31 +- 0.07 but adds no weight.
- **Not detected:** Hb, He II and [S II] (all below 2.1 sigma). The Hb limits bracket the paper's global > 11:
  weaker at P1, stronger at P3. The He II limits are weaker than the global < 0.026.
- **Marginal:** Ha at P3 reaches 3.5 sigma in one estimator only (0.6 and 1.7 sigma in the others; 2.7, 0.5 and
  1.4 sigma against the red null windows). It is not a detection.
- **[N II] 6583 at P3: probable, not secure.** The ratio is 0.062 / 0.045 / 0.068, that is 3.4 / 2.0 / 3.1 sigma
  with the adopted (blue-dominated) uncertainty and **4.1 / 2.8 / 4.7 sigma against the red null windows**. In the
  stacks the line is +0.052 R-equivalent in ON, -0.028 in OFF and 0.112 +- 0.018 in ON-OFF (the sign pattern of a
  ridge line, unlike Ha, which is positive in both). An independent re-derivation with separate code finds
  0.066 +- 0.011 (bootstrap), present in both tiles (0.078 +- 0.019 and 0.060 +- 0.014) and on both sides of the
  spine, with [N II] 6548 at 0.017 +- 0.009 (1/3 expected), against an rms of 0.010 at 73 line-free red
  positions. Its velocity, +21 +- 6 km/s, is 2 sigma from [O III]. In energy units this is about 0.11 R-equivalent
  at 5008 A (0.15 R at 6585 A). The data cannot tell filament [N II] from WIM [N II] that is brighter on the ridge
  fibers than on the off-ridge fibers of the same petal. Quote it as a 3-5 sigma feature and keep < 0.11 as the
  2 sigma limit; at P1 the same line is 0.04-0.05 (2.2-2.5 sigma, matched nulls).
- **Red lines.** Ha is positive at 0.12-0.19 R (5008 A energy equivalent) in **both** the ON and the OFF stacks,
  that is, relative to the sky fibers everywhere on these petals. It cancels in ON-OFF and in the regression
  intercept, but smooth WIM that varies across a petal cancels only partly, and the split null shows Ha
  fluctuating at twice its formal error. The Ha, [N II] and [S II] entries are upper limits on emission that
  follows the [O III] ridge.

## 8. Consistency with the paper

| Paper value | Here | Status |
|---|---|---|
| Table 2: +11.9 +- 1.9 (61) near P1 | +11.92 +- 1.95 (61) | reproduced |
| Table 2: +9.2 +- 1.0 (128) near P3 | +9.21 +- 1.03 (128) | reproduced |
| Southern-ridge median 1.62 R | 1.616 R (same global fit) | reproduced |
| Model sky brightness 2.1 R (P1), 1.5 R (P3) (`verify_desi.md`) | 2.07, 1.54 R | reproduced |
| Residuals within 0.1 deg: 1.18 +- 0.05 (17), 0.74 +- 0.04 (23) R | identical | reproduced |
| Acceptance median 0.566 at 5008 A | 0.566 | reproduced |
| Acceptance systematic +-8 % | +-20 % near P1 and P3 (unmeasured SV1 tiles) | **larger here** |
| 5007/4959 = 3.04 | 2.87-2.97 (P1), 3.08-3.29 (P3), +-0.15-0.19 | consistent |
| [O II]/[O III] = 1.11-1.16 (Dec < 19.5) | 1.00-1.05 (P1), 0.89-0.99 (P3) | consistent (about 1 sigma low) |
| [O III]/Hb > 11 (2 sigma) | > 8.5 (P1), > 18 (P3) | consistent |
| He II/[O III] < 0.026 | < 0.08 (P1), < 0.045 (P3) | weaker locally |
| Adopted V_helio +13.3 +- 1.0 | P1 +11.9 to +12.3; P3 +8.0 to +9.2 | P1 about 1 km/s lower; P3 4-5 km/s lower |

## 9. Caveats

- **Not sky fibers.** Science-target coadds only; the petal sky model itself was never inspected.
- **Petal-relative residuals.** Every flux is relative to the petal's sky-fiber mean, which is modeled with the
  mean over all fiberassign positions because sky-fiber positions are not public. Near P3 the ridge covers most
  of petal 9 (petal mean 0.60-0.62 of the peak), so the petal term is as large as the signal.
- **Velocities of residuals are amplified** by T/(T - <T>_petal), about 2 here. Table 2 values are means over
  fibers with real velocity structure of 5-6 km/s inside a region; the ON-OFF values are the cleaner local numbers.
- **S/N cut at a fixed velocity.** The S/N > 4 selection uses the flux at +14 km/s. At the 45 km/s line sigma
  the amplitude loss for a 20 km/s offset is 5 % (exp(-dv^2/4 sigma^2)), so the bias is small; cuts from 3 to 8
  move the means by <= 1.7 km/s (P1) and <= 1.1 km/s (P3). A Monte Carlo of the selection (JSON
  `selection_bias_simulation`; white noise, many faint fibers) gives a pull of the selected mean towards +14 km/s
  of +0.2 km/s for a true velocity of +9, -0.3 km/s for +24 and -1.1 km/s for +34 km/s: the cut cannot hide a
  line at the W26 velocities. In the same simulation the velocity pulls have rms 1.03-1.05, so the formal errors
  are right for pure noise and the excess scatter of the data (reduced chi^2 2.2-2.5) is not a fitting artifact.
- **Redshift-model subtraction.** Redrock mismatch adds noise at random observed wavelengths; the target's own
  lines are masked (+-8 A). Limiting to QSOs and z > 0.5 galaxies changes the Table 2 means by < 0.5 km/s.
- **Fixed line width.** Fluxes assume sigma = 0.759 A. The stacked line is unresolved, so the bias is small, but
  the LSF comes from 430 spectra and is good to 2-3 %.
- **Acceptance.** Unmeasured for the two SV1 tiles (section 6.4): +-20 % on the absolute brightness near P1 and
  P3. The same fills enter 46 % of the fibers behind the paper's southern-ridge median.
- **WIM cancellation.** Only emission that follows the [O III] ridge is measured. Smooth WIM cancels, WIM
  structure within a petal does not.
- **Small numbers and few petals.** Each region is covered by two tiles and two to six tile-petals. Bootstraps
  over fibers ignore errors shared by a petal; the null windows (scatter 0.98-1.12 of the bootstrap) suggest
  this is not a large effect, but they cannot test it at line wavelengths. Only 32 % of the eligible coadds near
  P1 were fetched.
- **Not recomputed here:** the global frame test and the stellar zero points (section 5.4).

## 10. Figure (`paper/figures/paper_appendix_desi.pdf`)

Rows 1 and 2 are P1 and P3; row 3 compares the two regions.
- (a), (e) Fibers with a good [O III] fit, colored by petal-relative brightness (R). Dotted lines: E and W
  strands. Dashed circle: 0.5 deg. Square: region centre. Rings: the three fibers of the next panel.
- (b), (f) Continuum-subtracted residual spectra (R per A) of the brightest, the median and the faintest
  selected fiber, with the tied-doublet fit at the fitted velocity. Labels: S/N, offset from the E spine, V_helio.
- (c), (g) Per-fiber V_helio against S/N. Blue: Table 2 selection. Gray: 2.5 < S/N <= 4. Black line: weighted
  mean. Gray band: adopted +13.3 +- 1.0 km/s. Dashed red: W26 Table 4. Dotted: S/N = 4.
- (d), (h) Stacks around Hb, 4959 and 5007: ON-OFF with its fit, ON and OFF.
- (i) Binned residual against offset from the E spine, with the forward-model prediction.
- (j), (k) ON-OFF stacks at [O II] 3726, 3729 and at Ha, [N II], with fits.
- (l) Ratios to [O III] 5007 from the regression with intercept (filled) and the ON-OFF stack (open), with the
  adopted uncertainties; dotted line: 1/2.98.
