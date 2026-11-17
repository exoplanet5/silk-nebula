# SDSS and LAMOST spectra of the Silk Nebula field (topic `spectra_sdss_lamost`)

Access date for all online sources: **2026-11-16**. All wavelengths are vacuum. Velocities are heliocentric (V_helio) unless
marked V_LSR. Intensities are in Rayleigh (R, photon units). Line ratios are photon ratios: multiply by λ_den/λ_num for energy
ratios (×0.97 for [O III]/Hβ, ×1.31 for [O III]5007/Hα).
Labels: **measured** = fitted here from public spectra; **derived** = computed from measured values plus stated assumptions;
**literature** = taken from a paper or a teammate's file (the file is named).

## 1. Summary

1. **Absolute [O III] from SDSS sky fibres (measured).** I used 9,092 SDSS SKY-fibre spectra (7,958 DR17 SDSS-I/II/BOSS/eBOSS + 1,134
   SDSS-V DR19), grouped into 285 plate–spectrograph groups; 256 groups passed the sky-brightness and fit-quality cuts. The [O III]
   strand shows up directly in the per-fibre map. Fitted strand amplitudes above the local background are:
   - **W strand (Dec 15–20°): 1.40 ± 0.04 R**
   - **E strand, Dec 14–15.5°: 1.49 ± 0.18 R**; Dec 15.5–17°: 1.17 ± 0.12 R
   - north extension, Dec 20–23°: 0.78 ± 0.04 R; Dec 23–29°: 0.43 ± 0.03 R

   The off-strand background is 0.3–0.4 R, with a ±0.2 R zero-point systematic. Within 0.3° of P1, P2 and P3 the mean [O III] is
   1.42 ± 0.14, 1.80 ± 0.20 and 1.70 ± 0.37 R. The brightest single fibres there reach 1.9, 3.5 and 2.3 R.
   The paper's imaging surface brightness converts to 10–13 R (paper-review teammate). The SDSS values are **≈ 5–7× fainter**.
2. **[O III]/Hβ (measured).** In the south (Dec < 20°: 141 ON fibres, 26 groups), the ON-minus-OFF strand spectrum gives:
   - [O III]5007 = 0.79 ± 0.02 R
   - Hβ = 0.017 ± 0.020 R, so **[O III]/Hβ > 13.8 (2σ)**
   - [O III]/Hα = 4.7 ± 1.0

   On the W strand alone [O III]/Hβ = 22 ± 11 (> 10.8 at 2σ). The north (Dec ≥ 20°) is different: [O III]/Hβ = 7.4 ± 2.5,
   [N II]6583/Hα = 1.85 ± 0.36 and [S II]/Hα = 1.77 ± 0.37. The LAMOST DR12 LRS stellar-residual stack independently gives
   ON−OFF [O III] = 0.88 ± 0.13 R with Hβ = −0.08 ± 0.09 R, so [O III]/Hβ > 4.9 (2σ).
3. **[O III] velocity (measured).**
   - South strand: **V_helio = +12.5 ± 2.5 (stat) ± 0.7 (zero point) km/s**, which is V_LSR ≈ +0.2 km/s.
   - W strand: +11.7 ± 2.4 km/s. North: +10.6 ± 4.1 km/s.
   - SDSS-I/II alone gives +13.7 ± 3.6 and BOSS alone +11.8 ± 3.4. The delete-one-group jackknife range is +11.7…+13.5
     (σ_jk = 2.0).
   - The paper gives +23.8…+26.8 km/s, so it is **12–14 km/s higher (≈ 3σ)**.
   - SDSS [O III] matches the paper's own red lines (+11.7…+19.7) and WHAM Hα (V_helio ≈ +12). I do not reproduce the
     [O III] velocity offset that the paper uses to separate the filament from the DIG and to get its kinematic distance.
4. **MRS-N nebular spectra are not public.**
   - The DR12 v2.0 catalogue holds only the 83 star spectra of plate `MN074209N160100N01` (2024-02-01): 80 F-Std and 3 Obj.
   - The 2,684 nebular fibres were used as the plate's "sky" fibres (header NSKIES = 165 per spectrograph).
   - The same LAMOST group states in arXiv:2510.08950 that MRS-N fibre spectra "have not yet been publicly released".
   - So I could not re-measure the paper's lines independently.
5. **Absolute scale of the paper's relative intensities (derived).** The red-arm scale tied to WHAM is 11.8 ± 1.3 R per relative
   unit, i.e. the OH 6553 intensity. The blue-arm [N I]5198/OH6553 ratio is 0.72 in the UVES atlas (Hanuschik 2003) and 0.14 at APO
   (median, 16–84 % range 0.06–0.38). The paper's [O III] at P1–P3 then comes out at **2.8–3.2 R (APO ratio) or 14–16 R (UVES
   ratio)**. The SDSS direct values favour the low end.

## 2. Data and sources

| Source | What | URL / endpoint | Version |
|---|---|---|---|
| SDSS SkyServer DR17 SQL | SpecObjAll and PlateX in RA 108.5–124.5, Dec 4.5–29.5 (box padded by one plate radius): 112,792 fibres, 9,403 SKY, 138 plate-MJDs | `https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch` | DR17 |
| SDSS SAS DR17 spec-lite | 7,958 SKY-fibre spectra (SDSS-I/II run2d 26/103/104; BOSS/eBOSS v5_13_2) | `https://data.sdss.org/sas/dr17/sdss/spectro/redux/<run2d>/spectra/lite/<plate>/` and `https://data.sdss.org/sas/dr17/eboss/spectro/redux/v5_13_2/spectra/lite/<plate>/` | DR17 |
| SDSS-V DR19 | spAll SKY fibres (SkyServer DR19 SqlSearch; query text not kept in a script) and 1,134 spec-lite files in 10 fields | `https://dr19.sdss.org/sas/dr19/spectro/boss/redux/v6_1_3/spectra/lite/<field>/<mjd>/` (data.sdss5.org redirects to data.sdss.org, where DR19 paths 404) | DR19, v6_1_3 |
| LAMOST DR12 v2.0 (international release, end of Sept 2026) | LRS catalogue (191,196 spectra, RA 112.5–120.5, Dec 10.5–29); MRS catalogue (61,237 spectra / 40 plans, RA 111–122, Dec 6–28); 4,000 LRS + 1,871 MRS strand-sample spectra + 83 MN-plate spectra | SQL `https://www.lamost.org/dr12/v2.0/sql/q`; FITS `https://www.lamost.org/dr12/v2.0/spectrum/fits/<obsid>`, `.../medspectrum/fits/<obsid>`; docs `https://www.lamost.org/dr12/v2.0/catalogue`, `.../doc/mr-data-production-description` | DR12 v2.0 |
| UVES sky atlas | Hanuschik 2003, A&A 407, 1157, tables 6–7 | `https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2bA/407/1157/table6` (and table7). vizier.cds.unistra.fr returned "Page unavailable" | VizieR J/A+A/407/1157 |
| MRS-N data status | Rosette MRS-N paper, section "Data Accessibility" | `https://arxiv.org/html/2510.08950v1` | v1 |
| Solar template | CALSPEC `sun_mod_001.fits` (Kurucz model, vacuum) in `data/raw/calspec_sun/` (STScI CALSPEC archive; download log not kept) | STScI CALSPEC | — |
| Paper | Wu et al. 2026, arXiv:2610.13329 (Tables 1–4; Sect. II.2, III) | `sources/arxiv/2610.13329v1.*` | v1 |
| Teammate files (read only) | WHAM nearest beams and V_LSR−V_helio: `data/derived/spectra_ism_wham_results.json`; imaging SB→R: `notes/paper_review.md` §4.10; filament trace v1.0: `data/derived/silk_trace.json` | — | — |

Downloads total ≈ 3.0 GB: SDSS spec-lite 1.5 GB + 0.42 GB, LAMOST 1.0 GB, catalogues 0.1 GB. The largest single file is the
73 MB MRS box catalogue.

## 3. Methods

### 3.1 SDSS sky fibres (absolute sky spectra)
- **Why spec-lite instead of spFrame/spCFrame.** For a SKY fibre, the spec-lite COADD `flux` column is the residual after
  pipeline sky subtraction. The `sky` column is the spectrograph's B-spline sky model at that fibre. So `flux + sky` is the
  fibre's own total sky, and the group mean of `flux + sky` is the sky model that spFrame would give. This route also keeps
  per-fibre spatial resolution (2″ or 3″ fibres spread over a 3° plate), at about 1/20 of the download size.
- **Extraction** (`spectra_sdss_lamost_sdss_extract.py`). Windows 4800–5070, 5175–5225, 5560–5600, 6280–6320 and 6490–6790 Å are
  taken on the native 10⁻⁴ dex grid (69 km/s per pixel) into `data/derived/spectra_sdss_lamost_sdss_skyfibres.npz`.
- **Group-level absolute fit** (`spectra_sdss_lamost_sdss_fit.py`). The mean total sky of each plate–MJD–spectrograph group is
  fitted, per window, with a linear model:
  - solar template (CALSPEC Kurucz, broadened to the `wdisp` LSF) × quadratic, for zodiacal light, scattered light and starlight
  - plus a broad 4984 Å bump (high-pressure-sodium light pollution)
  - plus OH lines (red window)
  - plus Gaussians at the nebular lines for V = +20 km/s (19 km/s intrinsic σ); [O III]4959 is tied to 5007/2.98
  - the solar-template velocity is chosen on a grid

  Null lines at 5026 and 4880 Å give the template systematic: [O III]-window null median −0.20 R (as calibrated), Hβ null +0.03 R.
- **Per-fibre values.** Each fibre's residual (`flux`) is fitted with a local baseline plus the same Gaussians. Absolute fibre
  intensity = group value + Δ_fibre − ⟨Δ⟩_group.
- **Diffuse-light correction.** SDSS fluxes are calibrated for point sources, so a uniform source is over-estimated by 1/f_psf.
  I multiply by f_psf, the point-source fraction in a 2″ or 3″ fibre for a Gaussian PSF at the plate's seeing (median 0.78,
  systematic ≈ 20 %).
- **Calibration check (measured).**
  - Corrected sky continuum at 5040 Å: median 21.75 AB mag/arcsec², which matches a dark APO sky (V ≈ 21.7–21.9).
  - [O I]5577: median 286 R, a normal value.
  - [N I]5198 component: 1.03 R (16–84 %: 0.53–2.57). OH 6553: 7.8 R (5.2–11.8).
- **Topocentric sign (fixed in this session).** The spec-lite `HELIO_RV` equals −v_corr, where v_corr is the astropy heliocentric
  correction. I checked this on 12 plates (|difference| < 1.2 km/s). The [O I]5577 centroid follows λ₀(1 − HELIO_RV/c): a fit
  of offset against HELIO_RV gives slope −2.03 and r = −0.986 under the opposite assumption.
  - The earlier version of the fit placed airglow, OH lines and the solar template with the wrong sign.
  - This affected only the group-level airglow and OH modelling. It did not affect the nebular-line residuals or the ON/OFF
    differences.
  - All products were regenerated after the fix. Old and new strand amplitudes agree.
  - The LAMOST header `HELIO_RV` follows the same convention: +8.10 km/s for the MN plate, against astropy −8.1 km/s.
- **Efficiency fix.** `np.load` of the npz is now read into memory once. Re-reading the NpzFile per fibre had made the fit take
  more than 15 minutes; it now takes 8 s.

### 3.2 Strand contrast, ON−OFF stacks, velocities
- **Geometry.** The geometry teammate's trace v1.0 (`silk_trace.json`) gives an E spine (Dec 12–28.5) and a W strand
  (Dec 15–20). I use signed distances to the nearest strand (`spectra_sdss_lamost_trace.py`).
- **Strand fit** (`spectra_sdss_lamost_sdss_maps.py`). I_f = b_group + A_E·G(d_E) + A_W·G(d_W), with Gaussian σ = 0.12° and
  5σ clipping. The per-group offsets absorb all plate-level systematics.
- **ON−OFF spectra** (`..._sdss_stack.py`, `..._sdss_substacks.py`). ON = fibres within 0.15° of a strand; OFF = fibres 0.45–2.0°
  away on the same plate–spectrograph. Group differences are combined with inverse-variance weights and 5σ per-pixel clipping.
  Gaussians are fitted with free common velocity for the [O III] and Hα+[N II] sets.
- **Wavelength zero point.** I fit the [O I]5577 centroid (blue camera) and [O I]6300/6364 (red camera) in every group's total
  sky:
  - blue median −0.2 km/s (16–84 %: −3.8…+2.1)
  - red median −5.3 km/s (−9.6…−1.7)

  Each stack is corrected by the ON-weighted mean of its groups.
- **Null region.** Fibres 1.0–1.15° from the strand compared with 1.3–2.5° give [O III] 0.010 ± 0.014 R and Hβ 0.003 ± 0.015 R.
  They also give Hα 0.080 ± 0.024 R and [S II] 0.066 R. So the low-ionization excesses carry a **±0.08 R systematic**.
- **Robustness of the velocity** (`..._sdss_vel_jackknife.py`): instrument split and delete-one-group jackknife.

### 3.3 LAMOST DR12 v2.0 stellar residuals
- **How the residual arises.** LAMOST subtracts a super-sky from each spectrograph's sky fibres: about 20 in normal plates, about
  165 nebular fibres on the MN plate. A star therefore carries I(star) − ⟨I⟩ over that spectrograph's sky fibres.
- **Method** (`spectra_sdss_lamost_lamost_resid.py`). Each spectrum is continuum-normalised. A data-driven template is built from
  the median of the 25 most similar OFF stars in the stellar rest frame. I fit template × linear + Gaussians, convert EW to R
  using PS1 g/r (LRS) or Gaia G (MRS) and the fibre aperture fraction (3.3″), and stack.
- **Samples.** Random subsets of STAR/snrg > 10: LRS has 2,500 ON stars (|d| < 0.2°, Dec 12–28.2) and 1,500 OFF (0.6–2.5°). MRS has
  971 ON (< 0.25°) and 900 OFF on 3 plates. The selection code was not kept as a script. I re-checked the criteria: recomputed
  distances match exactly and the S/N distribution matches the parent sample.
- **Fits.** Stack fits are in `spectra_sdss_lamost_lamost_stackfit.py`. Fixed-width per-plate MRS [O III] fits are in
  `spectra_sdss_lamost_mrs_onoff_velocity.py`.
- **MN plate.** Its 83 public stars were analysed with templates from the other plates (`spectra_sdss_lamost_mnplate_sample.py`,
  mode `mrsmn`).

### 3.4 Absolute scale of the paper's Table 2 (`spectra_sdss_lamost_paper_compare.py`)
- **Red arm.** The paper's relative Hα is compared with the WHAM-SS total Hα of the nearest beams (teammate file). The scale is
  10.6, 11.9, 11.0 and 13.6 R per unit for P1–P4, mean **11.8 ± 1.3 R**, which equals I(OH 6553) on that night.
- **Blue arm.** The reference line is [N I] 5197.9, which R ≈ 7500 resolves from 5200.3. The ratio I(5198)/I(OH 6553) is:
  - UVES atlas: 8.63 R / 11.97 R = 0.72. Here VizieR `Flux` is the peak flux density, so integrated = Flux × FWHM × 1.0645. This
    gives 313 R for [O I]5577, a sensible value.
  - APO (this work): median 0.144 (0.059–0.383 for 16–84 %).

## 4. Results

### 4.1 [O III] surface brightness

| Quantity | Value | Kind |
|---|---|---|
| W strand amplitude, Dec 15–20 | 1.40 ± 0.04 R (1.46 ± 0.05 with only Dec 15–20 data) | measured |
| E strand amplitude by Dec bin: 14–15.5 / 15.5–17 / 17–18.5 / 18.5–20 / 20–23 / 23–29 | 1.49 ± 0.18 / 1.17 ± 0.12 / 0.28 ± 0.08 / 0.23 ± 0.08 / 0.78 ± 0.04 / 0.43 ± 0.03 R | measured |
| ON−OFF, Dec 15–20: W strand vs E strand | 1.19 ± 0.03 R vs 0.14 ± 0.04 R | measured |
| Off-strand diffuse [O III] (profiles at 1.5–3°) | 0.32–0.38 R (± 0.2 R zero point) | measured |
| Mean within 0.3° of P1 / P2 / P3 / P4 | 1.42 ± 0.14 / 1.80 ± 0.20 / 1.70 ± 0.37 / 0.98 ± 0.24 R (P4 from a single plate) | measured |
| Paper imaging SB → R (P1–P3) | 12.8 / 9.9 / 10.1 R | literature (paper-review teammate) |
| Paper MRS-N relative × APO sky ratio × 11.8 R | 2.8 / 3.1 / 3.2 R (16–84 %: 1.2–8.6) | derived |
| Paper MRS-N relative × UVES ratio × 11.8 R | 14.2 / 15.5 / 16.1 R | derived |

Notes on the table:
- The strand amplitude and P-region averages are 1.4–1.8 R. Adding the 0.35 R background, the peak total is ≈ 1.8–2.2 R, and
  single fibres reach ≤ 3.5 R.
- The imaging-based 10–13 R would need a calibration error of about ×5. To reconcile the paper's own spectra with SDSS, the [N I]
  5198 line at Xinglong on 2024-02-01 must have been ≈ 0.5–1 R. That is within the APO distribution.
- In SDSS the W strand (through P2) carries most of the [O III] at Dec 15–20. This is useful for the geometry teammate.

### 4.2 Line ratios on the strand (ON−OFF, diffuse-corrected R)

| Region | [O III]5007 | Hβ | Hα | [N II]6583 | [S II]6716+31 | [O III]/Hβ | [O III]/Hα | [N II]/Hα | [S II]/Hα |
|---|---|---|---|---|---|---|---|---|---|
| South, Dec < 20 (141 fibres, 26 groups) | 0.79 ± 0.02 | 0.017 ± 0.020 | 0.17 ± 0.03 | 0.16 ± 0.03 | 0.10 ± 0.03 | **> 13.8 (2σ)** | 4.7 ± 1.0 | 0.96 ± 0.25 | 0.62 ± 0.23 |
| W strand, Dec 15–20 (55, 20) | 1.19 ± 0.03 | 0.055 ± 0.028 | 0.29 ± 0.05 | 0.24 ± 0.04 | 0.19 ± 0.05 | 22 ± 11 (> 10.8) | 4.2 ± 0.7 | 0.82 ± 0.19 | 0.67 ± 0.20 |
| North, Dec ≥ 20 (124, 35) | 0.43 ± 0.02 | 0.058 ± 0.019 | 0.18 ± 0.03 | 0.33 ± 0.03 | 0.31 ± 0.04 | 7.4 ± 2.5 | 2.4 ± 0.4 | 1.85 ± 0.36 | 1.77 ± 0.37 |
| All (267, 58) | 0.58 ± 0.02 | 0.040 ± 0.014 | 0.17 ± 0.02 | 0.23 ± 0.02 | 0.19 ± 0.03 | 14.7 ± 5.1 (> 8.6) | 3.4 ± 0.5 | 1.36 ± 0.21 | 1.11 ± 0.20 |
| LAMOST LRS ON−OFF (2,500/1,500 stars) | 0.88 ± 0.13 | −0.08 ± 0.09 | 0.17 ± 0.08 | 0.30 ± 0.03 | 0.25 ± 0.04 | > 4.9 (2σ) | — | (stellar Hα residuals) | — |

Notes on the table:
- Hα, [N II] and [S II] errors exclude the ±0.08 R null-test systematic, so the low-ionization ratios are uncertain by a further
  ~50 %.
- [O III]5007/4959 is 2.66 ± 0.19 (south) and 3.07 ± 0.23 (W strand), consistent with the atomic 2.98.
- [S II]6716/6731 is 1.44 ± 0.40 (all). This points to the low-density regime but constrains n_e only weakly. The paper's
  1.03–1.19 refers to the total sightline, which is dominated by the WIM.
- If the Hα excess is recombination, the expected Hβ in the south is 0.17/3.86 (case B in photons) = 0.044 R. The measured
  0.017 ± 0.020 R is within 1.4σ of that. So the strand's small low-ionization excess is consistent with ordinary recombination
  gas, but the **[O III] has no Balmer counterpart**.

### 4.3 Velocities (heliocentric; V_LSR = V_helio − 12.2 km/s)

| Measurement | V_helio (km/s) | Kind |
|---|---|---|
| SDSS [O III] strand, Dec < 20 | **+12.5 ± 2.5 stat ± 0.7 zp** (jackknife σ 2.0; SDSS-I/II +13.7 ± 3.6, BOSS +11.8 ± 3.4) | measured |
| SDSS [O III] W strand / north | +11.7 ± 2.4 / +10.6 ± 4.1 | measured |
| SDSS Hα+[N II] strand (all) | +24.9 ± 5.8 after red zp (−5.8); +19.1 without the red zp | measured (red zp less certain) |
| Paper [O III] P1/P2/P3 | +26.8 / +26.8 / +23.8 (± 3–4) | literature |
| Paper Hα P1–P4; [N II], [S II] | +16.3…+19.7; +11.7…+13.3 | literature |
| WHAM Hα, 4 nearest beams | +11.9…+12.7 (V_LSR −0.6…+0.4) | literature (ISM teammate file) |
| LAMOST DR12 MRS super-sky [O III] (OFF-star dips; fixed width) | +3.8 ± 12.8, +8.1 ± 5.5, +3.9 ± 4.5 (3 plates, v_corr from −17 to +17) | measured, low confidence |
| LAMOST MRS ON−OFF per plate | −55.6 ± 7.0, +23.4 ± 4.4, −0.1 ± 6.4 (combined +7.7 ± 4.6) | measured; systematics dominated |
| LAMOST MN plate (paper's plate), 77 stars | +0.46 ± 0.07 R at −3.7 ± 4.4 km/s; the positive sign points to template contamination | measured; not interpretable |

SDSS-minus-paper [O III] difference: 13.3 ± 3.3 km/s if the paper's three regions are independent errors, ≈ 3σ conservatively.

## 5. What the spectra imply for the SNR verdict

**Measured facts.**
- The southern Silk strand emits [O III] at ~1–2 R with [O III]/Hβ > 14 and [O III]/Hα ≈ 4–5.
- Its [O III] velocity equals the local WIM velocity (V_LSR ≈ 0).
- The weak low-ionization excess (~0.2 R ± 0.08) has [N II]/Hα ≈ 1 and [S II]/Hα ≈ 0.6.
- North of Dec 20°, [O III] weakens and [O III]/Hβ falls to ~7, with [N II]/Hα and [S II]/Hα of 1.5–2.

**Hypotheses, ordered by how well the data support them:**
1. **Collisionally excited [O III] at a hot-gas boundary, with suppressed recombination (favoured for the south).**
   - [O III]/Hβ > 14 is beyond photoionized WIM or DIG ([O III]/Hβ ≲ 1) and normal H II regions.
   - It is also beyond complete, steady radiative shocks. Shock-only MAPPINGS grids (Allen et al. 2008, ApJS 178, 20) give
     [O III]/Hβ typically ≲ 10 for v_s ~ 100–300 km/s; this is quoted from memory of the grids and was not recomputed here.
   - It is the signature of an *incomplete* shock: a front younger than its cooling time or truncated before the Balmer
     recombination zone, as in Cygnus Loop [O III] filaments (Raymond et al. 1988, ApJ 324, 869).
   - It also matches a conductive or turbulent-mixing interface between ~10⁶ K gas and cooler gas.
   - Either way a hot, X-ray-emitting volume should lie next to the strand. This supports a "shock/bubble boundary" in a stricter
     sense than the paper's [S II]/Hα criterion, which the paper-review teammate shows is satisfied by the WIM itself.
   - V_LSR ≈ 0 is what an edge-on, limb-brightened sheet shows. So the kinematics are consistent with a shell limb but give **no
     kinematic distance**.
2. **A different excitation regime in the north (tentative).** [O III]/Hβ ≈ 7 with strong [N II] and [S II] looks like a more
   complete radiative shock, perhaps where the front meets denser gas. It could also be an unrelated low-ionization structure on
   the same sightlines, such as the northern Hα patch or the 17° NSNS arc. This needs deeper spectra. A ±0.08 R systematic in
   Hα and [S II] is known.
3. **Photoionization by a hot star or a PN-like source (disfavoured).** It would need [O III]/Hβ ≳ 14 over a ~10° sheet with
   negligible Balmer emission. No ionizing source is required by the data.
4. **[O III] kinematically distinct from the DIG (not supported).** The paper's +10–15 km/s [O III]−Hα offset is not reproduced
   (SDSS gives 0 ± 3 km/s relative to WHAM). Possible causes include:
   - a blue-arm zero-point problem in the MRS-N secondary wavelength calibration, which has only the [N I]5198/5200 doublet as a
     sky line in 4950–5350 Å
   - a heliocentric-correction slip in one arm: 2|v_corr| = 16.2 km/s for this plate

   These are hypotheses and cannot be tested without the MRS-N nebular spectra. The paper's red-arm velocities agree with WHAM.

**Consequences for the paper's numbers.**
- At ~1.5–2 R the emission measure and density arguments of the paper-review teammate relax by ×5–7.
- The energy and age estimates rest on a kinematic distance from V_LSR([O III]) = +11–16 km/s. With V_LSR ≈ 0 near the
  anticentre, that distance is unconstrained, and so are E ≈ 2×10⁵¹ erg and t ≈ 6.8×10⁴ yr.

## 6. Figures

| File | Caption |
|---|---|
| `figures/spectra_sdss_lamost_sdss_oiii_map.png` | Absolute [O III]5007 per SDSS/BOSS/SDSS-V sky fibre (diffuse-corrected R) over the full box and the LAMOST field, with trace v1.0, the Hα patch and P1–P4. |
| `figures/spectra_sdss_lamost_sdss_strand_dec.png` | Within-plate strand amplitudes (Gaussian σ = 0.12°) of [O III], Hβ, Hα, [N II] and [S II] by Dec bin; diamonds = W strand. |
| `figures/spectra_sdss_lamost_sdss_profiles.png` | Absolute line intensity against distance from the nearest strand (all Dec), with sampling histogram. Mixes plates; absolute Hβ and Hα levels include template and geocoronal contributions; use the ON−OFF numbers. |
| `figures/spectra_sdss_lamost_sdss_diffstack.png` | ON-strand minus same-plate OFF sky-fibre spectrum (267 fibres, 58 groups): [O III] 0.58 R, Hβ 0.04 R. |
| `figures/spectra_sdss_lamost_sdss_substacks.png` | The same split into Dec < 20° (south), Dec ≥ 20° (north) and a null region. |
| `figures/spectra_sdss_lamost_paper_compare.png` | [O III] at P1–P3 by method (imaging, MRS-N relative × sky-line scales, SDSS), and V_helio from the paper, WHAM, SDSS and LAMOST MRS. |
| `figures/spectra_sdss_lamost_lamost_stacks.png` | LAMOST DR12 LRS and MRS stellar-residual stacks ON and OFF strand. LRS [O III] ON−OFF 0.88 ± 0.13 R; MRS dominated by sky-subtraction mismatch. |
| `figures/spectra_sdss_lamost_quicklook.png` | Exploratory: mean total sky of sky fibres for plates 4492 (116.5, +17.6), 4497 (116.2, +15.8), 11102 (117.2, +17.0), 4843 and 2940 (off-field), showing [O III] emission, Hβ emission in the solar absorption, and the 4984 Å HPS bump. |

## 7. Files

- **Scripts** (`scripts/`): `spectra_sdss_lamost_{query_sdss, get_sdss_sky, get_dr19_sky, sdss_extract, sdss_fit, sdss_maps,
  sdss_stack, sdss_substacks, sdss_vel_jackknife, get_lamost, lamost_resid, lamost_stackfit, mrs_onoff_velocity,
  mnplate_sample, paper_compare, explore, lines, trace}.py`.
- **Run order:** query → get → extract → fit → maps → stack → substacks → vel_jackknife; then lamost_resid lrs|mrs|mrsmn →
  lamost_stackfit → mrs_onoff_velocity → paper_compare.
- **Derived data** (`data/derived/spectra_sdss_lamost_*`):
  - `sdss_skyfibres.npz`: fibre windows
  - `sdss_groups.csv`: 285 groups, absolute fits, airglow
  - `sdss_fibres.csv`: 8,732 fibres
  - `sdss_onoff.json`: strand fits
  - `sdss_profiles.csv`
  - `sdss_diffstack.{json,npz}`
  - `sdss_substacks.json`
  - `sdss_zeropoints.csv`
  - `sdss_vel_jackknife.json`
  - `sdss_calib.json`
  - `sdss_globalstack.npy`
  - `lrs|mrs|mrsmn_resid_{stars.csv,stacks.json}`
  - `lamost_stackfit.json`
  - `mrs_onoff_velocity.json`
  - `mnplate_sample.csv`
  - `paper_compare.json`
- **Raw data:** `data/raw/sdss_dr17_skyserver/`, `sdss_dr17_speclite/`, `sdss_dr19_skyserver/`, `sdss_dr19_speclite/`,
  `lamost_dr12/`, `lamost_dr12_docs/`, `lamost_mrsn_literature/`, `vizier_hanuschik2003/`.

## 8. Limitations
- **Diffuse-light factor f_psf.** It is a Gaussian-PSF approximation with a ≈ 20 % systematic. The sky-continuum check supports
  it. Without it, every SDSS intensity rises by ×1.28 (median), which does not change the conclusions.
- **Absolute zero points.** Group-level absolute values depend on the solar-template continuum model. The [O III] zero point is
  ±0.2 R. The absolute Hβ (~1.6 R) and Hα (~5–7 R, including geocorona) levels are not reliable and are not used. All ratios use
  ON−OFF differences.
- **Sampling.** ON/OFF classes depend on trace v1.0. The SDSS strand sampling is sparse (52–141 ON fibres per region). Fibres
  do not sample P1–P3 at the paper's positions exactly.
- **Resolution.** At R ≈ 2000 (σ_LSF ≈ 69 km/s) line widths cannot be measured. The centroids assume symmetric LSFs. The red-camera
  zero point (−5.3 km/s) is less certain than the blue one.
- **LAMOST stellar residuals.** These are limited by template mismatch: Hβ and Hα sit in stellar absorption lines, and MRS shows
  ±1 R/Å features at ±50–80 km/s. MRS per-plate velocities scatter by ±25 km/s and are not used for conclusions. The MN-plate
  residual (+0.46 R) has the wrong sign for a super-sky dip and is treated as an artefact.
- **Sky-line ratios.** The ratio converting the paper's blue-arm intensities varies by ×5 (UVES 2000–01 vs APO medians). Only a
  range can be given.

## 9. Data gaps
- **MRS-N nebular spectra and per-fibre line catalogue** for `MN074209N160100N01`: not public. DR12 v2.0 has no nebular
  catalogue (catalogue page lists only stellar LRS/MRS catalogues), and arXiv:2510.08950 says they "have not yet been publicly
  released". So I could not independently re-measure the paper's [O III] velocity or widths.
- **No LAMOST sky-fibre spectra are released.** DR12 v2.0 fibre types in the box are only `Obj` and `F-Std`.
- **`dr10.lamost.org`** is unreachable. Not needed, since DR12 v2.0 at www.lamost.org was used.
- **SDSS spFrame/spCFrame** were not downloaded; the spec-lite sky-fibre route gives the same sky model. The SDSS object-spectrum
  residual method (i) was not run because the sky fibres measure the nebula directly. SPARCL was not used.
- **SDSS-V DR19.** 19 of 26 field-MJD groups were rejected as bright-sky (moon or twilight) or for poor [O III] fits, against
  10 of 259 DR17 groups.
- **Hanuschik 2003 on vizier.cds.unistra.fr** returned "Page unavailable"; the vizier.cfa.harvard.edu mirror was used.
- **No Hβ in MRS-N.** The paper notes this. SDSS now supplies the [O III]/Hβ limit.

## 10. Notes for teammates
- **Geometry.** In SDSS the W strand (through P2) is far brighter than the E strand at Dec 15–20 (1.19 vs 0.14 R ON−OFF).
  [O III] continues north of Dec 20 at 0.4–0.8 R.
- **Paper review.** Use ~1.5–2 R for [O III] on the strand, not 10–13 R. [O III]/Hβ > 14 (2σ). The [O III] velocity is
  +12.5 ± 2.6 km/s, which removes the [O III]−Hα offset.
- **Anyone using SDSS sky positions.** `HELIO_RV` in SDSS and LAMOST headers is −v_corr (astropy), so topocentric features sit at
  λ₀(1 − HELIO_RV/c).
