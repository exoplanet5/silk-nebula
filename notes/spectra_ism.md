# spectra_ism: ionized-gas kinematics and interstellar absorption toward the Silk Nebula

Topic tag `spectra_ism`. Access date for all sources is 2026-11-16. Every network request (URL, parameters, HTTP status, bytes, UTC
time) is logged in `data/raw/spectra_ism_fetch_log.jsonl`. Velocities are heliocentric (V_helio) unless marked LSR. V_LSR uses
the standard solar motion (20 km/s toward 18h, +30°), the convention of WHAM and Wu et al. 2026. In this direction
V_LSR − V_helio = −11.9…−12.6 km/s at P1–P4.

## 1. Questions and short answers

| Question | Answer (measured unless stated) |
|---|---|
| Does WHAM show Hα at the [O III] velocity (V_LSR ≈ +13.6)? | **No.** On-filament minus off-filament WHAM stack (23 matched pairs): Hα fixed at V_LSR +13.6 km/s with an [O III]-like width = **0.010 ± 0.026 R** (bootstrap), so < 0.063 R at 2σ. The beams have an NSNS [O III] excess of 0.66 ± 0.09 R, which gives **Hα/[O III] < 0.10 (2σ)** at the [O III] velocity. |
| Which velocity does WHAM Hα have? | One main component at V_LSR = +1.0 ± 0.2, −0.4 ± 0.1, +0.9 ± 0.4 and +0.3 ± 0.2 km/s (P1–P4), with intrinsic FWHM 30–32 km/s and 2.9–4.4 R. This matches the LAMOST [N II]/[S II] (V_LSR −0.9…+1.7) and is the same on and off the filament: it is the WIM/DIG. |
| LAMOST red lines in Rayleighs? | Scale k = 12.1 ± 1.5 R per sky-line unit, from WHAM Hα ÷ LAMOST relative Hα; the four regions agree to 13 %. P1: Hα 3.6 ± 0.5, [N II]6583 2.9 ± 0.4, [S II]6716+6731 ≈ 3.3 R. The other regions are in §3.1. |
| Interstellar absorber at +25 km/s? | **None was found with a clear detection.** The Ca II / Na I components cluster at +6…+15 km/s (the local gas, the same as WHAM and the red lines) and at +18…+20 km/s (LIC/Aur clouds, already present at 63 pc). The only possible components at +23…+28 km/s are weak, single-line and tentative (11 CMi at 92 pc; HD 66552 at 121 pc; HD 61038 at 149 pc). |
| Shock or high-velocity gas? | **One candidate.** HD 61717 (Gaia d = 272 pc, 0.97° W of the spine at Dec +22.3) has a narrow Ca II component at **V_helio = −17.0 km/s** (V_LSR −28), about 30 km/s blueward of the local gas. It has N(Ca II)/N(Na I) = **22 ± 2**, a lower limit because Ca II is mildly saturated. This is a grain-depleted (Routly–Spitzer) signature typical of shocked or intermediate-velocity gas. It is absent toward HD 64705 (280 pc, 2.6° E) and HD 61038 (149 pc, 1.9° W). |
| UV hot gas? | COS G130M of WD J074735.98+210635.83 (Gaia d = 99.5 ± 0.5 pc, **0.16° from the spine**): Si IV and N V are **not detected**. At V_helio +15…+40 the 3σ limits are log N(Si IV) < 12.2 and log N(N V) < 13.0. Only normal local gas (C II, Si II, S II, C II*) appears, at V ≈ +5…+15. No hot interface lies within 100 pc on this sightline. |
| Faraday rotation across the filament? | In the north half (spine Dec +20…+28) the NVSS RMs change sign across the filament: east minus west (2–6° from the spine) = **+36.8 ± 4.6 rad m⁻²**. 70 % of RMs are negative on the west side and 16 % on the east. The step lies at +2.0° (−0.6…+2.0°) E of the spine. In the south half (Dec +12…+20) there is no step: +5.3 ± 7.3 rad m⁻². |

## 2. Methods and sources

### 2.1 WHAM Hα (scripts/spectra_ism_wham.py, spectra_ism_wham_figures.py)
- **Data:** WHAM-SS DR1 v161116 (Haffner et al. 2003, 2016-11-16 build, FITS dated 2017-09-12). The file is
  `data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits` (91,200,960 B), downloaded by a teammate.
  - It is byte-identical to the file whampy uses, https://uwmadison.box.com/shared/static/4kccrw2bgad7muss3z2po7rezklenxxz.fits (https://github.com/Deech08/whampy, `skySurvey.py`).
  - The old page www.astro.wisc.edu/wham-site/ returns 404.
  - DR1 is Hα only. I found no public WHAM [S II], [N II] or [O III] data for this field (l ≈ 204°, b ≈ +19°). The pointed studies Haffner 1999 and Madsen 2006 target other directions. I did not verify this exhaustively.
- **Beams:** 1550 beams in RA 100–133°, Dec −4…+38°. Each beam has 0–3 Gaussian components, chosen by BIC. Intrinsic widths are the quadrature difference with the 12 km/s instrumental FWHM.
- **Filament coordinates:** taken from the geometry teammate's spine (`data/derived/silk_trace.json` v1.0).
- **On − off stack:** beams with |d_perp| < 0.5° minus the mean of east and west beams at 1–3°, matched within 0.8° in Dec. There are 23 pairs at Dec +12.8…+27.6. Errors come from bootstrapping over pairs.
- **Regression:** each channel regressed on the beam-averaged NSNS [O III] (teammate map `imaging_optical_nsns_maps.npz`, in R), plus a quadratic background surface. 288 beams; the northern Hα patch is excluded.
- **Velocity conversion:** LAMOST velocities were converted to LSR (`vlsr_recomputed_from_vhelio` in the JSON). The paper's own V_LSR column differs from my recomputation by up to 0.8 km/s for [O III] and up to 2.8 km/s for the red lines. For example, P4 Hα is +7.3 recomputed versus +4.5 in Table 4.

### 2.2 LISM cloud velocities
- Redfield & Linsky (2008) kinematic calculator, http://lism.wesleyan.edu/cgi-bin/dynlism.cgi (POST, Galactic coordinates). Raw HTML is in `data/raw/spectra_ism_lism/`.
- Toward P1 the sightline traverses the LIC (V_helio = +19.1 ± 1.1) and the Gem cloud (+33.9 ± 1.1), and passes near Aur (+20.2 ± 1.1). P2, P3, the north end and HD 61717 give values within 1.5 km/s of these.
- The [O III] velocity (+23.8…+26.8) therefore falls between two local (< 15 pc) clouds. Any absorption near +20 or +34 km/s is expected even toward the nearest stars.

### 2.3 Literature absorption catalogues (scripts/spectra_ism_absorption_catalogs.py, spectra_ism_vizier_parse.py)
- **Mirror:** VizieR CfA, https://vizier.cfa.harvard.edu/viz-bin/asu-tsv. vizier.cds.unistra.fr asu-tsv returned a CDS "page unavailable" page.
- **Catalogues used:**
  - J/A+A/510/A54 (Welsh+2010, Na I + Ca II, 1851 rows)
  - J/A+A/411/447 (Lallement+2003). Its W(D2) column is labelled "0.1 nm" but holds mÅ.
  - J/A+A/346/785 (Sfeir+1999)
  - J/A+A/507/833 (Megier+2009)
  - J/MNRAS/367/1478 (Hunter+2006)
  - J/ApJS/199/8 (Gudennavar+2012 UV/optical column compilation)
- **Not on VizieR:** J/ApJS/139/439, J/ApJ/673/283 (Redfield & Linsky component tables), J/ApJ/798/35, J/A+A/674/A40, J/ApJS/106/533 and J/ApJS/93/533 all returned "Table or Catalog not found".
- **Distances:** SIMBAD parallaxes from CfA SIMBAD TAP (`simbad_plx_abs_catalog_stars.csv`).
- **Result:** 22 catalogue stars lie within 12° of the spine, all within 300 pc. Output is `data/derived/spectra_ism_abs_catalog_stars.csv`. The catalogues give totals only, not component velocities.

### 2.4 Archival high-resolution spectra
- **Inventory** (`spectra_ism_archives.py`):
  - SIMBAD TAP: 1370 O/B/A stars with V < 10 in RA 104–128°, Dec +4…+34°.
  - ESO ObsCore (https://archive.eso.org/tap_obs/sync): 9547 phase-3 UVES/FEROS/HARPS/ESPRESSO spectra in the box.
  - ELODIE and SOPHIE region listings (http://atlas.obs-hp.fr/).
  - Result: 34 early-type stars with d < 1.5 kpc have spectra. All 6 within 3° of the spine were downloaded.
- **Downloads** (`spectra_ism_get_spectra.py`, `data/raw/spectra_ism_spectra/manifest.csv`): 40 products for 29 stars, 97 MB.
  - 19 ESO products: https://dataportal.eso.org/dataPortal/file/<ADP id>
  - 11 SOPHIE s1d
  - 10 ELODIE s1d
- **Frames** (`spectra_ism_frame_test.py`): checked empirically by cross-correlating telluric O₂ γ-band lines against a topocentric UVES spectrum (`data/derived/spectra_ism_frame_test.csv`).
  - UVES phase-3: topocentric (header ESO QC VRAD BARYCOR applied).
  - HARPS and SOPHIE s1d: barycentric. The telluric shift equals BERV within 0.1 km/s.
  - **ELODIE s1d: topocentric**. The shift is 0.0 ± 0.4 km/s for BERV from −29 to +20 km/s.
  - 11 CMi: BERV computed with astropy, because the header keyword is absent.
- **Telluric correction of Na D** (`spectra_ism_telluric.py`): median topocentric optical-depth template from all *other* spectra, each smoothed to the target's resolution, scaled outside the target's interstellar window and divided out.
- **Fits** (`spectra_ism_absorption_measure.py`):
  - Continuum from a running median (70 km/s window), so narrow lines survive and broad stellar lines divide out.
  - 0–3 Gaussian components in −60…+80 km/s, chosen by BIC.
  - Forced fit at +25 km/s with b = 5 km/s.
  - AOD columns N = 1.13×10¹⁷ W(mÅ)/(f λ²).
- **Vetting** (`spectra_ism_absorption_curate.py`): every profile was inspected. Lines are flagged good/ok/tentative/bad using SIMBAD RV and v sin i (`simbad_{basic,rot}_spectra_stars.csv`). Components within 10 km/s of the stellar RV in sharp-lined stars are flagged "possible stellar blend": HD 53588 (RV +16.4), HD 64705 (+31.5) and HD 67959.
- **Unusable:** HD 67959, BR CMi (composite spectrum), BD+17 1832, HD 57930 and HD 59605.
- **Outputs:** `spectra_ism_abs_components.csv` (all fits), `spectra_ism_abs_curated.csv/.json` (vetted), `spectra_ism_abs_stars.csv`, `spectra_ism_abs_telluric.csv` and `spectra_ism_abs_profiles.npz`.

### 2.5 UV (scripts/spectra_ism_mast_uv.py, spectra_ism_cos_wd.py)
- **Query:** MAST CAOM ObsCore TAP, https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync. vao.stsci.edu/CAOMTAP returned 404.
  - An unfiltered query hit the 100,000-row cap because of solar-system targets.
  - The final queries ran FUSE/IUE/HUT and HST COS/STIS/GHRS/FOS separately, with λ_min < 320 nm, and dropped moving targets by name.
  - Result: 492 fixed targets, listed in `data/derived/spectra_ism_mast_uv_targets.csv`.
- **Analysis:** COS G130M x1dsum `lfac2c010` (program 17420), https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/lfac2c010_x1dsum.fits. Gaia DR3 for the WD came from VizieR I/355/gaiadr3.

### 2.6 Rotation and dispersion measures (scripts/spectra_ism_rm_psr.py)
- **Data:** Taylor et al. (2009) NVSS RMs (J/ApJ/702/1230; 1494 sources in a 20° cone) and the ATNF catalogue snapshot on VizieR (B/psr).
- **Method:** a quadratic foreground surface is fitted. Cross-filament profiles, east/west steps and a bootstrapped on−off comparison are computed in filament coordinates.

## 3. Results

### 3.1 WHAM Hα: kinematics and the [O III] velocity
- **Line widths and intensities.** WHAM's main Hα component is at the LAMOST [N II]/[S II] velocity and width, not at the LAMOST Hα velocity.
  - WHAM P1–P4: V_LSR −0.4…+1.0 km/s, FWHM_int 29.6–31.9 km/s, I = 3.17 / 4.41 / 2.86 / 3.40 R (nearest beams).
  - LAMOST red lines, recomputed: [N II] V_LSR −0.9…+0.9 and Hα +4.2…+7.3.
  - The LAMOST Hα centroid is therefore **4–7 km/s redder than WHAM Hα and than LAMOST's own [N II]/[S II]**.
  - Off-filament P4 behaves the same as P1–P3. This supports the paper's statement that the red lines are foreground/background DIG.
- **Hypothesis (not tested on spectra): the LAMOST Hα offset is a sky-subtraction artefact.**
  - At Xinglong on 2024-02-01 (transit ~15:15 UT) the heliocentric correction is −8.3 km/s. Geocoronal Hα therefore lies at V_helio ≈ −8, V_LSR ≈ −20 km/s.
  - Over-subtracting a sky spectrum that contains geocoronal Hα and the blue half of the DIG would push the residual Hα centroid redward.
  - Testing this needs the per-fibre MRS-N spectra. dr10.lamost.org is not reachable.
  - The paper-review teammate raised the same point independently in `notes/paper_review.md` (row 18).
- **No Hα at the [O III] velocity.**
  - On−off stack: 0.010 ± 0.026 R at V_LSR +13.6 (FWHM_obs 47.8 km/s).
  - Window integrals of the difference spectrum: +9…+20 km/s gives −0.003 ± 0.012 R; +20…+40 gives −0.001 ± 0.009 R.
  - Regression slope (Hα in +9…+20 km/s per R of NSNS [O III]): 0.020 ± 0.037.
  - With a 0.66 R [O III] beam excess this gives **Hα/[O III] < 0.10 (2σ)**.
  - If the paper's 10–13 R [O III] is right rather than the NSNS calibration (§6), the limit is about 4× tighter (< 0.03).
- **Implication for the emission physics (derived).** A recombination-line ratio this low cannot come from photoionized gas at ~10⁴ K, which gives [O III]/Hα of a few at most.
  - It needs collisionally excited O²⁺ in gas hotter than ~2–3×10⁴ K that has not yet cooled and recombined. That points to an incomplete or non-radiative-to-radiative transition shock, or a conductive/turbulent mixing layer on a hot bubble.
  - This is the regime of the standard shock models (Raymond 1979, ApJS 39, 1; Hartigan, Raymond & Hartmann 1987, ApJ 316, 323). It favours a young or truncated shock front over a classic radiative SNR filament or plain DIG.
  - Emission measure of the [O III] gas at ~10⁴ K: EM < 0.17 pc cm⁻⁶ (beam average, T₄ = 1), or about 0.5 pc cm⁻⁶ at the filament's true width. For path lengths of 1–10 pc this means n_e < 0.7–0.2 cm⁻³ and ΔDM < 0.7–2 pc cm⁻³ in that phase.
- **A faint residual at a different velocity.** A free Gaussian fit to the on−off spectrum finds +0.089 ± 0.028 R at V_LSR −9.8 ± 2.4 km/s (V_helio ≈ +2), FWHM_obs 23 km/s. It is 3σ and not at the [O III] velocity. The cross-filament profiles in −20…+9 km/s show a +0.26 ± 0.09 R bump about 0.75° west of the spine. Both are possible but unconfirmed. The beam sampling allows only ~1–2 beams across.
- **Absolute LAMOST intensities (derived, ±13 % scale):**

| | Hα (R) | [N II]6583 (R) | [S II]6716 (R) | [S II]6731 (R) | WHAM Hα at the position (R) |
|---|---|---|---|---|---|
| P1 | 3.6 ± 0.5 | 2.9 ± 0.4 | 1.7 ± 0.3 | 1.6 ± 0.3 | 3.46 |
| P2 | 4.5 ± 0.6 | 3.2 ± 0.5 | 2.1 ± 0.4 | 1.8 ± 0.3 | 4.10 |
| P3 | 3.2 ± 0.5 | 2.4 ± 0.4 | 1.3 ± 0.3 | 1.2 ± 0.3 | 3.01 |
| P4 | 3.0 ± 0.5 | 2.4 ± 0.4 | 1.2 ± 0.3 | 1.2 ± 0.3 | 3.61 |

  - These are WIM intensities. The scale assumes that the 6553 Å sky-line normalisation is common to all fibres and that Hα has little sub-degree structure.
  - The implied [S II] of ~2.4–3.9 R is far above the teammate's NSNS [S II] (≤ 0.35 R on-ridge medians). Either the NSNS narrowband [S II] loses the smooth DIG level, or the LAMOST [S II]/Hα of 0.8–0.9 is biased. This should be cross-checked; I did not resolve it.

### 3.2 Interstellar Ca II / Na I absorption

Vetted components for the most relevant sightlines (V_helio km/s; EW mÅ):

| Star | d (pc) | sep / side | Ca II K components | Na I D2 components | Comment |
|---|---|---|---|---|---|
| 8 Cnc (HARPS, SOPHIE) | 63 | 4.7° E | +19.4 (1.9 ± 0.5); +19.1 (1.5 ± 0.5) | telluric residuals | LIC at 63 pc |
| 11 CMi (UVES) | 92 | 1.2° E (S end) | +18.1 (2.3 ± 0.3); **+25.5 ± 0.5 (1.0 ± 0.3)** | tentative | the +25.5 line is 3.8σ, Ca K only (H not covered at S/N) |
| HD 66552 (SOPHIE) | 121 | 3.9° E | +9.0 (3.1 ± 0.7); +23.5 ± 1.0 (4.2 ± 1.0, FWHM 14) | +29.9 (1.7 ± 0.5) | tentative, single line |
| HD 61038 (SOPHIE) | 149 | 1.9° W | +6.2 (4.9 ± 0.8) | +6.4 (4.8 ± 0.6), +18.4 (9.1 ± 0.9), +27.9 (3.0 ± 0.5, not in D1) | D1 confirms +8, +19 |
| HD 68332 (UVES) | 217 | 6.3° E | +13.1 (7.9 ± 1.5) | +14.1 (9.4 ± 0.9), +18.9 (2.8 ± 0.5) | |
| HD 63774 (UVES) | 220 | 6.2° E | +15.2 (7.2 ± 2.9), +19.6 (4.0 ± 3.3) | +15.4 (18.7 ± 2.0), +19.3 (4.2 ± 1.7) | |
| HD 63432 (UVES) | 250 | 2.7° E (N end) | stellar residual | −0.1 (4.8 ± 0.8), +9.9 (14 ± 1.5) | |
| HD 52124 (UVES) | 264 | 10.9° W | +8.7 (24.7 ± 2.2), +18.5 (7.9 ± 3.0) | +6.9 (12.4 ± 1.2), +14.1 (18.8 ± 1.3) | control |
| **HD 61717 (UVES)** | **272** | **0.97° W (Dec 22.3)** | **−17.0 ± 0.1 (23.7 ± 0.3)**, +9.7 (12.7 ± 2.6), +15.0 (12.7 ± 2.9) | **−17 (2.7 ± 0.3, forced)**, +4.5 (3.0 ± 0.7), +13.1 (21.2 ± 1.1) | high-velocity Ca-rich component |
| HD 64705 (SOPHIE) | 280 | 2.6° E | +9.6 (20.6 ± 1.1); +27.5 stellar (sharp-lined B8, RV* +31.5) | +9.8 (17.2 ± 0.7); +26.4 stellar? | the +26 lines are probably photospheric |
| HD 55869 (SOPHIE) | 423 | 4.8° W | +1.4 (25 ± 1), +11.2 (35 ± 2), +19.2 (8.1 ± 1.1) | +2.2 (10.4), +11.2 (67, saturated), +19.1 (31.6) | |

- **Three velocity groups.** The main local gas sits at +6…+15 km/s (V_LSR −6…+3); WHAM Hα and the LAMOST [N II]/[S II] have the same velocity. The LIC/Aur group sits at +18…+20 km/s and is already present at 63–92 pc. The third group is the HD 61717 component.
- **Nothing coherent at +25 km/s.** Forced fits at +25 km/s (b = 5 km/s) toward HD 61717, the best on-filament background star, give EW(K) = 1.4 ± 0.3 and EW(D2) = 2.2 ± 0.4 mÅ. Both are dominated by the wings of the +13/+15 components, so log N(Ca II, +25) ≲ 10.3.
  - The tentative +23…+28 lines appear in single lines toward 92–149 pc stars, at 1–4σ.
  - This could mean one of three things:
    - The [O III] gas has no associated neutral/low-ion column above ~10¹⁰ cm⁻² Ca II.
    - It is fully ionized (expected for T > 2×10⁴ K, §3.1).
    - It lies beyond 272 pc.
  - Absorption alone cannot bracket the [O III]-gas distance.
- **The HD 61717 component.** Measured values:
  - V_helio = −17.05 km/s (V_LSR −27.9). The line is narrow and barely resolved: b ≈ 1.3 km/s, which implies T ≲ 10⁴ K.
  - log N(Ca II)_AOD = 11.47 (K/H = 1.74, mildly saturated) and log N(Na I) = 10.14. N(Ca II)/N(Na I) = 22 ± 2 (lower limit).
  - Compare the ambient +13 km/s gas on the same sightline: ratio ≈ 3.
  - Interpretation (hypothesis): ratios > 10 at 20–30 km/s from the local rest frame are the classic Routly–Spitzer signature of grain-processed gas behind a shock.
  - Geometry: HD 61717 lies 0.97° on the *concave* (west) side of the arc. The geometry teammate's curvature centre is at RA 94.5°, radius 21°. A shell of radius R seen 1° inside its limb projects only sqrt(1 − (20/21)²) ≈ 0.3 of its expansion velocity. A near-side shell moving at ~100 km/s and lying in front of 272 pc would produce about −30 km/s relative to ambient.
  - Caveats: this is one sightline. HD 64705 (2.6° E, 280 pc) and HD 61038 (1.9° W, 149 pc) do not show it. It could also be an unrelated intermediate-velocity cloudlet.
- **Literature totals along the north part of the spine.**
  - HD 61219 (143 pc, 0.83° W, Dec 24.6) and HD 65856 (133 pc, 4.2° E, Dec 24.3) have unusually strong Ca II: EW(K) = 38.9 and 26.4 mÅ, log N(Ca II) = 11.6 and 12.05.
  - Their Na I is weak: 6.2 mÅ and < 2.5 mÅ D2, so N(Ca II)/N(Na I) ≈ 6 and > 37.
  - HD 61219, HD 65856 and HD 61717 are the three north-half stars at 130–280 pc with usable Ca II and Na I, and all three show Ca-rich gas.
  - For comparison, the +9.6 km/s component toward HD 64705 has a normal N(Ca II)/N(Na I) ≈ 2.7.
  - Component velocities for HD 61219 and HD 65856 are not available. A high-resolution spectrum of HD 61219 would test whether the −17 km/s gas lies within 143 pc.
- **Gaia DR3 DIB 862 nm:** not obtained (§6).

### 3.3 UV
- **Measured toward WD J074735.98+210635.83** (G = 15.79, BP−RP = −0.26, ϖ = 10.05 ± 0.05 mas, d = 99.5 ± 0.5 pc):
  - Si II 1260, C II 1334 and S II 1253 are detected at local velocities (τ-weighted centroids +9.6 and +16.9 km/s for C II and Si II).
  - C II* 1335.7 is detected at −20…+20 km/s: log N = 12.99 ± 0.12.
  - At +15…+40 km/s, the 3σ limits are log N(C II*) < 12.8, log N(Si IV) < 12.2 and log N(N V) < 13.0.
  - O I 1302 is swamped by geocoronal airglow.
  - The local columns are normal; no hot or interface gas lies within 100 pc on a sightline 0.16° from the [O III] spine.
- **Archive inventory** (no measurements made; `spectra_ism_mast_uv_targets.csv`):
  - HD 60848: O8V:pe, Gaia d ≈ 1.9 kpc, 2.3° W of the spine at Dec 17.1. 122 IUE SWP/LWR spectra; this is the best UV background sightline (follow-up).
  - IUE spectra of HD 61913 (0.95°), 11 CMi/HD 62832 (1.2°) and RE 0751+14 = PQ Gem (1.2°).
  - FUSE spectra of σ Gem (38 pc) and Pollux (10 pc), both LISM only.
  - Low-resolution COS G140L spectra of QSOs at 0.5–0.6° (J074815.44+220059.5, 2MASS J07452815+191952). These are not usable for velocity work.
  - COS G160M spectra of green-pea galaxies GP0747+2336 and GP0751+1638.
- **Not done:** O VI surveys (Savage & Lehner 2006; Barstow-type white-dwarf samples) were not checked for stars in this field.

### 3.4 Faraday rotation and pulsars
- **Foreground and on−off:** The 20° cone around (116°, +20°) has median RM = +21.6 rad m⁻² and a robust residual scatter of 24 rad m⁻² about the quadratic foreground. On-filament (|d| < 0.5°, n = 21) minus off (1–4°, n = 129) median residual = **+19.2 ± 8.1 rad m⁻²** (2.4σ).
- **North/south asymmetry:**
  - North (spine Dec +20…+28): E(2–6°) − W(2–6°) = +36.8 ± 4.6 rad m⁻². The step lies at +2.0° E (16–84 % range −0.6…+2.0°).
  - The west side, which is the concave side and also holds the Hα patch, HD 61717 and the Ca-rich stars, is dominated by negative RMs: 70 % negative versus 16 % on the east.
  - South (Dec +12…+20): no step (+5.3 ± 7.3).
  - Interpretation (hypothesis): a magnetised shell whose interior reverses B∥ in the north would produce this. A large-scale foreground RM gradient can do the same, so I do not claim the step as a filament signature.
- **Pulsars (B/psr):**
  - PSR B0656+14: DM = 13.98 pc cm⁻³, RM = +23.5 ± 0.4 rad m⁻², d = 288 pc, giving <n_e> = 0.049 cm⁻³ and <B∥> = 2.1 μG.
  - B0823+26: DM 19.48, RM +5.38, d = 357 pc, giving <n_e> = 0.055 and <B∥> = 0.34 μG.
  - J0751+1807: the only pulsar projected near the filament (2.5° from the centre, 0.95° E of the spine at Dec 18.1). DM = 30.25, timing parallax 1.6 ± 0.8 mas (d ≈ 0.4–1.3 kpc), no RM.
  - Its DM is consistent with a smooth <n_e> of 0.03–0.05 cm⁻³. The filament's ≤ 1–2 pc cm⁻³ (§3.1) is not separable.

## 4. Synthesis: what the spectra say about the SNR / DIG / bubble verdict
1. **The red lines are WIM/DIG.** WHAM confirms the paper's DIG statement independently: velocity, width and intensity are the same on and off the filament. No Hα, and hence no 10⁴ K recombining gas, is associated with the [O III] at a level of 10 % of [O III] (2σ). The paper's [S II]/Hα > 0.4 "shock" criterion and its BPT placement are properties of this background gas. They are not evidence for a shock in the Silk filament.
2. **The [O III] gas must be hot.** The required Hα/[O III] < 0.1 implies collisional excitation at T ≳ 2–3×10⁴ K without a cooled recombination zone: a non-equilibrium (incomplete) shock or an evaporative/mixing interface. This is consistent with a faint, fast shock in a low-density medium, for example the edge of an old, large SNR or superbubble such as the Monogem Ring/Gemini–Monoceros complex. It is inconsistent with photoionized DIG. The [O III] width (FWHM 43–49 km/s versus 30 km/s for the WIM) points the same way.
3. **The kinematic distance is fragile.** The [O III] V_helio of +25 km/s is 13 km/s from the local WIM/neutral gas (+12). It lies between the LIC (+19) and Gem (+34) clouds. At l ≈ 204°, b ≈ +19° differential rotation gives only ≈ 10 km/s per kpc (A sin 2l cos²b). A 13 km/s offset would therefore formally imply ~1.3–1.6 kpc. That is the size of peculiar motion an expanding shell produces, so neither the 160 pc nor the ~1.5 kpc kinematic distance is physically meaningful. The paper-review teammate gives a full treatment.
4. **Distance constraints from this work (measured):**
   - No neutral/low-ion counterpart at the [O III] velocity toward stars at 63–272 pc within 1–5° of the spine. Taken alone this does not constrain the distance, given point 2.
   - No hot ions (Si IV, N V) within 100 pc at 0.16° from the spine (WD J0747+2106).
   - A shock-like, Ca-rich, −30 km/s (relative) component lies in front of 272 pc, 1° inside the arc. If it is associated, the shell is nearer than 272 pc.
   - The data are mildly suggestive of a d ≈ 100–270 pc shell, consistent with the Monogem Ring/PSR B0656+14 at 288 pc only if the arc lies on the near side. This is a hypothesis.
5. **New hypothesis to test.** The Silk filament is the [O III]-bright, still-hot (non-radiative-to-radiative transition) edge of a nearby shell whose interior lies to the west. Tests:
   - Ca II/Na I in more stars at 100–300 pc just west of the spine (expect high-velocity, Ca-rich components).
   - Absorption velocities toward HD 61219 (143 pc).
   - Deep Hα/Hβ/[O III]λ4363 spectroscopy on the ridge (expect high T_e).
   - O VI / C IV toward HD 60848 (IUE archive, 1.9 kpc).
   - eROSITA/ROSAT soft X-rays interior to the arc (imaging teammate).

## 5. Figures
- `figures/spectra_ism_wham_maps.png`: WHAM-SS DR1 Hα around the filament (total Hα, main-component V_LSR, −60…−20 km/s channel, beam-averaged NSNS [O III]), with the spine, envelope and P1–P4.
- `figures/spectra_ism_wham_spectra.png`: WHAM spectra nearest P1–P4 with LAMOST line velocities, the on−off difference spectrum with the [O III]-velocity expectation, and cross-filament Hα profiles by velocity window.
- `figures/spectra_ism_rm_psr.png`: NVSS rotation measures and pulsars around the filament, a smoothed RM map, and north/south cross-filament RM profiles.
- `figures/spectra_ism_abs_profiles.png`: Ca II K and telluric-corrected Na I D2 profiles of 12 vetted early-type stars sorted by distance, with [O III], red-line, LIC and Gem velocities.
- `figures/spectra_ism_abs_veldist.png`: vetted absorption components (velocity versus distance), and literature plus this-work N(Na I) and N(Ca II) versus distance, coded by angle from the spine.
- `figures/spectra_ism_cos_wd.png`: COS G130M AOD profiles (Si II, C II, C II*, Si III, Si IV, N V) toward WD J0747+2106 at 99.5 pc, 0.16° from the spine.

## 6. Limitations
- **WHAM beam.** The 1° beam is about 3.6× the strand FWHM, so filament-confined Hα is diluted. The limits are beam averages.
- **[O III] calibration.** The Hα/[O III] limit uses the teammate's approximate NSNS [O III] calibration. The paper's 10–13 R would make the limit tighter.
- **LAMOST absolute scale.** The scale assumes a common sky-line normalisation and a smooth Hα field.
- **Absorption-fit method.** The running-median continuum can bias EWs of components broader than ~25 km/s, and it leaves stellar cores in sharp-lined stars (flagged).
- **Telluric correction.** The empirical correction is imperfect for ELODIE (R = 42,000) and HARPS 8 Cnc; those Na D results are tentative.
- **Column densities.** Columns are AOD or optically-thin values, so saturated components (for example HD 55869 Na D2 +11 km/s) are lower limits.
- **Instrumental FWHM.** ELODIE and SOPHIE resolution was assumed, not measured (42,000 and 75,000).
- **Distances.** These are SIMBAD (mostly Gaia DR3) parallaxes, without binary or extinction checks.
- **COS spectrum.** The S/N is about 9 per pixel. The WD's own photospheric lines were not modelled.
- **RM step.** The RM analysis cannot separate a filament-related step from a large-scale foreground structure.

## 7. Data gaps
- **Gaia DR3 DIB 862 nm:** not obtained.
  - ESA Gaia TAP (gea.esac.esa.int) failed with an SSL error through the proxy.
  - NOIRLab Data Lab returned 502, and a proxy error for the box query, on gaia_dr3.astrophysical_parameters.
  - VizieR I/355/paramp ("light") has no DIB columns, and J/A+A/674/A40 is not on VizieR.
- **Redfield & Linsky LISM component database:** not on VizieR; only the kinematic calculator was used.
- **Component velocities for the literature catalogues:** Welsh+2010 and Lallement+2003 list totals only. The original references were not retrieved.
- **LAMOST MRS-N per-fibre nebular spectra:** dr10.lamost.org is unreachable, so the geocoronal/sky-subtraction hypothesis is untested.
- **WHAM [S II]/[N II]/[O III] pointed or survey data** for this field: not found.
- **Other archives:** Keck KOA, PolarBase and FEROS-only targets were not queried. ESO FEROS appears in the ObsCore inventory, but no FEROS product was downloaded.
- **UV measurements not made:** IUE high-dispersion spectra of HD 60848, HD 61913 and 11 CMi were not measured. O VI literature for white dwarfs in the field was not checked.
- **Pulsar RM:** J0751+1807 has none.
