# imaging_optical: optical and narrowband precovery search for the Silk Nebula

Topic tag `imaging_optical`. Work date 2026-11-16. Every URL below was accessed on 2026-11-16 unless stated otherwise.

The search covered:
- Hα and narrowband surveys: Finkbeiner 2003, WHAM-SS, SHASSA, VTSS and NSNS.
- Photographic plates: DSS1 and DSS2.
- CCD broadband imaging: Pan-STARRS1 DR2 stacks and SDSS DR18.
- Public amateur images on AstroBin dated before Nov 2023.

The question in each survey is whether the [O III] filament, or the northern Hα/[S II] patch, was already recorded before the SWONS detection in Nov 2023.

NSNS DR0.2 is the reference template throughout. It is a calibrated, linear map in Rayleighs, but its data date from Dec 2022 to Apr 2025, so it is discovery-era data and not precovery.

## 1. Verdicts

The verdict table has three columns: survey, the epoch of its data, and the verdict with its numbers.

The template amplitude **A** is the broadband response measured against the NSNS DR0.2 [O III] map, divided by the response expected if the NSNS map were exactly calibrated. So A = 1 means "as bright as NSNS says", and A ≈ 2.7–6.2 corresponds to the paper's 10–13 R (§5).

| Survey | Epoch | Verdict |
|---|---|---|
| **SDSS DR18 drift-scan sky levels (g)** | 2002-12 to 2007-02 | **Detected (statistical precovery).** A = **1.35 ± 0.19** from 18 runs crossing the ridge at Dec 12–20.5; 17 of 18 runs are positive. Null tests are consistent with zero: r vs [O III] 0.08 ± 0.36; g vs Hα −0.27 ± 0.25. Shifting the template off the ridge gives a null scatter of 0.13 ± 0.25, so A = 1.17 at zero shift is a ≈ 4σ detection against that empirical null. |
| Pan-STARRS1 DR2 stacks (g, r) | 2010–2014 | **Not detected, and insensitive by construction.** Strips: g A = 0.35 ± 0.29 (P3), −0.39 ± 0.33 (P1), −0.25 ± 0.51 (P2), −0.05 ± 1.22 (Dec 22.3); weighted mean −0.02 ± 0.20. The stack sky subtraction removes structure larger than a few arcmin: the Medusa PN (Abell 21, 10′, roughly 140 R [O III] in NSNS) shows **−0.004…−0.027 nmgy/″² observed vs +0.25 predicted** in g. |
| DSS2 J (POSS-II IIIaJ+GG385) | 1987–94 | **Not detected.** A = 4.3 ± 3.2, so A < 14 at 3σ (Dec 12–20). That is roughly < 40 R at the ridge peak, or < 4.7 % of the J sky. Plate background structure is 0.4–1.2 % of sky per 2′ bin. |
| DSS1 O (POSS-I 103aO) | 1950–55 | **Not detected and not useful.** [O III] 5007 lies outside the band and only 4959 is inside. A = 10.6 ± 6.3, corresponding to 0.025 % of sky per R. |
| DSS2 F and DSS1 E (red plates) | 1950–96 | Used as controls, since [O III] is out of band. No [O III]-shaped signal: F −7.3 ± 11.8; E 2.5 ± 2.9, with χ²/dof = 6.4 from plate structure. |
| SHASSA Hα (continuum-subtracted) | 1997–2000 | **Not detected** at Dec 12.5–15 (the southern [O III]-only part). On minus off: −0.14 ± 0.34 R (Dec 12.5–14), −0.05 ± 0.28 R (Dec 14–15), +0.06 ± 0.25 R (P3). The 3σ upper limit is < 0.8–1.0 R. |
| WHAM-SS DR1 Hα (1° beams) | 1997–98 | No excess along the south. In V_LSR −10…+30: −0.02 ± 0.36 R at Dec 15–17. There is a marginal excess where NSNS shows Hα: +0.59 ± 0.34 R (Dec 17–19) and +0.87 ± 0.68 R (Dec 19–21). The single beam on the north Hα patch shows +1.5 ± 0.65 R (V_LSR −10…+30) and +2.3 ± 0.7 R (V_LSR −80…+80). |
| Finkbeiner 2003 Hα composite | WHAM/SHASSA/VTSS | Not detected in the south: −0.14 to +0.31 R. A +0.8 to +0.9 R excess appears at Dec 17–21, coincident with NSNS Hα. There the map is WHAM-limited at 1° resolution and the E−W systematics are 0.3–1.6 R. |
| VTSS Hα | 1995–2000s | **No coverage.** SkyView returns "No data found" for the 17° × 24° field, which lies at b ≈ +19°. |
| AstroBin public images, pre-2023-11 | ≤ 2023-10 | **No narrowband image covers the ridge.** That holds for 507 solved images tagged "OIII" and 154 tagged "narrowband". 169 wide-field (≥ 3°) solved images do cover it, but they are broadband nightscapes or meteor shots. The expected contrast there is 0.3–0.5 % of sky, so they cannot show the filament. Six were inspected: none shows it. |
| NSNS DR0.2 (reference) | 2022-12-12 … 2025-04-04 | Detected, but not precovery. [O III] on-ridge excess (\|Δ\| ≤ 0.25°): P1 2.04, P2 2.06, P3 2.18 R; peaks 2.95 / 3.64 / 3.37 R. |

**Bottom line.**
- The earliest record of the Silk Nebula's [O III] emission I can find is in the **SDSS imaging sky levels of 2002–2007**. The Silk raises the g-band sky by about 0.4 % on the ridge (~27.9 AB mag arcsec⁻² at P1).
- That was invisible in any single frame. It is recoverable only statistically, with the NSNS map as template.
- No photographic, Hα-survey or pre-2023 amateur image shows it.

## 2. Expected broadband and photographic signal

Script: `scripts/imaging_optical_expected.py`. Output: `data/derived/imaging_optical_expected_signal.json`.

Method:
- A line of intensity I (in R) gives a flux of I × 10⁶/4π photons s⁻¹ cm⁻² sr⁻¹.
- It is converted to AB f_ν through the band as f_ν = F · T(λ₀) · λ₀ / (c ∫T dλ/λ).
- Filter curves: PS1 from Tonry et al. 2012, Table 3 (https://ipp.ifa.hawaii.edu/ps1.filters/apj425122t3_mrt.txt); SDSS 2010 g, r from speclite (https://github.com/desihub/speclite).
- POSS photographic bands are approximated by top-hats: J 3850–5400, F 5900–6900, O 3300–5000, E 6200–6700 Å.

The surface brightness of 1 R in each band (in mag arcsec⁻²) is:

| Band | 1 R [O III] 5007 (+0.33 R 4959) | 1 R Hα | Sky assumed |
|---|---|---|---|
| PS1 g / SDSS g | 29.39 / 29.36 | — | 21.9 / 22.0 |
| PS1 r / SDSS r | — | 29.49 / 29.31 | 20.9 / 21.0 |
| POSS-II J / POSS-I O | 29.68 / 31.41 (4959 only) | — | 22.3 / 22.4 |
| POSS-II F / POSS-I E | — | 29.15 / 28.39 | 21.0 / 20.9 |

Checks of the brief's numbers:
- **The "12 R ≈ 27 AB mag arcsec⁻² in g" figure is verified.** 12 R of 5007 alone gives 27.0. Adding 4959 (4 R) gives 26.69 in PS1 g and 26.66 in SDSS g. That is 1.2–1.4 % of a 21.9–22.0 mag sky.
- **10 R of [O III] is 1.12 % of the POSS-II J sky** (0.112 % per R).
- The paper's Table 3 values of 23.60 / 23.85 / 23.87 mag arcsec⁻² convert to 12.7 / 10.1 / 9.9 R, assuming a 6.0 nm top-hat filter (SWONS, paper §II).

## 3. Data, methods and provenance

### 3.1 NSNS DR0.2 (template and reference)

Script: `scripts/imaging_optical_nsns.py`.

Data:
- HiPS: https://simg.de/nebulae3/dr0_2/{oiii,halpha,sii,vc,sn_*}/, order 3 (51.5″ pixels).
- Linear intensity in Rayleighs; CC-BY-NC-SA 4.0; S. Ziegenbalg.
- HiPS `properties`: `t_min` = MJD 59926 (2022-12-12), `t_max` = MJD 60769 (2025-04-04); release date 2025-06-07.
- Fetch script: `imaging_optical_fetch_nsns.py`; logs in `data/raw/nsns_hips/fetch_log_*.json`.

Method:
- Mosaicked to a 1′ TAN grid.
- The ridge was traced row by row on a 3′ × 12′ matched-filtered map. A quartic RA(Dec) fit gives `data/derived/imaging_optical_nsns_ridge.json` (Dec 9.9–27.6).
- My ridge agrees with the geometry teammate's spine v1.0 to within 0.01° at Dec 12–15. It lies 0.18° west of the spine at Dec 15–18 and 0.42° west at Dec 18–21, where it follows the western, Hα-bright strand. The SDSS and DSS fits use the NSNS map itself as template, so they are not sensitive to this choice.

NSNS results (`imaging_optical_nsns_measurements.json`):
- **[O III]** excess on ridge: south_bright 1.60 R, mid (Dec 17–19) 1.21 R, northern arc 0.58–0.62 R.
- **Hα** excess on ridge: P1 +0.23, P2 +0.44, P3 −0.01; mid +1.70 R.
- **[S II]** excess: P1 +0.16, P2 +0.15, P3 0.00; mid +0.69 R.
- **North Hα patch** (geometry v1.0 polygon): Hα +3.76 R, [S II] +2.88 R, [O III] +0.40 R, so [S II]/Hα ≈ 0.77.
- **P4:** [O III] −0.11 R, Hα −0.23 R, [S II] −0.25 R relative to its surroundings. It is a genuine nebula-free region.
- Uncertainties: statistical 0.01–0.09 R. The systematic error from E/W background asymmetry is 0.05–0.15 R for [O III] and 0.15–0.65 R for Hα.

### 3.2 Hα surveys

Script: `scripts/imaging_optical_halpha.py`.

Data:
- **Finkbeiner 2003:** `Halpha_fwhm06_1024.fits`, v1.1 (header dated 2003-07-27), Nside 1024, NESTED, Galactic, in R. Its size (50 342 400 B) matches https://faun.rc.fas.harvard.edu/dfink/skymaps/halpha/data/v1_1/healpix/Halpha_fwhm06_1024.fits. LAMBDA only hosts the Nside-512 file (`lambda_halpha_fwhm06_0512.fits`). The exact download URL was not logged.
- **WHAM-SS DR1:** `wham-ss-DR1-v161116-170912.fits` (VERSION 161116, file dated 2017-09-12; Haffner et al. 2003, ApJS 149, 405). The download URL was not logged. The WHAM site (www.astro.wisc.edu/wham-site) currently answers HTTP 202 or 404 to scripted requests.
- **SHASSA** via SkyView (https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl, surveys "SHASSA H/C/CC/Sm"): 116.5, +10.5, 13° × 11°, 48″ pixels. "Sm" is the smoothed continuum-subtracted map in dR. Exact query URLs are in `data/raw/skyview/SHASSA_*/*.json`.
- **VTSS** ("VTSS Ha", "VTSS Bg"): SkyView returns a 1-pixel "No data found" file.

Method:
- On-ridge = \|across\| ≤ 0.25° (maps) or ≤ 0.5° (WHAM beam centres).
- Off-ridge = 0.8–2° (maps) or 1–3° (WHAM).
- Errors are the off-region scatter and the E−W asymmetry.
- WHAM spectra are averaged over the beams in each region.

WHAM result: the mean on − off spectrum at Dec 13–19.5 peaks at V_LSR ≈ −2 km/s with 0.013 R (km/s)⁻¹, about 0.3 R in total. That velocity is the LAMOST red-line velocity, not the [O III] velocity (+11…+16).

SHASSA result: SHASSA gives no Hα at the 3σ < 0.8–1.0 R level on 4′ scales where NSNS sees 1.2–2.2 R of [O III]. So the southern filament has [O III]/Hα > 1.5–2.8 (3σ, NSNS scale), or > 2–3.8 on the SDSS-corrected [O III] scale.

### 3.3 DSS photographic plates

Scripts: `scripts/imaging_optical_fetch_dss.py`, `imaging_optical_dss.py`, `imaging_optical_dss_fig.py`.

Data:
- SkyView "DSS1 Blue/Red" and "DSS2 Blue/Red": 2° × 2° tiles at 6″ pixels, Clip sampler, 9 rows (Dec 11–27) × 2 columns (ridge ± 1°).
- Query URL and plate IDs per tile are in `data/raw/skyview/DSS*/…json`. Plates: DSS2 J xj428–xj703; F xp428–xp703; DSS1 O xo309–xo483; E xe309–xe483.
- Star masks: **Gaia DR3** G < 17.5 from the NOIRLab Data Lab TAP (https://datalab.noirlab.edu/tap/sync, table `gaia_dr3.gaia_source`; queries in `data/raw/gaia_dr3_datalab/query_log.json`). Mask radius r = 8″ + 110″·10^(−0.2(G−6)), plus 3σ outliers.
- Photometric calibration: PS1 DR2 mean photometry (https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv) of isolated 18.5–20.5 mag stars, fitting the density "volume" against flux.

Processing:
- 1′ block medians of the unmasked pixels.
- For each row, the E and W tiles share **one** zero point. A quadratic background in (across, along) is fitted at 0.7–2.2° from the ridge.
- The contrast (D − D_bg)/D_sky is taken as the fractional sky excess. This assumes a linear density response. The star-based calibration and the literature sky brightness agree to 20 % for J and O, and to a factor of 2.5–3.5 for F and E. Amplitudes are rescaled to the literature sky.
- The NSNS [O III] template is fitted with errors from 17 off-ridge template shifts.
- An earlier per-tile zero point produced a step at the tile seam near the ridge. It is superseded; see `_notes` in `imaging_optical_dss_summary.json`.

Result: the plates are dominated by 0.5–2 % emulsion and background structure on 0.3–1° scales. No [O III]-shaped excess appears.

### 3.4 Pan-STARRS1 DR2 stacks

Scripts: `scripts/imaging_optical_fetch_ps1.py`, `imaging_optical_ps1.py`, `imaging_optical_ps1_fig.py`.

Data:
- Stack images "stk.<b>.unconv" from https://ps1images.stsci.edu/cgi-bin/ps1filenames.py and fitscut.cgi.
- Four cross-filament strips at Dec 14.72 (P3), 16.60 (P1), 17.74 (P2) and 22.3, each 3.4° long. They are built from 26 overlapping 12′ cutouts per band, block-averaged to 2.5″.
- 208 cutouts plus 2 cutouts of the Medusa (Abell 21) test field, 76 MB in total. URLs are in `data/raw/ps1_dr2/cutouts/*.json`.
- Zero point ZP = 25 + 2.5 log EXPTIME. Median exposures are 600–900 s (g) and 1100–1760 s (r).

Method: Gaia masks, 1′ medians, 2′ cross-ridge bins and an NSNS template fit.

Result:
- The 2′ bin scatter is (1.8–4.1) × 10⁻³ nmgy arcsec⁻², i.e. 28.0–29.0 mag arcsec⁻².
- That is below the expected SDSS-scaled ridge peak of about 8 × 10⁻³. **Yet nothing is seen.**
- The Medusa test shows why: the PS1 stack sky model removes, and slightly over-subtracts, emission on scales of a few arcmin and larger.
- So the PS1 DR2 stacks cannot constrain the Silk. Only a re-reduction from warps or raw exposures with a large-scale background could.

### 3.5 SDSS DR18 imaging: drift-scan sky levels

Scripts: `scripts/imaging_optical_fetch_sdss.py` (fields), `imaging_optical_sdss_sky.py`, `imaging_optical_sdss_profile.py`, `imaging_optical_sdss_placebo.py`, `imaging_optical_sdss_placebo_fig.py`.

Data:
- SkyServer DR18 `Field` table: 20 435 fields at RA 109–124, Dec 6–30, 70 runs.
- Columns `sky_g` and `sky_r` give the median sky of each 13.5′ × 9.8′ frame in nmgy arcsec⁻².
- Query: https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch (full URL in `data/raw/sdss_dr18/fields_box.json`).
- No frame files were needed. SAS frames were not downloaded.

Why the sky level works as a photometer:
- Degree-scale emission is part of the "sky".
- Legacy and SEGUE stripes run at PA ≈ 47° here, obliquely across the N–S filament.
- The six camcols observe simultaneously, 0.42° apart perpendicular to the scan, so they cross the ridge at different times.

Model per run:
- sky[f, c] = S[f] + g_c + h_c·f + A·T[f, c].
- S[f] is a free common mode per field, i.e. per time. It absorbs airglow, moon and twilight changes.
- g_c and h_c are per-camcol offsets and drifts.
- T is the NSNS map at the field centre, 11′-median-filtered and smoothed with σ = 4′.

Error and selection choices:
- Errors are max(formal × √2, scatter of RA-shifted placebo templates).
- Fields within 0.25° of 53 Gaia stars with G < 5, or within 0.12° of stars with G < 7, are excluded (307 fields).
- Expected signal for an exactly calibrated NSNS map: 1.80 × 10⁻³ nmgy arcsec⁻² per R in g.

Robustness of A (A/A_exp, 18 core runs):

| Variant | A/A_exp |
|---|---|
| Raw template | 1.27 ± 0.18 |
| Median-filtered | 1.20 ± 0.18 |
| Ridge-only (tapered at 1.2°) | 1.41 ± 0.19 |
| Plus star exclusion (adopted) | **1.35 ± 0.19**, χ²/dof 0.41 |
| 2002–2004 runs (13) | 1.40 ± 0.21 |
| 2004–2007 runs (5) | 1.11 ± 0.43 |
| 4° window, formal errors | 1.17 ± 0.14 |

Shift test (`imaging_optical_sdss_placebo_w4.0.json`, figure `imaging_optical_sdss_shift.png`):
- The amplitude peaks at zero shift (1.17–1.26 at −0.2…0°) with negative side lobes at ±1.2°.
- That reproduces the shape of an injected NSNS-template signal, with a best-fit scale of 1.52 and residual rms 0.26.
- Templates shifted by \|s\| ≥ 2° give 0.13 ± 0.25.

All Dec ≥ 21 crossings (13 runs, mostly 2001–02) are noisy, with sky up to 2.6 nmgy″⁻² and field rms 2–4× higher. They give no constraint: A = −4 ± 3, with χ² dominated by a few runs. The r-band Hα test is not constraining (−0.19 ± 0.65).

Conversion to intensity:
- A ≈ 1.3 ± 0.2 (stat) ± 0.2 (method spread) scales the NSNS [O III] intensities, if all of the [O III]-correlated g excess is [O III] 4959+5007. Hβ and continuum would make this an upper bound.
- Ridge medians: P1 2.8, P2 2.8, P3 2.9 R. Peaks at 1′: 4.0, 4.9, 4.5 R (± ~20 %).
- At P1 that is g ≈ 27.9 mag arcsec⁻², 0.4 % of the sky.

### 3.6 Amateur archives

Scripts: `scripts/imaging_optical_astrobin_search.py`, `imaging_optical_astrobin_filter.py`.

Search:
- AstroBin public search API: https://www.astrobin.com/api/v2/images/image-search/
- Filters: `coord_ra_min/max` 100–133, `coord_dec_min/max` 0–36, `field_radius_min`, `date_published_max=2023-11-01`, and text.
- Raw pages are in `data/raw/astrobin_search/`. Candidate tables are `data/derived/imaging_optical_astrobin_candidates_<tag>.csv`.
- A field "covers" the ridge if NSNS ridge points (Dec 12–27) lie inside both 0.75 × fieldRadius and the solved RA/Dec box.

| Query | Images returned (plate-solved) | Cover the ridge |
|---|---|---|
| "OIII", r ≥ 1° | 617 (507) | **0** |
| "narrowband", r ≥ 2° | 166 (154) | **0** |
| "Halpha", r ≥ 2° | 9 (8) | **0** |
| "HOO", r ≥ 1° | 116 (107) | **0** |
| "Monogem", r ≥ 1° | 1 | **0** |
| any, r ≥ 3° | 2683 (1569) | 169: all broadband nightscapes, meteor or comet shots, constellation portraits |

Images inspected (HD thumbnails in `data/raw/astrobin_thumbs/`, log `download_log.json`):
- `7e7t1w`, B. Darvish, 2022-11, 12-panel L-eXtreme mosaic: ends at about RA 110°.
- `3dxohy`, Dong Liang, 2022-08, 100 h "Greater Orion": the Silk would sit at the extreme edge; nothing visible.
- `4p3q62`, 2022-10, an RGB field centred at 116.9, +21.4 on the filament. It is dominated by lunar halo gradients; no filament.
- `kgykw2`, Minh Lết, 2021-02, L-eXtreme Barnard's Loop field: its east edge is near the Rosette and Seagull (RA ≈ 105°), so the Silk lies outside.
- `tmjev3`, A. Wegdan, 2022-12, broadband ultra-wide Barnard's Loop: the Silk is outside or at the edge; nothing visible.
- `y0z9jq`, astroclausi, 2020-02, broadband DSLR "Gemini widefield" around Castor and Pollux, covering the northern arc: nothing visible.

None show the Silk. This is expected: the ridge is ≈ 0.4 % of the sky in a green channel, below the gradients in such images.

NSNS-era material and other pre-2023 sources:
- The NSNS "Milky Way from Monoceros to Gemini" mosaic (https://simg.de/nebulae3/mon-gem.html; JPGs downloaded to `data/raw/nsns_mongem/`) is 21° × 40°, centred 6h46 +5°. It has Hα and blue-continuum channels; the blue channel includes [O III] and Hβ. Its east edge is at about RA 111–112°, so **it does not reach the ridge** (RA ≈ 116.5°).
- Ziegenbalg states the arc appeared in "the preview images of the first clear night after the camera array reached its full configuration" (≈ Jan 2024; `sources/web/discovery_nsns_gem01.md`). DR0.2 does not give per-pixel epochs, so whether earlier 2022–23 NSNS [O III] frames covered the field cannot be established.
- No public pre-2023 Drechsler/Strottner/Schaeffer image of the field was found. That history is covered by the discovery teammate.

## 4. What the precovery results mean

1. **Intensity calibration.**
   - The SDSS sky-level photometry is an independent broadband calibration, 15–20 years older than the discovery.
   - It agrees with NSNS DR0.2 to a factor of 1.3 ± 0.3 on scales of 10′ and more.
   - It is **2.6–4.5× below the paper's 10–13 R** at P1–P3 (~24.9–25.3 instead of 23.6–23.9 mag arcsec⁻² in a 6 nm [O III] band). A paper-level amplitude (≥ 2.7) is excluded at roughly 4–7σ, depending on the error model.
   - The paper notes that its P1 photometry is contaminated by starlight.
   - The 10–13 R could still apply to arcminute-scale bright strands if the filling factor on 10′ scales is ≈ 1/3. Neither the SDSS nor the NSNS (51″) maps can test that.
   - Consequences for the paper's density, emission measure and energy estimates belong to the paper-review teammate.
2. **Temporal stability.** The [O III] filament was already present in 2002–2007. The two SDSS epoch subsets agree: A = 1.40 ± 0.21 (2002–04) and 1.11 ± 0.43 (2004–07). The SDSS-era to NSNS-era (2022–25) ratio, 1.3 ± 0.3, is degenerate with the NSNS absolute calibration, so it is not a variability measurement. All that can be said is that no order-of-magnitude change occurred over about 20 years, as expected for ages ≥ 10⁴ yr.
3. **[O III]/Hα.**
   - Southern filament: SHASSA (1997–2000) shows no Hα at < 0.8–1.0 R (3σ, 4′).
   - WHAM at Dec 15–17 gives −0.02 ± 0.36 R over V_LSR −10…+30 and 0.00 ± 0.10 R over +5…+25 km/s, the window around the [O III] velocity (1° beams).
   - Together these give [O III]/Hα ≳ 2–4.
   - Where Hα exists (Dec 17–21 and the north patch), it sits at V_LSR ≈ −2 km/s, the red-line velocity. That supports the paper's finding that the red lines differ kinematically from the [O III].
   - The north patch is a real Hα+[S II] structure, also seen in WHAM (+1.5 to +2.3 R in one beam) and in Finkbeiner (+0.9 R).
4. **Implications for other surveys.**
   - Sky-subtracted stacks (PS1 DR2, and by extension Legacy Surveys and DES) are blind to degree-scale emission at this level.
   - Survey sky estimates, or calibration-frame sky levels, can recover it.
   - The same `Field.sky_g` method could map the northern arc with deeper or more runs, or be applied to the DESI/BASS/MzLS sky levels if those become accessible.

## 5. Paper-scaling factor used in the figures

- Paper: 10.0–12.7 R (P1–P3).
- NSNS: on-ridge medians 2.04–2.18 R; 1′ peaks 2.95–3.64 R.
- Ratio: 2.7–6.2. The lower end, 2.7, is drawn as "paper, min." in the figures.

## 6. Figures

| File | Caption |
|---|---|
| `figures/imaging_optical_halpha_wide.png` | Finkbeiner Hα, WHAM-SS (−10…+30 km/s), NSNS Hα and NSNS [O III] over RA 108–125, Dec 5–29, with spine v1.0, the W strand, the Hα patch and P1–P4. The [O III] arc has no Hα counterpart except at Dec 17–21. |
| `figures/imaging_optical_shassa_south.png` | SHASSA continuum-subtracted Hα (4′) at the south end with NSNS [O III] contours at 1 and 2 R; cross-ridge Hα profiles (Dec 12.5–14, 14–15) show no Hα at the [O III] ridge. |
| `figures/imaging_optical_wham_spectra.png` | Mean WHAM-SS Hα spectra on and off the ridge (Dec 13–19.5) and of the north-patch beam; the on − off residual peaks at V_LSR ≈ −2 km/s. |
| `figures/imaging_optical_dss_strips.png` | NSNS [O III] beside processed DSS2 J, DSS2 F, DSS1 O and DSS1 E contrast mosaics (star-masked, plate-flattened, ±1.5 % of sky); the mean Dec 12–20 profile compared with SDSS-scaled and paper-scaled [O III]; 3σ limits. |
| `figures/imaging_optical_ps1_profiles.png` | PS1 g/r cross-ridge profiles at P3, P1, P2 and Dec 22.3 compared with SDSS-scaled [O III]; the Medusa PN test shows the PS1 sky subtraction removes extended emission. |
| `figures/imaging_optical_sdss_sky.png` | Fields of the 18 SDSS runs used over NSNS [O III]; the SDSS g sky residual vs distance from the ridge with the fitted (×1.35), ×1 and paper-scaled templates; per-run amplitudes vs crossing Dec, including null tests. |
| `figures/imaging_optical_sdss_shift.png` | SDSS amplitude vs RA shift of the template (observed vs injected NSNS × 1 and × 1.52, paper × 2.7), with the ±1σ null band from \|shift\| ≥ 2°. |

## 7. Limitations and data gaps

**Method limitations**
- **SDSS** is a template-correlation detection, not an image. It assumes the NSNS morphology. A, the amplitude ratio, carries the NSNS-template smoothing (≈ 8 % systematic) and the NSNS large-scale (> 3°) filtering.
  - Its errors are not fully Gaussian: per-run χ²/dof is 0.4 with placebo errors and 3.7 with formal errors. The empirical shift-null scatter (0.25) is the most conservative error.
  - Only Dec 12–20.5 is constrained. The northern arc (Dec 21–28) is not.
  - The `Field.sky_g` estimator is assumed to respond linearly to diffuse light. The template correlation and null tests support this, but it was not separately validated, for example by an injection into real frames.
- **DSS** contrast-to-intensity conversion assumes a linear density response (a factor-of-~2 systematic for blue plates and about 3 for red). Plate-to-plate structure dominates. Only a 4° wide band around the ridge was fetched.
- **PS1:** only the stacks were tested. Warps or raw exposures with large-scale backgrounds were not attempted; that would take a large download.
- The **cross-ridge coordinate** uses my NSNS ridge, which is offset 0.2–0.4° west of spine v1.0 at Dec 15–21. That is inside all windows; the templates are the maps themselves.

**Data not obtained**
- **Legacy Surveys** (legacysurvey.org, NERSC): unreachable; not tested. LS DR10 coadds are sky-subtracted in any case.
- **SDSS frame files:** not downloaded. The test used the Field-table sky values.
- **VTSS:** no coverage.
- **IPHAS/VPHAS+:** not applicable (\|b\| < 5–10°).
- **WHAM [O III]:** no published WHAM [O III] pointing list for this field could be checked (no catalogue access). Madsen et al. 2006 (ApJ 652, 401) observed [O III] only along selected sightlines. Gap.
- **AstroBin:**
  - Only metadata search and six HD thumbnails were used; full-resolution files were not downloaded.
  - Images without a plate solution, or with other wording, can be missed.
  - Pre-2023 NSNS frames and private amateur data are not public.
- **Provenance:** exact download URLs for the Finkbeiner and WHAM files were not logged; the file identities are given in §3.2.

**Downloads:** about 0.9 GB in total, all under `data/raw/{skyview,lambda_finkbeiner_halpha,wham_ss_dr1,nsns_hips,gaia_dr3_datalab,ps1_dr2,sdss_dr18,filter_curves,astrobin_search,astrobin_thumbs,nsns_mongem}`.

## 8. Derived files

| File | Contents |
|---|---|
| `data/derived/imaging_optical_verdicts.json` | Per-survey verdicts and numbers |
| `imaging_optical_sdss_sky.json` | Per-run fits for all variants and null tests |
| `imaging_optical_sdss_profile.json` | Model-light residual profile |
| `imaging_optical_sdss_placebo_w4.0.json` | Shift and injection curve |
| `imaging_optical_dss_summary.json` | Row-level DSS amplitudes; supersedes the amplitudes in `imaging_optical_dss_profiles.json` |
| `imaging_optical_dss_binned_<band>.npz` | Binned DSS tiles |
| `imaging_optical_ps1_profiles.json`, `imaging_optical_ps1_summary.json` | PS1 strip profiles and amplitudes |
| `imaging_optical_halpha_measurements.json`, `imaging_optical_wham_beams.csv` | Hα survey measurements and WHAM beams |
| `imaging_optical_nsns_measurements.json`, `imaging_optical_nsns_ridge.json`, `imaging_optical_nsns_maps.npz` | NSNS reference measurements, ridge and maps |
| `imaging_optical_expected_signal.json` | Expected signal per band |
| `imaging_optical_astrobin_candidates_*.csv` | AstroBin candidate tables |
