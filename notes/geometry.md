# Silk Nebula: geometry and filament trace (topic tag: `geometry`)

Date: 2026-11-16. All coordinates are ICRS/J2000 in degrees unless stated. Galactic (l, b) come from astropy.
Main deliverable: `data/derived/silk_trace.json` (`"version": "1.0"`, `"preliminary": false`). A rough preliminary
version (`0.1-prelim`, built from Table 1 and the 1°×10° extent only) was written at about 16:08 and then overwritten.

## 1. What other teammates should use

| Need | Use |
|---|---|
| Spine of the [O III] arc (ordered S→N) | `silk_trace.json → spine` (66 points, Dec 12.00–28.25, 0.25° step; `ra, dec, sigma_deg, source, hwhm_deg, l, b, flag`) |
| Second bright strand (P2 lies on it) | `strands.W.points` (Dec 15.0–19.75) |
| Faint far-west strand | `strands.W2.points` (Dec 21.5–24.0) |
| "On-filament" selection polygon | `envelope_polygon` (E+W+W2 ridges ± their HWHM) |
| Northern Hα/[S II] patch | `halpha_patch` (half-maximum polygon), `halpha_patch_info.outer_polygon_30pct` (larger, 30% of peak) |
| LAMOST plate | `lamost_field` = HIP 37797 (116.25248, +15.89500), r = 2.5° |
| P1–P4 | `regions` (Table 1 values, plus the outlines as actually drawn in Fig 1b, registered) |
| Digitized LAMOST [O III] fibres | `data/derived/geometry_lamost_oiii_points_digitized.csv` |
| Flat table of all strand points | `data/derived/geometry_spine_points.csv` |

`scripts/silk_geometry.py::load_trace()` returns the JSON as a dict.

## 2. Methods

### 2.1 Paper figures (Wu et al. 2026, arXiv:2610.13329v1, `sources/arxiv/2610.13329v1.pdf`)
* **Vector markers, not pixels.** Figs 2, 4 and 5 are matplotlib vector graphics. I used PyMuPDF `get_drawings()`
  (`scripts/geometry_extract_pdf_vectors.py`) to pull out the exact marker centres and fill colours:
  2684 red fibre circles in Fig 2 (the paper's "2,684 usable spectra"), 173 [O III] triangles in Fig 4, and
  2335 V_LSR points in Fig 5 (2174 in one scatter layer plus 161 in a second layer drawn last).
* **Axes.** I fitted the tick-mark positions (also vector) linearly. The residuals are ≤0.0001°.
  Fig 5 has decimal-degree labels: RA runs 119.516→113.000 and Dec 11.313→20.960. The RA labels of Figs 2 and 4 are
  truncated "7hMMm" strings. I therefore registered their RA axes by matching their markers to the Fig 5 points
  (Fig 2) and to the Fig 2 fibres (Fig 4) (`scripts/geometry_register_figs.py`):
  * Fig 2 has the same limits as Fig 5. 2335 of its fibres match Fig 5 points within 2″, with a median separation of 0.26″.
  * Fig 4: RA 118.800→114.600 and Dec 12.80→20.00. The tick labels resolve to 7h55.20m, 7h51.00m, 7h46.80m,
    7h42.60m and 7h38.40m. All 173 triangles fall on Fig 2 fibres, with rms 0.11″ in RA and 0.04″ in Dec.
* **Colours to values.** Fills are exact `jet` LUT colours (maximum mismatch 0.003), so I inverted them to values:
  * Fig 4 relative [O III] intensity: colourbar 8.46–77.33, quantisation ±0.13.
  * Fig 5 Hα V_LSR: colourbar 0–15 km/s, quantisation ±0.03 km/s. 149 points sit at the bottom end (≤0) and 31 at the
    top end (≥15). These are flagged as censored.
* **Fig 1 (SWONS deep Hα+[O III] composite).** The embedded high-resolution image (img-001, 1508×5274 px) has no axes.
  `scripts/geometry_register_fig1.py` works as follows:
  * Point sources are detected with a high-pass filter.
  * A coarse scale/offset vote is matched against Tycho-2 (VT<10).
  * An order-3 polynomial is then fitted on the gnomonic plane.
  * Result: **227 Tycho-2 stars, median residual 2.8″, rms 5.2″, 90th percentile 7.1″**, at 9.25″/px with normal parity.

  The overlays drawn on Fig 1b (white circle, P1–P4 outlines) were mapped from page coordinates into this WCS.
* **Fig 4 background ([O III] filter image).** I registered it independently to Tycho-2: 144 stars, order 2,
  rms 19″, median 10″.

### 2.2 NSNS images (S. Ziegenbalg, Northern Sky Narrowband Survey)
* Source page: https://simg.de/nebulae3/gem-01.html. The viewer data files `gem-01-ohs.js` and `gem-01c-ohs.js` give
  explicit WCS parameters. They were saved as text in `data/raw/nsns_simg/*.js.txt` and **not executed**.
  * Narrow view: STG projection, centre (116.0, +18.5), 10″/px at full resolution.
  * Wide view: STG, centre (83.0, +13.5), 10″/px.
* Images used:
  * `gem-01-ohs-2.jpg`: 1575×2575 px at 20″/px, channels R=[O III], G=Hα, B=[S II].
  * Six level-4 tiles of the wide view, `gem-01c-ohs/4-{0,1,2}-{0,1}.jpg`: 40″/px, continuum-free apart from the
    brightest stars.
* I checked the WCS against Tycho-2 (`scripts/geometry_nsns.py`):
  * Narrow: 330 stars, median residual 3.6″, offset (−0.005, +0.05) px.
  * Wide: 152 stars, median 6.5″, offset 0.33 px (applied).

### 2.3 Ridge (spine) measurement
* All four registered images were resampled onto one 0.02° RA/Dec grid (`scripts/geometry_build_maps.py`,
  `data/derived/geometry_maps.npz`). Stars were suppressed with median and grey-opening filters.
  The tracers used were:
  * SWONS colour index B−R;
  * the inverted Fig 4 grey image;
  * NSNS narrow [O III];
  * NSNS wide [O III].
* In 0.1° Dec bins, each RA profile was background-subtracted (1.6° running median) and fitted with one Gaussian per
  strand present at that Dec, plus a linear term (`scripts/geometry_trace_fit.py`).
  * Initial positions came from visual priors. North of Dec 22 and south of Dec 13.5, these priors were taken from the
    automatic peak catalogue `geometry_ridge_peaks.csv`.
  * Centres were constrained to ±0.25° of the prior. A second pass used spline-smoothed first-pass tracks as priors.
* Fits were accepted when amplitude S/N >5, the centre was not at a bound, and the error was <0.1°.
* **Combination** (`scripts/geometry_build_trace.py`):
  * In 0.25° Dec bins, take each image's median centre, then average the images.
  * `sigma_deg` = √(source-to-source scatter² + median fit error² + 0.005²), expressed as a great-circle angle.
    The scatter is the full standard deviation, not divided by √n, because the differences between images are systematic.
  * Where only NSNS-wide covers (Dec ≥26.25), the scatter term is set to 0.02°.
  * Ends were trimmed to the contiguous part with median S/N ≥8.
  * W-strand points north of Dec 19.75 were dropped, because the images disagree there by more than 0.1°.
* **Curvature.** Weighted small-circle fits on the sphere (residual = angular distance minus radius), with 400 bootstrap
  resamples. A great-circle fit and a model-free sagitta are given for comparison.
* **Hα patch.** I took the half-maximum region (bg + 0.5·(peak − bg)) of the NSNS narrow Hα channel after 0.06° Gaussian
  smoothing, in the box RA 114.9–116.4, Dec 19.2–22.4, keeping the connected component that contains the peak.
  For comparison I also give a 30%-of-peak contour, the same procedure on [S II], and on the SWONS red channel.

## 3. Sources and access (all 2026-11-16)

| What | URL / location | Notes |
|---|---|---|
| Wu et al. 2026, arXiv:2610.13329v1 (AJ doi:10.3847/1538-3881/aeb11c) | `sources/arxiv/2610.13329v1.pdf` | embedded images extracted to `data/raw/paper_figures/` |
| HIP 37797 = HD 62523 | https://simbad.cfa.harvard.edu/simbad/sim-id?Ident=HIP+37797&output.format=ASCII → `data/raw/simbad/hip37797_simbad_cfa_mirror.txt` | 07 45 00.596 +15 53 42.01 (Gaia DR3 3165676224716087168; V=7.87, A0). simbad.cds.unistra.fr returned a CDS "page unavailable" 404 page |
| Tycho-2 (Høg et al. 2000), VizieR I/259/tyc2 | https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=I/259/tyc2&-c=116.5+16.5&-c.bd=12.5x25&-out.max=200000&-out=RAmdeg,DEmdeg,BTmag,VTmag,TYC1,TYC2,TYC3&VTmag=<11 → `data/raw/vizier_tycho2/` | 8,140 stars. The TAP endpoint failed on both mirrors (an HTML error at unistra, "no table found" at CfA); ASU worked |
| NSNS O III arc page, viewer WCS and images | https://simg.de/nebulae3/gem-01.html, …/gem-01-ohs.js, …/gem-01-ohs-2.jpg, …/gem-01c-ohs.js, …/gem-01c-ohs/4-{d}-{p}.jpg, https://simg.de/viewer.js (read as text only to learn the tile naming) → `data/raw/nsns_simg/` | NSNS annotations O200.0, O200.0a and O200.0b were copied into `silk_trace.json → literature_geometry` |
| Fesen et al. 2026, arXiv:2608.26666v2 | teammate copy `data/raw/discovery_arxiv/fesen2026_2608.26666v2.pdf` (text read only) | §6.2.1 describes a "nearly vertical 15 degree long group of [O III] filaments … along the Monogem's eastern and northern limbs" (the Silk arc) and Hα+[S II] filaments "at the northern end of the remnant's [O III] filaments" (their Fig 38, NSNS 7.8°×12° images). The eROSITA shell they adopt is 24° in diameter, centred at (107.08, +13.5) (Knies et al. 2024) |

## 4. Results

### 4.1 The paper's figures (measured)
* **The SWONS composite (Fig 1) covers RA 114.57–118.86 and Dec 12.82–26.29 (3.9° × 13.5°).** This contradicts the
  caption's "field of view of 1°×10°", which only describes the nebula. The bright [O III] filament runs the full
  height of the image.
* **The Fig 1b white circle** is centred at (116.227, +15.934), 0.045° from HIP 37797. Its radius about HIP 37797 is
  2.37–2.49° (median 2.41°).
* **The Fig 2 fibres** extend to 2.495° from HIP 37797, with 95% inside 2.414°. This confirms the plate centre and the
  2.5° radius.
* **Fig 4 shows 173 [O III] fibres**, whereas the text says 279 spectra show relatively strong [O III]. They span
  Dec 13.51–18.35. Their relative intensity has a median of 34.7 (16–84%: 25.7–46.8).
* **Fig 5 (Hα V_LSR).** 157 of the 173 [O III] fibres are in the second, top scatter layer. The digitized Hα V_LSR has
  a median of 5.07 km/s (std 2.70) on that layer, against 5.42 km/s (std 3.33) for all other fibres. In other words,
  the Hα velocity on the filament does not differ measurably from the surroundings, consistent with the paper's
  red-line/background statement.
* **The background images in Figs 2, 4 and 5 are not registered to the plotted axes.** For Fig 4, positions implied by
  the axes differ from star-based positions by a median of 3.1′ and up to 10.4′. Measured on the Fig 4 image, the RA
  scale is off by about 5% (the residual slope is ~190″ per degree of RA). Fibre positions are plotted correctly, but
  their apparent alignment with the filament in those figures is only good to about 0.1°.
* **The P1–P4 outlines drawn on Fig 1b do not all match Table 1.**

  | Region | Offset of Table 1 position from drawn outline | Drawn outline |
  |---|---|---|
  | P2 | 0.04° | — |
  | P3 | 0.19° | — |
  | P1 | 0.48° | ellipse centre (116.606, +17.082), Dec 16.52–17.64; the Table 1 position sits on its southern tip |
  | P4 | 0.42° | box centre (115.445, +15.829); the Table 1 position (115.20, +16.18) lies outside the drawn box |

  On the filament, the Table 1 positions sit as follows: P1 on the E strand (offset −0.005°), P2 on the W strand
  (−0.04°), P3 0.115° west of the merged ridge, and P4 1.05–1.35° west of the strands, so P4 is valid as background.

### 4.2 Morphology (measured)
* **Main (E) strand = spine.**
  * Traced from Dec 12.00 (116.54) to Dec 28.25 (114.36), i.e. (l, b) = (208.15, +17.50) → (191.37, +21.85).
  * Path length is 16.94° and the end-to-end chord is 16.38°, consistent with NSNS's "~17°".
  * Per-point σ has a median of 0.022° and a maximum of 0.20°, at the junction.
  * Gaussian HWHM has a median of **0.138°** (16–84%: 0.076–0.191°), i.e. FWHM ≈0.28°. This depends on stretch and
    filtering, so treat it as ±30%.
  * Fainter emission continues to Dec ~11.3 in the south and ~28.9 in the north (NSNS wide peaks, S/N ~7), but is not
    in the spine.
* **Junction.** Between Dec 17.5 and 17.75 the E ridge steps 0.22° east (116.61 → 116.83), and it stays faint
  (S/N 6–12) up to Dec ~19. These points are flagged `junction`. The southern E sub-strand and the long northern arc
  may be separate, overlapping filaments.
* **W strand.**
  * Splits from E near Dec 15 and stays bright to Dec ~19.75, curving west (116.49 → 115.84). Length 4.84°.
  * Its separation from E is 0.08–1.08°.
  * It ends at the Hα/[S II] patch. P2 lies on it.
* **W2.** A faint strand at RA 115.1–115.35, Dec 21.5–24.0 (NSNS).
* **[O III] envelope.** Taking E+W+W2 ± their HWHM, the half-width has a median of 0.223° and a maximum of 0.83°
  (Dec ~19–20). The paper's "~1°" width corresponds to this envelope, not to a single filament.
* **NSNS O200.0b** (117.72, +12.52, "very faint [O III] filament"). A ridge at RA ~117.55–117.65, Dec 11.7–12.9 is
  detected in both NSNS views at 8–12σ prominence. This corroborates it; it is not traced in the JSON.

### 4.3 Curvature (derived)
* **Model-free.** The spine's maximum deviation from the great circle through its ends is **1.34°**, at Dec 19.75,
  bowing east. This sagitta implies an equivalent radius of about 25.8°.
* **Small-circle fit to the full E spine.** Centre **(94.50, +16.76)** with bootstrap 16–84% ranges RA 93.42–95.39 and
  Dec 16.52–16.96, i.e. (l, b) = (194.26, +0.44). **R = 21.33°** (20.48–22.35). The formal errors, scaled by χ²/dof,
  are 0.80°, 0.18° and 0.77°.
* **The circle is a poor model.** χ²/dof = 3744/63 and the rms residual is 0.150°. A great circle fits far worse
  (χ² = 40222 for 64 dof).

  | Segment fitted | Centre (RA, Dec) | R (16–84%) |
  |---|---|---|
  | North (Dec ≥17.75) | (94.18, +16.22) | 21.8° (19.6–25.5) |
  | South (Dec ≤17.5) | (69.35, +12.56) | 45.9° (22.5–71.4) |
  | W strand | (90.71, +13.40) | 25.0° (18.4–40.7) |

  The southern part is nearly straight, and the curvature is concentrated in the north. The arc is mildly S-shaped,
  with the kink at the junction. Bootstrap ranges understate the model error.
* **Literature centres measured against the E spine** (angular distance mean ± rms about the mean):

  | Candidate centre | Distance to spine | Comment |
  |---|---|---|
  | Geminga today (98.476, +17.770) | **17.56° ± 0.22°** (range 17.28–18.37°) | 3.9° from the fitted centre; nearly concentric. This is geometric only, since Geminga has moved ~17° since birth (NSNS) |
  | Monogem eROSITA centre (107.08, +13.5; 24° shell) | 9.35° (Dec 12) → 15.98° (Dec 28), rms 2.19° | not concentric. The arc lies inside the 12° rim south of Dec ≈21.25 and outside it further north |
  | PSR B0656+14 | rms 1.61° | — |
  | NSNS "C0" sphere centre (73.11, +6.57) | 44.09° ± 0.53° | — |

### 4.4 Hα patch (measured)
* **Half-maximum polygon** (NSNS Hα): centroid (115.474, +20.536), (l, b) = (199.46, +20.00). Major × minor axes
  1.45° × 0.31° at PA 163°, area 0.31 deg². Bounding box RA 115.23–115.79, Dec 19.83–21.30.
* **Checks.**
  * The [S II] half-maximum region overlaps the Hα one with IoU 0.65.
  * The 30%-of-peak contour is 2.24° × 1.28° (1.18 deg²).
  * The SWONS red channel gives centroid (115.33, +20.86) and 1.81° × 0.63°, with IoU 0.38. That image is truncated at
    RA ≈114.7 and contains continuum.
* **Comparison with literature.** The paper quotes 0.5° × 1°. NSNS O200.0a is 3.06° × 2.0° centred at (115.67, +21.07).
* **Size depends on the threshold.** The core is elongated along the arc (PA ~163°). It starts just north-west of
  the end of the W strand (115.84, +19.75), and its centroid lies 1.28° west of the E spine at the same Dec.

### 4.5 LAMOST [O III] fibres vs strands (derived)
* 87.3% of the 173 digitized [O III] fibres lie within 0.15° of the E or W ridge; the median |offset| is 0.061°.
  98 are nearest E and 75 nearest W.
* Intensity-weighted mean RA per 0.5° Dec bin follows E south of Dec 15 and lies between E and W for Dec 15–17.5.
  At Dec 18 it is on W (116.20, with W at 116.21). See `silk_trace.json → lamost_oiii_fibres_vs_strands`.

## 5. Figures
* `figures/geometry_trace.png`.
  * Left panel: the NSNS wide [O III] map with the E spine, the W and W2 strands, the [O III] envelope, the Hα patch,
    the LAMOST field, the small-circle fit with an arrow to its centre, NSNS O200.0b, and P1–P4.
  * Right panel: the registered SWONS composite over the LAMOST field, with the digitized [O III] fibres coloured by
    relative intensity, the strands, P1–P4 from Table 1 (×), the P1–P4 outlines as drawn in Fig 1b (dashed), and
    HIP 37797.
  * Both panels have an RA/Dec grid with galactic l (dashed) and b (dotted) lines.

## 6. Files
* Scripts (run in this order): `scripts/geometry_extract_pdf_vectors.py`, `geometry_register_figs.py`,
  `geometry_register_fig1.py`, `geometry_resample_fig1.py` (optional quick look), `geometry_build_maps.py`,
  `geometry_ridges.py`, `geometry_trace_fit.py`, `geometry_build_trace.py`, `geometry_make_figure.py`.
  Helper modules: `geometry_imreg.py`, `geometry_fig1_tools.py`, `geometry_nsns.py`.
  Also `geometry_preliminary_trace.py`, which wrote the early version.
* Derived data (`data/derived/`):
  * main outputs: `silk_trace.json`, `geometry_spine_points.csv`, `geometry_strand_fits.csv`, `geometry_ridge_peaks.csv`;
  * digitized paper data: `geometry_lamost_oiii_points_digitized.csv`, `geometry_lamost_fibres_fig2_digitized.csv`,
    `geometry_lamost_halpha_vlsr_fig5_digitized.csv`, the `geometry_fig{2,4,5}_*_pagecoords.csv` files;
  * registration: `geometry_pdf_axes.json`, `geometry_fig_registration.json`, `geometry_fig1_wcs.json`,
    `geometry_fig1_matched_stars.csv`, `geometry_fig1b_drawn_regions.json`;
  * maps: `geometry_maps.npz` with `geometry_maps_meta.json`, `geometry_fig1_resampled.npz`,
    `geometry_halpha_mask_nsns.npy`.
* Raw data (`data/raw/`): `paper_figures/`, `simbad/`, `vizier_tycho2/`, `nsns_simg/`. Total size is about 25 MB.

## 7. Limitations
* **Image quality.** All images are JPEG/PNG renderings with unknown non-linear stretches; none are calibrated FITS.
  Ridge centres are robust to monotonic stretches. Widths (HWHM), envelope and Hα-patch sizes are not, so treat them
  as ±30% and threshold-dependent.
* **Tracers differ.** The SWONS tracer (B−R colour index) is not pure [O III]: blue continuum and the Hα structure
  affect it. The Fig 4 background registration is only 19″ rms. The NSNS narrow view includes some continuum.
* **Strand identification is interpretive** at the E/W split (Dec 15–17.5) and at the junction (Dec 17.5–19).
  The `junction` flag and the large σ (0.14–0.20°) at Dec 18.5–18.75 reflect this.
* **Curvature.** The single-circle centre and radius describe a non-circular, S-shaped arc. Bootstrap ranges understate
  the model uncertainty, so use the section fits and the sagitta to gauge it.
* **Digitized LAMOST data.** These reproduce only what the paper plotted. The intensities are uncalibrated relative
  values, V_LSR is censored at 0 and 15 km/s, and only 173 of the 279 [O III] spectra mentioned in the text appear
  in Fig 4.
* **Single-source coverage.** North of Dec 26.25 only NSNS wide (40″/px) covers the arc, so σ there assumes 0.02°
  scatter between images.

## 8. Data gaps
* AstroBin plate-solve annotations and the SWONS original images were not retrieved. The geometry rests on the paper's
  embedded SWONS composite plus NSNS.
* Fesen et al. 2026 figures (Figs 31 and 38) were not separately registered; only their text was used. Their Fig 38 is
  NSNS data, which I used directly.
* SIMBAD's main CDS site (simbad.cds.unistra.fr) was unavailable, so I used the CfA mirror. VizieR TAP failed, so I
  used ASU on vizier.cfa.harvard.edu.
* The original LAMOST MRS-N fibre table, with true positions and line fluxes, was not used here (dr10.lamost.org is
  unreachable; the spectra teammate covers LAMOST). The digitized positions should be replaced by catalogue positions
  if they become available.
