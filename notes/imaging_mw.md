# imaging_mw: multiwavelength imaging counterparts of the Silk Nebula (UV, X-ray, IR/dust, HI, radio, gamma)

Topic tag `imaging_mw`, written 2026-11-16. All numbers below are **measured** by me from public survey data unless I mark them as *literature* or *assumption*.

**Geometry used.** `data/derived/silk_trace.json` **v1.0** (geometry teammate, file time 16:33), read at run time through `silk_geometry.load_trace()`:
- E-strand spine: Dec +12.0 to +28.25 (l 191.4–208.2, b +17.5…+21.8).
- Filament HWHM 0.138 deg, so FWHM f = 0.276 deg.
- The W/W2 strands and the H-alpha patch are overlaid on the maps.

**Offset convention.** The perpendicular offset from the spine is **positive to the EAST**. East is the convex, outer side. West is the concave side: it faces the spine's small-circle curvature centre (l,b) = (194.3, +0.5) and PSR B0656+14 (Monogem pulsar, l,b = 201.1, +8.3).

**Bug found and fixed.** A sign bug in the first pass of my offset code (east/west swapped) was caught by checking that P4 (west of the spine) gets a negative offset. All products were regenerated afterwards: P4 = −1.36 deg; P1 = −0.005; P2 = −0.46; P3 = −0.12.

---

## 1. Bottom line

1. **Soft X-rays: the Silk is the outer edge of the Monogem Ring X-ray emission (counterpart, morphological).**
   - In eRASS1 (0.2–0.6 keV) and ROSAT R12 (1/4 keV), the region inside (west of) the [O III] spine is brighter than outside:
     - eRASS1 0.2–0.6 keV: W/E = 1.32 (+1.59 cts s⁻¹ deg⁻², 2.4 σ against the parallel-shift null).
     - R12: W/E = 1.31 (+230 × 10⁻⁶ c s⁻¹ arcmin⁻²; 1.4 σ, as the map is noisier).
   - Over 11 segments of 1.5 deg along the spine, the X-ray half-drop edge lies **1.03 ± 0.20 deg inside** the spine (eRASS 0.2–0.6 keV) and **1.31 ± 0.22 deg inside** (R12).
   - Radial profiles around PSR B0656+14 in the two sectors holding the southern Silk (PA 0–40°):
     - the 0.2–0.6 keV ring peaks at r = 10.1 deg;
     - its outer half-drop is at r = 11.3–11.4 deg;
     - the Silk itself is at r = 11.2–11.6 deg.
   - So the X-ray emission fades out at the filament.
   - There is **no X-ray filament** on the Silk itself. On-filament contrast in all eRASS1 bands is |C| ≤ 1.3 σ.
   - Absorption cannot explain the interior excess: the interior has only 1.3 × 10¹⁹ cm⁻² less HI (GALFA), which changes 0.2–0.6 keV transmission by ≲ 7–15 %, against an observed +32 %.
   - This quantifies the qualitative statement of Fesen et al. (2026, arXiv:2608.26666, §6.2): "these [O III] filaments lie at the easternmost extent of the Monogem's soft X-ray emission".
2. **GALEX FUV: a weak excess at most, well below what a fast radiative shock should give (key test).**
   - On-filament FUV excess, measured four ways:
     - +47 ± 31 CU (per-tile local contrast, whole spine);
     - +75 ± 23 CU (per-tile, Dec 12–18);
     - +88 ± 30 CU (Murthy 2014 maps, whole spine);
     - +154 ± 46 CU (Murthy, south).
   - Errors are the rms of nulls with the trace shifted onto blank sky. The 2σ upper limits are 110–250 CU, i.e. ≲ 0.4–0.9 R for a 1549 Å line.
   - Expected FUV excess for a 120–140 km/s incomplete radiative shock making I(5007) = 10–13 R. The ratios are from Raymond et al. (1988) as quoted by Kim et al. (2010); the extinction correction is included:
     - **2,700–4,400 CU** with unscattered C IV;
     - **900–1,500 CU** with C IV suppressed 3.5× by resonance scattering;
     - **340–500 CU** with C IV suppressed 20×;
     - **210–300 CU** floor from O III] 1663 + He II 1640 alone.
   - The measurement lies **below the O III]-only floor** in three of four estimates, and **10–40× below** the predictions that include C IV.
   - So either the [O III] emission is not from a ≳ 100 km/s radiative shock, or I(5007) is several times below 10–13 R, or C IV is almost completely scattered or absorbed.
   - The NUV shows a comparable narrow excess (+64 ± 16 CU per tile; +113 ± 20 CU Murthy, 5.7 σ). The FUV/NUV colour of the excess (0.7–1.2) is not clearly bluer than the diffuse background (0.71). So the excess cannot be uniquely attributed to line emission, and some may be residual tile systematics.
3. **HI: a cavity inside and a wall ~0.9 deg outside the filament (ambiguous/suggestive).**
   - GALFA-HI DR2 at V_LSR +11…+16 km/s (the [O III] velocity) and +10…+18 km/s shows a low-column floor near and inside the spine.
   - It rises to a wall **0.86 ± 0.15 deg east (outside)** of the spine. The rise is found in 11/11 segments, with wall/floor = 1.8 (1.3–2.6).
   - The −1…+6 km/s gas (the H-alpha/[N II]/[S II] velocity) shows no comparable structure (wall found in only 6/11 segments).
   - There is **no HI filament on the Silk**: |C| ≤ 0.6 σ, and 2σ ΔN(HI) < 1.3 × 10¹⁹ cm⁻² at +11…+16 km/s.
   - The position–velocity diagrams show no clean expanding-shell arc.
4. **IR/dust: no counterpart.**
   - IRIS 100/60 µm, Planck 857/545 GHz and SFD show no filament.
   - 2σ: ΔI₁₀₀ < 0.19 MJy sr⁻¹ and ΔE(B−V) < 0.004 mag, i.e. N_H(dust) < 2.3 × 10¹⁹ cm⁻² using 5.8 × 10²¹ cm⁻² mag⁻¹.
   - The filament is **not** an ionized skin on a cirrus cloud edge. The local 100 µm level is low (3.2 MJy sr⁻¹), and cirrus lies 3–4 deg west and 3 deg east.
5. **Radio: no counterpart in total intensity; Stockert 1.4 GHz is ambiguous.**
   - Haslam 408 MHz: C = −0.07 ± 0.62 K, so Σ₁GHz < 1.3 × 10⁻²² W m⁻² Hz⁻¹ sr⁻¹ (2σ, beam-diluted; α = 0.5).
   - WMAP K and Planck 30 GHz give similar Σ limits.
   - NVSS and TGSS show nothing, but they resolve out the 17′-wide filament, so their non-detection (as in the paper) is weak evidence.
   - Stockert 1420 MHz shows +19 ± 5 mK on the spine (4.0 σ against nulls). If real, that is Σ₁GHz ≈ 3–5 × 10⁻²³. It sits on a W>E brightness step (+68 mK, 4.6 σ) also seen in Planck 30 GHz polarized intensity (4.2 σ). At 35′ resolution I cannot separate it from Galactic structure.
6. **H-alpha (Finkbeiner 2003, 6′):** no filament; 2σ < 0.43 R (beam-corrected). The diffuse H-alpha is 12 % brighter inside (1.2 σ).
7. **Gamma rays (Fermi-LAT SkyView bands 3–5) and the optional HAWC/Fermi TeV-halo geometry:** no filament contrast (|σ| ≤ 1.7). Only positions of PSR B0656+14 and Geminga are marked; no halo modelling was done.

**What this means for the paper's verdict** ("possible SNR; DIG or bubble not excluded").
- The multiwavelength morphology strongly favours the Silk being a **boundary of the Monogem Ring hot-gas region**, which matches Fesen et al. 2026: hot X-ray gas inside, the [O III] filament at the X-ray edge, more HI just outside.
- This argues against an isolated 28-pc SNR at 160 pc with its own interior. The Silk's concave side faces the Monogem interior, not an interior of its own.
- The GALEX limit argues against the paper's ~119 km/s *radiative-shock* excitation at the quoted 10–13 R. Two pictures are favoured instead:
  - (a) a slow or non-radiative/incomplete shock;
  - (b) photoionization or conduction at the hot-gas/neutral interface (e.g. a turbulent mixing layer), possibly with I(5007) lower than quoted.
- These inferences are mine and need spectroscopic and shock-model confirmation.

---

## 2. Per-band verdicts

C is the on-filament contrast: |offset| ≤ max(0.15°, beam/2) minus a local reference annulus on both sides, over the whole spine, with the null mean subtracted. (§1 quotes the raw Murthy FUV values, +88 and +154 CU; the null-mean-subtracted values are +97 and +154 CU.) σ is the rms of 14 nulls with the trace shifted by ±1.5…±4.5 deg. UL₂ is max(C,0) + 2σ (beam-diluted). Step is the mean over −2.5…−0.5 (west/inside) minus +0.5…+2.5 (east/outside), with significance relative to the null mean and rms. The full table is in `data/derived/imaging_mw_band_summary.csv` and `imaging_mw_verdicts.json`.

| Band (data) | Res. | Filament contrast C (σ) | 2σ UL / intrinsic | W−E step (σ) | Verdict |
|---|---|---|---|---|---|
| GALEX FUV, per-tile skybg | ~3′ | +47 CU (null rms 31; south +75 ± 23) | 110 CU ≈ 0.41 R at 1549 Å | — | weak excess / **limit far below radiative-shock prediction** |
| GALEX FUV, Murthy 2014 | 6′ | +97 CU (3.3σ); south +154 (3.3σ) | 157 / 178 CU | +57 (−0.7σ) | tentative weak excess |
| GALEX NUV, Murthy / tiles | 6′/3′ | +113 CU (5.7σ) / +64 ± 16 | 153 / 174 CU | +88 (−0.1σ) | weak excess, nature unclear (colour ≈ background) |
| ROSAT R12 (1/4 keV) | ~0.4° | −43 (−1.2σ) | 69 ×10⁻⁶ c/s/′² | **+230; W/E 1.31 (1.4σ)**; edge 1.3 ± 0.2° inside | **boundary counterpart; no filament** |
| ROSAT R45 (3/4 keV) | ~0.4° | +0.9 (+0.2σ) | 12 | +13; W/E 1.15 (1.6σ) | weak boundary, none on filament |
| ROSAT R67 (1.5 keV) | ~0.4° | −2.3 (−0.5σ) | 9.9 | +2.8 (0.3σ) | none |
| eRASS1 0.2–0.6 keV | 21′ (smoothed) | +0.05 (+0.3σ) | 0.33 cts/s/deg² | **+1.59; W/E 1.32 (2.4σ)**; edge 1.03 ± 0.20° inside | **boundary counterpart (Monogem rim); no filament** |
| eRASS1 0.6–1.0 keV | 21′ | +0.08 (+0.5σ) | 0.41 | +0.39; W/E 1.11 (1.4σ) | weak boundary |
| eRASS1 1.0–2.3 keV | 21′ | −0.02 (−0.2σ) | 0.19 | +0.05 (0.3σ) | none |
| eRASS1 O VII / O VIII narrow | 42′ | −1.3σ / +0.5σ | 0.0012 / 0.0012 cts/s/deg² | W/E 1.26 (2.0σ) / 1.14 (1.3σ) | weak boundary, none on filament |
| IRIS 100 µm | 4.3′ | +0.016 (+0.2σ) | 0.19 MJy/sr | −0.06 (−1.8σ) | **none** |
| IRIS 60 µm / Planck 857 / 545 | 4–5′ | ≤ +0.5σ | 0.04 / 0.19 / 0.06 MJy/sr | — | none |
| SFD E(B−V) | 6′ | +0.3σ | 0.004 mag | — | none (N_H < 2.3 ×10¹⁹) |
| Planck 353 GHz PI | ~15′ | +0.7σ | 8.6 µK | — | none |
| GALFA-HI +11…+16 km/s | 4′ | +0.2σ | 1.3 ×10¹⁹ cm⁻² | −1.2 ×10¹⁹ (−2.6σ) | **ambiguous: cavity + wall 0.86 ± 0.15° outside** |
| GALFA-HI −1…+6 km/s | 4′ | +0.6σ | 1.5 ×10¹⁹ | +1.2 ×10¹⁹ | none specific |
| GALFA-HI all (−30…+40) / HI4PI | 4′/16′ | ≤ 0.2σ | 4.1 / 3.9 ×10¹⁹ | interior −1.3 ×10¹⁹ | none on filament |
| Hα Finkbeiner 2003 | 6′ | +0.13 R (+0.9σ) | 0.43 R | +0.29 R (1.2σ) | none |
| Haslam 408 MHz | 51′ | −0.07 K (−0.1σ) | 1.2 K → Σ₁GHz < 1.3 ×10⁻²² | +1.3 K (−1.4σ) | none + limit |
| Stockert 1420 MHz | 35′ | **+19.5 mK (+4.0σ)** | 29 mK → Σ₁GHz ≲ 5 ×10⁻²³ | +68 mK (+4.6σ) | **ambiguous** |
| DRAO 22 MHz | ~1.5° | +0.5σ | — Σ₁GHz < 2.7 ×10⁻²² | — | none (weak) |
| Planck 30 GHz I / PI | 32′ | −0.2σ / +1.5σ | Σ₁GHz < 1.3 ×10⁻²² | PI step +4.2σ | none / PI ambiguous |
| WMAP K | 53′ | −0.1σ | Σ₁GHz < 1.0 ×10⁻²² | — | none |
| GB6 4.85 GHz | 3.5′ (filtered) | +0.2σ | 0.23 mJy/beam | — | none (scales > 20′ filtered) |
| NVSS 1.4 GHz | 45″ | −0.3σ (rms 0.34 mJy/beam at 1′) | 0.02 mJy/beam (3′ avg) | — | none, but **resolves out** the 17′ filament |
| TGSS 150 MHz | 25″ | +0.5σ (rms 2.2 mJy/beam at 1′) | 0.23 mJy/beam | — | none, interferometric |
| Fermi-LAT 0.3–1, 1–3, 3–300 GeV | ≥0.3° | −1.0, −1.5, −0.2σ | — | — | none |

Region medians (P1–P3 minus P4) for every band are in `imaging_mw_profiles.json` (`regions`). Examples:
- Murthy FUV: P1 − P4 = +279 CU and P3 − P4 = +268 CU. P4 is 1.36 deg west, and the FUV level there is affected by tile structure.
- eRASS1 0.2–0.6 keV: P1 − P4 = −1.4 cts/s/deg². P4 lies in the bright interior.

---

## 3. Methods

**Common grids** (`scripts/imaging_mw_common.py`):
- EQ: plate carrée, RA 108–125 × Dec +4…+30, 0.05° pixels.
- GAL: plate carrée, l 185–225 × b 0…+35, 0.1° pixels.
- All maps are reprojected with `reproject` (bilinear; adaptive for NVSS/TGSS; HEALPix for LAB).

**Trace-relative coordinates.** For each pixel I compute the signed perpendicular distance to the densified E spine (tangent plane) and the along-spine coordinate s. Pixels beyond the spine ends are excluded.

**Profiles, contrasts, steps and nulls** (`imaging_mw_profiles.py`):
- Binned medians at 0.1° (0.2° for beams ≥ 0.3°), by segment: whole spine, south (Dec 12–18), north (Dec 18–28), and ±0.6° along the spine around P1–P3.
- Contrast C and step S are defined in §2.
- Nulls repeat C and S with the trace shifted perpendicular by ±1.5…±4.5 deg. This samples unrelated sky with the same map systematics and large-scale gradients. Significance = (value − null mean)/null rms.
- These null errors are much larger, and more honest, than pixel-noise errors.
- Beam dilution of a Gaussian filament is D = f/√(f² + b²). Intrinsic limits are UL₂/D (`imaging_mw_summary.py`).

**X-ray cleaning** (`imaging_mw_xray_clean.py`):
- eRASS1 rate maps and RASS3 intensity maps: pixels > 4 robust σ above a 21′ running median are masked (2–5 % of pixels), then normalised Gaussian smoothing is applied (σ = 9′, FWHM 21′ on EQ; σ = 9′ on GAL).
- The O VII/O VIII narrow maps are sparse counts/exposure (40 % genuine zeros, which I keep). They are only smoothed (σ = 18′).
- RASS diffuse maps (Snowden et al. 1997; point sources already removed) are summed to R12, R45 and R67. Units are assumed to be 10⁻⁶ c s⁻¹ arcmin⁻²; the SkyView header has no BUNIT.

**Edges** (`imaging_mw_edges.py`):
- Along-spine X-ray edge = the first point east of the interior peak (searched −3.5…+0.5°) where the profile drops below the mid-level between that peak and the exterior median (+0.5…+3°), per 1.5° segment.
- HI wall = the first point (searched −1…+3°) where N(HI) rises above the mid-level between the floor and the wall peak.
- Monogem ring: radial profiles around PSR B0656+14 in 20° PA sectors (0.25° bins), with the outer half-drop measured beyond the ring peak.
- Caveat: the wall finder always returns a position (no formal null was run for the edge offsets). Their significance rests on the coherent ±0.15–0.2° scatter across 11 segments and on the maps.

**GALEX** (`imaging_mw_galex_tiles.py`, `imaging_mw_galex_tests.py`):
- 267 GR6/GR7 coadd pipeline skybg maps (127 FUV, 140 NUV; AIS median 218 s, plus MIS) with centres in RA 112.5–120.5, Dec 10.5–22.5.
- skybg is the pipeline's source-free sigma-clipped background, so it is star-masked at pipeline level.
- Maps are block-averaged to 1′ (≥ 80 % valid) within r ≤ 0.5° of the tile centre.
- Conversion to continuum units (CU) uses 1.40×10⁻¹⁵ (FUV) and 2.06×10⁻¹⁶ (NUV) erg s⁻¹ cm⁻² Å⁻¹ per cps, at λ 1528/2271 Å with 1.5″ pixels. That gives 2.04×10⁶ and 4.45×10⁵ CU per cps/pixel. The median FUV level is 879 CU (Murthy: 769 CU).
- Two analyses:
  - (i) a joint robust fit of tile offsets plus a common profile;
  - (ii) the per-tile local contrast (on |off| ≤ 0.15° minus 0.30–0.55° in the same tile), after dividing out the mean radial tile pattern (±1.5 % FUV, ±5 % NUV), with null shifts ±0.75…±3°.
- The Murthy (2014) HLSP Aitoff maps (0.1°) are analysed like the other bands.

**FUV prediction** (`imaging_mw_fig_galex.py`; all assumptions in `data/derived/imaging_mw_fuv_prediction.json`):
- I(5007) = 10–13 R. This is the paper's Table 3 conversion, reproduced by the paper-review teammate (`data/derived/paper_checks.json`: 9.9–12.8 R for AB, 60 Å).
- *Literature:* for 120–140 km/s incomplete shocks, O III] 1661,1666 / [O III] 4959,5007 ≈ 0.3 and C IV/O III] ≈ 10–13 (Raymond et al. 1988 models, as quoted by Kim et al. 2010, ApJ 722, 388, §4). I copied the paper to `data/raw/imaging_mw_literature/kim2010_1008.3240.pdf`.
- He II/O III] = 0.4–0.5 (same source).
- C IV resonance-scattering suppression of 1, 3.5 or 20. The 3.5 is from Sankrit & Wood 2001 as summarised in a web search, not read in full; the 20 is from Kim et al. 2010, region R2.
- *Assumptions:* W_eff = 300 Å (±30 %); relative FUV response 1.0 / 0.65 / 0.55 at 1549 / 1640 / 1663 Å, approximated from the GALEX effective-area curve (Morrissey et al. 2007), not from a tabulated file; A_FUV = 8.0 E(B−V) with E(B−V) = 0.04 (SFD along the spine, a full-column upper bound).
- These conversion assumptions carry a factor ≲ 2 uncertainty. That does not close the 10–40× gap to the C IV predictions, but it does matter for the comparison with the O III]-only floor.

**HI** (`imaging_mw_hi_build.py`, `imaging_mw_fig_hi.py`):
- 9 GALFA-HI DR2 *wide* cutouts (0.736 km/s channels, −29.8…+40.1 km/s, 4′ beam) on the native 1′ RA/Dec lattice, mosaicked by index arithmetic.
- Block-averaged to 3′. N(HI) = 1.823×10¹⁸ ∫T_B dv (optically thin).
- LAB 10 km/s HEALPix slices are used for the wide view.

---

## 4. Sources (all accessed 2026-11-16; per-file URL, bytes and SHA-256 are in each `data/raw/<archive>/MANIFEST.jsonl`)

| Archive / product | URL | Release / ref. | Files (size) |
|---|---|---|---|
| SkyView (RASS bkg R1–R7, RASS3 Int/Cnt, IRIS 12–100, SFD, Planck 857/545/353 PI&PA/030 I&PI, HI4PI & EBHIS N(HI), Haslam 408, Stockert 1420, GB6, DRAO 22, GEETEE 35, WMAP K, Finkbeiner Hα, Dame CO, Fermi 3–5, NVSS, TGSS) | https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl | SkyView "current"; survey references per SkyView pages (e.g. Snowden+97 ApJ 485,125) | `data/raw/skyview_imaging_mw/` 59 FITS (33 MB) |
| eROSITA-DE DR1 eRASS1 half-sky band maps 0.2–0.4, 0.4–0.6, 0.6–1.0, 1.0–2.3, 2.3–5 keV (3′, AIT) + O VII/O VIII narrow band (ZEA) | https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/HalfSkyMaps_dr1/maps_zheng/ | DR1 (Zheng et al. 2024, A&A 681, A77) | `data/raw/erosita_dr1/` 7 files (554 MB) |
| GALEX diffuse UV maps (Murthy 2014) | https://archive.stsci.edu/missions/hlsp/uv-bkgd/aitoff/ | HLSP uv-bkgd v1, GR6+GR7 (Murthy 2014 ApJS 213, 32) | `data/raw/mast_uvbkgd/` 2 files (104 MB) |
| GALEX coadd skybg maps | https://galex.stsci.edu/data/GR6/pipe/… (URLs from MAST CAOM TAP https://mast.stsci.edu/vo-tap/api/v0.1/caom) | GR6/GR7 pipeline | `data/raw/galex_mast/skybg/` 267 files (730 MB) + query tables |
| GALFA-HI DR2 wide cutouts | https://galfa-dr2-cutout-dptnqjtpkq-uc.a.run.app/cutout (service at https://www.galfa-hi.space, redirected from purcell.ssl.berkeley.edu; system of record Harvard Dataverse doi:10.7910/DVN/X4PQ1Z) | DR2 (Peek et al. 2018 ApJS 234, 2) | `data/raw/galfa_hi_dr2/` 9 tiles (229 MB) |
| LAB HI HEALPix slices | https://lambda.gsfc.nasa.gov/data/foregrounds/HI/lab_healpix.tar.bz2 | LAB (Kalberla et al. 2005), Land & Slosar repackaging | `data/raw/lambda_imaging_mw/` (141 MB tar; 7 slices extracted) |
| Kim et al. 2010 (G65.3+5.7 FUV) | https://arxiv.org/pdf/1008.3240 | ApJ 722, 388 | `data/raw/imaging_mw_literature/` |
| Literature read from teammates' captures | Fesen et al. 2026 arXiv:2608.26666; Knies et al. 2024 arXiv:2401.17289 (`data/raw/context_arxiv/`) | — | — |

*Literature context used:* Knies et al. 2024 report a Monogem kT = 0.14 ± 0.03 keV (0.18 keV in the east), age ≈ 1.2×10⁵ yr at D ≈ 300 pc, and a second candidate SNR G205.6+12.4 (hotter, ~5×10⁴ yr; radius 47 pc at 300 pc, about 9°). The Silk spine lies 11.1–16.5° (median 12.9°) from PSR B0656+14, and 5.6–16.5° (median 10.1°) from (205.6, +12.4). Its southern part (r 11.2–11.6° from the pulsar) coincides with the measured Monogem X-ray outer edge.

---

## 5. Figures (all `figures/imaging_mw_*.png`, 200 dpi, `silk_style`)

- `imaging_mw_eq_panels.png`: 12 bands on the RA/Dec grid (GALEX FUV/NUV Murthy, R12, R45, eRASS1 0.2–0.6 and 0.6–1.0 keV, IRIS 100, Hα, GALFA N(HI) at the [O III] and red-line velocities, Stockert 1420, Haslam 408), with the spine, W strands, Hα patch, P1–P4 and LAMOST field. The X-ray-bright region ends at the spine.
- `imaging_mw_gal_panels.png`: wide Galactic view (l 185–225, b 0–35) in 12 bands with PSR B0656+14, Geminga and the spine's small-circle fit. The Silk traces the northern/eastern rim of the Monogem Ring in R12 and eRASS1.
- `imaging_mw_monogem_edge.png`: eRASS1 0.2–0.6 keV with PSR-centred sectors, the radial X-ray profiles in the Silk sectors (Silk at the outer half-drop), and the X-ray-edge and HI-wall offsets along the spine.
- `imaging_mw_profiles.png`: cross-filament profiles (whole spine and P1/P2/P3, with the P4 level) for 14 bands plus the GALEX tile joint-fit FUV/NUV profiles.
- `imaging_mw_galex.png`: GALEX FUV/NUV skybg mosaics; per-tile contrasts against nulls; joint-fit profiles; and **expected radiative-shock FUV excess against measured (±2σ)**.
- `imaging_mw_hi_channels.png`: GALFA-HI channel maps from −10 to +30 km/s in 4 km/s bins, plus position–velocity diagrams across the spine (south and north) with the [O III] and red-line velocity bands.

---

## 6. Limitations and caveats

- **Trace dependence.** Every profile uses trace v1.0 (E strand). The W strand lies 0.2–0.5° west (P2 is −0.46°). The brightest narrow GALEX FUV band (RA ~116.0–116.4, Dec 14–18.5) lies nearer the W strand, so my on-region (±0.15°) partly misses it. Re-running `imaging_mw_galex_tiles.py`, `imaging_mw_profiles.py` and `imaging_mw_edges.py` picks up any trace update automatically.
- **GALEX systematics.** Within-tile flat-field and vignetting residuals and bright-star artefacts are comparable to the signal (null rms 16–46 CU). Airglow differences between visits are removed only as constant tile offsets. FUV coverage is patchy north of Dec 18. The FUV prediction depends on I(5007); the paper's magnitude system and background subtraction are unstated (see the paper-review note). If I(5007) were ~3 R, the O III]-only floor (~60–90 CU) would match the measured excess.
- **Edge offsets** come from a simple half-drop finder with no formal null. The X-ray edge positions scatter ±0.7° from segment to segment. The radial (PSR-centred) and perpendicular (spine-centred) estimates differ (Silk at the edge versus 1.0° outside the half-drop) because their geometries and smoothing differ. Both agree that the X-ray brightness declines across the ~1–1.5° just inside the filament and is at the exterior level by the filament.
- **X-ray depth.** eRASS1 is only 1 of 8 planned survey passes, so the narrow-band O VII/O VIII maps are count-limited. No spectral fitting was done. kT and absorption are taken from Knies et al. 2024 (literature).
- **Radio.** Single-dish beams (35′–53′) dilute a 17′ filament by 0.3–0.45. The Stockert +19 mK on-spine signal and the W/E step (also in Planck 30 GHz PI) may be Galactic structure or scanning effects. NVSS and TGSS are blind to ≳15′ structure.
- **HI.** N(HI) is optically thin. The GALFA eastern column (RA > 119.3°, Dec > 10.5°) is missing, so east-side profiles beyond about +2.5° in the south and +2.7° in the north rely on fewer pixels. The HI wall at +0.86° is inside the covered area.
- **H-alpha (Finkbeiner 2003)** has survey seams near Dec +15. Its resolution (6′) is adequate for a > 0.4 R filament, but it is not a deep test (the optical teammate covers deep Hα).

## 7. Data gaps

- **AKARI FIS** (WIDE-S, N160): SkyView returned HTML error pages for three requests. The files are renamed `*.INVALID_html_error.txt`. Not obtained.
- **WISE 12 µm WSSA** (Meisner & Finkbeiner): the host https://wise.skymaps.info responded, but I did not retrieve the tiles (time). IRIS 12 µm is dominated by zodiacal residual stripes here, so no useful 12 µm test was made.
- **GALFA-HI DR2**: three cutout tiles (RA 119.33–125, Dec 10.5–30) were refused with HTTP 429 (rate limit) on 8 attempts. A full single-request cutout returned HTTP 500.
- **HI4PI velocity cubes**: CDS (cdsarc) is unreachable, and the CfA VizieR mirror does not serve the FTP cubes. The MPIfR/AIfA profile servers returned 403. I used GALFA-HI and LAB instead.
- **Murthy 2014 2′ tables** (2.2 + 2.7 GB text) were not downloaded (size); I used the 0.1° Aitoff maps.
- **MAST JSON API** (mast.stsci.edu/api/v0/invoke) was blocked by the proxy; MAST TAP worked.
- **eRASS1 event files and tile images** were not needed, since the half-sky maps cover the field. A spectral or edge-hardness analysis with event data is a natural next step.
- **HAWC/Fermi TeV-halo geometry**: only source positions are marked; no halo maps were retrieved.
- **Planck/WMAP polarization angles and rotation measures** across the filament: not analysed beyond PI contrast.

## 8. Files

- Scripts (`scripts/`):
  - `imaging_mw_fetch.py`, `imaging_mw_fetch_galfa.py`, `imaging_mw_fetch_galfa_retry.py`, `imaging_mw_fetch_galex.py`, `imaging_mw_fetch_core.py`, `imaging_mw_galex_query.py`
  - `imaging_mw_common.py`, `imaging_mw_build_grids.py`, `imaging_mw_build_core.py`, `imaging_mw_xray_clean.py`, `imaging_mw_hi_build.py`
  - `imaging_mw_galex_tiles.py`, `imaging_mw_galex_tests.py`, `imaging_mw_profiles.py`, `imaging_mw_edges.py`, `imaging_mw_summary.py`
  - `imaging_mw_fig_panels.py`, `imaging_mw_fig_profiles.py`, `imaging_mw_fig_galex.py`, `imaging_mw_fig_monogem.py`, `imaging_mw_fig_hi.py`, `imaging_mw_quicklook.py` (diagnostic)
- Derived (`data/derived/`):
  - `imaging_mw_grids/` (all reprojected maps plus `grid_catalog.json`, the GALFA cube, GALEX mosaics)
  - `imaging_mw_profiles.json`, `imaging_mw_band_summary.csv`, `imaging_mw_verdicts.json`, `imaging_mw_edges.json`
  - `imaging_mw_galex_blocks.npz`, `imaging_mw_galex_profiles.json`, `imaging_mw_galex_tests.json`, `imaging_mw_fuv_prediction.json`
  - `imaging_mw_hi_pv.npz`, `imaging_mw_eq_tracecoords.npz`
- Re-run order: build_grids → build_core → xray_clean → hi_build → galex_tiles → galex_tests → profiles → edges → summary → fig_*.
