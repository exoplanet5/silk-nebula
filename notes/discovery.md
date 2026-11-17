# Silk Nebula: discovery history and observational characteristics

*Topic tag: `discovery`. Compiled 2026-11-16. All web sources were accessed on 2026-11-16 unless stated.*
*Machine-readable timeline: `data/derived/discovery_timeline.json` (39 events, each with a reliability flag). Geometry and brightness numbers: `data/derived/discovery_extents.json`.*

## 0. Summary

1. **At least three groups detected the filament independently in late 2023 and early 2024.** A fourth claim is hearsay.
   - **SWONS** (Chinese team; includes NAOC/LAMOST staff): detection in Nov (or Dec) 2023.
   - **S. Ziegenbalg / NSNS** (Germany): detection around Jan 2024 (self-reported).
   - **NHZ / Schaeffer–Fesen–Drechsler team**: imaged it as the eastern limb of the Monogem Ring; their own detection date is undocumented.
   - **S. Koposov in DESI data**: hearsay only.

   The AJ version of Wu et al. (published 2026-11-04) now says so explicitly: "*The SWONS team was not the sole discoverer of the Silk Nebula … each of these teams independently discovered the Silk Nebula as well.*" This sentence is **absent from arXiv:2610.13329v1**.
2. **Earliest public record.**
   - **First public image:** AstroBin `zh47r9`, posted by SWONS member Song-Yu Jia on **2024-01-14** ("Huge Undiscovered [OIII] Structure", 1°×10°, no coordinates).
   - **First public post with position and multi-line data:** Ziegenbalg's Cloudy Nights thread of **2024-06-11** (14° arc).
   - **First use of the name in English:** the QHYCCD press note, **2024-07-12**. In Chinese the name had already appeared in *China National Astronomy* 2024 No. 4.
   - **First journal paper:** **Ziegenbalg 2025, RNAAS 9, 227** (published 2025-08-28; "17°-large O III arc in Gemini"). This is 14 months before Wu et al. on arXiv.
3. **Quoted sizes (1°×10°, 1.5°×15°, 14°, 17°) are mostly differences in coverage, not disagreements about the object.**
   - The geometry teammate's spine (silk_trace v1.0) spans Dec 12.0–28.25. Its arc length is **16.9°** and its chord 16.4°.
   - Wu et al.'s own Fig. 1a, labelled "1°×10°", actually covers Dec 12.8–26.3. The spine is detected in it over **13.8°**.
   - The bright southern part up to the Hα patch is 9.8°; this is the "10°".
   - At the paper's 160 pc, 16.9° corresponds to 47 pc, not 28 pc.
4. **Surface-brightness estimates differ by a factor of about 6.**
   - Wu et al.: 23.6–23.9 mag/arcsec², which they quote as 10–13 R. I can reproduce that conversion only by assuming a 6 nm top-hat filter and AB magnitudes.
   - QHYCCD press note: 24.3 mag/arcsec², which is 3.3–6.7 R depending on the assumed bandwidth.
   - Ziegenbalg's NSNS estimates: 1.5–2.0 R.
5. **The sources disagree about the relation to the Monogem Ring.**
   - **Part of Monogem:** NHZ/Schaeffer (2024) and Fesen et al. 2026. Fesen et al. call it the "nearly vertical 15 degree long group of [O III] filaments" on the eastern/northern limb, shown in their Figs. 31 and 38. They never use the names "Silk" or "SWONS".
   - **Unlikely to be part of Monogem; possibly a Geminga SNR fragment:** Ziegenbalg (2025).
   - **Unnamed:** QHYCCD says only that it "intersects a hypothesized massive SNR".
   - **Independent ~160 pc SNR, yet "considered part of the Monogem Ring":** Wu et al. 2026. The paper does not reconcile these two statements.
   - **Geometric check against the eROSITA 24° circle** (measured spine): the southern 9.3° lies 2.4–2.9° inside the rim, and the northern 7.7° lies outside it (+4.0° at Dec 28). The spine crosses the rim at (116.71, +21.18). Its small-circle curvature centre (94.5, +16.8; r = 21.3°) is **12.6° from the Monogem centre** and 3.9° from Geminga. Caveat: this fit is driven by the northern half; see §6.
6. **None of the SNRs or candidates discussed by Fesen et al. is spatially adjacent to the Silk.**
   - Nearest: G190.4+12.5 is 9.3° from the spine; the edge gap is 6.2°.
   - G191.4+11.1: 10.6°. G206.6+6.1: 11.4°. G203.1+6.6: 11.9°.
   - G206.9+2.3: 15.2°. G195.9+2.2: 19.0°. G205.7−1.7: 19.3°.
   - The Silk is the easternmost structure in the Monogem field.
7. **DESI.**
   - Saydjari et al. 2026 (arXiv:2608.28740) maps only the dust-correlated Hα and [S II] in about 7.3° pixels. It does not map [O III] and never mentions Gemini, the Silk or Koposov.
   - The "Koposov found it in DESI data" statement exists only second-hand: S. Jackson, AstroBin `zwfplw`, 2025-05-02, quoting V. Belokurov.

---

## 1. Methods

- **Pages read in full** and saved as captures in `sources/web/discovery_*.md` (raw downloads in `data/raw/discovery_web/`): the QHYCCD press note, the AstroBin forum topic and posts, the AstroBin image records, the NSNS pages, the Cloudy Nights thread and the Chinese web pages.
  - Fetched with curl through the proxy where possible.
  - AstroBin pages return HTTP 403 to curl, so they were read through AstroBin's public JSON API (`/api/v2/forum/…`, `/api/v2/images/image/?hash=`, `/api/v2/nestedcomments/…`).
  - Cloudflare-protected pages were read with Tavily or Exa extracts.
- **Dating.** Dates come from machine timestamps wherever possible: API `published`/`created` fields, git commit times from the GitHub API for `stzie/simg.de-discoveries`, RSS `pubDate`, Crossref and DataCite received/accepted/issued dates, and arXiv submission history.
- **Papers.**
  - Fesen et al. 2026 (v2 PDF; v1 text in `data/raw/context_arxiv/`) and Saydjari et al. 2026 were read with `pdftotext`, and figure pages were inspected visually.
  - The AJ text of Wu et al. was captured from iopscience and diffed against arXiv v1 (`scripts/discovery_diff_aj_arxiv.py` → `data/derived/discovery_aj_vs_arxiv_diff.txt`).
- **Geometry** (`scripts/discovery_geometry.py`):
  - Ziegenbalg's catalogue ellipses (centre, major/minor axis, PA) were converted to major-axis end points.
  - The arc length, Monogem-rim crossings and neighbour distances were computed from the geometry teammate's measured spine (`data/derived/silk_trace.json` v1.0; I did not edit it).
  - Surface brightness (mag/arcsec²) was converted to Rayleigh under stated assumptions (§5).
- **Figures** (`scripts/discovery_figures.py`, silk_style, 200 dpi). Every figure was inspected.

## 2. Sources (URLs; accessed 2026-11-16)

| Source | What it gives | Capture |
|---|---|---|
| Wu et al. 2026, AJ 172, 298, doi:10.3847/1538-3881/aeb11c; arXiv:2610.13329v1 (2026-10-09) | SWONS discovery claim; LAMOST MRS-N spectra; 1°×10°; 23.6–23.9 mag/arcsec²; SNR? verdict | `sources/arxiv/…`, `sources/web/discovery_aj_vs_arxiv.md` |
| QHYCCD press note, https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/ (2024-07-12) | 1.5°×15°, 24.3 mag/arcsec², equipment, "hypothesized SNR" | `discovery_qhyccd_pr.md`; images in `data/raw/discovery_web/qhyccd_*.jpg/png` |
| AstroBin forum topic 126958 (https://app.astrobin.com/forum/topic/126958), post 166798 | QHYCCD repost (Yilia Yu); Ziegenbalg priority note; "part of the monogem" (T. Schaeffer); link to zh47r9 | `discovery_astrobin_forum_126958.md` |
| AstroBin zh47r9 (Song-Yu Jia, 2024-01-14) https://www.astrobin.com/zh47r9/ | Earliest public image | `discovery_astrobin_images.md`, `data/raw/discovery_web/astrobin_zh47r9_hd.jpg` |
| AstroBin zwfplw (S. Jackson / Hunting Horn, 2025-05-02) https://www.astrobin.com/zwfplw/ | 15°, 48 h RASA 11 [O III] mosaic; "3 independent discoverers" incl. Koposov (hearsay) | same |
| AstroBin y32271 (R. Raupach, 2025-02-05) https://www.astrobin.com/y32271/ | "DNA Nebula" [O III] filaments at 06:33 +18:35. **Not the Silk** (≈18° W) | same |
| AstroBin 1d8ivk (M. Drechsler) — linked in the QHYCCD note | M31 [O III] arc SDSO-1; illustrates the Drechsler team's [O III] work. **Unrelated to the Silk** | same |
| AstroBin forum 221269 (2026-02-28) | "AstroBin first?" discussion | same |
| Ziegenbalg NSNS page https://simg.de/nebulae3/gem-01.html; index https://simg.de/nebulae3/index.html; instruments https://simg.de/nebulae3/instruments.html; news/RSS | NSNS discovery, line intensities, Monogem/Geminga discussion, data-release dates | `discovery_nsns_gem01.md` |
| GitHub stzie/simg.de-discoveries (commits 2e1a2ac0, f8273e1d) | Dated catalogue entries K1a → O201.5 → O200.0 | same |
| Ziegenbalg 2025, RNAAS 9, 227, doi:10.3847/2515-5172/adfec7 | First journal description (17°×3°) | same |
| Zenodo snapshot doi:10.5281/zenodo.17215566 (DataCite metadata only) | Persistent copy of the gem-01 page (2025-09-27) | same |
| Cloudy Nights topic 925409 (pages 1–2) https://www.cloudynights.com/forums/topic/925409-possible-discovery-of-a-large-emission-nebula-in-gemini-and-a-tiny-one/ | Public priority and Monogem debate, 2024–2025 | `discovery_cloudynights_925409.md` |
| Fesen et al. 2026, arXiv:2608.26666 v1/v2 | Monogem optical filaments; the 15° group | `discovery_fesen2026_excerpts.md` |
| Saydjari et al. 2026, arXiv:2608.28740 | DESI sky-fibre DGL maps | `discovery_saydjari2026_desi_skyfibres.md` |
| G.-Y. Sun blog https://sunguoyou.lamost.org/sichou.html; Xinmin/Shangguan https://www.shobserver.com/news/detail?id=1214625; ifeng https://feng.ifeng.com/c/8y63sUUKpcN | Chinese-language dates and naming; press coverage | `discovery_chinese_sources.md` |
| Search and access log | Queries, null results, blocked sites | `discovery_search_log.md` |

## 3. Dated timeline (condensed; full version with URLs in the JSON)

Reliability: **H** = primary record with a machine timestamp; **M** = self-reported or approximate; **L** = hearsay.

| Date | Event | Rel. |
|---|---|---|
| 1981 | Monogem Ring recognised as an X-ray ring (HEAO-1; Nousek+81) | H (context) |
| 2007-08 | Kim+07 FIMS/SPEAR FUV: C IV/C III region "R1" (l ≈ 201–208°) roughly coincides with the southern arc; no O III in FUV | M |
| 2003–2012 | GALEX FUV said to show a thin filament along the arc (Ziegenbalg; a CN user) — **not verified here** | L |
| 2018-11 / 2022-12 | NSNS starts; first NSNS [O III] camera (2nd 2023-04; 3rd = full array 2024-01) | H |
| early 2023 | SWONS survey starts (50 mm f/0.95 + QHY1920M, 6 nm [O III], ~26°×16°) | M |
| **2023-11** (paper) / 2023-12 (G.-Y. Sun) / "late 2023" (QHYCCD) | SWONS detects the N–S [O III] filament in the Gemini winter field | M |
| ~2024-01 | NSNS: "arc detected in the preview images of the first clear night after the camera array reached its full configuration" (HTML-commented text) | M |
| 2024-01-10/11 | SWONS follow-up, Epsilon-180ED + QHY600PH, [O III] and Hα (AstroBin log, partial) | M |
| **2024-01-14** | **AstroBin zh47r9**: first public image, 1°×10°, no coordinates. Only comment: "where in the sky is it?" | H |
| 2024-02-01 | LAMOST MRS-N plate on HIP 37797 (RA 116.2525, Dec +15.8950), 3×900 s | H |
| 2024-04 (No. 4) | Name 丝绸星云 and a short description in *China National Astronomy* (not seen) | M |
| **2024-06-11** | Ziegenbalg (stefan822) posts on Cloudy Nights: 14° arc, "I would say it is a supernova remnant" | H |
| 2024-06-12 | T. Schaeffer (NHZ): "You stumbled upon a part of the Monogem SNR". Ziegenbalg: "not a new discovery, just the first (public) optical picture of a part of the Monogem ring" | H |
| 2024-06-15/18 | Ziegenbalg reverses: the arc lies outside the 18.4° primary ring and he will treat it as a new discovery | H |
| 2024-06-28 | NSNS page "Detection of an OIII arc in Gemini" (RSS, 14°); catalogue K1a "Probably a supernova remnant" | H |
| **2024-07-12** | QHYCCD press note and AstroBin forum 126958 (name "Silk Nebula", 1.5°×15°, 24.3 mag/arcsec²). The same day a user points out that Ziegenbalg has the same structure | H |
| 2024-07-13 | T. Schaeffer: "it's part of the monogem" | H |
| 2024-10-01 | NSNS DR0.1 | H |
| 2025-03-04 | M. Bartels reports a visual detection of the "OIII arc in Gemini" | M |
| 2025-05-02 | S. Jackson AstroBin zwfplw "15° long Silk Nebula"; Koposov/DESI claim | H (post) / L (claim) |
| 2025-05-22 | NSNS DR0.2 ([O III], [S II]); faint southern extension O201.5b | H |
| 2025-06-11 | Ziegenbalg DR0.2 update: 17°, northern part outside Monogem, Geminga hypothesis; catalogue O201.5 → **O200.0** | H |
| **2025-08-28** | **Ziegenbalg 2025, RNAAS 9, 227**: first journal paper on the object (rec. 08-01, acc. 08-06) | H |
| 2025-09-27 | Zenodo snapshot of the NSNS page | H |
| 2026-07-30 | Wu et al. received by AJ | H |
| 2026-08-27 (v2 09-07) | Fesen et al. arXiv:2608.26666: the 15° [O III] group is assigned to the Monogem limb | H |
| 2026-08-28 | Saydjari et al. DESI sky-fibre paper (no Silk content) | H |
| 2026-10-05 / 10-09 | Wu et al. accepted / arXiv v1 | H |
| **2026-11-04** | AJ 172, 298 published, with the new co-discoverer acknowledgment | H |
| 2026-11-16 | Chinese press (Xinmin/Shangguan, ifeng, Sina): "discovered by Chinese amateurs in early 2023", SNR presented as established; other discoverers not mentioned | H |

## 4. Discoverers, precovery and credit

**Independent detections** (ordered by documented date):

1. **SWONS** (Jia Songyu, Zhang Chao, Wu Chaojian, Zhang Weiguo, Qiu Hongyun, Sun Guoyou, per G.-Y. Sun's blog).
   - Detected Nov or Dec 2023 in a >5000 deg² [O III] mosaic.
   - Confirmed with the same 50 mm system (+1 h; the press note instead says a 1.5 h, 9-panel run with the 180 mm reflector) and then 36 h [O III] + 24 h Hα with a 180 mm f/2.8 reflector and QHY600M.
   - Two members (Zhang Chao, Wu Chaojian) are NAOC LAMOST MRS-N scientists and the paper's corresponding authors.
   - Zhang Chao is already credited on the January 2024 AstroBin image (with Jia, Zhang Weiguo, Qiu and Sun). Wu Chaojian appears in the May 2024 member list.
   - So SWONS was an amateur–professional collaboration from the start. This explains why LAMOST pointed at the field only about 2.5 months after the detection.
2. **Ziegenbalg (NSNS, Pulsnitz).**
   - Detection in about January 2024 (self-reported, exact night unknown).
   - First public posting with position, orientation and [O III]/Hα/[S II] data on 2024-06-11.
   - First peer-visible publication: RNAAS, 2025-08-28.
3. **NHZ / Schaeffer–Fesen–Drechsler team.**
   - In June 2024 T. Schaeffer said NHZ was building a full Monogem image and that "one already exists, but from another astrophotographer who didn't make it public" (unverifiable).
   - Their 216 h [O III] Monogem mosaic (Fesen+26 Table 4: Schaeffer, Björk, Body, Kottary, Arora, Sparkman) is undated. It was published on arXiv on 2026-08-27.
4. **Koposov (DESI)**: hearsay only (§7).

**"Precovery".**
- No optical [O III] image of the filament earlier than the SWONS November 2023 data was found. Ziegenbalg (CN, 2024-06-11) notes that the Sivan (1974) Hα atlas shows nothing at (l, b) = (200°, 20°).
- In the professional sense, the earlier records are at other wavelengths:
  - FIMS/SPEAR far-UV C III/C IV emission from region R1 (Kim et al. 2007, ApJ 665, L139). O III was not seen there, and the resolution is degrees.
  - GALEX FUV, where the filament is reportedly visible as a thin streak (unverified here; for the UV imaging teammate).
- NSNS [O III] cameras have operated since December 2022, so NSNS data of the field may predate the first SWONS detection. Ziegenbalg does not claim that.
- The Nov 2022 NSNS "Monoceros to Gemini" mosaic (Hα + broadband continuum, centre 6h46 +5°, 30°×28°) reaches only the edge of the Silk field and has no [O III] channel. It was not inspected; that belongs to the imaging-precovery teammate.

**Credit disputes and inconsistencies.**
- The arXiv v1 Conclusion says "discovered and named by the SWONS team" and "detected for the first time to its full extent". It credits no other discoverer, although Ziegenbalg's RNAAS preceded it by 14 months.
- The AJ version keeps those sentences but adds the acknowledgment of the Fesen/Drechsler/Schaeffer team and of Ziegenbalg, and adds Ziegenbalg 2025 to the references.
- Fesen et al. 2026 never mention SWONS, the name "Silk" or LAMOST. They credit the NSNS images ("invaluable") and other amateurs, and assign the filament to the Monogem Ring.
- In June 2024, Cloudy Nights users (T. Schaeffer; "yuzameh") argued that the arc was "nothing new" because it is part of the known Monogem SNR. Ziegenbalg's position: an emission nebula never published in visible light is a discovery even if it belongs to a known superstructure.
- AstroBin forum 221269 (2026) called the 2025 Hunting Horn image the only Silk image on AstroBin. That ignores SWONS's January 2024 post, which has a different title.
- The 2026-11-16 Chinese press dates the detection to "early 2023". That contradicts the paper (November 2023) and the G.-Y. Sun blog (December 2023).

## 5. Observational characteristics (as reported, with checks)

### 5.1 Names and designations
- **"Silk Nebula" / 丝绸星云** (SWONS; in print by *China National Astronomy* 2024 No. 4; in English in the 2024-07-12 press note).
- **"Huge Undiscovered [OIII] Structure"** (AstroBin, 2024-01-14).
- **NSNS K1a / K1b / K1c** (2024-06-28) → **O201.5** (2025-05-22) → **O200.0** (with O200.0a Hα+[S II]-only part, O200.0b faint southern extension, and the eastern filaments O205.0, ex-K1c) on 2025-06-11.
- **"OIII arc in Gemini"** (Ziegenbalg; Fesen et al. quote it as "the Gemini OIII Arc").
- Fesen et al. have no designation of their own: "NE filaments of the Monogem Ring SNR" (Fig. 38).
- No SIMBAD or Green-catalogue identifier was found.

### 5.2 Morphology
All sources agree on the basic picture:
- A long, nearly N–S [O III] ribbon of multiple parallel strands.
- Brightest and most structured in the south, where it is [O III]-only, with no Hα or [S II] counterpart.
- A faint arcuate northern continuation.
- A compact Hα + [S II] feature at about Dec +20.5.

Ziegenbalg's three-zone description:
- South: [O III] only.
- Centre (≈2.5–3°): only Hα + [S II].
- North: all three lines.

Fainter filaments east of the arc that "seem to cross" it (O205.0) are assigned by Ziegenbalg to Monogem.

Measured spine (geometry teammate, silk_trace v1.0):
- Main (E) strand: HWHM 0.138° (16–84 %: 0.076–0.191°).
- E+W envelope: median half-width 0.22°, maximum 0.83°.
- A bright W strand splits off near Dec +15 and ends at the Hα patch.
- Hα patch, half-maximum contour: centroid (115.47, +20.54), 1.45°×0.31°. The 30 % contour is 2.2°×1.3°.

So the quoted widths (1°, 1.5°, 2.5–3°) describe the whole multi-strand envelope. A single strand is only about 0.3° FWHM.

### 5.3 Quoted angular extents and how they reconcile

| Source (date) | Length × width | Notes |
|---|---|---|
| SWONS AstroBin (2024-01-14), Wu+ 2026 | 1° × 10° | Wu+ Fig. 1a footprint actually spans Dec 12.8–26.3; spine detected over **13.8°** in it |
| Ziegenbalg CN/NSNS K1a (2024-06) | 14° (14.3° × 2.4°, PA 175°) | Ends (116.62, +11.52) – (115.30, +25.76) |
| QHYCCD PR (2024-07-12) | 1.5° × 15° | — |
| Jackson AstroBin (2025-05) | 15° | Field 18.9° × 3.3° centred (116.02, +16.55) |
| NSNS O200.0 (2025-06), RNAAS (2025-08) | 17° × 2.5–3° (16.9° × 3.0°, PA 173°) | Ends (116.50, +11.54) – (114.28, +28.34); centre 07:41:48 +19:56 (RNAAS gives +19°54′) |
| Fesen+ 2026 | "nearly vertical 15 degree long group" | — |
| **Measured spine v1.0** | **16.94° arc (16.38° chord), Dec 12.0–28.25** | Bright S part (Dec 12.0–21.5) 9.8°; faint N arc 7.2°. NSNS wide image alone reaches Dec 28.25; NSNS narrow 25.5; SWONS 26.0 |

Physical length (arc, ±~5 % for curvature and end uncertainty): **47 pc at 160 pc**, 77 pc at 260 pc (Ziegenbalg's Geminga-SNR geometry), 89 pc at 300 pc (Monogem distance). The 10° extent gives 28 pc at 160 pc, which is the paper's "28 pc". The paper's Fig. 1a and the NSNS data both show more than its quoted 10°.

**The paper's "centre" (116.68, +16.42) is not the centre of the object.** It sits in the LAMOST/P1 region, 3.8° south of the spine midpoint.

### 5.4 Brightness

| Source | Quoted | → Rayleigh [O III] 5007 (AB, all in-band flux = line) |
|---|---|---|
| Wu+ Table 3, P1 / P2 / P3 | 23.60 / 23.87 / 23.85 mag arcsec⁻² ("~10–13 R") | 12.7 / 9.9 / 10.1 R for a 6 nm band; 10.6 / 8.3 / 8.4 (5 nm); 6.4 / 5.0 / 5.1 (3 nm) |
| QHYCCD PR, brightest region (V-calibrated) | 24.3 mag arcsec⁻² | 6.7 (6 nm) / 5.6 (5 nm) / 3.3 (3 nm) R |
| NSNS (Ziegenbalg; estimated transmissions, background-limited) | K1a 2.04 R; K1a′ 1.53 R | — |
| Fesen+2024 [O III] detection limit (as quoted by Wu+) | 0.5–5 R | — |

- My conversion reproduces Wu et al.'s 10–13 R only for a 6 nm (SWONS filter FWHM) top-hat bandwidth, an AB zero point, and no continuum in the band. The deep follow-up used a 5 nm filter (QHYCCD note).
- Wu et al. themselves note that P1, faintest in the image, is "brightest" in the photometry because of stellar contamination.
- The independent NSNS estimate is about 5–6× lower.
- **The absolute [O III] surface brightness is therefore uncertain by a factor of about 5 (≈1.5–13 R).** The paper's comparison "well above the 0.5–5 R Fesen limit" is not robust.
- Neither the LAMOST MRS-N (no flux calibration) nor the press-note number resolves this. A flux-calibrated spectrum (SDSS/DESI sky fibres) or a calibrated [O III] image would.

Other NSNS line intensities (R, Hα / [O III] / [S II]): K1a 0.02 / 2.04 / 0.06; K1a′ 1.68 / 1.53 / 1.14; K1b 2.20 / 0.63 / 2.27; K1c 0.24 / 0.72 / 1.54. K1b gives [S II]/Hα ≈ 1.0, as Fesen et al. note ("I(Hα)/I([S II]) ratio around 1.0"). Wu et al.'s LAMOST red-arm ratios are [S II]/Hα 0.83–0.93 on the nebula and 0.78 in background P4.

## 6. Relation to the Monogem Ring as stated by each source, and a geometric check

| Source | Statement |
|---|---|
| SWONS AstroBin (2024-01) | "previously unknown huge structure"; no Monogem mention |
| QHYCCD PR (2024-07) | "located on the outermost edge of a star-forming region. It also intersects with a hypothesized massive supernova remnant … since the 1980s … only low-resolution imaging data in X-rays"; nature undetermined |
| Ziegenbalg (CN 2024-06-11 → 06-18) | First guesses an SNR. Then accepts "first (public) optical picture of a part of the Monogem ring". Then: the arc lies outside the 18.4° primary ring, its implied circle centre is far from PSR B0656+14, and the implied diameter is about 60° |
| T. Schaeffer / NHZ (2024-06-12, 07-13) | "part of the Monogem SNR … 29° … full loop"; "it's part of the monogem" |
| Ziegenbalg NSNS page + RNAAS (2025) | Partially overlaps Monogem. The northern part "clearly lies outside", so "unlikely that the arc is part of the Monogem Ring SNR". The eastern faint filaments are probably Monogem. He notes that the line ratios correlate with position relative to the Monogem rim ([O III]-only inside, Hα+[S II] on the boundary, all lines outside). Speculative: a fragment of a Geminga-progenitor SNR at about 260 pc (~80 pc long) |
| Fesen et al. 2026 | The 15° [O III] group lies along Monogem's eastern and northern limbs, "at the easternmost extent of the Monogem's soft X-ray emission". The NE Hα/[S II] filaments mark where the shock meets a denser cloud; [O III] traces shocks faster than about 100 km/s. Mentions Ziegenbalg's Geminga alternative |
| Wu et al. 2026 | Possible independent SNR at about 160 pc (R ≈ 30 pc, 6.8×10⁴ yr, 2×10⁵¹ erg), **and** "has also been detected in Fesen's project … considered to be part of the Monogem Ring". The paper does not reconcile this with Monogem's ~300 pc distance, 8–12×10⁴ yr age and ~65 pc radius |
| R. Raupach (y32271) | Different [O III] filaments at about 06:33 +18:35, proposed as a Monogem optical component (not the Silk) |

**Geometric check** (measured spine vs Knies+2024 circles as digitised in Ziegenbalg's viewer links):
- **Single-SNR model (centre 107.08, +13.50, r = 11.97°).**
  - Spine distance from the rim: −2.62° at Dec 12, −2.43° at Dec 16, −0.50° at Dec 20, +1.45° at Dec 24, +4.01° at Dec 28 (negative = inside).
  - Crossing at (116.71, +21.18); 9.25° of arc inside, 7.68° outside.
  - So the filament is not parallel to the X-ray rim. It is inclined to it, cutting from about 2.6° inside to 4° outside.
- **Two-SNR model.**
  - The spine never enters the western "Monogem" circle (≥ 0.98° outside).
  - It crosses the G205.6+12.4 circle at (116.99, +19.01), with 7.2° of arc inside.
- **Curvature** (teammate's small-circle fit to the full spine): centre (94.50, +16.76), r = 21.3° (bootstrap 20.5–22.4°), (l, b) = (194.3, +0.4).
  - Distances from that centre: Monogem centre 12.6°, PSR B0656+14 10.4°, Geminga 3.9°, Ziegenbalg's Geminga-SNR centre 6.5°.
  - **Caveat:** the southern half alone (Dec ≤ 20) is nearly straight (fitted r ≈ 124°), and the full-spine fit has χ²/dof ≈ 59 (rms 0.15°). The curvature centre is set by the northern arc and should not be over-interpreted.
  - What the geometry does show is that the filament is not concentric with Monogem (it would need r ≈ 12°, centred 12.6° away), and that its northern 7–8° lie outside any published Monogem outline. This supports Ziegenbalg's objection, without proving his Geminga scenario.

## 7. DESI / Koposov

- **Only source:** S. Jackson, AstroBin zwfplw description (2025-05-02): "Prof Vasily Belokurov told me about how Prof Sergey Koposov found it by ingenious analysis of DESI data." No date, method, data product, talk or paper was found.
- **Saydjari et al. 2026 (arXiv:2608.28740, 10.8 M DESI Y3 sky fibres)** contains no emission-line map of the Silk.
  - Its spatial maps are the dust-correlated (IRAS-correlated) DGL, in NSIDE = 8 (~7.3°) pixels, for Hα and [S II] only.
  - Doubly ionised lines are "low signal to noise", so mapping them is "left to future work".
  - A dust-uncorrelated [O III] filament would not appear in these products.
- The "Koposov in DESI data" claim therefore stays unverified (reliability L). A direct test with the raw DESI DR1 sky-fibre spectra in the Silk box belongs to the DESI teammate.
- Precedent: Hewett et al. (2003) found the >2°-wide PN around PG 1034+001 from zero-velocity [O III] in SDSS spectra, including sky fibres (ING Newsletter 8, https://www.ing.iac.es/PR/newsletter/news8/science2.html).

## 8. AJ vs arXiv v1 (details in `sources/web/discovery_aj_vs_arxiv.md`)

- AJ dates: received 2026-07-30, accepted 2026-10-05, published 2026-11-04. AJ 172, 298, CC-BY.
- **Added in AJ:**
  - The acknowledgment naming the Fesen/Drechsler/Schaeffer team and Ziegenbalg as independent discoverers.
  - References to Ziegenbalg 2025 (RNAAS 9, 227) and Lumbreras-Calle et al. 2025.
  - The v1 "Slane et al. 2013" citation is corrected to Yamaguchi et al. 2014.
- **Unchanged:** all numerical results, the 1°×10° extent, and the "discovered and named by the SWONS team" wording in the Conclusion.
- The Fesen et al. citation (arXiv:2608.26666, posted 2026-08-27) must have been added after submission.
- OpenAlex shows 0 citing works as of 2026-11-16.

## 9. Figures

| File | Caption |
|---|---|
| `figures/discovery_timeline.png` | Dated public records 2023–2026 in six lanes (SWONS/LAMOST, Ziegenbalg/NSNS, Fesen/NHZ, other amateurs, DESI, press). Marker shape shows reliability (● timestamped, ■ self-reported, ✕ hearsay). |
| `figures/discovery_extents_map.png` | RA/Dec map comparing the NSNS 2024 (14°) and 2025 (17°) catalogue ellipses with sub-parts, the measured [O III] spine (silk_trace v1.0) with the W strand and Hα patch, P1–P4, the LAMOST plate, Monogem circles (single and two-SNR, Knies+2024), Fesen+2026 SNRs/candidates and pulsars. ◆ marks where the spine crosses the X-ray rims. |
| `figures/discovery_quoted_extents.png` | Left: quoted angular length vs publication date (marker size ∝ quoted width; right axis = length at 160 pc). Right: [O III] surface brightness in Rayleigh for the quoted mag/arcsec² values under 6/5/3 nm bandwidth assumptions, the NSNS estimates and the Fesen+2024 limit. |

## 10. Limitations

- Event dates for detections (as opposed to postings) are self-reported: SWONS Nov vs Dec 2023, NSNS "first clear night after Jan 2024". Neither the SWONS discovery frames nor the NSNS preview frames are public.
- I did not see *China National Astronomy* 2024 No. 4, so the first printed use of the name is dated only to the issue.
- The AJ text was read from an HTML extract. The rendering of tables and equations is incomplete, and the location of the Ziegenbalg 2025 in-text citation was not found.
- The AstroBin zwfplw footprint in the map is approximate: the rotation convention of AstroBin's `orientation` was assumed.
- Monogem circles are those Ziegenbalg digitised from Knies+2024 (and Fesen's 24° circle), not re-derived from eROSITA data. The X-ray rim is not a sharp circle.
- The Rayleigh conversions assume an AB zero point, a top-hat bandwidth and zero continuum, so they are indicative only.
- Neighbour separations use Fesen+2026 Table 1 centres (to the nearest arcmin) or the Galactic names (G203.1+6.6, G206.9+2.3).
- **Fesen+2026 inconsistencies:** G190.5+25.3 in the abstract vs G190.5−25.3 in Table 1 (RA 04:40, i.e. b < 0, far from the Silk), and G191.3+11.1 vs G191.4+11.1.
- **Silk–Monogem crossing points:** the trace is preliminary in places. Single-source points north of Dec 26.25 come from the NSNS wide image only.

## 11. Data gaps

- **NASA ADS:** not queried (no API token). Citation and full-text searches used OpenAlex, Crossref and arXiv instead.
- **Semantic Scholar:** HTTP 429.
- **zenodo.org:** HTTP 403 ("unusual traffic"). Snapshot metadata came from DataCite; files not downloaded.
- **AstroBin web pages:** HTTP 403 to curl; JSON API and extracts used instead.
- **vizier.cds.unistra.fr:** ASU queries return 404; the CfA mirror was used.
- **Koposov/Belokurov:** no primary statement found.
- **English-language news** (Sky & Telescope, Universe Today, phys.org, space.com): none found as of 2026-11-16.
- **NAOC press release:** not found.
- **Mel Bartels' visual observing report:** not found.
- **NSNS full-resolution FITS/HiPS and the Hunting Horn FITS mosaic** (Google Drive): not downloaded. These are for the imaging teammates.
- **GALEX FUV and FIMS/SPEAR coincidence:** asserted by Ziegenbalg, not checked here (UV teammate).

## 12. Files produced (topic `discovery`)

- **Notes and data:** `notes/discovery.md`; `data/derived/discovery_timeline.json`; `data/derived/discovery_extents.json`; `data/derived/discovery_aj_vs_arxiv_diff.txt`.
- **Scripts:** `scripts/discovery_geometry.py`, `scripts/discovery_figures.py`, `scripts/discovery_make_captures.py`, `scripts/discovery_diff_aj_arxiv.py`.
- **Figures:** `figures/discovery_timeline.png`, `figures/discovery_extents_map.png`, `figures/discovery_quoted_extents.png`.
- **Captures:** `sources/web/discovery_qhyccd_pr.md`, `discovery_astrobin_forum_126958.md`, `discovery_astrobin_images.md`, `discovery_nsns_gem01.md`, `discovery_cloudynights_925409.md`, `discovery_fesen2026_excerpts.md`, `discovery_saydjari2026_desi_skyfibres.md`, `discovery_aj_vs_arxiv.md`, `discovery_chinese_sources.md`, `discovery_search_log.md`.
- **Raw downloads:**
  - `data/raw/discovery_web/`: HTML, JSON API responses, QHYCCD and AstroBin images, GitHub diffs.
  - `data/raw/discovery_arxiv/`: Fesen v2 PDF/HTML, Saydjari PDF, arXiv abstract pages.
  - `data/raw/discovery_biblio/`: Crossref, OpenAlex, DataCite, arXiv API, AJ text extract.
  - `data/raw/discovery_vizier/`: HIP 37797.
